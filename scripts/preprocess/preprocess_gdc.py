"""Bounded DR46 ETL: samples x genes, ZSTD Parquet, no specific gene analysis."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,setup,write_json,now,hashes
from scripts.utils.gdc import RAW,PROCESSED,local_path,init_dirs,runtime_path,manifest_path
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb

VERSION='gdc-DR46-v1'
IDMAP={'SampleID':'sample_id','SampleBarcode':'sample','AliquotID':'aliquot_id',
       'AliquotBarcode':'aliquot','CaseID':'case_id','CaseBarcode':'case',
       'ProjectID':'project_id','SampleType':'sample_type','FileID':'file_id'}


def selections(kind):
    return list({r['file_id']:r for r in json.loads((RAW/f'metadata/{kind}_selected.json').read_text(encoding='utf-8'))}.values())


def verified_inputs(kind):
    rows=selections(kind)
    manifest=manifest_path()
    if not manifest.exists():
        raise FileNotFoundError('Complete download manifest is not yet available')
    with manifest.open(encoding='utf-8',newline='') as f:
        verified={r['file_id']:r for r in csv.DictReader(f) if r['status']=='verified'}
    failures=[]
    for row in rows:
        path=local_path(kind,row)
        if row['file_id'] not in verified or not path.exists() or path.stat().st_size!=row['file_size']:
            failures.append(row['file_id'])
    if failures:
        raise ValueError(f'{kind}: {len(failures)} unverified/missing inputs. Complete download first.')
    stamp=hashlib.sha256(''.join(r['file_id']+verified[r['file_id']]['SHA256'] for r in rows).encode()).hexdigest()
    return rows,stamp


def ready(dataset,stamp):
    path=PROCESSED/f'{dataset}.parquet'
    meta=path.with_suffix('.provenance.json')
    if not path.exists() or not meta.exists(): return False
    obj=json.loads(meta.read_text(encoding='utf-8'))
    return obj.get('source_set_sha256')==stamp and obj.get('processor_version')==VERSION and hashes(path)[0]==obj.get('SHA256')


def register(dataset,stamp,extra):
    path=PROCESSED/f'{dataset}.parquet'
    meta={'created_at':now(),'release':'46.0','processor_version':VERSION,
          'source_set_sha256':stamp,'SHA256':hashes(path)[0],'bytes':path.stat().st_size,**extra}
    write_json(path.with_suffix('.provenance.json'),meta)
    return meta


def safe_temp_remove(path):
    # Only generated memory-map temporary files within this fixed project directory.
    path=Path(path).resolve()
    permitted=(PROCESSED/'_work').resolve()
    if not path.is_relative_to(permitted) or path.suffix!='.mmap':
        raise ValueError('Unsafe temporary cleanup target')
    path.unlink(missing_ok=True)


def write_wide(dataset,rows,genes,values,stamp,transform=None,workflow=False):
    if ready(dataset,stamp):
        logging.info('Verified derived dataset; skip %s',dataset)
        return json.loads((PROCESSED/f'{dataset}.provenance.json').read_text(encoding='utf-8'))
    ids=list(IDMAP)
    if workflow: ids+=['Workflow']
    schema=pa.schema([(c,pa.string()) for c in ids]+[(g,pa.float64()) for g in genes.gene_id],
             metadata={b'release':b'GDC DR46.0',b'orientation':b'samples/aliquots x genes',
                       b'gene_id_version':b'original Ensembl IDs; symbols in genes sidecar',
                       b'source_set_sha256':stamp.encode()})
    target=PROCESSED/f'{dataset}.parquet'
    temp=target.with_suffix('.parquet.part')
    missing=0
    with pq.ParquetWriter(temp,schema,compression='zstd',compression_level=6) as writer:
        for begin in range(0,len(rows),1024):
            end=min(len(rows),begin+1024)
            block=np.asarray(values[begin:end,:len(genes)])
            if transform: block=transform(block)
            missing+=int(np.isnan(block).sum())
            arrays=[pa.array([row.get(IDMAP[c],None) for row in rows[begin:end]],type=pa.string()) for c in IDMAP]
            if workflow: arrays.append(pa.array([r['workflow'] for r in rows[begin:end]],type=pa.string()))
            arrays.extend(pa.array(block[:,j],from_pandas=True) for j in range(len(genes)))
            table=pa.Table.from_arrays(arrays,schema=schema)
            writer.write_table(table,row_group_size=1024)
            del table,arrays,block
            logging.info('%s Parquet rows %s/%s',dataset,end,len(rows))
    temp.replace(target)
    pq.write_table(pa.Table.from_pandas(genes,preserve_index=False),PROCESSED/f'{dataset}.genes.parquet',compression='zstd')
    return register(dataset,stamp,{'rows':len(rows),'genes':len(genes),'missing_cells':missing,
          'compression':'ZSTD','orientation':'samples/aliquots x genes','parsed_all_inputs':True})


def rna_frame(path):
    frame=pd.read_csv(path,sep='\t',comment='#',dtype={'gene_id':str,'gene_name':str,'gene_type':str})
    required={'gene_id','gene_name','gene_type','unstranded','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded'}
    if not required.issubset(frame): raise ValueError('Unexpected STAR schema')
    genes=frame.loc[frame.gene_id.str.startswith('ENSG')].copy()
    summary=frame.loc[~frame.gene_id.str.startswith('ENSG')].copy()
    if genes.gene_id.duplicated().any(): raise ValueError('Duplicate STAR gene IDs')
    return genes,summary


def expression():
    rows,stamp=verified_inputs('rna')
    names={'unstranded':'TCGA_STAR_UnstrandedCounts','tpm_unstranded':'TCGA_STAR_TPM',
           'fpkm_unstranded':'TCGA_STAR_FPKM','fpkm_uq_unstranded':'TCGA_STAR_FPKMUQ'}
    if all(ready(n,stamp) for n in list(names.values())+['TCGA_STAR_log2TPMplus1']):
        logging.info('All expression datasets already verified')
        return
    first,_=rna_frame(local_path('rna',rows[0]))
    genes=first[['gene_id','gene_name','gene_type']].reset_index(drop=True)
    baseline=genes.gene_id.to_numpy()
    work=PROCESSED/'_work';work.mkdir(parents=True,exist_ok=True)
    mmaps={field:np.memmap(work/f'{field}.mmap',mode='w+',dtype=np.float64,shape=(len(rows),len(genes))) for field in names}
    stats=[]
    for index,row in enumerate(rows):
        frame,summary=rna_frame(local_path('rna',row))
        if not np.array_equal(frame.gene_id.to_numpy(),baseline):
            raise ValueError(f'Inconsistent STAR annotation/gene order in {row["file_id"]}; no genes silently dropped')
        for field,values in mmaps.items():
            vector=frame[field].to_numpy(dtype=np.float64)
            if not np.isfinite(vector).all() or (vector<0).any():
                raise ValueError(f'STAR nonfinite/negative metric {field} in {row["file_id"]}')
            if field=='unstranded' and not np.equal(vector,np.floor(vector)).all():
                raise ValueError('Non-integer raw count')
            values[index,:]=vector
        stats.append({'FileID':row['file_id'],**{r.gene_id:float(r.unstranded) for r in summary.itertuples()}})
        if (index+1)%250==0: logging.info('STAR parsed %s/%s files',index+1,len(rows))
    for array in mmaps.values(): array.flush()
    pq.write_table(pa.Table.from_pylist(stats),PROCESSED/'TCGA_STAR_Library_Summary.parquet',compression='zstd')
    result={}
    for field,dataset in names.items():
        result[dataset]=write_wide(dataset,rows,genes,mmaps[field],stamp)
        if field=='tpm_unstranded':
            result['TCGA_STAR_log2TPMplus1']=write_wide('TCGA_STAR_log2TPMplus1',rows,genes,mmaps[field],stamp,
                                                     transform=lambda x:np.log2(x+1))
        filename=Path(mmaps[field].filename)
        mmaps[field]._mmap.close()
        safe_temp_remove(filename)
    write_json(runtime_path('gdc_expression_qc.json'),{'created_at':now(),'files':len(rows),
        'genes':len(genes),'unique_samples':len({r['sample_id'] for r in rows}),
        'unique_aliquots':len({r['aliquot_id'] for r in rows}),'datasets':result})


def cn_frame(path):
    frame=pd.read_csv(path,sep='\t',dtype={'gene_id':str,'gene_name':str},comment='#')
    if not {'gene_id','gene_name','copy_number'}.issubset(frame): raise ValueError('Unexpected gene CN schema')
    if frame.gene_id.duplicated().any(): raise ValueError('Duplicate CN gene IDs')
    return frame


def copy_number():
    rows,stamp=verified_inputs('cn')
    if ready('TCGA_GeneLevel_CN',stamp): return
    work=PROCESSED/'_work';work.mkdir(parents=True,exist_ok=True)
    capacity=65000
    path=work/'cn.mmap'
    matrix=np.memmap(path,mode='w+',dtype=np.float64,shape=(len(rows),capacity))
    matrix[:]=np.nan
    indices,annotations={},[]
    baseline=None
    baseline_order=None
    for i,row in enumerate(rows):
        frame=cn_frame(local_path('cn',row))
        if baseline is not None and np.array_equal(frame.gene_id.to_numpy(),baseline):
            order=baseline_order
        else:
            order=[]
            for item in frame[['gene_id','gene_name','chromosome','start','end']].itertuples(index=False,name=None):
                id=item[0]
                if id not in indices:
                    if len(indices)>=capacity: raise ValueError('Gene union exceeds budgeted capacity; re-estimate space before expanding')
                    indices[id]=len(indices)
                    annotations.append(dict(zip(['gene_id','gene_name','chromosome','start','end'],item)))
                order.append(indices[id])
            order=np.asarray(order)
            if baseline is None:
                baseline=frame.gene_id.to_numpy()
                baseline_order=order
        vector=pd.to_numeric(frame.copy_number,errors='raise').to_numpy(dtype=np.float64)
        if np.isinf(vector).any() or (vector[~np.isnan(vector)]<0).any(): raise ValueError('Invalid copy-number values')
        matrix[i,order]=vector
        if (i+1)%250==0: logging.info('CN parsed %s/%s files',i+1,len(rows))
    matrix.flush()
    genes=pd.DataFrame(annotations)
    result=write_wide('TCGA_GeneLevel_CN',rows,genes,matrix,stamp,workflow=True)
    matrix._mmap.close();safe_temp_remove(path)
    write_json(runtime_path('gdc_cn_qc.json'),{'created_at':now(),'files':len(rows),'genes':len(genes),
        'unique_samples':len({r['sample_id'] for r in rows}),'datasets':{'TCGA_GeneLevel_CN':result}})


def source_columns(path):
    opener=gzip.open if path.name.endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8-sig') as f:
        for line in f:
            if not line.startswith('#'):
                return line.rstrip('\r\n').split('\t')
    raise ValueError('No source header')


def long_table(kind,dataset):
    rows,stamp=verified_inputs(kind)
    if ready(dataset,stamp): return
    original=list(dict.fromkeys(c for row in rows for c in source_columns(local_path(kind,row))))
    standardized=['project_id','case_id','case_barcode','sample_id','sample_barcode','aliquot_id','aliquot_barcode','file_id','workflow']
    if kind=='segments': standardized+=['chromosome','start','end','copy_number','major_copy_number','minor_copy_number','purity','ploidy','gene_cn_workflow','workflow_match']
    else: standardized+=['gene','chromosome','position','ref','alt','variant_classification','variant_type','protein_change','VAF','identity_mapping_method']
    identity_map={};identity_audit=[]
    if kind=='mutation':
        bio=pq.read_table(PROCESSED/'TCGA_Biospecimen.parquet',columns=['project_id','case_id','case_submitter_id','sample_id','sample_submitter_id','sample_type','aliquot_id','aliquot_submitter_id']).to_pylist()
        for entry in bio:
            if entry['aliquot_id']:
                if entry['aliquot_id'] in identity_map and identity_map[entry['aliquot_id']]!=entry:raise ValueError('Ambiguous biospecimen UUID')
                identity_map[entry['aliquot_id']]=entry
    columns=list(dict.fromkeys(standardized+original))
    numeric={'start':pa.int64(),'end':pa.int64(),'position':pa.int64(),'copy_number':pa.float64(),
             'major_copy_number':pa.float64(),'minor_copy_number':pa.float64(),'purity':pa.float64(),'ploidy':pa.float64(),'VAF':pa.float64()}
    schema=pa.schema([(c,numeric.get(c,pa.string())) for c in columns],metadata={b'release':b'GDC DR46.0'})
    target=PROCESSED/f'{dataset}.parquet';temp=target.with_suffix('.parquet.part')
    records=0
    with pq.ParquetWriter(temp,schema,compression='zstd',compression_level=6) as writer:
        for number,row in enumerate(rows,1):
            path=local_path(kind,row)
            for frame in pd.read_csv(path,sep='\t',comment='#',dtype=str,chunksize=16384,low_memory=False):
                if frame.empty: continue
                data={c:frame[c].tolist() if c in frame else [None]*len(frame) for c in original}
                const={'project_id':row.get('project_id'),'case_id':row.get('case_id'),'case_barcode':row.get('case'),
                       'sample_id':row.get('sample_id'),'sample_barcode':row.get('sample'),'aliquot_id':row.get('aliquot_id'),
                       'aliquot_barcode':row.get('aliquot'),'file_id':row['file_id'],'workflow':row['workflow']}
                data.update({k:[v]*len(frame) for k,v in const.items()})
                if kind=='segments':
                    if 'GDC_Aliquot' in frame and not (frame.GDC_Aliquot==row['aliquot_id']).all():
                        raise ValueError('Segment source aliquot differs from selected metadata')
                    mapping={'chromosome':'Chromosome','start':'Start','end':'End','copy_number':'Copy_Number',
                             'major_copy_number':'Major_Copy_Number','minor_copy_number':'Minor_Copy_Number'}
                    for dest,src in mapping.items():
                        data[dest]=frame[src].tolist() if src in frame else [None]*len(frame)
                    for c in ['purity','ploidy']: data[c]=[None]*len(frame)
                    data['gene_cn_workflow']=[row.get('gene_cn_workflow')]*len(frame)
                    data['workflow_match']=[str(row.get('workflow_match'))]*len(frame)
                else:
                    mapping={'gene':'Hugo_Symbol','chromosome':'Chromosome','position':'Start_Position',
                             'ref':'Reference_Allele','alt':'Tumor_Seq_Allele2','variant_classification':'Variant_Classification',
                             'variant_type':'Variant_Type','protein_change':'HGVSp_Short'}
                    for dest,src in mapping.items(): data[dest]=frame[src].tolist() if src in frame else [None]*len(frame)
                    if 'Tumor_Sample_UUID' not in frame or 'Tumor_Sample_Barcode' not in frame:
                        raise ValueError('MAF lacks source tumor identity columns')
                    mapped={}
                    for source_uuid,source_barcode in frame[['Tumor_Sample_UUID','Tumor_Sample_Barcode']].drop_duplicates().itertuples(index=False,name=None):
                        source=identity_map.get(source_uuid)
                        if not source:raise ValueError(f'MAF tumor UUID absent from full biospecimen hierarchy: {source_uuid}')
                        if source['aliquot_submitter_id']!=source_barcode or source['case_id']!=row['case_id']:
                            raise ValueError(f'MAF barcode/case mismatch with source tumor biospecimen: {source_uuid}')
                        mapped[source_uuid]=source
                        if source_uuid!=row.get('aliquot_id'):
                            identity_audit.append({'file_id':row['file_id'],'api_linked_aliquot':row.get('aliquot_id'),'source_tumor_aliquot':source_uuid,'source_tumor_barcode':source_barcode,'case_id':row['case_id'],'method':'Source Tumor_Sample_UUID resolved against full official case biospecimen hierarchy'})
                    fields={'project_id':'project_id','case_id':'case_id','case_barcode':'case_submitter_id','sample_id':'sample_id','sample_barcode':'sample_submitter_id','aliquot_id':'aliquot_id','aliquot_barcode':'aliquot_submitter_id'}
                    for dest,src in fields.items():data[dest]=[mapped[v][src] for v in frame.Tumor_Sample_UUID]
                    data['identity_mapping_method']=['source MAF tumor UUID + barcode; full official case biospecimen hierarchy']*len(frame)
                    if 't_depth' in frame and 't_alt_count' in frame:
                        depth=pd.to_numeric(frame.t_depth,errors='raise')
                        alt=pd.to_numeric(frame.t_alt_count,errors='raise')
                        data['VAF']=(alt/depth.where(depth>0)).tolist()
                    else: data['VAF']=[None]*len(frame)
                arrays=[]
                for c in columns:
                    values=data.get(c,[None]*len(frame))
                    if c in numeric:
                        series=pd.to_numeric(pd.Series(values),errors='raise')
                        arrays.append(pa.array(series,type=numeric[c],from_pandas=True))
                    else:
                        arrays.append(pa.array(values,type=pa.string(),from_pandas=True))
                writer.write_table(pa.Table.from_arrays(arrays,schema=schema),row_group_size=16384)
                records+=len(frame)
            if number%250==0: logging.info('%s %s/%s files',dataset,number,len(rows))
    temp.replace(target)
    result=register(dataset,stamp,{'rows':records,'files':len(rows),'columns':len(columns),'parsed_all_inputs':True,'compression':'ZSTD'})
    write_json(runtime_path(f'gdc_{kind}_qc.json'),result)
    if kind=='mutation':write_json(runtime_path('gdc_maf_identity_qc.json'),{'created_at':now(),'file_metadata_discrepancies':identity_audit,'note':'File API-linked entities are retained as input provenance. Canonical row tumor identity comes from source MAF UUID/barcode and the full official biospecimen hierarchy; never from a normal-only file association.'})


def database():
    target=ROOT/'data/processed/tcga/tcga_current.duckdb'
    with duckdb.connect(str(target)) as con:
        con.execute("SET memory_limit='2GB'")
        for path in PROCESSED.glob('TCGA_*.parquet'):
            if path.name.endswith('.genes.parquet'): continue
            value=path.as_posix().replace("'","''")
            name=path.stem.replace('"','""')
            con.execute(f'CREATE OR REPLACE VIEW "{name}" AS SELECT * FROM read_parquet(\'{value}\')')
    required=['TCGA_STAR_TPM','TCGA_STAR_log2TPMplus1','TCGA_GeneLevel_CN','TCGA_CN_Segments',
              'TCGA_Masked_Somatic_Mutation','TCGA_Clinical','TCGA_Biospecimen','TCGA_Sample_Map']
    complete=all((PROCESSED/f'{name}.parquet').exists() for name in required)
    layers=json.loads((ROOT/'config/tcga_layers.json').read_text(encoding='utf-8'))
    layers.update(current_status='processed_unvalidated' if complete else 'incomplete',last_preprocessed_at=now())
    write_json(runtime_path('current_layer_status.json'),layers)
    return complete


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--kinds',nargs='+',default=['rna','cn','segments','mutation'],choices=['rna','cn','segments','mutation'])
    args=parser.parse_args()
    setup('preprocess_gdc');init_dirs()
    for kind in args.kinds:
        if kind=='rna': expression()
        elif kind=='cn': copy_number()
        elif kind=='segments': long_table(kind,'TCGA_CN_Segments')
        elif kind=='mutation': long_table(kind,'TCGA_Masked_Somatic_Mutation')
    print('Current database complete:',database())


if __name__=='__main__':
    main()
