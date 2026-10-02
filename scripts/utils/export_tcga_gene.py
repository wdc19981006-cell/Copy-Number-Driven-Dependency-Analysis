"""Export only one requested gene; keep current/reference provenance explicit."""
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT
from scripts.data_access import get_tcga_cn,get_tcga_gistic,get_tcga_expression,get_tcga_metadata

def main():
 p=argparse.ArgumentParser();p.add_argument('--gene',required=True);p.add_argument('--dataset',choices=['cn','gistic','expression'],required=True);p.add_argument('--layer',choices=['current','reference'],required=True);p.add_argument('--output',required=True);a=p.parse_args()
 functions={'cn':(get_tcga_cn,'CopyNumber'),'gistic':(get_tcga_gistic,'GISTIC'),'expression':(get_tcga_expression,'Expression')}
 fn,column=functions[a.dataset];frame=fn(a.gene,layer=a.layer).rename(columns={column:'Value'})
 if a.layer=='reference':
  meta=get_tcga_metadata(layer='reference')[['SampleID','CancerType','SampleType','TumorNormal']]
  frame=frame.merge(meta,on='SampleID',how='left',validate='one_to_one')
 else:frame['CancerType']=frame.ProjectID.str.replace('TCGA-','',regex=False)
 target=Path(a.output).resolve()
 if not target.is_relative_to((ROOT/'results').resolve()):raise ValueError('Output must be within project results')
 target.parent.mkdir(parents=True,exist_ok=True);frame.to_csv(target,index=False)

if __name__=='__main__':main()
