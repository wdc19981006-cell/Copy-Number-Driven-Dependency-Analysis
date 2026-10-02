"""GDC-specific mapping and workflow safety checks using small synthetic metadata."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.download.discover_gdc import select_cn,select_segments
from scripts.utils.gdc import target_entities


def row(workflow,id,sample='sample-1',aliquot='aliquot-1'):
    return {'workflow':workflow,'file_id':id,'file_name':id+'.tsv','file_size':10,'md5sum':'x',
            'sample_id':sample,'sample':'TCGA-AA-0001-01A','case_id':'case-1','case':'TCGA-AA-0001',
            'project_id':'TCGA-BRCA','aliquot_id':aliquot,'aliquot':aliquot,'sample_type':'Primary Tumor','state':'released'}


class SelectionTests(unittest.TestCase):
    def test_workflow_priority_is_per_sample_not_per_patient(self):
        rows=[row('ASCAT2','a'),row('ASCAT3','b'),row('ABSOLUTE LiftOver','c'),
              row('ASCAT3','d','sample-2')]
        selected=select_cn(rows)
        self.assertEqual(len(selected),2)
        self.assertEqual(selected[0]['selected_workflow'],'ABSOLUTE LiftOver')
        self.assertEqual(selected[0]['available_workflows'],'ABSOLUTE LiftOver;ASCAT3;ASCAT2')

    def test_unknown_workflow_stops(self):
        with self.assertRaisesRegex(ValueError,'Unknown CN workflow'):
            select_cn([row('New Pipeline','a')])

    def test_segment_exact_match_then_documented_fallback(self):
        genes=select_cn([row('ABSOLUTE LiftOver','g')])
        selected,audit=select_segments([row('ASCAT2','s2'),row('ASCAT3','s3')],genes)
        self.assertEqual(selected[0]['workflow'],'ASCAT3')
        self.assertFalse(selected[0]['workflow_match'])
        self.assertIn('not published',audit[0]['segment_selection_reason'])
        genes=select_cn([row('ASCAT2','g')])
        selected,_=select_segments([row('ASCAT2','s2'),row('ASCAT3','s3')],genes)
        self.assertEqual(selected[0]['workflow'],'ASCAT2')

    def test_filename_aliquot_identifies_tumor_not_matched_normal(self):
        def sample(id,barcode,type,aliquot):
            return {'sample_id':id,'submitter_id':barcode,'sample_type':type,
                    'portions':[{'analytes':[{'aliquots':[{'aliquot_id':aliquot,'submitter_id':aliquot}]}]}]}
        file={'file_name':'TCGA-BRCA.tumor-uuid.gene_level_copy_number.tsv','data_type':'Gene Level Copy Number',
              'cases':[{'case_id':'case','submitter_id':'TCGA-AA-0001','project':{'project_id':'TCGA-BRCA'},
                        'samples':[sample('normal','TCGA-AA-0001-10A','Blood Derived Normal','normal-uuid'),
                                   sample('tumor','TCGA-AA-0001-01A','Primary Tumor','tumor-uuid')]}]}
        mapped=target_entities(file)
        self.assertEqual(len(mapped),1)
        self.assertEqual(mapped[0]['sample_id'],'tumor')
        self.assertEqual(mapped[0]['aliquot_id'],'tumor-uuid')


if __name__=='__main__':
    unittest.main()
