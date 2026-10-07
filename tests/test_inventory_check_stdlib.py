import copy
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src' / 'provenance'))
from inventory_check import inspect_inventory, nonfinite_paths

class InventoryChecks(unittest.TestCase):
    def doc(self):
        return {'tools':[{'tool':'a','status':'ok','result':{'v':1}}], 'n_tools_ok':1,
                'n_tools_attempted':1, 'code_revision':'x','data_sha256':'x',
                'package_versions':{'a':'1'},'executed_at':'x'}
    def test_clean_is_not_certification(self):
        r=inspect_inventory(self.doc());self.assertEqual(r['issues'],[]);self.assertFalse(r['certified'])
    def test_nested_nonfinite(self):
        self.assertEqual(nonfinite_paths({'a':[float('inf'),float('nan')]}),['$.a[0]','$.a[1]'])
    def test_count_mismatch(self):
        d=self.doc();d['n_tools_ok']=9;self.assertTrue(any('ok count' in x for x in inspect_inventory(d)['issues']))
    def test_duplicate(self):
        d=self.doc();d['tools']*=2;self.assertIn('duplicate tool names',inspect_inventory(d)['issues'])
    def test_missing_provenance(self):
        d=self.doc();del d['data_sha256'];self.assertTrue(any('data_sha256' in x for x in inspect_inventory(d)['issues']))
    def test_ok_no_result(self):
        d=self.doc();d['tools'][0].pop('result');self.assertTrue(any('without result' in x for x in inspect_inventory(d)['issues']))
    def test_malformed(self):
        self.assertFalse(inspect_inventory({'tools':'bad'})['certified'])
    def test_unknown_status(self):
        d=self.doc();d['tools'][0]['status']='installed';self.assertTrue(any('unknown status' in x for x in inspect_inventory(d)['issues']))

if __name__=='__main__':unittest.main()
