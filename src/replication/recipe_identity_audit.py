"""Static recipe-identity audit. No labels, fits, thresholds or test data are changed."""
import ast,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
SOURCES={'observed':ROOT/'src/replication/twocohort_selection.py','permutation':ROOT/'experiments/permutation_arm.py'}
def audit():
 records={}
 for label,p in SOURCES.items():
  text=p.read_text();tree=ast.parse(text)
  calls=sorted({getattr(n.func,'id',getattr(n.func,'attr','')) for n in ast.walk(tree) if isinstance(n,ast.Call)})
  records[label]={'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'calls':calls,'symbol_collapse_present':'p2s' in text and 'best' in text,'filter':'de_genes' if 'de_genes(' in text else 'tstats_chunked' if 'tstats_chunked(' in text else 'unknown'}
 mismatch=[]
 for field in ['symbol_collapse_present','filter']:
  if records['observed'][field]!=records['permutation'][field]:mismatch.append(field)
 return {'scope':'static source inspection, not scientific null validation','recipes':records,'matched_recipe':not mismatch,'mismatches':mismatch,'manual_review_needed':['Compare numerical filtering equivalence and RNG/fold construction','Observed source docstring says L1 but LogisticRegression omits penalty and therefore uses library default; inspect installed version and effective params','Do not infer a matched permutation p-value from unmatched feature spaces'],'no_training_or_test_evaluation':True}
if __name__=='__main__':print(json.dumps(audit(),indent=2))
