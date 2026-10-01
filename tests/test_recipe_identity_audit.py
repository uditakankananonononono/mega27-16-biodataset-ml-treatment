import importlib.util
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'src/replication/recipe_identity_audit.py'
s=importlib.util.spec_from_file_location('recipe_audit',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_current_recipes_are_not_claimed_identical():
 a=m.audit();assert not a['matched_recipe'];assert set(a['mismatches'])=={'symbol_collapse_present','filter'};assert a['no_training_or_test_evaluation']
def test_sources_are_hash_bound():
 assert all(len(x['sha256'])==64 for x in m.audit()['recipes'].values())
