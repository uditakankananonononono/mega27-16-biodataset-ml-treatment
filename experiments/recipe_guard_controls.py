"""Synthetic-only recipe identity guard and locked control evaluation."""
import hashlib,json,os
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
ROOT=Path(__file__).resolve().parents[1]
FIELDS=('feature_space','collapse','filter','penalty','scaling','fold_policy','rng_policy')
CONTRACT=dict(feature_space='40 synthetic independent normal features',collapse='none',filter='none',penalty='l2 C=.5 liblinear max_iter500',scaling='within-dataset zscore',fold_policy='5fold stratified shuffled seed0 rebuilt from current labels',rng_policy='bootstrap0; permutation1000+j')
def require_match(observed,null):
    bad=[k for k in FIELDS if k not in observed or k not in null or observed[k]!=null[k]]
    if bad: raise ValueError('Unmatched recipe: '+','.join(bad))
    return hashlib.sha256(json.dumps(observed,sort_keys=True).encode()).hexdigest()
def stable_count(X,y):
    X=(X-X.mean(0))/(X.std(0)+1e-8)
    folds=list(StratifiedKFold(5,shuffle=True,random_state=0).split(X,y))
    rng=np.random.RandomState(0); coefs=[]
    for b in range(40):
        tr=folds[b%5][0];sub=rng.choice(tr,len(tr),replace=True)
        clf=LogisticRegression(penalty='l2',C=.5,solver='liblinear',max_iter=500).fit(X[sub],y[sub])
        coefs.append(clf.coef_[0])
    signs=np.sign(np.array(coefs))
    # Majority sign fraction: maximum of counts for -1,0,+1. No fitted thresholds.
    consistency=np.maximum.reduce([(signs==s).mean(0) for s in (-1,0,1)])
    return int((consistency>=.975).sum())
def generate(i):
    rng=np.random.RandomState(7100+i);X=rng.normal(size=(160,40))
    if i<5:
        score=X[:,:8].mean(1)+rng.normal(0,.10,160);y=(score>np.median(score)).astype(int)
    else:y=rng.permutation(np.r_[np.ones(80,dtype=int),np.zeros(80,dtype=int)])
    return X,y
def main():
    out=ROOT/'results/recipe_guard_synthetic';out.mkdir(parents=True,exist_ok=True)
    rejected=[]
    for field in FIELDS:
        wrong=dict(CONTRACT);wrong[field]='deliberate mismatch'
        try:require_match(CONTRACT,wrong)
        except ValueError:rejected.append(field)
    recipe_hash=require_match(CONTRACT,dict(CONTRACT))
    lock_hash=hashlib.sha256((ROOT/'PREREG_RECIPE_GUARD_SYNTHETIC.md').read_bytes()).hexdigest()
    rows=[]
    for i in range(10):
        path=out/f'case{i}.json'
        if path.exists():raise RuntimeError(f'Refusing case reselection or overwrite: {path}')
        X,y=generate(i);obs=stable_count(X,y)
        null=[stable_count(X,np.random.RandomState(1000+j).permutation(y)) for j in range(39)]
        p=(1+sum(v>=obs for v in null))/40
        row=dict(case=i,kind='positive' if i<5 else 'null',observed=obs,null_counts=null,p=p,recipe_sha256=recipe_hash,protocol_sha256=lock_hash)
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(row,indent=2));os.replace(tmp,path)
        rows.append(row);print(json.dumps({k:row[k] for k in ('case','kind','observed','p')}),flush=True)
    pos=sum(r['p']<=.05 for r in rows[:5]);false=sum(r['p']<=.05 for r in rows[5:])
    result=dict(protocol_sha256=lock_hash,recipe_sha256=recipe_hash,rejected_mismatch_fields=rejected,identical_contract_accepted=True,positive_detections=pos,null_detections=false,gate=len(rejected)==7 and pos>=4 and false<=1,caveat='Synthetic independent features only; not type-I calibration, real-cohort validation, clinical novelty or benchmark superiority')
    (out/'VERDICT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()
