import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('final_extension',Path(__file__).parents[1]/'scripts/final_extension.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

@pytest.mark.parametrize('group',list(mod.GROUPS))
def test_feature_groups_and_heldout_data_are_respected(group):
    rng=np.random.default_rng(20)
    x=np.column_stack([rng.normal(size=(60,10)),np.zeros((60,10))]);y=np.full(60,-1);y[[0,4,9,18,44]]=1
    mask=np.arange(60)<40;opt={'prior':.003,'l2':.01}
    a=mod.fit(x,y,mask,'nnpu',group,opt,15)
    changed=x.copy();changed[~mask,:10]=1e10;labels=y.copy();labels[~mask]=1
    b=mod.fit(changed,labels,mask,'nnpu',group,opt,15)
    np.testing.assert_array_equal(mod.score(a,x),mod.score(b,x))
    assert len(a['w'])==len(mod.GROUPS[group])

def test_mutation_model_does_not_use_other_assays():
    rng=np.random.default_rng(5)
    x=np.column_stack([rng.normal(size=(60,10)),np.zeros((60,10))]);y=np.full(60,-1);y[:6]=1;mask=np.arange(60)<40
    changed=x.copy();changed[:,[0,1,2,3,4,5,7,8,9]]=rng.normal(100,20,size=(60,9))
    opt={'prior':.003,'l2':.01}
    a=mod.fit(x,y,mask,'nnpu','mutation',opt,15);b=mod.fit(changed,y,mask,'nnpu','mutation',opt,15)
    np.testing.assert_array_equal(mod.score(a,x),mod.score(b,changed))

def test_tie_weight_preserves_requested_candidate_budget():
    weights=mod.top_weight(np.array([0.,1.,1.,1.,2.]),3)
    np.testing.assert_allclose(weights,[0,2/3,2/3,2/3,1])
    assert weights.sum()==pytest.approx(3)
