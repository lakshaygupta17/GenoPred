import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('final_robustness',Path(__file__).parents[1]/'scripts/final_robustness.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

@pytest.mark.parametrize('model,opt',[('nnpu',{'prior':.003,'l2':.01}),('gradient_boosting',{'level':0})])
def test_heldout_features_and_labels_cannot_change_fit(model,opt):
    rng=np.random.default_rng(5)
    x=np.column_stack([rng.normal(size=(60,10)),np.zeros((60,10))]);y=np.full(60,-1);y[[0,3,8,15,45,58]]=1
    mask=np.arange(60)<40;names=[str(i) for i in range(20)]
    a=mod.fit(x,y,mask,names,model,opt,41)
    changed=x.copy();changed[~mask,:10]=1e12
    changed_y=y.copy();changed_y[~mask]=1
    b=mod.fit(changed,changed_y,mask,names,model,opt,41)
    np.testing.assert_array_equal(mod.score(a,x,names),mod.score(b,x,names))

def test_constant_scores_do_not_get_arbitrary_top_k_credit():
    y=np.full(600,-1);y[[0,12,45,599]]=1
    m=mod.metric(y,np.zeros(600))
    assert m['recall_100']==pytest.approx(100/600)
    assert m['recall_500']==pytest.approx(500/600)
    assert m['ap']==pytest.approx(4/600)
