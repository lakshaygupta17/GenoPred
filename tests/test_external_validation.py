import numpy as np
from genopred_v2.external_validation import holm, ranks, top_ties, collapse_donors


def test_holm_restores_order_and_monotonicity():
    np.testing.assert_allclose(holm([.04,.001,.03,.8]), [.09,.004,.09,.8])


def test_negative_fitness_has_positive_expected_direction():
    score = np.array([1.,2.,3.,4.]); fitness = np.array([0.,-.1,-.5,-1.])
    assert np.isclose(ranks(score) @ ranks(-fitness), 1.)
    assert np.isclose(ranks(score) @ ranks(fitness), -1.)


def test_tied_cutoff_is_inclusive():
    assert top_ties(np.array([0.,1.,1.,2.]), 2).tolist() == [False,True,True,True]


def test_donor_collapsing_prevents_duplicate_weight():
    values = np.array([[1.,3.,9.], [0.,2.,4.]])
    np.testing.assert_allclose(collapse_donors(values, ['a','a','b']), [[2.,9.],[1.,4.]])


def test_rank_correlation_handles_score_ties():
    x = ranks(np.array([[1.,1.],[1.,2.],[3.,3.]]))
    np.testing.assert_allclose((x*x).sum(axis=0), [1.,1.])
