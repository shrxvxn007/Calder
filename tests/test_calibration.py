import numpy as np
import pytest

from calder.models.calibration import MonotoneCalibrator


def test_small_calibration_set_falls_back_without_fabrication():
    p = np.array([0.1, 0.9])
    c = MonotoneCalibrator()
    with pytest.warns(RuntimeWarning, match="raw probabilities"):
        c.fit(p, np.array([0, 1]))
    assert not c.is_fitted
    np.testing.assert_array_equal(c.predict(p), p)


def test_isotonic_is_monotone_and_clipped():
    p = np.linspace(0.1, 0.9, 20)
    labels = (p > 0.6).astype(int)
    c = MonotoneCalibrator().fit(p, labels)
    assert c.is_fitted
    out = c.predict(np.linspace(0, 1, 101))
    assert np.all(np.diff(out) >= 0)
    assert np.min(out) == 0 and np.max(out) == 1


def test_invalid_calibration_labels_and_single_class_fallback():
    with pytest.raises(ValueError):
        MonotoneCalibrator().fit(np.array([0.5]), np.array([3]))
    with pytest.warns(RuntimeWarning):
        c = MonotoneCalibrator().fit(np.linspace(0, 1, 20), np.ones(20))
    assert not c.is_fitted
