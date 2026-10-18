"""Regression checks for reproducible point-source selections."""

from argparse import Namespace
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest

from v2dl5.configuration import configuration
from v2dl5.data import Data


def test_name_only_target_does_not_inherit_crab(tmp_path):
    path = tmp_path / "config.yml"
    path.write_text('on_region:\n  target: "M 87"\n')
    config = configuration(Namespace(config=str(path), run_list=None, output_dir=str(tmp_path)))
    assert config["on_region"].get("lon") is None
    assert config["on_region"].get("lat") is None


def test_nested_target_is_rejected(tmp_path):
    path = tmp_path / "config.yml"
    path.write_text('datasets:\n  on_region:\n    target: "M 87"\n')
    with pytest.raises(ValueError, match="on_region"):
        configuration(Namespace(config=str(path), run_list=None, output_dir=str(tmp_path)))


@pytest.mark.parametrize("content", ["7\n7\n", "7 missing VPM\n", "# empty\n"])
def test_invalid_run_lists_are_rejected(tmp_path, content):
    path = tmp_path / "runs.txt"
    path.write_text(content)
    data = Data.__new__(Data)
    with pytest.raises(ValueError, match=r"[Rr]un|Duplicate|Invalid"):
        data._from_run_list(path)


def test_energy_dependent_radius_is_rejected():
    data = Data.__new__(Data)
    data.get_observations = lambda: [
        SimpleNamespace(rad_max=SimpleNamespace(quantity=np.array([[0.1], [0.2]]) * u.deg))
    ]
    with pytest.raises(ValueError, match="constant"):
        data.get_on_region_radius()


def test_reduced_spectrum_reload_preserves_fit_and_safe_masks(tmp_path):
    from gammapy.datasets import Datasets, SpectrumDatasetOnOff
    from gammapy.maps import MapAxis, RegionGeom

    from v2dl5.analysis import read_datasets

    axis = MapAxis.from_energy_bounds("0.1 TeV", "10 TeV", nbin=10)
    geom = RegionGeom.create("icrs;circle(83.63,22.01,0.1)", axes=[axis])
    dataset = SpectrumDatasetOnOff.create(
        geom, energy_axis_true=axis.copy(name="energy_true"), name="7"
    )
    dataset.mask_safe.data[:] = True
    dataset.mask_fit = geom.energy_mask(energy_min=0.3 * u.TeV, energy_max=3 * u.TeV)
    dataset.acceptance.data[:] = 1
    dataset.acceptance_off.data[:] = 5
    dataset.counts.data[:] = 3
    dataset.counts_off.data[:] = 10
    Datasets([dataset]).write(tmp_path / "datasets.yaml")

    loaded = read_datasets(tmp_path / "datasets.yaml")[0]
    np.testing.assert_array_equal(loaded.mask.data, dataset.mask.data)
    np.testing.assert_array_equal(loaded.counts_off.data, dataset.counts_off.data)
    np.testing.assert_array_equal(loaded.acceptance_off.data, dataset.acceptance_off.data)
    edges = axis.edges
    active = loaded.mask.data[:, 0, 0]
    assert np.all(edges[:-1][active] >= 0.3 * u.TeV)
    assert np.all(edges[1:][active] <= 3 * u.TeV)


def test_bti_filters_events_and_livetime_consistently():
    import logging

    from astropy.table import Table
    from astropy.time import Time
    from gammapy.data import GTI, EventList, Observation

    reference = Time("2023-01-01", scale="tt")
    gti = GTI.create([0] * u.s, [10] * u.s, reference_time=reference)
    events = EventList(
        Table(
            {
                "TIME": reference + np.array([1, 3, 5, 7, 9]) * u.s,
                "RA": np.zeros(5) * u.deg,
                "DEC": np.zeros(5) * u.deg,
                "ENERGY": np.ones(5) * u.TeV,
            }
        )
    )
    obs = Observation(obs_id=7, events=events, gti=gti)
    data = Data.__new__(Data)
    data._bti = [{"run": 7, "bti_start": 2.5, "bti_length": 3}]
    data._logger = logging.getLogger(__name__)
    selected = [obs]
    data._apply_bti(selected)
    obs = selected[0]
    assert obs.gti.time_sum.to_value("s") == pytest.approx(7)
    np.testing.assert_allclose((obs.events.time - reference).to_value("s"), [1, 7, 9], atol=1e-6)


def test_weak_source_has_finite_nonnegative_flux_limit():
    from gammapy.datasets import Datasets, SpectrumDatasetOnOff
    from gammapy.maps import MapAxis, RegionGeom

    from v2dl5.analysis import Analysis

    axis = MapAxis.from_energy_bounds("0.1 TeV", "20 TeV", nbin=20)
    geom = RegionGeom.create("icrs;circle(100,20,0.1)", axes=[axis])
    dataset = SpectrumDatasetOnOff.create(
        geom, energy_axis_true=axis.copy(name="energy_true"), name="weak"
    )
    dataset.exposure.data[:] = 1e8
    dataset.acceptance.data[:] = 1
    dataset.acceptance_off.data[:] = 5
    dataset.counts.data[:] = 1
    dataset.counts_off.data[:] = 20
    dataset.mask_safe.data[:] = True
    dataset.mask_fit = geom.energy_mask(energy_min=1 * u.TeV, energy_max=20 * u.TeV)
    analysis = Analysis(args_dict={"output_dir": ".", "integral_flux": {
        "energy": {"min": "2 TeV", "max": "10 TeV"}, "index": 2.5}})
    analysis._integral_flux(Datasets([dataset]))
    result = analysis.integral_flux.to_table(sed_type="flux")[0]
    assert result["success"]
    assert 0 <= result["flux"] < result["flux_ul"]
    assert np.isfinite(result["flux_ul"])
