""""
Main analysis class.
"""

import logging

import numpy as np
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from IPython.display import display
from regions import CircleSkyRegion
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
    SpectrumDatasetOnOff,
)
from gammapy.estimators import FluxPointsEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    ExpCutoffPowerLawSpectralModel,
    SkyModel,
    create_crab_spectral_model,
)


class Analysis:
    """
    Analysis class.

    Parameters
    ----------
    configuration : str
        Analysis configuration
    output_dir : str
        Output directory
    target : SkyCoord
        Target coordinates
    data : Data
        Datastore object

    """

    def __init__(
            self,
            configuration=None,
            output_dir=None,
            target=None,
            data=None):

        self._logger = logging.getLogger(__name__)

        self._configuration = configuration
        self._output_dir = output_dir
        self._target = target
        self._data = data

        self.on_region = None
        self.exclusion_regions = []
        self.exclusion_mask = None
        self.target_exclusion_radius = 0.5 * u.deg
        self.datasets = None

    def run(self):
        """
        Run analysis.

        """

        self._define_target_region()
        self._define_exclusion_regions()
        self._data_reduction()
        self._spectral_fits(Datasets(self.datasets).stack_reduce())

    def _define_target_region(self):
        """
        Define target region.

        TODO - get on radius from fits file.

        """

        on_region_radius = Angle("0.08944272 deg")
        self.on_region = CircleSkyRegion(center=self._target, radius=on_region_radius)
        self._logger.info("Target region: ra=%.2f deg, dec=%.2f deg, radius=%.3f deg",
                          self.on_region.center.ra.deg,
                          self.on_region.center.dec.deg,
                          self.on_region.radius.deg
                          )

    def _define_exclusion_regions(self):
        """
        Define exclusion regions.

        """

        # on region
        self.exclusion_regions.append(
            CircleSkyRegion(
                center=self._target,
                radius=self.target_exclusion_radius
            )
        )

        # bright stars
        # TODO

        # exclusion mask
        geom = WcsGeom.create(
            npix=(150, 150), binsz=0.05, 
            skydir=self._target.galactic, proj="TAN",
            frame="icrs"
        )

        self.exclusion_mask = ~geom.region_mask(self.exclusion_regions)
        self._logger.info(f"Number of exclusion regions: {len(self.exclusion_regions)}")

    def _data_reduction(self):
        """
        Reduce data.

        """

        energy_axis = MapAxis.from_energy_bounds(
            0.1, 40, nbin=10, per_decade=True, unit="TeV", name="energy"
        )
        energy_axis_true = MapAxis.from_energy_bounds(
            0.05, 100, nbin=20, per_decade=True, unit="TeV", name="energy_true"
        )

        geom = RegionGeom.create(region=self.on_region, axes=[energy_axis])
        dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

        dataset_maker = SpectrumDatasetMaker(
            containment_correction=False, selection=["counts", "exposure", "edisp"]
        )
        bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=self.exclusion_mask)
        safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)
        self._logger.info("Mask applied: %s, aeff_percent = %d",
                          safe_mask_masker.methods, safe_mask_masker.aeff_percent)

        self.datasets = Datasets()

        for obs_id, observation in zip(self._data.runs, self._data.get_observations()):
            dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
            dataset_on_off = bkg_maker.run(dataset, observation)
            dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
            self.datasets.append(dataset_on_off)

        self._print_results(self.datasets.info_table(cumulative=False))

    def _print_results(self, info_table):
        """
        Print results per run.

        Parameters
        ----------
        info_table : Table
            Table with results per run

        """

        for col in info_table.colnames:
            if info_table[col].dtype.kind == 'f':
                info_table[col].format = '{:.2f}'

        print(
            info_table[
                "name", "livetime", "ontime", "counts", "counts_off", "background",
                "alpha", "excess", "sqrt_ts"
                ]
        )

    def _spectral_fits(self, _datasets):
        """
        Perform spectral fits.

        """

        _datasets.models = self._define_spectral_models()

        _fit = Fit()
        result_joint = _fit.run(datasets=[_datasets])
        print(result_joint)
        display(result_joint.models.to_parameters_table())

    def _define_spectral_models(self, model):
        """"
        Spectral models

        """

        _spectral_model=None
        if model == "pl":
            _spectral_model = PowerLawSpectralModel(
                amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
                index=2,
                reference=1 * u.TeV,
            )
        elif model == "ecpl":
            _spectral_model = ExpCutoffPowerLawSpectralModel(
                amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
                index=2,
                lambda_=0.1 * u.Unit("TeV-1"),
                reference=1 * u.TeV,
            )

        return [
            SkyModel(spectral_model=_spectral_model, name=model),
        ]
