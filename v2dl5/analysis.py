""""
Main analysis class.
"""

import logging

import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from gammapy.datasets import Datasets, SpectrumDataset, FluxPointsDataset
from gammapy.estimators import FluxPointsEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
)
from IPython.display import display
from regions import CircleSkyRegion

import v2dl5.plot as v2dl5_plot


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

    def __init__(self, args_dict=None, target=None, data=None):
        self._logger = logging.getLogger(__name__)

        self.args_dict = args_dict
        self._output_dir = args_dict.get("output_dir", None)
        self._target = target
        self._data = data

        self.on_region = None
        self.exclusion_regions = []
        self.exclusion_mask = None
        self.target_exclusion_radius = 0.5 * u.deg
        self.datasets = None
        self.spectral_model = None
        self.flux_points = None

    def run(self):
        """
        Run analysis.

        """

        self._define_target_region()
        self._define_exclusion_regions()
        self._data_reduction()
        if self.args_dict["datasets"]["stack"]:
            _data_sets = Datasets(self.datasets).stack_reduce()
        else:
            _data_sets = self.datasets
        self._define_spectral_models(model=self.args_dict["fit"]["model"])
        self._spectral_fits(datasets=_data_sets)
        self._flux_points(_data_sets)

    def plot(self):
        """
        Plot all results.

        """

        for dataset in self.datasets:
            v2dl5_plot.plot_fit(dataset, self._output_dir)

        v2dl5_plot.plot_flux_points(self.flux_points, self._output_dir)
        v2dl5_plot.plot_sed(
            FluxPointsDataset(
                data=self.flux_points,
                models=self.spectral_model.copy()
            ),
            self._output_dir
        )
        return

    def write(self):
        """
        Write all results.

        """

        print("Model: ", self.spectral_model)
        print("Fluxpoints: ", self.flux_points)
        if self.flux_points is not None:
            FluxPointsDataset(
                data=self.flux_points,
                models=[self.spectral_model]
            ).write("flux_points.fits", overwrite=True)


    def _define_target_region(self):
        """
        Define target region.

        TODO - get on radius from fits file.

        """

        on_region_radius = Angle("0.08944272 deg")
        self.on_region = CircleSkyRegion(center=self._target, radius=on_region_radius)
        self._logger.info(
            "Target region: ra=%.2f deg, dec=%.2f deg, radius=%.3f deg",
            self.on_region.center.ra.deg,
            self.on_region.center.dec.deg,
            self.on_region.radius.deg,
        )

    def _define_exclusion_regions(self):
        """
        Define exclusion regions.

        """

        # on region
        self.exclusion_regions.append(
            CircleSkyRegion(center=self._target, radius=self.target_exclusion_radius)
        )

        # bright stars
        # TODO

        # exclusion mask
        geom = WcsGeom.create(
            npix=(150, 150), binsz=0.05, skydir=self._target.galactic, proj="TAN", frame="icrs"
        )

        self.exclusion_mask = ~geom.region_mask(self.exclusion_regions)
        self._logger.info("Number of exclusion regions: %d", len(self.exclusion_regions))

    def _data_reduction(self):
        """
        Reduce data using the reflected region maker.

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
        self._logger.info(
            "Mask applied: %s, aeff_percent = %d",
            safe_mask_masker.methods,
            safe_mask_masker.aeff_percent,
        )

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
            if info_table[col].dtype.kind == "f":
                info_table[col].format = "{:.2f}"

        print(
            info_table[
                "name",
                "livetime",
                "ontime",
                "counts",
                "counts_off",
                "background",
                "alpha",
                "excess",
                "sqrt_ts",
            ]
        )

    def _spectral_fits(self, datasets=None):
        """
        Perform spectral fits.

        """

        datasets.models = self.spectral_model

        _fit = Fit()
        result_joint = _fit.run(datasets=datasets)
        print(result_joint)
        display(result_joint.models.to_parameters_table())

    def _define_spectral_models(self, model=None):
        """ "
        Spectral models

        """

        _spectral_model = None
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

        self.spectral_model = SkyModel(spectral_model=_spectral_model, name=model)

        return [self.spectral_model]

    def _flux_points(self, datasets):
        """
        Calculate flux points.

        """

        e_min, e_max = 0.2, 30
        energy_edges = np.geomspace(e_min, e_max, 11) * u.TeV

        self._logger.info("Estimating flux points")
        fpe = FluxPointsEstimator(energy_edges=energy_edges, selection_optional="all")
        self.flux_points = fpe.run(datasets=datasets)

        display(self.flux_points.to_table(sed_type="dnde", formatted=True))
        display(self.flux_points.to_table(sed_type="e2dnde", formatted=True))
