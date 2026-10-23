"""Main analysis class."""

import hashlib
import json
import logging
import math
import uuid
from pathlib import Path

import astropy.units as u
import numpy as np
import yaml
from astropy.io import fits
from astropy.table import Table
from gammapy.datasets import Datasets, SpectrumDataset, SpectrumDatasetOnOff
from gammapy.estimators import FluxPointsEstimator, LightCurveEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
)

import v2dl5.plot as v2dl5_plot
import v2dl5.time as v2dl5_time


def _cache_inputs(provenance_inputs, datastore_path, observation_files):
    """Select observation, datastore index, and source inputs for cache keys."""
    index_files = {str(path.resolve()) for path in Path(datastore_path).glob("*-index.fits*")}
    return {
        key: value
        for key, value in provenance_inputs.items()
        if key in observation_files or key in index_files or key.endswith(".py")
    }


class Analysis:
    """
    Analysis class.

    Parameters
    ----------
    configuration : str
        Analysis configuration
    output_dir : str
        Output directory
    sky_regions: SkyRegions
        Sky regions
    data : Data
        Datastore object

    """

    def __init__(self, args_dict=None, sky_regions=None, v2dl5_data=None):
        self._logger = logging.getLogger(__name__)

        self.args_dict = args_dict
        self._output_dir = args_dict.get("output_dir", None)
        self.sky_regions = sky_regions
        self.v2dl5_data = v2dl5_data

        self.datasets = None
        self.spectral_model = None
        self.fit_results = None
        self.flux_points = None
        self.light_curves = {}
        self.integral_flux = None
        self.run_rows = []

    def run(self):
        """Run analysis."""
        self.stage = "reduction_failed"
        self._data_reduction()
        self.write()  # Preserve expensive reduction before inference or plots.
        data_sets = self.datasets
        if self.args_dict["datasets"]["stack"]:
            data_sets = Datasets([self.datasets.stack_reduce()])
        self._define_spectral_models(model=self.args_dict["fit"])
        self.stage = "fit_failed"
        try:
            self._spectral_fits(datasets=data_sets)
        except (ValueError, RuntimeError) as exc:
            self._logger.warning("Spectral fit failed: %s", exc)
        fit_success = self.fit_results is not None and self.fit_results.success
        self.write()
        self.stage = "flux_failed"
        self._integral_flux(data_sets)
        self.write()
        if fit_success:
            if self.args_dict["flux_points"].get("enabled", True):
                self._flux_points(data_sets)
            if self.args_dict["light_curve"].get("enabled", True):
                self.light_curves = self._analyse_light_curves(data_sets)
        self.write()
        if not fit_success:
            self.stage = "fit_failed"
            raise ValueError("Spectral fit failed; reduction and fixed-index flux retained")

    def _integral_flux(self, datasets):
        """Estimate a fixed-index PL integral flux, TS and one-sided 95% limit."""
        config = self.args_dict.get("integral_flux")
        if config is None:
            return
        low, high = (u.Quantity(config["energy"][key]) for key in ("min", "max"))
        axis = datasets[0].counts.geom.axes["energy"]
        edges = axis.edges
        low, high = axis.round(u.Quantity([low, high]))
        covered = np.zeros(len(edges) - 1, dtype=bool)
        for dataset in datasets:
            covered |= dataset.mask.data[:, 0, 0]
        bins = np.flatnonzero(covered)
        if not len(bins) or low < edges[bins[0]] or high > edges[bins[-1] + 1]:
            raise ValueError("Integral flux band extends outside usable data")
        band = (edges[:-1] < high) & (edges[1:] > low)
        if not np.all(covered[band]):
            raise ValueError("Integral flux band has unsupported energy bins")
        model = SkyModel(
            spectral_model=PowerLawSpectralModel(
                index=config.get("index", 2.5),
                amplitude="1e-12 cm-2 s-1 TeV-1",
                reference="1 TeV",
            ),
            name="reference_pl",
        )
        model.spectral_model.index.frozen = True
        model.spectral_model.amplitude.min = 0
        previous = datasets.models
        try:
            datasets.models = model
            estimator = FluxPointsEstimator(
                energy_edges=u.Quantity([low, high]),
                source="reference_pl",
                selection_optional=["ul", "errn-errp"],
                norm={"min": 0},
                n_sigma_ul=config.get("n_sigma_ul", 1.6448536269514722),
            )
            self.integral_flux = estimator.run(datasets)
            table = self.integral_flux.to_table(sed_type="flux")
            likelihood_table = self.integral_flux.to_table(sed_type="likelihood")
            success = np.asarray(table["success"], dtype=bool)
            finite_ul = np.isfinite(table["flux_ul"])
            failed_rows = np.flatnonzero(~success | ~finite_ul)
            if len(failed_rows):
                for row_index in failed_rows:
                    diagnostic = {
                        name: str(likelihood_table[name][row_index])
                        for name in likelihood_table.colnames
                    }
                    diagnostic["flux_ul"] = str(table["flux_ul"][row_index])
                    self._logger.error(
                        "Fixed-index flux/upper-limit result rejected for row %d "
                        "(success=%s, finite_flux_ul=%s): %s",
                        row_index,
                        success[row_index],
                        finite_ul[row_index],
                        diagnostic,
                    )
                reasons = []
                if not np.all(success):
                    reasons.append("fit unsuccessful")
                if not np.all(finite_ul):
                    reasons.append("upper limit is non-finite")
                raise ValueError(
                    "Fixed-index flux/upper-limit estimation failed: "
                    + ", ".join(reasons)
                    + "; see v2dl5.log for estimator values"
                )
        finally:
            datasets.models = previous

    def plot(self):
        """Plot all results."""
        _plot_dir = Path(f"{self._output_dir}/plots")

        plotter = v2dl5_plot.Plot(
            v2dl5_data=self.v2dl5_data,
            data_set=self.datasets,
            on_region=self.sky_regions.on_region,
            output_dir=_plot_dir,
        )
        plotter.plot_event_histograms()
        plotter.plot_irfs()
        plotter.plot_maps(exclusion_mask=self.sky_regions.exclusion_mask)
        plotter.plot_source_statistics()
        if (
            self.fit_results is not None
            and self.fit_results.success
            and self.flux_points is not None
        ):
            plotter.plot_spectra(
                flux_points=self.flux_points,
                model=self.spectral_model,
                y_min=self.args_dict["flux_points"].get("flux", {}).get("min"),
                y_max=self.args_dict["flux_points"].get("flux", {}).get("max"),
            )
        else:
            self._logger.warning("Skipping spectral plots because fit failed")
        plotter.plot_light_curves(self.light_curves)

    def write(self):
        """Write results to files."""
        data_dir = Path(self._output_dir) / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        self.datasets.write(data_dir / "datasets.yaml", overwrite=True)
        self.datasets.info_table(cumulative=False).write(
            data_dir / "runs.ecsv", format="ascii.ecsv", overwrite=True
        )
        Table(rows=self.run_rows, names=("obs_id", "status", "reason")).write(
            data_dir / "selection.ecsv", format="ascii.ecsv", overwrite=True
        )
        if self.integral_flux is not None:
            table = self.integral_flux.to_table(sed_type="flux", format="gadf-sed")
            table.meta.update(
                UL_TYPE="one-sided profile likelihood",
                ONE_SIDED_CONF=0.5
                * (
                    1
                    + math.erf(
                        self.args_dict["integral_flux"].get("n_sigma_ul", 1.6448536269514722)
                        / math.sqrt(2)
                    )
                ),
                REFERENCE_INDEX=self.args_dict["integral_flux"].get("index", 2.5),
            )
            table.write(data_dir / "integral_flux.ecsv", format="ascii.ecsv", overwrite=True)
        if self.fit_results is not None:
            self.fit_results.models.write(
                data_dir / "models.yaml", overwrite=True, write_covariance=True
            )
            self._write_yaml(
                {
                    "success": bool(self.fit_results.success),
                    "stat": float(self.fit_results.total_stat),
                },
                "fit.yaml",
            )
        self._write_datasets(self.flux_points, "flux_points.ecsv", "gadf-sed")
        for light_curve in self.light_curves.values():
            title_with_underscores = light_curve["title"].replace(" ", "_")
            self._write_datasets(
                light_curve["light_curve"],
                f"light_curve_{title_with_underscores}.ecsv",
                "lightcurve",
                "flux",
            )
        if self.spectral_model:
            self._write_yaml(self.spectral_model.to_dict(), "spectral_model.yaml")

    def _write_datasets(self, datasets, filename, file_format=None, sed_type=None):
        """
        Write datasets to disk.

        Parameters
        ----------
        datasets : Datasets
            Datasets
        filename : str
            Filename
        file_format: str
            Format specification (gammapy)

        """
        if datasets is None:
            return

        _data_dir = Path(self._output_dir) / "data"
        _data_dir.mkdir(parents=True, exist_ok=True)
        _out_file = _data_dir / filename
        self._logger.info(f"Writing datasets to {_out_file} ({file_format}, {sed_type})")
        if file_format is not None:
            if sed_type is not None:
                datasets.write(_out_file, overwrite=True, format=file_format, sed_type=sed_type)
            else:
                datasets.write(_out_file, overwrite=True, format=file_format)
        else:
            datasets.write(_out_file, overwrite=True)

    def _write_yaml(self, data_dict, filename):
        """
        Write model to disk.

        Parameters
        ----------
        data_dict : dict
            Data dictionary
        filename : str
            Filename

        """
        _data_dir = Path(f"{self._output_dir}/data")
        _data_dir.mkdir(parents=True, exist_ok=True)
        _out_file = f"{_data_dir}/{filename}"
        self._logger.info("Writing dataset to %s", _out_file)
        with open(_out_file, "w", encoding="utf-8") as outfile:
            yaml.dump(data_dict, outfile, default_flow_style=False)

    def _data_reduction(self):
        """Reduce data using the reflected region maker."""
        energy_axis = self._get_energy_axis(name="energy")
        energy_axis_true = self._get_energy_axis(name="energy_true")

        geom = RegionGeom.create(region=self.sky_regions.on_region, axes=[energy_axis])

        dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

        dataset_maker = SpectrumDatasetMaker(
            containment_correction=False, selection=["counts", "exposure", "edisp"]
        )
        bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=self.sky_regions.exclusion_mask)
        safe_mask_masker = SafeMaskMaker(
            methods=self.args_dict["datasets"]["safe_mask"]["methods"],
            **self.args_dict["datasets"]["safe_mask"]["parameters"],
        )
        self._logger.info(
            "Mask applied: %s, aeff_percent = %d",
            safe_mask_masker.methods,
            safe_mask_masker.aeff_percent,
        )

        self.datasets = Datasets()

        for observation in self.v2dl5_data.get_observations():
            obs_id = observation.obs_id
            cache_path = self._cache_path(obs_id)
            if cache_path is not None and cache_path.exists():
                with fits.open(cache_path) as hdus:
                    if any(hdu.verify_checksum() != 1 for hdu in hdus):
                        raise ValueError(f"Invalid cached dataset checksum: {cache_path}")
                cached = SpectrumDatasetOnOff.read(cache_path, format="gadf", checksum=True)
                self.datasets.append(cached)
                self.run_rows.append((int(obs_id), "accepted", "cached reduction"))
                continue
            dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
            dataset_on_off = bkg_maker.run(dataset, observation)
            if dataset_on_off is None or dataset_on_off.counts_off is None:
                self.run_rows.append((int(obs_id), "rejected", "No reflected OFF regions"))
                continue
            dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
            fit_range = self.args_dict["fit"]["fit_range"]
            dataset_on_off.mask_fit = dataset_on_off.counts.geom.energy_mask(
                energy_min=u.Quantity(fit_range["min"]), energy_max=u.Quantity(fit_range["max"])
            )
            valid = (
                np.isfinite(dataset_on_off.acceptance_off.data)
                & (dataset_on_off.acceptance_off.data > 0)
                & np.isfinite(dataset_on_off.acceptance.data)
                & (dataset_on_off.acceptance.data > 0)
            )
            exposure = dataset_on_off.exposure.data[:, 0, 0]
            response = dataset_on_off.edisp.get_edisp_kernel().pdf_matrix
            supported = np.sum(response * exposure[:, None], axis=0)
            valid &= (np.isfinite(supported) & (supported > 0))[:, None, None]
            dataset_on_off.mask_safe.data &= valid
            if not np.any(dataset_on_off.mask.data):
                self.run_rows.append((int(obs_id), "rejected", "No usable energy bins"))
                continue
            self.run_rows.append((int(obs_id), "accepted", ""))
            self.datasets.append(dataset_on_off)
            if cache_path is not None:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                temporary = cache_path.with_name(f"{uuid.uuid4()}.fits")
                dataset_on_off.write(temporary, format="gadf", checksum=True)
                temporary.replace(cache_path)

        if not self.datasets:
            Path(self._output_dir, "data").mkdir(parents=True, exist_ok=True)
            Table(rows=self.run_rows, names=("obs_id", "status", "reason")).write(
                Path(self._output_dir, "data", "selection.ecsv"),
                format="ascii.ecsv",
                overwrite=True,
            )
            self.stage = "no_usable_runs"
            raise ValueError("No usable runs after reduction")

        self._logger.info("Run-wise results:")
        self._print_results(self.datasets.info_table(cumulative=False))
        self._logger.info("Cumulative results:")
        self._print_results(self.datasets.info_table(cumulative=True))

    def _cache_path(self, obs_id):
        """Key optional reduction reuse by input contents and scientific settings."""
        cache = self.args_dict["datasets"].get("cache_dir")
        if not cache:
            return None
        provenance = json.loads(Path(self._output_dir, "provenance.json").read_text())
        store = self.v2dl5_data.get_data_store()
        base = Path(self.args_dict["observations"]["datastore"])
        files = {
            str((base / str(row["FILE_DIR"]) / str(row["FILE_NAME"])).resolve())
            for row in store.hdu_table
            if int(row["OBS_ID"]) == int(obs_id)
        }
        datasets = {
            key: value
            for key, value in self.args_dict["datasets"].items()
            if key not in ("cache_dir", "stack")
        }
        payload = {
            "inputs": _cache_inputs(provenance["inputs"], base, files),
            "versions": provenance["versions"],
            "datasets": datasets,
            "on_region": self.args_dict["on_region"],
            "fit_range": self.args_dict["fit"]["fit_range"],
            "bti": [item for item in (self.args_dict.get("bti") or []) if item["run"] == obs_id],
        }
        checksum = hashlib.sha256(yaml.safe_dump(payload, sort_keys=True).encode())
        checksum.update(self.sky_regions.exclusion_mask.data.tobytes())
        checksum.update(str(self.sky_regions.exclusion_mask.geom).encode())
        return Path(cache) / f"{obs_id}-{checksum.hexdigest()}.fits"

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
        print()

    def _spectral_fits(self, datasets=None):
        """Spectral fitting."""
        datasets.models = self.spectral_model

        _fit = Fit()
        self.fit_results = _fit.run(datasets=datasets)
        self._logger.info(self.fit_results)
        self.fit_results.models.to_parameters_table().pprint()

    def _define_spectral_models(self, model):
        """Spectral models."""
        _spectral_model = None
        if model.get("model", "pl") == "pl":
            _spectral_model = PowerLawSpectralModel(
                amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
                index=model.get("index", 2.0),
                reference=u.Quantity(model.get("reference_energy", "1.0 TeV")),
            )
        elif model.get("model") == "ecpl":
            _spectral_model = ExpCutoffPowerLawSpectralModel(
                amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
                index=model.get("index", 2.0),
                lambda_=u.Quantity(model.get("lambda", "0.1 TeV-1")),
                reference=u.Quantity(model.get("reference_energy", "1.0 TeV")),
            )
        else:
            raise ValueError(f"Unknown spectral model: {model.get('model')}")

        _spectral_model.amplitude.min = 0
        self.spectral_model = SkyModel(
            spectral_model=_spectral_model, name=model.get("model", "pl")
        )

    def _flux_points(self, datasets):
        """Calculate flux points."""
        energy_edges = (
            np.geomspace(
                u.Quantity(self.args_dict["flux_points"]["energy"]["min"]).to("TeV").value,
                u.Quantity(self.args_dict["flux_points"]["energy"]["max"]).to("TeV").value,
                self.args_dict["flux_points"]["energy"]["nbins"] + 1,
            )
            * u.TeV
        )

        self._logger.info("Estimating flux points")
        fpe = FluxPointsEstimator(energy_edges=energy_edges, selection_optional="all")
        self.flux_points = fpe.run(datasets=datasets)

        self.flux_points.to_table(sed_type="dnde", formatted=True).pprint()
        self.flux_points.to_table(sed_type="dnde", formatted=True).pprint_all()
        self.flux_points.to_table(sed_type="e2dnde", formatted=True).pprint()
        self.flux_points.to_table(sed_type="e2dnde", formatted=True).pprint_all()

    def _analyse_light_curves(self, data_sets):
        """
        Analyse light curves for different time interval.

        Parameters
        ----------
        data_sets : Datasets
            Datasets

        Returns
        -------
        light_curves : dict
            Light curves

        """
        light_curves = {}

        light_curves["per_observation"] = {
            "light_curve": None,
            "title": "per observation",
            "time_intervals": None,
        }
        light_curves["per_night"] = {
            "light_curve": None,
            "title": "per night",
            "time_intervals": (
                None
                if self.args_dict["datasets"]["stack"]
                else v2dl5_time.get_list_of_nights(
                    data_sets, time_zone=self.args_dict["light_curve"]["time_zone"]
                )
            ),
        }
        for time_bin_file in self.args_dict["light_curve"].get("time_bin_files", []):
            time_intervals = v2dl5_time.get_time_bins_from_file(time_bin_file)
            title = Path(time_bin_file).stem.replace("light_curve_", "")
            self._logger.info(f"Time intervals from {time_bin_file}: {time_intervals}")
            light_curves[title] = {
                "light_curve": None,
                "title": title,
                "time_intervals": time_intervals,
            }

        for light_curve in light_curves.values():
            self._logger.info(light_curve["title"])
            light_curve["light_curve"] = self._light_curve(data_sets, light_curve["time_intervals"])

        return light_curves

    def _light_curve(self, datasets, time_intervals=None):
        """
        Calculate light curve.

        Parameters
        ----------
        datasets : Datasets
            Datasets
        time_intervals : `~astropy.time.Time`
            Time intervals

        Returns
        -------
        _light_curve : `~gammapy.estimators.LightCurve`
            Light curve

        """
        self._logger.info(f"Estimating light curve {time_intervals}")

        lc_maker_1d = LightCurveEstimator(
            energy_edges=[
                u.Quantity(self.args_dict["light_curve"]["energy"]["min"]).to("TeV").value,
                u.Quantity(self.args_dict["light_curve"]["energy"]["max"]).to("TeV").value,
            ]
            * u.TeV,
            time_intervals=time_intervals,
            reoptimize=False,
            selection_optional=["all"],
        )
        _light_curve = lc_maker_1d.run(datasets)

        _table = _light_curve.to_table(sed_type="flux", format="lightcurve")
        print(
            _table[
                "time_min",
                "time_max",
                "e_min",
                "e_max",
                "flux",
                "flux_err",
                "flux_ul",
                "sqrt_ts",
            ]
        )
        _table.pprint_all()

        return _light_curve

    def _get_energy_axis(self, name="energy"):
        """Get energy axis."""
        _axes_dict = self.args_dict["datasets"]["geom"]["axes"][name]

        return MapAxis.from_energy_bounds(
            u.Quantity(_axes_dict["min"]).to("TeV").value,
            u.Quantity(_axes_dict["max"]).to("TeV").value,
            nbin=_axes_dict["nbins"],
            per_decade=_axes_dict.get("per_decade", True),
            unit="TeV",
            name=name,
        )


def read_datasets(filename):
    """Reload reduced GADF spectra preserving fit masks (Gammapy from_dict resets them)."""
    filename = Path(filename)
    content = yaml.safe_load(filename.read_text())
    return Datasets(
        [
            SpectrumDatasetOnOff.read(filename.parent / item["filename"], format="gadf")
            for item in content["datasets"]
        ]
    )
