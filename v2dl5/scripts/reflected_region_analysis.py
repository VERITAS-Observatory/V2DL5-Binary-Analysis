#!/usr/bin/python3
"""
Reflected region analysis using gammapy.

Example using a run list:

.. code-block:: console

    $ python reflected_region_analysis.py \
        --run_list region_run_list.txt \
        --output_dir 1D_reflected_region_output \
        --config examples/1D_reflected_region_config.yaml

"""

import argparse
import hashlib
import json
import logging
import platform
from importlib.metadata import version
from pathlib import Path

import yaml

import v2dl5.analysis
import v2dl5.configuration
import v2dl5.data
import v2dl5.sky_regions

logger = logging.getLogger(__name__)


def _observation_inputs(data, store):
    """Hash selected observation files, recording unavailable files for preflight."""
    selected_ids = set(data.runs)
    paths = {
        store / str(row["FILE_DIR"]) / str(row["FILE_NAME"])
        for row in data.get_data_store().hdu_table
        if int(row["OBS_ID"]) in selected_ids
    }
    inputs = {}
    for path in paths:
        try:
            checksum = hashlib.sha256(path.read_bytes()).hexdigest()
        except OSError:
            checksum = None
        inputs[str(path.resolve())] = checksum
    return inputs


def _check_resume(previous, fingerprint, output):
    """Require an identical completed attempt with intact saved products."""
    if (
        previous.get("fingerprint") != fingerprint
        or previous.get("status") not in ("complete", "partial")
    ):
        raise FileExistsError(f"Attempt cannot be reused: {output}; use a new output directory")
    products = previous.get("products")
    if not products:
        raise FileExistsError(f"Missing product manifest: {output}; use a new output directory")
    for name, checksum in products.items():
        path = output / name
        try:
            valid = (
                path.resolve().is_relative_to(output.resolve())
                and hashlib.sha256(path.read_bytes()).hexdigest() == checksum
            )
        except OSError:
            valid = False
        if not valid:
            raise FileExistsError(f"Invalid saved product: {path}; use a new output directory")


def _fingerprint(inputs, config, plots):
    """Create a stable fingerprint for the analysis inputs and runtime."""
    versions = {
        name: version(name)
        for name in (
            "gammapy",
            "astropy",
            "numpy",
            "regions",
            "scipy",
            "iminuit",
            "PyYAML",
        )
    }
    versions["python"] = platform.python_version()
    return hashlib.sha256(
        json.dumps(
            {
                "inputs": inputs,
                "config": config,
                "plots": plots,
                "versions": versions,
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()


def _parse():
    """
    Parse command line arguments.

    Returns
    -------
    dict
        Command line parameters.

    """
    parser = argparse.ArgumentParser(description="Reflected region analysis.")

    parser.add_argument(
        "--run_list",
        type=str,
        required=False,
        help="List of runs to be processed.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to configuration file.",
    )
    parser.add_argument("--resume", action="store_true", help="Reuse an exact completed result")
    parser.add_argument("--no-plots", action="store_true", help="Save numerical products only")
    return parser.parse_args()


def main():
    """Reflected region analysis."""
    logging.root.setLevel(logging.INFO)

    args = _parse()
    output = Path(args.output_dir)
    # Configuration errors also get a machine-readable disposition.
    output.mkdir(parents=True, exist_ok=True)
    status_file = output / "status.json"
    if status_file.exists() and not args.resume:
        raise FileExistsError(f"Existing attempt: {output}; use a new output directory")
    previous = json.loads(status_file.read_text()) if status_file.exists() else None
    analysis = None
    stage = "invalid_config"
    status = {"status": "running", "stage": stage}

    def save_status():
        temporary = output / "status.json.tmp"
        temporary.write_text(json.dumps(status, indent=2, allow_nan=False) + "\n")
        temporary.replace(status_file)

    lock = output / ".running"
    # mkdir is atomic, so only one process can own this output directory.
    # Remove this directory manually before resuming after a crashed attempt.
    lock.mkdir()
    try:
        args_dict = v2dl5.configuration.configuration(args=args)
        sky_regions = v2dl5.sky_regions.SkyRegions(args_dict=args_dict)
        args_dict["on_region"].update(
            frame="icrs",
            lon=f"{sky_regions.target.icrs.ra.deg} deg",
            lat=f"{sky_regions.target.icrs.dec.deg} deg",
        )
        stage = "missing_data"
        data = v2dl5.data.Data(args_dict=args_dict, target=sky_regions.target)
        inputs = {}
        paths = [Path(args.config), Path(__file__)]
        if args.run_list:
            paths.append(Path(args.run_list))
        store = Path(args_dict["observations"]["datastore"])

        stage = "missing_data"
        inputs.update(_observation_inputs(data, store))
        paths.extend(store.glob("*-index.fits*"))
        paths.extend(Path(v2dl5.analysis.__file__).parent.rglob("*.py"))
        from importlib.resources import files

        catalogue = args_dict["datasets"]["exclusion_region"].get("star_file")
        if catalogue:
            path = Path(catalogue)
            paths.append(
                path if path.exists() else Path(str(files("v2dl5").joinpath("data", path.name)))
            )
        regions = args_dict["datasets"]["exclusion_region"].get("regions_file")
        if regions:
            paths.append(Path(regions))
        paths.extend(Path(path) for path in args_dict["light_curve"].get("time_bin_files", []))
        for path in set(paths):
            inputs[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
        fingerprint = _fingerprint(inputs, args_dict, plots=not args.no_plots)
        if previous is not None:
            _check_resume(previous, fingerprint, output)
            logger.info("Reusing completed result with identical inputs: %s", output)
            return
        # Unusable runs are skipped (and listed in availability.ecsv).
        data.preflight(output / "availability.ecsv")
        status["fingerprint"] = fingerprint
        save_status()
        (output / "config.resolved.yaml").write_text(yaml.safe_dump(args_dict))
        (output / "provenance.json").write_text(
            json.dumps(
                {
                    "inputs": inputs,
                    "python": platform.python_version(),
                    "versions": {
                        name: version(name)
                        for name in (
                            "v2dl5",
                            "gammapy",
                            "astropy",
                            "numpy",
                            "regions",
                            "scipy",
                            "iminuit",
                            "PyYAML",
                        )
                    },
                },
                indent=2,
            )
            + "\n"
        )
        stage = "invalid_config"
        sky_regions.update_regions(
            args_dict=args_dict,
            on_region_radius=data.get_on_region_radius(),
            max_wobble_distance=data.get_max_wobble_distance(
                fov=args_dict["datasets"]["exclusion_region"]["fov"]
            ),
        )
        (output / "config.resolved.yaml").write_text(yaml.safe_dump(args_dict))
        sky_regions.exclusion_mask.write(output / "exclusion_mask.fits", overwrite=True)
        analysis = v2dl5.analysis.Analysis(
            args_dict=args_dict,
            sky_regions=sky_regions,
            v2dl5_data=data,
        )
        stage = "analysis_failed"
        analysis.run()
        analysis.write()
        stage = "plot_failed"
        if not args.no_plots:
            analysis.plot()
        rejected = [row for row in analysis.run_rows if row[1] != "accepted"]
        products = {
            str(path.relative_to(output)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in output.rglob("*")
            if path.is_file() and path.name not in ("status.json", "status.json.tmp", "v2dl5.log")
        }
        status.update(
            products=products,
            stage="finished",
            status="partial" if rejected else "complete",
            requested=len(data.runs),
            accepted=len(analysis.datasets),
            rejected=len(rejected),
        )
        if rejected:
            status["reason"] = "Some requested runs were rejected; see data/selection.ecsv"
        save_status()
    except FileExistsError:
        raise  # Preserve the original attempt, including its completion marker.
    except Exception as exc:
        if status.get("status") != "partial":
            status.update(
                status=getattr(analysis, "stage", stage) if stage == "analysis_failed" else stage,
                stage=stage,
                reason=str(exc),
            )
        save_status()
        raise
    finally:
        lock.rmdir()


if __name__ == "__main__":
    main()
