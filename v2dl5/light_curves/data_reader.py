import csv
import logging

import astropy.units as u
import yaml

import v2dl5.light_curves.orbital_period as orbit


class LightCurveDataReader:
    """
    Data reader for light-curve analysis.

    Allows to read gamma-ray, optical or x-ray data. Enrich data
    with information required for orbital analysis.

    Parameters
    ----------
    configuration_file: str
        Configuration file with data description.
    binary: dict
        Binary parameters.

    """

    def __init__(self, configuration_file, binary=None):

        self._logger = logging.getLogger(__name__)

        self.data_dict = {}
        try:
            with open(configuration_file, "r") as file:
                _yml_in = yaml.safe_load(file)
                self.config = _yml_in["data"]
        except FileNotFoundError as err:
            self._logger.error("Configuration file not found: %s", configuration_file)
            raise err
        except KeyError as err:
            self._logger.error("Data key not found in configuration file %s", configuration_file)
            raise err

        self.binary = binary

    def read_data(self):
        """
        Read data from files as described in configuration file.

        Returns
        -------
        dict
            Data dictionary with data for each instrument.

        """

        for data_config in self.config:
            self.data_dict[data_config["instrument"]] = self._read_fluxes_from_file(data_config)
            self._add_orbital_parameters(self.data_dict[data_config["instrument"]])

    def _read_fluxes_from_file(self, data_config):
        """
        Read flux from file.

        """

        try:
            if data_config["file_name"].endswith(".csv") or data_config["file_name"].endswith(
                ".ecsv"
            ):
                return self._read_fluxes_from_csv_file(data_config)
        except KeyError:
            self._logger.error("File name not found in configuration file")
            raise KeyError

    def _read_fluxes_from_csv_file(self, data_config, TimeMinMax=True, MJD_min=-1.0, MJD_max=-1.0):
        """
        Read gamma-ray fluxes from csv file (open gamma-ray format)
        - ignore all lines with '#'
        - accept a random number of spaces as delimiter

        Parameters:
        -----------
        data_config: dict
            configuration dictionary
        TimeMinMax: bool
            flag to read time_min and time_max
        MJD_min: float
            MJD min value for MJD cut
        MJD_max: float
            MJD max value for MJD cut

        """
        fp = open(data_config["file_name"])
        rdr = csv.DictReader(
            filter(lambda row: row[0] != "#", fp), delimiter=" ", skipinitialspace=True
        )
        f = {}
        f["time_min"] = []
        f["time_max"] = []
        f["flux"] = []
        f["flux_err"] = []
        for row in rdr:
            if TimeMinMax:
                t_min = float(row["time_min"])
                t_max = float(row["time_max"])
            else:
                t_min = float(row["time"])
                t_max = float(row["time"]) + 0.1
            if MJD_min > 0 and t_min < MJD_min:
                continue
            if MJD_max > 0 and t_max > MJD_max:
                continue

            f["time_min"].append(t_min)
            f["time_max"].append(t_max)

            f["flux"].append(float(row["flux"]))
            if "flux_err" in row:
                f["flux_err"].append(float(row["flux_err"]))
            elif "flux_up" in row:
                f["flux_err"].append(0.5 * abs(float(row["flux_up"]) - float(row["flux_down"])))
        fp.close()

        f["MJD"] = [0.5 * (a + b) for a, b in zip(f["time_min"], f["time_max"])]
        f["MJD_err"] = [0.5 * (b - a) for a, b in zip(f["time_min"], f["time_max"])]

        return f

    def _add_orbital_parameters(self, data):
        """
        Add orbital phase to light-curve data data.

        Parameters
        ----------
        data: dict
            Data dictionary with light-curve data.

        """

        orbital_period = self.binary["orbital_period"]
        mjd_0 = self.binary["mjd_0"]

        data["phase"] = [
            orbit.get_orbital_phase(mjd, orbital_period, mjd_0, True) for mjd in data["MJD"]
        ]
        data["phaseN"] = [
            orbit.get_orbital_phase(mjd, orbital_period, mjd_0, False) for mjd in data["MJD"]
        ]
        data["phase_err_low"] = [
            orbit.get_orbital_phase_range(a, b, c, False, orbital_period, mjd_0)
            for a, b, c in zip(data["time_min"], data["time_max"], data["phase"])
        ]
        data["phase_err_hig"] = [
            orbit.get_orbital_phase_range(a, b, c, True, orbital_period, mjd_0)
            for a, b, c in zip(data["time_min"], data["time_max"], data["phase"])
        ]
        data["phase_err"] = [data["phase_err_low"], data["phase_err_hig"]]

    def convert_photon_to_energy_flux(C_v, C_e, E_0, gamma):
        """
        Convert photon to energy flux

        """

        f = (-1.0 * gamma + 1) / (-1.0 * gamma + 2)
        # conversion to erg
        f = f * (E_0.to(u.erg)).value
        return [v * f for v in C_v], [e * f for e in C_e]
