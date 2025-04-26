"""Binary light curve plotting."""

import logging
import math

import matplotlib.pyplot as plt
import numpy as np

import v2dl5.plotting.utilities as plotting_utilities


class BinaryLightCurvePlotter:
    """
    Binary light curve plotter.

    Parameters
    ----------
    data : dict
        Data dictionary with light-curve data.
    config : dict
        Configuration dictionary.
    binary : dict
        Binary parameters.

    """

    def __init__(self, data, config, binary):
        self._logger = logging.getLogger(__name__)

        self.data = data
        self.config = config
        self.binary = binary

    def plot_flux_vs_time(
        self,
        time_axis="MJD",
        y_axis="flux",
        mjd_min=None,
        mjd_max=None,
        orbit_number=None,
        phase_min=None,
        phase_max=None,
        file_type=".pdf",
        figure_dir="./figures/",
        axes=None,
        fontsize=None,
        markersize=None,
    ):
        """
        Plot data (e.g., 'flux') vs time (MJD or orbital phase).

        Parameters
        ----------
        time_axis: str
            Time axis (MJD or orbital phase).
        y_axis: str
            Y-axis data (e.g. 'flux').
        mjd_min: float
            Minimum MJD value.
        mjd_max: float
            Maximum MJD value.
        phase_min: float
            Minimum phase value.
        phase_max: float
            Maximum phase value.
        orbit_number: int
            Orbit number to plot.
        file_type: str
            File type for plots (e.g., '.pdf', '.png').
        figure_dir: str
            Directory for saving figures.
        axes: matplotlib.axes.Axes
            Axes object to plot on.
        fontsize: int
            Font size for labels.
        markersize: int
            Marker size for points.
        """
        self._logger.info(
            f"Plotting {y_axis} vs {time_axis}(MJD {mjd_min}, {mjd_max}, orbit id {orbit_number})"
        )

        ax = axes if axes else plotting_utilities.paper_figures(None, None)

        for idx, (instrument, data) in enumerate(self.data.items()):
            if not self.plot_this_instrument(self.config[idx], y_axis):
                continue
            color, marker = self.get_marker_and_color(idx)

            x, y, e, x_ul, y_ul = self._get_light_curve_in_mjd_limits(
                data,
                y_axis,
                time_axis,
                mjd_min,
                mjd_max,
                orbit_number,
                phase_min,
                phase_max,
                self.config[idx].get("significance_min", None),
            )
            plt.errorbar(
                x,
                y,
                e,
                None,
                label=(
                    instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                marker=marker,
                linestyle="none",
                fillstyle="full",
                linewidth=plotting_utilities.get_line_width(),
                markersize=(
                    plotting_utilities.get_marker_size() if markersize is None else markersize
                ),
            )
            if len(y_ul) > 0:
                plt.errorbar(
                    x_ul,
                    y_ul,
                    yerr=0.1 * max(y_ul),
                    color=color,
                    fmt="_",
                    linestyle="none",
                    fillstyle="none",
                    uplims=True,
                    linewidth=plotting_utilities.get_line_width(),
                    markersize=plotting_utilities.get_marker_size() * 0.5,
                )

        ax.axhline(0, color="lightgray", linestyle="--")
        plt.xlabel(self._get_time_axis_label(time_axis), fontsize=fontsize)
        if mjd_min is not None and mjd_max is not None:
            ax.set_xlim([mjd_min, mjd_max])
        plt.ylabel(self.config[0].get(y_axis + "_axis_label", ""), fontsize=fontsize)
        if axes is None:
            plt.legend()
            plotting_utilities.print_figure(
                f"Light-Curve-{self.binary['name']}-{y_axis}-vs-{time_axis.replace(' ', '-')}",
                file_type=file_type,
                figure_dir=figure_dir,
            )

    def _get_number_columns_and_rows(self, number):
        """Get number of columns and rows for plotting."""
        n_columns = 4
        fig_size = (10, 4)
        if number == 9:
            n_columns = 3
            fig_size = (10, 8)
        if number > 9:
            n_columns = 6
            fig_size = (12, 6)
        if number == 16:
            n_columns = 4
            fig_size = (10, 8)
        if number > 16:
            n_columns = 8
            fig_size = (16, 10)
        return n_columns, math.ceil(number / n_columns), fig_size

    def plot_flux_vs_phase_for_individual_orbits(
        self, instrument, y_axis="flux", file_type=".pdf", figure_dir="./figures/"
    ):
        """Plot flux vs phase with one plot per orbit."""
        self._logger.info(f"Plotting {y_axis} vs phase for individual orbits")
        time_axis = "orbital phase"

        data = self.data[instrument]
        if y_axis not in data:
            self._logger.warning(f"Y-axis {y_axis} not found in data for {instrument}")
            return
        orbits = sorted(set(data["orbit_number"]))
        self._logger.info(f"Orbits for {instrument} (total number {len(orbits)}): {orbits}")

        y_min, y_max = self._global_value_extrema(data[y_axis], data.get(y_axis + "_err", 0.0))
        fontsize = 6
        n_columns, n_rows, figsize = self._get_number_columns_and_rows(len(orbits))
        self._logger.info(f"Number of columns: {n_columns}, number of rows: {n_rows}")
        plt.figure(figsize=figsize)

        for i, orbit_id in enumerate(orbits):
            axes = plt.subplot(n_rows, n_columns, i + 1)
            axes.set_xlim([0.0, 1.0])
            axes.set_ylim([y_min, y_max])
            self.plot_flux_vs_time(
                time_axis=time_axis,
                y_axis=y_axis,
                mjd_min=None,
                mjd_max=None,
                orbit_number=orbit_id,
                axes=axes,
                fontsize=fontsize,
                markersize=plotting_utilities.get_marker_size() * 0.5,
            )
            plt.rc("xtick", labelsize=fontsize)
            plt.rc("ytick", labelsize=fontsize)
            plt.text(
                0.05,
                0.95,
                f"Orbit {orbit_id}",
                transform=plt.gca().transAxes,
                fontsize=fontsize,
                verticalalignment="top",
            )

        plt.tight_layout()
        plotting_utilities.print_figure(
            f"Light-Curve-Orbits-{self.binary['name']}-{y_axis}-vs-{time_axis.replace(' ', '-')}",
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def plot_flux_vs_orbit_number(
        self, instrument, y_axis="flux", phase_bins=10, file_type=".pdf", figure_dir="./figures/"
    ):
        """Plot flux vs orbit number."""
        self._logger.info(f"Plotting {y_axis} vs orbit number")
        time_axis = "orbit number"

        data = self.data[instrument]
        if y_axis not in data:
            self._logger.warning(f"Y-axis {y_axis} not found in data for {instrument}")
            return
        orbits = sorted(set(data["orbit_number"]))
        self._logger.info(f"Orbits for {instrument} (total number {len(orbits)}): {orbits}")

        y_min, y_max = self._global_value_extrema(data[y_axis], data.get(y_axis + "_err", 0.0))
        fontsize = 6
        n_columns, n_rows, figsize = self._get_number_columns_and_rows(phase_bins)
        self._logger.info(f"Number of columns: {n_columns}, number of rows: {n_rows}")
        plt.figure(figsize=figsize)

        for i in range(phase_bins):
            axes = plt.subplot(n_rows, n_columns, i + 1)
            axes.set_xlim([min(orbits), max(orbits)])
            axes.set_ylim([y_min, y_max])

            phase_min = i * (1.0 / phase_bins)
            phase_max = (i + 1) * (1.0 / phase_bins)

            self.plot_flux_vs_time(
                time_axis=time_axis,
                y_axis=y_axis,
                mjd_min=None,
                mjd_max=None,
                orbit_number=None,
                phase_min=phase_min,
                phase_max=phase_max,
                axes=axes,
                fontsize=fontsize,
                markersize=plotting_utilities.get_marker_size() * 0.5,
            )
            plt.rc("xtick", labelsize=fontsize)
            plt.rc("ytick", labelsize=fontsize)
            plt.text(
                0.05,
                0.95,
                f"Orbit phase {phase_min:.2f} - {phase_max:.2f}",
                transform=plt.gca().transAxes,
                fontsize=fontsize,
                verticalalignment="top",
            )

        plt.tight_layout()
        plotting_utilities.print_figure(
            (
                f"Light-Curve-Orbit-Number-{self.binary['name']}-"
                f"{y_axis}-vs-{time_axis.replace(' ', '-')}"
            ),
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def _global_value_extrema(self, y, y_err=0.0):
        """Absolute min and max of y values."""
        y = np.array(y)
        y_err = np.array(y_err)

        y = y[~np.isnan(y) & np.isfinite(y)]
        y_err = y_err[~np.isnan(y_err) & np.isfinite(y_err)]

        if len(y) == 0:
            return -1, 1

        return (
            min([min(y) - max(y_err), 0.0]),
            max(y + y_err) * 1.2,
        )

    def _get_time_axis_label(self, time_axis):
        """Return time axis label."""
        if time_axis == "orbital phase":
            return "phase"
        if time_axis == "orbit number":
            return "orbit number"
        return "Modified Julian Date (MJD)"

    def _get_light_curve_in_mjd_limits(
        self,
        data,
        y_key,
        time_axis,
        mjd_min=None,
        mjd_max=None,
        orbit_number=None,
        phase_min=None,
        phase_max=None,
        min_significance=None,
    ):
        """
        Get light curve restricted in MJD or for a certain orbit number.

        Parameters
        ----------
        data : dict
            Light-curve data.
        y_key : str
            Key for y-axis data (e.g. 'flux').
        time_axis : str
            Time axis (MJD or orbital phase)
        mjd_min : float
            Minimum MJD value.
        mjd_max : float
            Maximum MJD value.
        orbit_number : int
            Select orbit number.
        phase_min : float
            Minimum phase value.
        phase_max : float
            Maximum phase value.
        min_significance : float
            Minimum significance for values.

        Returns
        -------
        list
            Light-curve data as lists.
        """
        x, y, e, x_ul, y_ul = [], [], [], [], []
        if y_key not in data:
            self._logger.warning(f"Y-axis {y_key} not found in data")
            return x, y, e, x_ul, y_ul
        mjd = data["MJD"]
        for i, t in enumerate(mjd):
            if (
                min_significance is not None
                and data.get("significance", [None] * len(mjd))[i] < min_significance
            ):
                continue
            if (mjd_min is not None and t < mjd_min) or (mjd_max is not None and t > mjd_max):
                continue
            if orbit_number is not None and data["orbit_number"][i] != orbit_number:
                continue
            if (phase_min is not None and data["phase"][i] < phase_min) or (
                phase_max is not None and data["phase"][i] > phase_max
            ):
                continue
            if time_axis == "MJD":
                w_x = data["MJD"][i]
            elif time_axis == "orbital phase":
                w_x = data["phase"][i]
            elif time_axis == "orbit number":
                w_x = data["orbit_number"][i]
            else:
                raise ValueError(f"Unknown time axis: {time_axis}")
            y_key_ul = y_key + "_ul"
            y_key_err = y_key + "_err"
            _ul = data.get(y_key_ul, [None] * len(mjd))[i]
            if _ul is not None and _ul > 0:
                x_ul.append(w_x)
                y_ul.append(data[y_key_ul][i])
            elif y_key_ul not in data or data[y_key_ul][i] < 0:
                x.append(w_x)
                y.append(data[y_key][i])
                if y_key_err not in data:
                    e.append(0.0)
                else:
                    e.append(data[y_key_err][i])
        return x, y, e, x_ul, y_ul

    def plot_this_instrument(self, config, y_axis):
        """Return if this instrument/axis should be plotted."""
        plot_this = config.get("plot_instrument", True)
        if y_axis not in config.get("plot_axis", []):
            plot_this = False
        return plot_this

    def get_marker_and_color(self, idx):
        """Return marker and color."""
        colors = plotting_utilities.get_color_list(len(self.data))
        color = (
            colors[idx]
            if self.config[idx].get("marker_color") is None
            else self.config[idx]["marker_color"]
        )
        marker = (
            plotting_utilities.get_marker_list()[idx]
            if self.config[idx].get("marker_type") is None
            else self.config[idx]["marker_type"]
        )
        return color, marker

    def plot_live_time_vs_phase_bin(
        self,
        instrument,
        phase_bins=10,
        file_type=".pdf",
        figure_dir="./figures/",
    ):
        """Plot live time histogram vs phase bin."""
        self._logger.info("Plotting live time histogram vs phase bin")

        ax = plotting_utilities.paper_figures(None, None)
        ax.set_xlim([0, 1])

        for idx, (instrument, data) in enumerate(self.data.items()):
            if self.config[idx].get("plot_live_time_histogram", False) is False:
                continue
            color, _ = self.get_marker_and_color(idx)

            x, y, _, _, _ = self._get_light_curve_in_mjd_limits(
                data, "live_time", "orbital phase",
            )
            if len(x) == 0:
                continue
            x = np.array(x)
            y = np.array(y)
            bin_edges = np.linspace(0, 1, phase_bins + 1)
            bin_width = bin_edges[1] - bin_edges[0]
            bin_heights, _ = np.histogram(x, bins=bin_edges, weights=y)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            plt.bar(
                bin_centers,
                bin_heights,
                width=bin_width*0.9,
                label=(
                    instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                alpha=0.7,
            )

        plt.xlabel("Phase")
        plt.ylabel("Live Time (h)")

        plt.legend()
        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-Live-Time-vs-Phase-Bin",
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def plot_index_vs_flux(
        self,
        instrument,
        file_type=".pdf",
        figure_dir="./figures/",
    ):
        """Plot spectral index vs flux."""
        self._logger.info("Plotting spectral index vs flux")

        _ = plotting_utilities.paper_figures(None, None)

        for idx, (instrument, data) in enumerate(self.data.items()):
            if self.config[idx].get("plot_flux_vs_index", False) is False:
                continue
            color, _ = self.get_marker_and_color(idx)

            _, x, x_e, _, _ = self._get_light_curve_in_mjd_limits(
                data, "flux", "orbital phase",
            )
            _, y, y_e, _, _ = self._get_light_curve_in_mjd_limits(
                data, "index", "orbital phase",
            )
            if len(x) == 0 or len(y) == 0:
                continue
            x = np.array(x)
            y = np.array(y)
            plt.errorbar(
                x,
                y,
                xerr=x_e,
                yerr=y_e,
                label=(
                    instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                linestyle="none",
                marker='o',
                alpha=0.7,
            )

        plt.xlabel(self.config[0].get("flux_axis_label", ""))
        plt.ylabel(self.config[0].get("index_axis_label", ""))

        plt.legend()
        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-Flux-vs-Index",
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def plot_distribution(
        self,
        instrument,
        y_axis="flux",
        file_type=".pdf",
        figure_dir="./figures/",
    ):
        """Plot distribution of y-axis data."""
        self._logger.info(f"Plotting {y_axis} distribution")

        _ = plotting_utilities.paper_figures(None, None)

        for idx, (instrument, data) in enumerate(self.data.items()):
            if self.config[idx].get("plot_1d_distribution", False) is False:
                continue
            color, _ = self.get_marker_and_color(idx)

            _, y, _, _, _ = self._get_light_curve_in_mjd_limits(
                data, y_axis, "orbital phase",
            )
            if len(y) == 0:
                continue
            y = np.array(y)
            plt.hist(
                y,
                bins=25,
                label=(
                    instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                alpha=0.7,
            )

        plt.xlabel(self.config[0].get(y_axis + "_axis_label", ""))
        plt.ylabel("Counts")

        plt.legend()
        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-{y_axis}-Distribution",
            file_type=file_type,
            figure_dir=figure_dir,
        )





# Temporary stuff - probably not needed
#
# def plotLightCurve_fluxvsPhase_inOrbits(
#         fDataDict,
#         PlotInstruments,
#         F,
#         lc_spline_bin_centers,
#         lc_spline_sv,
#         lc_spline_sv_err,
#         orbital_period_BL=315.,
#         orbital_periodBins=20,
#         plot_variable=None):
#     """plot flux vs orbital phase (separated in orbits, all in one plot)
#     """
#
#     MJD_firstOrbit, N_orbits = getNumberOfOrbits(
#         fDataDict, PlotInstruments,
#         orbital_period_BL, False)
#
#     lightCurvePlottingUtilities.paper_figures(4, 4)
#     colors = lightCurvePlottingUtilities.getColorList(N_orbits)
#     if N_orbits < 6:
#         markers = lightCurvePlottingUtilities.getMarkerList()
#     else:
#         markers = ['o']*N_orbits
#
#     if len(PlotInstruments) < 1:
#         return
#
#     # plot average and interpolated light curves
#     glabel = PlotInstruments[0] + " (average)"
#     if lc_spline_bin_centers and len(lc_spline_bin_centers) > 0:
#         plt.plot(
#             lc_spline_bin_centers,
#             lc_spline_sv,
#             color='tab:gray',
#             linestyle='--',
#             linewidth=lightCurvePlottingUtilities.getLineWidth())
#         plt.fill_between(lc_spline_bin_centers,
#                         lc_spline_sv - lc_spline_sv_err,
#                         lc_spline_sv + lc_spline_sv_err,
#                         color='tab:gray', linestyle='--',
#                         linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                         alpha=0.3)
#
#     # plot one light curve per orbital period
#     for i in range(0, N_orbits):
#         x = []
#         y = []
#         ex = []
#         mj = []
#         orbit_min = MJD_firstOrbit + i * orbital_period_BL
#         orbit_max = MJD_firstOrbit + (i + 1) * orbital_period_BL - 1.
#         for j in range(len(PlotInstruments)):
#             i_plotValue, i_plotError = lightCurvePlottingUtilities.get_plotting_variable(
#                 plot_variable, j)
#
#             for p in range(len(fDataDict[PlotInstruments[j]]['phaseN'])):
#                 if fDataDict[PlotInstruments[j]]['MJD'][p] >= orbit_min \
#                         and fDataDict[PlotInstruments[j]]['MJD'][p] < orbit_max + 1.:
#
#                     x.append(fDataDict[PlotInstruments[j]]['phase'][p])
#                     y.append(fDataDict[PlotInstruments[j]][i_plotValue][p])
#                     ex.append(fDataDict[PlotInstruments[j]][i_plotError][p])
#                     mj.append(fDataDict[PlotInstruments[j]]['MJD'][p])
#
#         if len(mj):
#             OrbitPhStr = "MJD %d - %d" % (
#                 orbit_min, orbit_max)
#             plt.errorbar(
#                 x,
#                 y,
#                 ex,
#                 None,
#                 color=colors[i],
#                 marker=markers[i],
#                 linestyle='none',
#                 label=OrbitPhStr,
#                 linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                 markersize=lightCurvePlottingUtilities.getMarkerSize())
#
#     plt.xlabel(
#         lightCurvePlottingUtilities.get_orbital_phase_axis_string(orbital_period_BL) )
#     plt.ylabel(
#         lightCurvePlottingUtilities.getFluxAxisString(
#         PlotInstruments[0],plot_variable) )
#     if PlotInstruments[0].find('Optical') < 0:
#         plt.axhline(y=0, linestyle=':')
#     if getPrintInstrumentName(PlotInstruments).find( "fwhm" ) < 0 and \
# getPrintInstrumentName(PlotInstruments).find( "ew" ) < 0:
#         plt.legend(prop={'size': 10}, framealpha=0.1)
#
#     lightCurvePlottingUtilities.printFigure(
#         getPrintInstrumentName(PlotInstruments) +
#         "-HESSJ0632p057-LC-phaseFolded-%dd-Orbits" %
#         orbital_period_BL)
#
#
#
# def plotAverageLightCurve_fluxvsPhase(
#         fDataDict,
#         PlotInstruments,
#         orbital_period_BL=315.,
#         orbital_periodBins=20):
#     """
#     plot average flux vs orbital phase
#
#     calculates also average light curves
#     """
#
#     print("Plot phase binned averaged light curve:")
#     print("\t orbital phase (%.1f d)" % orbital_period_BL)
#     print("\t number of phase bins (%d)" % orbital_periodBins)
#     print("\t calculating average light curve for ", PlotInstruments)
#     if len(PlotInstruments) < 1:
#         return
#
#     # copy all data into one set of arrays
#     # averaged light curve is calculated from
#     # light curve bins of all data
#     i_MJD = []
#     i_flux = []
#     i_flux_err = []
#     for I in PlotInstruments:
#         i_MJD.extend(fDataDict[I]['MJD'])
#         i_flux.extend(fDataDict[I]['flux'])
#         i_flux_err.extend(fDataDict[I]['flux_err'])
#
#     lightCurvePlottingUtilities.paper_figures(4, 4)
#
#     # calculate phase binned average light curve
#     lc_bincenters, lc_mean, lc_std = lightCurveAverageing.calculateAverageLightCurve(
#         i_MJD, i_flux, i_flux_err, orbital_period_BL, orbital_periodBins)
#
#     # bin width
#     lc_binw = np.repeat(0.5 / orbital_periodBins, len(lc_bincenters))
#     plt.errorbar(lc_bincenters, lc_mean, lc_std, lc_binw,
#                  color='b', linestyle='none',
#                  linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                  markersize=lightCurvePlottingUtilities.getMarkerSize(),
#                  marker='o', label='average')
#
#     # calculate cubic spline smoothed average light curve
#     F, lc_spline_bin_centers, lc_spline_sv = \
#         lightCurveAverageing.smoothCubeSplineLightCurve(
#             i_MJD, i_flux, i_flux_err,
#             orbital_period_BL, orbital_periodBins,
#             True, 0.)
#     # get errors on the same
#     F_temp, lc_spline_bin_centers, lc_spline_sv_err = \
#         lightCurveAverageing.smoothCubeSplineLightCurve(
#             i_MJD, i_flux, i_flux_err,
#             orbital_period_BL, orbital_periodBins,
#             True, 1.)
#     lc_spline_sv_err = list(
#         np.array(lc_spline_sv_err) -
#         np.array(lc_spline_sv))
#
#     plt.plot(lc_spline_bin_centers, lc_spline_sv, color='r', linestyle='--',
#              linewidth=lightCurvePlottingUtilities.getLineWidth(),
#              label='spline average')
#
#     # data set used for averaging
#     wP = []
#     wP = [lightCurveAnalysisorbital_period.getOrbitalPhase(
#         x, orbital_period_BL) for x in i_MJD]
#     plt.errorbar(wP, i_flux, i_flux_err, None,
#                  color='g', linestyle='none',
#                  linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                  markersize=lightCurvePlottingUtilities.getMarkerSize(),
#                  marker='o', label='average')
#
#     plt.fill_between(lc_spline_bin_centers,
#                      lc_spline_sv - lc_spline_sv_err,
#                      lc_spline_sv + lc_spline_sv_err,
#                      color='r', linestyle='--',
#                      linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                      alpha=0.3)
#
#     plt.xlabel(
#         lightCurvePlottingUtilities.get_orbital_phase_axis_string(orbital_period_BL))
#     plt.ylabel(lightCurvePlottingUtilities.getFluxAxisString(PlotInstruments[0]))
#     plt.legend()
#     plt.axhline(y=lc_spline_sv[0], linestyle=':')
#
#     lightCurvePlottingUtilities.printFigure(
#         getPrintInstrumentName(PlotInstruments) +
#         "-HESSJ0632p057-LC-phaseFolded-%dd-Average" %
#         orbital_period_BL)
#
#     return F, lc_spline_bin_centers, lc_spline_sv, lc_spline_sv_err
#
#
# def plotLightCurve_fluxvsMJD_XandGray(
#         fDataDict, mjd_min, mjd_max,
#         icrc2019Plots=False,
#         yaxis_min = -0.9, yaxis_max=7.5,
#         fColorDict=None,
#         plot_title=None,
#         convert_erg=False):
#     """plot flux vs MJD with two different y-axis
#
#     note: fixed limits for the y-axis
#     """
#
#     lightCurvePlottingUtilities.paper_figures(4, 4, 1, False)
#     # quick and dirty fix to get consistent colors and markers
#     if len(fColorDict) == 5:
#         colors = lightCurvePlottingUtilities.getColorList(3)
#         colors[-1] = 'black'
#         colors.append( 'gray' )
#     else:
#         colors = lightCurvePlottingUtilities.getColorList(len(fColorDict))
#     markers = lightCurvePlottingUtilities.getMarkerList(icrc2019Plots)
#     markers[3]='+'
#     markers[4]='x'
#
#     fig, ax1 = plt.subplots()
#     ax1_y = ax1.twinx()
#     ax1_x = ax1.twiny()
#     if convert_erg:
#         ax1.set_ylim(ymin=yaxis_min,ymax=yaxis_max*1.8)
#     else:
#         ax1.set_ylim(ymin=yaxis_min,ymax=yaxis_max)
#     ax1_y.set_ylim(ymin=yaxis_min,ymax=yaxis_max)
#
#     c = 0
#     for key, fData in fDataDict.items():
#         print('Plotting %s (%d data points)' % (key, len(fData)))
#
#         # scale all data by 1.e-12
#         # (scale is added to the legend text)
#         y = np.asarray(fData['flux']) / 1.e-12
#         y_err = np.asarray(fData['flux_err']) / 1.e-12
#
#         if key in fColorDict:
#             c = fColorDict[key]
#
#         if len(y) == 0:
#             c += 1
#             continue
#
#         # left labelled data
#         if "VERITAS" in key or \
#                 "HESS" in key or \
#                 "MAGIC" in key:
#             ln1 = ax1.errorbar(
#                 fData['MJD'],
#                 y,
#                 y_err,
#                 fData['MJD_err'],
#                 color=colors[c],
#                 marker=markers[c],
#                 label=key,
#                 linestyle='none',
#                 fillstyle='full',
#                 linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                 markersize=lightCurvePlottingUtilities.getMarkerSize())
#             ax1_y.plot(np.nan, '-r',
#                      color=colors[c], marker=markers[c], label=key,
#                      linestyle='none', fillstyle='full',
#                      linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                      markersize=lightCurvePlottingUtilities.getMarkerSize())
#         # right labelled data
#         else:
#             pLabel=key
#             if key.find('XRT')>=0:
#                 pLabel="$\it{Swift}$-XRT"
#             ln2 = ax1_y.errorbar(
#                 fData['MJD'],
#                 y,
#                 y_err,
#                 fData['MJD_err'],
#                 color=colors[c],
#                 marker=markers[c],
#                 label=pLabel,
#                 linestyle='none',
#                 fillstyle='full',
#                 linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                 markersize=lightCurvePlottingUtilities.getMarkerSize())
#         c += 1
#
#     ax1.axis(xmin=mjd_min,xmax=mjd_max)
#     ax1.locator_params(axis='x', tight=True, nbins=4)
#     ax1.set_xlabel('Modified Julian Day (MJD)')
#     # left label: assume always gamma rays
#     ax1.set_ylabel(lightCurvePlottingUtilities.getFluxAxisString("VERITAS",
#                                                                  None,
#                                                                  "10$^{-12} \\times$",
#                                                                  convert_erg))
#
#     if plot_title:
#         ty=0.9*yaxis_max
#         if convert_erg:
#             ty*=1.8
#         ax1.text(
#               mjd_min+0.75*(mjd_max-mjd_min),
#               ty,
#               plot_title,
#               fontsize=10,
#               weight='bold',
#               )
#
#     # right label: assume always X-rays
#     if 'Swift XRT' in fDataDict or 'NuSTAR' in fDataDict:
#         ax1_y.set_ylabel(lightCurvePlottingUtilities.getFluxAxisString(
#   "Swift XRT", None, "10$^{-12} \\times$"))
#
#     ax1_y.legend(loc=2)
#     if yaxis_min < 0.:
#         plt.axhline(y=0, linestyle=':')
#     lightCurvePlottingUtilities.align_yaxis(ax1, 0, ax1_y, 0)
#
#     # orbital phase axis
#     # (note: use default orbital period here
#     # as defined in lightCurveAnalysisorbital_period
#     phase_ticks = []
#     label_format = '%.2f'
#     for t in ax1_x.get_xticks():
#         mjd=mjd_min+t*(mjd_max-mjd_min)
#         phase_ticks.append(label_format % (
#             lightCurveAnalysisorbital_period.getOrbitalPhase(mjd),))
#     # remove first and last label from plotting
#     if len(phase_ticks)>0:
#         phase_ticks[0]=None
#         phase_ticks[-1]=None
#     ax1_x.set_xticklabels(phase_ticks)
#     ax1_x.set_xlabel('orbital phase')
#
#     lightCurvePlottingUtilities.printFigure('XG-HESSJ0632p057-LC')
#
