/*
 * Tools for VERITAS spectral analysis with Eventdisplay
 *
 * settings very similar to 2022 HESS J0632+057 paper
 *
 */

R__LOAD_LIBRARY($EVNDISPSYS/lib/libVAnaSum.so);

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

using namespace std;

/*
 * read configuration parameters from file
 */
map<string,double> read_config(string config_file) {
    map<string,double> params;
    ifstream fin(config_file.c_str());
    if (!fin.is_open()) {
        cerr << "Error opening config file: " << config_file << endl;
        exit( -1 );
    }

    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string key;
        double value;
        getline(ss, key, ':');
        ss >> value;

        // Remove whitespace from key
        key.erase(remove_if(key.begin(), key.end(), ::isspace), key.end());
        params[key] = value;
        cout << "Key: " << key << ", Value: " << value << endl;
    }
    return params;
}

/*
 * spectral analysis
 *
 * Typical settings for energy binning are 0.2 or 0.3 (in equal interval on the log-energy axis)
*/
void spectral_analysis(
    string anasumfile = "",
    string config_file = "spectral_plotting_config.txt",
    string csv_output_file = "")
{
    map<string,double> config = read_config(config_file);

    VEnergySpectrum a(anasumfile);
    a.setSignificanceParameters(1., 1., 0.95, 17, 4);
    a.setSpectralFitFluxNormalisationEnergy(config["DECORRELATIONENERGY_TEV"]);
    a.setEnergyBinning(config["ENERGYBINNING"]);
    a.setPlottingYaxis(config["FLUX_MIN"], config["FLUX_MAX"]);
    a.setPlottingEnergyRangeLinear(config["ENERGY_MIN"], config["ENERGY_MAX"]);
    TCanvas *c = a.plot();

    double iEMax_lin_TeV = a.getUpperEdgeofLastFilledEnergyBin( 0., 1. );
    double iEMin_lin_TeV = a.getLowerEdgeofFirstFilledEnergyBin( 0., 1. );

    a.setSpectralFitRangeLin(iEMin_lin_TeV, iEMax_lin_TeV);
    // power law fit
    a.setSpectralFitFunction(0);
    TF1 *fFitFunction = a.fitEnergySpectrum();

    // plotting of spectra plus fit
    a.setSignificanceParameters(1., 1., 0.95, 17, 4);
//    a.setPlottingUpperLimits(false);
    c = a.plot();
    TGraphAsymmErrors* i_cl = a.getEnergySpectrumGraph();
    if (fFitFunction != 0)
    {
        fFitFunction->Draw("same");
    }
    a.plotFitValues();
    if (c)
    {
        c->Print((csv_output_file + ".pdf").c_str());
    }

    // write to spectral points
    a.writeSpectralPointsToCSVFile((csv_output_file + ".ecsv").c_str());
}
