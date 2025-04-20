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
 read a list of runs from a text file
*/
vector< int > read_run_list(string run_list) {
    vector< int > runs;
    if( run_list.size() == 0 )
    {
        return runs;
    }
    ifstream fin(run_list.c_str());
    if (!fin.is_open()) {
        cerr << "Error opening run list file: " << run_list << endl;
        exit( -1 );
    }

    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        int run;
        ss >> run;
        runs.push_back(run);
    }
    return runs;
}

/*
 * spectral analysis
 *
*/
void spectral_analysis(
    string anasumfile = "",
    string config_file = "spectral_plotting_config.txt",
    string output_file = "",
    string run_list_file = "")
{
    map<string,double> config = read_config(config_file);
    vector< int > run_list = read_run_list(run_list_file);

    VEnergySpectrum a(anasumfile);
    a.setSignificanceParameters(1., 1., 0.95, 17, 4);
    a.setSpectralFitFluxNormalisationEnergy(config["DECORRELATIONENERGY_TEV"]);
    a.setEnergyBinning(config["ENERGYBINNING"]);
    a.setPlottingYaxis(config["FLUX_MIN"], config["FLUX_MAX"]);
    a.setPlottingEnergyRangeLinear(config["ENERGY_MIN"], config["ENERGY_MAX"]);
    if(run_list.size() > 0)
    {
        a.combineRuns(run_list);
    }
    TCanvas *c = a.plot();

    double iEMax_lin_TeV = a.getUpperEdgeofLastFilledEnergyBin( config["EXCESS_EVENTS_FIT_MIN"], config["SIGNIFICANCE_FIT_MIN"] );
    double iEMin_lin_TeV = a.getLowerEdgeofFirstFilledEnergyBin( config["EXCESS_EVENTS_FIT_MIN"], config["SIGNIFICANCE_FIT_MIN"] );

    a.setSpectralFitRangeLin(iEMin_lin_TeV, iEMax_lin_TeV);
    // power law fit
    a.setSpectralFitFunction(0);
    TF1 *fFitFunction = a.fitEnergySpectrum();

    // plotting of spectra plus fit
    a.setSignificanceParameters(config["SIGNIFICANCE_PLOT_MIN"], config["EXCESS_EVENTS_PLOT_MIN"], 0.95, 17, 4);
    a.setPlottingUpperLimits((int)config["PLOT_UPPER_LIMITS"]);
    c = a.plot();
    if( (int)config["PLOT_EVENT_NUMBERS"])
    {
        a.plotEventNumbers();
    }
    TGraphAsymmErrors* i_cl = a.getEnergySpectrumGraph();
    if (fFitFunction != 0)
    {
        fFitFunction->Draw("same");
    }
    a.plotFitValues();
    if (c)
    {
        c->Print((output_file + ".pdf").c_str());
    }

    // write spectral points
    a.writeSpectralPointsToCSVFile((output_file + ".ecsv").c_str());
}
