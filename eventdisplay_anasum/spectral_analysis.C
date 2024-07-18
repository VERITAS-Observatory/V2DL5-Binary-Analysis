/*
 * Tools for VERITAS spectral analysis with Eventdisplay
 *
 * settings very similar to 2022 HESS J0632+057 paper
 *
 */

R__LOAD_LIBRARY($EVNDISPSYS/lib/libVAnaSum.so);

#include <string>

using namespace std;

/*
 * spectral analysis
 *
 * Typical settings for energy binning are 0.2 or 0.3 (in equal interval on the log-energy axis)
*/
void spectral_analysis(
    string anasumfile = "",
    double iEnergyBinning = 0.2,
    string csv_output_file = "",
    double iDeCorrelationEnergy_TeV = 0.5)
{
    VEnergySpectrum a( anasumfile );
    a.setSignificanceParameters( 1., 1., 0.95, 17, 4 );
    a.setSpectralFitFluxNormalisationEnergy( iDeCorrelationEnergy_TeV );
    a.setEnergyBinning( iEnergyBinning );
    a.setPlottingYaxis( 1.e-15, 1e-10 );
    TCanvas *c = a.plot();

    double iEMax_lin_TeV = a.getUpperEdgeofLastFilledEnergyBin( 0., 1. );
    double iEMin_lin_TeV = a.getLowerEdgeofFirstFilledEnergyBin( 0., 1. );

    a.setSpectralFitRangeLin( iEMin_lin_TeV, iEMax_lin_TeV );
    // power law fit
    a.setSpectralFitFunction( 0 );
    TF1 *fFitFunction = a.fitEnergySpectrum();

    // plotting of spectra plus fit
    a.setSignificanceParameters( 1., 1., 0.95, 17, 4 );
//    a.setPlottingUpperLimits( false );
    c = a.plot();
    TGraphAsymmErrors* i_cl = a.getEnergySpectrumGraph();
    if( fFitFunction != 0 )
    {
        fFitFunction->Draw( "same" );
    }
    a.plotFitValues();
    if( c )
    {
        c->Print( (csv_output_file+".pdf").c_str());
    }

    // write to spectral points
    a.writeSpectralPointsToCSVFile( (csv_output_file+".ecsv").c_str());
}
