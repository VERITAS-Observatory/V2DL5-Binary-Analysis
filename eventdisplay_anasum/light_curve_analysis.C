/*
 * Tools for VERITAS light curve analysis with Eventdisplay
 *
 *  - light_curve_analysis : print a light curve line (ECSV format) for each period
 *  - TODO !! checkLightCurveIntervals - check that all runs are in MJD intervals
 *
 */

R__LOAD_LIBRARY($EVNDISPSYS/lib/libVAnaSum.so);

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace std;

/*
 * Read MJD min/max for flux calculation from file.
 *
 * File with MJD intervals: text file with <mjd start> <mjd stop>
 *
 */
vector<pair<double, double>> read_MJD_from_file(string iMJDIntervalFile)
{
     vector<pair<double, double>> MJD_min_max;

     ifstream is;
     is.open( iMJDIntervalFile.c_str(), ifstream::in );
     if( !is )
     {
          cout << "Error: MJD interval file " << iMJDIntervalFile << " not found" << endl;
          return MJD_min_max;
     }
     string is_line;
     while( getline( is, is_line ) )
     {
         if( is_line.size() > 0 )
         {
             istringstream is_stream( is_line );
              string t_min;
              string t_max;
              is_stream >> t_min;
              is_stream >> t_max;
              MJD_min_max.emplace_back( atof( t_min.c_str() ), atof( t_max.c_str() ) );
         }
     }
     is.close();
     return MJD_min_max;
}

/*
 * Read MJD start and stop time for each run.
 */
vector< pair<double, double> > read_runlist_from_file(string iAnaSumFile)
{
     vector< pair<double, double> > runlist;

     TFile *f = new TFile( iAnaSumFile.c_str() );
     if( f->IsZombie() )
     {
         return runlist;
     }
     TTree *t = (TTree*)f->Get("total_1/stereo/tRunSummary" );
     if( !t )
     {
         return runlist;
     }
     int runOn;
     double MJDrunstart;
     double MJDrunstop;
     t->SetBranchAddress( "runOn", &runOn );
     t->SetBranchAddress( "MJDrunstart", &MJDrunstart );
     t->SetBranchAddress( "MJDrunstop", &MJDrunstop );

     // add one 'safety' minute at run start and end
     double add_one_minute = 1./60./24.;

     for( int i = 0; i < t->GetEntries(); i++ )
     {
         t->GetEntry( i );
         if( runOn > 0 )
         {
            runlist.emplace_back( MJDrunstart-add_one_minute, MJDrunstop+add_one_minute );
         }
     }
     return runlist;
}

/*
 * print light curves for given MJD interval file
 *
 * **major functions, hardwired values are agreed**
 *
 * File with MJD intervals: text file with <mjd start> <mjd stop>
 *
 */
void light_curve_analysis(
        string iAnaSumFile = "anasum/anasum.combined.root",
        string iMJDIntervalFile = "HESSJ0632p057.MJD.v3.txt" )
{
     // energy threshold in TeV
     double i_fixed_Emin = 0.35;
     // assumed spectral index:
     double i_fixed_Index = -2.6;

     cout << "# Light curve generation ";
     cout << "E_min > " << i_fixed_Emin << " TeV; ";
     cout << "Spectral index " << i_fixed_Index << endl;

     vector<pair<double, double>> min_max;
     size_t n_intervals = 0;
     if( iMJDIntervalFile.size() > 0 && iMJDIntervalFile.find("RUNWISE") == string::npos )
     {
          min_max = read_MJD_from_file(iMJDIntervalFile);
     }
     else
     {
          min_max = read_runlist_from_file(iAnaSumFile);
     }
     for(unsigned int i = 0; i < min_max.size(); i++ )
     {
          cout << "MJD interval " << i << ": ";
          cout << setprecision(12) << min_max[i].first << " - " << min_max[i].second << endl;
     }

     for(unsigned int i = 0; i < min_max.size(); i++ )
     {
          VFluxCalculation a(
               iAnaSumFile.c_str(),
               1, -1, -1,
               min_max[i].first,
               min_max[i].second
          );
          a.setSignificanceParameters( -5., -10. );
          a.setSpectralParameters( i_fixed_Emin, 1., i_fixed_Index );
          a.calculateIntegralFlux( i_fixed_Emin );
          a.printECSVLine();
     }
}
