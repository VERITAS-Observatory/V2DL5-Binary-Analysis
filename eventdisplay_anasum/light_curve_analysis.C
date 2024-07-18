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
#include <vector>

using namespace std;

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

     ifstream is;
     is.open( iMJDIntervalFile.c_str(), ifstream::in );
     if( !is )
     {
          cout << "Error: MJD interval file " << iMJDIntervalFile << " not found" << endl;
          return;
     }
     string is_line;
     while( getline( is, is_line ) )
     {
         if( is_line.size() > 0 )
         {
             istringstream is_stream( is_line );
              string t;
              is_stream >> t;
              double MJD_min = atof( t.c_str() );
              is_stream >> t;
              double MJD_max = atof( t.c_str() );

              VFluxCalculation a( iAnaSumFile.c_str(), 1, 0, 1000000, MJD_min, MJD_max );
              a.setSignificanceParameters( -5., -10. );
              a.setSpectralParameters( i_fixed_Emin, 1., i_fixed_Index );
              a.calculateIntegralFlux( i_fixed_Emin );
              a.printECSVLine();
         }
     }

     is.close();
}

/*
 * check that all runs are in one of the intervals in MJD
 * used for the light curve generation
 *
 */
void checkLightCurveIntervals( string iAnaSumFile, string iMJDIntervalFile, bool iWriteRunLists = false )
{
     // read time ranges in MJD
     vector< double > iMJD_min;
     vector< double > iMJD_max;
     ifstream is;
     is.open( iMJDIntervalFile.c_str(), ifstream::in );
     if( !is )
     {
          cout << "Error: MJD interval file not found" << endl;
          return;
     }
     string is_line;
     while( getline( is, is_line ) )
     {
         if( is_line.size() > 0 )

         {
             istringstream is_stream( is_line );
             string td;
             is_stream >> td;
             iMJD_min.push_back( atof( td.c_str() ) );
             is_stream >> td;
             iMJD_max.push_back( atof( td.c_str() ) );
         }
     }
     cout << "Total number of intervals: " << iMJD_min.size() << endl;
     vector< int > iNruns( iMJD_min.size(), 0. );
     vector< double > iT( iMJD_min.size(), 0. );

     // runs not in any interval
     vector< int > iMissingRun;


     // open anasum file
     TFile *f = new TFile( iAnaSumFile.c_str() );
     if( f->IsZombie() )
     {
         return;
     }
     TTree *t = (TTree*)f->Get("total_1/stereo/tRunSummary" );
     if( !t )
     {
         return;
     }
     int runOn;
     double MJDOn;
     double tOn;
     t->SetBranchAddress( "runOn", &runOn );
     t->SetBranchAddress( "MJDOn", &MJDOn );
     t->SetBranchAddress( "tOn", &tOn );

     // list of runs per time bin
     vector< vector< int > > runsPerTimeBin;

     // loop over all time invervalls
     for( unsigned n = 0; n < iMJD_min.size(); n++ )
     {
         bool bFound = false;

         vector< int > iRunList;
         // loop over all runs
         for( int i = 0; i < t->GetEntries(); i++ )
         {
             t->GetEntry( i );

             if( runOn < 0 ) continue;

             if( MJDOn > iMJD_min[n] && MJDOn < iMJD_max[n]  )
             {
                iNruns[n]++;
                iT[n] += tOn;
                bFound = true;
                iRunList.push_back( runOn );
             }
         }
         if( !bFound )
         {
              iMissingRun.push_back( runOn );
              cout << "Missing run " << runOn << "\t" << MJDOn << endl;
         }
         runsPerTimeBin.push_back( iRunList );

     }

     cout << "Total number of missing runs: " << iMissingRun.size() << endl;

     for( unsigned int ii = 0; ii < iMissingRun.size(); ii++ )
     {
          cout << "Missing run " << iMissingRun[ii] << endl;
     }

     cout << "Summary: " << endl;
     for( int i = 0; i < iMJD_min.size(); i++ )
     {
          cout << "Interval " << iMJD_min[i] << ", " << iMJD_max[i] << ": ";
          cout << iNruns[i] << " runs, ";
          cout << iT[i]/3600. << " h" << endl;
     }

    //////////////////////////////
    // write run lists
    if( iWriteRunLists )
    {
         for( unsigned int n = 0; n < runsPerTimeBin.size(); n++ )
         {
              ostringstream iFileName;
              iFileName << "MJD" << iMJD_min[n] << "-MJD" << iMJD_max[n] << ".dat";
              ofstream iRunListFile;
              iRunListFile.open( iFileName.str().c_str() );
              cout << "Writing run list " << iFileName.str();
              cout << " with " << runsPerTimeBin[n].size() << " runs" << endl;
              for( unsigned int m = 0; m < runsPerTimeBin[n].size(); m++ )
              {
                   iRunListFile << runsPerTimeBin[n][m] << endl;
              }
              iRunListFile.close();
         }
    }

    //////////////////////////////
    // plot delta t vs MJD
    TGraph *iG_DT = new TGraph( 1 );
    iG_DT->SetTitle( "" );
    for( int i = 0; i < iMJD_min.size(); i++ )
    {
        iG_DT->SetPoint( i, 0.5*(iMJD_min[i]+iMJD_max[i]), (iMJD_max[i]-iMJD_min[i]) );
    }
    iG_DT->SetMarkerStyle( 21 );

    TCanvas *c1 = new TCanvas( "cDT", "delta T vs MJD", 10, 10, 600, 400 );
    c1->Draw();
    iG_DT->Draw( "ap" );
    iG_DT->GetHistogram()->SetXTitle( "MJD" );
    iG_DT->GetHistogram()->SetYTitle( "delta MJD" );
}
