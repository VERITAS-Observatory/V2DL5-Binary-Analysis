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
     vector<pair<double, double>> mjd_min_max;

     ifstream is;
     is.open( iMJDIntervalFile.c_str(), ifstream::in );
     if( !is )
     {
          cout << "Error: MJD interval file " << iMJDIntervalFile << " not found" << endl;
          return mjd_min_max;
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
              mjd_min_max.emplace_back( atof( t_min.c_str() ), atof( t_max.c_str() ) );
         }
     }
     is.close();
     return mjd_min_max;
}

/*
 * check that all runs are in one of the intervals in MJD
 * used for the light curve generation
 *
 */
void runlist_from_time_bins( string iAnaSumFile, string iMJDIntervalFile, bool iWriteRunLists = false )
{
     vector<pair<double, double>> min_max = read_MJD_from_file(iMJDIntervalFile);
     cout << "Total number of intervals: " << min_max.size() << endl;
     vector< int > iNruns( min_max.size(), 0. );
     vector< double > iT( min_max.size(), 0. );

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

     // loop over all time intervals
     for( unsigned n = 0; n < min_max.size(); n++ )
     {
         vector< int > iRunList;
         // loop over all runs
         for( int i = 0; i < t->GetEntries(); i++ )
         {
             t->GetEntry( i );

             if( runOn < 0 ) continue;

             if( MJDOn > min_max[n].first && MJDOn < min_max[n].second )
             {
                iNruns[n]++;
                iT[n] += tOn;
                iRunList.push_back( runOn );
             }
         }
         runsPerTimeBin.push_back( iRunList );
     }
     // check for missing runs
     for( int i = 0; i < t->GetEntries(); i++ )
     {
         t->GetEntry( i );
         if( runOn < 0 )
         {
             continue;
         }
         bool bFound = false;
         for(unsigned int t = 0; t < runsPerTimeBin.size(); t++ )
         {
             for(unsigned s = 0; s < runsPerTimeBin[t].size(); s++ )
             {
                 if( runOn == runsPerTimeBin[t][s] )
                 {
                     bFound = true;
                     break;
                 }
             }
         }
         if( !bFound )
         {
              iMissingRun.push_back( runOn );
              cout << "Missing run " << runOn << "\t" << MJDOn << endl;
         }
    }




     cout << "Total number of missing runs: " << iMissingRun.size() << endl;

     for( unsigned int ii = 0; ii < iMissingRun.size(); ii++ )
     {
          cout << "Missing run " << iMissingRun[ii] << endl;
     }

     cout << "Summary: " << endl;
     for( int i = 0; i < min_max.size(); i++ )
     {
          cout << "Interval " << min_max[i].first << ", " << min_max[i].second << ": ";
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
              iFileName << "MJD" << min_max[n].first;
              iFileName << "-MJD" << min_max[n].second << ".dat";
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
    for( int i = 0; i < min_max.size(); i++ )
    {
        iG_DT->SetPoint(
           i,
           0.5*(min_max[i].first+min_max[i].second),
          (min_max[i].second-min_max[i].first) );
    }
    iG_DT->SetMarkerStyle( 21 );

    TCanvas *c1 = new TCanvas( "cDT", "delta T vs MJD", 10, 10, 600, 400 );
    c1->Draw();
    iG_DT->Draw( "ap" );
    iG_DT->GetHistogram()->SetXTitle( "MJD" );
    iG_DT->GetHistogram()->SetYTitle( "delta MJD" );
    c1->SaveAs("delta_t_vs_MJD.png");
}
