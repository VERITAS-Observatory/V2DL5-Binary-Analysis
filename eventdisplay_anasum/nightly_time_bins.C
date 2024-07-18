/*
 * Tools for VERITAS light curve analysis with Eventdisplay
 * Prints a list of nightly time bins and the number of runs in each bin.
 *
 */

R__LOAD_LIBRARY($EVNDISPSYS/lib/libVAnaSum.so);

#include <iostream>
#include <string>
#include <set>

using namespace std;

/*
 * Print nightly time bins and the number of runs in each bin.
 * Used for the light curve generation
 *
 */
void nightly_time_bins( string iAnaSumFile )
{
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
     t->SetBranchAddress( "runOn", &runOn );
     t->SetBranchAddress( "MJDOn", &MJDOn );

     // loop over all runs
     set< unsigned int > MJD;
     for( int i = 0; i < t->GetEntries(); i++ )
     {
          t->GetEntry( i );
          if( runOn < 0 ) continue;
          MJD.insert( (int)MJDOn);
     }
     for (auto it = MJD.begin(); it != MJD.end(); ++it)
     {
          std::cout << *it << " " << *it + 1 << endl;
     }
}
