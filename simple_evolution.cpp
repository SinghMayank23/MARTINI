#include <iostream>
#include <fstream>
#include <math.h>
#include "MARTINI.h"
#include "TH1.h"
#include "TH2.h"
//#define PI 3.14159265359
using namespace std;


int main(int argc, char* argv[])
{
    // setup from the command line arguments
    string nevents    = argv[1]; // number of events
    string srun       = argv[2]; // subrun number
    string save_loc   = argv[3]; // save location
    string setup_file = argv[4]; // setup file: settings of this MARTINI run
    // number of events, jet pT info
    double jetPTHatBins[] = {  3., 5.  , 7.  , 9.  , 11. , 13. , 15. , 
                              17., 19. , 21. , 23. , 25. , 27. , 29. , 
                              32., 35. , 40. , 45. , 50. , 55. , 60. , 
                              65., 70. , 75. , 80. , 85. , 90. , 100., 
                             110., 120., 130., 140., 150., 160., 170., 
                             180., 190., 200., 225., 250., 0.};//makes for 40 subruns
    int numEvents      = stoi(nevents);


    // Prepare MARTINI
    MARTINI martini;
    martini.readFile(setup_file);
    int Subrun = atoi(srun.c_str());
    cout << "\n\n\t\tsubrun number : " << Subrun <<endl;
    cout << "\t\t\t --> pTMin: "<<jetPTHatBins[Subrun]<<" ==> pTMax: "<< jetPTHatBins[Subrun+1]<<endl;
    martini.readString("General:JetPTMin = " + to_string(jetPTHatBins[Subrun]));
    martini.readString("General:JetPTMax = " + to_string(jetPTHatBins[Subrun+1]));

    martini.init(0);
    martini.settings.listChanged(); // list the changed parameters. 
    
    vector<Parton> * plist = new vector<Parton>; //pointer to the parton list vector
    vector<Source> * slist = NULL; // pointer to the source list, null for now

    int event_counter = 0; 
    int mt; // maximal time steps
    double maxTime = martini.returnMaxTime();
    double dtfm = martini.returnDtfm();// dt in femtometers 
    mt = static_cast<int>(maxTime/dtfm+0.0001); //max number of steps
    int counter;
    int hadronized = 0;
    while( event_counter < numEvents )
    {
        plist->clear();
        martini.generateEvent(plist);
        Event fullEvent;
        counter = 0;
        if (martini.returnEvolution() == 1)  // evolve in medium if settings want that
        {
            for(int i=0; i<mt; i++) // loop over all time steps 
            {
                counter = martini.evolve(plist, slist, counter, i);
                counter+=1;
            }
        }
        // now do fragmentation and bin neutral pions and charged hadrons
        fullEvent.clear();
        martini.pythia.event.clear();

        if( martini.fragmentation(plist) ) // if you successfully fragment
        {
            hadronized += 1;
            // bin the particles
        }
            
        event_counter++;
    }
    // Cross Section Magic:
    double sigmaTot = martini.pythia.totalCrossSection();
    double sigmaEl = martini.pythia.elasticCrossSection();
    double sigmaInel = sigmaTot - sigmaEl;
    double sigmaGen = martini.pythia.info.sigmaGen();
    double ratio_photon = sigmaGen/(((double)numEvents)*sigmaInel);
    double ratio_hadron = sigmaGen/(((double)hadronized)*sigmaInel);
    // Scale your results by the relevant (cross section and num events) ratio
    // Save histograms to file

    cout << "numEvents : "<< numEvents << ", hadronized: "<< hadronized << " --> " << hadronized/(double)numEvents << endl;
    return 0;
};
