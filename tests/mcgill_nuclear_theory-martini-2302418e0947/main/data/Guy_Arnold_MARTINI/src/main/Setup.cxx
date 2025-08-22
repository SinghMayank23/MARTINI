
#ifdef HAVE_CONFIG_H
#include <config.h> //include to have access to VERSION
#endif
#include "Setup.h"

/// Setup reads in parameters and sets up output files

using namespace std;

Setup::Setup(int argc, char** argv, string f, const int size, const int rank)
{

/** ------------------------------------------------------------------------
    open and read parameter file:
    ------------------------------------------------------------------------ */

  inputFile = f;
  
  for(int i=1; i<argc; i++) 
    {
      if(!strcmp(argv[i],"-f")) inputFile = argv[i+1];
    }
  
  ifstr.open(inputFile.c_str(),ios::in);
  if(!ifstr) 
    {
      cerr << "[Setup::Setup]: ERROR: Unable to open file " << inputFile << endl;
      exit(1);
    }
 
  fileNameList = new FileNames(); 
  readparam = new ReadParameter(inputFile);
 
  readparam->read();
  ifstr.close();

  readparam->inputParam("output_file_name_1",tempFileName[1]);
  readparam->inputParam("output_file_name_2",tempFileName[2]);
  readparam->inputParam("output_file_name_3",tempFileName[3]);
  readparam->inputParam("output_file_name_4",tempFileName[4]);
  readparam->inputParam("trqq_file",tempFileName[5]);
  readparam->inputParam("trqg_file",tempFileName[6]);
  readparam->inputParam("trgq_file",tempFileName[7]);
  readparam->inputParam("trgg_file",tempFileName[8]);
  readparam->inputParam("trrad_file",tempFileName[9]);
  readparam->inputParam("temperature",T);
  readparam->inputParam("alpha_s",alphas);
  readparam->inputParam("timestep",dt);
  readparam->inputParam("energy_resolution",dE);
  readparam->inputParam("Nf",Nf);
  readparam->inputParam("initial_energy",initialE);
  readparam->inputParam("initial_width",width);
  readparam->inputParam("max_time",maxTime);
  readparam->inputParam("max_energy",maxEnergy);
  readparam->inputParam("smooth_transition_rates",strate);
  readparam->inputParam("bethe_heitler",BetheHeitler);
  readparam->inputParam("bdmps",BDMPS);
  readparam->inputParam("do_rad",doRad);
  readparam->inputParam("do_col",doCol);
  readparam->inputParam("make_new_radiative_gamma",newRadGamma);
  readparam->inputParam("measurements_per_fermi",mperfermi);

  // replace dt by dt/hbarc, so that the used dt is actually in GeV^-1.
  dt=dt/0.197327053;
  
  for(int i=1;i<=9;i++) //set max to maximal file number given above...
    {
      fileNameList->SetFN(tempFileName[i]);
    }
  
/** ------------------------------------------------------------------------
    initialize new instance of constants class and fill in the ones read in:
    ------------------------------------------------------------------------ */
  
  constants = new Constants();



/** ------------------------------------------------------------------------
    write information on used parameters in output file:
    ------------------------------------------------------------------------ */

  if(rank==0)
    {
      cout << "All information on used parameters etc. can be found in: " << tempFileName[1] << "." << endl;
      
      // This is the first time output is written into this file,
      // hence use ios::out and not ios::app to overwrite any exisiting data:  
      
      fstream fout(fileNameList->GetFN(1).c_str(),ios::out); 
      fout << "[Setup::Setup]: This is AMY energy loss Version " << VERSION << endl;
      fout << "                Names in [] indicate which function is writing the following output." << endl;
      fout << "                This is helpful for development and error tracking." << endl;    
      fout << "                _________________________________________" << endl;
      fout << "                Output is written to the following files:" << endl;
      fout << "                1) " << fileNameList->GetFN(1).c_str() << " (this file)" << endl;
      fout << "                2) " << fileNameList->GetFN(2).c_str() << endl;
      fout << "                3) " << fileNameList->GetFN(3).c_str() << endl;
      fout << "                4) " << fileNameList->GetFN(4).c_str() << endl;
//--> add information on additional output files here!

      fout << "                _________________________________________" << endl;
      fout << "                Used input files with transition rates:  " << endl;
      fout << "                1) " << "./rates/" << fileNameList->GetFN(5).c_str() << endl;
      fout << "                2) " << "./rates/" << fileNameList->GetFN(6).c_str() << endl;
      fout << "                3) " << "./rates/" << fileNameList->GetFN(7).c_str() << endl;
      fout << "                4) " << "./rates/" << fileNameList->GetFN(8).c_str() << endl;
      fout << "                5) " << "./rates/" << fileNameList->GetFN(9).c_str() << endl;
      fout << "                _________________________________________" << endl;
      fout << endl;
      fout.close();
      
    }
}

Setup::~Setup()
{
    delete readparam;
    delete fileNameList;
}

FileNames* Setup::getFileNameList()
{
  return fileNameList;
}

Control* Setup::createControl()
{    
  Control* c = new Control(fileNameList, constants, dt, dE, T, alphas, Nf, initialE, width, maxTime, maxEnergy, strate, BetheHeitler, BDMPS, doRad, doCol, newRadGamma, mperfermi, setupsSize, setupsRank);
  return c;
}
