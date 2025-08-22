
/**************************************************************************
 *                                                                        *
 *    AMY Energy Loss                                                     *
 *                                                                        *
 *    Bjoern Schenke, McGill University 2008, schenke@physics.mcgill.ca   *
 *                                                                        *
 **************************************************************************/

/** ------------------------------------------------------------------------
    precompiler checks for config.h and includes it:                        
    ------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/** ------------------------------------------------------------------------
    further header files to include:                                        
    ------------------------------------------------------------------------*/

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string>

//#include "mpi.h"
#include "main/Setup.h"

using namespace std;

/** ------------------------------------------------------------------------
    the main program:
    ------------------------------------------------------------------------ */

int main(int argc, char* argv[]) 
{

/** ------------------------------------------------------------------------
    initialize MPI for parallel computing (not included yet)
    ------------------------------------------------------------------------ */

//  MPI::Init(argc, argv);
//  int rank = MPI::COMM_WORLD.Get_rank();
//  int size = MPI::COMM_WORLD.Get_size();

/** ------------------------------------------------------------------------
    let people know what's running:
    ------------------------------------------------------------------------ */

      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;    
      cout << "  energy loss - Version " << VERSION << endl;
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
    
/** ------------------------------------------------------------------------
    create setup - object, which will automatically read in the parameters:
    ------------------------------------------------------------------------ */

      Setup* setup;
      setup = new Setup(argc,argv,"param.inp");

/** ------------------------------------------------------------------------
    create control - object, using setup that also passes parameters to it:
    ------------------------------------------------------------------------ */

      Control* control = setup->createControl();

/** ------------------------------------------------------------------------
    create and initialize particle list according to the input parameters:
    ------------------------------------------------------------------------ */

      FileNames* fileNameList = setup->getFileNameList();

/** ------------------------------------------------------------------------
    run evolution:
    ------------------------------------------------------------------------ */
      
      control->evolve();
}

