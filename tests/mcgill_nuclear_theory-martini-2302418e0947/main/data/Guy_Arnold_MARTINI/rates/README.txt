set the values in Constants.h to the ones that were used to generate the files with the transition rates.
Transition rates are generated with transitionRate.nb in ~/cvs/eloss/mathematica/transitionRate.nb
use "Export transition rates dGamma/domegadt(E,T,omega) on log grid (Log_e) method A" and B in there.
That runs for some time.
Note that when you change the coupling you have to generate new transition rate files, since the coupling cannot be scaled out.
However, the temperature scales out. I usually use 400 MeV in the data files - just set the value to the temperature you used when generating you file to the same on in Constants.h and everything should be fine.

now used data files are
logEntr...
logEnDtr..

and to compare 2q*_{at minimum} and 0.5q*_{at minimum}
use 
logEnq2tr..
logEnq05tr..
respectively

the T200 ones were just used for plotting a nice transition rate in the paper for T=200 MeV
