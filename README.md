# MARTINI #
MARTINI is a Modular Algorithm for Realistic Treatment of Heavy Ion Collisions.
Article: [Phys.Rev.C 80 (2009) 054913] 
## General ##
The software implements the AMY-McGill energy loss formalism for hard jets as they 
traverse a QGP medium. It expects two extrnal software packages, LHAPDF and PYTHIA.
PYTHIA used by MARTINI is modified: it has some functions implemented by Chanwook Park
(currently post doc at Wayne State-June 2021) and Bjoern Schekne (currently Staff Scientist
at Brookhaven National Lab) for calculation of jet cross section and other modifications for
nuclear PDF calculations. These packages are not included in this git repository but they 
should be available on the Beluga system for McGill-Nuclear Theory students. ROOT(by CERN)
is also required as a package for some histogramming functionality. If one wants to perform
jet analysis, the FASTJET3 is the most common choice which should also be provided by the 
user(s) themself(ves). 

This version is the result of Chanwook Park's PhD work and includes the addition of 
recoil partons, a NLO running of the coupling constant as well as formation time for 
suppression of radiation. Rouzbeh Yazdi has also added a perturbative calculation
 functionality for Conversion(collinear) and Bremsstrahlung photons. RY will also add functionality for
calculation of conversion photons beyond the collinear approximation.

Slightly longer term view is to include Dileptons from jet conversions (by RY) as well as
conversion photons from gluon jets. 

## Rates ##
MARTINI gives a choice of radiation rate sets:
    * LO-AMY with a small q cut
    * LO-AMY [cite AMY papers]
    * NLO-AMY NLO kernel [cite the NLO paper]
    * NP-AMY Non-perturbative kernel [cite the NP rate paper]

## Required Packages ##
    * PYTHIA 8.209 (local on Beluga, modified version)
    * LHAPDF5 for Parton Distribution Functions
    * ROOT (requirement removed for now. )
    * FASTJET3 (if jet analysis is the goal: should enter only in your own "main.cpp" file)

## Code Standard Suggestions ##
    * Use 4 space tabs 
    * code in C++11 standard 
    * place a space after statments like `if` etc
    * Curly brackets must be defined on a new line. 
### Contact ###
For help or more information, the current code master is Rouzbeh Yazdi and can be contacted
at rouzbeh.modarresi-yazdi@mail.mcgill.ca
