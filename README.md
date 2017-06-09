# sigomk1 (short for sigma, omega, rho - Mk. 1):
# ==============================================================================================================
sigomk1: Codes for calculating the equation of state for dense nuclear matter in a relativistic mean field (RMF) formalism
The main file is SigOmRho.cc, with the headers SigOmRhoFunc.hh (holds the relevant functions to be solved/integrated), 
SigOmRhoInt.hh (holds a seperate integrator for each recurring integral) and CONSTANTSandCONVERSIONS.hh (my usual header that
contains the relevant conversions for working with natural units or those of GeV^-1 or fm (distance, time) and GeV or fm^-1 
(energy)). This code only has the minimal particle content of a nucleon doublet and electrons (for charge neutrality). Future
versions of the code will have muons (no taus, as they are too unstable and are hardly populated) and hyperons.[1] There are also
early iterations of my codes for DM capture and for calculating the effect of DM on the EoS of dense nuclear matter. [2]

This code was written following the analytic calculations in Compact Stars (second edition) by N. K. Glendenning.

* [1] For these updated EoS codes, see https://github.com/pilambdaepsilon/dminns
* [2] The updated versions of these DM related codes can be found in https://github.com/pilambdaepsilon/dminns
