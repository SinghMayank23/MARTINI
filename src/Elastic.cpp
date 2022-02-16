// Elastic.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the elastic processes for energy loss and momentum broadening

#include "Elastic.h"

double Elastic::totalRate(double p, double T, double alpha_s, int Nf, int process)
{
    // compute the total transition rate in GeV (integrated over k, omega and q and angles) for fixed E=p and temperature T
    // using parametrization of numerically computed integral - then interpolate to get right value for used alpha_s 
    // processes: 1:qq, 2:gq, 3:qg, 4:gg
    // IMPORTANT: all computed values below are for a minimal omega of 0.05*T, so this is the cutoff to use in the calculation also!
    // also Nf=3 was used ... scales out though - see below
    // p is taken to be in GeV
    
    double alpha0 = 0.15;
    double deltaAlpha = 0.03;
    double iAlpha = floor((alpha_s-alpha0)/deltaAlpha+0.001);
    double alphaFrac = (alpha_s-alpha0)/deltaAlpha - iAlpha;
    double coefficient[2][6];
    double rateLower, rateUpper, rate;
    if (process == 1 || process == 2)  // qq or gq (scattering with a thermal quark)
    {
        if ( alpha_s >= 0.15 && alpha_s < 0.18 )
        {
            coefficient[0][0] = 0.18172488396136807;  coefficient[1][0] = 0.224596478395945;
            coefficient[0][1] = 0.6004740049060965;   coefficient[1][1] = 1.0874259848101948;
            coefficient[0][2] = 0.36559627257898347;  coefficient[1][2] = 0.6436398538984057;
            coefficient[0][3] = 0.10607576568373664;  coefficient[1][3] = 0.11585154613692052;
            coefficient[0][4] = 0.004322466954618182; coefficient[1][4] = -0.001719701730785056;
            coefficient[0][5] = 0.04731599462749122;  coefficient[1][5] = 0.06734745496415469;
        }
        else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
        {
            coefficient[0][0] = 0.224596478395945;    coefficient[1][0] = 0.2686436092048326;
            coefficient[0][1] = 1.0874259848101948;   coefficient[1][1] = 1.7286136256785387;
            coefficient[0][2] = 0.6436398538984057;   coefficient[1][2] = 0.9826325498183079;
            coefficient[0][3] = 0.11585154613692052;  coefficient[1][3] = 0.13136670133029682;
            coefficient[0][4] = -0.001719701730785056; coefficient[1][4] = -0.004876376882437649;
            coefficient[0][5] = 0.06734745496415469;  coefficient[1][5] = 0.09140316977554151;
        }
        else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
        {
            coefficient[0][0] = 0.2686436092048326;   coefficient[1][0] = 0.3137234778163784;
            coefficient[0][1] = 1.7286136256785387;   coefficient[1][1] = 2.445764079999846;
            coefficient[0][2] = 0.9826325498183079;   coefficient[1][2] = 1.3083241146035964;
            coefficient[0][3] = 0.13136670133029682;  coefficient[1][3] = 0.18341717903923757;
            coefficient[0][4] = -0.004876376882437649; coefficient[1][4] = 0.006098371807040589;
            coefficient[0][5] = 0.09140316977554151;  coefficient[1][5] = 0.12054238276023879;
        }
        else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
        {
            coefficient[0][0] = 0.3137234778163784;   coefficient[1][0] = 0.3597255453974444;
            coefficient[0][1] = 2.445764079999846;    coefficient[1][1] = 3.140669321831845;
            coefficient[0][2] = 1.3083241146035964;   coefficient[1][2] = 1.535549334026633;
            coefficient[0][3] = 0.18341717903923757;  coefficient[1][3] = 0.30505450230754705;
            coefficient[0][4] = 0.006098371807040589; coefficient[1][4] = 0.04285103618362223;
            coefficient[0][5] = 0.12054238276023879;  coefficient[1][5] = 0.1558288379712527;
        }
        else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
        {
            coefficient[0][0] = 0.3597255453974444;   coefficient[1][0] = 0.40656130602563223;
            coefficient[0][1] = 3.140669321831845;    coefficient[1][1] = 3.713430971987352;
            coefficient[0][2] = 1.535549334026633;    coefficient[1][2] = 1.5818298058630476;
            coefficient[0][3] = 0.30505450230754705;  coefficient[1][3] = 0.5269042544852683;
            coefficient[0][4] = 0.04285103618362223;  coefficient[1][4] = 0.11594975218839362;
            coefficient[0][5] = 0.1558288379712527;   coefficient[1][5] = 0.1982063104156748;
        }
        else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
        {
            coefficient[0][0] = 0.40656130602563223;  coefficient[1][0] = 0.45415805200862863;
            coefficient[0][1] = 3.713430971987352;    coefficient[1][1] = 4.0758813206143785;
            coefficient[0][2] = 1.5818298058630476;   coefficient[1][2] = 1.3775134184861555;
            coefficient[0][3] = 0.5269042544852683;   coefficient[1][3] = 0.873527536823307;
            coefficient[0][4] = 0.11594975218839362;  coefficient[1][4] = 0.23371456949506658;
            coefficient[0][5] = 0.1982063104156748;   coefficient[1][5] = 0.24840524848507203;
        }
        else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
        {
            coefficient[0][0] = 0.45415805200862863;  coefficient[1][0] = 0.5024541413891354;
            coefficient[0][1] = 4.0758813206143785;   coefficient[1][1] = 4.159425815179756;
            coefficient[0][2] = 1.3775134184861555;   coefficient[1][2] = 0.8719749565879445;
            coefficient[0][3] = 0.873527536823307;    coefficient[1][3] = 1.3606690530660879;
            coefficient[0][4] = 0.23371456949506658;  coefficient[1][4] = 0.4010658149846402;
            coefficient[0][5] = 0.24840524848507203;  coefficient[1][5] = 0.3067901992139913;
        }
        else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
        {
            coefficient[0][0] = 0.5024541413891354;   coefficient[1][0] = 0.5513999693402064;
            coefficient[0][1] = 4.159425815179756;    coefficient[1][1] = 3.893153859527746;
            coefficient[0][2] = 0.8719749565879445;   coefficient[1][2] = 0.009578762778659829;
            coefficient[0][3] = 1.3606690530660879;   coefficient[1][3] = 2.0095157488463244;
            coefficient[0][4] = 0.4010658149846402;   coefficient[1][4] = 0.6260756501912864;
            coefficient[0][5] = 0.3067901992139913;   coefficient[1][5] = 0.37424991045026396;
        }
        else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
        {
            coefficient[0][0] = 0.5513999693402064;   coefficient[1][0] = 0.600941593540798;
            coefficient[0][1] = 3.893153859527746;    coefficient[1][1] = 3.293344337592684;
            coefficient[0][2] = 0.009578762778659829; coefficient[1][2] = -1.1764805445298645;
            coefficient[0][3] = 2.0095157488463244;   coefficient[1][3] = 2.792180001243466;
            coefficient[0][4] = 0.6260756501912864;   coefficient[1][4] = 0.8949534049225013;
            coefficient[0][5] = 0.37424991045026396;  coefficient[1][5] = 0.44878529934031575;
        }
    
        rateLower = T*(coefficient[0][0] + coefficient[0][1]/pow(p/T,4.) - 
               coefficient[0][2]/pow((p/T),3.) - coefficient[0][3]/pow((p/T),2.) + 
               coefficient[0][4]/pow((p/T),1.5) - coefficient[0][5]/(p/T));
        rateUpper = T*(coefficient[1][0] + coefficient[1][1]/pow(p/T,4.) - 
               coefficient[1][2]/pow((p/T),3.) - coefficient[1][3]/pow((p/T),2.) + 
               coefficient[1][4]/pow((p/T),1.5) - coefficient[1][5]/(p/T));
        
        rate = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
        
        rate *= Nf/3.; // adjust number of flavors

        if (process == 2) rate *= 9./4.;
    }
    else if (process == 3 || process == 4)  // qg or gg (scattering with a thermal gluon)
    {
        if ( alpha_s >= 0.15 && alpha_s < 0.18 )
	    {
	          coefficient[0][0] = 0.9364689080337059;    coefficient[1][0] = 1.1485486950080581;
	          coefficient[0][1] = 2.626076478553979;     coefficient[1][1] = 4.993647646894147;
	          coefficient[0][2] = 2.1171556605834274;    coefficient[1][2] = 3.7295251994302876;
	          coefficient[0][3] = 0.13123339226210134;   coefficient[1][3] = -0.0017620287506503757 ;
	          coefficient[0][4] = 0.02875811664147147;   coefficient[1][4] = 0.010598257485913224;
	          coefficient[0][5] = 0.27736469898722244;   coefficient[1][5] = 0.3949856219367327;
	    }
        else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
	    {
	          coefficient[0][0] = 1.1485486950080581;    coefficient[1][0] = 1.3645568637616001;
	          coefficient[0][1] = 4.993647646894147;     coefficient[1][1] = 8.174225869366722;
	          coefficient[0][2] = 3.7295251994302876;    coefficient[1][2] = 5.732101892684938;
	          coefficient[0][3] = -0.0017620287506503757; coefficient[1][3] = -0.1416811579957863;
	          coefficient[0][4] = 0.010598257485913224;  coefficient[1][4] = 0.011703596451947428;
	          coefficient[0][5] = 0.3949856219367327;    coefficient[1][5] = 0.5354757997870718;
	    }
        else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
	    {
	          coefficient[0][0] = 1.3645568637616001;    coefficient[1][0] = 1.5839378568555678;
	          coefficient[0][1] = 8.174225869366722;     coefficient[1][1] = 11.785897000063443;
	          coefficient[0][2] = 5.732101892684938;     coefficient[1][2] = 7.758388282689373;
	          coefficient[0][3] = -0.1416811579957863;    coefficient[1][3] = -0.13163385415183002;
	          coefficient[0][4] = 0.011703596451947428;  coefficient[1][4] = 0.09016386041913003;
	          coefficient[0][5] = 0.5354757997870718;    coefficient[1][5] = 0.7042577279136836;
	    }
        else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
	    {
	        coefficient[0][0] = 1.5839378568555678;    coefficient[1][0] = 1.8062676019060235;
	        coefficient[0][1] = 11.785897000063443;    coefficient[1][1] = 15.344112642069764;
	        coefficient[0][2] = 7.758388282689373;     coefficient[1][2] =  9.384190917330093;
	        coefficient[0][3] = -0.13163385415183002;   coefficient[1][3] = 0.19709400976261568;
	        coefficient[0][4] = 0.09016386041913003;   coefficient[1][4] = 0.30577623140224813;
	        coefficient[0][5] = 0.7042577279136836;    coefficient[1][5] = 0.9066501895009754;
	    }
        else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
	    {
	        coefficient[0][0] = 1.8062676019060235;    coefficient[1][0] = 2.0312125903238236;
	        coefficient[0][1] = 15.344112642069764;    coefficient[1][1] = 18.36844006721506;
	        coefficient[0][2] = 9.384190917330093;     coefficient[1][2] = 10.209988454804193;
	        coefficient[0][3] = 0.19709400976261568;   coefficient[1][3] = 0.9957025988944573;
	        coefficient[0][4] = 0.30577623140224813;   coefficient[1][4] =  0.7109302867706849;
	        coefficient[0][5] = 0.9066501895009754;    coefficient[1][5] = 1.1472148515742653;
	    }
        else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
	    {
	        coefficient[0][0] = 2.0312125903238236;    coefficient[1][0] = 2.258502734110078;
	        coefficient[0][1] = 18.36844006721506;     coefficient[1][1] = 20.43444928479894;
	        coefficient[0][2] = 10.209988454804193;    coefficient[1][2] = 9.896928897847518 ;
	        coefficient[0][3] = 0.9957025988944573;    coefficient[1][3] = 2.3867073785159003;
	        coefficient[0][4] = 0.7109302867706849;    coefficient[1][4] = 1.3473328178504662;
	        coefficient[0][5] = 1.1472148515742653;    coefficient[1][5] = 1.429497460496924;
	    }
        else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
	    {
	        coefficient[0][0] = 2.258502734110078;     coefficient[1][0] = 2.4879110920956653;
	        coefficient[0][1] = 20.43444928479894;     coefficient[1][1] = 21.220550462966102;
	        coefficient[0][2] = 9.896928897847518;     coefficient[1][2] = 8.20639681844989;
	        coefficient[0][3] = 2.3867073785159003;    coefficient[1][3] = 4.445222616370339;
	        coefficient[0][4] = 1.3473328178504662;    coefficient[1][4] = 2.2381176005506016;
	        coefficient[0][5] = 1.429497460496924;     coefficient[1][5] = 1.7550164762706189;
	    }
        else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
	    {
	        coefficient[0][0] = 2.4879110920956653 ;   coefficient[1][0] = 2.7192501243929903;
	        coefficient[0][1] = 21.220550462966102;    coefficient[1][1] = 20.470583876561985;
	        coefficient[0][2] = 8.20639681844989 ;     coefficient[1][2] = 4.954737209403953;
	        coefficient[0][3] = 4.445222616370339;     coefficient[1][3] = 7.227667929705693;
	        coefficient[0][4] = 2.2381176005506016;    coefficient[1][4] = 3.401378906197122;
	        coefficient[0][5] = 1.7550164762706189;    coefficient[1][5] = 2.1251383942923474;
	    }
        else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
	    {
	        coefficient[0][0] = 2.7192501243929903 ;   coefficient[1][0] = 2.9523522354248817;
	        coefficient[0][1] = 20.470583876561985;    coefficient[1][1] = 18.027772799078463;
	        coefficient[0][2] = 4.954737209403953;     coefficient[1][2] = 0.050298242947981846;
	        coefficient[0][3] = 7.227667929705693;     coefficient[1][3] = 10.747352232336384;
	        coefficient[0][4] = 3.401378906197122;     coefficient[1][4] = 4.8378133911595285;
	        coefficient[0][5] = 2.1251383942923474;    coefficient[1][5] = 2.5391647730624003;
	    }

        rateLower = T*(coefficient[0][0] + coefficient[0][1]/pow(p/T,4.) - 
		     coefficient[0][2]/pow((p/T),3.) - coefficient[0][3]/pow((p/T),2.) + 
		     coefficient[0][4]/pow((p/T),1.5) - coefficient[0][5]/(p/T));
        rateUpper = T*(coefficient[1][0] + coefficient[1][1]/pow(p/T,4.) - 
		     coefficient[1][2]/pow((p/T),3.) - coefficient[1][3]/pow((p/T),2.) + 
		     coefficient[1][4]/pow((p/T),1.5) - coefficient[1][5]/(p/T));
      
        rate = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
        rate /= 3.; // historic reasons
      
        if (process == 4) rate *= 9./4.;
    }

  return rate;
}

double Elastic::totalRatePos(double p, double T, double alpha_s, int Nf, int process)
{
  // compute the total transition rate in GeV(integrated over k, omega and q and angles) for fixed E=p and temperature T
  // using parametrization of numerically computed integral - then interpolate to get right value for used alpha_s 
  // processes: 1:qq, 2:gq, 3:qg, 4:gg
  // IMPORTANT: all computed values below are for a minimal omega of 0.05*T, so this is the cutoff to use in the calculation also!
  // also Nf=3 was used ... scales out though


  double alpha0 = 0.15;
  double deltaAlpha = 0.03;
  double iAlpha = floor((alpha_s-alpha0)/deltaAlpha+0.001);
  double alphaFrac = (alpha_s-alpha0)/deltaAlpha - iAlpha;
  double coefficient[2][6];
  double rateLower, rateUpper, rate;


  if (process == 1 || process == 2) // qq or gq (scattering with a thermal quark)
    {
      if ( alpha_s >= 0.15 && alpha_s < 0.18 )
	{
	  coefficient[0][0] = 0.12199410313320332;  coefficient[1][0] = 0.15243607717720586;
	  coefficient[0][1] = 0.23732051765097376;  coefficient[1][1] = 0.5403120875137825;
	  coefficient[0][2] = -0.03285419708803458; coefficient[1][2] = 0.06440920730334501;
	  coefficient[0][3] = 0.2255419254079952;   coefficient[1][3] = 0.2881594349535524;
	  coefficient[0][4] = 0.03991522899907729;  coefficient[1][4] = 0.04948438583750772;
	  coefficient[0][5] = 0.05022641428394594;  coefficient[1][5] = 0.07152523367501308;
	}
      else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
	{
	  coefficient[0][0] = 0.15243607717720586;  coefficient[1][0] = 0.15243607717720586;
	  coefficient[0][1] = 0.5403120875137825;   coefficient[1][1] = 0.5403120875137825;
	  coefficient[0][2] = 0.06440920730334501;  coefficient[1][2] = 0.06440920730334501;
	  coefficient[0][3] = 0.2881594349535524;   coefficient[1][3] = 0.2881594349535524;
	  coefficient[0][4] = 0.04948438583750772;  coefficient[1][4] = 0.04948438583750772;
	  coefficient[0][5] = 0.07152523367501308;  coefficient[1][5] = 0.07152523367501308;
	}
      else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
	{
	  coefficient[0][0] = 0.15243607717720586;  coefficient[1][0] = 0.21661000995329158;
	  coefficient[0][1] = 0.5403120875137825;   coefficient[1][1] = 1.4087570376612657;
	  coefficient[0][2] = 0.06440920730334501;  coefficient[1][2] = 0.2713885880193171;
	  coefficient[0][3] = 0.2881594349535524;   coefficient[1][3] = 0.48681971936565244;
	  coefficient[0][4] = 0.04948438583750772;  coefficient[1][4] = 0.09567346780679847;
	  coefficient[0][5] = 0.07152523367501308;  coefficient[1][5] = 0.12780677622585393;
	}
      else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
	{
	  coefficient[0][0] = 0.21661000995329158;  coefficient[1][0] = 0.2501007467879627;
	  coefficient[0][1] = 1.4087570376612657;   coefficient[1][1] = 1.8034683081244214;
	  coefficient[0][2] = 0.2713885880193171;   coefficient[1][2] = 0.228092470920281;
	  coefficient[0][3] = 0.48681971936565244;  coefficient[1][3] = 0.6841577896561725;
	  coefficient[0][4] = 0.09567346780679847;  coefficient[1][4] = 0.15430793601338547;
	  coefficient[0][5] = 0.12780677622585393;  coefficient[1][5] = 0.1648297331159989;
	}
      else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
	{
	  coefficient[0][0] = 0.2501007467879627;   coefficient[1][0] = 0.28440720063047276;
	  coefficient[0][1] = 1.8034683081244214;   coefficient[1][1] = 2.0448244620634055;
	  coefficient[0][2] = 0.228092470920281;    coefficient[1][2] = -0.018574547528236382;
	  coefficient[0][3] = 0.6841577896561725;   coefficient[1][3] = 0.9863974758613413;
	  coefficient[0][4] = 0.15430793601338547;  coefficient[1][4] = 0.2503738253300167;
	  coefficient[0][5] = 0.1648297331159989;   coefficient[1][5] = 0.2090067594645225;
	}
      else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
	{
	  coefficient[0][0] = 0.28440720063047276;  coefficient[1][0] = 0.31945943548344036;
	  coefficient[0][1] = 2.0448244620634055;   coefficient[1][1] = 2.0482495934952256;
	  coefficient[0][2] = -0.018574547528236382;coefficient[1][2] = -0.5350999123662686;
	  coefficient[0][3] = 0.9863974758613413;   coefficient[1][3] = 1.4169725257394696;
	  coefficient[0][4] = 0.2503738253300167;   coefficient[1][4] = 0.3918202096574105;
	  coefficient[0][5] = 0.2090067594645225;   coefficient[1][5] = 0.26103455441873036;
	}
      else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
	{
	  coefficient[0][0] = 0.31945943548344036;  coefficient[1][0] = 0.35519799231686516;
	  coefficient[0][1] = 2.0482495934952256;   coefficient[1][1] = 1.7485135425544152;
	  coefficient[0][2] = -0.5350999123662686;  coefficient[1][2] = -1.3692232011881413;
	  coefficient[0][3] = 1.4169725257394696;   coefficient[1][3] = 1.9906086576701993;
	  coefficient[0][4] = 0.3918202096574105;   coefficient[1][4] = 0.5832315715098879;
	  coefficient[0][5] = 0.26103455441873036;  coefficient[1][5] = 0.32124694953933486;
	}
      else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
	{
	  coefficient[0][0] = 0.35519799231686516;  coefficient[1][0] = 0.39157507493019383;
	  coefficient[0][1] = 1.7485135425544152;   coefficient[1][1] = 1.0778995684787331;
	  coefficient[0][2] = -1.3692232011881413;  coefficient[1][2] = -2.5738838613236457;
	  coefficient[0][3] = 1.9906086576701993;   coefficient[1][3] = 2.727543221296746;
	  coefficient[0][4] = 0.5832315715098879;   coefficient[1][4] = 0.8323699786704292;
	  coefficient[0][5] = 0.32124694953933486;  coefficient[1][5] = 0.3905055907877247;
	}
      else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
	{
	  coefficient[0][0] = 0.39157507493019383;  coefficient[1][0] = 0.4285382777192131;
	  coefficient[0][1] = 1.0778995684787331;   coefficient[1][1] = 0.05505396151716547;
	  coefficient[0][2] = -2.5738838613236457;  coefficient[1][2] = -4.113979132685303;
	  coefficient[0][3] = 2.727543221296746;    coefficient[1][3] = 3.5992808060371506;
	  coefficient[0][4] = 0.8323699786704292;   coefficient[1][4] = 1.1252568207814462;
	  coefficient[0][5] = 0.3905055907877247;   coefficient[1][5] = 0.4667953957378259;
	}

      rateLower = T*(coefficient[0][0] + coefficient[0][1]/pow(p/T,4.) - 
		     coefficient[0][2]/pow((p/T),3.) - coefficient[0][3]/pow((p/T),2.) + 
		     coefficient[0][4]/pow((p/T),1.5) - coefficient[0][5]/(p/T));
      rateUpper = T*(coefficient[1][0] + coefficient[1][1]/pow(p/T,4.) - 
		     coefficient[1][2]/pow((p/T),3.) - coefficient[1][3]/pow((p/T),2.) + 
		     coefficient[1][4]/pow((p/T),1.5) - coefficient[1][5]/(p/T));
      
      rate = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
      
      rate *= Nf/3.; // adjust number of flavors

      if (process == 2) rate *= 9./4.;
    }
  else if (process == 3 || process == 4)  // qg or gg (scattering with a thermal gluon)
    {
      if ( alpha_s >= 0.15 && alpha_s < 0.18 )
	{
	  coefficient[0][0] = 0.6197775378922895;    coefficient[1][0] = 0.7680959463632293;
	  coefficient[0][1] = 1.5268694134079064;    coefficient[1][1] = 3.282164035377037;
	  coefficient[0][2] = 0.6939337312845367;    coefficient[1][2] = 1.6359849897319092;
	  coefficient[0][3] = 0.5967602676773388;    coefficient[1][3] =  0.6770046238563808;
	  coefficient[0][4] = 0.17320784052297564;   coefficient[1][4] = 0.22074166337990309;
	  coefficient[0][5] = 0.28964614117694565;   coefficient[1][5] = 0.4128184793199476;
	}
      else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
	{
	  coefficient[0][0] = 0.7680959463632293;    coefficient[1][0] = 0.9206225398305536;
	  coefficient[0][1] = 3.282164035377037;     coefficient[1][1] = 5.690562370150853;
	  coefficient[0][2] = 1.6359849897319092;    coefficient[1][2] = 2.8341906487774318;
	  coefficient[0][3] = 0.6770046238563808;    coefficient[1][3] = 0.7900156706763937;
	  coefficient[0][4] = 0.22074166337990309;   coefficient[1][4] = 0.2995126102416747;
	  coefficient[0][5] = 0.4128184793199476;    coefficient[1][5] = 0.5598645426609049;
	}
      else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
	{
	  coefficient[0][0] = 0.9206225398305536;    coefficient[1][0] = 1.0767954081327265;
	  coefficient[0][1] = 5.690562370150853;     coefficient[1][1] = 8.378841394880034;
	  coefficient[0][2] = 2.8341906487774318;    coefficient[1][2] = 3.9338968631891396;
	  coefficient[0][3] = 0.7900156706763937;    coefficient[1][3] = 1.0874771229885156;
	  coefficient[0][4] = 0.2995126102416747;    coefficient[1][4] = 0.46570985770548107;
	  coefficient[0][5] = 0.5598645426609049;    coefficient[1][5] = 0.7360069767362173;
	}
      else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
	{
	  coefficient[0][0] = 1.0767954081327265;    coefficient[1][0] = 1.2361819653856791;
	  coefficient[0][1] = 8.378841394880034;     coefficient[1][1] = 10.877148035367144;
	  coefficient[0][2] = 3.9338968631891396;    coefficient[1][2] = 4.526191560392149;
	  coefficient[0][3] = 1.0874771229885156;    coefficient[1][3] = 1.731930015138816;
	  coefficient[0][4] = 0.46570985770548107;   coefficient[1][4] = 0.7769917594310469;
	  coefficient[0][5] = 0.7360069767362173;    coefficient[1][5] = 0.9463662091275489;
	}
      else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
	{
	  coefficient[0][0] = 1.2361819653856791;    coefficient[1][0] = 1.3984393292278847;
	  coefficient[0][1] = 10.877148035367144;    coefficient[1][1] = 12.72181515837248;
	  coefficient[0][2] = 4.526191560392149;     coefficient[1][2] = 4.227297031355039;
	  coefficient[0][3] = 1.731930015138816;     coefficient[1][3] = 2.868526983329731;
	  coefficient[0][4] = 0.7769917594310469;    coefficient[1][4] = 1.2836917844304823;
	  coefficient[0][5] = 0.9463662091275489;    coefficient[1][5] = 1.1953148369630755;
	}
      else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
	{
	  coefficient[0][0] = 1.3984393292278847;    coefficient[1][0] = 1.5632880021613935;
	  coefficient[0][1] = 12.72181515837248;     coefficient[1][1] = 13.502896915302873;
	  coefficient[0][2] = 4.227297031355039;     coefficient[1][2] = 2.7113406243010467;
	  coefficient[0][3] = 2.868526983329731;     coefficient[1][3] = 4.615035662049938;
	  coefficient[0][4] = 1.2836917844304823;    coefficient[1][4] = 2.0259357821768784;
	  coefficient[0][5] = 1.1953148369630755;    coefficient[1][5] = 1.486253368704046;
	}
      else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
	{
	  coefficient[0][0] = 1.5632880021613935;    coefficient[1][0] = 1.730492163581557;
	  coefficient[0][1] = 13.502896915302873;    coefficient[1][1] = 12.913294655478987;
	  coefficient[0][2] = 2.7113406243010467;    coefficient[1][2] = -0.2477159937428581;
	  coefficient[0][3] = 4.615035662049938;     coefficient[1][3] = 7.042004003229154;
	  coefficient[0][4] = 2.0259357821768784;    coefficient[1][4] = 3.0253452576771465;
	  coefficient[0][5] = 1.486253368704046;     coefficient[1][5] = 1.8205651561017433;
	}
      else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
	{
	  coefficient[0][0] = 1.730492163581557;     coefficient[1][0] = 1.8998560359992867;
	  coefficient[0][1] = 12.913294655478987;    coefficient[1][1] = 10.708892844334745;
	  coefficient[0][2] = -0.2477159937428581;   coefficient[1][2] = -4.823210983922782;
	  coefficient[0][3] = 7.042004003229154;     coefficient[1][3] = 10.202109059054063;
	  coefficient[0][4] = 3.0253452576771465;    coefficient[1][4] = 4.298747764427364;
	  coefficient[0][5] = 1.8205651561017433;    coefficient[1][5] = 2.199497022778097;
	}
      else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
	{
	  coefficient[0][0] = 1.8998560359992867;    coefficient[1][0] = 2.071204284004704;
	  coefficient[0][1] = 10.708892844334745;    coefficient[1][1] = 6.741738604119316;
	  coefficient[0][2] = -4.823210983922782;    coefficient[1][2] = -11.099716230158746;
	  coefficient[0][3] = 10.202109059054063;    coefficient[1][3] = 14.106488110189458;
	  coefficient[0][4] = 4.298747764427364;     coefficient[1][4] = 5.846203546614067;
	  coefficient[0][5] = 2.199497022778097;     coefficient[1][5] = 2.62230136903594;
	}

      rateLower = T*(coefficient[0][0] + coefficient[0][1]/pow(p/T,4.) - 
		     coefficient[0][2]/pow((p/T),3.) - coefficient[0][3]/pow((p/T),2.) + 
		     coefficient[0][4]/pow((p/T),1.5) - coefficient[0][5]/(p/T));
      rateUpper = T*(coefficient[1][0] + coefficient[1][1]/pow(p/T,4.) - 
		     coefficient[1][2]/pow((p/T),3.) - coefficient[1][3]/pow((p/T),2.) + 
		     coefficient[1][4]/pow((p/T),1.5) - coefficient[1][5]/(p/T));
      
      rate = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
      rate /= 3.; // historic reasons
      
      if (process == 4) rate *= 9./4.;
    }
  //rate /= hbarc; // to get rate in 1/fm.
  //  cout.precision(8);
  //cout << "T=" << T << ", E=" << p << ", alpha_s=" << alpha_s << ", rate=" << rate << endl; 
  //cout << " lower rate=" << rateLower << endl; 
  //cout << " upper rate=" << rateUpper << endl; 
  return rate;
}

double Elastic::totalRateNeg(double p, double T, double alpha_s, int Nf, int process)
{
  // compute the total transition rate in GeV (integrated over k, omega and q and angles) for fixed E=p and temperature T
  // using parametrization of numerically computed integral - then interpolate to get right value for used alpha_s 
  // processes: 1:qq, 2:gq, 3:qg, 4:gg
  // IMPORTANT: all computed values below are for a minimal omega of 0.05*T, so this is the cutoff to use in the calculation also!
  // also Nf=3 was used ... scales out though


  double alpha0 = 0.15;
  double deltaAlpha = 0.03;
  double iAlpha = floor((alpha_s-alpha0)/deltaAlpha+0.001);
  double alphaFrac = (alpha_s-alpha0)/deltaAlpha - iAlpha;
  double coefficient[2][6];
  double rateLower, rateUpper, rate;
  if (process == 1 || process == 2) // qq or gq (scattering with a thermal quark)
    {
      if ( alpha_s >= 0.15 && alpha_s < 0.18 )
	{
	  coefficient[0][0] = 0.059730780828164666;  coefficient[1][0] = 0.07216040121873951;
	  coefficient[0][1] = 0.3631534872548789;    coefficient[1][1] = 0.5471138972952214;
	  coefficient[0][2] = 0.39845046966687;      coefficient[1][2] = 0.5792306465939813;
	  coefficient[0][3] = 0.11946615972422633;   coefficient[1][3] = 0.1723078888161528;
	  coefficient[0][4] = 0.03559276204445307;   coefficient[1][4] = 0.05120408756812135;
	  coefficient[0][5] = 0.00291041965645416;   coefficient[1][5] = 0.0041777787108426695;
	}
      else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
	{
	  coefficient[0][0] = 0.07216040121873951;   coefficient[1][0] = 0.0846236909779996;
	  coefficient[0][1] = 0.5471138972952214;    coefficient[1][1] = 0.7725791286875564;
	  coefficient[0][2] = 0.5792306465939813;    coefficient[1][2] = 0.7931123494736929;
	  coefficient[0][3] = 0.1723078888161528;    coefficient[1][3] = 0.23406373724706608;
	  coefficient[0][4] = 0.05120408756812135;   coefficient[1][4] = 0.06935459958589639;
	  coefficient[0][5] = 0.0041777787108426695; coefficient[1][5] = 0.005644055718614478;
	}
      else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
	{
	  coefficient[0][0] = 0.0846236909779996;    coefficient[1][0] = 0.09711346786308672;
	  coefficient[0][1] = 0.7725791286875564;    coefficient[1][1] = 1.0370070423372528;
	  coefficient[0][2] = 0.7931123494736929;    coefficient[1][2] = 1.036935526583188;
	  coefficient[0][3] = 0.23406373724706608;   coefficient[1][3] = 0.3034025403259155;
	  coefficient[0][4] = 0.06935459958589639;   coefficient[1][4] = 0.08957509599955729;
	  coefficient[0][5] = 0.005644055718614478;  coefficient[1][5] = 0.007264393465593115;
	}
      else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
	{
	  coefficient[0][0] = 0.09711346786308672;   coefficient[1][0] = 0.10962479860948156;
	  coefficient[0][1] = 1.0370070423372528;    coefficient[1][1] = 1.3372010137066646;
	  coefficient[0][2] = 1.036935526583188;     coefficient[1][2] = 1.307456863105879;
	  coefficient[0][3] = 0.3034025403259155;    coefficient[1][3] = 0.37910328734850873;
	  coefficient[0][4] = 0.08957509599955729;   coefficient[1][4] = 0.111456899829735;
	  coefficient[0][5] = 0.007264393465593115;  coefficient[1][5] = 0.009000895144744121;
	}
      else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
	{
	  coefficient[0][0] = 0.10962479860948156;   coefficient[1][0] = 0.1221541053951596;
	  coefficient[0][1] = 1.3372010137066646;    coefficient[1][1] = 1.6686065099273535;
	  coefficient[0][2] = 1.307456863105879;     coefficient[1][2] = 1.600404353394210;
	  coefficient[0][3] = 0.37910328734850873;   coefficient[1][3] = 0.4594932213772782;
	  coefficient[0][4] = 0.111456899829735;     coefficient[1][4] = 0.13442407314203592;
	  coefficient[0][5] = 0.009000895144744121;  coefficient[1][5] = 0.010800449048880756;
	}
      else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
	{
	  coefficient[0][0] = 0.1221541053951596;    coefficient[1][0] = 0.13469861652518803;
	  coefficient[0][1] = 1.6686065099273535;    coefficient[1][1] = 2.0276317271182074;
	  coefficient[0][2] = 1.600404353394210;     coefficient[1][2] = 1.912613330851788;
	  coefficient[0][3] = 0.4594932213772782;    coefficient[1][3] = 0.5434449889160747;
	  coefficient[0][4] = 0.13442407314203592;   coefficient[1][4] = 0.15810564016236883;
	  coefficient[0][5] = 0.010800449048880756;  coefficient[1][5] = 0.012629305933671075;
	}
      else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
	{
	  coefficient[0][0] = 0.13469861652518803;   coefficient[1][0] = 0.14725614907227047;
	  coefficient[0][1] = 2.0276317271182074;    coefficient[1][1] = 2.4109122726272654;
	  coefficient[0][2] = 1.912613330851788;     coefficient[1][2] = 2.241198157777867;
	  coefficient[0][3] = 0.5434449889160747;    coefficient[1][3] = 0.6299396046048817;
	  coefficient[0][4] = 0.15810564016236883;   coefficient[1][4] = 0.18216575652552597;
	  coefficient[0][5] = 0.012629305933671075;  coefficient[1][5] = 0.014456750325370632;
	}
      else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
	{
	  coefficient[0][0] = 0.14725614907227047;   coefficient[1][0] = 0.15982489441001274;
	  coefficient[0][1] = 2.4109122726272654;    coefficient[1][1] = 2.815254291049982;
	  coefficient[0][2] = 2.241198157777867;     coefficient[1][2] = 2.583462624103292;
	  coefficient[0][3] = 0.6299396046048817;    coefficient[1][3] = 0.7180274724508857;
	  coefficient[0][4] = 0.18216575652552597;   coefficient[1][4] = 0.20629432847931367;
	  coefficient[0][5] = 0.014456750325370632;  coefficient[1][5] = 0.01625568033747704;
	}
      else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
	{
	  coefficient[0][0] = 0.15982489441001274;   coefficient[1][0] = 0.17240331582158486;
	  coefficient[0][1] = 2.815254291049982;     coefficient[1][1] = 3.238290376079149;
	  coefficient[0][2] = 2.583462624103292;     coefficient[1][2] = 2.9374985881586273;
	  coefficient[0][3] = 0.7180274724508857;    coefficient[1][3] = 0.8071008047950518;
	  coefficient[0][4] = 0.20629432847931367;   coefficient[1][4] = 0.23030341585944009;
	  coefficient[0][5] = 0.01625568033747704;   coefficient[1][5] = 0.018010096397556033;
	}

      rateLower = T*(coefficient[0][0] + coefficient[0][1]/pow(p/T,4.) - 
		     coefficient[0][2]/pow((p/T),3.) + coefficient[0][3]/pow((p/T),2.) - 
		     coefficient[0][4]/pow((p/T),1.5) + coefficient[0][5]/(p/T));
      rateUpper = T*(coefficient[1][0] + coefficient[1][1]/pow(p/T,4.) - 
		     coefficient[1][2]/pow((p/T),3.) + coefficient[1][3]/pow((p/T),2.) - 
		     coefficient[1][4]/pow((p/T),1.5) + coefficient[1][5]/(p/T));
      
      rate = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
      
      rate *= Nf/3.; // adjust number of flavors

      if (process == 2) rate *= 9./4.;
    }
  else if (process == 3 || process == 4) // qg or gg (scattering with a thermal gluon)
    {
      if ( alpha_s >= 0.15 && alpha_s < 0.18 )
	{
	  coefficient[0][0] = 0.3166913701414167;    coefficient[1][0] = 0.3804527486448292;
	  coefficient[0][1] = 1.0992070651449564;    coefficient[1][1] = 1.7114836115114735;
	  coefficient[0][2] = 1.4232219292986843;    coefficient[1][2] = 2.093540209692791;
	  coefficient[0][3] = 0.4655268754156213;    coefficient[1][3] = 0.6787666526042345;
	  coefficient[0][4] = 0.1444497238817506;    coefficient[1][4] = 0.21014340589291994;
	  coefficient[0][5] = 0.012281442189758532;  coefficient[1][5] = 0.017832857383112792;
	}
      else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
	{
	  coefficient[0][0] = 0.3804527486448292;    coefficient[1][0] = 0.44393432393104637;
	  coefficient[0][1] = 1.7114836115114735;    coefficient[1][1] = 2.483663499207573;
	  coefficient[0][2] = 2.093540209692791;     coefficient[1][2] = 2.8979112438999044;
	  coefficient[0][3] = 0.6787666526042345;    coefficient[1][3] = 0.9316968286688833;
	  coefficient[0][4] = 0.21014340589291994;   coefficient[1][4] = 0.28780901378857465;
	  coefficient[0][5] = 0.017832857383112792;  coefficient[1][5] = 0.02438874287373154;
	}
      else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
	{
	  coefficient[0][0] = 0.44393432393104637;   coefficient[1][0] = 0.5071424487228405;
	  coefficient[0][1] = 2.483663499207573;     coefficient[1][1] = 3.4070556051784515;
	  coefficient[0][2] = 2.8979112438999044;    coefficient[1][2] = 3.824491419496227;
	  coefficient[0][3] = 0.9316968286688833;    coefficient[1][3] = 1.2191109771387096;
	  coefficient[0][4] = 0.28780901378857465;   coefficient[1][4] = 0.3755459972857442;
	  coefficient[0][5] = 0.02438874287373154;   coefficient[1][5] = 0.03174924882247299;
	}
      else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
	{
	  coefficient[0][0] = 0.5071424487228405;    coefficient[1][0] = 0.5700856365203443;
	  coefficient[0][1] = 3.4070556051784515;    coefficient[1][1] = 4.466964606692036;
	  coefficient[0][2] = 3.824491419496227;     coefficient[1][2] = 4.857999356928031;
	  coefficient[0][3] = 1.2191109771387096;    coefficient[1][3] = 1.5348360053714125;
	  coefficient[0][4] = 0.3755459972857442;    coefficient[1][4] = 0.471215528026891;
	  coefficient[0][5] = 0.03174924882247299;   coefficient[1][5] = 0.03971601962636114;
	}
      else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
	{
	  coefficient[0][0] = 0.5700856365203443;    coefficient[1][0] = 0.6327732610959403;
	  coefficient[0][1] = 4.466964606692036;     coefficient[1][1] = 5.646624908846933;
	  coefficient[0][2] = 4.857999356928031;     coefficient[1][2] = 5.982691423451806;
	  coefficient[0][3] = 1.5348360053714125;    coefficient[1][3] = 1.8728243844356356;
	  coefficient[0][4] = 0.471215528026891;     coefficient[1][4] = 0.572761497659723;
	  coefficient[0][5] = 0.03971601962636114;   coefficient[1][5] = 0.04809998538877525;
	}
      else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
	{
	  coefficient[0][0] = 0.6327732610959403;    coefficient[1][0] = 0.6952147319486842;
	  coefficient[0][1] = 5.646624908846933;     coefficient[1][1] = 6.931552369487635;
	  coefficient[0][2] = 5.982691423451806;     coefficient[1][2] = 7.185588273540373;
	  coefficient[0][3] = 1.8728243844356356;    coefficient[1][3] = 2.228328283532209;
	  coefficient[0][4] = 0.572761497659723;     coefficient[1][4] = 0.6786029643259804;
	  coefficient[0][5] = 0.04809998538877525;   coefficient[1][5] = 0.056755908207122875;
	}
      else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
	{
	  coefficient[0][0] = 0.6952147319486842;    coefficient[1][0] = 0.7574189285141091;
	  coefficient[0][1] = 6.931552369487635;    coefficient[1][1] = 8.307255807497631;
	  coefficient[0][2] = 7.185588273540373;    coefficient[1][2] = 8.454112812202247;
	  coefficient[0][3] = 2.228328283532209;     coefficient[1][3] = 2.596781386863294;
	  coefficient[0][4] = 0.6786029643259804;    coefficient[1][4] = 0.7872276571283385;
	  coefficient[0][5] = 0.056755908207122875;     coefficient[1][5] = 0.06554867983133447;
	}
      else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
	{
	  coefficient[0][0] = 0.7574189285141091;     coefficient[1][0] = 0.8193940883937045;
	  coefficient[0][1] = 8.307255807497631;    coefficient[1][1] = 9.761691032241623;
	  coefficient[0][2] = 8.454112812202247;    coefficient[1][2] = 9.777948193339808;
	  coefficient[0][3] = 2.596781386863294;     coefficient[1][3] = 2.9744411293541457;
	  coefficient[0][4] = 0.7872276571283385;    coefficient[1][4] = 0.8973688582323887;
	  coefficient[0][5] = 0.06554867983133447;    coefficient[1][5] = 0.07435862848596686;
	}
      else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
	{
	  coefficient[0][0] = 0.8193940883937045;    coefficient[1][0] = 0.8811479514201789;
	  coefficient[0][1] = 9.761691032241623;    coefficient[1][1] = 11.286034194965852;
	  coefficient[0][2] = 9.777948193339808;     coefficient[1][2] = 11.15001447311135;
	  coefficient[0][3] = 2.9744411293541457;    coefficient[1][3] = 3.3591358778545803;
	  coefficient[0][4] = 0.8973688582323887;     coefficient[1][4] = 1.0083901554550654;
	  coefficient[0][5] = 0.07435862848596686;     coefficient[1][5] = 0.08313659597360733;
	}

      rateLower = T*(coefficient[0][0] + coefficient[0][1]/pow(p/T,4.) - 
		     coefficient[0][2]/pow((p/T),3.) + coefficient[0][3]/pow((p/T),2.) - 
		     coefficient[0][4]/pow((p/T),1.5) + coefficient[0][5]/(p/T));
      rateUpper = T*(coefficient[1][0] + coefficient[1][1]/pow(p/T,4.) - 
		     coefficient[1][2]/pow((p/T),3.) + coefficient[1][3]/pow((p/T),2.) - 
		     coefficient[1][4]/pow((p/T),1.5) + coefficient[1][5]/(p/T));
      
      rate = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
      rate /= 3.; // historic reasons
      
      if (process == 4) rate *= 9./4.;
    }
  //rate /= hbarc; // to get rate in 1/fm.
  //  cout.precision(8);
  //cout << "T=" << T << ", E=" << p << ", alpha_s=" << alpha_s << ", rate=" << rate << endl; 
  //cout << " lower rate=" << rateLower << endl; 
  //cout << " upper rate=" << rateUpper << endl; 
  return rate;
}

// gets the rate from the table using interpolation, done in Import
double Elastic::function(double p, double omega, double alpha_s, Import *import, int process)
{
  return import->getElasticRate(p, omega, alpha_s, process);
}

// calculates the area under the envelope function when using the rejection method
// (integrals had been solved analytically before)
double Elastic::area (double y, double alpha_s, int posNegSwitch, int process)
{
    if (process==1 || process==2)
    {
        if (posNegSwitch==1)
	        return(0.0333333*alpha_s*(7.+ 4.*alpha_s)*(2.99573+log(y)));
        else 
	        return(-0.133333*alpha_s*(1.75+alpha_s)*log(-0.0833333*y));
    }
    else if (process==3 || process==4)
    {
        if (posNegSwitch==1)
	        return(0.05*alpha_s*(7.+ 4.*alpha_s)*(2.99573+log(y)));
        else 
	        return(-0.2*alpha_s*(1.75+alpha_s)*log(-0.0833333*y));
    }
    else
        return 0.0;
}

double Elastic::findValuesRejection(double u, double T, double alpha_s, int Nf,
		Random *random, Import *import, int process)
{
    double f;
    double Pos, Neg;
    double x;
    double y;
    Pos = totalRatePos(u,T,alpha_s, Nf, process);
    Neg = totalRateNeg(u,T,alpha_s, Nf, process);
    
    cout.precision(6);
    // this switch will hold the decision whether k is positive or negative: 0=negative, 1=positive
    int posNegSwitch = 1; 

    if (process == 1 || process == 2) // for qq or gq
    {
        // decide whether omega shall be positive or negative (if x (uniform on [0,1]) < area(k<0)/area(all k) then k<0)
        if (random->genrand64_real1()<Neg/(Neg+Pos)) posNegSwitch = 0;
        // generate random number according to envelope distribution f(y)
        if( posNegSwitch == 1 ) // if omega > 0
        {
            do
            {
                y = exp((-1.41428*pow(10.,9.)*alpha_s - 8.08158*pow(10.,8.)*alpha_s*alpha_s +
      	                2.02327*pow(10.,9.)*random->genrand64_real1()*area(u,alpha_s,posNegSwitch,process))
      	                /(alpha_s*(4.72097*pow(10.,8.) + 2.6977*pow(10.,8.)*alpha_s)));
                // here random->genrand64_real1()*area(u.,posNegSwitch) 
                // is a uniform random number on [0, area under f(x)]
                x = random->genrand64_real1();
                // x is uniform on [0,1]
                //i++;
                if (function(u,y,alpha_s,import,process)>(((alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y))))))
                    cout << "[Elastic::findValuesRejection]:WARNING, envelope too low!" << endl;
            } while( x > (function(u,y,alpha_s,import,process))/((alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y))))); 
            // reject if x is larger than the ratio p(y)/f(y), f(y)=(alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y)))
            f=y;
        }
        else // if omega < 0
        {
            do
            {
                y = -12. * exp(-(30.* random->genrand64_real1()*area(-0.05,alpha_s,posNegSwitch,process))/(alpha_s*(7.+ 4.* alpha_s)));
                // here random->genrand64_real1()*area(u+12.,posNegSwitch) 
                // is a uniform random number on [0, area under f(x)]
                x = random->genrand64_real1();
                // x is uniform on [0,1]
                //i++;
                //if (i>200)
                // {
                // cout << "(omega<0,process 1 or 2) with p=" << u << ", we" << " are having " << i << " rejections." << endl;
                // }
                if (function(u,y,alpha_s,import,process)>(((alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y))))))
                    cout << "[Elastic::findValuesRejection]:WARNING, envelope too low!" << endl;
             } while( x > (function(u,y,alpha_s,import,process))/((alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y))))); 
	        // reject if x is larger than the ratio p(y)/f(y)
	        f=y;
        }
    }
    else if (process == 3 || process == 4) // for qg or gg
    {
        // decide whether omega shall be positive or negative (if x (uniform on [0,1]) < area(k<0)/area(all k) then k<0)
        if (random->genrand64_real1()<Neg/(Neg+Pos)) posNegSwitch = 0;
        // generate random number according to envelope distribution f(y)
        if( posNegSwitch == 1 ) // if omega > 0
	    {
	        do
	        {
	            y = exp((-2.32591*pow(10.,17.)*alpha_s - 1.32909*pow(10.,17.)*alpha_s*alpha_s +
		             2.2183*pow(10.,17.)*random->genrand64_real1()*area(u,alpha_s,posNegSwitch,process))
		            /(alpha_s*(7.76406*pow(10.,16.) + 4.43661*pow(10.,16.)*alpha_s)));
	            // here random->genrand64_real1()*area(u,posNegSwitch) 
	            // is a uniform random number on [0, area under f(x)]
	            x = random->genrand64_real1();
	            // x is uniform on [0,1]
	            if (function(u,y,alpha_s,import,process)>(1.5*(alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y)))))
	                cout << "WARNING, envelope too low, for omega>0!" << endl;
	        } while( x > (function(u,y,alpha_s,import,process))/(1.5*(alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y))))); 
	        // reject if x is larger than the ratio p(y)/f(y), f(y)=1.5*(alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y)))
	        f=y;
	    }
        else // if omega < 0
	    {
            do
            {
	            y = -12.*exp(-(2.81475*pow(10.,15.)*random->genrand64_real1()*area(-0.05,alpha_s,posNegSwitch,process))
		             /(alpha_s*(9.85162*pow(10.,14.) + 5.6295*pow(10.,14.)*alpha_s)));
	            // here random->genrand64_real1()*area(u,posNegSwitch) 
	            // is a uniform random number on [0, area under f(x)]
	            x = random->genrand64_real1();
	            // x is uniform on [0,1]
	            if (function(u,y,alpha_s,import,process)>(1.5*(alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y)))))
		            cout << "WARNING, envelope too low, for omega<0!" << endl;
            } while( x > (function(u,y,alpha_s,import,process))/(1.5*(alpha_s/0.15)*((0.035 + alpha_s*0.02)/(sqrt(y*y))))); 
	        // reject if x is larger than the ratio p(y)/f(y) 
	        f=y;
	    }
    }
    return f;
}

double Elastic::functionOmegaQ(double p, double omega, double q, double alpha_s, Import *import, int process)
{
  return import->getElasticRateOmegaQ(p, omega, q, alpha_s, process);
}

// calculates the area under the envelope function when using the rejection method
// (integrals had been solved analytically before)
double Elastic::areaOmegaQ (double y, double omega, double alpha_s, int process)
{
  double A, B;
  if (process==1 || process==2)
    {
      A = (0.7+alpha_s-0.3)*0.0014*(1000. + 40./sqrt(omega*omega) + 10.5* pow(omega,4.))*alpha_s;
      B = 2.*sqrt(omega*omega) + 0.01;
    }
  else if (process==3 || process==4)
    {
      A = (0.7+alpha_s-0.3)*0.0022*(1000. + 40./sqrt(omega*omega) + 10.5* pow(omega,4.))*alpha_s;
      B = 2.*sqrt(omega*omega) + 0.002;
    }      
  return (0.5*A*(atan(y*y/sqrt(B))-atan(omega*omega/sqrt(B))))/sqrt(B);
}

double Elastic::areaOmegaQ2 (double y, double omega, double alpha_s, int process)
{
  double A;
  if (process==1 || process==2)
    {
      A =  0.016*exp(0.5*omega)*alpha_s;
      if (omega > 22.) A = 0.0014*exp(0.5*omega)*alpha_s;
    }
  else if (process==3 || process==4)
    {
      A =  0.024*exp(0.5*omega)*alpha_s;
      if (omega > 22.) A = 0.0021*exp(0.5*omega)*alpha_s;
    }      
  return A*(2.*exp(-0.5*omega)-2.*exp(-0.5*y));
}

// use rejection method for all -3 < omega < 4, and Metropolis for the rarer omegas outside that range 
double Elastic::findValuesRejectionOmegaQ(double p, double omega, double T, double alpha_s, int Nf,
					  Random *random, Import *import, int process)
{
    double x, y;
    double A, B;
    cout.precision(6);
    int i = 0;

    if ( omega < 10. && omega >-3 )
    {
        //cout << "doing small omega=" << omega << endl;
        if (process == 1 || process == 2) // for qq or gq
        {
            A = (0.7+alpha_s)*0.0014*(1000. + 40./sqrt(omega*omega) + 10.5* pow(omega,4.))*alpha_s;
            B = 2.*sqrt(omega*omega) + 0.01;
        }
        else if (process == 3 || process == 4) // for qg or gg
        {
            A = (0.7+alpha_s)*0.0022*(1000. + 40./sqrt(omega*omega) + 10.5* pow(omega,4.))*alpha_s;
            B = 2.*sqrt(omega*omega) + 0.002;
        }
        // generate random number according to envelope distribution f(y)
        do
	    {
	        y = pow(B,0.25)*sqrt(tan((2.*sqrt(B)*random->genrand64_real1()*areaOmegaQ(p,omega,alpha_s,process)
			        	    +A*atan(omega*omega/sqrt(B)))/A));
	        //cout << "y=" << y << endl;
	        // is a uniform random number on [0, area under f(x)]
	        x = random->genrand64_real1();
	        // x is uniform on [0,1]
	        if (functionOmegaQ(p,omega,y,alpha_s,import,process)>(A*y/(pow(y,4.)+B)))
	            cout << "WARNING, envelope too low determining q, for omega=" << omega << " and process=" << process << endl;
	        i++;
        } while( x > (functionOmegaQ(p,omega,y,alpha_s,import,process))/(A*y/(pow(y,4.)+B)) || y > p ); 
        // reject if x is larger than the ratio p(y)/f(y), f(y)=(A*y/(pow(y,4.)+B))
        if (i>2000)
	    {
            cout << defaultfloat << " with omega=" << omega << ", we" << " had " << i << " rejections." << endl;
	    }
    }
    //else if (omega < 10. && omega > -3)
    //{
    //    //cout << "doing large omega=" << omega << endl;
    //    if (process == 1 || process == 2) // for qq or gq
    //    {
    //        A = 0.016*exp(0.5*omega)*alpha_s;
    //        if (omega > 10.) A = 0.01*exp(0.5*omega)*alpha_s;
    //        if (omega > 22.) A = 0.0014*exp(0.5*omega)*alpha_s;
    //    }
    //    else if (process == 3 || process == 4) // for qg or gg
    //    {
    //        A = 0.024*exp(0.5*omega)*alpha_s;
    //        if (omega > 10.) A = 0.016*exp(0.5*omega)*alpha_s;
    //        if (omega > 22.) A = 0.0021*exp(0.5*omega)*alpha_s;
    //    }
    //    // generate random number according to envelope distribution f(y)
    //    do
    //    {
    //        y = -2.*log(-(0.5*(-2.*A*exp(-0.5*omega) + random->genrand64_real1()*areaOmegaQ2(p,omega,alpha_s,process)))/A);
    //        //cout << "y=" << y << endl;
    //        // is a uniform random number on [0, area under f(x)]
    //        x = random->genrand64_real1();
    //        // x is uniform on [0,1]
    //        if (functionOmegaQ(p,omega,y,alpha_s,import,process)>(A*exp(-y*0.5)))
    //        {
    //            cout << "WARNING, omega > 6, envelope too low for omega=" << omega << ", and q=" << y << endl;
    //            cout << functionOmegaQ(p,omega,y,alpha_s,import,process) << ", env=" << (A*exp(-y*0.5)) << endl;
    //        }
    //        i++;
    //        if (i>10000)
    //        {
    //            cout << "stuck, i=" << i << ", omega=" << omega << endl;
    //            cout << functionOmegaQ(p,omega,y,alpha_s,import,process) << ", env=" << (A*exp(-y*0.5)) << endl;
    //            cout << "y=" << y << endl << endl;
    //        }
    //    } while( x > (functionOmegaQ(p,omega,y,alpha_s,import,process))/(A*exp(-y*0.5)) || y > p ); 
    //    // reject if x is larger than the ratio p(y)/f(y), f(y)=(A*exp(-y*0.5))
    //    if (i>2000)
    //    {
    //        cout << defaultfloat << " with omega=" << omega << ", we" << " had " << i << " rejections." << endl;
    //    }
    //}
    else // for large omega
    {
        //cout << "doing large omega Metropolis, omega=" << omega << endl;
        cout.precision(12);
        double g = 0, g_new = 0;
        double ratio;
        double y_new;
        //int rejections=0;
        //int acceptances=0;
        //int rejections1 = 0;
        // the ranges in which the variables u and phi need to be sampled
        const double y_min = sqrt(omega*omega);
        const double y_max = p;
        // randomly select initial values of q=y, such that
        do
	    {
	        y = y_min + random->genrand64_real1()*(y_max - y_min);
	        g = functionOmegaQ(p,omega,y,alpha_s,import,process);
	        //cout << "omega=" << omega << ", y=" << y << ", g=" << g << endl;
            //rejections1 +=1 ;
	    } while(g < 1E-35); // RY: g==0. before. now just check to make sure it's "slightly" above zero
        //cout<< "had "<<rejections1<<" rejections in first loop"<<endl; 
        // number of steps in the Markov chain
        const int n_steps = 500;
        
        for ( int i=0; i < n_steps; i++ )
        {        
            do
            {
                y_new = y_min + random->genrand64_real1()*(y_max - y_min); 

            }
            while( y_new < y_min || y_new > y_max); // check that the new value is in range
        
            g_new = functionOmegaQ(p,omega,y_new,alpha_s,import,process); // calculate the function at the proposed point
            ratio = g_new / g;                                            // ratio of g(y_new) / g(y)
            if ( ratio > 1.0 || random->genrand64_real1() < ratio )         // accept if g(y_new) > g(y) 
            // or with probability g(y_new) / g(y)
            {
                y = y_new;
                g = g_new;
                //acceptances+=1;
            }
            //else rejections+=1;
	    }
        //cout << "Metropolis rejections = " << rejections << endl;
        //cout << "Metropolis acceptances = " << acceptances << endl;
    }
    return y;
}

Vec4 Elastic::getNewMomentum(Vec4 vecpRest, double omega, double q, Random *random)
{
    // everything in GeV - make sure the input (vecpRest, omega, q) is in GeV already ( no more temperature )
    Vec4 newvecp, newvecpTemp;
    Vec4 vecet;
    Vec4 vecqt;
    Vec4 vecql;
    double qt;
    double ql;
    double cosTheta_pq;
    double pRest=vecpRest.pAbs();
    double phi;
    double M[3][3]; //rotation matrix
    double u;
    Vec4 r;
  
    //if (omega == q)
    if ( fabs(omega - q) < 1E-14)
    {
        newvecp=vecpRest*(pRest-omega)/pRest;
        return newvecp;
    }
  
    cosTheta_pq = (-omega*omega+2.*pRest*omega+q*q)/(2.*pRest*q);
    if (q<sqrt(omega*omega))
        cout << "WARNING: q<omega=" << q << "<" << omega << endl;
    qt = q*sqrt(1.-cosTheta_pq*cosTheta_pq);                     // transverse momentum transfer
    ql = q*cosTheta_pq;
    /*
    if(cosTheta_pq>-1.)
      {
        cout << "Theta_pq=" << acos(cosTheta_pq)/PI*180. << endl;
        cout << "omega=" << omega << endl;
        cout << "q=" << q << endl;
        cout << endl;
      }
    */
    //cout << "qt=" << qt << endl;
    //cout << "ql=" << ql << endl << " " << endl;
    
    //cout << "sqrt(ql*ql+qt*qt)=" << sqrt(ql*ql+qt*qt) << endl;
    //cout << "vecpRest=" << vecpRest;
    //cout << "|vecpRest|=" << vecpRest.pAbs() << endl;
    
    if ( vecpRest.py()*vecpRest.py() > vecpRest.px()*vecpRest.px() )
    {
        vecet.px(0.);
        vecet.py(-vecpRest.pz());
        vecet.pz(vecpRest.py());
    }
    else
    {
        vecet.px(vecpRest.pz());
        vecet.py(0.);
        vecet.pz(-vecpRest.px());
    }
  
    vecet/=vecet.pAbs();                            // normalize \vec{et} to one
    //cout << "|vecet|=" << vecet.pAbs() << endl;
    
    vecqt = qt*vecet;                               // this is the transverse transferred momentum vector
    //cout << "vecqt=" << vecqt;
  
    vecql = vecpRest/pRest*ql;
    //cout << "vecql=" << vecql << endl;
  
    newvecpTemp = vecpRest - vecqt;                 // change transverse momentum
    
    newvecpTemp-=vecql;                             // change longitudinal momentum
  
    //newvecpTemp*=(pRest - omega*T)/vecpRest.pAbs();     // new quark momentum
      
    //cout << "new newvecpTemp=" << newvecpTemp;
  
    //cout << "|newvecpTemp|=" << newvecpTemp.pAbs() << endl;
    //cout << "p-omega=" << pRest-omega << endl;
    // now rotate the new vector about the old one by a random angle phi.
  
    phi=2.*PI*random->genrand64_real1();
    r = vecpRest/vecpRest.pAbs();
    //cout << "r=" << r;
    u = 1.-cos(phi);
    //cout << "phi=" << phi << endl;
   
    // define the rotation matrix for rotations around pvecRest
    
    M[0][0]=r.px()*r.px()*u+cos(phi);
    M[1][0]=r.px()*r.py()*u-r.pz()*sin(phi);
    M[2][0]=r.px()*r.pz()*u+r.py()*sin(phi);
  
    M[0][1]=r.py()*r.px()*u+r.pz()*sin(phi);
    M[1][1]=r.py()*r.py()*u+cos(phi);
    M[2][1]=r.py()*r.pz()*u-r.px()*sin(phi);
  
    M[0][2]=r.pz()*r.px()*u-r.py()*sin(phi);
    M[1][2]=r.pz()*r.py()*u+r.px()*sin(phi);
    M[2][2]=r.pz()*r.pz()*u+cos(phi);
  
    newvecp.px(M[0][0]*newvecpTemp.px()+M[0][1]*newvecpTemp.py()+M[0][2]*newvecpTemp.pz());
    newvecp.py(M[1][0]*newvecpTemp.px()+M[1][1]*newvecpTemp.py()+M[1][2]*newvecpTemp.pz());
    newvecp.pz(M[2][0]*newvecpTemp.px()+M[2][1]*newvecpTemp.py()+M[2][2]*newvecpTemp.pz());
    newvecp.e(sqrt(newvecp.px()*newvecp.px()+newvecp.py()*newvecp.py()+newvecp.pz()*newvecp.pz()));
  
    //newvecp.px(newvecpTemp.px());
    //newvecp.py(newvecpTemp.py());
    //newvecp.pz(newvecpTemp.pz());
    
    //cout << "|newvecpTemp|=" << newvecpTemp.pAbs() << endl;
    //cout << "|newvecp|=" << newvecp.pAbs() << endl << " " << endl;
    //cout << "newvecp=" << newvecp << endl;
    return newvecp;
}

Vec4 Elastic::getRecoilMomentum(Vec4 vecq, double temp, int kind, Random *random)
{
  double q, omega;
  double k, k_min;
  double cosTh, sinTh;
  Vec4 u1, u2;    // rotation axis
  double t1, t2;  // 1-cos(angle)
  double M1[3][3], M2[3][3]; //rotation matrix
  double phi;

  Vec4 veckTemp, veck;
  Vec4 vecRecoil;


  q = vecq.pAbs();
  omega = vecq.e();

  // minimum momenum of thermal parton k that makes recoil parton on-shell
  if (fabs(omega) > q) omega *= q/fabs(omega);  // momentum transfer is always space-like
  k_min = (q-omega)/2.;

  k = random->thermal2(k_min, temp, kind);
  cosTh = (2*k*omega - q*q + omega*omega)/(2*k*q);

  // choose unit vector perpendicular to vecq with which vecq is rotated by theta
  u1.px(vecq.py());
  u1.py(-vecq.px());
  u1.pz(0.);
  if(u1.pAbs() > 0.)
    u1 /= u1.pAbs();
  else
    u1.px(1.);
  u1.e(u1.pAbs());
  //cout << "u1 = " << u1 << " u1*vecq = " << u1.px()*vecq.px() + u1.py()*vecq.py() + u1.pz()*vecq.pz() << endl;;

  t1 = 1.-cosTh;
  sinTh = sqrt(1.-cosTh*cosTh);
 
  // define the rotation matrix for theta rotations 
  M1[0][0]=u1.px()*u1.px()*t1+cosTh;
  M1[1][0]=u1.px()*u1.py()*t1-u1.pz()*sinTh;
  M1[2][0]=u1.px()*u1.pz()*t1+u1.py()*sinTh;

  M1[0][1]=u1.py()*u1.px()*t1+u1.pz()*sinTh;
  M1[1][1]=u1.py()*u1.py()*t1+cosTh;
  M1[2][1]=u1.py()*u1.pz()*t1-u1.px()*sinTh;

  M1[0][2]=u1.pz()*u1.px()*t1-u1.py()*sinTh;
  M1[1][2]=u1.pz()*u1.py()*t1+u1.px()*sinTh;
  M1[2][2]=u1.pz()*u1.pz()*t1+cosTh;

  //cout << "cosTh = " << cosTh << " sinTh = " << sinTh << " 1-cosTh = " << t1 << endl;
  //cout << "M1[i][0]  = " << M1[0][0] << " " << M1[1][0] << " " << M1[2][0] << endl;
  //cout << "M1[i][1]  = " << M1[0][1] << " " << M1[1][1] << " " << M1[2][1] << endl;
  //cout << "M1[i][2]  = " << M1[0][2] << " " << M1[1][2] << " " << M1[2][2] << endl;

  // get the momentum rotated by theta
  veckTemp.px(M1[0][0]*vecq.px()+M1[0][1]*vecq.py()+M1[0][2]*vecq.pz());
  veckTemp.py(M1[1][0]*vecq.px()+M1[1][1]*vecq.py()+M1[1][2]*vecq.pz());
  veckTemp.pz(M1[2][0]*vecq.px()+M1[2][1]*vecq.py()+M1[2][2]*vecq.pz());
  veckTemp *= k/veckTemp.pAbs();
  veckTemp.e(veckTemp.pAbs());

  //cout << "[getRecoilMomentum]::veckTemp = " << veckTemp;
  // do this again for rnadom azimuthal angle.  
  phi=2.*M_PI*random->genrand64_real1();
  u2 = vecq/vecq.pAbs();
  t2 = 1.-cos(phi);
  
  M2[0][0]=u2.px()*u2.px()*t2+cos(phi);
  M2[1][0]=u2.px()*u2.py()*t2-u2.pz()*sin(phi);
  M2[2][0]=u2.px()*u2.pz()*t2+u2.py()*sin(phi);

  M2[0][1]=u2.py()*u2.px()*t2+u2.pz()*sin(phi);
  M2[1][1]=u2.py()*u2.py()*t2+cos(phi);
  M2[2][1]=u2.py()*u2.pz()*t2-u2.px()*sin(phi);

  M2[0][2]=u2.pz()*u2.px()*t2-u2.py()*sin(phi);
  M2[1][2]=u2.pz()*u2.py()*t2+u2.px()*sin(phi);
  M2[2][2]=u2.pz()*u2.pz()*t2+cos(phi);

  veck.px(M2[0][0]*veckTemp.px()+M2[0][1]*veckTemp.py()+M2[0][2]*veckTemp.pz());
  veck.py(M2[1][0]*veckTemp.px()+M2[1][1]*veckTemp.py()+M2[1][2]*veckTemp.pz());
  veck.pz(M2[2][0]*veckTemp.px()+M2[2][1]*veckTemp.py()+M2[2][2]*veckTemp.pz());
  veck.e(veck.pAbs());

  //cout << "[getRecoilMomentum]::vecq = " << vecq;
  //cout << "[getRecoilMomentum]::veck = " << veck;
  vecRecoil = veck + vecq;

  return vecRecoil;
}

Vec4 Elastic::getThermalMomentum(Vec4 vecq, double temp, int kind, Random *random)
{
    double q, omega;
    double k, k_min;
    double cosTh, sinTh;
    Vec4 u1, u2;    // rotation axis
    double t1, t2;  // 1-cos(angle)
    double M1[3][3], M2[3][3]; //rotation matrix
    double phi;
  
    Vec4 veckTemp, veck;
  
    q = vecq.pAbs();
    omega = vecq.e();
  
    // minimum momenum of thermal parton k that makes recoil parton on-shell
    if ( fabs(omega) > q)
        omega *= q/fabs(omega);  // momentum transfer is always space-like
    k_min = (q-omega)/2.;
  
    k = random->thermal2(k_min, temp, kind);
    cosTh = (2*k*omega - q*q + omega*omega)/(2*k*q);
  
    // choose unit vector perpendicular to vecq with which vecq is rotated by theta
    u1.px(vecq.py());
    u1.py(-vecq.px());
    u1.pz(0.);
    if(u1.pAbs() > 0.)
        u1 /= u1.pAbs();
    else
        u1.px(1.);
    u1.e(u1.pAbs());
    //cout << "u1 = " << u1 << " u1*vecq = " << u1.px()*vecq.px() + u1.py()*vecq.py() + u1.pz()*vecq.pz() << endl;;
  
    t1 = 1.-cosTh;
    sinTh = sqrt(1.-cosTh*cosTh);
   
    // define the rotation matrix for theta rotations 
    M1[0][0]=u1.px()*u1.px()*t1+cosTh;
    M1[1][0]=u1.px()*u1.py()*t1-u1.pz()*sinTh;
    M1[2][0]=u1.px()*u1.pz()*t1+u1.py()*sinTh;
  
    M1[0][1]=u1.py()*u1.px()*t1+u1.pz()*sinTh;
    M1[1][1]=u1.py()*u1.py()*t1+cosTh;
    M1[2][1]=u1.py()*u1.pz()*t1-u1.px()*sinTh;
  
    M1[0][2]=u1.pz()*u1.px()*t1-u1.py()*sinTh;
    M1[1][2]=u1.pz()*u1.py()*t1+u1.px()*sinTh;
    M1[2][2]=u1.pz()*u1.pz()*t1+cosTh;
  
    //cout << "cosTh = " << cosTh << " sinTh = " << sinTh << " 1-cosTh = " << t1 << endl;
    //cout << "M1[i][0]  = " << M1[0][0] << " " << M1[1][0] << " " << M1[2][0] << endl;
    //cout << "M1[i][1]  = " << M1[0][1] << " " << M1[1][1] << " " << M1[2][1] << endl;
    //cout << "M1[i][2]  = " << M1[0][2] << " " << M1[1][2] << " " << M1[2][2] << endl;
  
    // get the momentum rotated by theta
    veckTemp.px(M1[0][0]*vecq.px()+M1[0][1]*vecq.py()+M1[0][2]*vecq.pz());
    veckTemp.py(M1[1][0]*vecq.px()+M1[1][1]*vecq.py()+M1[1][2]*vecq.pz());
    veckTemp.pz(M1[2][0]*vecq.px()+M1[2][1]*vecq.py()+M1[2][2]*vecq.pz());
    veckTemp *= k/veckTemp.pAbs();
    veckTemp.e(veckTemp.pAbs());
  
    //cout << "[getRecoilMomentum]::veckTemp = " << veckTemp;
    // do this again for rnadom azimuthal angle.  
    phi=2.*M_PI*random->genrand64_real1();
    u2 = vecq/vecq.pAbs();
    t2 = 1.-cos(phi);
    
    M2[0][0]=u2.px()*u2.px()*t2+cos(phi);
    M2[1][0]=u2.px()*u2.py()*t2-u2.pz()*sin(phi);
    M2[2][0]=u2.px()*u2.pz()*t2+u2.py()*sin(phi);
  
    M2[0][1]=u2.py()*u2.px()*t2+u2.pz()*sin(phi);
    M2[1][1]=u2.py()*u2.py()*t2+cos(phi);
    M2[2][1]=u2.py()*u2.pz()*t2-u2.px()*sin(phi);
  
    M2[0][2]=u2.pz()*u2.px()*t2-u2.py()*sin(phi);
    M2[1][2]=u2.pz()*u2.py()*t2+u2.px()*sin(phi);
    M2[2][2]=u2.pz()*u2.pz()*t2+cos(phi);
  
    veck.px(M2[0][0]*veckTemp.px()+M2[0][1]*veckTemp.py()+M2[0][2]*veckTemp.pz());
    veck.py(M2[1][0]*veckTemp.px()+M2[1][1]*veckTemp.py()+M2[1][2]*veckTemp.pz());
    veck.pz(M2[2][0]*veckTemp.px()+M2[2][1]*veckTemp.py()+M2[2][2]*veckTemp.pz());
    veck.e(veck.pAbs());
  
    //cout << "[getRecoilMomentum]::vecq = " << vecq;
    //cout << "[getRecoilMomentum]::veck = " << veck;
  
    return veck;
}
