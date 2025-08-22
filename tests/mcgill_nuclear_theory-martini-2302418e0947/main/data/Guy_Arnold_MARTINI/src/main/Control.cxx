#include "Control.h"
#include "Constants.h"

using namespace std;

Control::Control(const FileNames* fileList, const Constants* constants, double idt, double idE, double iT, const double ialphas, int iNf, const double iinitialE, const double iwidth, const double imaxTime, const double imaxEnergy, const int istrate, const int iBetheHeitler, const int iBDMPS, const int idoRad, const int idoCol, const int inewRadGamma, const int imperfermi, const int isize, const int irank)
{
  
  dt = idt;
  dE = idE;
  T = iT;
  alphas = ialphas;
  Nf = iNf;
  initialE = iinitialE;
  width = iwidth;
  maxTime = imaxTime;
  maxEnergy = imaxEnergy;
  strate = istrate;
  BetheHeitler = iBetheHeitler;
  BDMPS = iBDMPS;
  doRad = idoRad;
  doCol = idoCol;
  newRadGamma = inewRadGamma;
  mperfermi = imperfermi;
  size=isize;
  rank=irank;

  LogEmax=constants->LogEmax;
  LogEmin=constants->LogEmin;
  LogStepE=constants->LogStepE;
  LogOmegaMin=constants->LogOmegaMin;
  LogStepOmega=constants->LogStepOmega;
  Tfile=constants->Tfile;
  stepT=constants->stepT;


  controlsFileNameList=fileList;
  
  int esize = static_cast<int>(maxEnergy/dE+1+0.0001);
  cout << "esize=" << esize << endl;
  Pq = new double [(esize+1)];
  Pg = new double [(esize+1)];
  Pem = new double [(esize+1)];
  PqA = new double [(esize+1)];
  PgA = new double [(esize+1)];
 
  if(strate==1)
    {
      string path="../rates/";
      string file[4];
      ifstream fin[4];
      
      //open files with transition rates to read in:
      for(int i=0; i<=3; i++)
	{
	  file[i] = path+controlsFileNameList->GetFN(i+5).c_str();
	  fin[i].open(file[i].c_str(),ios::in);     
	  if(!fin[i])
	    {
	      cerr << "[Control::Control]: ERROR: Unable to open file " << file[i] << endl;
	      exit(1);
	    }
	}
        
      double temp[4];
      int Esize=static_cast<int>((LogEmax-LogEmin)/LogStepE+1);
      //int Tsize=static_cast<int>((Tmax-Tmin)/stepT+1);
      int omegaSize=static_cast<int>((2*(LogEmax-LogOmegaMin))/LogStepOmega+3);
      vector<double> trqqvec;
      vector<double> trqgvec;
      vector<double> trgqvec;
      vector<double> trggvec;

      trqq = new double[(Esize)*(omegaSize)];
      trqg = new double[(Esize)*(omegaSize)];
      trgq = new double[(Esize)*(omegaSize)];
      trgg = new double[(Esize)*(omegaSize)];
      
      int count=1;
      
      for(int j=0;j<4;j++)
	{
	  while ( !fin[j].eof() )
	    {

	      fin[j] >> temp[j];
	      if(!(count%4)) 
		{
		  if(j==0) trqqvec.push_back( temp[0] ) ;
		  if(j==1) trqgvec.push_back( temp[1] ) ;
		  if(j==2) trgqvec.push_back( temp[2] ) ;
		  if(j==3) trggvec.push_back( temp[3] ) ;
		}
	      count++;
	    }
	}

      count=0;


      //ofstream fout(controlsFileNameList->GetFN(9).c_str(),ios::out);     
      for(int iE=0;iE<Esize;iE++)
	  for(int iOmega=0;iOmega<omegaSize;iOmega++)
	    {
	      trqq[iE*((omegaSize))+iOmega]=trqqvec[count];
	      trqg[iE*((omegaSize))+iOmega]=trqgvec[count];
	      trgq[iE*((omegaSize))+iOmega]=trgqvec[count];
	      trgg[iE*((omegaSize))+iOmega]=trggvec[count];
	      count++;
	      /*cout << pow(10,LogEmin+iE*LogStepE) << " " 
		   << Tmin+iT*stepT << " " 
		   << " "
		   << trgq[iE*(Tsize*omegaSize)+iT*omegaSize+iOmega] << endl;*/ // test ouput
	    }
      //fout.close();
	     
      //test interpolation:
      
      double En=3;
      double T=0.4;
      int Tsize=1;
      double omega=-0.2;
      int iE=floor((log(En)-LogEmin)/LogStepE);
      int iT=0;
      if(iE<0) iE=0;
      
      int iOmega;
      if(omega>0) iOmega=floor((log(omega)-LogOmegaMin+(LogEmax-LogOmegaMin))/LogStepOmega+1);
      if(omega<0) iOmega=floor((-log(-omega)+LogOmegaMin+(LogEmax-LogOmegaMin))/LogStepOmega);
   
      
      double fracO;
      if(omega>0) 
	fracO=(omega-(exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin))))
	  /(exp(((iOmega)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin))
	    -exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin)));
      
      if(omega<0)
	fracO=(omega-(-exp((-(iOmega)*LogStepOmega+LogEmax))))
	  /((exp((-(iOmega)*LogStepOmega+LogEmax)))-exp((-(iOmega+1)*LogStepOmega+LogEmax)));
      
      double fracE=(En-exp((LogEmin+iE*LogStepE)))/(exp((LogEmin+(iE+1)*LogStepE))-exp((LogEmin+iE*LogStepE)));
      
      double trval = ((1-fracE)*((1-fracO)*trqq[iE*(omegaSize)+iOmega]+fracO*trqq[iE*(omegaSize)+(iOmega+1)])
		      +fracE*((1-fracO)*trqq[(iE+1)*(omegaSize)+iOmega]+fracO*trqq[(iE+1)*(omegaSize)+(iOmega+1)]));
      if(omega>0)
	{
	  cout << "omega=" << omega << ", iOmega=" << iOmega 
	       << ", value at iOmega=" << exp(((iOmega-1)*LogStepOmega-(LogEmax-LogOmegaMin)+LogOmegaMin)) 
	       << ", trval=" << trval << ", fracE=" << fracE 
	       << ", fracO=" << fracO << ", value at iOmega: " 
	       << trqq[iE*(omegaSize)+iOmega]
	       << ", value at iOmega+1: " <<  trqq[(iE)*(omegaSize)+(iOmega+1)] 
	       << ", value at iOmega: " 
	       << trqq[(iE+1)*(omegaSize)+iOmega]
	       << ", value at iOmega+1: " <<  trqq[(iE+1)*(omegaSize)+(iOmega+1)] 
	       << ", E at iE=" << exp((LogEmin+iE*LogStepE))
	       << ",LogStepE=" << LogStepE << endl;
	}
      else
	{
	  cout << "omega=" << omega << ", iOmega=" << iOmega 
	       << ", value at iOmega=" << -exp((-(iOmega)*LogStepOmega+LogEmax)) 
	       << ", trval=" << trval << ", fracE=" << fracE 
	       << ", fracO=" << fracO << ", value at iOmega: " 
	       << trqq[iE*(omegaSize)+iOmega] 
	       << ", value at iOmega+1: " <<  trqq[iE*(omegaSize)+(iOmega+1)] << endl;
	       }      
    }//end if(strate==1)


  if(rank==0)
    {
      ofstream fout(controlsFileNameList->GetFN(1).c_str(),ios::app);     
      fout << endl;
      fout << "[Control::Control]:" << endl;
      fout << "                maximal time = " << maxTime << " fm" << endl;
      fout.close();
    }

  timeEvolution = new Evolution(controlsFileNameList, constants, dt, dE, T, alphas, Nf, initialE, width, maxEnergy, BetheHeitler, BDMPS, doRad, doCol, newRadGamma);
  
  timeEvolution->dEdtOutput();

  timeEvolution->initP(Pq, Pg, Pem);
}

void Control::evolve()
{
  int it;
  int mt = static_cast<int>(maxTime/hbarc/dt+0.0001); // compute number of steps
  int esize = maxEnergy/dE+1+0.0001;
  ofstream fout2(controlsFileNameList->GetFN(4).c_str(),ios::out);     
  ofstream fout(controlsFileNameList->GetFN(3).c_str(),ios::out);     
  fout << "E Pq Pg" << endl;
  
  for(it=0;it<=mt;it++)
    {
      if((mperfermi*it*dt*hbarc-floor(mperfermi*it*dt*hbarc))<mperfermi*dt*hbarc-0.00000001)
	{
	  cout << "writing table at time " << it*dt*hbarc << " fm." << endl;
	  for(int i=1;i<esize;i++)
	    {
	      fout << it*dt*hbarc << " " << i*dE << " " << Pq[i] << " " << Pg[i] << endl;
	    }
	  fout << " " << endl;
	}

      if(it%10==0) cout << "time = " << it*dt*hbarc << " fm" << endl;
      if(strate==1)
	{
	  timeEvolution->runsmooth(Pq,Pg,Pem,PqA,PgA,trqq,trqg,trgq,trgg);
	}
      else 
	{
	  timeEvolution->run(Pq,Pg,Pem,PqA,PgA);
	}

      
      double sum=0;
      double mean=0;
      double mean2=0;
      
      for(int i=50;i<esize;i++)
	{
	  sum+=(Pq[i]+Pg[i])*dE;
	  mean+=i*dE*(Pq[i])*dE;
	  mean2+=i*dE*(Pq[i]+Pg[i])*dE;
	}
      cout << (it+1)*dt*hbarc << " " << sum << endl;
      fout2 << (it+1)*dt*hbarc << " " << mean << " " << mean2 << endl;

    }
  
  fout.close();
  fout2.close();

}
