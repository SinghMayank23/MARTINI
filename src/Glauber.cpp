// Glauber.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines for sampling the initial geometry according to the Glauber model

#include "Glauber.h"

void Glauber::lookup_nucleus_info(Nucleus *nucleus, string nucleus_name)
{
    char *temp_name = new char[50];
    strcpy(temp_name, nucleus_name.c_str());
    nucleus->name = temp_name;

    nucleus->rho_WS = 0.15;   //default rho. this WILL change to the right value that gives integral(rho) = A
    if(nucleus_name.compare("p") == 0)
    {
        nucleus->A = 1;
        nucleus->Z = 1;
        nucleus->w_WS = 0.0;
        nucleus->R_WS = 1.0;
        nucleus->a_WS = 1.0;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else if(nucleus_name.compare("C") == 0)
    {
        nucleus->A = 12;
        nucleus->Z = 6;
        nucleus->w_WS = 1.403;
        nucleus->R_WS = 2.44;
        nucleus->a_WS = 1.635;
        nucleus->AnumFunc = 1; //Anum2HO;
        nucleus->AnumFuncIntegrand = 1; //Anum2HOInt;
        nucleus->DensityFunc = 1; //NuInt2HO;
    }
    else if(nucleus_name.compare("O") == 0)
    {
        nucleus->A = 16;
        nucleus->Z = 8;
        nucleus->R_WS = 2.608;
        nucleus->w_WS = -0.051;
        nucleus->a_WS = 0.513;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else if(nucleus_name.compare("Al") == 0)
    {
        nucleus->A = 27;
        nucleus->Z = 13;
        nucleus->R_WS = 3.07;
        nucleus->w_WS = 0.0;
        nucleus->a_WS = 0.519;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else if(nucleus_name.compare("Cu") == 0)
    {
        nucleus->A = 63;
        nucleus->Z = 29;
        nucleus->R_WS = 4.163;
        nucleus->w_WS = 0.0;
        nucleus->a_WS = 0.606;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else if(nucleus_name.compare("Au") == 0)
    {
        nucleus->A = 197;
        nucleus->Z = 79;
        nucleus->w_WS = 0.0;
        nucleus->R_WS = 6.38;
        nucleus->a_WS = 0.505;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else if(nucleus_name.compare("Pb") == 0)
    {
        nucleus->A = 208;
        nucleus->Z = 82;
        nucleus->w_WS = 0.0;
        nucleus->R_WS = 6.62;
        nucleus->a_WS = 0.546;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else if(nucleus_name.compare("U") == 0)
    {
        nucleus->A = 238;
        nucleus->Z = 92;
        nucleus->w_WS = 0.0;
        nucleus->R_WS = 6.874;
        nucleus->a_WS = 0.556;
        nucleus->AnumFunc = 3; //Anum3Fermi;
        nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
        nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
    else
    {
        fprintf(stderr, "Unknown_nucleus: %s \n", nucleus_name.c_str());
        fprintf(stderr, "Exiting...\n");
        exit(0);
    }
}

void Glauber::FindNucleusData(Nucleus *nucleus, char *target, char *file_name)
{
 char *s, *name, c, *func_name, *tmp_name;
 static int ind;
 string tmp_str;
 FILE *input;
 FILE *tmp_file;
 int bytes_read;

 s = char_malloc(120);
 name = char_malloc(120);

 input = fopen(file_name, "r");
 tmp_str ="tmp.dat";

 tmp_name = (char*)malloc((tmp_str.length()+10) * sizeof(char));
 memset((void*)tmp_name,0,sizeof(tmp_name));
 strcpy(tmp_name,tmp_str.c_str());

 tmp_file = fopen(tmp_name, "w");

 bytes_read=fscanf(input, "%s", s);

 ind = 0;
 while(strcmp(s, "EndOfData") != 0)
  {
   bytes_read=fscanf(input, "%s", name);
   if(strcmp(name, target) == 0)
    {
     ind++;
     fprintf(tmp_file, "Name %s\n", name);
     c = getc(input);
     while(c != 'N')
      {
       fprintf(tmp_file, "%c", c);
       c = getc(input);
      }/* while */
     break;
    }/* if target is found */
   bytes_read=fscanf(input, "%s", s);
  }/* while end of the file is not encountered */

 fprintf(tmp_file, "%s", "EndOfData");
 
 fclose(input);
 fclose(tmp_file);

 nucleus->rho_WS = 0.15; /* default rho.  this WILL change to the right 
			 value that gives integral(rho) = A */
 
 nucleus->name = StringFind(tmp_name, "Name");
 nucleus->A = DFind(tmp_name, "A");
 nucleus->Z = DFind(tmp_name, "Z");
 if(ind != 0)
  {
   nucleus->w_WS = DFind(tmp_name, "w_WS");
   nucleus->a_WS = DFind(tmp_name, "a_WS");
   nucleus->R_WS = DFind(tmp_name, "R_WS");
   func_name = StringFind(tmp_name, "density_func");
  }
 else
  {
   nucleus->w_WS = 0.0; 
   nucleus->a_WS = 0.53;
   nucleus->R_WS = 1.15*pow(nucleus->A, 1.0/3.0);
   func_name = "3Fermi";
  }
 
 if(strcmp(func_name, "2HO")==0) 
  {
    nucleus->AnumFunc = 1; //Anum2HO;
    nucleus->AnumFuncIntegrand = 1; //Anum2HOInt;
    nucleus->DensityFunc = 1; //NuInt2HO;
  }
 else if(strcmp(func_name, "3Gauss")==0) 
  {
    nucleus->AnumFunc = 2; //Anum3Gauss;
    nucleus->AnumFuncIntegrand = 2; //Anum3GaussInt;
    nucleus->DensityFunc = 2; //NuInt3Gauss;
  }
 else if(strcmp(func_name, "3Fermi")==0) 
  {
    nucleus->AnumFunc = 3; //Anum3Fermi;
    nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
    nucleus->DensityFunc = 3; //NuInt3Fermi;
  }
 char_free(s);
 char_free(name);
 remove(tmp_name);

}/* FindNucleusData */


void Glauber::PrintLexusData()
{
 fprintf(stdout, "LexusData.SigmaNN = %e\n",  LexusData.SigmaNN);
 fprintf(stdout, "LexusData.InterMax = %d\n", LexusData.InterMax);
}/* PrintLexusData */


void Glauber::PrintNucleusData(Nucleus *nucleus)
{
 fprintf(stdout, "Nucleus Name: %s\n", nucleus->name);
 
 //fprintf(stderr, "Nucleus.A = %e\n", nucleus->A);
 //fprintf(stderr, "Nucleus.Z = %e\n", nucleus->Z);
 //fprintf(stderr, "Nucleus.w_WS = %e\n", nucleus->w_WS);
 //fprintf(stderr, "Nucleus.a_WS = %e\n", nucleus->a_WS);
 //fprintf(stderr, "Nucleus.R_WS = %e\n", nucleus->R_WS);

}

void Glauber::PrintCollisionSystem(Data *nucleus)
{
  fprintf(stdout, "Projectile Name: %s\n", nucleus->Projectile.name);
  fprintf(stdout, "Target Name: %s\n", nucleus->Target.name);
}


int Glauber::LinearFindXorg(double x, double *Vx, int ymax)
{
/* finds the first of the 4 points, x is between the second and the third */
 
 int x_org;
 double nx;

 nx = ymax*(x - Vx[0])/(Vx[ymax] - Vx[0]);
 
 x_org = (int) nx;
 x_org -= 1;

 if( x_org <= 0 ) return 0;
 else if(x_org >= ymax - 3) return ymax - 3;
 else return x_org;

}/* Linear Find Xorg */


double Glauber::FourPtInterpolate(double x, double *Vx, double *Vy, double h, int x_org, int ymax)
{
 /* interpolating points are x_org, x_org+1, x_org+2, x_org+3 */
 /* cubic polynomial approximation */

 double a, b, c, d, f;

 MakeCoeff(&a, &b, &c, &d,  Vy, Vx, h, x_org);
 
 f = a*pow(x - Vx[x_org], 3.);
 f += b*pow(x - Vx[x_org], 2.);
 f += c*(x - Vx[x_org]);
 f += d;
 
 return f;
}/* FourPtInterpolate */


void Glauber::MakeCoeff(double *a, double *b, double *c, double *d, 
			double *Vy, double *Vx, double h, int x_org)
{
 double f0, f1, f2, f3;

 f0 = Vy[x_org];
 f1 = Vy[x_org+1];
 f2 = Vy[x_org+2];
 f3 = Vy[x_org+3];

 *a =  (-f0 + 3.0*f1 - 3.0*f2 + f3)/(6.0*h*h*h);
 
 *b =  (2.0*f0 - 5.0*f1 + 4.0*f2 - f3)/(2.0*h*h);

 *c =  (-11.0*f0 + 18.0*f1 - 9.0*f2 + 2.0*f3)/(6.0*h);

 *d = f0;

}/* MakeCoeff */

int Glauber::FindXorg(double x, double *Vx, int ymax)
{
 int i, x_org;

 i = 0;
 while(Vx[i] < x) i++;
 
 x_org = i - 2;
 
 if( x_org <= 1 ) return 1;
 else if(x_org >= ymax - 3) return ymax - 3;
 else return x_org;

}/* Find Xorg */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::VInterpolate(double x, double *Vx, double *Vy, int ymax)
{
 int x_org;
 double h;

 if( (x < Vx[0])||(x > Vx[ymax]) )
  {
   fprintf(stderr, 
           "VInterpolate: x = %le is outside the range (%le, %le).\n", 
	    x,Vx[0],Vx[ymax]);
   fprintf(stderr, "This can't happen.  Exiting...\n");
   exit(0);
  }

/* we only deal with evenly spaced Vx */
/* x_org is the first of the 4 points */

 x_org = LinearFindXorg(x, Vx, ymax);

 h = (Vx[ymax] - Vx[0])/ymax;

 return FourPtInterpolate(x, Vx, Vy, h, x_org, ymax);

}/* VInterpolate */

double *Glauber::MakeVx(double down, double up, int maxi_num)
{
 static double dx, *vx;
 int i;

 vx = vector_malloc(maxi_num + 1);
 dx = (up - down)/maxi_num;
 
 for(i=0; i<=maxi_num; i++)
  {
   vx[i] = dx*i;
  }

 return vx;

}/* MakeVx */


double *Glauber::MakeVy(char *s, double *vx, int maxi_num)
{
 int i, di;
 static double *vy;
 static char *st, *dst;
 FILE *data_file;

 if(maxi_num > 200) di = 100;
 if(maxi_num <= 200) di = 20;
 
 vy = vector_malloc(maxi_num + 1);

 st = char_malloc(120);
 dst = char_malloc(120);
 strcpy(dst, ".dat");

 st = strcpy(st, s);
 st = strcat(st, dst);

 data_file = fopen(st, "a");
 fprintf(data_file, "EndOfData\n");
 
 for(i=0; i<=maxi_num; i++)
  {
   vy[i] = NuInS(vx[i]);
   if(i % di == 0)
    {
     fprintf(stderr, "%s[%d] = %le\n", s, i, vy[i]); 
    }
   fprintf(data_file, "%20.16e  %20.16e\n", vx[i], vy[i]);
  }

 fclose(data_file);
 char_free(st);
 char_free(dst);
 
 return vy;
}/* MakeVy */


double *Glauber::ReadInVx(char *file_name, int maxi_num, int quiet)
{
 static double x, *vx;
 int i;
 FILE *input;
 static char *s, *sx;
 int bytes_read;
 s = char_malloc(120);
 sx = char_malloc(120);
 
 vx = vector_malloc(maxi_num + 1);

 if(quiet == 1)
  {
   fprintf(stdout, "Reading in Vx from %s ...\n", file_name);
   }
 
 input = fopen(file_name, "r");
 bytes_read=fscanf(input, "%s", s);
 while(strcmp(s, "EndOfData") != 0)
  {
   bytes_read=fscanf(input, "%s", sx);
   bytes_read=fscanf(input, "%s", s);
  }

 for(i=0; i<=maxi_num; i++)
  {
   bytes_read=fscanf(input, "%lf", &x);
   vx[i] = x;
   bytes_read=fscanf(input, "%lf", &x);
  }
 fclose(input);

 char_free(sx);
 char_free(s);
 return vx;

}/* ReadInVx */


double *Glauber::ReadInVy(char *file_name, int maxi_num, int quiet)
{
 static double y, *vy;
 int i;
 FILE *input;
 static char *s, *sy;
 int bytes_read;
 s = char_malloc(120);
 sy = char_malloc(120);
 
 vy = vector_malloc(maxi_num + 1);

 if(quiet == 1)
 {
  fprintf(stdout, "Reading in Vy from %s ...\n", file_name);
 }
 
 input = fopen(file_name, "r");
 bytes_read=fscanf(input, "%s", s);
 while(strcmp(s, "EndOfData") != 0)
  {
   bytes_read=fscanf(input, "%s", sy);
   bytes_read=fscanf(input, "%s", s);
  }

 for(i=0; i<=maxi_num; i++)
  {
   bytes_read=fscanf(input, "%lf", &y);
   bytes_read=fscanf(input, "%lf", &y);
   vy[i] = y;
  }
 fclose(input);
 
 char_free(s);
 char_free(sy);
 return vy;

}/* ReadInVy */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::InterNuPInSP(double s)
{
 double y;
 static int ind = 0;
 static double up, down; 
 static int maxi_num; 
 static double *vx, *vy;
 FILE *output;
 static char *st, *dst;
 ind++;

 if(LexusData.Projectile.A == 1.0) return 0.0;

 CalcRho(&(LexusData.Projectile));

 up = 2.0*LexusData.SCutOff;
 down = 0.0; 
 maxi_num = LexusData.InterMax; 

 string path = "";
 const char* MARTINIPATH = "MARTINIPATH";
 char* envPath = getenv(MARTINIPATH);
 
 if (envPath != 0 && *envPath != '\0') 
   {
     int i = 0;
     while (*(envPath+i) != '\0') path += *(envPath+(i++));
   }
 else path = ".";
 string pathAndFile = path+"/main/data/NuPInSP";
 
 char *paf;
 paf = (char*)malloc((pathAndFile.length()+10) * sizeof(char));
 memset((void*)paf,0,sizeof(paf));
 strcpy(paf,pathAndFile.c_str());

 st = char_malloc(120);
 dst = char_malloc(120);
 strcpy(dst, ".dat");

 st = strcpy(st, paf);
 st = strcat(st, dst);

 if(ind == 1)
  {
    if(IsFile(st))
      {
	vx = ReadInVx(st, maxi_num, 1);
	vy = ReadInVy(st, maxi_num, 1);
      }
    else
      {
	output = fopen(st, "w");
	fprintf(output, "Name %s\n", LexusData.Projectile.name);
	fprintf(output, "SigmaNN %e\n", LexusData.SigmaNN);
	fclose(output);
	
	vx = MakeVx(down, up, maxi_num);
	vy = MakeVy(paf, vx, maxi_num);
      }
  }/* if ind */
 
 free(paf);
 char_free(st);
 char_free(dst);

 if(s > up) return 0.0;
 else{
   y = VInterpolate(s, vx, vy, maxi_num);
   if( y < 0.0 ) return 0.0; 
   else return y;
  }
}/* InterNuPInSP */

double Glauber::InterNuTInST(double s)
{
 double y;
 static int ind = 0;
 static double up, down; 
 static int maxi_num; 
 static double *vx, *vy;
 static char *st, *dst;
 FILE *output;

 ind++;
 if(LexusData.Target.A == 1.0) return 0.0;

 string path = "";
 const char* MARTINIPATH = "MARTINIPATH";
 char* envPath = getenv(MARTINIPATH);
 
 if (envPath != 0 && *envPath != '\0') 
   {
     int i = 0;
     while (*(envPath+i) != '\0') path += *(envPath+(i++));
   }
 else path = ".";
 string pathAndFile = path+"/main/data/NuTInST";
 
 char *paf;
 paf = (char*)malloc((pathAndFile.length()+10) * sizeof(char));
 memset((void*)paf,0,sizeof(paf));
 strcpy(paf,pathAndFile.c_str());

 st = char_malloc(120);
 dst = char_malloc(120);
 strcpy(dst, ".dat");


 st = strcpy(st, paf);
 st = strcat(st, dst);
 
 if(ind == 1) 
  {
   CalcRho(&(LexusData.Target));
   
   up = 2.0*LexusData.SCutOff;
   down = 0.0; 
   maxi_num = LexusData.InterMax; 
     
   if(IsFile(st))
    {
     vx = ReadInVx(st, maxi_num, 1);
     vy = ReadInVy(st, maxi_num, 1);
    }
   else
    {
     output = fopen(st, "w");
     fprintf(output, "Name %s\n", LexusData.Target.name);
     fprintf(output, "SigmaNN %e\n", LexusData.SigmaNN);
     fclose(output);
  
     vx = MakeVx(down, up, maxi_num);
     vy = MakeVy(paf, vx, maxi_num);
    }
  }/* if ind */
 
 free(paf);
 char_free(st);
 char_free(dst);

 if(s > up) return 0.0;
 else{
   y = VInterpolate(s, vx, vy, maxi_num);
   if( y < 0.0 ) return 0.0; 
   else return y;
  }
}/* InterNuTInST */

void Glauber::CalcRho(Nucleus *nucleus)
{
/* to pass to AnumIntegrand */ 
 
 Nuc_WS = nucleus;
 
 double R_WS = nucleus->R_WS;   
 double f = 1.0;
 
 if (nucleus->AnumFunc==1)
   f = Anum2HO(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==2)
   f = Anum3Gauss(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==3)
   f = Anum3Fermi(R_WS)/(nucleus->rho_WS);
 nucleus->rho_WS = (nucleus->A)/f;
}/* CalcRho */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::NuInS(double s)
{
 double y;
 int count;
 int id;

/* to pass to the DensityFunc's */
 NuInS_S = s;

 id = Nuc_WS->DensityFunc;

 count = 0;
 y = integral(id, 0.0, 1.0, TOL, &count); 
 
 return y;
}

double Glauber::Anum3Fermi(double R_WS)
{
 int count=0;
 double up, down, a_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 rho = Nuc_WS->rho_WS;

/* to pass to Anumintegrand */
 AnumR = R_WS/a_WS;
 
 down = 0.0;
 up = 1.0;
 
 f = integral(4, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum3Fermi */


double Glauber::Anum3FermiInt(double xi) 
{
 double f;
 double r;
 double R_WS, w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0-tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 R_WS = AnumR;
 r = -log(xi);

 f = r*r;
 f *= ( 1.0 + w_WS*pow(r/R_WS, 2.) );
 f /= (xi + exp(-R_WS));
  
 return f;
}/* Anum3FermiInt */


double Glauber::NuInt3Fermi(double xi)
{
 double f;
 double c;
 double z, r, s;
 double w_WS, R_WS, a_WS, rho;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 
   
/* devide by a_WS, make life simpler */
   
   s = NuInS_S/a_WS;
   z = -log(xi);
   r = sqrt(s*s + z*z);
   R_WS /= a_WS;

   c = exp(-R_WS);
   
   f = 2.0*a_WS*rho*(LexusData.SigmaNN);
   f *= 1.0 + w_WS*pow(r/R_WS, 2.);
   f /= xi + c*exp(s*s/(r + z)); 
 
 return f;
}/* NuInt3Fermi */


/* %%%%%%% 3 Parameter Gauss %%%%%%%%%%%% */


double Glauber::Anum3Gauss(double R_WS)
{
 int count=0;
 double up, down, a_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 rho = Nuc_WS->rho_WS;

/* to pass to Anumintegrand */
 AnumR = R_WS/a_WS;
 
 down = 0.0;
 up = 1.0;
 
 f = integral(5, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum3Gauss */


double Glauber::Anum3GaussInt(double xi) 
{
 double y;
 double r_sqr;
 double R_WS, w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0 - tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 
 R_WS = AnumR;
 
 r_sqr = -log(xi);

 y = sqrt(r_sqr);
 
 y *= 1.0 + w_WS*r_sqr/pow(R_WS, 2.);

/* 2 comes from dr^2 = 2 rdr */

 y /= 2.0*(xi + exp(-R_WS*R_WS));
 
 return y;
}/* Anum3GaussInt */


double Glauber::NuInt3Gauss(double xi)
{
 double f;
 double c;
 double z_sqr, r_sqr, s;
 double w_WS, R_WS, a_WS, rho;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 

/* devide by a_WS, make life simpler */

   s = NuInS_S/a_WS;
   z_sqr = -log(xi);
   r_sqr = s*s + z_sqr;
   R_WS /= a_WS;

   c = exp(-R_WS*R_WS);

   f = a_WS*rho*(LexusData.SigmaNN);
   f *= 1.0 + w_WS*r_sqr/pow(R_WS, 2.);
   f /= sqrt(z_sqr)*(xi + c*exp(s*s)); 
 
 return f;
}/* NuInt3Gauss */



/* %%%%%%% 2 Parameter HO %%%%%%%%%%%% */


double Glauber::Anum2HO(double R_WS)
{
 int count=0;
 double up, down, a_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 rho = Nuc_WS->rho_WS;

 down = 0.0;
 up = 1.0;
 
 f = integral(6, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum2HO */


double Glauber::Anum2HOInt(double xi) 
{
 double y;
 double r_sqr, r;
 double w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0 - tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 
 r_sqr = -log(xi);
 
 r = sqrt(r_sqr);
 
/* 2 comes from dr^2 = 2 rdr */
 
 y = r + w_WS*r*r_sqr;
 y /= 2.0;
 
 return y;
}/* Anum2HOInt */


double Glauber::NuInt2HO(double xi)
{
 double f;
 double z_sqr, r_sqr, s;
 double w_WS, a_WS, rho;
 
 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 

/* devide by a_WS, make life simpler */
   
   s = NuInS_S/a_WS;
   z_sqr = -log(xi);
   r_sqr = s*s + z_sqr;

/* no need to divide by 2 here because -infty < z < infty and 
   we integrate only over positive z */
   
   if(z_sqr < 0.0) z_sqr = tiny;
   f = a_WS*rho*(LexusData.SigmaNN);
   f *= (1.0 + w_WS*r_sqr)*exp(-s*s)/sqrt(z_sqr);
 
 return f;
}/* NuInt2HO */




int Glauber::IsFile(char *file_name)
{
 FILE *temp;

 if( (temp = fopen(file_name,"r")) == NULL) return 0;
 else 
  {
   fclose(temp);
   return 1;
  }
}/* IsFile */

char *Glauber::StringFind(char *file_name, char *st)
{
  char *x = char_malloc(120);
  char *s = char_malloc(120);

  FILE *input;
  input = fopen(file_name,"r");

  int bytes_read=fscanf(input, "%s", s);
  int ind = 0;

  while(strcmp(s, "EndOfData") != 0)
  {
    bytes_read=fscanf(input, "%s", x);
    if(strcmp(s, st) == 0)
    {
      ind++;
      fclose(input);
      free(s);
      return x;
    }/* if right, return */
    free(s);
    s = char_malloc(120);
    bytes_read=fscanf(input, "%s", s);
  }/* while */
  free(s);
  fclose(input);

  if(ind == 0)
  {
    fprintf(stdout, "%s not found in %s.\n", st, file_name);
    printf("Enter %s = ", st);
    bytes_read=scanf("%s", x);
    printf("Rewriting %s...\n", file_name);
    ReWriteString(file_name, st, x);
    char_free(x);
    char_free(s);
    return x;
  }
  strcpy(x, "N/A");
  return(x);
}/* StringFind */

double Glauber::DFind(char *file_name, char *st)
{
 char *s;
 double x;
 
 s = StringFind(file_name, st);

 sscanf(s, "%lf", &x);
 return x;

}/* DFind */

int Glauber::integer(double x)
{
 return (int)(x+0.5);
}

char *Glauber::char_malloc(int n1)
{
    char *char_ptr;

    /* pointer to the n1 array */
    char_ptr = (char *) malloc (sizeof(char)*n1);

return char_ptr;
}

double *Glauber::vector_malloc(int n1)
{
 double *d1_ptr;

    /* pointer to the n1 array */
    d1_ptr = (double *) malloc (sizeof(double )*n1);

 return d1_ptr;
}

void Glauber::char_free(char *vec)
{
 free(vec);
}

double Glauber::integral (int id, double down, double up, double tol, int *count)
{
  double dx, y, g1[7];
  int i;
  
  if (down == up) y = 0.0;
  else
    {
      dx = (up-down)/6.0;
      for( i=0; i<7; i++) 
	{
	  if (id==1) g1[i] = NuInt2HO(down + i*dx);
	  else if (id==2) g1[i] = NuInt3Gauss(down + i*dx);
	  else if (id==3) g1[i] = NuInt3Fermi(down + i*dx);
	  else if (id==4) g1[i] = Anum3FermiInt(down + i*dx);
	  else if (id==5) g1[i] = Anum3GaussInt(down + i*dx);
	  else if (id==6) g1[i] = Anum2HOInt(down + i*dx);
	  else if (id==7) g1[i] = OLSIntegrand(down + i*dx);
	}
      *count = 7;
      y = qnc7(id, tol, down, dx, g1, 0.0, 0.0, count);
    }
  return y;
} /* end of integral */

double Glauber::qnc7(int id, double tol, double down, double dx, double *f_of, 
		     double pre_sum, double area, int *count)
{
  int i;
  double left_sum, right_sum, ans;
  static double w[] = 
    {41.0/140.0, 54.0/35.0, 27.0/140.0, 68.0/35.0, 27.0/140, 54.0/35.0,
       41.0/140.0};
  double fl[7];
  double fr[7];
  /*
    qnc7 calculates integral over left and right half of the given interval
    and branches
    to do so, first halve dx
    */
  
  dx /= 2.0;

  /*
    first calculate the left estimate
    f_of[] contains the evaluated values at down+i, 0< i <7
    store half distanced values for the left sum in fl[]
    */

  
  if (id==1)
    {
      fl[1] = NuInt2HO(down + dx);
      fl[3] = NuInt2HO(down + 3.0*dx);
      fl[5] = NuInt2HO(down + 5.0*dx);
    }
  else if (id==2)
    {
      fl[1] = NuInt3Gauss(down + dx);
      fl[3] = NuInt3Gauss(down + 3.0*dx);
      fl[5] = NuInt3Gauss(down + 5.0*dx);
    }
  else if (id==3)
    {
      fl[1] = NuInt3Fermi(down + dx);
      fl[3] = NuInt3Fermi(down + 3.0*dx);
      fl[5] = NuInt3Fermi(down + 5.0*dx);
    }
  else if (id==4)
    {
      fl[1] = Anum3FermiInt(down + dx);
      fl[3] = Anum3FermiInt(down + 3.0*dx);
      fl[5] = Anum3FermiInt(down + 5.0*dx);
    }  
  else if (id==5)
    {
      fl[1] = Anum3GaussInt(down + dx);
      fl[3] = Anum3GaussInt(down + 3.0*dx);
      fl[5] = Anum3GaussInt(down + 5.0*dx);
    }  
  else if (id==6)
    {
      fl[1] = Anum2HOInt(down + dx);
      fl[3] = Anum2HOInt(down + 3.0*dx);
      fl[5] = Anum2HOInt(down + 5.0*dx);
    }  
  else if (id==7)
    {
      fl[1] = OLSIntegrand(down + dx);
      fl[3] = OLSIntegrand(down + 3.0*dx);
      fl[5] = OLSIntegrand(down + 5.0*dx);
    }  
  
  fl[0] = f_of[0];
  fl[2] = f_of[1];
  fl[4] = f_of[2];
  fl[6] = f_of[3];
  
  *count += 3;

  left_sum = 0.0;
  for(i=0; i<7; i++) left_sum += w[i]*fl[i];
  left_sum *= dx;

/*printf("leftsum is %le\n", left_sum);*/

  /*
    like wise, the right sum is in fr[]
    */

  if (id==1)
    {
      fr[1] = NuInt2HO(down + 7.0*dx);
      fr[3] = NuInt2HO(down + 9.0*dx);
      fr[5] = NuInt2HO(down + 11.0*dx);
    }
  else if (id==2)
    {
      fr[1] = NuInt3Gauss(down + 7.0*dx);
      fr[3] = NuInt3Gauss(down + 9.0*dx);
      fr[5] = NuInt3Gauss(down + 11.0*dx);
    }
  else if (id==3)
    {
      fr[1] = NuInt3Fermi(down + 7.0*dx);
      fr[3] = NuInt3Fermi(down + 9.0*dx);
      fr[5] = NuInt3Fermi(down + 11.0*dx);
    }
  else if (id==4)
    {
      fr[1] = Anum3FermiInt(down + 7.0*dx);
      fr[3] = Anum3FermiInt(down + 9.0*dx);
      fr[5] = Anum3FermiInt(down + 11.0*dx);
    }
  else if (id==5)
    {
      fr[1] = Anum3GaussInt(down + 7.0*dx);
      fr[3] = Anum3GaussInt(down + 9.0*dx);
      fr[5] = Anum3GaussInt(down + 11.0*dx);
    }
  else if (id==6)
    {
      fr[1] = Anum2HOInt(down + 7.0*dx);
      fr[3] = Anum2HOInt(down + 9.0*dx);
      fr[5] = Anum2HOInt(down + 11.0*dx);
    }
  else if (id==7)
    {
      fr[1] = OLSIntegrand(down + 7.0*dx);
      fr[3] = OLSIntegrand(down + 9.0*dx);
      fr[5] = OLSIntegrand(down + 11.0*dx);
    }
  
  fr[0] = f_of[3];
  fr[2] = f_of[4];
  fr[4] = f_of[5];
  fr[6] = f_of[6];

  *count += 3;

  right_sum = 0.0;
  for(i=0; i<7; i++) right_sum += w[i]*fr[i];
  right_sum *= dx;

/*printf("rightsum is %le\n", right_sum);*/

  ans = left_sum + right_sum;

/*printf("ans is %le\n", ans);*/


  /* 
    up date total area subtract previously assigned area for this interval
    and add newly calculated area
    */
/*printf("pre_area is %le\n", area);*/

  area += -fabs(pre_sum) + fabs(left_sum) + fabs(right_sum);

  /* 
    printf("presum is %le\n", pre_sum);
    printf("area is %le\n", area);
    */
  /*
    branch if the refined sum is finer than the previous estimate
    */

  if( fabs(ans - pre_sum) > tol*fabs(area) && (*count < limit))
    {
      /*
	branch by calling the function itself
	by calling the qnc7 twice, we are branching
	since left hand side is being calculated first, until the condition
	is satisfied, the left branch keeps branching
	when finally the condition is met by one left-most interval, 
	qnc7 returns the right hand side of one up level,
	and the same process resumes 
	until the criterion is met by all the branched
	intervals,
	then qnc7 returns to the original right branch and resumes halving 
	until the condition is met by all intervals
	(funk, down, dx, f_of[7], pre_ans, ans)
	*/

      tol /= 1.414;
      left_sum = qnc7(id, tol, down, dx, fl, left_sum, area, count);
      right_sum = qnc7(id, tol, down+dx*6, dx, fr, right_sum, area, count);

      ans = left_sum + right_sum;

      /* printf("ans is %le\n", ans);*/

    }/* belongs to if*/

  return ans;
} /* end of qnc */

void Glauber::ReWriteString(char *file_name, char *st, char *x)
{
  int bytes_read;
  char tmp_file[100];
  FILE *output;
  
  /* tmp_file = tmpnam(); this does not work with gcc */
  
  //tmpnam(tmp_file); /* this works with both cc and gcc */
  
  bytes_read=mkstemp(tmp_file);
  
  output = fopen(tmp_file, "w");
  fprintf(output, "%s %s\n", st, x);
  fclose(output);
  
  FileCat(file_name, tmp_file);
  
  FileCopy(tmp_file, file_name);
  
  remove(tmp_file);
  
}/* ReWriteString */


void Glauber::FileCopy(char *in_file, char *out_file)
{
  FILE *output, *input;
  char c; 
  
  input = fopen(in_file, "r");
  output = fopen(out_file, "w");
  
  c = getc(input);
  while(c != EOF)
    {
      fprintf(output, "%c", c);
      c = getc(input);
    }
  
  fclose(input);
  fclose(output);
}/* FileCopy */


void Glauber::FileCat(char *in_file, char *out_file)
{
  FILE *output, *input;
  char c; 
  
  input = fopen(in_file, "r");
  output = fopen(out_file, "a");
  
  c = getc(input);
  while(c != EOF)
    {
      fprintf(output, "%c", c);
      c = getc(input);
    }
  
  fclose(input);
  fclose(output);
}/* FileCopy */


double Glauber::OLSIntegrand(double s)
{
  double sum, arg, x, r;
  int k, m;
  m = 20;
  sum = 0.0;
  for(k=1; k<=m; k++)
    {
      arg = M_PI*(2.0*k - 1.0)/(2.0*m);
      x = cos(arg);
      r = sqrt(s*s + b*b + 2.0*s*b*x);
      sum += InterNuTInST(r);
    }/* k */
  
  return s*sum*M_PI/(1.0*m)*InterNuPInSP(s);
  
}/* OLSIntegrand */


double Glauber::TAB()
{
  double f;
  int count = 0;
  f = integral(7, 0.0, LexusData.SCutOff, TOL, &count); // integrate OLSIntegrand(s)
  f *= 2.0/(LexusData.SigmaNN); //here TAB is the number of binary collisions, dimensionless 
                                //(1/fm^4 integrated over dr_T^2 (gets rid of 1/fm^2), divided by sigma (gets rid of the other))
  return f;
}/* TAB */


double Glauber::PAB(double x, double y)
{
  double result;
  double s1=sqrt(pow(x+b/2.,2.)+y*y);
  double s2=sqrt(pow(x-b/2.,2.)+y*y);
  result = InterNuPInSP(s1)*InterNuTInST(s2)/(currentTAB*LexusData.SigmaNN);
  if (glauberEnvelope) 
    result *= (1.4+1.4*tanh((sqrt(0.0001*pow(y,4.))))-1.4*tanh(sqrt(0.01*pow(x,4)))); 
  return result;
}/* PAB */


void Glauber::preInit(double SigmaNN, string Target, string Projectile, double inb, int imax)
{
  char* Target_Name;
  Target_Name = (char*)malloc((Target.length()+10) * sizeof(char));
  memset((void*)Target_Name,0,sizeof(Target_Name));
  strcpy(Target_Name,Target.c_str());

  char *Projectile_Name;
  Projectile_Name = (char*)malloc((Projectile.length()+10) * sizeof(char));
  memset((void*)Projectile_Name,0,sizeof(Projectile_Name));
  strcpy(Projectile_Name,Projectile.c_str());

  /*
  string path = "";
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
    }
  else path = ".";
  string pathAndFile = path+"/data/known_nuclei.dat";

  char *paf;
  paf = (char*)malloc((pathAndFile.length()+10) * sizeof(char));
  memset((void*)paf,0,sizeof(paf));
  strcpy(paf,pathAndFile.c_str());
  
  cout << "opening file " << paf << " for reading ... " << endl;
  
  if(!IsFile(paf))
    {
      fprintf(stderr, "No known_nuclei.dat.  Please provide one.\n");
      fprintf(stderr, "Exiting...\n");
      exit(0);
    }
  */

  /* energy unit is always GeV and length unit is fm */
  //FindNucleusData(&(LexusData.Target), Target_Name, paf);
  lookup_nucleus_info(&(LexusData.Target), Target_Name);
  Target_A = LexusData.Target.A;
  //PrintNucleusData(&(LexusData.Target));
  
  //FindNucleusData(&(LexusData.Projectile), Projectile_Name, paf);
  lookup_nucleus_info(&(LexusData.Projectile), Projectile_Name);
  //PrintNucleusData(&(LexusData.Projectile));
  Projectile_A = LexusData.Projectile.A;

  PrintCollisionSystem(&LexusData);

  currentA = max(Projectile_A, Target_A);
  //currentA = LexusData.Projectile.A;
  //cout << " currentA=" << currentA << endl;
  free(Target_Name);
  free(Projectile_Name);
}

void Glauber::init(double SigmaNN, string Target, string Projectile, double inb, int imax, bool glauberEnvelopeIn)
{
  glauberEnvelope = glauberEnvelopeIn;

  char* Target_Name;
  Target_Name = (char*)malloc((Target.length()+10) * sizeof(char));
  memset((void*)Target_Name,0,sizeof(Target_Name));
  strcpy(Target_Name,Target.c_str());

  char *Projectile_Name;
  Projectile_Name = (char*)malloc((Projectile.length()+10) * sizeof(char));
  memset((void*)Projectile_Name,0,sizeof(Projectile_Name));
  strcpy(Projectile_Name,Projectile.c_str());

  /*
  string path = "";
  const char* MARTINIPATH = "MARTINIPATH";
  char* envPath = getenv(MARTINIPATH);
  
  if (envPath != 0 && *envPath != '\0') 
    {
      int i = 0;
      while (*(envPath+i) != '\0') path += *(envPath+(i++));
    }
  else path = ".";
  string pathAndFile = path+"/data/known_nuclei.dat";

  char *paf;
  paf = (char*)malloc((pathAndFile.length()+10) * sizeof(char));
  memset((void*)paf,0,sizeof(paf));
  strcpy(paf,pathAndFile.c_str());
  
  cout << "opening file " << paf << " for reading ... " << endl;
  
  if(!IsFile(paf))
    {
      fprintf(stderr, "No known_nuclei.dat.  Please provide one.\n");
      fprintf(stderr, "Exiting...\n");
      exit(0);
    }
  */

  /* pp total cross-section : 40 mb = 4 fm**2 */
  /* LexusData.SigmaNN = 4.0; */
  
  /* energy unit is always GeV and length unit is fm */
  //FindNucleusData(&(LexusData.Target), Target_Name, paf);
  lookup_nucleus_info(&(LexusData.Target), Target_Name);
  Target_A = LexusData.Target.A;
  //PrintNucleusData(&(LexusData.Target));
  
  //FindNucleusData(&(LexusData.Projectile), Projectile_Name, paf);
  lookup_nucleus_info(&(LexusData.Projectile), Projectile_Name);
  //PrintNucleusData(&(LexusData.Projectile));
  Projectile_A = LexusData.Projectile.A;

  PrintCollisionSystem(&LexusData);

  LexusData.SigmaNN = 0.1*SigmaNN; // sigma in fm^2 
  
  currentA = max(Projectile_A, Target_A);
  
  fprintf(stderr, "SigmaNN = %le\n", LexusData.SigmaNN);
  
  /* Run Specific */
  
  LexusData.InterMax = imax;
  LexusData.SCutOff = 12.;
  PrintLexusData();
  
  b=inb; 
  //currentTAB=TAB();
  
  free(Target_Name);
  free(Projectile_Name);
}/* init */


ReturnValue Glauber::SamplePAB(Random *random) //slow 
{
  // samples the initial position of the nucleon-nucleon collision according to the distribution P_AB(x,y,b)
  // using a 2D-Metropolis algorithm
  ReturnValue returnVec;
  double x, y;
  double x_new, y_new;
  double g = 0, g_new = 0;
  double ratio;
  
  // the ranges in which the variables u and phi need to be sampled
  const double x_min = 0.0;
  const double x_max = LexusData.SCutOff;
  const double y_min = 0.0;
  const double y_max = LexusData.SCutOff;
  
  do
    {
      x = x_min + random->genrand64_real1()*(x_max - x_min);
      y = y_min + random->genrand64_real1()*(y_max - y_min);
      g = PAB(x,y);
    } while(g==0.0);
  
  cout.precision(6);

  // number of steps in the Markov chain
  const int n_steps =25;
  
  for(int i=0; i<n_steps; i++)
    {
      do
	{
	  x_new = x_min + random->genrand64_real1()*(x_max - x_min);  
	  // propose new x using a uniform distribution over the entire range
	} while( x_new < x_min || x_new > x_max);
      
      do
	{
	  y_new = y_min + random->genrand64_real1()*(y_max - y_min);
	  // propose new y using a uniform distribution over the entire range
	} while( y_new < y_min || y_new > y_max);  
      
      g_new = PAB(x_new,y_new);
      
      ratio = g_new / g;        // ratio of  g(x',y') / g(x,y)
      
      // accept if g(x',y') > g(x,y) or with probability g(x',y') / g(x,y)
      if ( ratio>=1.0 || random->genrand64_real1() < ratio )
	{
	  x = x_new;
	  y = y_new;
	  g = g_new;
	}
    }
  if(random->genrand64_real1() <= 0.5) returnVec.x = x;
  else returnVec.x = -x;
  if(random->genrand64_real1() <= 0.5) returnVec.y = y;
  else returnVec.y = -y;
  return returnVec; 
}

double Glauber::areaPAB(double x, double A)
{
  double f;
  //f = 0.5*height*gamma*(- log(gamma)+log(gamma+x*x));
  f = (3710.33-3710.33*exp((-0.05-0.005*b)*x*x))*A*(4.+b)/(10.+b);
  return f;
}

ReturnValue Glauber::SamplePABRejection(Random *random)
{
  ReturnValue returnVec;

  double r, x, y, tmp;
  double phi;
  cout.precision(10);
  double A=0.00015;
  /*
  for (int tp=0; tp<=1000; tp++)
    {
      double dx=0.01;
      cout << dx*tp << " " << PAB(dx*tp,0.) << endl;
    }
  
  cout << endl;
  for (int tp=0; tp<=1000; tp++)
    {
      double dx=0.01;
      cout << dx*tp << " " << PAB(0.,dx*tp) << endl;
    }
  exit(1);
  */
  if(currentA<177.)
    A*=197./(197.-currentA-20.);
  gamma = 13.;
  height = 0.022*(1.+b/10.); //PAB(0.,0.)/gamma;
  do
    {
      phi = 2.*M_PI*random->genrand64_real1();
      //r = sqrt(-1.+exp(2.*(areaPAB(LexusData.SCutOff)*random->genrand64_real1())/height/gamma))*sqrt(gamma);
      r = 4.47213595499958*sqrt(-log(1.+(-0.0026951787996341873-0.00026951787996341873*b)
				     *(areaPAB(LexusData.SCutOff,A)*random->genrand64_real1())/(A*(4.+b))))/sqrt(1.+0.1*b);
      // here random->genrand64_real1()*areaPAB(LexusData.SCutOff) 
      // is a uniform random number on [0, area under f(x)]
      tmp = random->genrand64_real1();
      // x is uniform on [0,1]
      //i++;
      x=r*cos(phi);
      y=r*sin(phi);
      //cout << "r=" << r  << " x=" << x << " y=" << y << endl;
      if ((r*PAB(x,y))>A*(1. + 0.25*b)*exp(5. + 0.05*(-1. - 0.1*b)*r*r)*r) 
	cout << "WARNING: PAB>envelope: " << "r*PAB=" << PAB(x,y)*r 
	     << ", f=" << A*(1. + 0.25*b)*exp(5. + 0.05*(-1. - 0.1*b)*r*r)*r << endl;
    } while( tmp > (r*PAB(x,y))/(A*(1. + 0.25*b)*exp(5. + 0.05*(-1. - 0.1*b)*r*r)*r)); 
  // reject if tmp is larger than the ratio p(y)/f(y), f(y)=height*gamma/(y*y+gamma)
  returnVec.x=x;
  returnVec.y=y;
  return returnVec; 
}

double Glauber::areaTA(double x, double A)
{
  double f;
  f = A*220.*(1.-exp(-0.025*x*x)); 
  return f;
}

ReturnValue Glauber::SampleTARejection(Random *random)
{
  ReturnValue returnVec;

  double r, x, y, tmp;
  double phi;
  double A=1.2*LexusData.SigmaNN/4.21325504715; // increase the envelope for larger sigma_inel (larger root_s) (was originally written
  // for root(s)=200 GeV, hence the cross section of 4.21325504715 fm^2 (=42.13 mb)
  cout.precision(10);
  do
    {
      phi = 2.*M_PI*random->genrand64_real1();
      r = 6.32456*sqrt(-log((-0.00454545*(-220.*A+areaTA(15.,A)*random->genrand64_real1()))/A));
      // here random->genrand64_real1()*areaTA(LexusData.SCutOff) 
      // is a uniform random number on [0, area under f(x)]
      tmp = random->genrand64_real1();
      // x is uniform on [0,1]
      if ( r*InterNuPInSP(r) > A*r*11.*exp(-r*r/40.) ) 
	cout << "WARNING: TA>envelope: " << "TA=" << r*InterNuPInSP(r) 
	     << ", f=" << A*r*11.*exp(-r*r/40.) << endl;
    } while( tmp > r*InterNuPInSP(r)/(A*r*11.*exp(-r*r/40.))); 
  // reject if tmp is larger than the ratio p(y)/f(y), f(y)=A*r*11.*exp(-r*r/40.))
  x=r*cos(phi);
  y=r*sin(phi);
  returnVec.x=x;
  returnVec.y=y;
  return returnVec; 
}
