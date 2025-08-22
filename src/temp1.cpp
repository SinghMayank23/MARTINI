
HydroInfo HydroSetup::getHydroValues(double x, double y, double z, double t)
{
    double tau, eta;
    if(use_tau_eta_coordinate == 1)
    {
        if(t*t>z*z)
        {
            tau = sqrt(t*t-z*z);
            eta = 0.5*log((t+z)/(t-z));
        }
        else
        {
            tau = 0.;
            eta = 0.;
        }
    }
    else
    {
        // if the medium is given in cartesian coordinates 
        // set tau and eta to t and z
        tau = t;
        eta = z;
    }

    int ieta = floor((hydroZmax+eta)/hydroDz+0.0001);
    if(hydroWhichHydro == 8 || hydroWhichHydro == 6)
        ieta = 0;

    int itau = floor((tau-hydroTau0)/hydroDtau+0.0001);       
    int ix = floor((hydroXmax+x)/hydroDx+0.0001);
    int iy = floor((hydroXmax+y)/hydroDx+0.0001);

    double xfrac = (x-( (double)ix*hydroDx-hydroXmax ))/hydroDx;
    double yfrac = (y-( (double)iy*hydroDx-hydroXmax ))/hydroDx;
    double etafrac = eta/hydroDz-(double)ieta+0.5*(double)ietamax;
    double taufrac = (tau-hydroTau0)/hydroDtau-(double)itau;

    HydroInfo info;
    if ( ix < 0 || ix >= ixmax ) 
    {
        cout << "[MARTINI:HydroSetup::getHydroValues]: "
             << "WARNING - x out of range x=" << x 
             << ", ix=" << ix << ", ixmax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta 
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau 
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info.T = 0.0;
        info.QGPfrac = 0.0;
        info.vx = 0.0;
        info.vy = 0.0;
        info.vz = 0.0;
        info.veta = 0.0;
        return(info);
    }
    if ( iy < 0 || iy >= ixmax )
    {
        cout << "[MARTINI:HydroSetup::getHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy 
             << ", iymax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta 
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau 
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info.T = 0.0;
        info.QGPfrac = 0.0;
        info.vx = 0.0;
        info.vy = 0.0;
        info.vz = 0.0;
        info.veta = 0.0;
        return(info);
    }
    if ( itau < 0 || itau >= itaumax )
    {
        cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << itaumax 
             << endl;
        cout << "[MARTINI:HydroSetup::getHydroValues]: tau= " << tau 
             << ", hydroTauMax = " << hydroTauMax 
             << ", hydroDtau = " << hydroDtau << endl;

        info.T = 0.0;
        info.QGPfrac = 0.0;
        info.vx = 0.0;
        info.vy = 0.0;
        info.vz = 0.0;
        info.veta = 0.0;
        return(info);
    }
    if ( ieta < 0 || ieta >= ietamax  ) 
    {
        cout << "[MARTINI:HydroSetup::getHydroValues]: WARNING - "
             << "eta out of range, ieta=" << ieta << ", ietamax=" << ietamax 
             << endl;
        
        info.T = 0.0;
        info.QGPfrac = 0.0;
        info.vx = 0.0;
        info.vy = 0.0;
        info.vz = 0.0;
        info.veta = 0.0;
        return(info);
    }
  
  //The array of positions on the 4-dimensional rectangle:
  int position[2][2][2][2];
  for(int ipx=0; ipx<2; ipx++)
  {
      int px;
      if(ipx==0 || ix==ixmax-1 ) 
          px = ix;
      else 
          px = ix+1;
      for(int ipy=0; ipy<2; ipy++)
      {
          int py;
          if(ipy==0 || iy==ixmax-1 )
              py = iy;
          else 
              py = iy+1;
          for(int ipeta=0; ipeta<2; ipeta++)
          {
              int peta;
              if(ipeta==0 || ieta==ietamax-1 ) 
                  peta = ieta;
              else 
                  peta = ieta+1;
              for(int iptau=0; iptau<2; iptau++)
              {
                  int ptau;
                  if(iptau==0 || itau==itaumax-1 ) 
                      ptau = itau;
                  else 
                      ptau = itau+1;
                  position[ipx][ipy][ipeta][iptau] = (
                                  px+ixmax*(py+ixmax*(peta+ietamax*ptau)));
              }
          }
      }
  }

  //And now, the interpolation:
  double T;
  double QGPfrac;
  double vx;
  double vy;
  double vz;
  double veta;
  T = QGPfrac = vx = vy = vz = veta = 0.;

  HydroCell HydroCell_temp1, HydroCell_temp2;
  for(int iptau = 0; iptau < 2; iptau++)
  {
      double taufactor;
      if(iptau == 0) 
          taufactor = 1. - taufrac;
      else 
          taufactor = taufrac;
      for(int ipeta = 0; ipeta < 2; ipeta++)
      {
          double etafactor;
          if(ipeta == 0) 
              etafactor = 1. - etafrac;
          else 
              etafactor = etafrac;
          for(int ipy = 0; ipy < 2; ipy++)
          {
              double yfactor;
              if(ipy == 0) 
                  yfactor = 1. - yfrac;
              else 
                  yfactor = yfrac;
                  
              double prefrac = yfactor*etafactor*taufactor;

              HydroCell_temp1 = (*lattice)[position[0][ipy][ipeta][iptau]];
              HydroCell_temp2 = (*lattice)[position[1][ipy][ipeta][iptau]];

              T += prefrac*( (1. - xfrac)*HydroCell_temp1.T
                              + xfrac*HydroCell_temp2.T );
              vx += prefrac*( (1. - xfrac)*HydroCell_temp1.vx
                              + xfrac*HydroCell_temp2.vx );
              vy += prefrac*( (1. - xfrac)*HydroCell_temp1.vy
                              + xfrac*HydroCell_temp2.vy );
              if (hydroWhichHydro != 8)
              {
                  QGPfrac += prefrac*( (1. - xfrac)*HydroCell_temp1.QGPfrac
                                       + xfrac*HydroCell_temp2.QGPfrac );
                  vz += prefrac*( (1. - xfrac)*HydroCell_temp1.vz
                                  + xfrac*HydroCell_temp2.vz );
              }
          }
      }
  }
  
  if (hydroWhichHydro == 8)   // for boost invariant medium
  {
      if (T < hydroTfinal)
          QGPfrac = 0.;
      else
          QGPfrac = 1.;
      vz = z/t;               // Bjorken flow in the lab frame 
      double gamma_L_inv = sqrt(1. - vz*vz);
      vx = vx*gamma_L_inv;    // convert vx and vy to lab frame
      vy = vy*gamma_L_inv;
      veta = 0.0;
  }

  info.T = T;
  info.QGPfrac = QGPfrac;
  info.vx = vx;
  info.vy = vy;
  info.vz = vz;
  info.veta = veta;
  
  return(info);
}

