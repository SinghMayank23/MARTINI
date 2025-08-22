

PreHydroInfo HydroSetup::getPreHydroValues(double x, double y, double z, double t)
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

    int ieta = 0.;

    int itau = floor((tau-prehydroTau0)/prehydroDtau+0.0001);       
    int ix = floor((prehydroXmax+x)/prehydroDx+0.0001);
    int iy = floor((prehydroXmax+y)/prehydroDx+0.0001);

    double xfrac = (x-( (double)ix*prehydroDx-prehydroXmax ))/prehydroDx;
    double yfrac = (y-( (double)iy*prehydroDx-prehydroXmax ))/prehydroDx;
    double taufrac = (tau-prehydroTau0)/prehydroDtau-(double)itau;
    double etafrac = 0.;

    PreHydroInfo info;
    if ( ix < 0 || ix >= iprexmax ) 
    {
        cout << "[MARTINI:HydroSetup::getPreHydroValues]: "
             << "WARNING - x out of range x=" << x 
             << ", ix=" << ix << ", ixmax=" << iprexmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta 
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau 
             << " itau=" << itau << " itaumax=" << ipretaumax << endl;

        info.T = 0.0;
        info.ux = 0.0;
        info.uy = 0.0;
        info.utau = 1.0;
        return(info);
    }
    if ( iy < 0 || iy >= iprexmax )
    {
        cout << "[MARTINI:HydroSetup::getPreHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy 
             << ", iymax=" << iprexmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta 
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau 
             << " itau=" << itau << " itaumax=" << ipretaumax << endl;

        info.T = 0.0;
        info.ux = 0.0;
        info.uy = 0.0;
        info.utau = 1.0;
        return(info);
    }
    if ( itau < 0 || itau >= ipretaumax )
    {
        cout << "[MARTINI:HydroSetup::getPreHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << ipretaumax 
             << endl;

        info.T = 0.0;
        info.ux = 0.0;
        info.uy = 0.0;
        info.utau = 1.0;
        return(info);
    }
  
  //The array of positions on the 4-dimensional rectangle:
  int position[2][2][2][2];
  for(int ipx=0; ipx<2; ipx++)
  {
      int px;
      if(ipx==0 || ix==iprexmax-1 ) 
          px = ix;
      else 
          px = ix+1;
      for(int ipy=0; ipy<2; ipy++)
      {
          int py;
          if(ipy==0 || iy==iprexmax-1 )
              py = iy;
          else 
              py = iy+1;
          for(int ipeta=0; ipeta<2; ipeta++)
          {
              int peta;
              if(ipeta==0 || ieta==ipreetamax-1 ) 
                  peta = ieta;
              else 
                  peta = ieta+1;
              for(int iptau=0; iptau<2; iptau++)
              {
                  int ptau;
                  if(iptau==0 || itau==ipretaumax-1 ) 
                      ptau = itau;
                  else 
                      ptau = itau+1;
                  position[ipx][ipy][ipeta][iptau] = (
                                  py+iprexmax*(px+iprexmax*(peta+ipreetamax*ptau)));
              }
          }
      }
  }

  //And now, the interpolation:
  double T;
  double ux;
  double uy;
  double utau;
  T = ux = uy = utau = 0.;

  PreHydroCell PreHydroCell_temp1, PreHydroCell_temp2;
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

              PreHydroCell_temp1 = (*prehydro_lattice)[position[0][ipy][ipeta][iptau]];
              PreHydroCell_temp2 = (*prehydro_lattice)[position[1][ipy][ipeta][iptau]];

              T += prefrac*( (1. - xfrac)*PreHydroCell_temp1.T
                              + xfrac*PreHydroCell_temp2.T );
              ux += prefrac*( (1. - xfrac)*PreHydroCell_temp1.ux
                              + xfrac*PreHydroCell_temp2.ux );
              uy += prefrac*( (1. - xfrac)*PreHydroCell_temp1.uy
                              + xfrac*PreHydroCell_temp2.uy );
              utau += prefrac*( (1. - xfrac)*PreHydroCell_temp1.utau
                              + xfrac*PreHydroCell_temp2.utau );
          }
      }
  }
  

  info.T = T;
  info.ux = ux;
  info.uy = uy;
  info.utau = utau;
  
  return(info);
}
