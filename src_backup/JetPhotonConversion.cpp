#include "JetPhotonConversion.h"

ConversionPhotonGrid::ConversionPhotonGrid(std::string grid_info_file, std::string grid_data_file)
{
    std::ifstream settings (grid_info_file, std::ios::in);
    double temp;
    std::string label;
    while (!settings.eof() && settings.is_open())
    {
      settings >> label >> temp;
      if (label == "jet_energy_grid_size")
        jet_energy_grid_size = int(temp);
      else if (label == "temperature_grid_size")
        temperature_grid_size = int(temp);
      else if (label == "photon_energy_grid_size")
        photon_energy_grid_size = int(temp);
      else if (label == "jet_energy_min")
        jet_energy_min = temp;
      else if (label == "jet_energy_max")
        jet_energy_max = temp;    
      else if (label == "temperature_min")
        temperature_min = temp;
      else if (label == "temperature_max")
        temperature_max = temp;
      else if (label == "photon_energy_min")
        photon_energy_min = temp;
      else if (label == "photon_energy_max")
        photon_energy_max = temp;
      else
      {
        // a tag that doesn't match the above is uselse. Next!
        continue;
      }
    }
    settings.close();

    dEjet = (jet_energy_max-jet_energy_min)/jet_energy_grid_size;
    dT = (temperature_max - temperature_min)/temperature_grid_size;
    dEgamma = (photon_energy_max - photon_energy_min)/photon_energy_grid_size;
    std::cout<<"jet_energy_grid_size: "   <<jet_energy_grid_size;
    std::cout<<", photon_energy_grid_size: "<<photon_energy_grid_size; 
    std::cout<<", temperature_grid_size: "  <<temperature_grid_size<<std::endl;
    // now read file and populate the grid
    FILE* fin = fopen(grid_data_file.c_str(), "r+b");
    //std::cout<<fin<<std::endl;
    //std::cout<<(fin==NULL)<<std::endl;
    if (fin == NULL) 
    {
        std::cout<<"\x1B[31m"<<"Failed to open file. Abort.";
        std::cout<<"Initialize the grid to zeros and return."<<"\033[0m"<<std::endl;
        grid = new double**[jet_energy_grid_size];
        for(int i = 0; i < jet_energy_grid_size; i++)
        {
            grid[i] = new double*[temperature_grid_size];
            for (int j = 0;  j < temperature_grid_size; j++)
            {
                grid[i][j] = new double[photon_energy_grid_size];
                for (int k = 0; k < photon_energy_grid_size; k++)
                    grid[i][j][k] = 0;
            }
        }
    }
    else
    {
        int status;
        int size = sizeof(double);
        grid = new double**[jet_energy_grid_size];
        for(int i = 0; i < jet_energy_grid_size; i++)
        {
            grid[i] = new double*[temperature_grid_size];
            for (int j = 0;  j < temperature_grid_size; j++)
            {
                grid[i][j] = new double[photon_energy_grid_size];
                for (int k = 0; k < photon_energy_grid_size; k++)
                {
                    double rate;
                    status = fread(&rate, size, 1, fin);
                    if (status)
                    {
                        //cout<<"i: "<<i<<" j: "<< j <<" k: "<< k<<" ==> rate: "<<rate<<endl;
                        grid[i][j][k] = rate;
                    }
                }
            }
        }
        fclose(fin);
    }
}

ConversionPhotonGrid::~ConversionPhotonGrid()
{
    for(int i = 0; i < jet_energy_grid_size; i++)
    {
        for (int j = 0;  j < temperature_grid_size; j++)
        {
            delete[] grid[i][j];
        }
        delete[] grid[i];
    }
}
//interpolate the values stored in the grid to get the rate for the given
// set of numbers. Trilinear interpolator.
double ConversionPhotonGrid::getConvRate_Trilinear(double eJet, double T, double eGamma)
{
    int ijet = floor((eJet - jet_energy_min)/dEjet);
    int iT = floor((T - temperature_min)/dT);
    int iEgamma = floor((eGamma - photon_energy_min)/dEgamma);
    //std::cout<<"I'm here in the trilinear converison function"<<std::endl;
    if (ijet+1 > jet_energy_grid_size || iT+1 > temperature_grid_size || iEgamma+1 > photon_energy_grid_size)
    {
        std::cout<<"\x1B[31m";
        std::cout<<"Return 0 for the rate as it is outside of the grid. No extrapolation for now.";
        std::cout<<"\033[0m"<<std::endl;
        return 0;
    }
    if (eGamma < photon_energy_min || eGamma > photon_energy_max || eJet < jet_energy_min || eJet > jet_energy_max || T > temperature_max || T < temperature_min )
    {
        //outside of the region of validity for our table. 
        return 0;
    }
    double E1 = jet_energy_min + ijet*dEjet;
    double E2 = jet_energy_min + (ijet+1)*dEjet;

    double T1 = temperature_min + iT*dT;
    double T2 = temperature_min + (iT+1)*dT;

    double W1 = photon_energy_min + iEgamma*dEgamma;
    double W2 = photon_energy_min + (iEgamma+1)*dEgamma;

    double volume = (E2-E1)*(T2-T1)*(W2-W1); // total volume surrounding our point.
    double w111 = (E2 - eJet)*(T2 - T)*(W2 - eGamma)/volume;
    double w112 = (E2 - eJet)*(T2 - T)*(eGamma - W1)/volume;
    double w121 = (E2 - eJet)*(T - T1)*(W2 - eGamma)/volume;
    double w122 = (E2 - eJet)*(T - T1)*(W2 - eGamma)/volume;
    double w211 = (eJet - E1)*(T2 - T)*(W2 - eGamma)/volume;
    double w212 = (eJet - E1)*(T2 - T)*(eGamma - W1)/volume;
    double w221 = (eJet - E1)*(T - T1)*(W2 - eGamma)/volume;
    double w222 = (eJet - E1)*(T - T1)*(eGamma - W1)/volume;

    double r111 = grid[ijet][iT][iEgamma];
    double r112 = grid[ijet][iT][iEgamma+1];
    double r121 = grid[ijet][iT+1][iEgamma];
    double r122 = grid[ijet][iT+1][iEgamma+1];
    double r211 = grid[ijet+1][iT][iEgamma];
    double r212 = grid[ijet+1][iT][iEgamma+1];
    double r221 = grid[ijet+1][iT+1][iEgamma];
    double r222 = grid[ijet+1][iT+1][iEgamma+1];
    return w111*r111 + w112*r112 + w121*r121 + w122*r122 + w211*r211 + w212*r212 + w221*r221 + w222*r222;

}

ConversionAngularGrid::ConversionAngularGrid(std::string grid_info_file, std::string grid_data_file)
{
    /* 
    * Constructor for the angular grid class.
    * read the gridinfo file to get information on the grid
    * use that information to read in the gridFile(binary) 
    */ 
    std::ifstream settings (grid_info_file, std::ios::in);
    double temp;
    std::string label;
    while (!settings.eof() && settings.is_open())
    {
      settings >> label >> temp;
      if (label == "jet_momentum_grid_size")
        jet_momentum_grid_size = int(temp);
      else if (label == "temperature_grid_size")
        temperature_grid_size = int(temp);
      else if (label == "transverse_k_grid_size")
        transverse_k_grid_size = int(temp);
      else if (label == "longitudinal_k_grid_size")
        longitudinal_k_grid_size = int(temp);
      else if (label == "jet_mom_min")
        jet_mom_min = double(temp);
      else if (label == "jet_mom_max")
        jet_mom_max = double(temp);    
      else if (label == "temperature_min")
        temperature_min = double(temp);
      else if (label == "temperature_max")
        temperature_max = double(temp);
      else if (label == "kT_min")
        kT_min = double(temp);
      else if (label == "kT_max")
        kT_max = double(temp);
      else if (label == "kZ_min")
        kZ_min = double(temp);
      else if (label == "kZ_max")
        kZ_max = double(temp);
      else
      {
        // a tag that doesn't match the above is uselse. Next!
        continue;
      }
    }
    settings.close();
    dPjet = 0.5;
    dT = (temperature_max - temperature_min)/(temperature_grid_size-1); // because I made a mistake in submitting jobs. table has
                                                                        // 50 for grid size but I submitted 51 jobs for temperature
                                                                        // so have to correct for it.
    dkT = (kT_max - kT_min)/transverse_k_grid_size;
    dkZ = (kZ_max - kZ_min)/longitudinal_k_grid_size;
    std::cout<<"jet_momentum_grid_size: "   <<jet_momentum_grid_size;
    std::cout<<", temperature_grid_size: "  <<temperature_grid_size;
    std::cout<<", transverse_k_grid_size: "  <<transverse_k_grid_size;
    std::cout<<", longitudinal_k_grid_size: "<<longitudinal_k_grid_size<<std::endl;
    // now read file and populate the grid
    FILE* fin = fopen(grid_data_file.c_str(), "r+b");
    //std::cout<<fin<<std::endl;
    //std::cout<<(fin==NULL)<<std::endl;
    if (fin == NULL) 
    {
        std::cout<<"\x1B[31m"<<"Failed to open file. Abort.";
        std::cout<<"Initialize the grid to zeros and return."<<"\033[0m"<<std::endl;
        grid = new double***[jet_momentum_grid_size];
        for(int i = 0; i < jet_momentum_grid_size; i++)
        {
            grid[i] = new double**[temperature_grid_size];
            for (int j = 0;  j < temperature_grid_size; j++)
            {
                grid[i][j] = new double*[transverse_k_grid_size];
                for (int k = 0; k < transverse_k_grid_size; k++)
                {
                    grid[i][j][k] = new double[longitudinal_k_grid_size];
                    for (int w = 0; w < longitudinal_k_grid_size; w++)
                    {
                        grid[i][j][k][w] = 0;
                    }
                }
            }
        }
    }
    else
    {
        int status;
        int size = sizeof(double);
        grid = new double***[jet_momentum_grid_size];
        std::cout.precision(5);
        for(int i = 0; i < jet_momentum_grid_size; i++)
        {
            grid[i] = new double**[temperature_grid_size];
            for (int j = 0;  j < temperature_grid_size; j++)
            {
                grid[i][j] = new double*[transverse_k_grid_size];
                for (int k = 0; k < transverse_k_grid_size; k++)
                {
                    grid[i][j][k] = new double[longitudinal_k_grid_size];
                    for(int w = 0; w < longitudinal_k_grid_size; w++)
                    { 
                        double rate;
                        status = fread(&rate, size, 1, fin);
                        if (status)
                        {
                            //std::cout<<"i: "<<i<<" j: "<< j <<" k: "<< k<<" ==> rate: "<<std::scientific<<rate<<endl;
                            grid[i][j][k][w] = rate;
                        }
                    }
                }
            }
        }
        fclose(fin);
    }
    //TEST: print what the maximum values of (pjet, T, kT, kZ) are
    // then let's pick a triplet (pjet, T, kT) and see if it makes 
    // sense
    std::cout<<"pjet max is "<<jet_mom_min + dPjet*jet_momentum_grid_size<<endl;
    std::cout<<"T_max is "<<temperature_min + dT*temperature_grid_size<<endl;
    std::cout<<"kT_max is "<<kT_min + transverse_k_grid_size*dkT<<endl;
    std::cout<<"kZ_max is "<<kZ_min + longitudinal_k_grid_size*dkZ<<endl;
    //now let's pick some indices: ipjet = 3, iT = 20, ikT = 40
    // these should correspond to (3.5, 0.36, 1) GeV values for the triplet
    //ofstream f("test.dat", std::fstream::out);
    //double test_kz = kZ_min;
    //std::cout<<"begin testing"<<endl;
    //for (int i = 0; i < longitudinal_k_grid_size; i++)
    //    f<<test_kz + i*dkZ<<","<<grid[3][20][40][i]<<endl;
    //f.close();
    //exit(0);
}
ConversionAngularGrid::~ConversionAngularGrid()
{
    for (int i = 0; i < jet_momentum_grid_size; i++)
    {
        for (int j = 0; j < temperature_grid_size; j++)
        {
            for (int k = 0; k < transverse_k_grid_size; k++)
                delete[] grid[i][j][k];
            delete[] grid[i][j];
        }
        delete[] grid[i];
    }
}
double ConversionAngularGrid::interpolate(double jet_momentum, double T, double kT, double kZ)
{
    //use the rate table and linear interpolation to get the rate of jet-photon
    //conversion at the given quadruplet (Pjet, T, kT, kZ)
    // maximum and minimum values are
    if (jet_momentum > jet_mom_max || jet_momentum < jet_mom_min ||
        T < temperature_min || T > temperature_max ||
        kT > kT_max || kT < kT_min || kZ > kZ_max || kZ < kZ_min)
        return 0;
    else
    {
        //we're within the limits of the grid.
        // get an estimate on where my values would be in the grid
        double dPjet = 0.5; // TODO: this one I hardcode for now
                            // currently I only have 16 of the 40 
                            // jet energies calculated.
        double dT = (temperature_max - temperature_min)/(temperature_grid_size-1);
        double dkT = (kT_max - kT_min)/transverse_k_grid_size;
        double dkZ = (kZ_max - kZ_min)/longitudinal_k_grid_size; 
        int ijet = floor((jet_momentum - jet_mom_min)/dPjet);
        int iT   = floor((T - temperature_min)/dT);
        int ikT  = floor(kT/dkT);
        int ikZ  = floor((kZ - kZ_min)/dkZ);
        if (ijet + 2 > 20 || ijet-1 < 0 || 
            iT + 2 > temperature_grid_size || iT - 1 < 0 ||
            ikT + 2 > transverse_k_grid_size || ikT - 1 < 0 ||
            ikZ + 2 > longitudinal_k_grid_size || ikZ - 1 < 0)
        {
            //std::cout<<"index out of bounds: "<<" iJet: "<<ijet<<", iT: "<<iT<<", ikT: "<<ikT<<", ikZ: "<<ikZ;
            //std::cout<<"pJet: "<<jet_momentum<<"kT: "<<kT<<"kZ: "<<kZ<<"T: "<<T<<" return 0."<<std::endl;
            return 0;
        }
        else
        {
            // std::cout<<" iJet: "<<ijet<<", iT: "<<iT<<", ikT: "<<ikT<<", ikZ: "<<ikZ<<endl;
            // Now generate the (E, T, kT, kZ) values that braket the provided quadruplet
            double p0 = jet_mom_min + (ijet-1)*dPjet;
            double p1 = jet_mom_min + (ijet+1)*dPjet;
            double T0 = temperature_min + (iT-1)*dT;
            double T1 = temperature_min + (iT+1)*dT;
            double kT0 = (ikT-1)*dkT;
            double kT1 = (ikT+1)*dkT;
            double kZ0 = kZ_min + (ikZ-1)*dkZ;
            double kZ1 = kZ_min + (ikZ+1)*dkZ;
            //std::cout<<"p0 = "<<p0<<", pjet = "<<jet_momentum<<", p1 = "<<p1<<std::endl;
            //std::cout<<"T0 = "<<T0<<", T = "<<T<<", T1 = "<<T1<<std::endl;
            //std::cout<<"kT0 = "<<kT0<<", kT = "<<kT<<", kT1 = "<<kT1<<std::endl;
            //std::cout<<"kZ0 = "<<kZ0<<", kz = "<<kZ<<", kZ1 = "<<kZ1<<std::endl;
            //now form the differences
            double p1p0 = (p1 - p0);
            double p1pJ = (p1 - jet_momentum);
            double pJp0 = (jet_momentum - p0);
            double T1T0 = (T1 - T0);
            double T1TJ = (T1 - T);
            double TJT0 = (T - T0);
            double kT1kT0 = (kT1 - kT0);
            double kT1kTJ = (kT1 - kT);
            double kTJkT0 = (kT - kT0);
            double kZ1kZ0 = (kZ1 - kZ0);
            double kZ1kZJ = (kZ1 - kZ);
            double kZJkZ0 = (kZ - kZ0);
            //denominator of the weights is shared among all of them:
            double denom = p1p0*T1T0*kT1kT0*kZ1kZ0;
            //weights multiplied by the relevant rate value:
            double N0 = T1TJ*kT1kTJ*kZ1kZJ*grid[ijet-1][iT-1][ikT-1][ikZ-1];
            double N1 = T1TJ*kT1kTJ*kZJkZ0*grid[ijet-1][iT-1][ikT-1][ikZ+1];
            double N2 = T1TJ*kTJkT0*kZJkZ0*grid[ijet-1][iT-1][ikT+1][ikZ+1];
            double N3 = TJT0*kTJkT0*kZJkZ0*grid[ijet-1][iT+1][ikT+1][ikZ+1];
            double N4 = TJT0*kT1kTJ*kZJkZ0*grid[ijet-1][iT+1][ikT-1][ikZ+1];
            double N5 = TJT0*kTJkT0*kZ1kZJ*grid[ijet-1][iT+1][ikT+1][ikZ-1];
            double N6 = TJT0*kT1kTJ*kZ1kZJ*grid[ijet-1][iT+1][ikT-1][ikZ-1];
            double N7 = T1TJ*kTJkT0*kZ1kZJ*grid[ijet-1][iT-1][ikT+1][ikZ-1];
            double N8 = T1TJ*kT1kTJ*kZ1kZJ*grid[ijet+1][iT-1][ikT-1][ikZ-1];
            double N9 = T1TJ*kT1kTJ*kZJkZ0*grid[ijet+1][iT-1][ikT-1][ikZ-1];
            double N10 = T1TJ*kTJkT0*kZJkZ0*grid[ijet+1][iT-1][ikT+1][ikZ+1];
            double N11 = TJT0*kTJkT0*kZJkZ0*grid[ijet+1][iT+1][ikT+1][ikZ+1];
            double N12 = TJT0*kT1kTJ*kZ1kZJ*grid[ijet+1][iT+1][ikT-1][ikZ-1];
            double N13 = TJT0*kTJkT0*kZ1kZJ*grid[ijet+1][iT+1][ikT+1][ikZ-1];
            double N14 = TJT0*kT1kTJ*kZJkZ0*grid[ijet+1][iT+1][ikT-1][ikZ+1];
            double N15 = T1TJ*kTJkT0*kZ1kZJ*grid[ijet+1][iT-1][ikT+1][ikZ-1];
            double rate = p1pJ*(N0+N1+N2+N3+N4+N5+N6+N7) + pJp0*(N8+N9+N10+N11+N12+N13+N14+N15);
            return rate/denom;
        }
    }
    return -1; // if you're here something has gone wrong that I didn't catch. Return a nonsensical number
}
/*
    double ConversionAngularGrid::integrand_kz(double jetx, double jety, double jetx, double kT, double phiGamma, double v1, double v2, double v3, double T)
        (jetx, jety, jetz) : jet momentum components
        (kT, phiGamma) : photon transverse momentum component and its azimuthal angle in the lab frame
        (vx, vy, vz) : 
*/
double ConversionAngularGrid::integrand_kz(double kz, double phiGamma, double jetx, double jety, double jetz, double kT, double vx, double vy, double vz, double T)
{
    double jet_momentum = sqrt(jetx*jetx + jety*jety + jetz*jetz);
    double cosTheta = jetz/jet_momentum;
    double cosPhi = jetx/sqrt(jetx*jetx + jety*jety);
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double sinPhi = sqrt(1. - cosPhi*cosPhi);
    double v2 = vx*vx + vy*vy + vz*vz;
    double gamma = 1./sqrt(1. - v2);
    double k, kx, ky;  // photon in the lab frame
    double kxPrime, kyPrime, kzPrime; // boosted frame (cell rest frame) photon components
    double qx, qy, qz; //photon momentum components in the boosted+rotated frame 
    kx = kT*cos(phiGamma);
    ky = kT*sin(phiGamma);

    k = sqrt(kx*kx + ky*ky + kz*kz);
    kxPrime = kx*(1. + (gamma-1.)*vx*vx/v2) + vx*(ky*vy + kz*vz)*(gamma-1.)/v2 - gamma*vx*k;
    kyPrime = ky*(1. + (gamma-1.)*vy*vy/v2) + vy*(kx*vx + kz*vz)*(gamma-1.)/v2 - gamma*vy*k;
    kzPrime = kz*(1. + (gamma-1.)*vz*vz/v2) + vz*(kx*vx + ky*vy)*(gamma-1.)/v2 - gamma*vz*k;
    // rotate the momenta
    qx = kxPrime*cosTheta*cosPhi + kyPrime*cosTheta*sinPhi -kzPrime*sinTheta;
    qy = -kxPrime*sinPhi + kyPrime*cosPhi;
    qz = kxPrime*cosPhi*sinTheta + kyPrime*sinTheta*sinPhi + kzPrime*cosTheta;
    // now calculate the rate at this combination. The return values contains the boost factor to
    // bring the calculation back to the lab frame as well as the kT factor from the jacobian of 
    // the spherical coordinates.
    //double boostFactor = gamma*(1.-(kx*vx-ky*vy-kz*vz)/(k*sqrt(v2)));
    // I don't boost because the integral  measure is in fact over the lab frame components
    // not cell frame, and I already boost this in MARTINI proper by gamma*(1 - u dot p (jet))
    // I have to divide by a 2 pi though, to account for the dphi in the cell frame that I trivially
    // integrated over.
    return kT*interpolate(jet_momentum, T, sqrt(qx*qx + qy*qy), qz)/(2*PI);
}

double ConversionAngularGrid::integrate_kz(double phiGamma, double eta, double jetx, double jety, double jetz, double kT, double vx, double vy, double vz, double T)
{
    double w[] = {5./9, 8./9, 5./9};
    double x[] = {-sqrt(3./5), 0, sqrt(3./5)};
    double result=0;
    double kZMin = -kT*sinh(eta);
    double kZMax = kT*sinh(eta);
    double deltaKz = (kZMax - kZMin)/10;
    double f1, f2, f3;
    double z1, z2, z3;
    double a = 0, b = 0;
    for(int i = 0; i < 10; i++)
    {
        a = b;
        b += deltaKz;
        z1 = (b-a)*x[0]/2 + (a+b)/2;
        z2 = (b-a)*x[1]/2 + (a+b)/2;
        z3 = (b-a)*x[2]/2 + (a+b)/2;
        f1 = integrand_kz(z1, phiGamma, jetx, jety, jetz, kT, vx, vy, vz, T);
        f2 = integrand_kz(z2, phiGamma, jetx, jety, jetz, kT, vx, vy, vz, T);
        f3 = integrand_kz(z3, phiGamma, jetx, jety, jetz, kT, vx, vy, vz, T);
        result += (w[0]*f1 + w[1]*f2 + w[2]*f3)*(b-a)/2;
    }
    return result;
}

double ConversionAngularGrid::convert(double jetx, double jety, double jetz, double kT, double eta, double vx, double vy, double vz, double T)
{
    /*
        convert(double jetx, double jety, double jetz, double kT, double eta, double v1, double v2, double v3, double T):
            -jetx, jety, jetz = jet momentum vector in the cell rest frame, needed for the rotation
            -kT, eta = photon transverse momentum and pseudorapidity
            -vx, vy, vz, T = cell velocity vector, cell temperature
        -Integrate over the lab frame phi_gamma and kz. 
        -Integrals will be Gaussian Quadrature rules. The kz iterval will be divided into 10 subintervals
         and a Gauss-Legendre 3 point rule would be implemented. The phi integral will also be 10 sub intervals
         with a 3 point Gauss-Legendre rule.
    */
    //weights and points of n=3 Gauss-Legendre quadrature rules
    double w[] = {5./9, 8./9, 5./9};
    double x[] = {-sqrt(3./5), 0, sqrt(3./5)};
    double f1, f2, f3;
    double z1, z2, z3;
    double a = 0, b = 0;
    double deltaPhi = (2*PI)/10;
    double result=0;
    for (int i = 0; i < 10; i++)
    {
        //compute the current limits of the phi integral: a_i, b_i
        a = b;
        b = a + deltaPhi;
        z1 = (b-a)*x[0]/2 + (a+b)/2;
        z2 = (b-a)*x[1]/2 + (a+b)/2;
        z3 = (b-a)*x[2]/2 + (a+b)/2;
        f1 = integrate_kz(z1, eta, jetx, jety, jetz, kT, vx, vy, vz, T);
        f2 = integrate_kz(z2, eta, jetx, jety, jetz, kT, vx, vy, vz, T);
        f3 = integrate_kz(z3, eta, jetx, jety, jetz, kT, vx, vy, vz, T); 
        result += (f1*w[0] + f2*w[1] + f3*w[2])*(b-a)/2;
    } 
    return result;
}
