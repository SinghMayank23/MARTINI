#include <fstream>
#include <iomanip>

#include "Background.h"

using namespace fastjet;
using namespace std;

Background::Background(double _pt_min_trk, double _et_min_trk, double _eta_max_trk,
                       double _r_size, string _subrun_id, bool _output)
{
  // set basic variables
  output = _output;
  pt_min_trk  = _pt_min_trk;
  et_min_trk  = _et_min_trk;
  eta_max_trk = _eta_max_trk;
  r_size = _r_size;
  subrun_id = _subrun_id;

  // initialize cells/grids
  et_cell.resize(n_grid_eta, vector<double>(n_grid_phi));
  et_cell_subtracted.resize(n_grid_eta, vector<double>(n_grid_phi));
  inside_jet.resize(n_grid_eta, vector<bool>(n_grid_phi));

  et_cell_avg.resize(n_grid_eta);
  et_cell_avgsq.resize(n_grid_eta);
  et_cell_std.resize(n_grid_eta);

  for(int i=0; i<n_grid_eta; i++)
  {
    double val_eta = ((double)i+0.5)*(2.*eta_max_trk)/n_grid_eta-eta_max_trk;
    eta_grid.push_back(val_eta);
  }

  for(int j=0; j<n_grid_phi; j++)
  {
    double val_phi = ((double)j+0.5)*(2.*M_PI)/n_grid_phi;
    phi_grid.push_back(val_phi);
  }

}

ClusterSequence Background::noise_pedestal_subtraction
    (vector<PseudoJet> hard_particles, vector<PseudoJet> background,
     double kappa)
{
  vector<PseudoJet> recoils, negative_particles;
  return noise_pedestal_subtraction
    (hard_particles, recoils, negative_particles, background, kappa);
}

ClusterSequence Background::noise_pedestal_subtraction
    (vector<PseudoJet> hard_particles, vector<PseudoJet> recoils,
     vector<PseudoJet> negative_particles, vector<PseudoJet> background,
     double kappa, bool recoil)
{
  // assign 0 to every cells
  for(auto &phi : et_cell) fill(phi.begin(), phi.end(), 0.);
  for(auto &phi : et_cell_subtracted) fill(phi.begin(), phi.end(), 0.);
  for(auto &phi : inside_jet) fill(phi.begin(), phi.end(), false);

  fill(et_cell_avg.begin(), et_cell_avg.end(), 0.);
  fill(et_cell_avgsq.begin(), et_cell_avgsq.end(), 0.);
  fill(et_cell_std.begin(), et_cell_std.end(), 0.);

  // compute ET at given cell
  compute_et_cell(hard_particles);
  compute_et_cell(background);
  if(recoil)
  {
    compute_et_cell(recoils);
    compute_et_cell(negative_particles, true);
  }

  //cout << "recoil:" << recoil << endl;
  //cout << "size of hard particles:" << hard_particles.size() << endl;
  //cout << "size of recoils:" << recoils.size() << endl;
  //cout << "size of negative particles:" << negative_particles.size() << endl;
  //cout << "size of background:" << background.size() << endl;

  for(int i=0; i<n_grid_eta; i++)
  {
    // compute average and standard deviation of ET at given rapidity ring
    for(int j=0; j<n_grid_phi; j++)
    {
      et_cell_avg[i] += et_cell[i][j]/n_grid_phi;
      et_cell_avgsq[i] += et_cell[i][j]*et_cell[i][j]/n_grid_phi;
    }
    et_cell_std[i] = sqrt(et_cell_avgsq[i] - et_cell_avg[i]*et_cell_avg[i]);

    // subtract background density from each cell
    for(int j=0; j<n_grid_phi; j++)
    {
      double subtracted_et = et_cell[i][j] - et_cell_avg[i] - kappa*et_cell_std[i];
      et_cell_subtracted[i][j] = (subtracted_et > 0.) ? subtracted_et : 0.;
    }
  }
  //cout << endl << "first run" << endl;
  //for(int i=0; i<n_grid_eta; i++)
  //  cout << eta_grid[i] << " " << et_cell_avg[i] << " " << et_cell_std[i] << endl;
  //cout << endl;

  vector<PseudoJet> temp_particles;
  for(int i=0; i<n_grid_eta; i++)
  {
    for(int j=0; j<n_grid_phi; j++)
    {
      if(et_cell_subtracted[i][j] < 1e-5) continue;

      double px = et_cell_subtracted[i][j]*cos(phi_grid[j]);
      double py = et_cell_subtracted[i][j]*sin(phi_grid[j]);
      double pz = et_cell_subtracted[i][j]*sinh(eta_grid[i]);
      double E  = et_cell_subtracted[i][j]*cosh(eta_grid[i]);
      temp_particles.push_back( PseudoJet(px, py, pz, E) );
    }
  }

  JetDefinition jet_def(antikt_algorithm, r_size);

  // run the clustering, extract the jets using temp_particles
  ClusterSequence temp_cs(temp_particles, jet_def);
  vector<PseudoJet> temp_jets = sorted_by_pt(temp_cs.inclusive_jets(et_min_jet));

  vector<PseudoJet> subtracted_jets;

  if(temp_jets.size() == 0)
  // if there are no jets, no further background estimation
  {
    subtracted_jets = sorted_by_pt(temp_cs.inclusive_jets());
    return temp_cs;
  }
  else
  // re-do estimating background without jets whose pt ~ et > et_min_Jet
  {
    // re-assign 0 to every cells
    for(auto &phi : et_cell_subtracted) fill(phi.begin(), phi.end(), 0);
    fill(et_cell_avg.begin(), et_cell_avg.end(), 0.);
    fill(et_cell_avgsq.begin(), et_cell_avgsq.end(), 0.);
    fill(et_cell_std.begin(), et_cell_std.end(), 0.);

    for(int i=0; i<n_grid_eta; i++)
    {
      int n_cell_phi = 0;

      et_cell_avg[i] = 0.;
      et_cell_avgsq[i] = 0.;
      et_cell_std[i] = 0.;

      for(int j=0; j<n_grid_phi; j++)
      {
        // iterate through all jets
        for(auto jet : temp_jets)
        {
          double delta_eta = eta_grid[i] - jet.eta();
          double delta_phi = phi_grid[j] - jet.phi();
          if(delta_phi > M_PI) delta_phi -= 2.*M_PI;
          else if (delta_phi < -M_PI) delta_phi += 2.*M_PI;
          double delta_r = sqrt(pow(delta_eta, 2.) + pow(delta_phi, 2.));

          //cout << i << " eta:" << eta_grid[i] << " etaJ:" << jet.eta() <<
          //             " dEta:" << delta_eta << " "
          //     << j << " phi:" << phi_grid[j] << " phiJ:" << jet.phi() <<
          //             " dphi:" << delta_phi << " "
          //     << "delta_r:" << delta_r << endl;

          // exclude cells within constructed jets
          if(delta_r < r_size) 
            inside_jet[i][j] = true;
        }

        // estimate avg and std of background excluding the cells incide jets
        if(!inside_jet[i][j])
        {
          et_cell_avg[i] += et_cell[i][j];
          et_cell_avgsq[i] += pow(et_cell[i][j], 2.);
          n_cell_phi++;
        }
      }

      //cout << "eta:"<< eta_grid[i] << " size of phi cell:" << n_cell_phi << endl;
      et_cell_avg[i] /= n_cell_phi;
      et_cell_avgsq[i] /= n_cell_phi;
      et_cell_std[i] = sqrt(et_cell_avgsq[i] - et_cell_avg[i]*et_cell_avg[i]);

      // subtract background density from each cell
      for(int j=0; j<n_grid_phi; j++)
      {
        double subtracted_et = et_cell[i][j] - et_cell_avg[i] - kappa*et_cell_std[i];
        et_cell_subtracted[i][j] = (subtracted_et > 0.) ? subtracted_et : 0.;
        //et_cell_subtracted[i][j] = subtracted_et;
      }
    }
    //cout << endl << "second run" << endl;
    //for(int i=0; i<n_grid_eta; i++)
    //  cout << eta_grid[i] << " " << et_cell_avg[i] << " " << et_cell_std[i] << endl;
    //cout << endl;

    // cluster jets again
    temp_particles.clear();
    for(int i=0; i<n_grid_eta; i++)
    {
      for(int j=0; j<n_grid_phi; j++)
      {
        if(et_cell_subtracted[i][j] < 1e-5) continue;

        double px = et_cell_subtracted[i][j]*cos(phi_grid[j]);
        double py = et_cell_subtracted[i][j]*sin(phi_grid[j]);
        double pz = et_cell_subtracted[i][j]*sinh(eta_grid[i]);
        double E  = et_cell_subtracted[i][j]*cosh(eta_grid[i]);
        temp_particles.push_back( PseudoJet(px, py, pz, E) );
      }
    }
    
    // run the clustering, extract the jets using temp_particles
    ClusterSequence cs(temp_particles, jet_def);
    subtracted_jets = sorted_by_pt(cs.inclusive_jets());

    if(output)
    {
      vector<PseudoJet> subtracted_jets2 = sorted_by_pt(cs.inclusive_jets(20.));

      // original jets
      ClusterSequence cs_temp(hard_particles, jet_def);
      vector<PseudoJet> original_jets = sorted_by_pt(cs_temp.inclusive_jets(20.));

      string file_name = "list_JES_" + subrun_id + ".dat";
      ofstream fout(file_name, ios::app);

      for(auto j_subt : subtracted_jets2)
      {
        int idx = -1;
        for(unsigned int j=0; j<original_jets.size(); j++)
        {
          double delta_phi = j_subt.phi() - original_jets[j].phi();
          if(delta_phi > M_PI) delta_phi -= 2.*M_PI;
          else if (delta_phi < -M_PI) delta_phi += 2.*M_PI;
          double delta_eta = j_subt.eta() - original_jets[j].eta();
          double delta_r = sqrt(pow(delta_phi, 2.) + pow(delta_eta, 2.));
          if(delta_r < r_size*1.5)
          {
            idx = j;
            break;
          }
        }
        if(idx == -1) continue;

        double pt_orig = 0.;
        double pt_bg = 0.;

        // ratio between original jet pt and bkgd-subtracted jet pt
        double JES = original_jets[idx].pt()/j_subt.pt();

        for(auto j : hard_particles)
        {
          double d_phi = j.phi() - j_subt.phi();
          double d_eta = j.eta() - j_subt.eta();
          if(d_phi > M_PI) d_phi -= 2.*M_PI;
          else if (d_phi < -M_PI) d_phi += 2.*M_PI;
          double delta_r = sqrt(pow(d_eta, 2.) + pow(d_phi, 2.));
          if(delta_r < r_size)
            pt_orig += j.pt();
        }
        for(auto j : background)
        {
          double d_phi = j.phi() - j_subt.phi();
          double d_eta = j.eta() - j_subt.eta();
          if(d_phi > M_PI) d_phi -= 2.*M_PI;
          else if (d_phi < -M_PI) d_phi += 2.*M_PI;
          double delta_r = sqrt(pow(d_eta, 2.) + pow(d_phi, 2.));
          if(delta_r < r_size)
            pt_bg += j.pt();
        }
        double pt_sum = pt_orig + pt_bg;
        double subtracted = pt_sum - j_subt.pt();
        fout << j_subt.pt() << " " << JES << " " << subtracted << endl;
      }
      fout.close();
    }
    return cs;
  }
}

void Background::compute_et_cell(vector<PseudoJet> particles, bool negative)
{
  for(auto p : particles)
  {
    double eta = p.eta();
    double phi = p.phi();
    double et = p.e()*p.pt()/sqrt(p.pt()*p.pt()+p.pz()*p.pz());

    int i_eta = (int)((eta+eta_max_trk)/(2.*eta_max_trk)*n_grid_eta);
    int i_phi = (int)(phi/(2.*M_PI)*n_grid_phi);
    if(i_eta >= 0 && i_eta < n_grid_eta && i_phi >= 0 && i_phi < n_grid_phi)
    {
      if(negative) et_cell[i_eta][i_phi] -= et;
      else         et_cell[i_eta][i_phi] += et;
    }
  }
}
