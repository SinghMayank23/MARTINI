#ifndef Background_H
#define Background_H

#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace std;

class Background
{
  private:
    const int n_grid_eta = 40;
    const int n_grid_phi = 72;

    const int et_min_jet = 30.;


    int num_CF_events;
    vector< vector<PseudoJet> > background;

    bool output;
    double pt_min_trk;
    double et_min_trk;
    double eta_max_trk;
    double r_size;

    string subrun_id;

    vector< vector<double> > et_cell;
    vector< vector<double> > et_cell_subtracted;
    vector< vector<bool> > inside_jet;
    vector<double> eta_grid, phi_grid;
    vector<double> et_cell_avg, et_cell_avgsq, et_cell_std;


  public:
    Background(double _pt_min_trk, double _et_min_trk, double _eta_max_trk,
               double _r_size, string _subrun_id, bool _output);
    ~Background() {};

    ClusterSequence noise_pedestal_subtraction
        (vector<PseudoJet> hard_particles, vector<PseudoJet> background, double kappa);
    ClusterSequence noise_pedestal_subtraction
        (vector<PseudoJet> hard_particles, vector<PseudoJet> recoils,
         vector<PseudoJet> negative_particles,vector<PseudoJet> background,
         double kappa, bool recoil = false);
    void compute_et_cell(vector<PseudoJet> particles, bool negative = false);
};

#endif
