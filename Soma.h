#pragma once

#include <vector>
#include <cmath>
#include <iostream>

#include "../include/Variables.h"
#include "../include/utils.h"

class Dendrite;

class Soma {
    std::vector<float> w; // weights
    float V = V_rest; // membrane potential, Na^+ in Volt
    float V_Dt = 0.0; // membrane potential
    float C = 0.0; // calcium, Ca^2+

    float r = 0.0;
    float NMDA_mult = 1.0; // councidence detection, and if NMDA is frequently high, pruning is less likely
    float NMDA_act = 0.0;

    int S_Dt = 0;
    int S_count = 0;

    std::vector<float> lastSpikeTimes;
    e_learningrule learningRule;
    std::vector<float> neuromodulators; // stores the active neuromodulators, and their strength
    std::vector<float> neuromod_mult; // predefined multiplier for neuromodulators

    int _homeo_counter = 0;

    bool refractory = false;
    int t_ref = refractory_length;
    int ref = 0;
    bool isInhibitory; // if it's inhibitory or excitatory

    //temp
    float It_tot_log;
public:
    int S = 0; // if the soma sent a signal or not
    explicit Soma(const std::vector<float>& weights, bool inhibitory, e_learningrule learningrule, std::vector<float> neuromodmultiplier)
        : w(weights), isInhibitory(inhibitory), learningRule(learningrule), neuromodulators(neuromodmultiplier.size(), 0.0f), neuromod_mult(neuromodmultiplier) {
        lastSpikeTimes.resize(w.size(), -1000.0f);
    }
    float modulation();
    void clampWeights();
    std::vector<float> homeostasis();
    void updateWeights(const int t);
    void step(const int t, const std::vector<int> S, const double delta_t, const std::vector<float> mod, float V_dend_sum, int num_root_dendrites, float I_dend_spike);
    // helpers
    float getItTot() const { return It_tot_log; }
    float getVoltage() const { return V; }
    float getCalcium() const { return C; }
    int getSpike() const { return S_Dt; }
    float getNMDA() const { return NMDA_mult; }
    float getWeight(int idx) const { return idx < w.size() ? w[idx] : 0; }
};