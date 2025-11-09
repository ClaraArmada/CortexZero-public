#include "Dendrite.h"

float Dendrite::modulation() {
    float modulator_effect = 0.0;
    for (size_t i = 0; i < neuromodulators.size(); ++i) {
        if (neuromodulators[i] >= neuromod_th) {
            modulator_effect += 1.0 - std::exp(-neuromod_mult[i] * neuromodulators[i]);
            modulator_effect = std::clamp(modulator_effect, 0.0f, 5.0f);
        }
    }
    if (C < C_low) {
        return 0.0;
    }
    else if (C > C_high) {
        float excess = (C - C_high) / C_high;
        return std::max(0.0f, 1 - excess * overactivity_penalty) * modulator_effect;
    }
    else {
        float norm_C = (C - C_low) / (C_high - C_low);
        return norm_C;
    }
}

void Dendrite::clampWeights()
{
    for (float& weight : w) {
        if (isInhibitory) {
            weight = std::clamp(weight, -10.0f, 0.0f);
        }
        else {
            weight = std::clamp(weight, 0.0f, 10.0f);
        }
    }
}

std::vector<float> Dendrite::homeostasis() {
    if (_homeo_counter >= homeo_interval) {
        _homeo_counter = 0;
        float eps = 1e-6;
        float scale = std::clamp(target_rate / (r + eps), 0.9f, 1.1f);
        for (int weight = 0; weight < w.size(); weight++) {
            w[weight] *= scale;
        }
    }
    clampWeights();
    return w;
}

void Dendrite::updateWeights(const int t)
{
    if (S_Dt == 1) {
        float mod = modulation();
        for (int i = 0; i < w.size(); i++) {
            float weightChange = 0.0;
            float delta_t = t - lastSpikeTimes[i];

            // Calcium gating function

            float f_Ca = 0;
            if (C >= C_high) {
                f_Ca = +((C - C_high) / (C_max - C_high));
            }
            else if (C <= C_low) {
                f_Ca = -((C_low - C) / std::max(C_low, 1e-6f));
            }

            if (delta_t > 0) {
                int pre_spike = 1;
                int post_spike = 1;
                switch (learningRule) {
                case e_learningrule::stdp:
                    if (delta_t > 0) {
                        weightChange = learning_rate * A_plus * std::exp(-delta_t / tau);
                    }
                    else {
                        weightChange = learning_rate * -A_minus * std::exp(delta_t / tau);
                    }
                    break;

                case e_learningrule::voltage_dependent: {
                    float V_pre;
                    if (delta_t < 5) {
                        V_pre = V_th_dendrite;
                    }
                    else {
                        V_pre = V_rest;
                    }

                    if (V > V_th_dendrite) {
                        weightChange = learning_rate * (V - V_th_dendrite) / 20.0;
                    }
                    else {
                        weightChange = learning_rate * -0.1;
                    }
                    break;
                }

                case e_learningrule::hebbian:
                    if (w[i] > 0) {
                        weightChange = learning_rate * 0.012 * pre_spike * post_spike;
                    }
                    break;

                case e_learningrule::antihebbian:
                    weightChange = learning_rate * 0.01 * -pre_spike * post_spike;
                    break;

                case e_learningrule::covariance: {
                    float mean_pre = 0.05;
                    float mean_post = r;
                    weightChange = learning_rate * (pre_spike - mean_pre) * (post_spike - mean_post);
                    break;
                }

                default:
                    std::cerr << "Invalid rule\n";
                    return;
                }
                w[i] += weightChange * mod * NMDA_mult * (1 + beta_C * f_Ca);
                Dendrite::clampWeights();
            }
        }
    }
}

void Dendrite::step(const int t, const std::vector<int> S, const double delta_t, const std::vector<float> mod)
{
    // last spike time update
    for (size_t i = 0; i < S.size(); ++i) {
        if (S[i] == 1) {
            lastSpikeTimes[i] = t;
        }
    }

    // neuromodulator add and decay
    size_t N = std::min(neuromodulators.size(), mod.size());
    for (size_t i = 0; i < N; i++) {
        neuromodulators[i] += mod[i]; // add incoming strength
        neuromodulators[i] = neuromodulators[i] * std::exp(-delta_t / neuromod_decay)
            + neuromod_rest * (1.0f - std::exp(-delta_t / neuromod_decay));
    }

    // sum of spike*weight
    float It_tot = 0.0f; // total sum

    for (size_t i = 0; i < w.size() && i < S.size(); ++i) {
        It_tot += w[i] * static_cast<float>(S[i]);
    }

    float It_exc = std::max(It_tot, 0.0f);

    // NMDA simulation
    float H_M;
    float H_I;
    if (V - V_depol > 0) {
        H_M = V - V_depol;
    }
    else {
        H_M = 0;
    }

    if (It_exc - NMDA_input_th > 0) {
        H_I = It_exc - NMDA_input_th;
    }
    else {
        H_I = 0;
    }
    float trigger = H_M * H_I;

    NMDA_act = NMDA_act * std::exp(-delta_t / T_NMDA) + A_NMDA * trigger;

    NMDA_mult = 1 + k * NMDA_act;

    NMDA_mult = std::min(NMDA_mult, NMDA_max);

    // calcium simulation (with NMDA influence)
    float depol = std::max(0.0f, V - V_rest);

    // NMDA-driven calcium influx (step 1 multiplier already computed elsewhere)
    float dCa_NMDA = alpha_NMDA * (NMDA_mult - 1.0f);

    // total calcium update
    if (S_Dt == 1) {
        C = alpha_C * C + V * depol + beta_S + dCa_NMDA;
        r = (1.0f - alpha_r) * r + alpha_r * 1.0f;
    }
    else {
        C = alpha_C * C + V * depol + dCa_NMDA;
        r = (1.0f - alpha_r) * r + alpha_r * 0.0f;
    }

    C = std::clamp(C, 0.0f, C_max);

    // I_eff
    float It_eff = It_tot * (NMDA_mult * It_eff_scale);

    //refractory period
    if (refractory) {
        V = V_reset;
        ref += 1;
        if (ref == t_ref) {
            refractory = false;
            ref = 0;
        }
    }
    // membrane potential
    else {
        float V_new = V + (delta_t / T_m) * (-(V - V_rest) + (R * It_eff));
        V_new = std::clamp(V_new, -85.0f, -30.0f);

        if (V_new >= V_th_dendrite) {
            S_Dt = 1;
            V = V_reset;
            refractory = true;
        }
        else {
            S_Dt = 0;
            V = V_new;
        }
    }

    //V += f_normal_random(0.0, 0.05);

    // update weights
    updateWeights(t);

    // L2 regularization
    if (S_count % 50 == 0) {
        float activity_error = r - target_rate;
        float decay_strength = 0.0f;

        if (activity_error > 0.005f) {
            decay_strength = 0.0005f * activity_error / target_rate;
        }
        else if (activity_error < -0.005f) {
            decay_strength = -0.0003f;
        }

        for (int i = 0; i < w.size(); i++) {
            if (w[i] > 0.1f) {
                w[i] -= decay_strength * w[i];
                w[i] = std::max(0.1f, w[i]);
            }
        }
    }

    // homeostasis
    _homeo_counter += 1;
    homeostasis();
}

bool Dendrite::isDistalNode() const {
    return isDistal;
}