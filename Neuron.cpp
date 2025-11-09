#include "Neuron.h"
#include "Soma.h"
#include "Dendrite.h"

Neuron::Neuron(std::vector<float> neuromodmult, const std::vector<float> weightRange, const int maxDendriticDepth, const int dendriticTreeCount, const int branchingwidth, bool inhibitory)
{
    std::vector<float> dendLR_P;
    std::vector<float> somaLR_P;
    if (inhibitory) dendLR_P = dendrititicRule_I_P; else dendLR_P = dendrititicRule_P;
    if (inhibitory) somaLR_P = somaticRule_I_P; else somaLR_P = somaticRule_P;
    float weightLoc = (weightRange[1] + weightRange[0]) / 2.0f;

    for (int root = 0; root < dendriticTreeCount; root++) {
        int treedepth = std::max(1, rand() % maxDendriticDepth + 1);

        // create root dendrite
        Dendrite* rootPtr = new Dendrite(
            f_random_weights(weightRange[0], weightRange[1], branchingwidth),
            inhibitory,
            f_sampleLearningRule(dendLR_P),
            neuromodmult
        );
        rootDendrites.push_back(rootPtr);

        // recursively generate branches
        buildTree(rootPtr, treedepth - 1, weightRange, branchingwidth, inhibitory, neuromodmult);
    }
    collectDistal();

    int total_dendritic_inputs = rootDendrites.size() + distalDendrites.size();

    std::vector<float> soma_weights;
    for (int i = 0; i < total_dendritic_inputs; i++) {
        soma_weights.push_back(5.0f);
    }

    soma = new Soma(
        soma_weights,
        inhibitory,
        f_sampleLearningRule(somaLR_P),
        neuromodmult);
}

void Neuron::buildTree(Dendrite* parent, int remainingDepth,
    const std::vector<float>& weightRange,
    int branchingwidth, bool inhibitory, std::vector<float> neuromodmult)
{
    std::vector<float> dendLR_P;
    if (inhibitory) dendLR_P = dendrititicRule_I_P; else dendLR_P = dendrititicRule_P;

    if (remainingDepth <= 0) return;

    float weightLoc = (weightRange[1] + weightRange[0]) / 2.0f;

    for (int b = 0; b < branchingwidth; b++) {
        Dendrite* child = new Dendrite(
            f_random_weights(weightRange[0], weightRange[1], branchingwidth),
            inhibitory,
            f_sampleLearningRule(dendLR_P),
            neuromodmult
        );

        parent->dendriticNodes.push_back(child);

        buildTree(child, remainingDepth - 1, weightRange, branchingwidth, inhibitory, neuromodmult);
    }
}

void Neuron::traverse(Dendrite* dend)
{
    if (dend->dendriticNodes.empty()) {
        distalDendrites.push_back(dend);
        dend->isDistal = true;
        return;
    }

    for (Dendrite* d : dend->dendriticNodes) {
        traverse(d);
    }
}

void Neuron::collectDistal()
{
    distalDendrites.clear();
    for (Dendrite* root : rootDendrites) {
        traverse(root);
    }
}

//std::vector<float> Neuron::homeostasis() {
//
//}

void Neuron::traverseStep(Dendrite* dend, const int t, const double delta_t, const std::vector<float>& mod)
{
    for (Dendrite* child : dend->dendriticNodes) {
        traverseStep(child, t, delta_t, mod);
    }

    if (dend->dendriticNodes.empty()) {
        return;
    }

    std::vector<int> inputs;
    inputs.reserve(dend->dendriticNodes.size());
    for (Dendrite* child : dend->dendriticNodes) {
        inputs.push_back(child->S);
    }

    dend->step(t, inputs, delta_t, mod);
}

void Neuron::step(const int t,
    const std::vector<std::vector<int>>& S,
    const double delta_t,
    const std::vector<float>& mod)
{
    for (int d = 0; d < distalDendrites.size(); d++) {
        if (d < S.size()) {
            distalDendrites[d]->step(t, S[d], delta_t, mod);
        }
    }

    for (Dendrite* root : rootDendrites) {
        traverseStep(root, t, delta_t, mod);
    }

    std::vector<int> S_soma;
    S_soma.reserve(rootDendrites.size() + distalDendrites.size());

    for (Dendrite* root : rootDendrites)
        S_soma.push_back(root->S);

    for (Dendrite* distal : distalDendrites)
        S_soma.push_back(distal->S);

    float I_dend_spike = 0.0f;
    for (Dendrite* d : rootDendrites)
        if (d->S) I_dend_spike += 25.0f;
    for (Dendrite* d : distalDendrites)
        if (d->S) I_dend_spike += 10.0f;

    float V_dend_sum = 0.0f;
    int count_for_avg = 0;
    for (Dendrite* d : rootDendrites) { V_dend_sum += d->getVoltage(); ++count_for_avg; }
    for (Dendrite* d : distalDendrites) { V_dend_sum += d->getVoltage(); ++count_for_avg; }

    if (count_for_avg == 0) count_for_avg = 1;

    last_I_dend_spike = I_dend_spike;
    last_V_dend_sum = V_dend_sum;

    float V_dend_avg = V_dend_sum / static_cast<float>(count_for_avg);
    soma->step(t, S_soma, delta_t, mod, V_dend_avg,
        static_cast<int>(rootDendrites.size()), I_dend_spike);
}


float Neuron::getDendriticVoltageSum() const {
    float sum = 0.0f;
    for (Dendrite* d : rootDendrites)
        sum += d->getVoltage();
    return sum;
}

// --- Helper functions ---
Soma* Neuron::getSoma() {
    return soma;
}

Dendrite* Neuron::getDistalDendrite(int idx) {
    return (idx >= 0 && idx < distalDendrites.size()) ? distalDendrites[idx] : nullptr;
}

int Neuron::getNumDistalDendrites() const {
    return distalDendrites.size();
}

int Neuron::getSomaSpike() {
    return soma->S;
}