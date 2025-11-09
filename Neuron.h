#pragma once

#include <vector>
#include <cmath>

#include "../include/Variables.h"
#include "../include/utils.h"

// Forward declarations to avoid circular include
class Dendrite;
class Soma;

class Soma;

class Neuron {
    std::vector<Dendrite*> trees; /* dendritic node roots that send signal to the soma,
    they are connected to various subnodes in a tree-like structure,
    child nodes can vary in number, so can the number of trees */
    Soma* soma; // Soma of the neuron
    //int Id; // unused for now, stores the id of the neuron
    int _homeo_counter = 0;
    // [DendriteA[Input1, Input2, ...], DendriteB[Input1, Input2, ...], ...]
public:
    std::vector<Dendrite*> distalDendrites;
    std::vector<Dendrite*> rootDendrites; // dendritic nodes that send signal to the soma, they can contain other dendrites
    float last_I_dend_spike = 0.0f;
    float last_V_dend_sum = 0.0f;
    explicit Neuron(std::vector<float> neuromodmult, const std::vector<float> weightRange = { 0.1f, 0.5f }, const int maxDendriticDepth = 3, const int dendriticTreeCount = 1, const int branchingwidth = 2, bool inhibitory = false);
    void buildTree(Dendrite* parent, int remainingDepth,
        const std::vector<float>& weightRange, int branchingwidth, bool inhibitory, std::vector<float> neuromodmult);
    void traverse(Dendrite* dend);
    void collectDistal();
    //std::vector<float> homeostasis();
    void traverseStep(Dendrite* dend, const int t, const double delta_t, const std::vector<float>& mod);
    void step(const int t, const std::vector<std::vector<int>>& S, const double delta_t, const std::vector<float>& mod);
    // helpers
    Soma* getSoma();
    Dendrite* getDistalDendrite(int idx);
    int getNumDistalDendrites() const;
    int getSomaSpike();
    float getDendriticVoltageSum() const;
    float getDendriticCurrent() const { return last_I_dend_spike; }
};