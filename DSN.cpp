#include "DSN.h"

//excitatory
std::vector<float> dendrititicRule_P = { 0.45, 0.35, 0.15, 0.03, 0.02 };
std::vector<float> somaticRule_P = { 0.05, 0.25, 0.25, 0.10, 0.35 };
//inhibitory
std::vector<float> dendrititicRule_I_P = { 0.05, 0.15, 0.10, 0.40, 0.30 };
std::vector<float> somaticRule_I_P = { 0.02, 0.10, 0.10, 0.35, 0.043 };