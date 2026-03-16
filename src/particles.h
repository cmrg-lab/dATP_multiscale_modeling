#ifndef PARTICLES_H
#define PARTICLES_H

#include "problemDefines.h"
#include <vector>

class initParticleArgs
{
public:
    std::vector< std::pair<float, float> > experimentalData;
    //---------------------------------------------
    // model reference parameters that we need to optimize
    //--------------------------------------------
    float gamma_B; // [unitless] - RU-RU cooperative coefficient
    float gamma_M; // [unitless] - XB-RU/RU-XB cooperative coefficient (Note: gamma_M = mu_B)
    float mu_M; // [unitless] - inter-RU XB-XB cooperative coefficient
    float k2_plus_ref; // [1/ms] - XB attachment rate (dATP), C -> M1
    float k3_plus;     // [1/ms] - power stroke forward rate (dATP), M1 -> M2
    float k4_plus_ref; // [1/ms] - XB detachment rate (dATP), M2 -> C
    float kB_plus_ref;  // [1/ms] - RU ON rate, B*/B -> C*/C  (p = plus)
    float kB_minus_ref; // [1/ms] - RU OFF rate, C*/C -> B*/B (m = minus)
    float lambda; // [unitless, 0-1] - scales Ca2+ unbinding from Ca-bound non-permissive states
    float kCa_plus_ref;  // [1/(uM*ms)] - Ca2+ binding rate to troponin C
    float kCa_minus_ref; // [1/ms] - Ca2+ unbinding rate from troponin C
    float percent_dATP;  // [0-1] - fraction of ATP replaced by dATP (0 = all ATP, 1 = all dATP)
    float k_force;       // [unitless] - force-feedback coefficient for SR-to-DRX transition rate
    float k_plus_SR_ref;  // [1/ms] - SRX (super-relaxed) to DRX (disordered-relaxed) transition rate
    float k_minus_SR_ref; // [1/ms] - DRX to SRX transition rate
    initParticleArgs(std::vector< std::pair<float, float> > experimentalData,
                     float gamma_B,
                     float gamma_M,
                     float mu_M,
                     float k2_plus_ref,
                     float k3_plus,
                     float k4_plus_ref,
                     float kB_plus_ref,
                     float kB_minus_ref,
                     float lambda,
                     float kCa_plus_ref,
                     float kCa_minus_ref,
                     float percent_dATP,
                     float k_force,
                     float k_plus_SR_ref,
                     float k_minus_SR_ref):experimentalData(experimentalData),
        gamma_B(gamma_B),
        gamma_M(gamma_M),
        mu_M(mu_M),
        k2_plus_ref(k2_plus_ref),
        k3_plus(k3_plus),
        k4_plus_ref(k4_plus_ref),
        kB_plus_ref(kB_plus_ref),
        kB_minus_ref(kB_minus_ref),
        lambda(lambda),
        kCa_plus_ref(kCa_plus_ref),
        kCa_minus_ref(kCa_minus_ref),
        percent_dATP(percent_dATP),
        k_force(k_force),
        k_plus_SR_ref(k_plus_SR_ref),
        k_minus_SR_ref(k_minus_SR_ref)
       {}


};

void init_particle(initParticleArgs & args);

#endif // PARTICLES_H
