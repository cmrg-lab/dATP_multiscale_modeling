//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |          Function Name                |                         //
//                    |          force_pCa_curve ()           |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
#include "force_pCa_curve.h"
#include "problemDefines.h"
//#include "experimentalDataHost.h"
#include "gpuErrchk.h"
#include <math.h>
#include <stdio.h>
#include <boost/thread.hpp>
#include <boost/atomic.hpp>
//------------------
// Functions used
//-----------------
#include "rates_trans_matrix.h"
#include "update_RUs.h"
#include "repeat_simul.h"
#include "setGPU.h"
//-------------------

static boost::atomic<int> totalThreadsFinishedMallocing(0);
static boost::mutex lock;

//------------------------------------------
// The force_pCa_Curve function definition:
//------------------------------------------

void force_pCa_curve(initParticleArgs & args,
                     unsigned long randSeed,
                     float * ForceArrays,
                     float * Fss,
                     float * McArrays,
                     float * CArrays,
                     float * BArrays,
                     float * SRArrays,
                     int cc
                    )
{

// grab a new GPU to balance load
int GPUid = getGPU();
setGPU(GPUid);
// select beginning of this loop's Force array
float * Force = &(ForceArrays[cc * MAX_TSTEPS]);
float * Mc = &(McArrays[cc * MAX_TSTEPS]);
float * C = &(CArrays[cc * MAX_TSTEPS]);
float * B  = &(BArrays[cc * MAX_TSTEPS]);
float * SR  = &(SRArrays[cc * MAX_TSTEPS]);
float * kB_plus;
gpuErrchk(cudaMallocManaged(&kB_plus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(kB_plus, 0, sizeof(float)*N_S*N_S));
float * kB_minus;
gpuErrchk(cudaMallocManaged(&kB_minus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(kB_minus, 0, sizeof(float)*N_S*N_S));
float * k2_plus_dATP;
gpuErrchk(cudaMallocManaged(&k2_plus_dATP, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k2_plus_dATP, 0, sizeof(float)*N_S*N_S));
float * k2_plus_ATP;
gpuErrchk(cudaMallocManaged(&k2_plus_ATP, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k2_plus_ATP, 0, sizeof(float)*N_S*N_S));
float * k2_minus;
gpuErrchk(cudaMallocManaged(&k2_minus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k2_minus, 0, sizeof(float)*N_S*N_S));
float * k4_plus_dATP;
gpuErrchk(cudaMallocManaged(&k4_plus_dATP, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k4_plus_dATP, 0, sizeof(float)*N_S*N_S));
float * k4_plus_ATP;
gpuErrchk(cudaMallocManaged(&k4_plus_ATP, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k4_plus_ATP, 0, sizeof(float)*N_S*N_S));
float * k4_minus;
gpuErrchk(cudaMallocManaged(&k4_minus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k4_minus, 0, sizeof(float)*N_S*N_S));

float gamma_B = args.gamma_B;
float gamma_M = args.gamma_M;
float mu_B = args.gamma_M;
float mu_M = args.mu_M;
float kB_plus_ref = args.kB_plus_ref;
float kB_minus_ref = args.kB_minus_ref;
float k2_plus_ref_dATP = args.k2_plus_ref;
float k3_plus_dATP   = args.k3_plus;
float k4_plus_ref_dATP = args.k4_plus_ref;
float percent_dATP = args.percent_dATP;
float kCa_plus_ref = args.kCa_plus_ref;
float kCa_minus_ref = args.kCa_minus_ref;
float k_force_dATP = args.k_force;
float k_plus_SR_ref_dATP = args.k_plus_SR_ref;
float k_minus_SR_ref = args.k_minus_SR_ref;

float k2_plus_ref_ATP = 0.0025; // [1/ms] - XB attachment rate (ATP), C -> M1
float k3_plus_ATP = 0.05;       // [1/ms] - power stroke forward rate (ATP), M1 -> M2
float k4_plus_ref_ATP = 0.135;  // [1/ms] - XB detachment rate (ATP), M2 -> C
float k_plus_SR_ref_ATP = 16;   // [1/ms] - SRX to DRX transition rate (ATP)
float k_force_ATP = 0.2;        // [unitless] - force-feedback coefficient for SR-to-DRX transition (ATP)

//-------------------------------
//   Set rates using the input arguments
//-------------------------------
float r = 1; // [unitless] - XB-XB cooperativity exponent (fixed)
float q = 1; // [unitless] - RU-RU cooperativity exponent (fixed)
float lambda = 0; // [unitless, 0-1] - scales Ca2+ unbinding from Ca-bound non-permissive states (0 = no effect)
// calculating rates for XB cycling - use Tanner 2007/ Daniel 1998/ Pate & Cooke 1989
float k2_minus_ref, k3_minus, k4_minus_ref;
float conc_ADP,conc_Pi, conc_ATP, x_preR, g_Ca, g_Cb, g_Mc, g_Md, delta_G_ATP, delta_G, k_xb, x_xb;
//float  A, B, C, D, M, N, P, x_b0;
//metabolite concentrations in cytosol
conc_ADP    = 30;        //uM, Dawson et al 1978/ Kushmerick et al 1969 (frog)
conc_ATP    = 3e3;        //uM
conc_Pi     = 3e3;         //uM

//thermodynamic parameters
//r_gas         = 8.314;      // Gas constant, J/mol*K
//tc            = 15;
//temp          = tc + 273;   //temperature in Kelvin


// other constants
float alpha = 0.28; // [unitless] - fraction of ATP hydrolysis energy allocated to M1 (Mc) state
float eta = 0.68;   // [unitless] - fraction of ATP hydrolysis energy allocated to M2 (Md) state
//A = 2000; 
//B = 100; // all from Tanner et al, 2007.
//C = 1;
//D = 1;
//M = 3600;
//N = 40;
//P = 20;
k_xb = 5; // [RT/nm^2] - XB elastic spring constant

delta_G_ATP = 13; // units = RT
delta_G = delta_G_ATP - log(conc_ATP/(conc_ADP*conc_Pi)); // units = RT

x_preR      = 0;      // [nm] XB distortion in pre-power-stroke (pre-rotated) state
x_xb        = 0.075; // [nm] XB distortion in post-power-stroke (Md) state
//x_b0        = eta * delta_G / k_xb; // xb distortion due to ATP hydrolysis



g_Cb    =  0                                    ;// [RT] free energy of XB state Cb (detached, reference = 0)
g_Mc    = alpha * delta_G + k_xb * (x_preR)     ;// [RT] free energy of XB state Mc (weakly-bound M1)
g_Md    = eta* delta_G + k_xb*pow(x_xb,2)       ;// [RT] free energy of XB state Md (strongly-bound M2)
g_Ca    =   g_Cb                                ;// [RT] free energy of XB state Ca (same as Cb)


// to get reverse values, keep in mind that rij/rji = e^(gi - gj)
// so, to find r21 = r12/e^(g1 - g2)

//kCa_plus_ref    = 0.09;
//kCa_minus_ref   = 0.113;                    //X_kCa_minus_ref_PSO[i];
//kB_minus_ref    = 0.327;                    //X_kB_minus_ref_PSO[i];
//k2_plus_ref     = A * pow(k_xb/2*M_PI,0.5)*exp(-k_xb*pow(x_preR-x_b0,2)/2); // from tanner 2007
k2_minus_ref    = k2_plus_ref_ATP/ exp(g_Cb - g_Mc); // [1/ms] using vals from optimization_0227 (k2_plus = 0.615440)
//k3_plus         = (B/pow(k_xb,.5))*(1-tanh(C*pow(k_xb,.5)*(x_xb-x_b0)))+D;        //X_k3_plus_PSO[i];
k3_minus        = k3_plus_ATP / exp(g_Mc - g_Md) ;// [1/ms] reverse power stroke, M2 -> M1
//k4_plus_ref     = pow(k_xb,0.5)*(pow(M*pow(x_xb,2),0.5)-N*x_xb)+ P;                 //X_k4_plus_PSO[i];
k4_minus_ref    = k4_plus_ref_ATP * exp(g_Ca - g_Md - delta_G); // [1/ms] XB reattachment, C -> M2

//-------------------------------------
// Call the transition rates function:
//-------------------------------------

rates_trans_matrix(N_S,
kB_plus_ref,
kB_minus_ref,
k2_plus_ref_dATP,
k2_plus_ref_ATP,
k2_minus_ref,
k4_plus_ref_dATP,
k4_plus_ref_ATP,
k4_minus_ref,
gamma_B,
gamma_M,
mu_B,
mu_M,
r,
q,
kB_plus,
kB_minus,
k2_plus_dATP,
k2_plus_ATP,
k2_minus,
k4_plus_dATP,
k4_plus_ATP,
k4_minus
);

    //-----------------------------------------------------
    // start Ca- loop i.e., to get the entire F-Ca curve:
    //-----------------------------------------------------

    float Ftemp        = 0.0;                      // is used to calculate the steady-state force at the end
    float Cal_conc     = pow(10.0f,-(args.experimentalData[cc].first-6)); // [uM] Ca2+ concentration (pCa=6 -> 1 uM)
    float kCa_plus     = Cal_conc*kCa_plus_ref;  // [1/ms] = [uM]*[1/(uM*ms)] - Ca2+ on rate
    float kCa_minus    = kCa_minus_ref;           // [1/ms] - Ca2+ off rate
    const int n_pCa = args.experimentalData.size();
    
    
    //---------------------------------
    // Call the repeat_simul function:
    //----------------------------------
    cudaStream_t s;
    cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
    totalThreadsFinishedMallocing++;
    while(totalThreadsFinishedMallocing < n_pCa){
        boost::thread::yield();
    }
repeat_simul<<<MAX_REPS/32, 32, 0, s>>>(lambda,
                                        randSeed,
                                        k4_plus_dATP,
                                        k4_plus_ATP,
                                        k4_minus,
                                        k3_plus_dATP,
                                        k3_plus_ATP,
                                        k3_minus,
                                        k2_plus_dATP,
                                        k2_plus_ATP,
                                        k2_minus,
                                        kB_plus,
                                        kB_minus,
                                        kCa_plus,
                                        kCa_minus,
                                        percent_dATP,
                                        k_force_dATP,
                                        k_force_ATP,
					k_plus_SR_ref_dATP,
					k_plus_SR_ref_ATP,
					k_minus_SR_ref,
                                        Force,
                                        Mc,
                                        C,
                                        B,
                                        SR,
                                        cc
                                           );

    gpuErrchk(cudaStreamSynchronize(s)); // wait for device to finish repeat_simul
    gpuErrchk(cudaStreamDestroy(s));
    //--------------------------------------------------------------------------------------
    // Calculate The Steady-State Force using Impluse using data from the last 0.5 sec
    // (i.e., just 100000 time steps) only using numerical trapaziodal integration
    //--------------------------------------------------------------------------------------

    for (int n = MAX_TSTEPS-100000; n < MAX_TSTEPS-1; n++)  // time marching
    {
        Ftemp = Ftemp+Force[n];
    }

    Fss[cc] = (Ftemp + (0.5f * Force[MAX_TSTEPS-100001]) + (0.5f * Force[MAX_TSTEPS-1])) / 100000.0f / MAX_REPS;    //Fss[cc] = 1;

    //--------------------------------

    // free allocated memory
gpuErrchk(cudaFree(kB_plus));
gpuErrchk(cudaFree(kB_minus));
gpuErrchk(cudaFree(k2_plus_dATP));
gpuErrchk(cudaFree(k2_plus_ATP));
gpuErrchk(cudaFree(k2_minus));
gpuErrchk(cudaFree(k4_plus_dATP));
gpuErrchk(cudaFree(k4_plus_ATP));
gpuErrchk(cudaFree(k4_minus));

    lock.lock();
    if(totalThreadsFinishedMallocing == n_pCa){
        totalThreadsFinishedMallocing = 0;
    }
    lock.unlock();

} // end main function
