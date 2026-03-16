//--------------------------------------------------------------------------------------//
//                     University of California San Diego                               //
//                         Dept. of Bioengineering                                      //
//                      Cardiac Mechanics Research Group                                //
//                         PI: Prof. Andrew McCulloch                                   //
//--------------------------------------------------------------------------------------//
// Authors: Kim McCabe, Yasser Aboelkassem,  & Andrew McCulloch                         //
// Year  :  08/2017    									//
// Updated: Abby Teitgen 05/2023                                                        //
//-----------------------------                                                         //
//          This code uses the "Markov Chain Monte Carlo" algorithm to solve            //
//           cardiac thin filament activation / crossbridge cycling problem.            //
//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |      The TF/XB 10 State Model         |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
#include "problemDefines.h"

#include "gpuErrchk.h"
#include "setGPU.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//--------------------------
// Function to be called
//--------------------------
#include "particles.h"
#include "csvReader.h"

//-------------------------
// Main body code
//------------------------
int main(int argc, const char *argv[])
{
if(argc < 2){
std::cerr << "Experimental data CSV needed as first argument." << std::endl;
return 1;
}
if(argc < 11){
std::cerr << "Parameter list needed." << std::endl;
return 1;
}
std::vector< std::pair<float, float> > experimentalData = csvReader(argv[1]);
srand(SEED); //Random-Seed initialization (must be outside any loop)
std::cout << "SEED: " << SEED << std::endl;
long long startTime = time(NULL);
// make sure new GPUs are used
initGPUSelection();

//---------------------------------------------
// Model reference parameters that we need to optimize
//--------------------------------------------
float gamma_B = atof(argv[2]); // [unitless] - RU-RU cooperative coefficient
float gamma_M = atof(argv[3]); // [unitless] - XB-RU/RU-XB cooperative coefficient (Note: gamma_M = mu_B)
float mu_M = atof(argv[4]); // [unitless] - inter-RU XB-XB cooperative coefficient
float k2_plus_ref = atof(argv[5]); // [1/ms] - XB attachment rate (dATP), C -> M1
float k3_plus = atof(argv[6]);     // [1/ms] - power stroke forward rate (dATP), M1 -> M2
float k4_plus_ref = atof(argv[7]); // [1/ms] - XB detachment rate (dATP), M2 -> C
float kB_plus_ref = atof(argv[8]); // [1/ms] - RU ON rate, B*/B -> C*/C  (p = plus)
float kB_minus_ref = atof(argv[9]); // [1/ms] - RU OFF rate, C*/C -> B*/B (m = minus)
float lambda = atof(argv[10]); // [unitless, 0-1] - scales Ca2+ unbinding from Ca-bound non-permissive states
float kCa_plus_ref = atof(argv[11]); // [1/(uM*ms)] - Ca2+ binding rate to troponin C
float kCa_minus_ref = atof(argv[12]); // [1/ms] - Ca2+ unbinding rate from troponin C
float percent_dATP = atof(argv[13]); // [0-1] - fraction of ATP replaced by dATP (0 = all ATP, 1 = all dATP)
float k_force = atof(argv[14]); // [unitless] - force-feedback coefficient for SR-to-DRX transition rate
float k_plus_SR_ref = atof(argv[15]); // [1/ms] - SRX (super-relaxed) to DRX (disordered-relaxed) transition rate
float k_minus_SR_ref = atof(argv[16]); // [1/ms] - DRX to SRX transition rate

initParticleArgs args = initParticleArgs(experimentalData,
gamma_B,
gamma_M,
mu_M,
k2_plus_ref,
k3_plus,
k4_plus_ref,
kB_plus_ref,
kB_minus_ref,
lambda,
kCa_plus_ref,
kCa_minus_ref,
percent_dATP,
k_force,
k_plus_SR_ref,
k_minus_SR_ref
);
init_particle(args);
std::cout << "One iteration runtime: " << (time(NULL)-startTime) << " second(s)" << std::endl;

return 0;

} // end main function
