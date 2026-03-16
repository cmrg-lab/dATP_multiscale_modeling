#ifndef PROBLEM_DEFINES_H
#define PROBLEM_DEFINES_H

#define N_RU 26 // Number of regulatory units (RUs) per thin filament
#define N_S 6 // Number of states per RU: B*, C*, B, C, M1 (Mc), M2 (Md)
#define MAX_REPS (640) // Max number used to repeat the simulation 32 blocks * 32 threads //150016
#define DT 5e-4f // fixed time step [ms] (0.5 us)
#define MAX_TSTEPS (4000001) // Max number of time steps; total simulation time = MAX_TSTEPS * DT = 2000 ms (2 s)
#define SEED (time(NULL))
#define RANDVAL (rand())
#endif // PROBLEM_DEFINES_H
