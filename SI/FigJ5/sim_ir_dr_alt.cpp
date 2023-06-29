#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_ir_dr_alt.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

// Sensitive strain parameters
#define bS 0.6          // birth rate sensitive
#define dS 0.01         // death rate sensitive
#define S_init 10       // initial number of sensitive pathogens
#define hS 0.3          // half concentration of sensitive strain

// Resistant strain parameters
#define bR 0.59         // birth rate resistant
#define dR 0.01         // death rate resistant
#define R_init 0        // initial number of resistant pathogens

// Immune response parameters
#define xi 0.075        // attack rate of immune cells
#define alpha 0.05      // production of immune cells
#define delta 0.05      // degradation rate of immune cells
#define I_init 2        // initial number of immune cells

// General parameters
#define nu 0.01         // mutation rate S -> R 
#define t0 0            // initial time

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

// Main code
int RUN(double,double,double);

int main(int argc,char *argv[])
{
    double tfin = atof(argv[2]);                // treatment time (in days) (if = 0 this is the infinite treatment scenario, i.e. no treatment end)    
    double c = atof(argv[1]);                   // antimicrobial concentration
    double hR = atof(argv[3]);                  // half concentration of resistant strain
    int success = 0;                            // survival counter
    
    int r_ind;                                  // repetition index
    int repeats = 10000;                        // number of repetitions
    
    for (r_ind=0; r_ind<repeats; r_ind++)
    {
        gsl_rng_set(r,r_ind);                   // setting the seed
        success += RUN(c,tfin,hR);              // updating survival counter (+1 if resistant strain survived)
    }
    
    double avg;                                 // empirical survival probability of R
    avg = (double)success/(double)repeats;
    
    ofstream file ("survIR_day_vary_c_factor_" + std::to_string(hR) + "_alt.txt", ios::app);   // file output, average
    file << c;
    file << ",";
    file << avg;  
    file << "\n";
    file.close();
      
    return(0);
}

// Stochastic simulation
int RUN(double c,double tfin,double hR)
{
    int S, R, I;                 // auxiliary variables (sensitive, resistant type, immune cells)
    double t;                    // auxiliary variable (time) 
    int return_value = 0;        // return value (0 if R=0, 1 if R>1 at tfin)
        
    // Initialization of the system
    t = (double)t0;
    S = S_init;
    R = R_init;
    I = I_init;
    
    // Pharmacodynamics as in Regoes 2004 Eq. 3 - parameter values for ciprofloxacin
    double rS, rR;
    rS = bS*(1.-tanh(15.*(0.-hS)));
    rR = bR*(1.-tanh(15.*(0.-hR)));

    // Transition vector
    int transitions[2];
    transitions[0] = 1;
    transitions[1] = -1;
    
    // Simulation 
    while (S+R<100 && S+R>0)
    {
        // Update
        int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
        int ind = 0;            // index for the transition to update
        
        // Transition rate vector
        double rates[6];
            
        // Treatment phase transition rates
        rates[0] = max(0.,rS*(double)S);                                   // birth sensitive
        rates[1] = dS *(double)S + xi* (double)S * (double)I;                // death sensitive
        rates[2] = max(0.,rR*(double)R);                                   // birth resistant
        rates[3] = dR *(double)R + xi* (double)R * (double)I;                // death resistant
        rates[4] = alpha * ((double)S + (double)R);                             // production immune cell
        rates[5] = delta * (double)I;                                           // degradation immune cell
             
                  
        // Draw two uniform random numbers for updating
        double rand1, dt;
        rand1 = gsl_ran_flat(r, 0.0, 1.0);   
        dt = gsl_ran_exponential(r, 1/accumulate(rates,rates+6,0.));  
        
        // Time update
        t += dt;
            
        // Population update
        while (update == 0)
        {
            if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+6,0.))
            {
                update = 1;
                   
                // Update sensitive                
                if (ind <= 1)
                {
                    if (ind%2 == 0)             //birth event
                    {
                        // Mutation event
                        if (gsl_ran_bernoulli(r, nu)==1)
                        {
                            R++;
                        }
                        else S++;
                    }
                    else S--;
                }
                    
                // Update resistant
                else if (ind <= 3)
                {
                    if (ind%2 == 0)              //birth event
                    {
                        R++;
                    }
                    else R--;
                }

                // Update immune cells
                else
                {
                    if (ind%2 == 0)             // production event
                    {
                        I++;
                    }
                    else I--;
                }
            }
                
            ind++;
        }       
        if (R>=100) return_value = 1;
    }

    // update treatment
    rS = bS*(1.-tanh(15.*(c-hS)));
    rR = bR*(1.-tanh(15.*(c-hR)));
    double tchange = t;

    while (t<=tfin + tchange)
    {
        // Update
        int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
        int ind = 0;            // index for the transition to update
        
        // Transition rate vector
        double rates[6];
            
        // Treatment phase transition rates
        rates[0] = max(0.,rS*(double)S);                                   // birth sensitive
        rates[1] = dS *(double)S + xi* (double)S * (double)I;                // death sensitive
        rates[2] = max(0.,rR*(double)R);                                   // birth resistant
        rates[3] = dR *(double)R + xi* (double)R * (double)I;                // death resistant
        rates[4] = alpha * ((double)S + (double)R);                             // production immune cell
        rates[5] = delta * (double)I;                                           // degradation immune cell
             
                  
        // Draw two uniform random numbers for updating
        double rand1, dt;
        rand1 = gsl_ran_flat(r, 0.0, 1.0);   
        dt = gsl_ran_exponential(r, 1/accumulate(rates,rates+6,0.));  
        
        // Time update
        t += dt;
            
        // Population update
        while (update == 0)
        {
            if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+6,0.))
            {
                update = 1;
                   
                // Update sensitive                
                if (ind <= 1)
                {
                    if (ind%2 == 0)             //birth event
                    {
                        // Mutation event
                        if (gsl_ran_bernoulli(r, nu)==1)
                        {
                            R++;
                        }
                        else S++;
                    }
                    else S--;
                }
                    
                // Update resistant
                else if (ind <= 3)
                {
                    if (ind%2 == 0)              //birth event
                    {
                        R++;
                    }
                    else R--;
                }

                // Update immune cells
                else
                {
                    if (ind%2 == 0)             // production event
                    {
                        I++;
                    }
                    else I--;
                }
            }
                
            ind++;
        }       
        if (R>=100) return_value = 1;
    }
    
    return(return_value);
} 
