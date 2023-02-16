#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_trajectories_ir_day.cpp `pkg-config --libs gsl` command for compiling

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
#define I_init 0        // initial number of immune cells

// General parameters
#define nu 0.01         // mutation rate S -> R 
#define t0 0            // initial time

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

// Main code
double RUN(double,double,double,int);

int main(int argc,char *argv[])
{
    double tfin = atof(argv[2]);                // treatment time (in days) (if = 0 this is the infinite treatment scenario, i.e. no treatment end)    
    double c = atof(argv[1]);                   // antimicrobial concentration
    double hR = atof(argv[3]);                  // half concentration of resistant strain
    int success = 0;                            // survival counter
    
    int r_ind;                                  // repetition index
    int repeats = 10;                           // number of repetitions
    
    while (success<repeats)
    {
        gsl_rng_set(r,r_ind);                     // setting the seed
        double time = RUN(c,tfin,hR,success);                  // simulation output
        if (time>0)                             // Condition for resistant extinction (and not sensitive!)
        {           
            success++;
            cout << success;
            cout << "\n";
        }
        r_ind++;                                // next simulation       
    }
    
    return(0);
}

// Stochastic simulation
double RUN(double c,double tfin,double hR,int ct)
{
    int S, R, I;                    // auxiliary variables (sensitive, resistant type)
    double t;                    // auxiliary variable (time)        
    double return_value;        // return value (extinction time)

    // File to save trajectory    
    ofstream file ("trajectory_" + std::to_string(ct) + "_highlow.txt");   

    // Initialization of the system
    t = (double)t0;
    S = S_init;
    R = R_init;
    I = I_init;
    
    file << S;
    file << ",";
    file << R;
    file << ",";
    file << I;
    file << ",";
    file << t;
    file << "\n";

    // Pharmacodynamics as in Regoes 2004 Eq. 3 - parameter values for ciprofloxacin
    double rS, rR;
    rS = bS*(1.-tanh(15.*(0.-hS)));
    rR = bR*(1.-tanh(15.*(0.-hR)));

    // Transition vector
    int transitions[2];
    transitions[0] = 1;
    transitions[1] = -1;
    
    // Simulation 
    while (S+R <= 100 && S+R >0)
    {
        // Update
        int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
        int ind = 0;            // index for the transition to update
        
        // Transition rate vector
        double rates[6];
        
        // Treatment phase transition rates
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

        file << S;
        file << ",";
        file << R;
        file << ",";
        file << I;
        file << ",";
        file << t;
        file << "\n";   
    }

    // update treatment
    double tchange = t;
    rS = bS*(1.-tanh(15.*(c-hS)));
    rR = bR*(1.-tanh(15.*(c-hR)));

    while (t <= tfin + tchange)
    {
        // Update
        int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
        int ind = 0;            // index for the transition to update
        
        // Transition rate vector
        double rates[6];
        
        // Treatment phase transition rates
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

        file << S;
        file << ",";
        file << R;
        file << ",";
        file << I;
        file << ",";
        file << t;
        file << "\n";   
    }

    // update end of treatment
    rS = bS*(1.-tanh(0.*(c-hS)));
    rR = bR*(1.-tanh(0.*(c-hR)));

    while (t <= 20.)
    {
        // Update
        int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
        int ind = 0;            // index for the transition to update
        
        // Transition rate vector
        double rates[6];
        
        // Treatment phase transition rates
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

        file << S;
        file << ",";
        file << R;
        file << ",";
        file << I;
        file << ",";
        file << t;
        file << "\n";   
    }
    
    return_value = t;

    file.close();

    return(return_value);
}
