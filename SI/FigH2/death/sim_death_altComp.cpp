#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_death_altComp.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

// Sensitive strain parameters
#define bS 2.5          // birth rate sensitive
#define dS 0.5          // death rate sensitive
#define cS 1.           // competition rate sensitive

// Resistant strain parameters
#define bR 2.25          // birth rate resistant
#define dR 0.5          // death rate resistant
#define cR 1.           // competition rate resistant
#define R_init 1        // initial value resistant

// General parameters
#define K 1000          // Carrying capacity scaling
#define nu 0            // mutation rate S -> R 
#define t0 0            // initial time

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

// Main code
int RUN(double,double,int,double);

int main(int argc,char *argv[])
{
    double tfin = atof(argv[2]);      // treatment time (in days) (if = 0 this is the infinite treatment scenario, i.e. no treatment end)    
    double c = atof(argv[1]);                          // antimicrobial concentration
    double fac = atof(argv[4]);
    int scenario = atoi(argv[3]);                 // 0 = bacteriostatic, 1 = bactericidal
    int success = 0;                             // survival counter
    
    int r_ind;                                   // repetition index
    int repeats = 100000;                         // number of repetitions
    
    for (r_ind=0; r_ind<repeats; r_ind++)
    {
        gsl_rng_set(r,r_ind);                    // setting the seed
        success += RUN(c,tfin,scenario,fac);                  // updating survival counter (+1 if resistant strain survived)
    }
    
    double avg;                                  // empirical survival probability of R
    avg = (double)success/(double)repeats;
    
    ofstream file ("surv_death_vary_c_fac_" + std::to_string((int)fac) + "_sc_" + std::to_string(scenario) + ".txt", ios::app);   // file output, average
    file << c;
    file << ",";
    file << avg;  
    file << "\n";
    file.close();
    
    return(0);
}

// Stochastic simulation
int RUN(double c,double tfin,int scenario,double fac)
{
    int S, R;                    // auxiliary variables (sensitive, resistant type)
    double t;                    // auxiliary variable (time) 
    int return_value = 0;        // return value (0 if R=0, 1 if R>1 at tfin)
    double gamma_alt = cS/dS;    // reparameterize competition to fit with original model
        
    // Initialization of the system
    t = (double)t0;
    S = (int)((double)K*(bS-dS)/cS);
    R = R_init;
    
    // Pharmacodynamics as in Regoes 2004 Eq. 3 - parameter values for ciprofloxacin
    double kappaS, kappaR, micS, micR, psimaxS, psimaxR, psiminS, psiminR;
    kappaS = 1.1;                            // steepness antimicrobial reaction for sensitive
    kappaR = kappaS;                            // steepness antimicrobial reaction for resistant
    micS = 0.017;                             // minimal inhibitory concentration for sensitive
    micR = micS *fac;                            // minimal inhibitory concentration for resistant
    psiminS = log(10)*(-6.5)*24.;                // minimal growth rate for any antimicrobial concentration sensitive
    psiminR = psiminS;                       // minimal growth rate for any antimicrobial concentration resistant
    psimaxS = bS-dS;                         // exponential growth rate sensitive
    psimaxR = bR-dR;                         // exponential growth rate resistant
    
    double aS, aR;
    aS = (psimaxS-psiminS)*pow(c/micS,kappaS) / (pow(c/micS,kappaS) - psiminS/psimaxS);
    aR = (psimaxR-psiminR)*pow(c/micR,kappaR) / (pow(c/micR,kappaR) - psiminR/psimaxR);

    // Transition vector
    int transitions[2];
    transitions[0] = 1;
    transitions[1] = -1;
    
    // Simulation 
    
    // Infinite treatment scenario
    if (tfin == 0)
    {
        while (R>0 && S>0)
        {
            // Update
            int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
            int ind = 0;            // index for the transition to update
            
            // Transition rate vector
            double rates[4];
            
            // Treatment phase transition rates
            // Bacteriostatic scenario
            if (scenario == 0)
            {
                rates[0] = max(0.,(bS-aS)*(double)S);        // birth sensitive
                rates[1] = dS*(1.+ gamma_alt*((double)S+(double)R)/(double)K) *(double)S;                // death sensitive
                rates[2] = max(0.,(bR-aR)*(double)R);        // birth resistant
                rates[3] = dR*(1.+ gamma_alt*((double)S+(double)R)/(double)K) *(double)R;                // death resistant
            }
            
            // Bactericidal scenario
            if (scenario == 1)
            {
                rates[0] = max(0.,(bS)*(double)S);              // birth sensitive
                rates[1] = (dS + aS)*(1. + gamma_alt*((double)S+(double)R)/(double)K) *(double)S;                // death sensitive
                rates[2] = max(0.,(bR)*(double)R);              // birth resistant
                rates[3] = (dR + aR)*(1.+ gamma_alt*((double)S+(double)R)/(double)K) *(double)R;                // death resistant
            }     
                  
            // Draw two uniform random numbers for updating
            double rand1, dt;
            rand1 = gsl_ran_flat(r, 0.0, 1.0);  
            dt = gsl_ran_exponential(r, 1/accumulate(rates,rates+4,0.));  
        
            // Time update
            t += dt;
            
            // Population update
            while (update == 0)
            {
                if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+4,0.))
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
                    if (ind >= 2)
                    {
                        if (ind%2 == 0)              //birth event
                        {
                            R++;
                        }
                        else R--;
                    }
                }
                
                ind++;
            }       
        }
        
        if (R > 0)
        {
            return_value = 1;
        }
    }
    
    // Finite treatment scenario
    else
    {
        while (t <= tfin && R>0)
        {
            // Update
            int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
            int ind = 0;            // index for the transition to update
            
            // Transition rate vector
            double rates[4];
            
            // Treatment phase transition rates
            // Bacteriostatic scenario
            if (scenario == 0)
            {
                rates[0] = max(0.,(bS-aS)*(double)S);        // birth sensitive
                rates[1] = dS*(1.+ gamma_alt*((double)S+(double)R)/(double)K) *(double)S;                // death sensitive
                rates[2] = max(0.,(bR-aR)*(double)R);        // birth resistant
                rates[3] = dR*(1.+ gamma_alt*((double)S+(double)R)/(double)K) *(double)R;                // death resistant
            }
            
            // Bactericidal scenario
            if (scenario == 1)
            {
                rates[0] = max(0.,(bS)*(double)S);              // birth sensitive
                rates[1] = (dS + aS)*(1. + gamma_alt*((double)S+(double)R)/(double)K) *(double)S;                // death sensitive
                rates[2] = max(0.,(bR)*(double)R);              // birth resistant
                rates[3] = (dR + aR)*(1.+ gamma_alt*((double)S+(double)R)/(double)K) *(double)R;                // death resistant
            }      
                  
            // Draw two uniform random numbers for updating
            double rand1, dt;
            rand1 = gsl_ran_flat(r, 0.0, 1.0);  
            dt = gsl_ran_exponential(r, 1./accumulate(rates,rates+4,0.));  
        
            // Time update
            t += dt;
            
            // Population update
            while (update == 0)
            {
                if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+4,0.))
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
                    if (ind >= 2)
                    {
                        if (ind%2 == 0)              //birth event
                        {
                            R++;
                        }
                        else R--;
                    }
                }
                
                ind++;
            }       
        }
        
        if (R > 0)
        {          
            return_value = 1;
        }
    }  
    
    return(return_value);
}
