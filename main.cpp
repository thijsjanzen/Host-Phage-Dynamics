//
//  main.cpp
//  Host_Gillespie
//
//  Created by Thijs Janzen on 16/08/16.
//  Copyright (c) 2016 Thijs Janzen. All rights reserved.
//

#include <iostream>
#include "randomc.h"
#include <vector>
#include <cmath>
#include <fstream>
#include "GetParams.h"
#include <unistd.h>  // for mac specific run code
#include <time.h>
#include <numeric>
#include <algorithm>

void macstart(const char * argv[]);  // forward declaration

struct Host{
    double Resource;
    double reproProb;
    double reproSize;

    Host() {
        Resource = 0.0;
        reproProb = 0.0;
    }

    void updateReproProb(double& totalProb, double R_r) {
        totalProb -= reproProb;
        if(Resource <= 0.0) {
            reproProb = 0.0;
        } else {
            reproProb = exp(-0.5 * (R_r - Resource));
        }
        totalProb += reproProb;
    }

    void uptake(double& R, double H, double alpha) {
        double uptake = alpha * R/(R+H);
        Resource += uptake;
        R -= uptake;
    }

    void maintenance(double c) {
        Resource -= c;
    }

    Host& operator=(const Host& other) {
        Resource = other.Resource;
        reproProb = other.reproProb;
        return *this;
    }
};


struct Host_infected {
    double Resource;
    int numberPhages;
    int lysisTime;

    Host_infected() {
        Resource = 0.0;
        numberPhages = 0.0;
        lysisTime = 1e6;
    }

    Host_infected& operator=(const Host_infected& other) {
        Resource = other.Resource;
        numberPhages = other.numberPhages;
        lysisTime = other.lysisTime;
        return *this;
    }
};




int pickEvent(const std::vector<long double>& v, const double& sum) {
    double r = uniform() * sum;
    for(int i = 0; i < v.size(); ++i) {
        r -= v[i];
        if(r <= 0.0) return i;
    }
    return (int)v.size() -1;
}

int pickIndiv(const std::vector<Host>& v, const double& maxProb) {

    while(1 == 1) {
        int index = random_number((int)v.size());
        double prob = v[index].reproProb / maxProb;
        if(uniform() < prob) {
            return index;
        }
    }

    return (int)(v.size() -1);
}

bool compFn(const Host& A, const Host& B) {
    return A.reproProb < B.reproProb;
}


void updateReproProb(const std::vector<Host>& v, double& maxProb) {
    auto it = std::max_element(v.begin(), v.end(), compFn);
    maxProb = (*it).reproProb;
}

template <typename T>
double calculateSD(const std::vector<T>& v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();


    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
    return stdev;
}



void doSimulation(const GetParams& Params, int repl)
{
    std::vector<Host> population;
    std::vector<Host_infected> population_infected;

    double totalReproProb = 0.0;
    double maxReproProb = 0.0;


    for(int i = 0; i < Params.initHostPopSize; ++i) {
        Host temp = Host();
        temp.reproSize = Params.reproductionSize;

        temp.updateReproProb(totalReproProb, Params.reproductionSize);

        if(temp.reproProb > maxReproProb) {
            maxReproProb = temp.reproProb;
        }
        totalReproProb += temp.reproProb;
        population.push_back(temp);
    }

    double t = 0.0;
    double R = 0.0;


    std::ofstream outFile("outFile.txt");

    double prevT = 0.0;

    double avgR = 0.0;
    for(auto it = population.begin(); it != population.end(); ++it) {
        avgR += (*it).Resource;
    }
    avgR = avgR / population.size();

    //   outFile << t << "\t" << population.size() << "\t" << R << "\t" << avgR << "\n";



    int numPhages = 0;

    while(t < Params.maxTime) {
        long double uptakeRate = population.size() * Params.alpha * R / (R+Params.H);
        long double maintenanceRate = population.size() * Params.maintenanceCost;
        long double reproRate  = totalReproProb;
        long double infectRate = population.size() * numPhages * Params.omega;
        long double phageDecay = numPhages * Params.gamma;
        long double Inf_Hosts_Rate = (int)population_infected.size() * Params.lambda / Params.beta;
        long double Inf_Hosts_maintenance = (int)population_infected.size() * Params.maintenanceCost;

        long double totalRate = uptakeRate + maintenanceRate + reproRate + infectRate + phageDecay + Inf_Hosts_Rate + Inf_Hosts_maintenance;
        double dt = Expon(totalRate);

        std::vector<long double> rates = {maintenanceRate, uptakeRate, reproRate,infectRate, phageDecay, Inf_Hosts_Rate,Inf_Hosts_maintenance};
        int event = pickEvent(rates, totalRate);

        //update flow effects
        R -= Params.D * R * dt;
        R += Params.Inflow * Params.D * dt;
        // remove individuals that flow out:
        double lossProb = Params.D * dt;
        if(lossProb > 1) {
            break;
        }

        if(lossProb > 0 && !population.empty()) {
            int loss = Binom((int)population.size(), lossProb);
            
            for(int i = 0; i < loss; ++i) {
                int index = (int)random_number((int)population.size());
                totalReproProb -= population[index].reproProb;

                if(population[index].reproProb == maxReproProb) {
                    population[index] = population.back();
                    population.pop_back();
                    updateReproProb(population,maxReproProb);
                } else {
                    population[index] = population.back();
                    population.pop_back();
                }
            }
        }

        //infected hosts:
        if(lossProb > 0 && !population_infected.empty()) {
            int loss = Binom((int)population_infected.size(), lossProb);

            for(int i = 0; i < loss; ++i) {
                int index = (int)random_number((int)population_infected.size());
                population_infected[index] = population_infected.back();
                population_infected.pop_back();
            }
        }

        if(numPhages > 0 && lossProb > 0) {
            int loss = Binom((int)numPhages, lossProb);
            numPhages -= loss;
            if(numPhages < 0) {
                numPhages = 0;
            }
        }

        if((int)population.size() <= 0) {
            break;
        }

        switch (event)
        {
            case 0: //maintenance
            {
                int index =random_number((int)population.size());
                population[index].maintenance(Params.maintenanceCost);

                if( population[index].reproProb == maxReproProb) {
                    population[index].updateReproProb(totalReproProb,Params.reproductionSize);
                    updateReproProb(population,maxReproProb);
                } else {
                    population[index].updateReproProb(totalReproProb,Params.reproductionSize);
                }

                if(population[index].Resource <= 0.0) {
                    population[index] = population.back();
                    population.pop_back();
                }
                break;
            }
            case 1: //uptake
            {
                int index = random_number((int)population.size());
                population[index].uptake(R, Params.H, Params.alpha);
                population[index].updateReproProb(totalReproProb, Params.reproductionSize);
                if(population[index].reproProb > maxReproProb) {
                    maxReproProb = population[index].reproProb;
                }
                break;
            }
            case 2:
            { //reproduction
                int reproIndiv = pickIndiv(population, maxReproProb);

                population[reproIndiv].Resource *= 0.5;

                if( population[reproIndiv].reproProb == maxReproProb) {
                    population[reproIndiv].updateReproProb(totalReproProb,Params.reproductionSize);
                    updateReproProb(population,maxReproProb);
                } else {
                    population[reproIndiv].updateReproProb(totalReproProb,Params.reproductionSize);
                }


                Host kid = population[reproIndiv];
                totalReproProb += kid.reproProb;

                population.push_back(kid);
                break;
            }
            case 3: //infection
            {
                int index = random_number((int)population.size());
                Host_infected new_infected_host;
                new_infected_host.Resource = population[index].Resource;
                population_infected.push_back(new_infected_host);

                population[index].updateReproProb(totalReproProb, Params.reproductionSize);
                totalReproProb -= population[index].reproProb;

                if(population[index].reproProb == maxReproProb) {
                    population[index] = population.back();
                    population.pop_back();
                    updateReproProb(population,maxReproProb);
                } else {
                    population[index] = population.back();
                    population.pop_back();
                }
                numPhages--; //one phage has infected and is lost from the population
                break;
            }
            case 4: //phage decay
            {
                numPhages--;
                break;
            }
            case 5: //update infected hosts
            {
                if(!population_infected.empty()) {
                    int index = random_number((int)population_infected.size());
                    if(population_infected[index].Resource >= Params.lambda) {
                        population_infected[index].Resource -= Params.lambda;
                        population_infected[index].numberPhages++;
                    } else {
                        //lysis!!!
                        numPhages += population_infected[index].numberPhages;
                        population_infected[index] = population_infected.back();
                        population_infected.pop_back();
                    }
                }
                break;
            }

            case 6: //maintentance infected hosts
            {
                if(!population_infected.empty()) {
                    int index = random_number((int)population_infected.size());
                    population_infected[index].Resource -= Params.maintenanceCost;
                    if(population_infected[index].Resource <= 0.0) {
                        //lysis!
                        numPhages += population_infected[index].numberPhages;
                        population_infected[index] = population_infected.back();
                        population_infected.pop_back();
                    }

                }
            }
        }


        t += dt;
        if(t >= Params.infectionTime && (t-dt) < Params.infectionTime) {
            numPhages = Params.initPhagePopSize;
        }

        double currentT = t;
        if(currentT - prevT >= 0.05) {
            double avgR = 0.0;
            for(auto it = population.begin(); it != population.end(); ++it) {
                avgR += (*it).Resource;
            }
            avgR = avgR / population.size();

            double threshold = Params.lambda / (1-Params.beta * Params.gamma);

            int N = (int)population.size();
            double a = (Params.lambda + Params.maintenanceCost) * (Params.gamma + N * Params.omega);
            double b = (N*Params.omega * (1 - Params.gamma * Params.beta));
            double threshold2 = a/b;

            outFile << Params.lambda << "\t" << Params.beta << "\t" << repl << "\t" << t << "\t" << population.size() << "\t" << numPhages << "\t" << R << "\t" << avgR  << "\t" << threshold << "\t" << population_infected.size() << "\t" << threshold2 << "\n";


            outFile.flush();
            prevT = currentT;
            std::cout << t << "\t" << population.size() << "\t" << numPhages << "\t" << avgR << "\t" << threshold2 << "\n";
        }


        if(population.size() <= 0.0) {
            break;
        }
    }

    avgR = 0.0;
    for(auto it = population.begin(); it != population.end(); ++it) {
        avgR += (*it).Resource;
    }
    avgR = avgR / population.size();

    return;
}


int main(int argc, const char * argv[]) {

    macstart(argv);
    GetParams P;
    P.readFromIni("config.ini");
    set_seed(P.seed);

    clock_t start2 = clock();
    doSimulation(P,0);
    clock_t end2 = clock();
    double elapsed_secs2 = double(end2 - start2) / CLOCKS_PER_SEC;
    std::cout << "that took: " << elapsed_secs2 << " seconds\n";

    return 0;
}

void macstart(const char * argv[])  {
    std::cout << "\n\n\n";
#ifdef __APPLE__
    char *dirsep = strrchr(argv[0], '/');
    if ( dirsep != NULL ) *dirsep = 0;
    int changeDir = chdir(argv[0]);
    std::cout << "Changing Dir: " << changeDir << "\n";
    std::string cwd = getcwd(NULL, 0);
    std::cout << cwd << "\n";
    std::cout << "Starting simulation\n";
#endif
}