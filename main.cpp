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
            reproProb = exp(-10.0 * (R_r - Resource));
        }
        totalProb += reproProb;
    }

    void uptake(double& R, double H, double alpha) {
        double uptake = alpha * R/(R+H);
        if(uptake > R) uptake = R;
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

    std::vector<double> growthTimes;

    Host_infected() {
        Resource = 0.0;
        numberPhages = 0;
        lysisTime = 1e6;
    }

    Host_infected& operator=(const Host_infected& other) {
        Resource = other.Resource;
        numberPhages = other.numberPhages;
        lysisTime = other.lysisTime;
        growthTimes = other.growthTimes;
        return *this;
    }
};



template <typename T>
int pickEvent(const std::vector<T>& v, const T& sum) {
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

template <typename T>
double calcMean(const std::vector<T>& v)
{
    if(v.empty()) {
        return 0;
    }

    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    return(mean);
}

void doSimulation(const GetParams& Params, int repl)
{
    std::vector<Host> population;
    std::vector<Host_infected> population_infected;

    double totalReproProb = 0.0;
    double maxReproProb = 0.0;




    //double R_star = Params.D * Params.H / (Params.alpha - Params.D);
   // double initCalcHostSize = (Params.Inflow - R_star) * (R_star + Params.H) * (Params.alpha - Params.D) / (Params.alpha * Params.H);

    for(int i = 0; i < Params.initHostPopSize; ++i) {
        Host temp = Host();
        temp.Resource = Params.reproductionSize / 2;
        temp.reproSize = Params.reproductionSize;

        temp.updateReproProb(totalReproProb, Params.reproductionSize);

        if(temp.reproProb > maxReproProb) {
            maxReproProb = temp.reproProb;
        }
    //    totalReproProb += temp.reproProb;
        population.push_back(temp);
    }

    double t = 0.0;
    double R = Params.Inflow * Params.D;


    //std::ofstream outFile("outFile.txt",std::ios::app);
    std::ofstream outFile("outFile.txt", std::ios::app);

    double prevT = 0.0;

    double avgR = 0.0;
    for(auto it = population.begin(); it != population.end(); ++it) {
        avgR += (*it).Resource;
    }
    avgR = avgR / population.size();

    //   outFile << t << "\t" << population.size() << "\t" << R << "\t" << avgR << "\n";

    static int counter = 1;
    counter++;

    int numPhages = 0;
    int burstCounter = 0;
    int summedBurstSize = 0;
    double doublingTime = -1;



    while(t < Params.maxTime) {
        long double maintenanceRate = population.size() * Params.maintenanceCost;
        long double uptakeRate = population.size() * Params.alpha * R / (R+Params.H);
        long double reproRate  = totalReproProb;

        long double infectRate = population.size() * numPhages * Params.omega;
        long double phageDecay = numPhages * Params.gamma;
        long double Inf_Hosts_Rate = (int)population_infected.size() * Params.beta; //beta is the rate of phage production, e.g. 1/beta is the average time until one phage is produced
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
            return;
        //    break;
        }

        int allHosts = (int)population.size() + (int)population_infected.size();
        if(lossProb > 0 && allHosts > 0 ) {
            int loss = Binom(allHosts, lossProb);
            for(int i = 0; i < loss; ++i) {
                std::vector<int> probs = {(int)population.size(), (int)population_infected.size() };
                int event = pickEvent(probs, allHosts);
                switch(event) {
                    case 0: {
                        if(!population.empty()) {
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
                        break;
                    }
                    case 1: {
                        if(!population_infected.empty()) {
                            int index = (int)random_number((int)population_infected.size());
                            population_infected[index] = population_infected.back();
                            population_infected.pop_back();
                        }
                        break;
                    }
                }
            }
        }

        if(numPhages > 0 && lossProb > 0) {
            int loss = Binom((int)numPhages, lossProb);
            numPhages -= loss;
            if(numPhages < 0) {
                numPhages = 0;
            }
        }

        if(t > Params.infectionTime || Params.infectionTime > Params.maxTime) {
            if(doublingTime > 0 && burstCounter > 0) {
                break;
            }

            if(population.size() <= 0.0) {
               break;
            }
        }


        switch (event)
        {
            case 0: //maintenance
            {
                if(population.empty()) break;
                int index = random_number((int)population.size());
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
                if(population.empty()) break;
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
                if(population.empty()) break;
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
                if(population.empty()) break;
                int index = random_number((int)population.size());
                Host_infected new_infected_host;

                new_infected_host.Resource = population[index].Resource;
                new_infected_host.growthTimes.push_back(t); //initial time
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
                        population_infected[index].growthTimes.push_back(t);
                    } else {

                        //lysis!!!
                        numPhages += population_infected[index].numberPhages;
                        burstCounter++;
                        summedBurstSize += population_infected[index].numberPhages;

                        /*std::ofstream timeFile("growthTimes.txt",std::ios::app);
                        for(auto it = population_infected[index].growthTimes.begin(); it != population_infected[index].growthTimes.end(); ++it) {
                            timeFile << (*it) << "\t";
                        }
                        timeFile << "\n";
                        timeFile.close();*/

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
                        burstCounter++;
                        summedBurstSize += population_infected[index].numberPhages;

                       /* std::ofstream timeFile("growthTimes.txt",std::ios::app);
                        for(auto it = population_infected[index].growthTimes.begin(); it != population_infected[index].growthTimes.end(); ++it) {
                            timeFile << (*it) << "\t";
                        }
                        timeFile << "\n";
                        timeFile.close();*/

                        population_infected[index] = population_infected.back();
                        population_infected.pop_back();
                    }
                }
                break;
            }
        }


        t += dt;
        if(t >= Params.infectionTime && (t-dt) < Params.infectionTime) {
            numPhages = Params.initPhagePopSize;
        }

        double currentT = t;
        if(currentT - prevT >= 0.001) {

            if(doublingTime < 0) {
                if(population.size() >= 4 * Params.initHostPopSize) {
                    doublingTime = t;
                }
            }

            prevT = currentT;
            avgR = 0.0;
            for(auto it = population.begin(); it != population.end(); ++it) {
                avgR += (*it).Resource;
            }
            avgR = avgR / population.size();

           /* std::ofstream trackFile("trackThings.txt",std::ios::app);

            trackFile << repl << "\t" << Params.Inflow << "\t" << Params.alpha << "\t" << Params.maintenanceCost << "\t" <<
                        population.size() << "\t" << numPhages << "\t" << R << "\t" << avgR << "\t" << t << "\n";
            trackFile.close();*/
        }
    }


    avgR = 0.0;
    for(auto it = population.begin(); it != population.end(); ++it) {
        avgR += (*it).Resource;
    }
    avgR = avgR / population.size();

    outFile  << repl << "\t" << Params.Inflow << "\t" << Params.alpha << "\t" << Params.maintenanceCost << "\t" <<
    population.size() << "\t" << numPhages << "\t" << 1.0 * summedBurstSize / burstCounter <<  "\t" << doublingTime << "\t" << avgR << "\n";

    return;
}


int main(int argc, const char * argv[]) {

    macstart(argv);
    GetParams P;
    P.readFromIni("config.ini");
    set_seed(P.seed);

    /*
    clock_t start2 = clock();
    doSimulation(P,0);
    clock_t end2 = clock();
    double elapsed_secs2 = double(end2 - start2) / CLOCKS_PER_SEC;
    std::cout << "that took: " << elapsed_secs2 << " seconds\n";
*/



    for(double alpha = 1.0; alpha <= 1.5; alpha+= 0.1) {
            P.alpha = alpha;

       // for(double c = 0.0; c < 0.1; c+= 0.01) {
       //         P.maintenanceCost = c;

            for(int r = 0; r < 5; ++r) {
                std::cout << P.Inflow << "\t" << P.alpha << "\t" << P.maintenanceCost << "\t" << r << "\n";
                doSimulation(P,r);
            }
        //}
     }

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
