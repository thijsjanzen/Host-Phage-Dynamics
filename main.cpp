//
//  main.cpp
//  Host-Phage
//
//  Created by Thijs Janzen on 30/05/16.
//  Copyright (c) 2016 Thijs Janzen. All rights reserved.
//

#include "./randomc.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <unistd.h>



template <class T>
double sum(const std::vector<T>& v) {
    double a = std::accumulate(v.begin(), v.end(), 0.0);
    return a;
}

class Host {
 public:
    double reproduction_limit;
    double R;  // resources
    double P;  // number of phages residing inside the host
    double lysisTimeOfPhage;  // type of phage that infected the host
    double remainingLysisTime;  // remaining time until lysis
    bool dead;
    bool infected;

    void updatePhages(double PhageGrowthRate,
                      std::vector< int >* Phages,
                      int* numberOfPhages);
    void maintain(double maintenance);
    void infection(std::vector< int >* Phages,
                   int numberOfHosts,
                   double infectProbability,
                   int* numberOfPhages);

    Host()  {
        R = 0;
        P = 0;
        remainingLysisTime = 1e6;
        lysisTimeOfPhage = -1;
        dead = false;
        infected = false;
    }

    Host(double initR, double maxR):
        R(initR), reproduction_limit(maxR) {
        P = 0;
        remainingLysisTime = 1e6;
        lysisTimeOfPhage = -1;
        dead = false;
        infected = false;
    }

    Host& operator=(const Host& other)  {
        R = other.R;
        P = other.P;
        remainingLysisTime = other.remainingLysisTime;
        lysisTimeOfPhage = other.lysisTimeOfPhage;
        dead = other.dead;
        infected = other.infected;
        reproduction_limit = other.reproduction_limit;
        return *this;
    }

    Host(const Host& other) {
        R = other.R;
        P = other.P;
        remainingLysisTime = other.remainingLysisTime;
        lysisTimeOfPhage = other.lysisTimeOfPhage;
        dead = other.dead;
        infected = other.infected;
        reproduction_limit = other.reproduction_limit;
    }
};

void Host::updatePhages(double PhageGrowthRate,
                        std::vector< int >* Phages,
                        int* numberOfPhages)    {
    // if there are no phages, don't update them
    if (!infected) return;

    remainingLysisTime--;
    if (R >= PhageGrowthRate)    {
        P += PhageGrowthRate;
        R -= PhageGrowthRate;
    }

    if (remainingLysisTime <= 0) {
        // lysis!
        dead = true;
        for (int i = 0; i < static_cast<int>(P); ++i) {
            int index = lysisTimeOfPhage;
            if (uniform() < 0.01) {
                if (uniform() < 0.5) {
                    index++;
                } else {
                    index--;
                }
            }
            (*Phages)[index]++;
            (*numberOfPhages)++;
        }
    }
    return;
}

void Host::maintain(double maintenance) {
    // float precision should be enough and is faster
    double loss = maintenance *  powf(static_cast<float>(1+R), 0.75f);
    R -= loss;
    if (R <= 0.0) {
        dead = true;
    }
    return;
}

void Host::infection(std::vector< int >* Phages,
                     int numberOfHosts,
                     double infectProbability,
                     int* numberOfPhages)  {
    // only non-infected individuals can become infected
    if (infected) return;
    // dead individuals can't become infected
    if (dead) return;

    double probInfection = 1.0 * (*numberOfPhages) *
                                 numberOfHosts *
                                 infectProbability;

    int lysisTime = -1;
    if (uniform() < probInfection)  {
        int r = 1;
        // we don't need to draw a random number
        // when there is a single phage left
        if (*numberOfPhages > 1) {
            r  += random_number(*numberOfPhages);
        }

        for (int i = 0; i < Phages->size(); ++i) {
            r -= (*Phages)[i];
            if (r <= 0) {
                lysisTime = i;
                break;
            }
        }
    }

    if (lysisTime > 0) {
        infected = true;
        lysisTimeOfPhage = lysisTime;
        remainingLysisTime = lysisTime;
        (*Phages)[lysisTime]--;
        numberOfPhages--;
    }
}

void update(int HostPopulationSize,
            const std::vector< Host >::iterator& it,
            std::vector< Host >* toAdd,
            double* Resources,
            double maintenance,
            double infectProbability,
            std::vector< int >* Phages,
            double PhageGrowthRate,
            int* numberOfPhages)    {
    if ((*it).dead) return;


    (*it).maintain(maintenance);
    if ((*it).dead) return;

    if ((*it).infected) {
        (*it).updatePhages(PhageGrowthRate, Phages, numberOfPhages);
        if ((*it).dead) return;
    }

    // also check if dead because he could have died
    // in the function maintain (yes, this is a bit silly)
    if (!(*it).infected && !(*it).dead) {
        double uptake   = 1.0 * (*Resources) / (*Resources + HostPopulationSize);
        (*it).R         += uptake;
        *Resources       -= uptake;

        // double prob = 1.0 * (*it).R/ (*it).reproduction_limit;
        // prob *= prob;
        double prob = expf(-0.5 * ((*it).reproduction_limit - (*it).R));
        if (uniform() < prob)    {
            (*it).R *= 0.5;
            Host copy = (*it);
            if (uniform() < 0.05) {
                copy.reproduction_limit += normal(0.0, 1);
                if (copy.reproduction_limit < 1) copy.reproduction_limit = 1;
            }
            toAdd->push_back(copy);
        }

        (*it).infection(Phages,
                        HostPopulationSize,
                        infectProbability,
                        numberOfPhages);
    }
    return;
}

// helper function to remove dead individuals from the population
bool isDead(const Host& H) {
    return H.dead;
}

void removeDeadAddKids(std::vector< Host >* HostPopulation,
                       std::vector< Host >* toAdd) {
    // now we first have to remove all the dead individuals
    // and then add the born individuals
    auto it = std::remove_if(HostPopulation->begin(), HostPopulation->end(),
                             isDead);
    HostPopulation->erase(it, HostPopulation->end());

    if (!toAdd->empty()) {
        HostPopulation->insert(HostPopulation->end(),
                               toAdd->begin(), toAdd->end());
    }
}

void decayPhages(std::vector< int >* Phages,
                 const double& phageDecay,
                 int* numberOfPhages) {
    for (auto it = Phages->begin(); it != Phages->end(); ++it) {
        if ((*it) >= 1) {
            int loss = Binom((*it), phageDecay);
            (*it) -= loss;
            (*numberOfPhages) -= loss;
            if ((*it) < 0) (*it) = 0;
        }
    }
}


void macstart(const char * argv[]);
void writeOutput(const std::vector< Host >& HostPopulation,
                 const std::vector< int > Phages,
                 int t,
                 int numberOfPhages);


int main(int argc, const char * argv[]) {
    // apple specific code making sure that
    // we are in the same directory as the executable
    macstart(argv);

    // initial host population size
    int HostN = 1e2;
    // size of host at which it reproduces
    double maxR = 20;
    double maintenance = 0.1;
    double PhageGrowthRate = 1;
    double infectionProbability = 1e-9;
    double phageDecay = 0.1;

    // influx of the chemostat
    double Inflow = 1e3;
    double Resources = 0;

    int maxT = 1e7;
    // time at which the phages are introduced
    int critT = 1e9;

    std::random_device rdev;
    unsigned int chosen_seed = rdev();
    // we only use the "standard" random number generator for std::shuffle
    std::mt19937 urng(chosen_seed);

    // all other random numbers are generated
    // using code from Agner Fog:  www.agner.org/random
    set_seed(chosen_seed);

    // population of hosts
    std::vector< Host > HostPopulation;

    // temporary vector containing newly born hosts
    std::vector< Host > toAdd;

    // vector of the histogram of Phages
    std::vector < int > Phages(maxR*4, 0);

    // we initialize a population of N individuals
    // that each have a minimum amount of resources already
    HostPopulation.resize(HostN, Host(maintenance * 10, maxR));

    // clear the output files
    std::ofstream outFile_temp("output.txt");
    outFile_temp.close();

    std::ofstream outFile_temp3("output_phages.txt");
    outFile_temp3.close();

    std::clock_t start = clock();

    int numberOfPhages = 0;
    // the maximum time is quite arbitrary
    for (int t = 0; t < maxT; t++) {
        // add resources
        Resources += Inflow;
        toAdd.clear();
        // allocate plenty of memory
        // (this vector will hold all new additions to the population)
        toAdd.reserve(HostPopulation.size());

        // shuffle all individuals
        std::shuffle(HostPopulation.begin(), HostPopulation.end(), urng);

        for (auto it = HostPopulation.begin();
                  it != HostPopulation.end(); ++it) {
            update(static_cast<int>(HostPopulation.size()),
                   it, &toAdd, &Resources,
                   maintenance, infectionProbability,
                   &Phages, PhageGrowthRate, &numberOfPhages);
        }

        removeDeadAddKids(&HostPopulation, &toAdd);

        decayPhages(&Phages, phageDecay, &numberOfPhages);

        if (t == critT) {
            int maxMaxr = 0;
            for (std::size_t i = 0; i < HostPopulation.size(); ++i) {
                if (HostPopulation[i].reproduction_limit > maxMaxr)
                    maxMaxr = HostPopulation[i].reproduction_limit;
            }

            for (int i = 1; i < maxMaxr; ++i)   {
                int relNum = static_cast<int>(HostPopulation.size() / maxMaxr);
                Phages[i] += relNum;
                numberOfPhages += relNum;
            }
        }

        if (t % 10 == 0) writeOutput(HostPopulation, Phages, t, numberOfPhages);

        if (HostPopulation.size() < 1) {
            break;
        }
    }

    HostPopulation.clear();
    Phages.clear();

    std::clock_t end = clock();

    std::cout << "Time taken: " << "\t"
              << static_cast<double>((end - start)/CLOCKS_PER_SEC)
              << "\n";

    return 0;
}

void writeOutput(const std::vector< Host >& HostPopulation,
                 const std::vector< int > Phages,
                 int t,
                 int numberOfPhages) {
    int numberOfHosts = static_cast<int>(HostPopulation.size());

    double meanResources = 0.0;
    double meanMaxR = 0.0;
    double minMaxr = 1e6;
    double maxMaxr = -1;

    for (auto it = HostPopulation.begin(); it != HostPopulation.end(); ++it) {
        meanResources += (*it).R;
        meanMaxR += (*it).reproduction_limit;
        if ((*it).reproduction_limit > maxMaxr)
            maxMaxr = (*it).reproduction_limit;
        if ((*it).reproduction_limit < minMaxr)
            minMaxr = (*it).reproduction_limit;
    }
    meanResources = 1.0 * meanResources / HostPopulation.size();
    meanMaxR = 1.0 * meanMaxR / HostPopulation.size();

    std::cout << t << "\t"
              << numberOfHosts << "\t"
              << numberOfPhages << "\t"
              << meanMaxR << "\n";

    std::ofstream outFile;
    outFile.open("output.txt", std::ios::app);
    if (outFile.is_open())   {
        outFile << t << "\t" << numberOfHosts << "\t"
                             << meanResources << "\t"
                             << numberOfPhages << "\t"
                             << meanMaxR << "\t"
                             << minMaxr << "\t"
                             << maxMaxr << "\n";
    }
    outFile.close();

    std::ofstream outFile_Phages("output_phages.txt", std::ios::app);
    outFile_Phages << t << "\t";
    for (int i = 0; i < Phages.size(); ++i) {
        outFile_Phages << Phages[i] << "\t";
    }
    outFile_Phages << "\n";
    outFile_Phages.close();
}

void macstart(const char * argv[])  {
    std::cout << "\n\n\n";
#ifdef __APPLE__
    {
        char *dirsep = strrchr(argv[0], '/');
        if (dirsep != NULL) *dirsep = 0;
        int changeDir = chdir(argv[0]);
        std::cout << "Changing Dir: " << changeDir << "\n";
        std::string cwd = getcwd(NULL, 0);
        std::cout << cwd << "\n";
        std::cout << "Starting simulation\n";
    }
#endif
}
