//
//  GetParams.h
//  Adapted from Hanno Hildenbrandt
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2016 Thijs Janzen. All rights reserved.
//

#ifndef FINITE_CHROMOSOME_ALWAYS_RECOM_GETPARAMS_H_
#define FINITE_CHROMOSOME_ALWAYS_RECOM_GETPARAMS_H_

#include <sstream>
#include <string>
#include <vector>


class GetParams {
 public:
    GetParams();
    void readFromIni(const char * filename);

    template <typename T>
    void readNameValuePair(std::stringstream& ss,
                           std::string iniName, T& value);

    // GLOBAL VARIABLES
    int seed;
    int initHostPopSize;
    // size of host at which it reproduces
    double reproductionSize;
    double maintenanceCost;
    double phageGrowthRate;
    double infectionProbability;
    double phageDecayRate;

    // influx of the chemostat
    double Inflow;
    double Resources;

    int maxTime;
    // time at which the phages are introduced
    int infectionTime;
};

#endif  // FINITE_CHROMOSOME_ALWAYS_RECOM_GETPARAMS_H_
