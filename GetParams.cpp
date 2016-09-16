//
//  GetParams.cpp
//  Adapted from Hanno Hildenbrandt, 2007
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2016 Thijs Janzen. All rights reserved.
//



#include "GetParams.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

GetParams::GetParams() {
    seed = 1;
    initHostPopSize = 1e2;
    initPhagePopSize = 0;
    // size of host at which it reproduces
    reproductionSize = 20;
    maintenanceCost = 0.1;

    // influx of the chemostat
    Inflow = 1e3;
    Resources = 0;

    maxTime = 1e7;
    // time at which the phages are introduced
    infectionTime = 1e9;
}

void GetParams::readFromIni(const char * filename ) {
    // locate file and tranfer text to stringstream
    std::ifstream ifs(filename);
    std::stringstream ss;

    // only for succesfully created ifstream: otherwise null-pointer?
    if (ifs) {
        // config.ini content is transferred to stringstream: easier to search?
        ss << ifs.rdbuf();
    } else {
        throw "Can't locate file";
    }

    while (ss.good()) {
        readNameValuePair(ss,  "seed", seed);
        readNameValuePair(ss,  "initHostPopSize", initHostPopSize);
        readNameValuePair(ss,  "initPhagePopSize", initPhagePopSize);
        readNameValuePair(ss,  "reproductionSize", reproductionSize);
        readNameValuePair(ss,  "maintenanceCost", maintenanceCost);
        readNameValuePair(ss,  "H", H);
        readNameValuePair(ss,  "Inflow", Inflow);
        readNameValuePair(ss,  "maxTime", maxTime);
        readNameValuePair(ss,  "infectionTime", infectionTime);
        readNameValuePair(ss,  "D",             D);
        readNameValuePair(ss,  "omega", omega);
        readNameValuePair(ss,  "gamma", gamma);
        readNameValuePair(ss,  "lambda", lambda);
        readNameValuePair(ss,  "beta", beta);
        readNameValuePair(ss,  "alpha", alpha);

    }
}

template <typename T>
void GetParams::readNameValuePair(std::stringstream& ss,
                                  std::string iniName, T& value ) {
    std::string name;
    char sign;

    // >> copy ss content to string until white space is encountered
    ss >> name;
    if (name != iniName )
        throw "expect parameter";
    ss >> sign;  // copies ss content to character
    if (sign != '=' )
        throw "text format of ini file is not compatible";
    ss >> value;
    std::cout << iniName << ": " << value << std::endl;
}
