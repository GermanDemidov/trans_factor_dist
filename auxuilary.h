//
//  auxuilary.h
//  trans_factors_distrib
//
//  Created by Герман Демидов on 22.04.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#ifndef __trans_factors_distrib__auxuilary__
#define __trans_factors_distrib__auxuilary__

#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <random>
#include <iostream>


std::string reverse_compliment(std::string);

double sm(std::vector<double>);
double sd(std::vector<double>);

double generate_next_time(double);
int generate_next_coordinate_change();

#endif /* defined(__trans_factors_distrib__auxuilary__) */
