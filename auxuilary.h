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
#include <map>
#include <tgmath.h>
#include <algorithm>
#include  <iterator>

class Precalculated_values {
public:
    Precalculated_values(double, double, int);
    double return_exponent_of_double(double);
    
private:
    double min_val_to_add = 0.0;
    void calculate_vector_of_exponents(double, double);
    std::vector<double> exponents_values;
    int multiplier_to_round = 1;

};


std::string reverse_compliment(std::string);

double sm(std::vector<double>);
double sd(std::vector<double>);

double generate_next_time(double);
int generate_next_coordinate_change();

double calculate_vector_of_exponents(double i, double j);

double fRand(double, double);

template<class T>
class Choose_element_from_array {
public:
    Choose_element_from_array() {
        auto const seed = std::random_device()();
        std::mt19937 eng(seed);
        engine = eng;
    }

private:
    std::mt19937 engine;
};

struct compare
{
    inline bool operator () (const std::pair<double, bool>& x, const std::pair<double, bool>& y) const
    {
        if(x.first != y.first)
            return x.first > y.first;
        else
            return x.second > y.second;
    }
};


std::default_random_engine & global_urng( );
void randomize( );
int pick_a_number( int from, int thru );

std::pair<int, int> reverse_interval_for_another_strand(long, std::pair<int, int>);

std::pair<int, int> broad_borders(int, std::pair<int, int>);

int return_normal_number(double mean, double sd);

#endif /* defined(__trans_factors_distrib__auxuilary__) */
