//
//  auxuilary.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 22.04.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "auxuilary.h"

Precalculated_values::Precalculated_values(double min_value_exp, double max_value_exp, int precision_exp) {
    min_val_to_add = min_value_exp;
    int multiplier_to_round_up = 1;
    for (int i = 0; i < precision_exp; ++i)
    {
        multiplier_to_round_up *= 10;
    }
    multiplier_to_round = multiplier_to_round_up;
    calculate_vector_of_exponents(min_value_exp, max_value_exp);
}

void Precalculated_values::calculate_vector_of_exponents(double min_val, double max_val) {
    int min_val_to_calc = (int)(min_val * multiplier_to_round);
    int max_val_to_calc = (int)(max_val * multiplier_to_round);
    for (int i = min_val_to_calc; i <= max_val_to_calc; i++) {
        double x = (1.0 * i) / multiplier_to_round;
        exponents_values.push_back(exp(x));
    }
};

double Precalculated_values::return_exponent_of_double(double x) {
    return exponents_values[(int)(abs(min_val_to_add) + x * multiplier_to_round)];
};

double sm(std::vector<double> sample) {
    double mean = 0.0;
    for (int i = 0; i < sample.size(); ++i) {
        mean += sample[i];
    }
    return mean / sample.size();
}

double sd(std::vector<double> sample) {
    double var = 0.0;
    double sum_of_squares = 0.0;
    double sum = 0.0;
    for (int i = 0; i < sample.size(); ++i) {
        sum += sample[i];
        sum_of_squares += sample[i] * sample[i];
    }
    var += sum_of_squares;
    var -= ((sum * sum) / sample.size());
    return sqrt(var / (sample.size() - 1));
}

std::string reverse_compliment(std::string DNA) {
    std::string reversed_and_complimented_dna = "";
    for (long i = DNA.size() - 1; i != -1; --i) {
        if (DNA.at(i) == 'A') reversed_and_complimented_dna += 'T';
        if (DNA.at(i) == 'T') reversed_and_complimented_dna += 'A';
        if (DNA.at(i) == 'G') reversed_and_complimented_dna += 'C';
        if (DNA.at(i) == 'C') reversed_and_complimented_dna += 'G';
    }
    return reversed_and_complimented_dna;
}

double generate_next_time(double sum_of_propensities) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    double lambda = 0.000005 * sum_of_propensities;
    double result = 0.0;
    std::exponential_distribution<double> distribution(lambda);
    result = distribution(gen);

    return result;
}

double fRand(double x, double y)
{
    double fMin = std::min(x,y);
    double fMax = std::max(x,y);
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


std::default_random_engine & global_urng( ) {
    static std::default_random_engine  u{};
    return u;
}
void randomize( ) {
    static std::random_device  rd{};
    global_urng().seed( rd() );
}

int pick_a_number( int from, int thru ) {
    randomize();
    static std::uniform_int_distribution<>  d{};
    using  parm_t  = decltype(d)::param_type;
    return d( global_urng(), parm_t{from, thru} );
}

std::pair<int, int> reverse_interval_for_another_strand(long length, std::pair<int, int> coords) {
    return std::make_pair(length - coords.second, length - coords.first);
}

std::pair<int, int> broad_borders(int broad_to, std::pair<int, int> coords) {
    return std::make_pair(coords.first - broad_to, coords.second + broad_to);    
}

int return_normal_number(double mean, double sd) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean, sd);
    double number = distribution(generator);
    return (int)number;
}





