//
//  auxuilary.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 22.04.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "auxuilary.h"

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
    
    double lambda = 1 / sum_of_propensities;
    double result = 0.0;
    std::exponential_distribution<double> distribution(lambda);
    result = distribution(gen);

    return result;
}

int generate_next_coordinate_change() {
    return 5 + (int)generate_next_time(4.0);
}
