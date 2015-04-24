//
//  parser.h
//  trans_factors_distrib
//
//  Created by Герман Демидов on 10.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#ifndef __trans_factors_distrib__parser__
#define __trans_factors_distrib__parser__

#include <stdio.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class Parser_fasta {
    
private:
    std::map<std::string, std::string> sequences;
    
public:
    Parser_fasta(std::string);
    std::string& get_sequence_by_name(std::string);
    std::map<std::string, std::string> return_sequences();
    
    void debug_print();
    
};

class Parser_pcm {
    
private:
    std::vector<char> nucleotides;
    std::map<std::string, std::vector<std::map<char, double>>> pcm;
    std::map<std::string, std::vector<std::map<char, double>>> pwm;
    std::vector<std::string> protein_names;

    
public:
    Parser_pcm(std::string);
    
    void calculate_pwm(std::map<char, double>&);
    std::vector<std::map<char, double>> return_profile_for_protein_pcm (std::string);
    std::vector<std::map<char, double>> return_profile_for_protein_pwm (std::string);
    std::vector<std::string> return_protein_names();

};

class Parser_dnase_acc {
private:
    std::vector<std::pair<int, int>> open_acc_intervals;

public:
    Parser_dnase_acc(std::string);
    bool is_in_interval(std::pair<int, int>);
    bool is_in_interval(int);
};



#endif /* defined(__trans_factors_distrib__parser__) */
