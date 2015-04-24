//
//  cell.h
//  trans_factors_distrib
//
//  Created by Герман Демидов on 11.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#ifndef __trans_factors_distrib__cell__
#define __trans_factors_distrib__cell__

#include <stdio.h>
#include <map>
#include <vector>
#include <utility>
#include <deque>
#include "auxuilary.h"
#include "transcription_factor.h"
#include "parser.h"
#include <ctime>
#include "transcription_factors_in_cell.h"
#include <queue>

class Cell {
public:
    Cell(std::string&, std::string&, int, std::vector<Transcription_factor>&, Parser_dnase_acc&, std::map<std::string, double>&,
         std::vector<Transcription_factor>&);
    void generate_next_event(Transcription_factors_in_cell&);
    void generate_TF_appearance_time();
    void find_specific_binding_sites(Transcription_factor, Parser_dnase_acc&, bool);
    
    void find_potential_strength(bool);
    
    void start_simulation();
    
private:
    int cell_id;
    const int total_number_of_TFs_to_bind = 100;
    
    // first - time, second - true if it is operation of binding, false if it is operation of sliding or unbinding
    std::priority_queue<std::pair<double, bool>> timeline;
    
    std::string forward_DNA;
    std::string reverse_DNA;
    std::vector<Transcription_factor> transcription_factors;
    // pre-calculated weights of binding
    std::map<std::string, std::vector<double> > weights_of_binding_of_all_tfs_forward;
    std::map<std::string, std::vector<double> > weights_of_binding_of_all_tfs_revcompl;
    
    // maps that indicate the specific sites occupancy
    std::map<std::string, std::map<int, bool> > specific_binding_sites_forward;
    std::map<std::string, std::map<int, bool> > specific_binding_sites_revcompl;
    
    // potential in pairs, 5 bp backward, 5 bp forward (just mean of binding strenths)
    std::map<std::string, std::vector<std::pair<double, double>> > potential_strength_forward;
    std::map<std::string, std::vector<std::pair<double, double>> > potential_strength_revcompl;

    // return special value if the index is out of vector coordinates
    double return_weight_of_binding(int, std::vector<double>);
    
    // deque for the appearance of TFs, their times (without specifying the type of TF)
    //std::deque<double> TF_appearance_time;
    
    // map with concentrations for each TF, that used as propensities
    std::map<std::string, double> concentrations;
    // pre-calculated sum of propensities
    double sum_of_concentrations = 0.0;
    
    // minimum of potentials for each TF
    std::map<std::string, double> minimums_for_normalizations_of_potentials;
    
    // events of binding and unbinding
    void bind_tf_to_dna(int, Transcription_factors_in_cell& );
    void unbind_tf_from_dna(int, Transcription_factors_in_cell& );
    
    // calculate next movement
    int move_tf_right_or_left(int, double, double, std::string& );
    
    // test for unbinding
    bool test_for_unbinding(double );

};

#endif /* defined(__trans_factors_distrib__cell__) */
