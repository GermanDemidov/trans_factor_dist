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
    
    void start_simulation(Transcription_factors_in_cell& tfs, Parser_dnase_acc&);
    
private:
    int cell_id;
    int total_number_of_TFs_to_bind = 100;
    int number_of_steps_before_stabilization = 100;
    int number_of_one_dim_slidings_in_step = 3;
    
    // first - time, second - true if it is operation of binding, false if it is operation of sliding or unbinding
    std::priority_queue<std::pair<double, bool>, std::vector<std::pair<double, bool>>, compare> timeline;
    
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

    // vector with occupied coordinates
    std::set<std::pair<int, int> > no_access_dna_forward;
    std::set<std::pair<int, int> > no_access_dna_revcompl;
    
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
    void bind_tf_to_dna(int, Transcription_factors_in_cell&, Parser_dnase_acc& dnase);
    void unbind_tf_from_dna(int, Transcription_factors_in_cell&);
    
    // calculate next movement
    int move_tf_right_or_left(int, double, double, std::string& );
    
    // test for unbinding
    bool test_for_unbinding(double );

    // test to not to fall in no access dna
    bool binding_site_is_free(bool forward_or_reverse, int coord, int len_of_tf);

    // the function that denotes one dimensional sliding
    void one_dimensional_slinding_of_TF(int i, Transcription_factors_in_cell&);
    
    // function that finds if tf is binded specifically and apply the result
    void test_for_specific_binding(int i, Transcription_factors_in_cell&);
    
    

};




#endif /* defined(__trans_factors_distrib__cell__) */
