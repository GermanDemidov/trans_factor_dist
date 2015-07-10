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
#include <fstream>
#include <iostream>
#include <queue>

class Cell {
public:
    Cell(std::string&, std::string&, int, std::vector<Transcription_factor>&, Parser_dnase_acc&, std::map<std::string, double>&,
         std::vector<Transcription_factor>&);
    void generate_next_event(Transcription_factors_in_cell&);
    void find_specific_binding_sites(Transcription_factor, Parser_dnase_acc&, bool);
    
    void find_potential_strength(bool);
    
    void start_simulation(Transcription_factors_in_cell tfs, Parser_dnase_acc);
    
    std::map<std::vector<bool>, bool> get_final_frequency_of_combinations();
    
    void null_everything(Parser_dnase_acc&);
    
    std::vector<bool> get_final_combo();
    std::vector<double> get_times_of_changes();
    
private:
    int cell_id;
    int total_number_of_TFs_to_bind;
    int number_of_steps_before_stabilization = 5000;
    int number_of_one_dim_slidings_in_step = 5;
    std::string repressor = ">hb";
    int length_of_repression = 100;
    
    void generate_next_event_with_binded(Transcription_factors_in_cell&);
    void generate_next_event_with_unbinded(Transcription_factors_in_cell&);
    
    std::map<std::string, int> sizes_of_tfs;
    void determine_sizes_of_tfs(Transcription_factor, Parser_dnase_acc&, bool);
    
    double time_before_stabilization = 5000.0;
    double upper_bound_for_time = 15000.0;
    int bound_for_number_of_specific_sites = 60;


    // first - time, second - 1 if it is operation with unbinded, 2 if operation with binded
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, compare> timeline;
    
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
    
    // map with concentrations for each TF, that used as propensities
    std::map<std::string, double> concentrations;
    // pre-calculated sum of propensities
    double sum_of_concentrations = 0.0;
    
    // minimum of potentials for each TF
    std::map<std::string, double> minimums_for_normalizations_of_potentials;
    
    // events of binding and unbinding
    bool bind_tf_to_dna(int, Transcription_factors_in_cell&, Parser_dnase_acc& dnase);
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
    bool test_for_specific_binding(int i, Transcription_factors_in_cell&);
    
    // output of the combinations of specific binding sites occupancies
    int number_of_specific_binding_sites = 0;
    std::map<std::string, std::map<int, int>> codes_of_specific_binding_sites_forward;
    std::map<std::string, std::map<int, int>> codes_of_specific_binding_sites_revcompl;
    std::vector<bool> specific_sites_binded_or_not;
    std::map<std::vector<bool>, bool> final_frequency_of_combinations;
    std::vector<double> time_for_unbinded_state;
    std::vector<double> time_for_binded_state;
    
    // counter of steps without any changes
    int counter_of_steps_without_changes = 0;
    
    std::vector<double> times_of_changes;

};




#endif /* defined(__trans_factors_distrib__cell__) */
