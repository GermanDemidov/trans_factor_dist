//
//  Transcription_factors_in_cell.h
//  trans_factors_distrib
//
//  Created by Герман Демидов on 24.04.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#ifndef __trans_factors_distrib__Transcription_factors_in_cell__
#define __trans_factors_distrib__Transcription_factors_in_cell__

#include <stdio.h>
#include <vector>
#include "transcription_factor.h"
#include "auxuilary.h"
#include <set>

class Transcription_factors_in_cell {
public:
    Transcription_factor& get_tf_by_index(int);

    Transcription_factors_in_cell(std::vector<std::string>, std::vector<Transcription_factor>&);
    int choose_next_unbinded_DNA_to_interact();
    int choose_next_binded_DNA_to_interact();
    int number_of_binded_dna();
    double generate_next_characteristic_time_for_binded_TFs();
    double generate_next_characteristic_time_for_unbinded_TFs();
    std::string get_type_of_TF(int i);
    int get_size_of_TF(int i);
    
    // change options of binding
    void bind_tf_to_dna(int, bool);
    void unbind_tf_from_dna(int);
    
    int num_of_binded_TFs();
    int num_of_unbinded_TFs();

private:
    std::vector<Transcription_factor> tfs;
    std::vector<std::string> protein_names;
    std::map<std::string, int> number_of_unbinded_TFs;
    std::set<int> indexes_of_binded_TFs;
    Choose_element_from_array<int> choose_for_int_vect;
    

};

#endif /* defined(__trans_factors_distrib__Transcription_factors_in_cell__) */
