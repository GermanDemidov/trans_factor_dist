//
//  cell.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 11.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "cell.h"

Cell::Cell(std::string& for_DNA, std::string& rev_DNA,
           int id_of_cell,
           std::vector<Transcription_factor>& types_of_TFs_for_preprocessing, Parser_dnase_acc& dnase, std::map<std::string, double>& conc, std::vector<Transcription_factor>& transcription_factor_instances){
    std::vector<std::string> protein_names;
    for (int i = 0; i < types_of_TFs_for_preprocessing.size(); ++i) {
        protein_names.push_back(types_of_TFs_for_preprocessing[i].type_name);
    }
    Transcription_factors_in_cell tfs(protein_names, transcription_factor_instances);
    forward_DNA = for_DNA;
    reverse_DNA = rev_DNA;
    cell_id = id_of_cell;
    concentrations = conc;
    
    for (std::map<std::string, double>::iterator it = concentrations.begin(); it != concentrations.end(); ++it) {
        sum_of_concentrations += it->second;
    }
    
    for (int i = 0; i < types_of_TFs_for_preprocessing.size(); ++i) {
        minimums_for_normalizations_of_potentials.insert(std::make_pair(types_of_TFs_for_preprocessing[i].type_name, 1000.0));
        find_specific_binding_sites(types_of_TFs_for_preprocessing[i], dnase, true);
        find_specific_binding_sites(types_of_TFs_for_preprocessing[i], dnase, false);
    }
        
    //find_potential_strength(true);
    //find_potential_strength(false);
    
    //generate_TF_appearance_time();
    
    
    



}

// Finds specific binding sites and weights of all sites
void Cell::find_specific_binding_sites(Transcription_factor tf, Parser_dnase_acc& dnase, bool forward) {
    std::vector<double> weights_of_binding;
    std::map<int, bool> specific_sites;
    if (forward) {
        for (int j = 0; j < forward_DNA.size() - tf.get_size(); ++j) {
            tf.change_coordinate_in_sequence(j);
            double weight_of_current_binding = tf.calculate_weight_of_binding(forward_DNA);
            weights_of_binding.push_back(weight_of_current_binding);
            if (weight_of_current_binding > 3.0 && dnase.is_in_interval(std::make_pair(j, j + tf.get_size()))) {
                specific_sites.insert(std::make_pair(j, false));
            }
        }
        weights_of_binding_of_all_tfs_forward.insert(std::make_pair(tf.type_name, weights_of_binding));
        specific_binding_sites_forward.insert(std::make_pair(tf.type_name, specific_sites));
    } else {
        for (int j = 0; j < reverse_DNA.size() - tf.get_size(); ++j) {
            tf.change_coordinate_in_sequence(j);
            double weight_of_current_binding = tf.calculate_weight_of_binding(reverse_DNA);
            weights_of_binding.push_back(weight_of_current_binding);
            if (weight_of_current_binding > 3.0 && dnase.is_in_interval(std::make_pair(reverse_DNA.size() - j - tf.get_size(),
                                                                                       reverse_DNA.size() - j))) {
                specific_sites.insert(std::make_pair(j, true));
            }
        }
        weights_of_binding_of_all_tfs_revcompl.insert(std::make_pair(tf.type_name, weights_of_binding));
        specific_binding_sites_revcompl.insert(std::make_pair(tf.type_name, specific_sites));
    }
}

double Cell::return_weight_of_binding(int i, std::vector<double> weights_of_binding) {
    double weight_of_other_DNA = -20.0;
    if (i >= 0 && i < weights_of_binding.size()) return weights_of_binding[i];
    else return weight_of_other_DNA;
}

void Cell::find_potential_strength(bool forward) {
    std::map<std::string, std::vector<double>> weights_of_binding_of_all_tfs;
    if (forward) weights_of_binding_of_all_tfs = weights_of_binding_of_all_tfs_forward;
    else weights_of_binding_of_all_tfs = weights_of_binding_of_all_tfs_revcompl;

    for (std::map<std::string, std::vector<double>>::iterator it = weights_of_binding_of_all_tfs.begin(); it != weights_of_binding_of_all_tfs.end(); ++it)
    {
        std::vector<std::pair<double, double>> weights_of_binding;
        std::string type_of_prot = it->first;
        
        std::cout << "\n" << type_of_prot << "\n";
        
        std::deque<double> deque_left;
        std::deque <double> deque_right;
        double value_to_the_left = 0.0;
        double value_to_the_right = 0.0;
        
        const int number_of_letters_to_define_potential = 6;
        
        int counter = 1;
        while (counter < number_of_letters_to_define_potential) {
            deque_left.push_back(return_weight_of_binding(-counter, it->second));
            value_to_the_left += return_weight_of_binding(-counter, it->second);
            deque_right.push_back(return_weight_of_binding(counter, it->second));
            value_to_the_right += return_weight_of_binding(counter, it->second);
            counter++;
        }
        weights_of_binding.push_back(std::make_pair(value_to_the_left, value_to_the_right));
        
        for (int i = 1; i < it->second.size(); ++i) {
            double value_to_decrease_left = deque_left.front();
            deque_left.pop_front();
            
            double value_to_decrease_right = deque_right.front();
            deque_right.pop_front();
            
            value_to_the_left += return_weight_of_binding(i - 1, it->second);
            deque_left.push_back(return_weight_of_binding(i - 1, it->second));
            value_to_the_left -= value_to_decrease_left;
            
            value_to_the_right += return_weight_of_binding(i + number_of_letters_to_define_potential, it->second);
            deque_right.push_back(return_weight_of_binding(i + number_of_letters_to_define_potential, it->second));

            value_to_the_right -= value_to_decrease_right;
            
            if (value_to_the_left / 5 < minimums_for_normalizations_of_potentials.at(it->first)) {
                minimums_for_normalizations_of_potentials[it->first] = value_to_the_left / 5;
            }
            
            weights_of_binding.push_back(std::make_pair(value_to_the_left / 5, value_to_the_right / 5));
        }
    }
}

void Cell::bind_tf_to_dna(int i, Transcription_factors_in_cell& tfs) {
    
}

void Cell::unbind_tf_from_dna(int i, Transcription_factors_in_cell& tfs) {
    
}

void Cell::generate_next_event(Transcription_factors_in_cell& tfs) {
    if (timeline.empty()) {
        timeline.push(std::make_pair(tfs.generate_next_characteristic_time_for_unbinded_TFs(), true));
    } else {
        if (timeline.top().second == false) {
            timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_unbinded_TFs(), false));
        } else {
            timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_unbinded_TFs(), true));
        }
    }
}

int Cell::move_tf_right_or_left(int i, double left, double right, std::string& type_name_of_prot) {
    double minimum = -1 * minimums_for_normalizations_of_potentials[type_name_of_prot] + 1;
    double a = -0.25;
    double b = 1.5;
    int answer = i;
    if (fRand(0.0, 2 * minimum + left + right) > minimum + left) {
        answer -= (int)(a * left + b);
    } else {
        answer += (int)(a * right + b);
    }
    return answer;
}

bool Cell::test_for_unbinding(double weight) {
    double a = -0.25;
    double b = -20.0;

    double dice = fRand(0.0, 1.0);
    double bound = 1.0 / (1.0 + exp(-a * (weight - b)));
    if (dice < bound) return true;
    else return false;
}

/*void Cell::generate_TF_appearance_time() {
    double current_time = 0.0;
    for (int i = 0; i < total_number_of_TFs_to_bind; ++i) {
        //std::cout << current_time << "\n";
        TF_appearance_time.push_back(current_time);
        current_time += generate_next_time(1 / sum_of_concentrations);
    }
}*/
