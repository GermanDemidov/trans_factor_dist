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
    total_number_of_TFs_to_bind = transcription_factor_instances.size();
    
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
    
    start_simulation(tfs, dnase);
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
            if (weight_of_current_binding > 6.0 && dnase.is_in_interval(std::make_pair(j, j + tf.get_size() - 1))) {
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
            if (weight_of_current_binding > 6.0 && dnase.is_in_interval(std::make_pair(reverse_DNA.size() - j - tf.get_size(),
                                                                                       reverse_DNA.size() - j - 1))) {
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

void Cell::bind_tf_to_dna(int i, Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
    int choose_strand_to_bind = pick_a_number(0,1);
    int coordinate = 0;
    if (choose_strand_to_bind == 0) {
        coordinate = pick_a_number(0, forward_DNA.size() - 1);
    } else coordinate = pick_a_number(0, reverse_DNA.size() - 1);
    if (choose_strand_to_bind == 0) {
        if (dnase.is_in_interval(std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1))) {
            if (binding_site_is_free(true, coordinate, tfs.get_size_of_TF(i))) {
                tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
                tfs.bind_tf_to_dna(i, true);
                no_access_dna_forward.insert(std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1));
            } else {std::cout << "\n\n\n\n\n\n\nBinding site was occupied\n\n\n\n\n\n\n";}
        } else if (dnase.is_in_interval(std::make_pair(reverse_DNA.size() - coordinate - tfs.get_size_of_TF(i),
                                                       reverse_DNA.size() - coordinate - 1))) {
            if (binding_site_is_free(false, coordinate, tfs.get_size_of_TF(i))) {

            tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
            tfs.bind_tf_to_dna(i, false);
            no_access_dna_revcompl.insert(std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1));

            } else {std::cout << "\n\n\n\n\n\n\nBinding site was occupied\n\n\n\n\n\n\n";}
        }
    }
}

void Cell::unbind_tf_from_dna(int i, Transcription_factors_in_cell& tfs) {
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        no_access_dna_forward.erase(std::make_pair(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), tfs.get_tf_by_index(i).get_coordinate_in_sequence() + tfs.get_tf_by_index(i).get_size() - 1));
    } else {
        no_access_dna_revcompl.erase(std::make_pair(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), tfs.get_tf_by_index(i).get_coordinate_in_sequence() + tfs.get_tf_by_index(i).get_size() - 1));
    }
    tfs.unbind_tf_from_dna(i);
}

void Cell::generate_next_event(Transcription_factors_in_cell& tfs) {
    if (timeline.empty()) {
        double first_event_time = tfs.generate_next_characteristic_time_for_unbinded_TFs();
        timeline.push(std::make_pair(first_event_time, true));
    } else {
        if (timeline.top().second == false && tfs.number_of_binded_dna() != 0) {
            timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_binded_TFs(), false));
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


bool Cell::binding_site_is_free(bool forward, int coord, int len_of_tf) {
    
    std::pair<int, int> pair_to_test = std::make_pair(coord, coord + len_of_tf - 1);
    int common_part = 0;
    std::set<std::pair<int, int>> no_access_dna;
    if (forward) {
        no_access_dna = no_access_dna_forward;
    } else {
        no_access_dna = no_access_dna_revcompl;
    }
    for (auto& coords : no_access_dna) {
        if (coords.first > pair_to_test.second) {
                break;
        } else {
                int first_diff = coords.first - pair_to_test.first;
                int second_diff = coords.second - pair_to_test.second;
                if (first_diff >= 0 && second_diff >= 0)
                    common_part = pair_to_test.second - coords.first;
                else if (first_diff >= 0 && second_diff <= 0)
                    return false;
                else if (first_diff <= 0 && second_diff >= 0)
                    return false;
                else if (first_diff <= 0 && second_diff <= 0)
                    common_part = coords.second - pair_to_test.first;
        }
        if (common_part > 0)
                    return false;
    }
    return true;
}

void Cell::start_simulation(Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
    int number_of_steps_without_change = 0;
    generate_next_event(tfs);
    
    while (number_of_steps_without_change < number_of_steps_before_stabilization) {
        std::pair<double, bool> next_event = timeline.top();
        std::cout << next_event.first << " " << next_event.second << "\n";
        if (next_event.second == true) {
            std::cout << number_of_steps_without_change << " next_event.second == true" << "\n";
            int next_tf = tfs.choose_next_unbinded_DNA_to_interact();
            bind_tf_to_dna(next_tf, tfs, dnase);
            if (tfs.number_of_binded_dna() == 0)
                timeline.push(std::make_pair(timeline.top().first + generate_next_time(0.01), false));
            test_for_specific_binding(next_tf, tfs);
        } else {
            std::cout << number_of_steps_without_change << " next_event.second == false" << "\n";
            int next_tf = tfs.choose_next_binded_DNA_to_interact();
            if (next_tf >= 0) {
                std::cout << "NEXT BINDED TF: " << next_tf << "\n";
                
                for (int i = 0; i < number_of_one_dim_slidings_in_step; ++i) {
                    one_dimensional_slinding_of_TF(next_tf, tfs);
                    test_for_specific_binding(next_tf, tfs);
                }
            }
        }
        generate_next_event(tfs);
        timeline.pop();
        number_of_steps_without_change++;
    }
}

void Cell::test_for_specific_binding(int i, Transcription_factors_in_cell& tfs) {
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        if (specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name].find(tfs.get_tf_by_index(i).get_coordinate_in_sequence()) != specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name].end()) {
            specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name][tfs.get_tf_by_index(i).get_coordinate_in_sequence()] = true;
            tfs.get_tf_by_index(i).set_binded_to_dna_specifically();
        }
    } else {
        if (specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name].find(tfs.get_tf_by_index(i).get_coordinate_in_sequence()) != specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name].end()) {
            specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name][tfs.get_tf_by_index(i).get_coordinate_in_sequence()] = true;
            tfs.get_tf_by_index(i).set_binded_to_dna_specifically();
        }
        
    }
}

void Cell::one_dimensional_slinding_of_TF(int i, Transcription_factors_in_cell& tfs) {
    std::cout << "\n\n\nWOW STARTED SLIDING!\n\n\n";
    int new_coord = 0;
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        new_coord = move_tf_right_or_left(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), potential_strength_forward[tfs.get_tf_by_index(i).type_name][i].first, potential_strength_forward[tfs.get_tf_by_index(i).type_name][i].second, tfs.get_tf_by_index(i).type_name);
        if (binding_site_is_free(true, new_coord, tfs.get_tf_by_index(i).get_size())) {
            tfs.get_tf_by_index(i).change_coordinate_in_sequence(new_coord);
        }
    } else {
        new_coord = move_tf_right_or_left(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][i].first, potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][i].second, tfs.get_tf_by_index(i).type_name);
        if (binding_site_is_free(false, new_coord, tfs.get_tf_by_index(i).get_size())) {
            tfs.get_tf_by_index(i).change_coordinate_in_sequence(new_coord);
        }
    }
}

/*void Cell::generate_TF_appearance_time() {
    double current_time = 0.0;
    for (int i = 0; i < total_number_of_TFs_to_bind; ++i) {
        //std::cout << current_time << "\n";
        TF_appearance_time.push_back(current_time);
        current_time += generate_next_time(1 / sum_of_concentrations);
    }
}*/
