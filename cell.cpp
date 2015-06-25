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
        
    find_potential_strength(true);
    find_potential_strength(false);
    
    start_simulation(tfs, dnase);
}

// Finds specific binding sites and weights of all sites
void Cell::find_specific_binding_sites(Transcription_factor tf, Parser_dnase_acc& dnase, bool forward) {
    std::vector<double> weights_of_binding;
    std::map<int, bool> specific_sites;
    std::map<int, int> codes_of_specific_sites;
    if (forward) {
        for (int j = 0; j < forward_DNA.size() - tf.get_size(); ++j) {
            tf.change_coordinate_in_sequence(j);
            double weight_of_current_binding = tf.calculate_weight_of_binding(forward_DNA);
            weights_of_binding.push_back(weight_of_current_binding);
            if (weight_of_current_binding > 6.0 && dnase.is_in_interval(std::make_pair(j, j + tf.get_size() - 1))) {
                specific_sites.insert(std::make_pair(j, false));
                number_of_specific_binding_sites++;
                codes_of_specific_sites.insert(std::make_pair(j, number_of_specific_binding_sites));
            }
        }
        weights_of_binding_of_all_tfs_forward.insert(std::make_pair(tf.type_name, weights_of_binding));
        specific_binding_sites_forward.insert(std::make_pair(tf.type_name, specific_sites));
        codes_of_specific_binding_sites_forward.insert(std::make_pair(tf.type_name, codes_of_specific_sites));
    } else {
        for (int j = 0; j < reverse_DNA.size() - tf.get_size(); ++j) {
            tf.change_coordinate_in_sequence(j);
            double weight_of_current_binding = tf.calculate_weight_of_binding(reverse_DNA);
            weights_of_binding.push_back(weight_of_current_binding);
            if (weight_of_current_binding > 6.0 && dnase.is_in_interval(std::make_pair(reverse_DNA.size() - j - tf.get_size(),
                                                                                       reverse_DNA.size() - j - 1))) {
                specific_sites.insert(std::make_pair(j, true));
                number_of_specific_binding_sites++;
                codes_of_specific_sites.insert(std::make_pair(j, number_of_specific_binding_sites));
            }
        }
        weights_of_binding_of_all_tfs_revcompl.insert(std::make_pair(tf.type_name, weights_of_binding));
        specific_binding_sites_revcompl.insert(std::make_pair(tf.type_name, specific_sites));
        codes_of_specific_binding_sites_forward.insert(std::make_pair(tf.type_name, codes_of_specific_sites));
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
        if (forward == true) {
            potential_strength_forward[type_of_prot] = weights_of_binding;
        } else {
            potential_strength_revcompl[type_of_prot] = weights_of_binding;
        }
    }
}

bool Cell::bind_tf_to_dna(int i, Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
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
                return true;
            } else {
                return false;
            }
        } else if (dnase.is_in_interval(std::make_pair(reverse_DNA.size() - coordinate - tfs.get_size_of_TF(i),
                                                       reverse_DNA.size() - coordinate - 1))) {
            if (binding_site_is_free(false, coordinate, tfs.get_size_of_TF(i))) {

            tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
            tfs.bind_tf_to_dna(i, false);
            no_access_dna_revcompl.insert(std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1));
                return true;
            } else {
                return false;
            }
        }
    }
    return false;
}

void Cell::unbind_tf_from_dna(int i, Transcription_factors_in_cell& tfs) {
    std::cout << "UNBINDING!\n";
    int coordinate_in_sequence = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        no_access_dna_forward.erase(std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_tf_by_index(i).get_size() - 1));
    } else {
        no_access_dna_revcompl.erase(std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_tf_by_index(i).get_size() - 1));
    }
    tfs.unbind_tf_from_dna(i);
}

void Cell::generate_next_event(Transcription_factors_in_cell& tfs) {
    if (timeline.empty()) {
        double first_event_time = tfs.generate_next_characteristic_time_for_unbinded_TFs();
        timeline.push(std::make_pair(first_event_time, 1));
    } else {
        if (timeline.top().second == 2 && tfs.number_of_binded_dna() != 0) {
            timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_binded_TFs(), 2));
        } else {
            timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_unbinded_TFs(), 1));
        }
    }
}

void Cell::generate_next_event_with_binded(Transcription_factors_in_cell& tfs) {
    timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_binded_TFs(), 2));
}

void Cell::generate_next_event_with_unbinded(Transcription_factors_in_cell& tfs) {
    timeline.push(std::make_pair(timeline.top().first + tfs.generate_next_characteristic_time_for_unbinded_TFs(), 1));
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
    std::cout << "new coord " << answer << "\n";
    return answer;
}

bool Cell::test_for_unbinding(double weight) {
    std::cout << "Testing for unbinding...\n";
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


bool Cell::test_for_specific_binding(int i, Transcription_factors_in_cell& tfs) {
    int coordinate_in_sequence = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
    std::string tf_type = tfs.get_tf_by_index(i).type_name;
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        if (specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name].find(coordinate_in_sequence) != specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name].end()) {
            specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name][tfs.get_tf_by_index(i).get_coordinate_in_sequence()] = true;
            tfs.get_tf_by_index(i).set_binded_to_dna_specifically();
            no_access_dna_forward.insert(std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_size_of_TF(i) - 1));
            int code_of_current_specific_site = codes_of_specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name][coordinate_in_sequence];
            specific_sites_binded_or_not[code_of_current_specific_site] = true;
            final_frequency_of_combinations[specific_sites_binded_or_not] = true;
            counter_of_steps_without_changes = 0;
            std::cout << "Binded specifically!\n";
            return true;
        }
    } else {
        if (specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name].find(coordinate_in_sequence) != specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name].end()) {
            specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name][tfs.get_tf_by_index(i).get_coordinate_in_sequence()] = true;
            tfs.get_tf_by_index(i).set_binded_to_dna_specifically();
            no_access_dna_revcompl.insert(std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_size_of_TF(i) - 1));
            int code_of_current_specific_site = codes_of_specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name][coordinate_in_sequence];
            specific_sites_binded_or_not[code_of_current_specific_site] = true;
            final_frequency_of_combinations[specific_sites_binded_or_not] = true;
            counter_of_steps_without_changes = 0;
            std::cout << "Binded specifically!\n";
            return true;
        }
    }
    std::cout << "Not a specific binding...\n";
    counter_of_steps_without_changes++;
    return false;
}

void Cell::one_dimensional_slinding_of_TF(int i, Transcription_factors_in_cell& tfs) {
    int new_coord = 0;
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        std::cout << tfs.get_tf_by_index(i).get_coordinate_in_sequence() << "\n";
        std::cout << potential_strength_forward[tfs.get_tf_by_index(i).type_name][i].first << "\n";
        std::cout << potential_strength_forward[tfs.get_tf_by_index(i).type_name][i].second << "\n";
        std::cout << tfs.get_tf_by_index(i).type_name << "\n";
        int current_coordinate = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
        new_coord = move_tf_right_or_left(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), potential_strength_forward[tfs.get_tf_by_index(i).type_name][i].first, potential_strength_forward[tfs.get_tf_by_index(i).type_name][i].second, tfs.get_tf_by_index(i).type_name);
        no_access_dna_forward.erase(no_access_dna_forward.find(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) ));
        if (binding_site_is_free(true, new_coord, tfs.get_tf_by_index(i).get_size())) {
            tfs.get_tf_by_index(i).change_coordinate_in_sequence(new_coord);
            no_access_dna_forward.insert(std::make_pair(new_coord, new_coord + tfs.get_size_of_TF(i) - 1) );
        } else {
            no_access_dna_forward.insert(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) );
            std::cout << "CAN NOT MOVE!\n";
        }
    } else {
        std::cout << tfs.get_tf_by_index(i).get_coordinate_in_sequence() << "\n";
        std::cout << potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][i].first << "\n";
        std::cout << potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][i].second << "\n";
        std::cout << tfs.get_tf_by_index(i).type_name << "\n";
        int current_coordinate = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
        new_coord = move_tf_right_or_left(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][i].first, potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][i].second, tfs.get_tf_by_index(i).type_name);
        no_access_dna_revcompl.erase(no_access_dna_revcompl.find(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) ));

        if (binding_site_is_free(false, new_coord, tfs.get_tf_by_index(i).get_size())) {
            tfs.get_tf_by_index(i).change_coordinate_in_sequence(new_coord);
            no_access_dna_revcompl.insert(std::make_pair(new_coord, new_coord + tfs.get_size_of_TF(i) - 1) );
        }  else {
            no_access_dna_revcompl.insert(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) );
            std::cout << "CAN NOT MOVE!\n";
        }
    }
}


std::map<std::vector<bool>, bool> Cell::get_final_frequency_of_combinations() {
    return final_frequency_of_combinations;
}



void Cell::start_simulation(Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
    int number_of_steps_without_change = 0;
    generate_next_event(tfs);
    specific_sites_binded_or_not.resize(number_of_specific_binding_sites);
    
    
    while (counter_of_steps_without_changes < number_of_steps_before_stabilization) {
        std::pair<double, int> next_event = timeline.top();
        std::cout << next_event.first << " " << next_event.second << "\n";
        if (next_event.second == 1) {
            std::cout << number_of_steps_without_change << " operation with unbinded" << "\n";
            int next_tf = tfs.choose_next_unbinded_DNA_to_interact();
            bool result_of_binding = bind_tf_to_dna(next_tf, tfs, dnase);
            if (result_of_binding == true) {
                generate_next_event_with_binded(tfs);
            }
            else generate_next_event_with_unbinded(tfs);
            test_for_specific_binding(next_tf, tfs);
        } else if (next_event.second == 2){
            std::cout << number_of_steps_without_change << " operation with binded" << "\n";
            int next_tf = tfs.choose_next_binded_DNA_to_interact();
            if (next_tf >= 0) {
                std::cout << "Next binded TF: " << next_tf << "\n";
                
                for (int i = 0; i < number_of_one_dim_slidings_in_step; ++i) {
                    std::cout << "Here should be 1D sliding\n";
                    one_dimensional_slinding_of_TF(next_tf, tfs);
                    test_for_specific_binding(next_tf, tfs);
                }
            }
            generate_next_event_with_binded(tfs);
        }
        timeline.pop();
        number_of_steps_without_change++;
    }
}


