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
    std::vector<std::pair<int, int> > dnase_acc = dnase.return_no_acc_intervals();
    for (int i = 0; i < dnase_acc.size(); i++) {
        no_access_dna_forward.insert(dnase_acc[i]);
        no_access_dna_revcompl.insert(std::make_pair((int)reverse_DNA.size() - dnase_acc[i].second - 1, (int)reverse_DNA.size() - dnase_acc[i].first));
    }
    
    cell_id = id_of_cell;
    concentrations = conc;
    total_number_of_TFs_to_bind = transcription_factor_instances.size();
    
    for (std::map<std::string, double>::iterator it = concentrations.begin(); it != concentrations.end(); ++it) {
        sum_of_concentrations += it->second;
    }
    
    for (int i = 0; i < types_of_TFs_for_preprocessing.size(); ++i) {
        sizes_of_tfs[types_of_TFs_for_preprocessing[i].type_name] =types_of_TFs_for_preprocessing[i].get_size();
        minimums_for_normalizations_of_potentials.insert(std::make_pair(types_of_TFs_for_preprocessing[i].type_name, 1000.0));
        find_specific_binding_sites(types_of_TFs_for_preprocessing[i], dnase, true);
        find_specific_binding_sites(types_of_TFs_for_preprocessing[i], dnase, false);
    }
    std::cout << "\nTotal number of specific binding sites on the current DNA " << number_of_specific_binding_sites << "\n";
    
    find_potential_strength(true);
    find_potential_strength(false);
    times_of_changes.push_back(0.0);
    
    start_simulation(tfs, dnase);
}

void Cell::null_everything(Parser_dnase_acc& dnase) {
    for (int i = 0; i != specific_sites_binded_or_not.size(); i++)
        specific_sites_binded_or_not[i] = false;
    time_for_binded_state.clear();
    times_of_changes.clear();
    times_of_changes.push_back(0.0);
    time_for_binded_state.push_back(0.0);
    time_for_unbinded_state.clear();
    time_for_unbinded_state.push_back(0.0);
    final_frequency_of_combinations.clear();
    final_frequency_of_combinations[specific_sites_binded_or_not] = true;
    no_access_dna_forward.clear();
    no_access_dna_revcompl.clear();
    std::vector<std::pair<int, int> > dnase_acc = dnase.return_no_acc_intervals();
    for (int i = 0; i < dnase_acc.size(); i++) {
        no_access_dna_forward.insert(dnase_acc[i]);
        no_access_dna_revcompl.insert(std::make_pair((int)reverse_DNA.size() - dnase_acc[i].second - 1, (int)reverse_DNA.size() - dnase_acc[i].first));
    }
    
    // tfs - ?
    counter_of_steps_without_changes = 0;
    while(!timeline.empty())
        timeline.pop();
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
            //std::cout << weight_of_current_binding << ",";
            weights_of_binding.push_back(weight_of_current_binding);
            if (weight_of_current_binding > 3.0 && dnase.is_in_interval(std::make_pair(j, j + tf.get_size() - 1))) {
                specific_sites.insert(std::make_pair(j, false));
                std::cout << "Forward\t" << tf.type_name << "\t" << j << "\t" << j + tf.get_size() - 1 << "\t"
                << reverse_DNA.size() - j - tf.get_size() + 1 << "\t" << reverse_DNA.size() - j << "\t"
                << weight_of_current_binding << "\t" << number_of_specific_binding_sites << "\n";
                codes_of_specific_sites.insert(std::make_pair(j, number_of_specific_binding_sites));
                number_of_specific_binding_sites++;
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
            if (weight_of_current_binding > 3.0 && dnase.is_in_interval(std::make_pair(reverse_DNA.size() - j - tf.get_size(),
                                                                                       reverse_DNA.size() - j - 1))) {
                specific_sites.insert(std::make_pair(j, true));
                std::cout << "Revcompl\t" << tf.type_name << "\t" << j << "\t" << j + tf.get_size() - 1 << "\t"
                <<  forward_DNA.size() - j - tf.get_size() + 1 << "\t" << forward_DNA.size() - j << "\t"
                <<  weight_of_current_binding << "\t" << number_of_specific_binding_sites << "\n";
                codes_of_specific_sites.insert(std::make_pair(j, number_of_specific_binding_sites));
                number_of_specific_binding_sites++;
            }
        }
        weights_of_binding_of_all_tfs_revcompl.insert(std::make_pair(tf.type_name, weights_of_binding));
        specific_binding_sites_revcompl.insert(std::make_pair(tf.type_name, specific_sites));
        codes_of_specific_binding_sites_revcompl.insert(std::make_pair(tf.type_name, codes_of_specific_sites));
    }
}

double Cell::return_weight_of_binding(int i, std::vector<double>& weights_of_binding) {
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
        weights_of_binding.push_back(std::make_pair(value_to_the_left / (number_of_letters_to_define_potential - 1) , value_to_the_right / (number_of_letters_to_define_potential - 1)));
        
        for (int i = 0; i < it->second.size() + sizes_of_tfs[type_of_prot] - 1; ++i) {
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
            
            if (value_to_the_left / (number_of_letters_to_define_potential - 1) < minimums_for_normalizations_of_potentials.at(it->first)) {
                minimums_for_normalizations_of_potentials[it->first] = value_to_the_left / (number_of_letters_to_define_potential - 1);
            }
            
            weights_of_binding.push_back(std::make_pair(value_to_the_left / (number_of_letters_to_define_potential - 1), value_to_the_right / (number_of_letters_to_define_potential - 1)));
        }
        /*for (int i = 0; i < weights_of_binding.size(); i++) {
         std::cout << weights_of_binding[i].first << weights_of_binding[i].second << "\n";
         }
         std::cout << "LENGTH " << weights_of_binding.size() << " " << forward << "\n";*/
        if (forward == true) {
            potential_strength_forward[type_of_prot] = weights_of_binding;
        } else {
            potential_strength_revcompl[type_of_prot] = weights_of_binding;
        }
    }
}

int Cell::choose_position_with_maximum_weight_inside_neighborhood(int coord, int size_of_tf, std::vector<double>& weights) {
    int new_coord = coord;
    double max_weight = -1000.0;
    for (int j = std::max(coord - size_of_tf, 0); j < std::min(coord + size_of_tf, (int)reverse_DNA.size()); j++) {
        if (weights[j] > max_weight) {
            max_weight = weights[j];
            new_coord = j;
        }
    }
    return new_coord;
}

/*bool Cell::bind_tf_to_dna(int i, Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
    int choose_strand_to_bind = pick_a_number(0,1);
    int coordinate = 0;
    if (choose_strand_to_bind == 0) {
        coordinate = pick_a_number(0, (int)forward_DNA.size() - 1);
    } else {
        coordinate = pick_a_number(0, (int)reverse_DNA.size() - 1);
    }
    std::pair<int, int> coordinates = std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1);
    if (choose_strand_to_bind == 0) {
        if (binding_site_is_free(true, coordinate, tfs.get_size_of_TF(i))) {
            tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
            tfs.bind_tf_to_dna(i, true, coordinate);
            no_access_dna_forward.insert(coordinates);
            no_access_dna_revcompl.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
            return true;
        } else {
            return false;
        }
    }
    else if (binding_site_is_free(false, coordinate, tfs.get_size_of_TF(i))) {
        //std::cout << "free...\n";
        
        tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
        tfs.bind_tf_to_dna(i, false, coordinate);
        no_access_dna_revcompl.insert(coordinates);
        no_access_dna_forward.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
        return true;
    } else {
        //std::cout << "not free..." << coordinate << "\n";
        return false;
    }
    return false;
}*/

bool Cell::bind_tf_to_dna(int i, Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
    int choose_strand_to_bind = pick_a_number(0,1);
    int coordinate = 0;
    if (choose_strand_to_bind == 0) {
        coordinate = pick_a_number(0, (int)forward_DNA.size() - 1);
    } else {
        coordinate = pick_a_number(0, (int)reverse_DNA.size() - 1);
    }
    std::vector<double> weights_of_binding;
    int size_of_tf = tfs.get_tf_by_index(i).get_size();
    int new_coordinate = 0;
    // maximum number of stochastic jumps before protein lands
    int counter_of_stochastic_jumps = 0;
    
    while (counter_of_stochastic_jumps < max_number_of_stochastic_jumps_before_landing) {
        if (choose_strand_to_bind == 0) {
            weights_of_binding = weights_of_binding_of_all_tfs_forward[(tfs.get_tf_by_index(i).type_name)];
            coordinate = choose_position_with_maximum_weight_inside_neighborhood(coordinate, size_of_tf, weights_of_binding);
        } else {
            weights_of_binding = weights_of_binding_of_all_tfs_revcompl[(tfs.get_tf_by_index(i).type_name)];
            coordinate = choose_position_with_maximum_weight_inside_neighborhood(coordinate, size_of_tf, weights_of_binding);
        }
        
        if (choose_strand_to_bind == 0) {
            if (binding_site_is_free(true, coordinate, tfs.get_size_of_TF(i))) {
                tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
                tfs.bind_tf_to_dna(i, true, coordinate);
                
                std::pair<int, int> coordinates = std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1);

                no_access_dna_forward.insert(coordinates);
                no_access_dna_revcompl.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
                return true;
            } else {
                counter_of_stochastic_jumps++;
                choose_strand_to_bind = pick_a_number(0,1);
                new_coordinate = return_normal_number(coordinate, sd_of_jumps);
                //std::cout << "Jump " << coordinate << " " << new_coordinate << " " << counter_of_stochastic_jumps << "\n";
                if (choose_strand_to_bind == 0) {
                    if (new_coordinate < 0 || new_coordinate > (int)forward_DNA.size() - 1) {
                        return false;
                    }
                } else {
                    if (new_coordinate < 0 || new_coordinate > (int)reverse_DNA.size() - 1) {
                        return false;
                    }
                }
                if (choose_strand_to_bind != 0) {
                    //std::cout << "Changed strand!\n";
                    coordinate = (int)reverse_DNA.size() - new_coordinate;
                } else {
                    coordinate = new_coordinate;
                }
            }
        }
        else {
            if (binding_site_is_free(false, coordinate, tfs.get_size_of_TF(i))) {
                tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
                tfs.bind_tf_to_dna(i, false, coordinate);
                
                std::pair<int, int> coordinates = std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1);

                no_access_dna_revcompl.insert(coordinates);
                no_access_dna_forward.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
                return true;
            } else {
                counter_of_stochastic_jumps++;
                choose_strand_to_bind = pick_a_number(0,1);
                new_coordinate = return_normal_number(coordinate, sd_of_jumps);
                //std::cout << "Jump " << coordinate << " " << new_coordinate << " " << counter_of_stochastic_jumps << "\n";
                if (choose_strand_to_bind == 0) {
                    if (new_coordinate < 0 || new_coordinate > (int)forward_DNA.size() - 1) {
                        return false;
                    }
                } else {
                    if (new_coordinate < 0 || new_coordinate > (int)reverse_DNA.size() - 1) {
                        return false;
                    }
                }
                if (choose_strand_to_bind != 1) {
                    //std::cout << "Changed strand!\n";
                    coordinate = (int)reverse_DNA.size() - new_coordinate;
                } else {
                    coordinate = new_coordinate;
                }
            }
        }
    }
    return false;
}

void Cell::unbind_tf_from_dna(int i, Transcription_factors_in_cell& tfs) {
    int coordinate_in_sequence = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
    std::pair<int, int> coordinates = std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_size_of_TF(i) - 1);

    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        no_access_dna_forward.erase(coordinates);
        no_access_dna_revcompl.erase(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
    } else {
        no_access_dna_revcompl.erase(coordinates);
        no_access_dna_forward.erase(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
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
    double a = -0.75;
    double b = 1;
    int answer = i;
    double first_square = exp(left / 10); //minimum + left;
    double second_square = exp(right / 10); //minimum + right;
    double bound_to_left = first_square / (first_square + second_square);
    double random_number = ((double) rand() / (RAND_MAX));
    
    if (random_number < bound_to_left) {
        int dist_to_move = (int)(a * left + b);
        if (dist_to_move < 1) {
            dist_to_move = 1;
        }
        answer -= dist_to_move;
    } else {
        int dist_to_move = (int)(a * right + b);
        if (dist_to_move < 1) {
            dist_to_move = 1;
        }
        answer += dist_to_move;
    }
    //std::cout << "new coord " << answer << "\n";
    return answer;
}

bool Cell::test_for_unbinding(double weight) {
    double a = -0.3;
    double b = -15.0;
    
    double dice = fRand(0.0, 1.0);
    double bound = 1.0 / (1.0 + exp(-a * (weight - b)));
    if (dice < bound) {
        return true;
    }
    else return false;
}


bool Cell::binding_site_is_free(bool forward, int coord, int len_of_tf) {
    
    std::pair<int, int> pair_to_test = std::make_pair(coord, coord + len_of_tf - 1);
    int common_part = 0;
    std::set<std::pair<int, int>>* no_access_dna;
    if (forward) {
        no_access_dna = &no_access_dna_forward;
    } else {
        no_access_dna = &no_access_dna_revcompl;
    }
    for (auto& coords : *no_access_dna) {
        common_part = std::min(coords.second, pair_to_test.second) - std::max(coords.first, pair_to_test.first);
        if (common_part >= 0)
            return (false);
    }
    
    return true;
}


bool Cell::test_for_specific_binding(int i, Transcription_factors_in_cell& tfs) {
    int coordinate_in_sequence = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
    int current_coordinate = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
    std::string tf_type = tfs.get_tf_by_index(i).type_name;
    bool flag_of_previous_binding = false;
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        std::pair<int, int> coordinates = std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_size_of_TF(i) - 1) ;
        if (no_access_dna_forward.find(coordinates ) != no_access_dna_forward.end()) {
            flag_of_previous_binding = true;
        }
        if (flag_of_previous_binding) {
            no_access_dna_forward.erase(coordinates );
            no_access_dna_revcompl.erase(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
        }
        if (specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name].find(coordinate_in_sequence) != specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name].end()) {
            specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name][tfs.get_tf_by_index(i).get_coordinate_in_sequence()] = true;
            tfs.get_tf_by_index(i).set_binded_to_dna_specifically();
            if (tf_type != repressor) {
                no_access_dna_forward.insert(broad_borders(1, coordinates));
                no_access_dna_revcompl.insert(broad_borders(1, reverse_interval_for_another_strand(reverse_DNA.size(), coordinates)));
            } else {
                if (
                    (binding_site_is_free(true, coordinate_in_sequence - length_of_repression, tfs.get_size_of_TF(i) + 2 * length_of_repression)) && 
                    (binding_site_is_free(false, reverse_DNA.size() - coordinate_in_sequence - tfs.get_size_of_TF(i) - length_of_repression,
                                          tfs.get_size_of_TF(i) + 2 * length_of_repression))
                    ) {
                        no_access_dna_forward.insert(std::make_pair(coordinate_in_sequence - length_of_repression, coordinate_in_sequence + tfs.get_size_of_TF(i) + length_of_repression));
                        no_access_dna_revcompl.insert( broad_borders(length_of_repression, reverse_interval_for_another_strand(reverse_DNA.size(), coordinates) ));
                    //std::cout << "REPRESSION HAPPENED FORWARD " << coordinate_in_sequence << "\n";
                } else {
                    no_access_dna_forward.insert(std::make_pair(coordinate_in_sequence - 1, coordinate_in_sequence + tfs.get_size_of_TF(i)));
                    no_access_dna_revcompl.insert(broad_borders(1, reverse_interval_for_another_strand(reverse_DNA.size(), coordinates)));
                    //std::cout << "REPRESSION NOT HAPPENED " << coordinate_in_sequence << "\n";
                }
            }
            int code_of_current_specific_site = codes_of_specific_binding_sites_forward[tfs.get_tf_by_index(i).type_name][coordinate_in_sequence];
            specific_sites_binded_or_not[code_of_current_specific_site] = true;
            final_frequency_of_combinations[specific_sites_binded_or_not] = true;
            counter_of_steps_without_changes = 0;
            //std::cout << tf_type << " " << coordinate_in_sequence << " Binded specifically!\n";
            return true;
        }
        if (flag_of_previous_binding) {
            no_access_dna_forward.insert(std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_size_of_TF(i) - 1) );
            no_access_dna_revcompl.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
        }
    } else {
        std::pair<int, int> coordinates = std::make_pair(coordinate_in_sequence, coordinate_in_sequence + tfs.get_size_of_TF(i) - 1) ;

        if (no_access_dna_revcompl.find(coordinates ) != no_access_dna_revcompl.end()) {
            flag_of_previous_binding = true;
        }
        if (flag_of_previous_binding) {
            no_access_dna_revcompl.erase(coordinates);
            no_access_dna_forward.erase(reverse_interval_for_another_strand(forward_DNA.size(), coordinates));
        }
        if (specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name].find(coordinate_in_sequence) != specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name].end()) {
            specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name][tfs.get_tf_by_index(i).get_coordinate_in_sequence()] = true;
            tfs.get_tf_by_index(i).set_binded_to_dna_specifically();
            if (tf_type != repressor) {
                no_access_dna_revcompl.insert(broad_borders(1, coordinates));
                no_access_dna_forward.insert(broad_borders(1, reverse_interval_for_another_strand(forward_DNA.size(), coordinates)));
            } else {
                if (
                    (binding_site_is_free(false, coordinate_in_sequence - length_of_repression, tfs.get_size_of_TF(i) + 2 * length_of_repression)) &&
                    (binding_site_is_free(true, forward_DNA.size() - coordinate_in_sequence - tfs.get_size_of_TF(i) - length_of_repression,
                                          tfs.get_size_of_TF(i) + 2 * length_of_repression))
                    ) {
                    no_access_dna_revcompl.insert(std::make_pair(coordinate_in_sequence - length_of_repression, coordinate_in_sequence + tfs.get_size_of_TF(i) + length_of_repression));
                    no_access_dna_forward.insert( broad_borders(length_of_repression, reverse_interval_for_another_strand(forward_DNA.size(), coordinates) ));
                    //std::cout << "REPRESSION HAPPENED REVCOMPL " << coordinate_in_sequence << "\n";
                } else {
                    no_access_dna_revcompl.insert(std::make_pair(coordinate_in_sequence - 1, coordinate_in_sequence + tfs.get_size_of_TF(i)));
                    no_access_dna_forward.insert(broad_borders(1, reverse_interval_for_another_strand(forward_DNA.size(), coordinates)));
                }
            }
            int code_of_current_specific_site = codes_of_specific_binding_sites_revcompl[tfs.get_tf_by_index(i).type_name][coordinate_in_sequence];
            specific_sites_binded_or_not[code_of_current_specific_site] = true;
            final_frequency_of_combinations[specific_sites_binded_or_not] = true;
            counter_of_steps_without_changes = 0;
            //std::cout << tf_type << " " << coordinate_in_sequence << " Binded specifically!\n";
            return true;
        }
        if (flag_of_previous_binding) {
            no_access_dna_revcompl.insert(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) );
            no_access_dna_forward.insert(reverse_interval_for_another_strand(forward_DNA.size(), coordinates));

        }
    }
    counter_of_steps_without_changes++;
    return false;
}

void Cell::one_dimensional_slinding_of_TF(int i, Transcription_factors_in_cell& tfs) {
    int new_coord = 0;
    if (tfs.get_tf_by_index(i).is_binded_to_forward()) {
        //std::cout << tfs.get_tf_by_index(i).type_name << "\n";
        int current_coordinate = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
        //std::cout << "Forward! \n";
        new_coord = move_tf_right_or_left(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), potential_strength_forward[tfs.get_tf_by_index(i).type_name][current_coordinate].first, potential_strength_forward[tfs.get_tf_by_index(i).type_name][current_coordinate].second, tfs.get_tf_by_index(i).type_name);
        if (new_coord > 0 && new_coord < forward_DNA.length()) {
            no_access_dna_forward.erase(no_access_dna_forward.find(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) ));
            no_access_dna_revcompl.erase(no_access_dna_revcompl.find(reverse_interval_for_another_strand(reverse_DNA.size(), std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) )));

            if (binding_site_is_free(true, new_coord, tfs.get_tf_by_index(i).get_size())) {
                tfs.get_tf_by_index(i).change_coordinate_in_sequence(new_coord);
                no_access_dna_forward.insert(std::make_pair(new_coord, new_coord + tfs.get_size_of_TF(i) - 1) );
                no_access_dna_revcompl.insert(reverse_interval_for_another_strand(reverse_DNA.size(), std::make_pair(new_coord, new_coord + tfs.get_size_of_TF(i) - 1) ));

            } else {
                no_access_dna_forward.insert(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) );
                no_access_dna_revcompl.insert(reverse_interval_for_another_strand(reverse_DNA.size(), std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1)) );

                //std::cout << "Can not move!\n";
            }
        }
    } else {
        int current_coordinate = tfs.get_tf_by_index(i).get_coordinate_in_sequence();
        //std::cout << tfs.get_tf_by_index(i).type_name << "\n";
        //std::cout << "Revcompl! \n";
        
        new_coord = move_tf_right_or_left(tfs.get_tf_by_index(i).get_coordinate_in_sequence(), potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][current_coordinate].first, potential_strength_revcompl[tfs.get_tf_by_index(i).type_name][current_coordinate].second, tfs.get_tf_by_index(i).type_name);
        if (new_coord > 0 && new_coord < reverse_DNA.length()) {
            
            no_access_dna_revcompl.erase(no_access_dna_revcompl.find(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) ));
            no_access_dna_forward.erase(no_access_dna_forward.find(reverse_interval_for_another_strand(forward_DNA.size(), std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) )));

            if (binding_site_is_free(false, new_coord, tfs.get_tf_by_index(i).get_size())) {
                tfs.get_tf_by_index(i).change_coordinate_in_sequence(new_coord);
                no_access_dna_revcompl.insert(std::make_pair(new_coord, new_coord + tfs.get_size_of_TF(i) - 1) );
                no_access_dna_forward.insert(reverse_interval_for_another_strand(forward_DNA.size(), std::make_pair(new_coord, new_coord + tfs.get_size_of_TF(i) - 1) ));

            }  else {
                no_access_dna_revcompl.insert(std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1) );
                no_access_dna_forward.insert(reverse_interval_for_another_strand(forward_DNA.size(), std::make_pair(current_coordinate, current_coordinate + tfs.get_size_of_TF(i) - 1)) );

                //std::cout << "Can not move!\n";
            }
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
    int total_number_of_occupied_specific_sites = 0;
    double time_of_previous_change = 0.0;
    double current_time = 0.0;
    
    while ((total_number_of_occupied_specific_sites > bound_for_number_of_specific_sites && current_time - time_of_previous_change < time_before_stabilization) || (total_number_of_occupied_specific_sites <= bound_for_number_of_specific_sites && current_time - time_of_previous_change < upper_bound_for_time) || total_number_of_occupied_specific_sites < 10) {
        //std::cout << "TIME LEFT " << time_before_stabilization - current_time + time_of_previous_change << "\n";
        std::pair<double, int> next_event = timeline.top();
        current_time = next_event.first;
        //std::cout << next_event.first << " " << next_event.second << "\n";
        if (timeline.empty()) {
            tfs.generate_next_characteristic_time_for_binded_TFs();
            tfs.generate_next_characteristic_time_for_unbinded_TFs();
        }
        if (next_event.second == 1 && tfs.num_of_unbinded_TFs() > 0) {
            //std::cout << number_of_steps_without_change << " operation with unbinded" << "\n";
            int next_tf = tfs.choose_next_unbinded_DNA_to_interact();
            if (next_tf >= 0) {
                bool result_of_binding = bind_tf_to_dna(next_tf, tfs, dnase);
                if (result_of_binding == true) {
                    generate_next_event_with_binded(tfs);
                    generate_next_event_with_unbinded(tfs);
                    tfs.get_tf_by_index(next_tf).time_of_binding = next_event.first;
                    time_for_unbinded_state.push_back(next_event.first - tfs.get_tf_by_index(next_tf).time_of_unbinding);
                }
                else generate_next_event_with_unbinded(tfs);
                bool result_of_specific_binding = test_for_specific_binding(next_tf, tfs);
                if (result_of_specific_binding) {
                    total_number_of_occupied_specific_sites++;
                    number_of_steps_without_change = 0;
                    time_of_previous_change = next_event.first;
                    times_of_changes.push_back(next_event.first);
                }
                
            }
        } else if (next_event.second == 2 && tfs.num_of_binded_TFs() > 0) {
            // std::cout << number_of_steps_without_change << " operation with binded" << "\n";
            int next_tf = tfs.choose_next_binded_DNA_to_interact();
            if (next_tf >= 0) {
                bool binded_specifically = tfs.get_tf_by_index(next_tf).is_binded_specifically();
                if (!binded_specifically) {
                    if (next_tf >= 0) {
                        bool result_of_specifi_binding = false;
                        for (int i = 0; i < number_of_one_dim_slidings_in_step; ++i) {
                            one_dimensional_slinding_of_TF(next_tf, tfs);
                            result_of_specifi_binding = test_for_specific_binding(next_tf, tfs);
                            if (result_of_specifi_binding) {
                                number_of_steps_without_change = 0;
                                total_number_of_occupied_specific_sites++;
                                time_of_previous_change = next_event.first;
                                times_of_changes.push_back(next_event.first);
                                break;
                            }
                        }
                        if (!result_of_specifi_binding) {
                            double weight_of_binding = 0.0;
                            std::string tf_type = tfs.get_tf_by_index(next_tf).type_name;
                            int current_coordinate = tfs.get_tf_by_index(next_tf).get_coordinate_in_sequence();
                            
                            if (tfs.get_tf_by_index(next_tf).is_binded_to_forward()) {
                                weight_of_binding = weights_of_binding_of_all_tfs_forward[tf_type][current_coordinate];
                            } else {
                                weight_of_binding = weights_of_binding_of_all_tfs_revcompl[tf_type][current_coordinate];
                            }
                            if (test_for_unbinding(weight_of_binding)) {
                                unbind_tf_from_dna(next_tf, tfs);
                                tfs.get_tf_by_index(next_tf).time_of_unbinding = next_event.first;
                                time_for_binded_state.push_back(next_event.first - tfs.get_tf_by_index(next_tf).time_of_binding);
                                generate_next_event_with_unbinded(tfs);
                            } else {
                                generate_next_event_with_binded(tfs);
                            }
                        }
                    }
                }
                else {
                    generate_next_event_with_binded(tfs);
                }
            }
        }
        timeline.pop();
        number_of_steps_without_change++;
    }
    
    std::cout << "\nAverage sliding length, by type: ";
    std::map<std::string, std::vector<double>> average_lengths = tfs.return_average_sliding_lengths_vectors();
    for (std::map<std::string, std::vector<double>>::iterator it = average_lengths.begin(); it != average_lengths.end(); ++it) {
        std::cout << "\nType " << it->first << " : ";
        double average = 0.0;
        for (int i = 0; i < it->second.size(); ++i) {
            average += it->second[i];
        }
        average /= it->second.size();
        std::cout << (int)average;
    }

    /*for (std::map<std::vector<bool>, bool>::iterator it = final_frequency_of_combinations.begin(); it != final_frequency_of_combinations.end(); ++it) {
        for (int i = 0; i < it->first.size(); i++) {
         std::cout << it->first[i];
         }
        std::cout << "\n";
    }*/
    double mean_for_unbinded = 0.0;
    for (int i = 0; i < time_for_unbinded_state.size(); i++) {
        mean_for_unbinded += time_for_unbinded_state[i];
    }
    mean_for_unbinded /= time_for_unbinded_state.size();
    double mean_for_binded = 0.0;
    for (int i = 0; i < time_for_binded_state.size(); i++) {
        mean_for_binded += time_for_binded_state[i];
    }
    mean_for_binded /= time_for_binded_state.size();
    std::cout << "\nAVERAGE TIME FOR UNBINDED " << mean_for_unbinded << "\n";
    std::cout << "AVERAGE TIME FOR BINDED " << mean_for_binded << "\n";
    
}

std::vector<bool> Cell::get_final_combo() {
    return specific_sites_binded_or_not;
}

std::vector<double> Cell::get_times_of_changes() {
    return times_of_changes;
}




