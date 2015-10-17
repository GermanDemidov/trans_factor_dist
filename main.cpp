//
//  main.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 10.02.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include <iostream>
#include <thread>
#include <future>
#include "parser.h"
#include "transcription_factor.h"
#include "auxuilary.h"
#include "cell.h"

class Solver {
    std::map<char, double> background_probabilities;
    std::map<std::string, double> concentrations;
public:
    Solver(std::string fasta_filename, std::string pcm_filename, std::string gene_to_study,
           std::string first_prot, std::string second_prot, std::string path_to_files, std::string file_to_out, int number_of_simulation) {
        // initialization of background frequencies
        background_probabilities.insert(std::make_pair('A', 0.28768819278776));
        background_probabilities.insert(std::make_pair('C', 0.21231180721224));
        background_probabilities.insert(std::make_pair('G', 0.21231180721224));
        background_probabilities.insert(std::make_pair('T', 0.28768819278776));
        // number of transcription factors
        const int number_of_transcription_factors = 100;
        
        number_of_simulation /= 4;
        
        double first_concentration = 5.0;
        double second_concentration = 5.0;
        // initialization of concentrations
        concentrations.insert(std::make_pair(first_prot, first_concentration));
        concentrations.insert(std::make_pair(second_prot, second_concentration));
        
        const int number_of_simulations = 1;
        
        // file with the data about euchromatin/
        //std::string euchromatin_input = gene_to_study + "_dnaseAccS05.btrack";
        std::string euchromatin_input = path_to_files + gene_to_study + "_dnaseAccS05.btrack";
        std::string output_file_name = path_to_files + file_to_out;
        std::string output_file_name_raw = path_to_files + "raw_" + file_to_out;
        
        std::ofstream outfile_raw;
        outfile_raw.open(output_file_name_raw, std::ios::app);
        
        std::map<int, std::vector<double> > times_of_steps_for_each;
        for (int i = 0; i < 1000; i++) {
            times_of_steps_for_each[i];
        }
        
        
        // parse FASTA-file, creates hash map "sequences" inside fasta_sequences
        // key - fasta id, value - fasta sequence
        Parser_fasta fasta_sequences(fasta_filename);
        
        // parse PCM-file, creates PWM - hash map
        // key - TF name, value - vector of doubles
        Parser_pcm pcm_for_TFs(pcm_filename);
        std::vector<std::string> protein_names = pcm_for_TFs.return_protein_names();
        pcm_for_TFs.calculate_pwm(background_probabilities);
        
        std::map<std::string, std::vector<int> > highest_binding_sites;
        std::map<std::string, std::string> seqs = fasta_sequences.return_sequences();
        for(std::map<std::string, std::string>::iterator iterator = seqs.begin(); iterator != seqs.end(); iterator++) {
            std::vector<int> initialization_vector;
            highest_binding_sites.insert(std::make_pair(iterator->first, initialization_vector));
        }
        
        Parser_dnase_acc euchromatin(euchromatin_input);
        std::vector<double> binding_weight;
        
        std::vector<Transcription_factor> transcription_factors_to_initialize;
        
        for (std::map<std::string, std::string>::iterator iterator = seqs.begin(); iterator != seqs.end(); iterator++) {
            bool flag_of_gene_to_study = false;
            for (int k = 0; k < gene_to_study.size(); ++k) {
                if (iterator->first.at(k + 1) == gene_to_study.at(k)) {
                    flag_of_gene_to_study = true;
                } else {
                    flag_of_gene_to_study = false;
                    break;
                }
            }
            std::vector<std::string> prots_to_study;
            
            if (flag_of_gene_to_study == true) {
                std::string sequence_for_gene_to_study = iterator->second;
                std::string promoter_sequence = sequence_for_gene_to_study.substr(0, 12000);
                std::string revcompl_promoter_sequence = reverse_compliment(sequence_for_gene_to_study.substr(0, 12000));
                prots_to_study.push_back(">" + first_prot);
                prots_to_study.push_back(">" + second_prot);
                
                
                for (int i = 0; i != protein_names.size(); ++i) {
                    if (protein_names[i] == ">" + first_prot || protein_names[i] == ">" + second_prot) {
                        std::cout << protein_names[i] << " ";
                        Transcription_factor tf(0, pcm_for_TFs.return_profile_for_protein_pwm(protein_names[i]), protein_names[i]);
                        transcription_factors_to_initialize.push_back(tf);
                    }
                }
                
                double sum_of_concentrations = 0.0;
                for (std::map<std::string, double>::iterator it = concentrations.begin(); it != concentrations.end(); ++it) {
                    sum_of_concentrations += it->second;
                }
                
                double part =(first_concentration / sum_of_concentrations);
                int bound_between_different_types = (int) (part * (1.0 * (double)number_of_transcription_factors));
                for (int k = 0; k < number_of_simulations; ++k) {
                    std::vector<Transcription_factor> tfs_instances;
                    for (int j = 0; j < number_of_transcription_factors; j++) {
                        if (j < bound_between_different_types) {
                            tfs_instances.push_back(Transcription_factor(j, pcm_for_TFs.return_profile_for_protein_pwm(prots_to_study[0]), prots_to_study[0]));
                        } else {
                            tfs_instances.push_back(Transcription_factor(j, pcm_for_TFs.return_profile_for_protein_pwm(prots_to_study[1]), prots_to_study[1]));
                        }
                    }
                    std::map<std::vector<bool>, int> final_combinations;
                    std::map<std::vector<bool>, bool> final_combo;
                    std::vector<bool> quasi_stable_combo;
                    
                    std::vector<double> times_of_steps;
                    std::vector<double> times_of_steps1;
                    std::vector<double> times_of_steps2;
                    std::vector<double> times_of_steps3;
                    
                    
                    
                    
                    Cell new_cell = *new Cell(promoter_sequence, revcompl_promoter_sequence, k,
                                              transcription_factors_to_initialize, euchromatin, concentrations, tfs_instances);
                    /*final_combo = new_cell.get_final_frequency_of_combinations();
                     for (std::map<std::vector<bool>, bool>::iterator it = final_combo.begin(); it != final_combo.end(); ++it){
                     std::vector<bool> final_comb = it->first;
                     outfile_raw.open(output_file_name_raw, std::ios::app);
                     for (int s = 0; s < final_comb.size(); s++) {
                     std::cout << final_comb[s];
                     outfile_raw << final_comb[s];
                     }
                     std::cout << "\n";
                     outfile_raw << "\n";
                     outfile_raw.close();
                     }*/
                    
                    quasi_stable_combo = new_cell.get_final_combo();
                    
                    final_combinations[quasi_stable_combo] = 1;
                    
                    outfile_raw.open(output_file_name_raw, std::ios::app);
                    outfile_raw.close();
                    
                    outfile_raw.open(output_file_name_raw, std::ios::app);
                    for (int s = 0; s < quasi_stable_combo.size(); s++) {
                        std::cout << quasi_stable_combo[s];
                        outfile_raw << quasi_stable_combo[s];
                    }
                    std::cout << "\n";
                    outfile_raw << "\n";
                    outfile_raw.close();
                    
                    
                    
                    times_of_steps = new_cell.get_times_of_changes();
                    if (times_of_steps.size() > 2) {
                        for (int k = 0; k < times_of_steps.size() - 1; k++) {
                            times_of_steps_for_each[k].push_back(times_of_steps[k+1] - times_of_steps[k]);
                        }
                    }
                    
                    
                    
                    new_cell.null_everything(euchromatin);
                    Cell new_cell1(new_cell);
                    Cell new_cell2(new_cell);
                    Cell new_cell3(new_cell);
                    
                    // block with parallelization
                    for (int l = 0; l < number_of_simulation; l++) {
                        std::cout << "\n\nSTEP " << l << "\n";
                        std::vector<std::string> protein_names;
                        for (int i = 0; i < transcription_factors_to_initialize.size(); ++i) {
                            protein_names.push_back(transcription_factors_to_initialize[i].type_name);
                        }
                        
                        new_cell.null_everything(euchromatin);
                        new_cell1.null_everything(euchromatin);
                        new_cell2.null_everything(euchromatin);
                        new_cell3.null_everything(euchromatin);
                        
                        Transcription_factors_in_cell tfs(protein_names, tfs_instances);
                        Transcription_factors_in_cell tfs1(protein_names, tfs_instances);
                        Transcription_factors_in_cell tfs2(protein_names, tfs_instances);
                        Transcription_factors_in_cell tfs3(protein_names, tfs_instances);
                        
                        tfs.null_everything(protein_names);
                        tfs1.null_everything(protein_names);
                        tfs2.null_everything(protein_names);
                        tfs3.null_everything(protein_names);
                        
                        std::thread thr(&Cell::start_simulation, &new_cell, std::ref(tfs), std::ref(euchromatin));
                        std::thread thr1(&Cell::start_simulation, &new_cell1, std::ref(tfs1), std::ref(euchromatin));
                        std::thread thr2(&Cell::start_simulation, &new_cell2, std::ref(tfs2), std::ref(euchromatin));
                        std::thread thr3(&Cell::start_simulation, &new_cell3, std::ref(tfs3), std::ref(euchromatin));
                        
                        thr.join();
                        thr1.join();
                        thr2.join();
                        thr3.join();
                        /*new_cell.start_simulation(tfs, euchromatin);
                         new_cell1.start_simulation(tfs1, euchromatin);
                         new_cell2.start_simulation(tfs2, euchromatin);
                         new_cell3.start_simulation(tfs3, euchromatin);*/
                        
                        // new_cell.start_simulation(tfs, euchromatin);
                        /*final_combo = new_cell.get_final_frequency_of_combinations();
                         std::map<std::vector<bool>, bool> final_combo1;
                         final_combo1 = new_cell1.get_final_frequency_of_combinations();
                         std::map<std::vector<bool>, bool> final_combo2;
                         final_combo2 = new_cell2.get_final_frequency_of_combinations();
                         std::map<std::vector<bool>, bool> final_combo3;
                         final_combo3 = new_cell3.get_final_frequency_of_combinations();
                         std::vector<std::map<std::vector<bool>, bool >> final_combos;
                         final_combos.push_back(final_combo);
                         final_combos.push_back(final_combo1);
                         final_combos.push_back(final_combo2);
                         final_combos.push_back(final_combo3);*/
                        
                        quasi_stable_combo = new_cell.get_final_combo();
                        std::vector<bool> quasi_stable_combo1;
                        quasi_stable_combo1 = new_cell1.get_final_combo();
                        std::vector<bool> quasi_stable_combo2;
                        quasi_stable_combo2 = new_cell2.get_final_combo();
                        std::vector<bool> quasi_stable_combo3;
                        quasi_stable_combo3 = new_cell3.get_final_combo();
                        std::vector<std::vector<bool> > final_combos_quasi_stable;
                        final_combos_quasi_stable.push_back(quasi_stable_combo);
                        final_combos_quasi_stable.push_back(quasi_stable_combo1);
                        final_combos_quasi_stable.push_back(quasi_stable_combo2);
                        final_combos_quasi_stable.push_back(quasi_stable_combo3);
                        
                        
                        for (int m = 0; m != final_combos_quasi_stable.size(); m++) {
                            if (final_combinations.find(final_combos_quasi_stable[m]) != final_combinations.end()) {
                                final_combinations[final_combos_quasi_stable[m]] += 1;
                            } else {
                                final_combinations[final_combos_quasi_stable[m]] = 1;
                            }
                            outfile_raw.open(output_file_name_raw, std::ios::app);
                            for (int s = 0; s < final_combos_quasi_stable[m].size(); s++) {
                                outfile_raw << final_combos_quasi_stable[m][s];
                            }
                            outfile_raw << "\n";
                            outfile_raw.close();
                        }
                        
                        
                        for (int l = 0; l < final_combos_quasi_stable.size(); l++) {
                            int number_of_specific_sites = 0;
                            /*for (int d = 0; d < final_combos[l].size(); d++) {
                             if (final_combos[l][d])
                             number_of_specific_sites++;
                             }*/
                            //std::cout << "Total amount of occupied sites: " << number_of_specific_sites << "\n";
                            /*if (final_combinations.find(final_combos[l]) != final_combinations.end()) {
                             final_combinations[final_combos[l]] += 1;
                             } else {
                             final_combinations[final_combos[l]] = 1;
                             }*/
                            /*for (std::map<std::vector<bool>, bool>::iterator it = final_combos[l].begin(); it != final_combos[l].end(); ++it){
                             outfile_raw.open(output_file_name_raw, std::ios::app);
                             for (int s = 0; s < it->first.size(); s++) {
                             outfile_raw << it->first[s];
                             }
                             outfile_raw << "\n";
                             outfile_raw.close();
                             }*/
                            times_of_steps = new_cell.get_times_of_changes();
                            if (times_of_steps.size() > 3) {
                                for (int k = 0; k < times_of_steps.size() - 1; k++) {
                                    times_of_steps_for_each[k].push_back(times_of_steps[k+1] - times_of_steps[k]);
                                }
                            }
                            times_of_steps1 = new_cell1.get_times_of_changes();
                            if (times_of_steps1.size() > 3) {
                                for (int k = 0; k < times_of_steps1.size() - 1; k++) {
                                    times_of_steps_for_each[k].push_back(times_of_steps1[k+1] - times_of_steps1[k]);
                                }
                            }
                            times_of_steps2 = new_cell2.get_times_of_changes();
                            if (times_of_steps2.size() > 3) {
                                for (int k = 0; k < times_of_steps2.size() - 1; k++) {
                                    times_of_steps_for_each[k].push_back(times_of_steps2[k+1] - times_of_steps2[k]);
                                }
                            }
                            times_of_steps3 = new_cell3.get_times_of_changes();
                            if (times_of_steps2.size() > 3) {
                                for (int k = 0; k < times_of_steps3.size() - 1; k++) {
                                    times_of_steps_for_each[k].push_back(times_of_steps3[k+1] - times_of_steps3[k]);
                                }
                            }
                        }
                    }
                    
                    std::ofstream outfile (output_file_name);
                    for (std::map<std::vector<bool>, int>::iterator it = final_combinations.begin(); it != final_combinations.end(); ++it) {
                        for (int i = 0; i < it->first.size(); i++) {
                            std::cout << it->first[i];
                            outfile << it->first[i];
                        }
                        std::cout << " " << it->second << "\n";
                        outfile << " " << it->second << "\n";
                    }
                }
            }
        }
        std::cout << "\n";
        /*int counter_for_steps = 0;
         for (std::map<int, std::vector<double>>::iterator it = times_of_steps_for_each.begin(); it != times_of_steps_for_each.end(); ++it) {
         if (!it->second.empty()) {
         std::cout << "step" << counter_for_steps << ",";
         //std::cout << '(';
         for (int i = 0; i < it->second.size(); i++) {
         std::cout << it->second[i];
         if (i < it->second.size() - 1) {
         std::cout <<  ",";
         }
         }
         //std::cout << ")";
         std::cout << "\n";
         counter_for_steps++;
         }
         }*/
        /*for(std::map<std::string, std::vector<int>>::iterator iterator = highest_binding_sites.begin(); iterator != highest_binding_sites.end(); iterator++) {
         std::cout << iterator->first << "\n";
         for (int i = 0; i < iterator->second.size(); i++) {
         std::cout << iterator->second[i] << " ";
         }
         std::cout << "\n\n";
         }*/
        /*int count_of_hops = 0;
         int count_of_long_hops = 0;
         std::vector<double> test;
         for (int i = 0; i < 1000000; i++) {
         double x = generate_next_coordinate_change();
         if (x > 15) count_of_hops++;
         if (x > 25) count_of_long_hops++;
         
         test.push_back(x);
         }
         std::cout << count_of_hops << "\n";
         std::cout << count_of_long_hops << "\n";
         std::cout << "\nMEAN: " << sm(test) << " SD: " << sd(test) << "\n";*/
        outfile_raw.close();
    }
    
};

int main(int argc, const char * argv[]) {
    // it is a working directory. You should replace this with your directory.
    std::string path_to_files = "/Users/german/Desktop/Gurskiy/project/trans_factors_distrib/";
    int number_of_simulations = 100;
    std::string pcm_file_input = path_to_files + argv[1];
    std::string fasta_file_input = path_to_files + argv[2];
    std::string gene_to_study = argv[3];
    std::string protein_to_study_first = argv[4]; // activator
    std::string protein_to_study_second = argv[5]; // represor
    std::string output_file_name = argv[6];
    std::cout << "Starting with PCM " << pcm_file_input << "\n" << "and input fasta " << fasta_file_input << "...\n";
    
    Solver solve(fasta_file_input, pcm_file_input, gene_to_study, protein_to_study_first, protein_to_study_second, path_to_files, output_file_name, number_of_simulations);
    
    return 0;
}
