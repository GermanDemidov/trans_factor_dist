//
//  main.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 10.02.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include <iostream>
#include "parser.h"
#include "transcription_factor.h"
#include "auxuilary.h"
#include "cell.h"

class Solver {
    std::map<char, double> background_probabilities;
    std::map<std::string, double> concentrations;
public:
    Solver(std::string fasta_filename, std::string pcm_filename, std::string gene_to_study,
           std::string first_prot, std::string second_prot) {
        // initialization of background frequencies
        background_probabilities.insert(std::make_pair('A', 0.28768819278776));
        background_probabilities.insert(std::make_pair('C', 0.21231180721224));
        background_probabilities.insert(std::make_pair('G', 0.21231180721224));
        background_probabilities.insert(std::make_pair('T', 0.28768819278776));
        
        // initialization of concentrations
        concentrations.insert(std::make_pair(first_prot, 10.7));
        concentrations.insert(std::make_pair(second_prot, 14.8));
        
        const int number_of_simulations = 100;
        
        // file with the data about euchromatin/
        std::string euchromatin_input = "/Users/german/Desktop/Gurskiy/project/trans_factors_distrib/trans_factors_distrib/dnase/" + gene_to_study + "_dnaseAccS05.btrack";

        
        // parse FASTA-file, creates hash map "sequences" inside fasta_sequences
        // key - fasta id, value - fasta sequence
        Parser_fasta fasta_sequences(fasta_filename);
        //fasta_sequences.debug_print();
        
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


                //std::cout << tf.motif[0].at('A') << "\n";
                //typedef std::map<std::string, std::string>::iterator it_type;
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
            if (flag_of_gene_to_study == true) {
                for (int i = 0; i != protein_names.size(); ++i) {
                    if (protein_names[i] == ">" + first_prot || protein_names[i] == ">" + second_prot) {
                        std::cout << protein_names[i] << " ";
                        Transcription_factor tf(0, pcm_for_TFs.return_profile_for_protein_pwm(protein_names[i]), protein_names[i]);
                        transcription_factors_to_initialize.push_back(tf);
                    }
                        std::string sequence_for_gene_to_study = iterator->second;
                        std::string promoter_sequence = sequence_for_gene_to_study.substr(0, 12000);
                        std::string revcompl_promoter_sequence = reverse_compliment(sequence_for_gene_to_study.substr(0, 12000));
                        std::cout << promoter_sequence << "\n";
                        std::cout << revcompl_promoter_sequence << "\n";
                        
                        for (int k = 0; k < number_of_simulations; ++k) {
                            Cell(promoter_sequence, revcompl_promoter_sequence, k,
                                 transcription_factors_to_initialize, euchromatin);
                        }
                        
                    /*highest_binding_sites[iterator->first].push_back(0);
                    std::string sequence_for_prot = iterator->second;
                    for (int j = 0; j != 12000; ++j) {
                        tf.change_coordinate_in_sequence(j);
                        double weight_of_current_binding = tf.calculate_weight_of_binding(sequence_for_prot);
                        binding_weight.push_back(weight_of_current_binding);
                        if (weight_of_current_binding > 3.0 && euchromatin.is_in_interval(std::make_pair(j, j + tf.get_size()))) {
                            long last_elem = highest_binding_sites[iterator->first].size();
                            highest_binding_sites[iterator->first][last_elem - 1]++;
                        }
                    }
                
                    std::string sequence_for_prot_rev = reverse_compliment(iterator->second);
                    for (int j = (int)sequence_for_prot.size() - 12000 - tf.get_size(); j != sequence_for_prot.size() - tf.get_size(); ++j) {
                        tf.change_coordinate_in_sequence(j);
                        double weight_of_current_binding = tf.calculate_weight_of_binding(sequence_for_prot_rev);
                        binding_weight.push_back(weight_of_current_binding);
                        if (weight_of_current_binding > 3.0 && euchromatin.is_in_interval(std::make_pair(j - 12000 - tf.get_size(), j - 12000))) {
                            long last_elem = highest_binding_sites[iterator->first].size();
                            highest_binding_sites[iterator->first][last_elem - 1]++;
                        }
                    }*/
                }
            }
        }
        std::cout << "\n";
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
        
    }
};

int main(int argc, const char * argv[]) {
    // insert code here...
    std::string path_to_files = "/Users/german/Desktop/Gurskiy/project/trans_factors_distrib/";
    std::string pcm_file_input = path_to_files + argv[1];
    std::string fasta_file_input = path_to_files + argv[2];
    std::string gene_to_study = argv[3];
    std::string protein_to_study_first = argv[4];
    std::string protein_to_study_second = argv[5];
    std::cout << "Starting with PCM " << pcm_file_input << "\n" << "and input fasta " << fasta_file_input << "...\n";
    
    Solver solve(fasta_file_input, pcm_file_input, gene_to_study, protein_to_study_first, protein_to_study_second);
    
    
    
    return 0;
}
