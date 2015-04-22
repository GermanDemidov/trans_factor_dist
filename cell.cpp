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
           std::vector<Transcription_factor>& types_of_TFs_for_preprocessing, Parser_dnase_acc& dnase) {
    
    forward_DNA = for_DNA;
    reverse_DNA = rev_DNA;
    cell_id = id_of_cell;
    
    for (int i = 0; i < types_of_TFs_for_preprocessing.size(); ++i) {
        
    }
    
}

void Cell::find_specific_binding_sites(Transcription_factor tf, Parser_dnase_acc& dnase, bool forward) {
    if (forward) {
        
    } else {
        
    }
}
