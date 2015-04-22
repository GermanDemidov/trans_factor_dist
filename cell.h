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
#include "transcription_factor.h"
#include "parser.h"

class Cell {
public:
    Cell(std::string&, std::string&, int, std::vector<Transcription_factor>&, Parser_dnase_acc&);
    void generate_next_event();
    void generate_TF_appearance_time();
    void find_specific_binding_sites(Transcription_factor, Parser_dnase_acc&, bool);
    
private:
    int cell_id;
    
    std::string forward_DNA;
    std::string reverse_DNA;
    std::vector<Transcription_factor> transcription_factors;
    std::map<std::string, std::vector<double> > weights_of_binding_of_all_tfs_forward;
    std::map<std::string, std::vector<double> > weights_of_binding_of_all_tfs_revcompl;
    
    std::map<std::string, std::map<int, bool> > specific_binding_sites_forward;
    std::map<std::string, std::map<int, bool> > specific_binding_sites_backward;
    
};

#endif /* defined(__trans_factors_distrib__cell__) */
