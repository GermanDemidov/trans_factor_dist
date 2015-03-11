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

class Cell {
public:
    Cell(std::map<std::string, std::string>, std::vector<int>);
    void generate_next_event();
    
private:
    std::map<std::string, std::string> DNA;
    std::vector<Transcription_factor> transcription_factors;
    
};

#endif /* defined(__trans_factors_distrib__cell__) */
