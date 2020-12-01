//
//  main.cpp
//  covar_network
//
//  Created by Wangfei MA on 7/5/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include <iostream>
#include <string>

#include "covar_networks.hpp"
#include "mole_network.hpp"
#include "snp_network.hpp"

int main(int argc, const char * argv[]) {
    
    string inputfile=string(argv[1]);
    
    string outputfile=string(argv[2]);
    
    float cutoff=0.4;
    if (argc>3) cutoff=atof(argv[3]);
    
    cout << "Starting Program" <<endl;
    
    covar_networks covar_networks;
    
    covar_networks.init(inputfile);
    
    vector<float> cutoffs{cutoff,cutoff,cutoff};
    covar_networks.run(cutoffs);
    
    covar_networks.output(outputfile);
    
    return 0;
}
