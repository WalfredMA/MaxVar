//
//  MDTag.cpp
//  testcpp
//
//  Created by Wangfei MA on 6/28/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "MDTag.hpp"

using namespace std;

bool MDTag::GrepMD(string& des, const char *mdtext) { 
    
    const char *bar_search=strstr(mdtext, "MD:Z:");
    
    if (bar_search == NULL ){
        
        return false;
    }
    
    char current_mdtext;
    
    int current_segment=0;
    
    int current_type=0;
    
    for (size_t i = (bar_search-mdtext+5); i < strlen(mdtext); ++i){
        
        current_mdtext=mdtext[i];
        if (current_mdtext=='^'){
            
            des.insert(des.end(),current_segment,'M');
            
            current_segment=0;
            
            current_type=1;
        }
        
        else if (current_mdtext>='0' && current_mdtext<='9'){
            
            current_type=0;
            
            current_segment*=10;
            
            current_segment+=(int)(current_mdtext-'0');
            
        }
        
        else{

            des.insert(des.end(),current_segment,'M');
            
            current_segment=0;
            
            if (current_type){
                
                des.push_back('D');
            }
            
            else {
                
                des.push_back(current_mdtext);
            }
        }
    }
    
    des.insert(des.end(),current_segment,'M');
    
    return true;
}


