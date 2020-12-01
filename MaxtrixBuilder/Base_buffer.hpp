//
//  Base_buffer.hpp
//  testcpp
//
//  Created by Wangfei MA on 8/26/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#ifndef Base_buffer_hpp
#define Base_buffer_hpp

#include <stdio.h>

#include "Config.h" 


class Buff_ctrl
{
    
public:
    
    Buff_ctrl(){};
    ~Buff_ctrl(){};
    
    void update(int init_index,int column, int barindex, char qual, short extra_data);
    
    void update(int init_index, int column, int barindex, char qual, std::intpair extra_data);
    
    void update(int init_index, int column, int barindex, char qual, int extra_data);
    
    void Release(int release_index);
    
    void ReleaseAllBases(int release_size);
    
    void resize(Base* move_to, int target_allele, int addsize);
    
    void move(int movesize);
    
    void transfer(int old_pos, int new_pos, int target_allele, std::vector<int>& moles);
    
    Base* Lock_Base(int index);
    
    void Mini_Lock(int index);
    
    void Mini_unLock(int index);
    
    
    
};


#endif /* Base_buffer_hpp */
