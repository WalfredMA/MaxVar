//
//  MLE.cpp
//  sample_analyzer
//
//  Created by Wangfei MA on 7/17/20.
//  Copyright © 2020 UCSF_Kwoklab. All rights reserved.
//

#include "MLE.hpp"

using namespace std;


bool MLE::log_Gaussian(vector<double> &frequencies, int binstart, float *coef_a, float *coef_b, float *coef_c)
{
    
    vector<float> values_y;
    
    vector<float> values_x;
    
    for (int i=0; i<frequencies.size(); i++)
    {
        
        values_x.push_back(i+binstart);
        
        values_y.push_back(log(frequencies[i]+1));
         
    }
    if (linear_fitting(values_x, values_y, coef_a, coef_b)<0.6 ){return 0;}

    return 1;
    
}

bool MLE::Poisson(vector<double> &frequencies, int binstart, float *result){
    
    
    size_t binsize=frequencies.size();
    
    long long highvalue=0;
    int highindex=0;
    
    for (int i=0;i<binsize;i++){
        
        if (frequencies[i]>highvalue){
            
            highvalue=frequencies[i];
            
            highindex=i;
        }
    }
    
    highindex+=binstart;
    
    vector<float> values_y;
    
    vector<float> values_x;
    
    for (int i=0; i<binsize; i++){
        
        values_x.push_back(i+binstart);
        
        values_y.push_back(log(frequencies[i]+1)+permutation_log(1, i+binstart));
        
    }
    
    float lambda;
    if (slope(values_x, values_y, &lambda )>0.8 && fabs((float)exp(lambda)-highindex)<0.1*highindex){
        
        (*result)=exp(lambda);
        
        return 1;
    }
    
    else {
        
        return 0;
    }
}


bool MLE::Gaussian(vector<double> &frequencies, int binstart, float *err, float *mean)
{
    
    size_t binsize=frequencies.size();
    
    vector<float> values_y;
    
    vector<float> values_x;
    
    for (int i=0; i<binsize; i++)
    {
        
        values_x.push_back(i);
        
        values_y.push_back(log(frequencies[i]+1));
        
    }
    
    float mode;
    
    if (quadratic_fitting(values_x, values_y, err, mean, &mode)<0.5 || *err>=0){return 0;}
    
    *err=-0.5/(*err);
    
    *mean=binstart+(*mean)*(*err);
    
    *err=sqrt(*err);
    
    return 1;
}



float MLE::quadratic_fitting(const vector<float>& x, const vector<float>& y, float *coef_a, float *coef_b, float *coef_c)
{
    
    size_t binsize = (y.size()<x.size())?y.size():x.size();
    
    if (binsize<3){return 0.0;}
    
    vector<float> increments;
    for(size_t i=0; i<binsize-1; i++)
    {
        
        increments.push_back((y[i+1]-y[i])/(x[i+1]-x[i]));
    }
    
    if (linear_fitting(x, increments,coef_a, coef_b)<0.8){return 0.0;}
    
    (*coef_b)-=*coef_a;
    (*coef_a)*=0.5;
    
    vector<float>residuals;
    for(size_t i=0; i<binsize; i++)
    {
        
        residuals.push_back(y[i]-((*coef_a)*x[i]*x[i]+(*coef_b)*x[i]));
        
    }
    
    float avgY=accumulate(y.begin(), y.begin()+binsize, 0.0)/binsize;
    
    *coef_c=accumulate(residuals.begin(), residuals.begin()+binsize, 0.0)/binsize;
    
    
    float regression=0;
    float Var_y=0;
    for(size_t i=0; i<binsize; i++)
    {
        
        regression+=(residuals[i]-*coef_c)*(residuals[i]-*coef_c);
        Var_y+=(y[i]-avgY)*(y[i]-avgY);
    }
    
    if (!Var_y) return 1.0;
    
    return 1.0-regression/Var_y;
    
}


float MLE::linear_fitting(const vector<float>& x, const vector<float>& y, float *slope, float *intercept)
{
    
    size_t binsize = (y.size()<x.size())?y.size():x.size();
    
    if (binsize<2){return 0;}
    
    double avgX = accumulate(x.begin(), x.begin()+binsize, 0.0) / binsize;
    double avgY = accumulate(y.begin(), y.begin()+binsize, 0.0) / binsize;
    
    
    double Var_xy= 0.0;
    double Var_xx = 0.0;
    double Var_yy = 0.0;
    
    for(size_t i=0; i<binsize; ++i)
    {
        Var_xy += (x[i] - avgX) * (y[i] - avgY);
        Var_xx += (x[i] - avgX) * (x[i] - avgX);
    }
    
    if (!Var_xx)
    {
        
        return 0.0;
    }
    
    else
    {
        
        *slope=Var_xy/Var_xx; 
        
        *intercept=avgY-avgX*(*slope);
        
        for(size_t i=0; i<binsize; i++)
        {
            
            Var_yy+=(y[i]-avgY)*(y[i]-avgY);
        }
        
        if (!Var_yy) return 1.0;
        
        
        return Var_xy*Var_xy/(Var_yy*Var_xx);
        
    }
    
}


float MLE::slope(const vector<float>& x, const vector<float>& y, float *slope )
{
    
    size_t binsize = (y.size()<x.size())?y.size():x.size();
    
    double avgX = accumulate(x.begin(), x.begin()+binsize, 0.0) / binsize;
    double avgY = accumulate(y.begin(), y.begin()+binsize, 0.0) / binsize;
    
    
    double Var_xy = 0.0;
    double Var_yy = 0.0;
    double Var_xx = 0.0;
    
    for(size_t i=0; i<binsize; ++i)
    {
        Var_xy += (x[i] - avgX) * (y[i] - avgY);
        Var_xx += (x[i] - avgX) * (x[i] - avgX);
    }
    
    if (!Var_xx)
    {
        
        return 0.0;
    }
    
    else
    {
        
        *slope=Var_xy/Var_xx;
        
        for(size_t i=0; i<binsize; i++)
        {
            
            Var_yy+=(y[i]-avgY)*(y[i]-avgY);
        }
        
        if (!Var_yy) return 1.0;
        
        return Var_xy*Var_xy/(Var_yy*Var_xx);
        
    }
}


float MLE::permutation_log(int start,int end)
{
    
    return factorial_log(end)-factorial_log(start);
}


float MLE::factorial_log(int end)
{
    
    static float results_storage[100];
    static bool results_index[100]={0};
    
    if (end<100 && results_index[end])
    {
        
        return results_storage[end];
        
    }
    
    float factorial=0.0;
    
    for (int i=1;i<=end;i++)
    {
        factorial+=log(i);
    }
    
    
    if (end<100)
    {
        
        results_index[end]=1;
        results_storage[end]= factorial;
        
    }
    
    return  factorial;
}

