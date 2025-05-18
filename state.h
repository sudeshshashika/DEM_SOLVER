#pragma once
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream>
#include <string>
#include <fstream>
#include "Particle.h"
#include "Cell.h"

class state{

public:
       
    std::string  part_input_file;
    double       time_end;              
    double       timestep_length;      
    int          part_out_freq;         
    std::string  part_out_name_base;    
    int          vtk_out_freq;          
    std::string  vtk_out_name_base;     
    int          cl_workgroup_1dsize;   
    int          cl_workgroup_3dsize_x; 
    int          cl_workgroup_3dsize_y; 
    int          cl_workgroup_3dsize_z;                 
    double       k_s;                   
    double       k_dn;                 
    double       x_min;                 
    double       y_min;                 
    double       z_min;                 
    double       x_max;                 
    double       y_max;                 
    double       z_max;                 
    double       x_n;                   
    double       y_n;                   
    double       z_n;                   
    int          reflect_x;             
    int          reflect_y;             
    int          reflect_z;     
    double       g_x;                   
    double       g_y;                   
    double       g_z;
    double       len_x;
    double       len_y;
    double       len_z;
    int          particleCount;          


void setState(std::string infile1Name){
    std::string temp; 
    std::fstream infile1(infile1Name);
    infile1 >> temp >> part_input_file;
    infile1 >> temp >> time_end;
    infile1 >> temp >> timestep_length;
    infile1 >> temp >> part_out_freq;
    infile1 >> temp >> part_out_name_base;   
    infile1 >> temp >> vtk_out_freq;     
    infile1 >> temp >> vtk_out_name_base;  
    infile1 >> temp >> cl_workgroup_1dsize;
    infile1 >> temp >> cl_workgroup_3dsize_x; 
    infile1 >> temp >> cl_workgroup_3dsize_y;
    infile1 >> temp >> cl_workgroup_3dsize_z; 
    infile1 >> temp >> k_s; 
    infile1 >> temp >> k_dn; 
    infile1 >> temp >> x_min;
    infile1 >> temp >> y_min;
    infile1 >> temp >> z_min; 
    infile1 >> temp >> x_max;
    infile1 >> temp >> y_max;
    infile1 >> temp >> z_max; 
    infile1 >> temp >> x_n;
    infile1 >> temp >> y_n; 
    infile1 >> temp >> z_n;
    infile1 >> temp >> reflect_x;
    infile1 >> temp >> reflect_y;
    infile1 >> temp >> reflect_z;
    infile1 >> temp >> g_x;
    infile1 >> temp >> g_y;
    infile1 >> temp >> g_z;
    len_x   = abs(x_max - x_min) / x_n;
    len_y   = abs(y_max - y_min) / y_n;
    len_z   = abs(z_max - z_min) / z_n;
    infile1.close();
}
void setParticleCount(std::string infile2Name){
    std::fstream infile2(infile2Name);
    infile2 >> particleCount;
    infile2.close();
    
}
void printState(){
       std::cout<< part_input_file<<std::endl;        
       std::cout<< time_end<<std::endl;
       std::cout<< timestep_length<<std::endl;
       std::cout<< part_out_freq<<std::endl;
       std::cout<< part_out_name_base   <<std::endl;
       std::cout<< vtk_out_freq     <<std::endl;
       std::cout<< vtk_out_name_base <<std::endl;
       std::cout<< cl_workgroup_1dsize<<std::endl;
       std::cout<< cl_workgroup_3dsize_x <<std::endl;
       std::cout<< cl_workgroup_3dsize_y<<std::endl;
       std::cout<< cl_workgroup_3dsize_z <<std::endl;
       std::cout<< k_s <<std::endl;
       std::cout<< k_dn<<std::endl;
       std::cout<< x_min<<std::endl;
       std::cout<< y_min<<std::endl;
       std::cout<< z_min<<std::endl;
       std::cout<< x_max<<std::endl;
       std::cout<< y_max<<std::endl;
       std::cout<< z_max <<std::endl;
       std::cout<< x_n<<std::endl;
       std::cout<< y_n <<std::endl;
       std::cout<< z_n<<std::endl;
       std::cout<< reflect_x<<std::endl;
       std::cout<< reflect_y<<std::endl;
       std::cout<< reflect_z<<std::endl;
       std::cout<< g_x<<std::endl;
       std::cout<< g_y<<std::endl;
       std::cout<< g_z<<std::endl;
       std::cout<< particleCount<<std::endl;
       std::cout<<len_x<<std::endl;
       std::cout<<len_y<<std::endl;
       std::cout<<len_z<<std::endl;
}

};

