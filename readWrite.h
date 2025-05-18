#pragma once
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>
#include "Particle.h"
#include "state.h"


void reading(const std::string &path, Particle* Particles){
    std::ifstream infile3{path};
    std::string temp{};
    int Number_of_particles;
    infile3 >> Number_of_particles;
    for(int i =0; i< Number_of_particles; ++i){
        infile3 >> temp; 
        std::istringstream iss{temp};
        if(iss >> Particles[i].m){
            iss >> Particles[i].m;
        }else{
            Particles[i].m = 2147483647;//INFINITY;
        }
        infile3 >> Particles[i].r >> Particles[i].d[0] >> Particles[i].d[1] >> Particles[i].d[2] 
                >> Particles[i].v[0] >> Particles[i].v[1] >> Particles[i].v[2] ;
    }
    infile3.close();
}

//writing to out file
void WriteOut(Particle* P, std::string path, int ITERATION, int part_out_freq, std::string part_out_name_base, int N){
    
    std::string part_out_PATH = path;                //change folder path
    std::string file_extension     = ".out";
    if( (ITERATION % part_out_freq)==0){
        std::ofstream outfile{part_out_PATH+part_out_name_base + std::to_string(ITERATION/part_out_freq)+ file_extension};
        outfile << N <<"\n";
        for(int i = 0; i < N; ++i){
            outfile<< std::fixed<< P[i].m << " " << P[i].d[0]  << " " << P[i].d[1]   << " " << P[i].d[2]
                                          << " " << P[i].v[0] << " " << P[i].v[1]  << " " << P[i].v[2] <<"\n";
        }
        outfile.close();   
    }
}
//writing to vtk file
void WriteVTK(Particle* P, std::string path, int ITERATION, int vtk_out_freq, std::string vtk_out_name_base, int N){

    std::string part_out_PATH = path;                //change folder path
    std::string file_extension     = ".vtk";

    if( (ITERATION % vtk_out_freq)==0){
        std::ofstream outfile{part_out_PATH+ vtk_out_name_base + std::to_string(ITERATION/vtk_out_freq)+ file_extension};
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "Particle data\n";
        outfile << "ASCII\n";
        outfile << "DATASET UNSTRUCTURED_GRID\n";
        outfile << "POINTS " << N << " double" <<"\n";

        for(int i = 0; i < N; ++i){
            outfile<< std::fixed<< P[i].d[0] << " " << P[i].d[1]  << " " << P[i].d[2] << "\n";
        }

        outfile << "\n"<<"\n";
        outfile << "CELLS " << "2 " << "4\n";
        outfile << "1 "<<"0"<<"\n";
        outfile << "1 "<<"1"<<"\n";
        outfile << "\n"<<"\n";
        outfile << "CELL_TYPES " << "2\n";
        outfile << "1"<<"\n";
        outfile << "1"<<"\n";
        outfile << "\n"<<"\n";
        outfile << "POINT_DATA " << N << "\n";
        outfile << "SCALARS mass double\n";
        outfile << "LOOKUP_TABLE default\n";

        for(int i = 0; i < N; ++i){
            outfile<< std::fixed<< P[i].m<< "\n";
        }
        outfile << "VECTORS v double\n";

        for(int i = 0; i < N; ++i){
            outfile<< std::fixed<< P[i].v[0] << " " << P[i].v[1]  << " " << P[i].v[2] << "\n";
        }

        outfile.close();
    }

}
