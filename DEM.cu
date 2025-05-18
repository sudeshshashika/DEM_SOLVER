/**main programme for the DEM SOLVER**/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "Particle.h"
#include "Cell.h"
#include "state.h"
#include "readWrite.h"
#include "kernels.h"
#include "device.h"

int main(){

    std::string path_par {"input/bounce.par"};
    std::string path_in  {"input/bounce.in"}; 
    std::string out_path {"output/"}; 
    double time{0};                                        
    int ITERATION{0};   
    state D; 
    device gpu;
    D.setState(path_par);
    D.setParticleCount(path_in);
    gpu.setThreadsPerBlock(1024);
    gpu.setBytes(D);
    gpu.setcNumBlocks();
    gpu.setNumBlocks(D);
    Particle* P = (Particle*)malloc(gpu.nBytes);
    Cell* C     = (Cell*)malloc(gpu.cBytes);   
    reading(path_in,P);                                         
    Particle * d_P = nullptr;
    Cell* d_C      = nullptr;  
    cudaMalloc(&d_C, gpu.cBytes);
    cudaMalloc(&d_P,gpu.nBytes); 
    auto START = std::chrono::system_clock::now();
    while( time <= D.time_end){
        WriteOut(P,out_path,ITERATION,D.part_out_freq,D.part_out_name_base,D.particleCount);
        WriteVTK(P,out_path,ITERATION,D.vtk_out_freq,D.vtk_out_name_base,D.particleCount);
        cudaMemcpy(d_P, P, gpu.nBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_C, C, gpu.cBytes, cudaMemcpyHostToDevice);
        RESET              <<< gpu.num_blocks, gpu.TPB >>> (d_P, D);
        cellInit           <<< gpu.cNumBlocks, gpu.TPB >>> (d_C, gpu);
        MoveAllParticles   <<< gpu.num_blocks, gpu.TPB >>> (d_C, d_P, D);
        ForceTotal         <<< gpu.num_blocks, gpu.TPB >>> (d_C, d_P, D, D.x_n, D.y_n, D.z_n);
        wallForce          <<< gpu.num_blocks, gpu.TPB >>> (d_P, D);
        UPDATE             <<< gpu.num_blocks, gpu.TPB >>> (d_P, D);
        boundarycondition  <<< gpu.num_blocks, gpu.TPB >>> (d_P, D);
        cudaMemcpy(P, d_P, gpu.nBytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(C, d_C, gpu.cBytes, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        time += D.timestep_length;
        ++ITERATION;
    }
    auto END   = std::chrono::system_clock::now();
    std::chrono::duration<double>elepsed_time = END-START;
    printf("RUN_TIME = %f sec. \n",elepsed_time.count());
    cudaFree(d_P);
    cudaFree(d_C);
    free(P);
    free(C);
    return 0;
}
