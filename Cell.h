/**This is for the cell list structure*/
#pragma once
#include <iostream>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

class Cell{

public: 
    int x;                                                       //cell id
    CUDA_CALLABLE_MEMBER Cell( int x_ = -1 ): x{x_}{}            // constructor
    CUDA_CALLABLE_MEMBER Cell( const Cell& src): x{src.x}{}      // copy constructor
    CUDA_CALLABLE_MEMBER Cell &operator=( const Cell& src){     // copy assignment
        if( this == &src)
            return *this;
        x = src.x;
        return *this;
    }
    
    CUDA_CALLABLE_MEMBER ~Cell(){}                             //destructor
};
