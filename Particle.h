#pragma once

#include <iostream>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

class Particle{
    /*friend std::ostream &operator<<(std::ostream& os, const Particle &obj){
        os << obj.m << "     " << obj.d[0] << " " << obj.d[1] << " " << obj.d[2] 
                    << "     " << obj.v[0] << " " << obj.v[1] << " " << obj.v[2];
        return os;
    }*/
public:
    double m;         //particle mass
    double r;         //radius
    double d[3];      //position vector
    double v[3];      // velocity vector
    double f_s[3];    //spring force
    double f_d[3];    // damping force
    int    id;
    CUDA_CALLABLE_MEMBER Particle(double m_=1.0,double r = 0, int x_ = 0): m{m_},r{r}, id{x_}{
        for(int i = 0; i < 3; ++i){
            d[i] = 0.0;
            v[i] = 0.0;
            f_s[i] = 0.0;
            f_d[i] = 0.0;
        }
    }
    CUDA_CALLABLE_MEMBER Particle(const Particle& src): m{src.m},r{src.r}, id{src.id}{
        for(int i =0; i < 3; ++i){
            d[i] = src.d[i];
            v[i] = src.v[i];
            f_s[i] = src.f_s[i];
            f_d[i] = src.f_d[i];
        }
    }
    CUDA_CALLABLE_MEMBER Particle& operator=(const Particle& src){
        if(this == &src)
            return *this;
        m  = src.m;
        r  = src.r;
        id = src.id;
        for(int i =0; i < 3; ++i){
            d[i] = src.d[i];
            v[i] = src.v[i];
            f_s[i] = src.f_s[i];
            f_d[i] = src.f_d[i];
        }
        return *this;           
    }
    CUDA_CALLABLE_MEMBER ~Particle(){}
};
