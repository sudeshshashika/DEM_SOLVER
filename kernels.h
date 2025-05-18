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
#include "state.h"
#include "device.h"

CUDA_CALLABLE_MEMBER void InteractionForce(Particle& A, Particle& B, double k_s, double k_dn){
    double x_ij[3]{0.0,0.0,0.0};         //distance vector
    double vn[3]{0.0,0.0,0.0};   //velocity difference
    double n[3]{0.0,0.0,0.0};    //unit vector
    double X_ij= 0.0;            // scaler distance
    for(int i = 0; i < 3; ++i){
        x_ij[i]  = A.d[i] - B.d[i];
        vn[i]    = A.v[i] - B.v[i];
        X_ij     = X_ij + (x_ij[i] * x_ij[i]);

    }
    double rab  = sqrt(X_ij);  //distance between particles
    for(int i = 0; i < 3; ++i){
        n[i] = x_ij[i]/rab;//B.d[i] / rab;    //normalizing the distance
    }    
    double p = A.r + B.r - rab;
    if(p >= 0.0){
        for(int i = 0; i < 3; ++i){
            A.f_s[i] = A.f_s[i] + k_s*p*n[i] - k_dn*vn[i];                          //force between two particles
            B.f_s[i] = B.f_s[i] - k_s*p*n[i] - k_dn*vn[i];                          // force on particle B is the opposite of A
        }
    } 
}
//Total Force Particle parallel cell list
__global__ void ForceTotal(Cell* c,Particle* p,state D,int x_n, int y_n, int z_n){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){

        if(D.reflect_x ==0 && D.reflect_y ==0 && D.reflect_z == 0){
            int X = (p[idx].d[0] -D.x_min )/ D.len_x;
            int Y = (p[idx].d[1] -D.y_min )/ D.len_y;
            int Z = (p[idx].d[2] -D.z_min) / D.len_z;
            for(int i = X-1; i <= X+1; ++i){
                for(int j = Y-1; j <= Y+1; ++j){
                    for(int k = Z-1; k <= Z+1; ++k){
                            int  ncell_ID = ((i+x_n)%x_n)*(y_n*z_n) + ((j+y_n)%y_n)*z_n + ((k+z_n)%z_n);  //calculation of neighbor cell index  
                            int HeadneighborParticle = c[ncell_ID].x;
                            while(HeadneighborParticle != -1){
                                if(idx != HeadneighborParticle){
                                    InteractionForce(p[idx],p[HeadneighborParticle],D.k_s,D.k_dn);
                                } 
                                HeadneighborParticle = p[HeadneighborParticle].id;
                            }
                    }
                }
            }     
        }else{        
            int X = (p[idx].d[0] -D.x_min )/ D.len_x;
            int Y = (p[idx].d[1] -D.y_min )/ D.len_y;
            int Z = (p[idx].d[2] -D.z_min)/  D.len_z;
            for(int i = X-1; i <= X+1; ++i){
                for(int j = Y-1; j <= Y+1; ++j){
                    for(int k = Z-1; k <= Z+1; ++k){
                        if(i>=0 && i<x_n && j>=0 && j<y_n && k>=0 && k<z_n){       //continue;}    //avoid being beyound boundary cells
                            int  ncell_ID = k + (z_n * j ) + (z_n * y_n ) * i;     //calculation of neighbor cell index  
                            int HeadneighborParticle = c[ncell_ID].x;
                            while(HeadneighborParticle != -1){
                                if(idx != HeadneighborParticle){
                                    InteractionForce(p[idx],p[HeadneighborParticle],D.k_s,D.k_dn);
                                } 
                                HeadneighborParticle = p[HeadneighborParticle].id;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void cellInit(Cell* c, device D){
    long long idx = blockIdx.x * blockDim.x  + threadIdx.x;
    if(idx < D.nCells){
        c[idx].x = -1.0;
    }
}
//moveParticle
__global__ void MoveAllParticles(Cell*c, Particle* P,state D){
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;   
    if(idx < D.particleCount){
        P[idx].id = idx;
        int x = (P[idx].d[0] -D.x_min )/ D.len_x;
        int y = (P[idx].d[1] -D.y_min )/ D.len_y;
        int z = (P[idx].d[2] -D.z_min )/  D.len_z;
        int index = x * (D.y_n * D.z_n) + y * D.z_n + z;
        int val = atomicExch(&c[index].x , idx); 
        P[idx].id = val;
    }  
}

CUDA_CALLABLE_MEMBER void wallInteraction(Particle& A, double x_limi,double y_limi,double z_limi,double k_s, double k_dn){
    if(A.m == 1.0){
        double x_ij[3]{0.0,0.0,0.0}; 			//distance vector
        double disb[3]{x_limi,y_limi,z_limi}; 	// wall coord
        double velb[3]{0.0,0.0,0.0}; 		    // wall velocity	
        double vn[3]{0.0,0.0,0.0};  			//velocity difference
        double n[3]{0.0,0.0,0.0};    			//unit vector
        double X_ij= 0.0;            			// scaler distance
        for(int i = 0; i < 3; ++i){
            x_ij[i]  = A.d[i] - disb[i];
            vn[i]    = A.v[i] - velb[i];
            X_ij     = X_ij + (x_ij[i] * x_ij[i]);

        }
        double rab  = sqrt(X_ij); 				 //distance between particles
        for(int i = 0; i < 3; ++i){
            n[i] = x_ij[i]/rab;				    //normalizing the distance
        }    
        double p = A.r + 0.0 - rab;
        if(p >= 0.0){
            for(int i = 0; i < 3; ++i){
                A.f_d[i] = A.f_d[i] + 5.0*k_s*p*n[i]- k_dn* vn[i]; //force from the wall
            }
        }
    } 
}
//boundary force from all particles
__global__ void wallForce(Particle* P,state D){
    
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        if(D.reflect_x ==1){
            if(P[idx].d[0] < D.x_min + P[idx].r + 1.0e-15){
                wallInteraction(P[idx], D.x_min, P[idx].d[1] ,P[idx].d[2],D.k_s, D.k_dn);

            }
            if(P[idx].d[0] > D.x_max- + P[idx].r + 1.0e-15){
                wallInteraction(P[idx],D.x_max, P[idx].d[1],P[idx].d[2],D.k_s, D.k_dn);
            }
        }
        if(D.reflect_y == 1){
            if(P[idx].d[1] < D.y_min + + P[idx].r + 1.0e-15){
                wallInteraction(P[idx],P[idx].d[0], D.y_min ,P[idx].d[2],D.k_s, D.k_dn);
        
            }
            if(P[idx].d[1] > D.y_max - P[idx].r + 1.0e-15){
                wallInteraction(P[idx],P[idx].d[0], D.y_max ,P[idx].d[2],D.k_s, D.k_dn);
            }  
        }
        if(D.reflect_z ==1){
            if(P[idx].d[2] < D.z_min+ P[idx].r + 1.0e-15){
                wallInteraction(P[idx],P[idx].d[0], P[idx].d[1], D.z_min,D.k_s, D.k_dn);
            }
            if(P[idx].d[2] > D.z_max- P[idx].r + 1.0e-15){
                wallInteraction(P[idx],P[idx].d[0], P[idx].d[1] , D.z_max,D.k_s, D.k_dn);
            } 
        }
    }
}

__global__ void UPDATE(Particle *P, state D){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        double gravity[3]{D.g_x,D.g_y,D.g_z};
        // if(P[idx].m == 1.0){ gravity[0]=g_x;gravity[1]=g_y;gravity[2]=g_z;};
        if(P[idx].m > 1.0) { gravity[0]=0;gravity[1]=0;gravity[2]=0;} 
        for(int dim = 0; dim< 3; ++dim){
            P[idx].v[dim] = P[idx].v[dim] + D.timestep_length  * (gravity[dim] + ((P[idx].f_s[dim] + P[idx].f_d[dim]) / P[idx].m)); 
            if(P[idx].m == 1.0){
                P[idx].d[dim] = P[idx].d[dim] + (D.timestep_length * P[idx].v[dim]);
            } 
        } 
    }  
}

__global__ void boundarycondition(Particle*  P, state D){
                                                                                                                  
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        if(D.reflect_x == 0){
            if(P[idx].d[0] < D.x_min + P[idx].r){
                P[idx].d[0] = D.x_max - P[idx].r;
            }
            if(P[idx].d[0] > D.x_max - P[idx].r){
                P[idx].d[0] = D.x_min + P[idx].r;
            }
        }
        if(D.reflect_y == 0){
            if(P[idx].d[1] < D.y_min + P[idx].r){
                P[idx].d[1] = D.y_max - P[idx].r;
            }
            if(P[idx].d[1] > D.y_max - P[idx].r){
                P[idx].d[1] = D.y_min + P[idx].r;
            }
        }
        if(D.reflect_z == 0){
            if(P[idx].d[2] < D.z_min + P[idx].r){
                P[idx].d[2] = D.z_max - P[idx].r;
            }
            if(P[idx].d[2] > D.z_max - P[idx].r){
                P[idx].d[2] = D.z_min + P[idx].r;
            }
        }
    }

}

__global__ void RESET(Particle* p, state D){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        for(int dim=0;dim<3;++dim){
            p[idx].f_s[dim] = 0.0;
            p[idx].f_d[dim] = 0.0;
        }
    }
}

