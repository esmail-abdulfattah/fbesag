
#include <assert.h>
#include <iostream>

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include "cgeneric.h"
#define Calloc(n_, type_) (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x)) 

double *inla_cgeneric_pbesag_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{
    bool debug = false;
    double *ret = NULL, *prec = NULL; //lprec = (theta ? theta[0] : NAN);
    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    int M = data->ints[3]->ints[2];   
    int npart = data->ints[2]->ints[0];

    if(theta){
        prec = Calloc(npart, double);
        for(int i=0; i<npart; i++) prec[i] = exp(theta[i]);
    }

    if(debug){

        std::cout << npart << std::endl;

        for(int i=0; i<10; i++) //VEC_CGENERIC_GRAPH
            std::cout << data->ints[3]->ints[i] << ", ";

        std::cout << std::endl;
        for(int i=0; i<13; i++) //misc
            std::cout << data->ints[4]->ints[i] << ", ";
        std::cout << std::endl; 
    
        std::cout << data->doubles[0]->doubles[0] << " -pc- " << std::endl;
        std::cout << data->doubles[1]->doubles[0] << " -sd- " << std::endl;

    }

    switch (cmd) {

        case INLA_CGENERIC_VOID:
        {
            assert(!(cmd == INLA_CGENERIC_VOID));
            break;
        }

        case INLA_CGENERIC_GRAPH:
        {       
            ret = Calloc(2 + 2*M, double);
            ret[0] = data->ints[3]->ints[1]; //It is N: dimesnion of the graph
            ret[1] = M; //length of ii vector which means number of non-zero
            for (int i = 0; i < M; i++) {
                ret[2 + i] = data->ints[3]->ints[3 + i];
                ret[2 + M + i] = data->ints[3]->ints[3 + M + i];
            }

            /*
            for (int i = 0; i < 2 + 2*M ; i++){
                std::cout << ret[i] << ", ";
            }
            exit(0);
            */
            
            break;
        }

        case INLA_CGENERIC_Q:
        {
            
            double scaled_cnst = data->doubles[1]->doubles[2];
            ret = Calloc(2 + M, double);
            ret[0] = -1;
            ret[1] = M;
            int count = 2, k = 0;
            for (int i = 0; i < N; i++) {

                int num_nei_i = data->ints[4]->ints[k];
                k++;
                double mas_prec = prec[data->ints[4]->ints[k]];
                k++;
                int size_neighbors = data->ints[4]->ints[k];
                k++;

                int save_count = count;
                ret[count] = 0;
                count++;
                double tmp_neighbors = 0.0, sum_tmp_neighbors = 0.0;
                for(int j=0; j<size_neighbors; j++){

                    tmp_neighbors = -0.5*prec[data->ints[4]->ints[k]];

                    k++;
                    int tick = data->ints[4]->ints[k];
                    k++;

                    if(tick>i){
                        ret[count] = scaled_cnst*(tmp_neighbors - 0.5*mas_prec);
                        count++;
                    }
                    
                    sum_tmp_neighbors += tmp_neighbors;
                }

                ret[save_count] = scaled_cnst*(-sum_tmp_neighbors + 0.5*num_nei_i*mas_prec) + 1e-5;
            }      

            /*
            for (int i = 0; i < 2 + M ; i++){
                std::cout << ret[i] << ", ";
            }
            exit(0);
            */
 
            break;
        }

        case INLA_CGENERIC_MU:
        {
            ret = Calloc(npart, double);
            ret[0] = npart;
            for(int i=0; i<npart;i++)
                ret[i] = 0;
            break;

        }
        
        case INLA_CGENERIC_INITIAL:
        {
            ret = Calloc(npart+1, double);
            ret[0] = npart;
            //for(int i=1; i<=npart;i++)
            //    ret[i] = 4.0;

            for(int i=1; i<=npart;i++)
                ret[i] = data->doubles[2]->doubles[i-1]; 

            //exit(0);
            break;
        }

        case INLA_CGENERIC_LOG_NORM_CONST:
        {
            return NULL;
            break;
        }

        case INLA_CGENERIC_LOG_PRIOR:
        {
            double lam = data->doubles[0]->doubles[0], mean_theta=0;
            double val1 = data->doubles[1]->doubles[0], val2 = data->doubles[1]->doubles[1];

            for(int i=0 ;i <npart; i++){
                mean_theta += theta[i];
                //std::cout << theta[i] << std::endl;
            }
                
            mean_theta = mean_theta/npart;
            
            double res1 = log(lam) - (lam)*exp(-0.5*mean_theta) -0.5*mean_theta, res2 = 0;
            for(int i=0; i < npart; i++){
                res2 += val1*(theta[i]-mean_theta)*(theta[i]-mean_theta);
                for(int j=0; j < i; j++){
                    res2 += 2*val2*(theta[i]-mean_theta)*(theta[j]-mean_theta);
                }
            }
            
            //std::cout << res1 << std::endl;
            //std::cout << -0.5*res2 << std::endl;

            ret = Calloc(1, double);
            ret[0] = res1 -0.5*res2;
            break;
        }

        case INLA_CGENERIC_QUIT:
        default:
            break;
    }

    return (ret);
}