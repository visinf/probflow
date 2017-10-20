#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>

#include "mex.h"

// Function based on EpicFlow implementation of Revaud et al. [36]
void sor_coupled_slow_but_readable(float *du, float *dv, const float *a11, const float *a12, const float *a22, const float *b1, const float *b2, const float *dpsis_horiz, const float *dpsis_vert, const int iterations, const float omega, int width, int height){
    int i,j,iter;

    float sigma_u,sigma_v,sum_dpsis,A11,A22,A12,B1,B2,det;
    for(iter=0; iter<iterations; iter++){
        for(j=0; j<height; j++){
	        for(i=0; i<width; i++){
	            sigma_u = 0.0f;
	            sigma_v = 0.0f;
	            sum_dpsis = 0.0f;
	            if(j>0){
		            sigma_u -= dpsis_vert[(j-1)*width+i] * du[(j-1)*width+i];
		            sigma_v -= dpsis_vert[(j-1)*width+i] * dv[(j-1)*width+i];
		            sum_dpsis += dpsis_vert[(j-1)*width+i];
		        }
	            if(i>0){
                    sigma_u -= dpsis_horiz[j*width+i-1] * du[j*width+i-1];
                    sigma_v -= dpsis_horiz[j*width+i-1] * dv[j*width+i-1];
                    sum_dpsis += dpsis_horiz[j*width+i-1];
		        }
	            if(j<(height-1)){
		            sigma_u -= dpsis_vert[j*width+i] * du[(j+1)*width+i];
		            sigma_v -= dpsis_vert[j*width+i] * dv[(j+1)*width+i];
		            sum_dpsis += dpsis_vert[j*width+i];
		        }
	            if(i<(width-1)){
		            sigma_u -= dpsis_horiz[j*width+i] * du[j*width+i+1];
		            sigma_v -= dpsis_horiz[j*width+i] * dv[j*width+i+1];
		            sum_dpsis += dpsis_horiz[j*width+i];
		        }
                A11 = a11[j*width+i] + sum_dpsis;
                A12 = a12[j*width+i];
                A22 = a22[j*width+i] + sum_dpsis;
                det = A11*A22 - A12*A12;
                B1 = b1[j*width+i] - sigma_u;
                B2 = b2[j*width+i] - sigma_v;
                du[j*width+i] = (1.0f-omega) * du[j*width+i] + omega * ( A22*B1-A12*B2)/det;
                dv[j*width+i] = (1.0f-omega) * dv[j*width+i] + omega * (-A12*B1+A11*B2)/det;
	        }
	    }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs != 11 || nlhs != 2)
    mexErrMsgTxt("wrong number of parameters");


    int m = mxGetM(prhs[0]);
    int n = mxGetN(prhs[0]);
    size_t s = n * m * sizeof(float);

    float *du= malloc(s);
    memcpy(du, mxGetData(prhs[0]), s);

    float *dv= malloc(s);
    memcpy(dv, mxGetData(prhs[1]), s);

    int niter= (int) mxGetScalar(prhs[9]);
    float omega= (float) mxGetScalar(prhs[10]);

    sor_coupled_slow_but_readable(du, dv, mxGetData(prhs[2]), mxGetData(prhs[3]), mxGetData(prhs[4]), mxGetData(prhs[5]), mxGetData(prhs[6]), mxGetData(prhs[7]), mxGetData(prhs[8]), niter, omega, m, n);

    plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
    float * duOut = (float *) mxGetData(plhs[0]);
    memcpy(duOut,du,s);
    plhs[1] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
    float * dvOut = (float *) mxGetData(plhs[1]);
    memcpy(dvOut,dv,s);

    free(du);
    free(dv);
}
