#include "Solver.h"

using namespace std ; 

vector<double>  solve3Diag(const vector <double> & lDiag, const vector <double> & Diag, const vector <double> & uDiag,
	                               const vector <double> & rHS) {

   // --- Initialize cuSPARSE
    cusparseHandle_t handle;    cusparseCreate(&handle);

    const int N     = 5;        // --- Size of the linear system

    // --- Lower diagonal, diagonal and upper diagonal of the system matrix
    double *h_ld = (double*)malloc(N * sizeof(double));
    double *h_d  = (double*)malloc(N * sizeof(double));
    double *h_ud = (double*)malloc(N * sizeof(double));

    h_ld[0]     = 0.;
    h_ud[N-1]   = 0.;
    for (int k = 0; k < N - 1; k++) {
        h_ld[k + 1] = -1.;
        h_ud[k]     = -1.;
    }
    for (int k = 0; k < N; k++) h_d[k] = 2.;

    double *d_ld;   cudaMalloc(&d_ld, N * sizeof(double));
    double *d_d;    cudaMalloc(&d_d,  N * sizeof(double));
    double *d_ud;   cudaMalloc(&d_ud, N * sizeof(double));

    cudaMemcpy(d_ld, h_ld, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_d,  h_d,  N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ud, h_ud, N * sizeof(double), cudaMemcpyHostToDevice);

    // --- Allocating and defining dense host and device data vectors
    double *h_x = (double *)malloc(N * sizeof(double)); 
    h_x[0] = 100.0;  h_x[1] = 200.0; h_x[2] = 400.0; h_x[3] = 500.0; h_x[4] = 300.0;

    double *d_x;       cudaMalloc(&d_x, N * sizeof(double));   
    cudaMemcpy(d_x, h_x, N * sizeof(double), cudaMemcpyHostToDevice);

    // --- Allocating the host and device side result vector
    double *h_y = (double *)malloc(N * sizeof(double)); 
    double *d_y;        cudaMalloc(&d_y, N * sizeof(double)); 

    cusparseDgtsv(handle, N, 1, d_ld, d_d, d_ud, d_x, N);

    cudaMemcpy(h_x, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
    for (int k=0; k<N; k++) printf("%f\n", h_x[k]);
	vector < double> ans ; 
	for (int k=0; k<N; k++) {
	   ans.push_back(h_x[k]); 

	}
	return ans ; 

}
 
