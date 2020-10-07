#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lapacke.h>

// BAND MATRIX CREATION AND USE IS NOT INCLUDED

void read_input(double *L, int *N, double *K, double *gamma, double *freqA, double *freqB, int *I_0);
void read_coefficients(double *mu, double *E);

int reindex(int N, int i) {
    int j = 0;
    if (i < N/2) {
        j = 2*i;
    } else {
        j = 2 * (N - i) - 1;
    }
    return j;
}

int main () {
    
    band_mat bmat;
    
    FILE *f = fopen("output.txt","w");
    
    double L;
    int N;
    double K;
    double gamma;
    double freqA;
    double freqB;
    int I_0;
  
    read_input(&L, &N, &K, &gamma, &freqA, &freqB, &I_0);

    long ncols;
    if (N < 3){
        ncols = 4 * N;
    }
    else{
        ncols = 2 * N;
    }
    long nbands_low = 4;		
    long nbands_up = 4;	
	double dx = L/N;
	double dx_sq = dx * dx;

    init_band_mat (&bmat, nbands_low, nbands_up, ncols);
    
    double *x = malloc(sizeof (double) * ncols);
    double *b = malloc(sizeof (double) * ncols);
    
    double *mu = malloc(sizeof (double) * N);
    double *E = malloc(sizeof (double) * (N+1));
    
    read_coefficients(mu, E);
	E[N] = E[0];
    
    for (int q = 0; q <= (I_0-1); q++) {
        
        double freqQ = freqA + q*(freqB - freqA)/I_0;

	// SET BOUNDARIES
	// Matrix is solved in a specific pattern, each edge being resolved after the other.
        
        for (int i = 0; i < ncols; i++) {
            double a_0 = (1/dx_sq)*((freqQ * freqQ * mu[i/2] * dx_sq) - E[i/2 + 1] - E[i/2]);
            double b_0 = (freqQ * gamma);
            double a_1 = (1/dx_sq)*E[i/2];
            double a_end = (1/dx_sq) * E[0] * cos (K * L);
            double b_end = (1/dx_sq) * E[0] * sin (K * L);
			
            if (i % 2 == 0) {
                b[reindex(ncols,i)] = -1.0*cos(K*i*dx/2);
                if (i > 0) {
                    setv(&bmat, reindex(ncols, i - 1), reindex(ncols, i), 0);
                    setv(&bmat, reindex(ncols, i), reindex(ncols, i - 1), 0);
                }
            }
            else{
                b[reindex(ncols,i)] = -1.0*sin(K*(i-1)*dx/2);
                setv(&bmat, reindex(ncols, i - 1), reindex(ncols, i), -1.0*b_0);
                setv(&bmat, reindex(ncols, i), reindex(ncols, i - 1), b_0);
            }

            if (i > 1) {
                setv(&bmat, reindex(ncols, i-2), reindex(ncols, i), a_1);
                setv(&bmat, reindex(ncols, i), reindex(ncols, i-2), a_1);
            }
            
            setv(&bmat, reindex(ncols, i), reindex(ncols, i), a_0);
            setv(&bmat, reindex(ncols, 0), reindex(ncols, ncols - 2), a_end);
            setv(&bmat, reindex(ncols, 0), reindex(ncols, ncols - 1), b_end);
            setv(&bmat, reindex(ncols, 1), reindex(ncols, ncols - 1), a_end);
            setv(&bmat, reindex(ncols, 1), reindex(ncols, ncols - 2), -1.0*b_end);
            setv(&bmat, reindex(ncols, ncols - 2), reindex(ncols, 0), a_end);
            setv(&bmat, reindex(ncols, ncols - 1), reindex(ncols, 0), b_end);
            setv(&bmat, reindex(ncols, ncols - 1), reindex(ncols, 1), a_end);
            setv(&bmat, reindex(ncols, ncols - 2), reindex(ncols, 1), -1.0*b_end);
        }
        solve_Ax_eq_b(&bmat, x, b);
        
        //Check for singularity. Function works when info = 0, else there is a singular matrix.

        if (solve_Ax_eq_b(&bmat, x, b) != 0){
            fprintf(f, "%g \n", NAN);
        }
        else{
            for (int i = 0; i < N; i++) {
                fprintf(f,"%g %g %g %g \n", freqQ, i*dx, x[reindex(ncols,2*i)], x[reindex(ncols,2*i+1)]);
            }
        }
    }
    free(mu);
    free(E);
    free(x);
    free(b);
    fclose(f);
}

void read_input(double *L, int *N, double *K, double *gamma, double *freqA, double *freqB, int *I_0) {
    FILE *input;
    
    if (!(input = fopen ("input.txt", "r"))){
      printf ("Error opening file\n");
      exit (1);
    }
    
    if (7 != fscanf (input, "%lf %d %lf %lf %lf %lf %d", L, N, K, gamma, freqA, freqB, I_0)){
      printf ("Error reading parameters from file\n");
      exit (1);
    }
    
    fclose (input);
}

void read_coefficients(double *mu, double *E){
	FILE *coefficients;
	
	int m = 0;
	double n; double o;
    
    if((coefficients=fopen("coefficients.txt","r")) == NULL) {
        printf("Error opening file\n");
        exit(1);
    }
	
	while(2 == fscanf(coefficients, "%lf %lf", &n, &o)){
		mu[m] = n; E[m] = o; 
		m++;
	}
	
    fclose(coefficients);
}