#include <stdlib.h>
#include <stdio.h>
#include <math.h>
    
//Only main section of code is include.

void read_input(double *L, int *N, double *t_F, double *t_D, double *K, double *A, double *B);

int main(void) {
    
    //output
    FILE *f = fopen("output.txt", "w");
    
    //inputs
    double L;
    int N;
    double t_F;
    double t_D;
    double K;
    double A;
    double B;

    // Read in from file; 
    read_input(&L, &N, &t_F, &t_D, &K, &A, &B);

    double dx = L/N;
    double dt;
    double ctime = 0;
    
    double *U, *U_next;  //y at current and next timestep
    double *V, *V_next;  //y at current and next timestep
    double *Z; //

    /* Allocate memory according to size of nx */
    U       = (double *) malloc(N * sizeof(double));
    U_next  = (double *) malloc(N * sizeof(double));
    V       = (double *) malloc(N * sizeof(double));
    V_next  = (double *) malloc(N * sizeof(double));
    Z       = (double *) malloc(N * sizeof(double));
    
    if (U==NULL||U_next==NULL||V==NULL||V_next==NULL||Z==NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    int j;
    double x;
    
    const double pi = 4.0 * atan(1.0);
    for (j = 0; j < N; j++){
        x = j*dx;
        double cini = cos(2*pi*x/L);
        double sini = sin(2*pi*x/L);
        U[j] = K + A*cini;
        V[j] = B*sini;
    }
    
    double next_time = 0;
    
    while (ctime < t_F){
        
	//tactic 1 for linearisation of PDE
        for (j = 0; j < N; j++){
            Z[j] = sqrt(pow(U[j],2)+pow(V[j],2));
            dt = 1/((1/pow(dx, 2)) + (2*Z[j] + 4*pow(Z[j],3)))*(dx*dx);
        }
        
        double dt0 = dt;
        int output = 0;
        if (ctime + dt0 > next_time){
            dt0 = next_time - ctime;
            output = 1;
        }
        
        // Set boundary
        double rcnst = dt/(dx*dx);
        
        U_next[0] = U[0] + rcnst*(U[1] - 2*U[0] + U[N-1]) + (dt*U[0])*(2 - (4*pow(U[0], 2)) - (4*pow(V[0], 2)));
        V_next[0] = V[0] + rcnst*(V[1] - 2*V[0] + V[N-1]) + (dt*V[0])*(2 - (4*pow(U[0], 2)) - (4*pow(V[0], 2)));
    
        for(j=1; j<N-1; j++){
            double derivU = (U[j+1] + U[j-1] - 2*U[j])/(dx*dx) + (U[j])*(2 - (4*pow(U[j], 2)) - (4*pow(V[j], 2)));
            double derivV = (V[j+1] + V[j-1] - 2*V[j])/(dx*dx) + (V[j])*(2 - (4*pow(U[j], 2)) - (4*pow(V[j], 2)));
            U_next[j] = U[j] + derivU * dt;
            V_next[j] = V[j] + derivV * dt;
        }
        
        // Set boundary end (N-1 is last element of array, as begins at 0)
        U_next[N-1] = U[N-1] + rcnst*(U[1] - 2*U[N-1] + U[N-2]) + (dt*U[N-1])*(2 - (4*pow(U[N - 1], 2)) - (4*pow(V[N - 1], 2)));
        V_next[N-1] = V[N-1] + rcnst*(V[1] - 2*V[N-1] + V[N-2]) + (dt*V[N-1])*(2 - (4*pow(U[N - 1], 2)) - (4*pow(V[N - 1], 2)));
        
        double *temp;
        
        temp = V_next;
        V_next = V; 
        V = temp;
        temp = U_next;
        U_next = U; 
        U = temp;
        
        // increment time
        ctime += dt0;
        if (output) {
          for (j=0; j<N; j++) {
              x = j*dx;
              fprintf(f, "%g %g %g %g\n",ctime,x,U[j],V[j]);
          }
          next_time += t_D;
        }
        
    }
    
    fclose(f);
    free(Z);
    free(U);
    free(U_next);
    free(V);
    free(V_next);
    return 0;
}

