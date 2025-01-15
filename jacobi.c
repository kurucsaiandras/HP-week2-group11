/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>

void
jacobi2(double ***f, double ***u, double ***u_2, int N, int iter_max, double tolerance) {
    double ***temp;
    double prev;
    double diff = INFINITY;
    double div = 1.0/(N*N*N);
    double delta_2 = 4.0 / (N*N);
    int iter = 0;
    while(iter < iter_max) {
        diff = 0.0;
        iter++;
        #pragma omp parallel for default(none) shared(f, u, u_2, delta_2, N)
        //#pragma omp parallel for
        for(int i = 1; i < N-1; i++) {
            for(int j = 1; j < N-1; j++) {
                for(int k = 1; k < N-1; k++) {
                    u_2[i][j][k] = (u[i-1][j][k] + u[i+1][j][k] +
                                    u[i][j-1][k] + u[i][j+1][k] +
                                    u[i][j][k-1] + u[i][j][k+1] +
                                    delta_2 * f[i][j][k]) / 6.0;

                }
            }
        }

        // Swap u and u_2
        temp = u;
        u = u_2;
        u_2 = temp;
    }
}



void jacobi(double ***f, double ***u, double ***u_2, int N, int iter_max, double tolerance) {
    double ***temp;
    double delta_2 = 4.0 / (N*N);
    int iter = 0;
    while(iter < iter_max) {
        iter++;

        #pragma omp parallel for schedule(static, 16) default(none) shared(f, u, u_2, delta_2, N)
        for(int i = 1; i < N-1; i++) {
            for(int j = 1; j < N-1; j++) {
                for(int k = 1; k < N-1; k++) {
                    u_2[i][j][k] = (u[i-1][j][k] + u[i+1][j][k] +
                                    u[i][j-1][k] + u[i][j+1][k] +
                                    u[i][j][k-1] + u[i][j][k+1] +
                                    delta_2 * f[i][j][k]) / 6.0;

                }
            }
        }

        // Swap u and u_2
        temp = u;
        u = u_2;
        u_2 = temp;
    }
}