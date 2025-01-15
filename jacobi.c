/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>

void
jacobi(double ***f, double ***u, double ***u_2, int N, int iter_max, double tolerance) {
    double ***temp;
    double prev;
    double diff = INFINITY;
    double div = 1.0/(N*N*N);
    double delta_2 = 4.0 / (N*N);
    int iter = 0;
    while(iter < iter_max && diff*div > tolerance) {
        diff = 0.0;
        iter++;
        for(int i = 1; i < N-1; i++) {
            for(int j = 1; j < N-1; j++) {
                for(int k = 1; k < N-1; k++) {
                    u_2[i][j][k] = (u[i-1][j][k] + u[i+1][j][k] +
                                    u[i][j-1][k] + u[i][j+1][k] +
                                    u[i][j][k-1] + u[i][j][k+1] +
                                    delta_2 * f[i][j][k]) / 6.0;
                    diff += fabs(u_2[i][j][k] - u[i][j][k]);
                }
            }
        }
        // Swap u and u_2
        temp = u;
        u = u_2;
        u_2 = temp;
        printf("diff: %f\n", diff*div);
    }
}
