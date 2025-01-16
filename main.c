/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"
#include <omp.h>

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;
    double 	***f = NULL;

    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }


    // allocate memory
    if ( (u = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( (f = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
    }

    #ifdef _JACOBI
    double 	***u_2 = NULL;
    if ( (u_2 = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u_2: allocation failed");
        exit(-1);
    }
    #endif

    double start = omp_get_wtime();

    // --------Initialization-----------
    init_u(u, N+2, start_T);
    init_f(f, N+2);
    #ifdef _JACOBI
    init_u(u_2, N+2, start_T);
    #endif

    #pragma omp parallel default(none) shared(f, u, u_2, N, iter_max, tolerance, start_T)
    {

    // --------Calculation-----------
    #ifdef _JACOBI
    jacobi(f, u, u_2, N+2, iter_max, tolerance);
    #endif

    #ifdef _GAUSS_SEIDEL
    gauss_seidel(f, u, N+2, iter_max, tolerance);
    #endif

    }

    double elapsed_time = omp_get_wtime() - start;

    // Calculate lups
    double Mlups = (double)N * N * N * iter_max / elapsed_time / 1e6;

    // Calculate memory footprint
    #ifdef _JACOBI
    int mem_footprint = N * N * N * 3 * 8;
    #endif
    #ifdef _GAUSS_SEIDEL
    int mem_footprint = N * N * N * 2 * 8;
    #endif

    printf("%d\t%d\t%f\t%f\n", omp_get_max_threads(), mem_footprint, elapsed_time, Mlups);

    // dump  results if wanted 
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N+2, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N+2, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free_3d(u);
    #ifdef _JACOBI
    free_3d(u_2);
    #endif
    free_3d(f);

    return(0);
}
