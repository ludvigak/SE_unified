#include "../SE_fgg.h"
#include <stdlib.h>
#include <stdio.h>

// gcc -Wall -g -O3 -DTHREE_PERIODIC -DFGG_SPLIT -std=c99 fgg_test.c ../SE_fgg.c -lm

#ifdef FGG_SPLIT
#define PRECOMP_FGG_EXPA 1
#else
#define PRECOMP_FGG_EXPA 0
#endif

#define FGG_PREORDER 1

int main(int argc, char* argv[])
{
    const int N = atoi(argv[1]);
    const int M = atoi(argv[2]);
    const int P = atoi(argv[3]);

    const double c = 2.31;
    const double h = 1.0/(double)M;

    /* parameter struct */
    SE_FGG_params params;
    SE_FGG_pack_params(&params, N, M, M, M, P, c, h);

    /* system of point charges */
    SE_state s;
    SE_init_system(&s, &params);

    /* storage for termporaries */
    SE_FGG_work work;
    SE_FGG_allocate_workspace(&work, &params, 1, 1);

    /* holds final result */
    double* H_per = SE_FGG_allocate_grid(&params); 

    /* done with initialization. now compute */

    /* possibly permute s in grid-index order */ 
    if(FGG_PREORDER)
    {
	SE_FGG_expand_all(&work, &s, &params);
	SE_FGG_reorder_system(&s, &work, &params);
    }

    /* compute Gaussian on P^3 grid */
    double t0 = SE_gettime();
    SE_FGG_base_gaussian(&work, &params);
    if(PRECOMP_FGG_EXPA)
    {
	/* pre-compute FGG expansion vectors, then grid */
	double t1 = SE_gettime();
	SE_FGG_expand_all(&work, &s, &params);
	t1 = SE_gettime()-t1;

	double t2 = SE_gettime();
	SE_FGG_grid_split_SSE_dispatch(&work, &s, &params);	
	t2 = SE_gettime()-t2;

	printf("[FG G split] %f\t%f\t flops: %e\n",t1,t2,(N/t2)*pow(P,3)*3 );
    }
    else
    {
	/* unified expansion-gridding */ 
	SE_FGG_grid(&work, &s, &params);
    }
    /* have computed grid function in non-periodic domain. wrap it to
       produce periodicity */
    SE_FGG_wrap_fcn(H_per,&work,&params);

    printf("[FGG GRID] %d %d %d \t%f\n", N, M, P, SE_gettime()-t0);

    // now the integration/interpolation part
    
    // assume zs (static gaussian already computed)
    double* phi = SE_FGG_allocate_vec(N);

    t0=SE_gettime();

    // extend periodic function H_per non-periodic grid
    SE_FGG_extend_fcn(&work, H_per, &params);

    if(PRECOMP_FGG_EXPA)
    {
	double t1 = SE_gettime();
	SE_FGG_int_split_SSE_dispatch(phi, &work, &params);	
	t1 = SE_gettime()-t1;
	printf("[FG I split] %f\t\t\t flops: %e\n",t1,(N/t1)*pow(P,3)*4 );
    }
    else
    {
	SE_FGG_int(phi, &work, &s, &params);	
    }
    printf("[FGG INT] %d %d %d \t%f\n", N, M, P, SE_gettime()-t0);

    /* clean up */ 
    SE_FGG_free_workspace(&work);
    SE_free_system(&s);
    free(H_per);
    free(phi);

    return 0;
}
