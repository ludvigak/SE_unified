#ifndef ROTLET_DIRECT_H
#define ROTLET_DIRECT_H

#define PI 3.141592653589793

typedef struct 
{
    double box[3];
    double xi;
    double rc;
} ewald_opts;


void rotlet_direct_rsrc(double*, const double*, int,
				const double*, const double*, int, 
				const ewald_opts);

#endif
