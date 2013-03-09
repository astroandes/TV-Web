#ifndef __FIELD_H__
#define __FIELD_H__

#include "density_struct.h"

#define USE_FFT 1
#define USE_FD  2

int ComputeDensityGradiant(density_grid *,int);
int ComputeDensityHessian(density_grid *,int);
int DensityField2DensityGrid(density_field,density_grid *);
int Compute_CIC_Density(float *,float *,int,int,long long,density_grid *,int,int,int,
			float,float,float,float,float,float);
int GaussSmooth(density_grid *,float, int);
int FilterSmooth(density_grid *grid,double (*filterfunc)(double,double,double,void *),const char *filtername,void *parms);
int ComputePotential(density_grid *grid);
int ComputeVelocityTorsion(density_grid *density);
#endif
