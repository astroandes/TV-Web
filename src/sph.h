#ifndef SPH_H_
#define SPH_H_

#include "nei_find.h"

#ifndef PI
#define PI 3.1415926535897932384
#endif

double SPH_Kernel(double,double);
double SPH_DiffKernel(double[],double,int,double);
double SPH_ComputeDensity(Nei_List *,int,double (*)(double,double),double,int);
double SPH_Estimateh(Nei_List *,int,double (*)(double,double),double *,int);

#endif
