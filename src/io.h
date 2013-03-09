#ifndef _H_IO
#define _H_IO

#include "struct.h"
#include "endian.h"

#define FLAG_POS  ((int)1<<0)
#define FLAG_VEL  ((int)1<<1)
#define FLAG_MASS ((int)1<<2)
#define FLAG_ID   ((int)1<<3)
#define FLAG_TYPE ((int)1<<4)
#define FLAG_GAS  ((int)1<<5)
#define FLAG_ALL  (((int)1<<6)-1)
#define FLAG_DM   ((int)FLAG_POS|FLAG_VEL|FLAG_MASS|FLAG_ID)
//if the file is not of the same endian type as the computer
//this is done automatically
#define FLAG_SWAPENDIAN ((int)1<<31)

#define SKIP(filename) {int skpdummy;fread(&skpdummy,sizeof(skpdummy),1,filename);}

int SaveEigenvectorGrid(char *fname, density_grid *density, int n_eigen);
int LoadEnvGrid(char *fname,density_grid *density);
ssize_t Mygetline (char **lineptr, int *n, FILE *stream);
int IsGadgetFile(char *);
int ReadGadget(char*,snapshot_data*,int);
int ReadSIMPLE2Gadget(char*,snapshot_data*,int);
int LoadSurveyCone(SurveyCone *,const char *);
int LoadTree(const char *,Grid **);
int LoadDensity(char *,density_field *);
int LoadDensityGrid(char *,density_grid *);
int LoadEigenvalueGrid(char *,density_grid *);
int SaveDensityGrid(char *,density_grid *);
int SaveEigenvalueGrid(char *,density_grid *, int);
int SaveTraceGrid(char *,density_grid *);
int SaveEnvGrid(char *fname,density_grid *density);
int WriteGadget(char *fname,snapshot_data *snap,int flags);
int ReadART(char *, char *, snapshot_data *, int, int);
int LoadVelocityGrid(char *fname, density_grid *density, int comp);
#endif
