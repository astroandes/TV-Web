#include <stdio.h>
#include <stdlib.h>
#define USAGE "./load_env.x filename"
/*
   USAGE: ./load_cic.x filename
   AUTHOR: Jaime Forero-Romero j.e.forero.romero@gmail.com
   DESCRIPTION:
   Loads a 3Dimensional grid into a 3D array.
   NOTES:
   n_x, n_y, n_z is the grid size
   x_0, y_0, z_0 is the posiion of the cell for the grid(i,j,k) element in kpc/h
   dx,dy,z is the grid size in kpc/h
   
   grid[n] runs over x,y,z in space, following fortran array conventions
   for memory handling
   
   i.e. a 3D point grid(i,j,k) corresponds to  grid[i + N_X * (j + N_Y * k)] where N_X is the grid dimension size in x and N_Y is the grid dimension size in y (i,j,k start at 0).

*/
#define FLOAT float
int main(int argc, char **argv){
  FILE *in;
  FLOAT *grid;
  int dumb;
  char line[30];
  long long i;
  long long n_total;
  int n_x, n_y, n_z;
  int n_nodes;
  float dx, dy, dz, x_0, y_0, z_0;
  FLOAT max_val, min_val;
  float pos_x, pos_y, pos_z;
  int pos_i, pos_j, pos_k;

  if(argc!=5){
    fprintf(stderr, "USAGE: %s\n", USAGE);
    exit(1);
  }
  
  pos_x = atof(argv[2]);
  pos_y = atof(argv[3]);
  pos_z = atof(argv[4]);

  if(!(in=fopen(argv[1], "r"))){
    fprintf(stderr, "Problem opening file %s\n", argv[1]);
    exit(1);
  }
  fread(&dumb,sizeof(int),1,in);
  fread(line,sizeof(char)*30,1,in);
  fread(&dumb,sizeof(int),1,in);

  fread(&dumb,sizeof(int),1,in);
  fread(&n_x,sizeof(int),1,in);    
  fread(&n_y,sizeof(int),1,in);    
  fread(&n_z,sizeof(int),1,in);    
  fread(&n_nodes,sizeof(int),1,in);    
  fread(&x_0,sizeof(float),1,in);    
  fread(&y_0,sizeof(float),1,in);    
  fread(&z_0,sizeof(float),1,in);    
  fread(&dx,sizeof(float),1,in);    
  fread(&dy,sizeof(float),1,in);    
  fread(&dz,sizeof(float),1,in);    
  fread(&dumb,sizeof(int),1,in);

  n_total = n_x * n_y * n_z;
  if(!(grid=malloc(n_total * sizeof(FLOAT)))){
    fprintf(stderr, "problem with array allocation\n");
    exit(1);
  }
  
  fread(&dumb,sizeof(int),1,in);
  fread(&(grid[0]),sizeof(FLOAT), n_nodes, in);
  fread(&dumb,sizeof(int),1,in);  

  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  fclose(in);


  min_val = 1.0E10;
  max_val = -1.0E10;

  pos_i = (int)(pos_x/dx);
  pos_j = (int)(pos_y/dy);
  pos_k = (int)(pos_z/dz);
  
  fprintf(stdout, "%d %d %d\n", pos_i, pos_j, pos_k);
  fprintf(stdout, "%g\n", grid[pos_i + n_x * (pos_j + n_y * pos_k)]);


  return 0;
}
