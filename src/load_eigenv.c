#include <stdio.h>
#include <stdlib.h>
#define USAGE "./load_eigenvec.x filename vecID"
/*
   USAGE: ./load_eigenvec.x filename
   AUTHOR: Jaime Forero-Romero j.e.forero.romero@gmail.com
   DESCRIPTION:
   Loads a 3Dimensional grid of eigenvectors into a 3D array.
   NOTES:
   n_x, n_y, n_z is the grid size
   x_0, y_0, z_0 is the posiion of the cell for the grid(i,j,k) element in kpc/h
   dx,dy,z is the grid size in kpc/h
   
   grid[n] runs over w,x,y,z following fortran array conventions
   for memory handling.
   w -> dimension of eigen vectors (runs over 0,1,2)
   x,y,z -> spatial dimension over the grid
   
   i.e. a 4D point grid(l,i,j,k) corresponds to component l of the vector at 
   grid[l + 3 * (i + N_X * (j + N_Y * k))] where N_X is the grid dimension size in x 
   and N_Y is the grid dimension size in y (i,j,k start at 0).

*/
#define FLOAT float
int main(int argc, char **argv){
  FILE *in;
  FLOAT *grid;
  int dumb;
  char line[30];
  long long i, j,k, n_total, vecid;  
  int n_x, n_y, n_z;
  int n_nodes;
  float dx, dy, dz, x_0, y_0, z_0;
  FLOAT max_val, min_val;
  long long min_i,min_j,min_k, ind;
  

  if(argc!=3){
    fprintf(stderr, "USAGE: %s\n", USAGE);
    exit(1);
  }
  

  if(!(in=fopen(argv[1], "r"))){
    fprintf(stderr, "Problem opening file %s\n", argv[1]);
    exit(1);
  }
  
  vecid = atoll(argv[2]);

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
  n_total = n_x *n_y * n_z;

  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  fprintf(stderr, "Ntotal %i\n", n_total);

  if(!(grid=malloc(3 * n_total * sizeof(float)))){
    fprintf(stderr, "problem with array allocation\n");
    exit(1);
  }
  
  fread(&dumb,sizeof(int),1,in);
  fread(&(grid[0]),sizeof(float), 3 * n_total, in);
  fread(&dumb,sizeof(int),1,in);  

  fclose(in);

  fprintf(stdout, "Eigenvector at point %lld [%g, %g, %g]\n", vecid, grid[vecid + 0], grid[vecid + 1], grid[vecid + 2]);

  exit(1);

  fprintf(stdout, "Eigenvector of grid point [0,0,0]->[%f,%f,%f]\n", grid[0], grid[1], grid[2]);
  fprintf(stdout, "Eigenvector of grid point [10,20,30]->[%f,%f,%f]\n", 
	  grid[0 + 3*(10 + n_x * (20 + n_y * 30))],
	  grid[1 + 3*(10 + n_x * (20 + n_y * 30))],
	  grid[2 + 3*(10 + n_x * (20 + n_y * 30))]);
  fprintf(stdout, "Eigenvector of grid point [%d,%d,%d]->[%f,%f,%f]\n", 
	  n_x-1 , n_y-1, n_z-1, 
	  grid[0 + 3*((n_x-1) + n_x * ((n_y-1) + n_y * (n_z-1)))],
	  grid[1 + 3*((n_x-1) + n_x * ((n_y-1) + n_y * (n_z-1)))],
	  grid[2 + 3*((n_x-1) + n_x * ((n_y-1) + n_y * (n_z-1)))]);
  

  min_val = 1.0E10;
  max_val = -1.0E10;
  for(i=0;i<n_x;i++){
    for(j=0;j<n_y;j++){
      for(k=0;k<n_z;k++){
	ind = 	2 + 3*(i + n_x * (j + n_y * k));
	if(grid[ind]<min_val){
	  min_val = grid[ind];
	  min_i = i;
	  min_j = j;
	  min_k = k;
	}
	if(grid[ind]>max_val){
	  max_val = grid[ind];
	}
      }
    }
  }
  fprintf(stdout, "min val along x component %g\n", min_val);
  fprintf(stdout, "min position: %i %i %i\n", min_i, min_j, min_k);
  ind = 3*(min_i + n_x * (min_j + n_y * min_k));
  fprintf(stdout, "full vector [%g,%g,%g]\n", grid[ind], grid[1+ind], grid[2+ind]);
  fprintf(stdout, "max val along x component %g\n", max_val);
  
  return 0;
}
