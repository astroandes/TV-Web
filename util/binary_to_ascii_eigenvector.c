#include <stdio.h>
#include <stdlib.h>
#define USAGE "./load_env.x filename"
/*
   USAGE: ./load_cic.x filename
   AUTHOR: Jaime Forero-Romero j.e.forero.romero@gmail.com
   DESCRIPTION:
   Loads a 4Dimensional grid into a 4D array.
   NOTES:
   n_x, n_y, n_z is the grid size
   x_0, y_0, z_0 is the posiion of the cell for the vector(i,j,k) element in kpc/h
   dx,dy,z is the grid size in kpc/h
   
   vector[n] represents a 3D vector associated with a point in in space, following fortran array conventions
   for memory handling
   
   i.e. component C of the vector at point i,j,k in the grid is

    grid(i,j,k) corresponds to  vector[C + 3*(i + N_X * (j + N_Y * k))] where N_X is the grid dimension size in x and N_Y is the grid dimension size in y (i,j,k start at 0).

*/
#define FLOAT float
int main(int argc, char **argv){
  FILE *in;
  FLOAT *grid;
  int dumb;
  char line[30];
  long long i,j,k;
  long long n_total;
  int n_x, n_y, n_z;
  long long n_nodes;
  float dx, dy, dz, x_0, y_0, z_0;
  FLOAT max_val, min_val;
  FLOAT vec_x, vec_y, vec_z;
  if(argc!=2){
    fprintf(stderr, "USAGE: %s\n", USAGE);
    exit(1);
  }
  

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
  fread(&n_nodes,sizeof(long long),1,in);    
  fread(&x_0,sizeof(float),1,in);    
  fread(&y_0,sizeof(float),1,in);    
  fread(&z_0,sizeof(float),1,in);    
  fread(&dx,sizeof(float),1,in);    
  fread(&dy,sizeof(float),1,in);    
  fread(&dz,sizeof(float),1,in);    
  fread(&dumb,sizeof(int),1,in);

  n_total = n_x * n_y * n_z;
  if(!(grid=malloc(n_total * 3 * sizeof(FLOAT)))){
    fprintf(stderr, "problem with array allocation\n");
    exit(1);
  }
  
  fread(&dumb,sizeof(int),1,in);
  fread(&(grid[0]),sizeof(FLOAT), n_nodes*3, in);
  fread(&dumb,sizeof(int),1,in);  

  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  fclose(in);


  for(i=0;i<n_x;i++){
    for(j=0;j<n_y;j++){
      for(k=0;k<n_z;k++){
	vec_x = grid[0 + 3*(i + n_x * (j + n_y * k))];
	vec_y = grid[1 + 3*(i + n_x * (j + n_y * k))];
	vec_z = grid[2 + 3*(i + n_x * (j + n_y * k))];
	fprintf(stdout, "%d %d %d %f %f %f\n", i, j, k, vec_x, vec_y, vec_z);
      }
    }
  }


  return 0;
}
