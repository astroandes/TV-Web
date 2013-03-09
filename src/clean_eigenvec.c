#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
FLOAT *load_eigenvec(char * filename, int *n_x, int *n_y, int *n_z, 
		     int *n_nodes, float *x_0, float *y_0, float *z_0, 
		     float *dx, float *dy, float *dz, char *line
		     );
void dump_eigenvec(FLOAT *grid, char * filename, int *n_x, int *n_y, int *n_z, 
		   int *n_nodes, float *x_0, float *y_0, float *z_0, 
		   float *dx, float *dy, float *dz, char *line
		   );
int check_vector(FLOAT *vec);
void cross_product(FLOAT *vec_in_A, FLOAT *vec_in_B, FLOAT *vec_out_C);
int main(int argc, char **argv){

  FLOAT *grid_1;
  FLOAT *grid_2;
  FLOAT *grid_3;

  int n_x, n_y, n_z;
  int n_nodes;
  float x_0, y_0, z_0, dx, dy, dz;
  char line[30];
  long long ind;
  long long i, j, k, l;
  FLOAT vec_1[3];
  FLOAT vec_2[3];
  FLOAT vec_3[3];
  FLOAT new_vec[3];
  n_x = 0;
  n_y = 0;
  n_z = 0;

  if(argc<3){
    fprintf(stderr, "USAGE: %s\n", USAGE);
    exit(1);
  }
  if(!strcmp(argv[1], argv[4])||!strcmp(argv[2],argv[5])||!strcmp(argv[3],argv[6])){
   fprintf(stdout, "writing to\n");
    for(i=1;i<argc;i++){
      fprintf(stdout, "\t %s\n", argv[i]);
    }
    exit(1);
  }

  grid_1 = load_eigenvec(argv[1], &n_x, &n_y, &n_z, &n_nodes, &x_0, &y_0, &z_0, &dx, &dy, &dz, line);
  grid_2 = load_eigenvec(argv[2], &n_x, &n_y, &n_z, &n_nodes, &x_0, &y_0, &z_0, &dx, &dy, &dz, line);
  grid_3 = load_eigenvec(argv[3], &n_x, &n_y, &n_z, &n_nodes, &x_0, &y_0, &z_0, &dx, &dy, &dz, line);

  fprintf(stdout, "start cleaning\n");
  for(i=0;i<n_x;i++){
    for(j=0;j<n_y;j++){      
      for(k=0;k<n_z;k++){
	for(l=0;l<3;l++){
	  ind = l + 3*(i + n_x * (j + n_y * k));
	  vec_1[l] = 	grid_1[ind];
	  vec_2[l] = 	grid_2[ind];
	  vec_3[l] = 	grid_3[ind];
	  new_vec[l] = 0.0;
	}
	if(check_vector(vec_1)){
	  fprintf(stdout, "Problem 1: %g %g %g\n", vec_1[0], vec_1[1], vec_1[2]);
	  cross_product(&(vec_2[0]), &(vec_3[0]), &(new_vec[0]));
	  fprintf(stdout, "Solution 1: %g %g %g\n", new_vec[0], new_vec[1], new_vec[2]);
	  for(l=0;l<3;l++){
	    ind = l + 3*(i + n_x * (j + n_y * k));
	    grid_1[ind] = new_vec[l];
	  }
	}
	if(check_vector(vec_2)){
	  fprintf(stdout, "Problem 2: %g %g %g\n", vec_2[0], vec_2[1], vec_2[2]);
	  cross_product(&(vec_3[0]), &(vec_1[0]), &(new_vec[0]));
	  fprintf(stdout, "Solution 2: %g %g %g\n", new_vec[0], new_vec[1], new_vec[2]);
	  for(l=0;l<3;l++){
	    ind = l + 3*(i + n_x * (j + n_y * k));
	    grid_2[ind] = new_vec[l];
	  }
	}
	if(check_vector(vec_3)){
	  fprintf(stdout, "Problem 3: %g %g %g\n", vec_3[0], vec_3[1], vec_3[2]);
	  cross_product(&(vec_1[0]), &(vec_2[0]), &(new_vec[0]));
	  fprintf(stdout, "Solution 3: %g %g %g\n", new_vec[0], new_vec[1], new_vec[2]);
	  for(l=0;l<3;l++){
	    ind = l + 3*(i + n_x * (j + n_y * k));
	    grid_3[ind] = new_vec[l];
	  }
	}
      }
    }
  }

  dump_eigenvec(grid_1, argv[4], &n_x, &n_y, &n_z, &n_nodes, &x_0, &y_0, &z_0, &dx, &dy, &dz, line);
  dump_eigenvec(grid_2, argv[5], &n_x, &n_y, &n_z, &n_nodes, &x_0, &y_0, &z_0, &dx, &dy, &dz, line);
  dump_eigenvec(grid_3, argv[6], &n_x, &n_y, &n_z, &n_nodes, &x_0, &y_0, &z_0, &dx, &dy, &dz, line);
  fprintf(stdout, "finished cleaning\n");
  return 0;
}

void cross_product(FLOAT *vec_in_A, FLOAT *vec_in_B, FLOAT *vec_out_C){
  vec_out_C[0] = (vec_in_A[1]*vec_in_B[2] - vec_in_A[2]*vec_in_B[1]);
  vec_out_C[1] = -(vec_in_A[0]*vec_in_B[2] - vec_in_A[2]*vec_in_B[0]);
  vec_out_C[2] = (vec_in_A[0]*vec_in_B[1] - vec_in_A[1]*vec_in_B[0]);
}

int check_vector(FLOAT *vec){
  int GOOD;
  int BAD;
  GOOD = 0;
  BAD = 1;
  //  if(vec[0]<-1.0||vec[1]<-1.0||vec[2]<-1.0||vec[0]>1.0 ||vec[1]>1.0||vec[2]>1.0){
  if(isnan(vec[0])||isnan(vec[1])||isnan(vec[2])){
    return BAD;
  }else{
    return GOOD;
  }
}

FLOAT *load_eigenvec(char * filename, int *n_x, int *n_y, int *n_z, 
		     int *n_nodes, float *x_0, float *y_0, float *z_0, 
		     float *dx, float *dy, float *dz, char *line
		     ){
  FILE *in;
  long long n_total;
  /*
  char line[30];

  int n_x, n_y, n_z;
  int n_nodes;
  float dx, dy, dz, x_0, y_0, z_0;
  */
  FLOAT *grid;
  int dumb;
  if(!(in=fopen(filename, "r"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  fprintf(stdout, "reading %s\n", filename);
  
  fread(&dumb,sizeof(int),1,in);
  fread(line,sizeof(char)*30,1,in);
  fread(&dumb,sizeof(int),1,in);

  fread(&dumb,sizeof(int),1,in);
  fread(n_x,sizeof(int),1,in);    
  fread(n_y,sizeof(int),1,in);    
  fread(n_z,sizeof(int),1,in);    
  fread(n_nodes,sizeof(int),1,in);    
  fread(x_0,sizeof(float),1,in);    
  fread(y_0,sizeof(float),1,in);    
  fread(z_0,sizeof(float),1,in);    
  fread(dx,sizeof(float),1,in);    
  fread(dy,sizeof(float),1,in);    
  fread(dz,sizeof(float),1,in);    
  fread(&dumb,sizeof(int),1,in);
  n_total = (*n_x)*( *n_y) * (*n_z);
  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", *n_x, *n_y, *n_z, *n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", *x_0, *y_0, *z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", *dx, *dy, *dz);
  /*  *nx = n_x; *ny = n_y; *nz = n_z;
  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  */
  fprintf(stderr, "Ntotal %lld\n", n_total);

  if(!(grid=malloc(3 * n_total * sizeof(float)))){
    fprintf(stderr, "problem with array allocation\n");
    exit(1);
  }
  
  fread(&dumb,sizeof(int),1,in);
  fread(&(grid[0]),sizeof(float), 3 * n_total, in);
  fread(&dumb,sizeof(int),1,in);  

  fclose(in);
  return grid;
}


void dump_eigenvec(FLOAT *grid, char * filename, int *n_x, int *n_y, int *n_z, 
		   int *n_nodes, float *x_0, float *y_0, float *z_0, 
		   float *dx, float *dy, float *dz, char *line
		   ){
  FILE *in;
  long long n_total;

  int dumb=4;
  if(!(in=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  fprintf(stdout, "writing %s\n", filename);
  
  fwrite(&dumb,sizeof(int),1,in);
  fwrite(line,sizeof(char)*30,1,in);
  fwrite(&dumb,sizeof(int),1,in);

  fwrite(&dumb,sizeof(int),1,in);
  fwrite(n_x,sizeof(int),1,in);    
  fwrite(n_y,sizeof(int),1,in);    
  fwrite(n_z,sizeof(int),1,in);    
  fwrite(n_nodes,sizeof(int),1,in);    
  fwrite(x_0,sizeof(float),1,in);    
  fwrite(y_0,sizeof(float),1,in);    
  fwrite(z_0,sizeof(float),1,in);    
  fwrite(dx,sizeof(float),1,in);    
  fwrite(dy,sizeof(float),1,in);    
  fwrite(dz,sizeof(float),1,in);    
  fwrite(&dumb,sizeof(int),1,in);
  n_total = (*n_x)*( *n_y) * (*n_z);
  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", *n_x, *n_y, *n_z, *n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", *x_0, *y_0, *z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", *dx, *dy, *dz);
  /*  *nx = n_x; *ny = n_y; *nz = n_z;
  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  */
  fprintf(stderr, "Ntotal %lld\n", n_total);

  
  fwrite(&dumb,sizeof(int),1,in);
  fwrite(&(grid[0]),sizeof(float), 3 * n_total, in);
  fwrite(&dumb,sizeof(int),1,in);  

  fclose(in);
}
