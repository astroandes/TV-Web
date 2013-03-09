#include <stdio.h>
#include <stdlib.h>
#define USAGE "./load_env.x filename"
/*
   USAGE: ./load_env.x filename
   AUTHOR: Jaime Forero-Romero j.e.forero.romero@gmail.com
   DESCRIPTION:
   Loads a 3Dimensional environment grid into a 3D array.
   NOTES:
   n_x, n_y, n_z is the grid size
   x_0, y_0, z_0 is the posiion of the cell for the grid(i,j,k) element in kpc/h
   dx,dy,z is the grid size in kpc/h
   
   grid[n] runs over x,y,z in space, following fortran array conventions
   for memory handling
   
   i.e. a 3D point grid(i,j,k) corresponds to  grid[i + N_X * (j + N_Y * k)] where N_X is the grid dimension size in x and N_Y is the grid dimension size in y (i,j,k start at 0).

   grid[n] can take one of 4 values: 0,1,2,3
   3-> voids, 2-> sheets, 1->filaments, 0->peaks
*/
int main(int argc, char **argv){
  FILE *in;
  int *grid;
  int dumb;
  char line[30];
  int i;
  int n_x, n_y, n_z;
  int n_nodes;
  float dx, dy, dz, x_0, y_0, z_0;
  int n_void, n_sheet, n_fil, n_peak, n_bound;

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
  fread(&n_nodes,sizeof(int),1,in);    
  fread(&x_0,sizeof(float),1,in);    
  fread(&y_0,sizeof(float),1,in);    
  fread(&z_0,sizeof(float),1,in);    
  fread(&dx,sizeof(float),1,in);    
  fread(&dy,sizeof(float),1,in);    
  fread(&dz,sizeof(float),1,in);    
  fread(&dumb,sizeof(int),1,in);
  
  if(!(grid=malloc(n_nodes * sizeof(int)))){
    fprintf(stderr, "problem with array allocation\n");
    exit(1);
  }
  
  fread(&dumb,sizeof(int),1,in);
  fread(grid,sizeof(int), n_nodes, in);
  fread(&dumb,sizeof(int),1,in);  

  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  /*
    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);
  */
  fclose(in);


 n_void = 0;
 n_sheet = 0;
 n_fil = 0;
 n_peak = 0;
 n_bound = 0;
 for(i=0;i<n_nodes;i++){
   if(grid[i]==3){
     n_peak++;
   }else if(grid[i]==2){
     n_fil++;
   }else if(grid[i]==1){
     n_sheet++;
   }else if(grid[i]==0){
     n_void++;
   }else if(grid[i]==-1){
     n_bound++;
   }
 }
 fprintf(stdout, "n_void %d void volume filling fraction %g\n", n_void, 1.0*n_void/n_nodes);
 fprintf(stdout, "n_sheet %d void volume filling fraction %g\n", n_sheet, 1.0*n_sheet/n_nodes);
 fprintf(stdout, "n_fil %d void volume filling fraction %g\n", n_fil, 1.0*n_fil/n_nodes);
 fprintf(stdout, "n_peak %d void volume filling fraction %g\n", n_peak, 1.0*n_peak/n_nodes);
 fprintf(stdout, "n_bound %d void volume filling fraction %g\n", n_bound, 1.0*n_bound/n_nodes);

  return 0;
}
