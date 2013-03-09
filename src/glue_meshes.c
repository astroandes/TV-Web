#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "io.h"
#include "field.h"
#include "READ_ART.h"
#include "tools.h"
#define USAGE "./glue.x <file_list_in> file_out"

int main(int argc, char **argv){
  int i,n_snap_in,n;
  density_grid *Grids;


  if(argc<3){
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
  }

  n_snap_in = argc - 2;
  if(!(Grids=malloc((n_snap_in+1) * sizeof(density_grid)))){
    fprintf(stderr, "Problem with the allocation of structure array\n");
    exit(1);
  }

  for(i=0;i<n_snap_in+1;i++){
    Grids[i].grid = NULL;
  }

  /*Initialize the final grid*/
  LoadDensityGrid(argv[1],&(Grids[n_snap_in]));
  for(n=0;n<Grids[n_snap_in].NNodes;n++){
    Grids[n_snap_in].grid[n] = 0.0;
  }  

  /*Load the grids
    loop over the grid items and add*/

  for(i=1;i<=n_snap_in;i++){
    LoadDensityGrid(argv[i],&(Grids[i-1]));
    fprintf(stdout, "Adding %d nodes\n", Grids[n_snap_in].NNodes);
    for(n=0;n<Grids[n_snap_in].NNodes;n++){
      Grids[n_snap_in].grid[n] += Grids[i-1].grid[n];
    }    
    free(Grids[i-1].grid);
  }
  
  SaveDensityGrid(argv[n_snap_in+1],&(Grids[n_snap_in]));

  return 0;
}
