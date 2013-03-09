#ifndef __MCUBES
#define __MCUBES

#include "struct.h"

#define MC_DELTA 2

#define MCFLAG_GETNORMALS     (1<<0) /* To compute the normals */
#define MCFLAG_GETVERTEXTREE  (1<<1) /* To compoute the point neighours */
#define MCFLAG_GETGRAD        (1<<2) /* To compute the value of the norm² of the gradiant for every vertex  */
#define MCFLAG_MARKFACES      (1<<3) /* for every face, record the number of the box it belonged to */
#define MCFLAG_CORRECT        (1<<4)
#define MCFLAG_GETVOLUME      (1<<5)

typedef struct MC_Object_str
{
    float *Vertex; //{x1,y1,z1,x2,y2,z2,...,Xn,Yn,Zn}
    int *Face; //3 Vertices in every face{f1.1,f1.2,f1.3,f2.1,...,fn.3}
    //Face is an index in Vertex, if the n^th point is the second vertex in the i^th face, Face[3*i+1]=3*n   
    
    char *Invalid;

    float *Normal; //the normal to the surface at a vertex (same structure as Point)
    float *Color; //color (same struct as Vertex {r1,g1,b1,...})
    float *alpha; //transparency

    int NFaces; //number of faces
    int NVertices; //number of vertices

    int *BoxIndex; //BoxIndex[i] is the index of the cell the ith face belonged to

    float *grad; //Value of the interpolated gradient of the field  at every vertex (x1,y1,z1,...)
    float *gradNorm;
    int *NNeighbours; //Number of direct neighbours to a given vertex
    int *VertexNeighbour[28]; //index of the n^th neighbours of vertex v (PointNeighbours[n][v])
    //vertices v,PointNeighbours[v][n] and PointNeighbours[v][n+1] form a face.
    double VolumeDown;
    double VolumeUp;
} MC_Object;

int InitMCubes(void);
int GetIsoSurface(density_grid *,float,MC_Object *,int);
int free_MC_Object(MC_Object **);

#endif
