#ifndef __TREE_STRUCT_H
#define __TREE_STRUCT_H

struct Tree_str 
{
    int N; //number of this point
    int Group; //which graph it belongs to
    int NbOthers; //number of previous points
    struct Tree_str **Other; //pointer to eventual other points linked to it
    struct Tree_str *Next; //pointer to the nearest point
    float d; //distance to the point pointed by Next (closest)
};

typedef struct Tree_str Tree;

typedef struct Grid_str
{
    int N; //Number of points

    int NDx; //Number of divisions on X axis
    int NDy;
    int NDz;
    
    double DMax; //maximum diameter of a sphere fully contained in one cell

    int NC; //Number of cells

    float XMax,XMin,YMax,YMin,ZMax,ZMin; //box coordinates
    float Sx,Sy,Sz; //size of the box

    float *x,*y,*z; //points coordinates

    int *WhichCell; //the cell in which a point is 
    int *NInCell; //number of points in a given cell
    int **InCell; //list of the points in a given cell
    
    int *CellXCoord; //X coordinate of cell N (in the grid)
    int *CellYCoord; //Y coordinate of cell N (in the grid)
    int *CellZCoord; //Z coordinate of cell N (in the grid)

    Tree *mst; //minimal spanning tree (an array and linked list at the same time, wonderfull no ?)

} Grid;

#endif
