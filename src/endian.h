#ifndef __ENDIAN_SWAP
#define __ENDIAN_SWAP

#include <stdio.h>

int fread_sw(void *,int,int,FILE *,int);

int swapI(int);
float swapF(float);
double swapD(double);
inline void Dswap2B(void*);
inline void Dswap4B(void*);
inline void Dswap8B(void*);
void Dswap2BArr(void*,int);
void Dswap4BArr(void*,int);
void Dswap8BArr(void*,int);

#endif
