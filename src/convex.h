#include <stdio.h>
#include <math.h>
#include <stdlib.h> /* for srand() */
#include <string.h>
#include <time.h> /* for random seeding */
#include <ctype.h>
#include <unistd.h>
#include <ginac/ginac.h>
#include <iostream>
#include <glpk.h>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
using namespace std;
using namespace GiNaC;


#define NR_END 1
#define FREE_ARG char*

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

long *booltoinds(bool *vec, int len, int *n1);
int isinconvhull(int **imat1, long *inds, int nlen, int mlen, int *vec);
int isinspan(int **imat1, long *inds, int nlen, int mlen, int *vec, bool pos);
int isinspan(int **imat1, bool *boolinds, int nlen, int mlen, int *vec, bool pos);
int isincone(int **imat1, long *inds, int nlen, int mlen, int *vec);
int isincone(int **imat1, bool *boolinds, int nlen, int mlen, int *vec);
int isinrspan(int **imat1, long *inds, int nlen, int mlen, int *vec, bool pos);
int isinrspan(int **imat1, bool *boolinds, int nlen, int mlen, int *vec, bool pos);
int isinrcone(int **imat1, long *inds, int nlen, int mlen, int *vec);
int isinrcone(int **imat1, bool *boolinds, int nlen, int mlen, int *vec);
int isvert(int **imat1, int *inds, int nlen, int mlen, int ind);
int getvertboolT(int **imat1, int nlen, int mlen, bool *verts);
int hasposkervec(int **imat1, int nlen, int mlen, bool strict);
ex LCmult(ex *exlist, int n);
void posintkervec(int **imat, int nlen, int mlen, int intvec[], bool q);
void intkervec(int **imat, int nlen, int mlen, int intvec[], bool q);
int hasposimvec(int **imat1, int nlen, int mlen);
int hasposrkervec(int **imat1, int nlen, int mlen, bool strict);
int hasposlimvec(int **imat1, int nlen, int mlen);
bool has_nonz_signed_col(int **imat, int n, int m);
int **conemat(int **A, int n, int m, bool *verts, int *m1);
int affrank(int **A, int n, int m);
int getfacestruct(int **Sl, int n, int m, bool *verts, bool ***allfacets, int debug);
int getfacestruct2(int **Sl, int n, int m, bool *verts, bool ***allfacets, int debug);
int getfacestruct1(int **Sl, int n, int m, bool *verts, bool ***allfacets, int *rk, int **rankvec, int *totverts, int debug);
int getextremeraybool(int **imat1, int nlen, int mlen, bool *verts);
int getextremerayboolT(int **imat1, int nlen, int mlen, bool *verts);
int **Minkowski2(int **Sl1, int n, int m1, bool *b1, int **Sl2, int m2, bool *b2, int *numverts);
int **MinkowskiN(int **Sl, int n, int m, bool **b, int *yc, int k, int *numverts);
int PolytopeVol(int **S, int **Sl, int n, int m, int debug);
int **minkerbasis(int **M, int n, int m, int *tot, int *rk, int *deg, int debug);
int **poskerbasis(int **Si, int n, int m, int *tot, int *rk);
