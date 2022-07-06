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

int genmod(int numv, int **mons, int q);
int ppp1(int **lst, double *cflst, long r, int numv, int polymindeg, int polymaxdeg, int preprocess, int ldegreemin, int ldegree, int strict, int debug);
int ppp(ex tmp1, int preprocess, int maxldegree, int strict, int *ldeg, int debug);
int pppfull(ex tmp, int numv, int preprocess, int maxldegree, int strict, int *ldeg, int debug);
int pppfull1(int **lst, double *cflst, long r, int numv, int polymindeg, int polymaxdeg, int preprocess, int ldegreemin, int ldegree, int strict, int debug);
int genmons(int numv, int maxdeg, int **mons);
int **genmons1(int numv, int mindeg, int maxdeg, int *tot);
int polyprod(int **lst1, int len1, double *cfs2, int **lst2, int len2, int numv, double ***lst1out, int ***lst2out);
int addoneintstra(int k, int *s, int len1, int ***t, int *ind);
int ceiling(int x, int y);
int symmattovec(int i, int j, int N);
void newconstraint(int ***constraints, int *numconstraints, int mon_ind);
void growconstraint(int **constraints, int ind, int blk, int row, int col);
double *genfloatvec(char *line, int *num);
int addpolyterm(int lstlen, int *newvec, double cf, int veclen, int ***veclst, double **cflst, int *pos);
int nearint(double p, int imin, int imax, double tol);
int **genevenpows(int numv, int mindeg, int maxdeg, int *tot);
int **genemptymonslist(int numv, int mindeg, int maxdeg, int *tot);


