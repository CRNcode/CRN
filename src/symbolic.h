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

const possymbol & get_possymbol(const string & s);
const symbol & get_symbol(const string & s);
int isapower(ex e);
ex getroot(ex e, int *powval);
ex gcdpoly(ex e);
void monovtomonoexp(int *k, int r, int *out, int numv);
void monovtomonoexp0(int *k, int *out, int numv);
lst monotolist(ex e);
void monotointlist(int *cf, ex e, int **lst, int *r);
void monotointlist(double *cf, ex e, int **lst, int *r);
void checkmonotolist();
void checkmonotointlist();
int monotointlist1(ex e, char **vars, int numv, int **lst, int *r);
void monotointlist0(int *p, ex e, int **lst);
void monotointlist0(double *p, ex e, int **lst);
int polytointlist(ex e, int ***lst, int **cflst, long *r);
int polytointlist(ex e, int ***lst, double **cflst, long *r);
void polytointlist0(ex e, int ***lst, int **cflst, long *r);
void polytointlist0(ex e, int ***lst, double **cflst, long *r);
int sumtointlist(ex e, int **lst);
int **monlisttoexplist(int **mons, int dim, long len, int numv, bool srt);
int **moncflisttoexplist(int **mons, int dim, int *cfs, long len, int numv, bool srt);
int **moncflisttoexplist(int **mons, int dim, double *cfs, long len, int numv, bool srt);
int **moncflisttoexplist0(int **mons, int *cfs, long len, int numv, bool srt);
int **moncflisttoexplist0(int **mons, double *cfs, long len, int numv, bool srt);
void monovars(ex e, char ***vars, int *numvars);
ex monosimp(ex e, char ***vars, int *numvars);
void monovarss(ex e, ex *vars, char **varnames, int numv);
int polyvars(ex e, char ***vars);
int polyvars(ex e, char ***vars, int *numv);
ex polysimp(ex e, char ***vars, int *numv);
int polyvars1(ex e, char ***vars);
void polyvarss(ex e, ex *vars, char **varnames, int numv);
void extolst(ex tmp, int numv, int ***explst, int **cflst, long *r);
void extolst(ex tmp, int numv, int ***explst, double **cflst, long *r);
void extolst1(ex tmp, int numv, int ***explst, int **cflst, long *r, int *mindegree, int *maxdegree);
void extolst1(ex tmp, int numv, int ***explst, double **cflst, long *r, int *mindegree, int *maxdegree);
ex intltomono(int cf, int *l, int n);
ex expltomono(int cf, int *l, int n);
void expltomonoprint(int cf, int *l, int n);
int polylistsimp(int **lst, int numv, long r, int ***lstnew, int mkhom, int tmppolymindeg, int tmppolymaxdeg, int *polymindeg, int *polymaxdeg);
ex lrexptopolysimp(int **l, int *cfs, long numinds, int n);
ex lrexptopolysimp(int **l, double *cfs, long numinds, int n);
ex lrexptopolysimp(int **l, double *cfs, double tol, long numinds, int n);
ex polyhom(ex e, ex var, int offset);
ex polyhomsimp(ex tmp1, int *numvin, int *numvout, int *deg, int mkhom, int dcfs, int debug);
ex polyfromstring(const char *str, int *numvin, int *numv, int *deg, int mkhom, int *dcfs, int debug);
ex polyfromfile(const char filename[], int *numvin, int *numv, int *deg, int mkhom, int *dcfs, int debug);
void polysfromfile(const char fname[], int ***mons, double ***cfs, int *numv, long *totmons, int *totpols, int debug);
matrix makexvec(int n);
matrix makexvec(int j, int n);
matrix makediagxmat(int k, int n, bool minus);
matrix makexvect(int j, int n);
matrix makexvecinv(int n);
matrix makexveccomp(int n);
