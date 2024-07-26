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

int ispospoly(ex detx, int numv, int *allflg, bool isfct, int doublecfs, int maxpppdeg, int *pppdegused, int debug);
int ispospolynonneg(ex detx, int numv, int *allflg, bool isfct, int doublecfs, int maxpppdeg, int *pppdegused, int debug);
int ispospolynonneg(ex detx, int numv, int *allflg, bool isfct, int doublecfs, int debug);
int isfullpospoly(ex detex, int numv, bool isfct, int maxpppdeg, int *pppdegused, int nonnegonly, int debug);
int NewtonPolyVertSgn1(ex tmp, int numv, int *allflg, int debug);
int J2pI_recurr(ex tmp, ex orig_poly, ex rems[], int numv, int *level, int *totiter, float *bestq, ex *best, int *lastsplitfail, int maxlevel, int maxcliques, int maxiter, int debug);
int heuristic_squares(ex detex, int numv, int debug);
int canbezero(ex detex, int numv, bool debug);
int testposneg(ex tmp, double min, double max, int tottries, double tol);
int AMGMsimp(ex detex, int numv, int debug);
int *sympt(int *vecmid, int *vecl, int cf, int n, bool pqswitch);
int splt(int a, int b);
float numsplits(long indx, int n, int **allt, int *cfs, long r, int **needcount, int **needexcl, int **needs, int **neededby, int debug);
int testforpos(int **lst, int *cfs, int polytot, int numv, double min, double max, int tottries, double tol);
int testforpos(int **lst, double *cfs, int polytot, int numv, double min, double max, int tottries, double tol);
int hasvertnonsq(int **lst, int numv, long r, int *allsq);
int NewtonPolyVertPosSq(ex tmp, int numv);
int NewtonPolyVertPosSq_a(int **lst, int *cflst, int numv, long r);
int NewtonPolyVertPosSq_a(int **lst, double *cflst, int numv, long r);
int vertsgns(int **lst, int *cflst, int numv, long r, int *allpos, int debug);
int vertsgns(int **lst, double *cflst, int numv, long r, int *allpos, int debug);
int *midpt(int *vecl, int *vecr, int cfl, int cfr, int n);
void gather(int **allt, int *cfs, long r, int n, int **neededby,int **needs,int *needexcl,long **lt,long **rt,long **neg,long **mid,long **all,int *numlt,int *numrt,int *numneg,int *nummid,int m0,int startind,int maxsz, int *numsofar, bool onlft, int *failflag, bool tight, bool excflag);
void sizesort4(int *cliquesz, long **lt, long **rt, int *numlt, int *numrt, int *absnegmax, int left, int right);
void lkeeponly(long **imat, int n, int m);
int allcliques1(int **allt, int *cfs, long r, int n, int **neededby, int **needs, int *needexcl, long ***left, long ***right, int **numleft, int **numright, int **absnegmax, int maxsz, int maxgetcliques, bool tight, bool onlycmplt, bool excflag, int debug);
ex makesquare(int **l, int *cfs, long *leftinds, int numleft, long *rightinds, int numright, int n, int maxnegcf, bool simp);
double devalpoly(int **lst, int *inds, int *cflst, int cfsindexed, int tot, int numv, int offset, double *vals);
double devalpoly(int **lst, int *inds, double *cflst, int cfsindexed, int tot, int numv, int offset, double *vals);
int notsq(int *lst, int numv);
ex facsq(int **l, long *leftinds, int *kleft, int numleft, long *rightinds, int *kright, int numright, int numv, int minfac);
int *intersect(int **l, long *leftinds, int numleft, long *rightinds, int numright, int numv);
int DetSgn(matrix J, int n, int numv, int maxpppdeg, int debug);
int DetSgnk(matrix J, int k, int n, int numv, int maxpppdeg, int debug);
int isPmatrix(matrix J, int n, int debug);
int isPmatrixstrict(int **J, int n);
int isposdef(int **J, int n);
int isposdef(matrix J, int n, int maxpppdeg, int *pppdegused, int debug);
int isP0matorth(matrix J, int n, bool skiptop, int maxpppdeg, int *pppdegused, int debug);
int isP0matorth(int **J, int n, bool skiptop, int debug);
int isP0matorth(matrix J, int n, int maxpppdeg, int debug);
int isP0matorth(int **J, int n, int debug);
int mixedtrace(matrix J, int n, int maxpppdeg, int *pppdegused, int debug);
int mixedtrace(int **J, int n, int debug);
int isQmatrix(matrix J, int n);
int signsym(int **J, int n);
int signsym(matrix J, int n, int maxpppdeg, int debug);
int AdComp2DetPos(matrix J, int n, int numv, int maxpppdeg, int *pppdegused, int debug);
int AdComp2isP0(matrix J, int n, int maxpppdeg, int debug);
int polytest(ex tmp, int numv, int deg, int dcfs, const char *filter, int maxpppdeg, int debug);
int polytest(const char *file, const char *filter, int maxpppdeg, int debug);
int polymixedvol(const char fname[], int debug);
