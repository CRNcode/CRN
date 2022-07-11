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


int sgn(int i);
int par(int k);
int gcd(int a, int b);
int gcd(int *a, unsigned long num);
void gcdlist(int *lstin, unsigned long num, int *lstout);
int isend(char c);
bool iseven(int *vec, int n);
int binltoint(bool *vec, int l);
void inttobinl(int k1, bool *vec, int l);
long lbinltoint(bool *vec, int l);
void linttobinl(unsigned long k1, bool *vec, int l);

//Files and strings, etc
char *readfileintostr(const char fname[]);
bool filehasstr(char *fname, char *str);
unsigned long numflines(const char *fname, int *maxline);
int gtline(FILE *fp, char s[], int lim);
int getline0(FILE *fp, char s[], int lim);
char *getlinefromstr(long *pos, char *block);
char *getfnameraw(char *fname);
char *getdname(char *fname);
char *strchop2(char *t, int n, int n1);
void chop(char *str, char **str1, char **str2, const char *str3);
char *lrtrim(char s[]);
int startswith(char *s, char *lstr);
void reverse(char s[]);
void reverse(int s[], int tot);
int endswith(char *s, char *rstr);
int chemgts2(char *s, char ***v, char sep);
int isonlyspace(char *s);
int iscomline(char s[]);
int getnumints(char *s);
int getnumints1(char *s);
int getnumsegs(char *s);
char *getnthint(char *s, int n);
char *getnthint1(char *s, int n);
char *getnthwd(char *s, int n);
int getint(char **s);
int *getintvec(char *str);
unsigned long *getulvec(char *str);
int ispureint(char *s);

//output
void printvec(int *ivec, int n);
void printvec(double *ivec, int n);
void printvec(bool *ivec, int n);
void printvec(long *ivec, int n);
void printvec(unsigned long *ivec, int n);
void printvec1(int *ivec, int n);
void printvec1(double *ivec, int n);
void printvec1(bool *ivec, int n);
void printvec1(long *ivec, int n);
void printmat(int **imat, int n, int m);
void printmat(bool **imat, int n, int m);
void printmat(matrix imat, int n, int m);
void printmat(char **cmat, int **imat, int n, int m);
void printmat1(int **imat, int n, int m);
void printmat1(bool **imat, int n, int m);
void printmat1(matrix imat, int n, int m);
void printmat1(char **cmat, int **imat, int n, int m);
void printmat1(matrix imat, int n, int m);
void printmaximamat(int **imat, int n, int m);
void printmaximamat(matrix mat, int n, int m);
void printmaximamat1(int **imat, int n, int m);
void printmaximamat1(matrix mat, int n, int m);
void printtens(int ***tens, int n, int m);
void fprintvec(FILE *fd, int *ivec, int n);
void fprintveca(FILE *fd, int *ivec, int n);
void printsubmat(int **imat, int *vec1, int *vec2, int k1, int k2);
void printsubmat(matrix imat, int *vec1, int *vec2, int n, int m);
void printindsubmat(int **imat, int *vec1, int *vec2, int k1, int k2);
void printmat1(char **cmat, int **imat, int n, int m);
void printmatpat(char **cmat, int **imat, int n, int m);
void listprint(int *cfs, int **mons, int dim, long len);

// Lists and arrays
int isinlist(int i, int ilst[], int tot);
int isinlist(long i, long ilst[], int tot);
int areequal(int *vec1, int *vec2, int n);
int areequal(double *vec1, double *vec2, int n);
bool issubvec(bool *vec1, bool *vec2, int n);
int issubvec(bool *vec1, bool **vec2, int strt, int finish, int n);
bool unordlistareeq(int *a, int na, int *b, int nb);
bool unordlistareeq(long *a, int na, long *b, int nb);

void inittozero(int *vec, int len);
void inittozero(bool *vec, int len);
void inittozero(double *vec, int len);
void inittozero(int **mat, int n, int m);
void inittoone(int *vec, int len);
void inittoone(bool *vec, int len);
void boolunion(bool *v, bool *v1, int n);
int addtovec(int *v, int *numv, int k);
int addv(int k, const char *s, char ***t, int *ind);
int addv1(int k, const char *s, char ***t);
int addv2(int k, char *s, char ***t);
int addnewtoarray(int ***outmat, int sz, int *xc, int len);
int addnewtoarray(long ***outmat, int sz, long *xc, int len);
void addnewto1Darray(int **outvec, int sz, int newint);
void addnewto1Darray(bool **outvec, int sz, bool newbool);
void addnewto1Darray(double **outvec, int sz, double newint);
void addnewto1Dlarray(long **outvec, int sz, long newint);
int addvec1(int lstlen, int *newvec, int veclen, int ***veclst, int check, int *pos);
int addvec1(int lstlen, double *newvec, int veclen, double ***veclst, int check, int *pos);
int addvec1(int lstlen, bool *newvec, int veclen, bool ***veclst, int check, int *pos);
int isinarray(char *v[], int numv, const char *s);
long isinarray(int **v, long numv, int *s, int n);
int isinarray(int **v, int numv, int *s, int n);
int isinarray(bool **v, int numv, bool *s, int n);
int isinarray(double **v, int numv, double *s, int n);
unsigned long isinarray(char *v[], unsigned long numv, const char *s);
int isinarrayoffset2(int **v, int numv, int *s, int n);
char **genarraydat(const char *datfile, int *num);
char **genarraydat(const char *datfile, unsigned long *num);


//allocation and freeing
char **charar(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
bool **bmatrix(long nrl, long nrh, long ncl, long nch);
void free_bmatrix(bool **m, long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_imat(int **A, int numA);
void free_bmat(bool **A, int numA);
ex **exmatrix(long nrl, long nrh, long ncl, long nch);
void free_exmatrix(ex **m, long nrl, long nrh, long ncl, long nch);
void ifree(int **imat, int n);
void ifree(double **imat, int n);
void ifree(int **imat, long n);
void ifree(double **imat, long n);
void lfree_i(long **imat, int n);
void ulfree_ul(unsigned long **imat, unsigned long n);
int freearraydat(char **array, int lim);
int freearraydat(char **array, unsigned long lim);

//search
long binisearch(int **v, long n, int *x, int vlen);
void growlist(int **lst, int *r);
void growlist(bool **lst, int *r);
void iswap(int *v, int i, int j);
void iswap(double *v, int i, int j);
void qsortt(int *v, int left, int right);
void qsort2(int *v, int *v1, int left, int right);
void iswapi(int **v, int i, int j);
void qsorti(int **v, int dim, long left, long right);
void qsorti2(int **v, int *cfs, int dim, long left, long right);
void qsorti2(int **v, double *cfs, int dim, long left, long right);
void ipswap(int **p, int i, int j);
void lswap(long *p, int i, int j);
void lpswap(long **p, int i, int j);

// permutations and combinations
int choose(int n, int k);
unsigned long factorial(int x);
unsigned long comb(int n, int k);
void firstcomb(int *vec, int n, int n1);
int nextcomb(int *vec, int n, int n1);
int **allperms(int n);
int **allperms1(int *vec, int n);
int **allcombsgen(int n, int n1);
int nextperm(int **vec, int **veclr, int *par, int n);
int firstcombfrom(int *vec, int n, int n1, int i_init);
void nextcombk(int *vec, int n, int n1, long k);
int nextnum(int *vec, int n, int base);
int nextnumb(int *vec, int n, int *base);
long base10(int *vec, int n, int base);
void basek(int *vec, int n, int base, long num);
void nextnumk(int *vec, int n, int base, long k);

//graphs
void strongconnect(int i, int **mat, int n, int *inds, int *lowlink, int *ind, int ***S, int *totS, int *Sc, int *Scpos);
int **Tarjan(int **mat, int n, int *totS);
void strongconnect_light(int i, int **mat, int n, int *inds, int *lowlink, int *ind, int *totS, int *Sc, int *Scpos);
int Tarjan_light(int **mat, int n);

//vectors and matrices
int **submat(int **imat, int n, int m, int *xc, int n1, int *yc, int m1);
int **submatfrombool(int **A, int n, int m, bool *vecx, bool *vecy, int *n1, int *m1);
int **submatfromcolbool(int **A, int n, int m, bool *vecy, int *m1);
int **submatfromrowbool(int **A, int n, int m, bool *vecx, int *n1);
int **submat1(int **imat, int n, int i1, int i2);
int **submatminus(int **mat, int n, int m);
int **submatgen(int **imat, int n, int m, int *i1, int *i2, int dim);

int unionlvec(long *a, int na, long *b, int nb, long **outvec);

void veccp(int *cflst,long r,int *cflst1);
void veccp(int *vecin,int r,int *vecout);
int *veccp(int *vecin,int r);
void veccp(double *vecin,int r,double *vecout);
void veccp(bool *vecin,int r,bool *vecout);
int **cpmat(int **mat, int n, int m);
void cpmat(int **mat, int **out, int n, int m);
matrix minusmat(matrix Jin, int n, int m);
void minusmat(int **a, int n, int m);
int **idmat(int n, bool minus);

int isgrtr(int *vec1, int *vec2, int n);
int isgrtrf(int *vec1, int *vec2, int n);

int *halve(int *a, int n);
void vecadd(int *a, int *b, int n);
void vecadd(double *a, double *b, int n);
void matadd(int **a, int **b, int n, int m);
void vecsum2(int *out, int *a, int *b, int numv);
void vecsum3(int *out, int *a, int *b, int *c, int numv);
int *subtract(int *a, int *b, int numv);
int **subtract(int **a, int **b, int n, int m);

int nonzentries(int **imat, int n, int m);
int nonzentries(bool *bvec, int n);
int nonzentries(int *ivec, int n);
void multAB(int **A, int **B, int **C, int n, int r, int m);
matrix multAB(matrix A, int **B, int n, int r, int m);
matrix multAB(int **A, matrix B, int n, int r, int m);
matrix multAB(matrix A, matrix B, int n, int r, int m);
matrix multAB(matrix S, matrix V, int n, int m);
matrix multABT(matrix S, matrix V, int n, int m);
matrix multABT(matrix S, int **V, int n, int m);
matrix multABT(int **S, matrix V, int n, int m);
int **multABT(int **S, int **V, int n, int m);
int gcd_vec(int *v, int n);
int reduce_vec(int *v, int n);
bool disjointcols(int **A, int n, int m, int *xc, int len);
bool disjointcols(int **A, int n, int m, int i, int j);
bool signcompatcols(int **A, int n, int m, int *xc, int len, bool minus);
bool signcompatcols(int **A, int n, int m, int i, int j, bool minus);
bool isreversecol(int **S, int n, int j, int k);
bool isequalcol(int **S, int n, int j, int k);
void sgn_mat(int **imat, int n, int m);
bool reduce_mat(int **imat, int n, int m, bool transpose);
bool hasposdiag(int **imat, int n);
int rowmax(int **imat, int n, int m, int rowk);
int colmax(int **imat, int n, int m, int colk);
int **scalecols(int **imat, int n, int m, int *vec);
matrix scalecols(matrix imat, int n, int m, matrix vec);
matrix scalecols(int **imat, int n, int m, matrix vec);
void scalerows(int **imat, int n, int m, int *vec);
matrix scalerows(int **imat, int n, int m, matrix vec);
matrix scalerows(matrix imat, int n, int m, matrix vec);
int **transposemat(int **imat, int n, int m);
void transposemat(int **imatin, int n, int m, int **imatout);
matrix transposemat(matrix imat, int n, int m);
void translatematbyrow(int **imatin, int n, int m, int rownum, int **imatout);
void translatematbycol(int **imatin, int n, int m, int colnum, int **imatout);
int **doublemat(int **imat, int n, int m);
int colsum(int **imat1, int nlen, int mlen, int j);
int rowsum(int **imat1, int nlen, int mlen, int j);
matrix imattoexmat(int **A, int n, int m);
matrix isubmattoexmat(int **A, int *vec1, int n1, int *vec2, int m1);
int rowdegree(int **imat, int n, int m);
int coldegree(int **imat, int n, int m);
int det(int **imat, int n);
ex det(matrix imat, int n);
int detsubmat(int **imat, int n, int m, int *i1, int *i2, int dim);
ex symbdetsubmat(matrix imat, int n, int m, int *i1, int *i2, int dim);
int matrank(int **A, int n, int m);
int submatrank(int **A, int *vec1, int n, int *vec2, int m);
int submatrank(int **A, int n, int m, bool *vecx, bool *vecy);
matrix AdComp2(matrix J, int ktot);
ex AdComp2Det(matrix J, int n, int debug);
ex getminorsum0(matrix M, int n, int k);
