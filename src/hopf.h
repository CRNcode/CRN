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

void elongate(int ***explst, int **cflst, long *r, int numv);
void homify(int **lst1, long r1, int **lst2, long r2, int numv);
ex ReJpiI(matrix J, int n, int numv, bool hom);
ex ImJpiI(matrix J, int n, int numv, int hom);
void ReImJpiI(ex tmp1, ex tmp2, int numv, int ***explst1, int **cflst1, long *r1, int ***explst2, int **cflst2, long *r2, bool hom);
void ReImJpiI(matrix J, int n, int numv, int ***explst1, int **cflst1, long *r1, int ***explst2, int **cflst2, long *r2, bool hom);
int BTpotential(matrix J, int n, int rk, int numv, int maxpppdegree, int *pppdegused, int debug);
ex lsts_to_poly(int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int numv);
ex remsqs(ex tmp, ex *sqs, int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int numv, int *used1, int *used2, int clq1, int clq2, bool q);
int J2pI_seq(ex tmp, ex *tmpout, int n, int numv, int *totiter, float *bestq, ex *best, int *lastsplitfail, int maxlevel, int maxcliques, int maxiter, bool q);
int J2pI_seqWrap(ex tmp, ex *sqs, ex *finalpos, int n, int numv, int *totiter, float *bestq, ex *best, int *lastsplitfail, int maxlevel, int maxcliques, int maxiter, bool q);
int admitsIpair2(matrix J, int numv, int maxpppdegree, int *pppdegused, int debug);
int admitsIpair3(matrix J, int numv, int maxpppdegree, int *pppdegused, int debug);
int admitsIpair4(matrix J, int numv, int maxpppdegree, int *pppdegused, int debug);
int J2pIwrap(matrix J, int nlen, int numv, int maxpppdegree, int *pppdegused, int debug);
void numhonspl(int ** explst, int *cflst, long r, int numv, int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int *out);
int honestsplits(long indx, int n, int **allt, int *cfs, long r, int **lst1, int *cfs1, int r1, int **lst2, int *cfs2, int r2, int *used1, int *used2, int **needcount1, int **needs1, int **neededby1, int **needcount2, int **needs2, int **neededby2, int **needexcl, bool q);
void modcliques(int k, int *used1, long r1, int *used2, long r2, int ul, int ur, int *clq1, int *clq2, bool q);
int negcnt(ex tmp);
bool isprod(int *a, int *b, int *prod, int numv);
int findpairsplit1(ex tmp, int **explst, long t1, long *u1, long *u2, int **explst1, int **explst2, int *cflst1, int *cflst2, long r1, long r2, int *used1, int *used2, int clq1, int clq2, int numv, bool best, bool q);
bool ueq2v(int *u, int *v, int n);
bool flgset(int newf, int a, int b);
int validpair1(int **explst, long t1, long t2, long *u1, long *u2, int **explst1, int **explst2, int *cflst1, int *cflst2, long r1, long r2, int *used1, int *used2, int numv, int newf);
int randseedpair(int **explst, int **explst1, int **explst2, int *cflst, int *cflst1, int *cflst2, long r, long r1, long r2, int *used1, int *used2, int numv, int **needs1, int **needs2, long *ul, long *ur, int flg);
int randseedpair_old(int **explst, int **explst1, int **explst2, int *cflst, int *cflst1, int *cflst2, long r, long r1, long r2, int *used1, int *used2, int numv, int **needs, int *needexcl, long *ul, long *ur);
void forcedclq(long indx, long r, int *needexcl, int **needs, int **neededby, int *used, int *tot);
int isfreehonest(int **explst, int *cflst, int numv, long r, int **explst1, int *cflst1, int r1, int **explst2, int *cflst2, int r2, int *used, int *used1, int *used2, int *used1tmp, int *used2tmp);
ex honestforced(int **explst, int *cflst, int numv, long r, int *needexcl, int **needs, int **neededby,int **explst1, int *cflst1, int r1, int **explst2, int *cflst2, int r2, int *used1, int *used2, int *used1tmp, int *used2tmp, bool hon, int firstbigall, bool q);
int BigSq1(ex tmp, ex tmpfull, ex *sqs, ex *finalpos, int n, int numv, int **explst1, int *cflst1, int *used1, long r1, int **explst2, int *cflst2, int *used2, long r2, int *clq1, int *clq2, int *level, bool q);
int BigSqWrap1(matrix J, int n, int numv, bool q);
int Hopfforbid(int **imatSi, int **imatSl, int nlen, int mlen, const char *filter, int effort, int *pppdegused, int debug);
int HopfForbid(char *str, const char *filter, int effort, int *pppdegused, int debug);
int HopfForbid(char *di6, int n, int m, const char *filter, int effort, int *pppdegused, int debug);
int HopfForbid(int **AM, int n, int m, const char *filter, int effort, int *pppdegused, int debug);
int HopfForbid(char *di6, int n, int m, int effort, int *pppdegused, int debug);
void ruleouthopf(const char *fname, const char *outfhopf, const char *outfnohopf, int n, int m, const char *filter, int effort, int debug);
void ruleouthopf(const char *fname, const char *outfhopf, const char *outfnohopf, int n, int m, int effort, int debug);
