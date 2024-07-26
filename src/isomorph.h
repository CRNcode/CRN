#include <stdio.h>
#include <math.h>
#include <stdlib.h> /* for srand() */
#include <string.h>
#include <time.h> /* for random seeding */
#include <ctype.h>
#include <ginac/ginac.h>
#include <glpk.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace GiNaC;


#define NR_END 1
#define FREE_ARG char*

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

int genallbicomplexes(int n, int N, int **l);
unsigned long genbimol(int n, int m, char *outfile, int open, int rev, const char procstr[], int debug);
unsigned long genbitrimol(int n, int m, char *outfile, int open, const char procstr[], int debug);
unsigned long genzeroonetritrimol(int n, int m, char *outfile, int open, const char procstr[], int debug);
unsigned long genbitetramol(int n, int m, char *outfile, int open, const char procstr[], int debug);
unsigned long genbimultimol(int n, int m, char *outfile, int open, const char procstr[], int molecularity, int debug);
bool goodorder(int **L, int **R, int n, int m);
int bin6toint(bool *vec);
void amtodig(bool *V, int l, char out[]);
int numshortg();
FILE *mergeandshort(const char fname[], const char ptn[], int numparts);
char *dynamic_isomorphMA(char *di6, int n, int m);
char *dynamic_isomorphGK(char *di6, int n, int m);
int dynamicshortfile(const char *infile, const char *outfile, const char *flag, int n, int m, bool q);
unsigned long buildsupCRN(char *infile, char *outfile, int n, int m, unsigned long *totin, bool allowopen, bool r, int rankchange, int debug);
bool CRNsubset(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5, int n, int m, int debug);
unsigned long getallisomorphs(char *isotablefile, char *infilepre, char *infileDI, char *infile, char *outfile, unsigned long *numin);
int CRNwithsources(const char *filein, const char *filesources, const char *fileout, int n, int m, int debug);
unsigned long genmultisrc(int n, int m, char *outfile, int molecularity, int debug);
