#include <stdio.h>
#include <math.h>
#include <stdlib.h> /* for srand() */
#include <string.h>
#include <time.h> /* for random seeding */
#include <ctype.h>
#include <unistd.h>
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

int checkrankLR(int **L, int **R, int n, int m);
int checkrankS(int **S, int n, int m);
int checkrank(int **AM, int n, int m);
int checkrank(char *di6, int n, int m);

int isSSD2(int **imat, int n, int m);
int isSSD1(int **imat, int n, int m, bool allm, int *rkbad, int q);
int isSSD(int **mat, int n, int m, int q);

int isWSD2(int **imat, int n, int m);
int isWSD1(int **imat, int n, int m, bool allm, int *rkbad, int q);
int isCSD(int **imat, int n, int m, int q);
int isWS(int **mat, int n);
int reversify(int **S, int **V, int n, int m, int ***Sout, int ***Vout, int *mout);
void simppair(int **imat, int **imat1, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1);
void WSpair(int **imat, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1);
int strmatisSNSSS(char *fname);
int strmatisSSD(char *fname, int q);
int strmatcheckS2minors(char *fname, char *fnameout);
int strmatisCSD(char *fname, int q);
int minorisWS(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k, unsigned long fk);
int minorhas0rc(int **imat, int n, int m, int *vec1, int *vec2, int k);
int doubleisWSD2(int **imat, int n, int m);
int doubleisWSD1(int **imat, int n, int m);
int doubleisWSD(int **imat, int n, int m, int q, bool rval);
int isWSD(int **mat, int n, int m, int q);
int mats_compat(int **imat1, int **imat2, int n, int m, bool allm, int q);
int **detsk(int **imat, int n, int m, int k);
void printsiphons(int **allsiphons, int totsiphons, char **chems);
void printsiphons(int **allsiphons, int totsiphons, int **allminsiphons, int totminsiphons, char **chems);
void printsiphons(int **allsiphons, int totsiphons);
void printsiphons(int **allsiphons, int totsiphons, int **allminsiphons, int totminsiphons);

bool structpersist(int **imat3, int nlen, int cols3, int **allminsiphons, int totminsiphons);
bool structpersist(int **Si, int **Sil, int n, int m, int debug);
bool structpersist(int **AM, int n, int m, int debug);
bool structpersist(char *di6, int n, int m, int debug);

long getpos(int *vec, int n, int r);
long getpos1(int *vec, int omit, int n, int r);
int minor1(int **imat, int *xcombs, int *ycombs, int n, int m, int k, int **dets);
int **detsk1(int **imat, int n, int m, int k, int **dets);
int ***allminors(int **imat, int n, int m);
void free_allminors(int ***t, int n, int m);
void printmaximaindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2);
int **simpmat(int **mat, int n, int m, int *n1, int *m1);
int **readmatrixfromstr(char *str, int *nlen, int *mlen);
int readmatpairfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ***mat2);
int mat_signpat_compat(int **imat1, int **imat2, int n, int m, int q);
int getreac(char *str, char ***leftchems, int **leftstoics, int *numleft, char ***rightchems, int **rightstoics, int *numright, int *rev);
int chemgts2a(char *s, char ***v, char sep);
int guessreacformat(char *fname);
int guessoutputformat(char *fname);
int getallreacs(char *str, int ***imat1, int ***imat2, int ***imat3, int ***imat4, int ***stoichl, int ***stoichr, char ***chems, bool *haszerocomplex, int *n, int *m, int *cols3, int *allrev, int *allgood, int debug);
int analysereacs(const char fname[], int q, bool htmlswitch, bool statswitch);
int fixedminorcompat(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k);
int allminorsigns(int **imat1, int **imat2, int n, int m, int q);
int genMAXMAreacs(char *fname, int **imat1, int **imat2, int n, int m);
int S2(int **S, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1);
int arecompat(int **imat1, int **imat2, int n, int m, int debug);
int arecompatsimp(int **imat1, int **imat2, bool maxonly, int n, int m, int *strict, int debug);
int S2a(int **S, int **V, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1);
int S2b(int **S, int **V, int n, int m, int ***imat1, int **valsvec, int **indsvec, int *base, int *basetot, int *n1, int *m1);
int S2arecompat(int **S, int **V, int n, int m, int q);
int S2isWSD(char *fname, int q, int maxonly);
int minorisSNSsing2(int **imat, int n, int m, int *vec1, int *vec2, int k);
int qualdetsubmat(int **imat, int n, int m, int *i1, int *i2, int dim);
int qualdetsubmat1(int **imat, int n, int m, int *i1, int *i2, int dim);
void printdotindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2, long lab);
int cuttails(int **mat, int n, int m, int *rw, int *col, int **n1, int **m1);
int addtwointstr(int k, int *s, int *s1, int len1, int len2, int ***t);
int S2detstocheck(char *fname, int q, int maxonly);
int findbad(int ** lists, int numlists, int *rowset, int latestcol, int level, int dim, bool **bmat, int **imat, matrix mmat, int n, int m, long *totbad, long *totall, int **hist);
ex intpairtopoly(int *l1, int cf1, int n1, int *l2, int cf2, int n2);
int **cpmat(int **mat, int n, int m);
int ReacIso(const char fname[], int q);
int SauroSplit(char *fname, char *prefix);
int *digraph6toam(const char d6[], int *n);
char *di6toreacstr(char *di6, int n, int m, int open);
int **SSltoAM(int **S, int **Sl, int n, int m, bool minus);
void AMtoSSl(int **AM, int n, int m, bool minus, int ***S, int ***Sl);
void di6toSSl(char *di6, int n, int m, bool minus, int ***S, int ***Sl);
int **di6toCRNam(char *di6, int n, int m, int *entries);
int **di6tosignpat(char *di6, int n);
char *di6toSauro(char *di6, int n, int m);
void di6toSauro2(char *di6, int n, int m, char *outstr);
void di6toSauro1(FILE *fd, char *di6, int n, int m);
unsigned long di6toSaurofile(const char *infile, const char *outfile, int n, int m);
unsigned long di6toreacfile(const char *infile, const char *outfile, int n, int m, int open);
unsigned long di6tosimpstrfile(const char *infile, const char *outfile, int n, int m);
unsigned long di6todotfile(const char *infile, const char *outfile, int n, int m);
unsigned long reactodi6file(const char *fname, const char *fout, const char *sepline, int *n, int *m);
unsigned long Saurotodi6file(const char *fname, const char *fout, int *n, int *m);
unsigned long simpstrtodi6file(const char *fname, const char *fout, int n, int m);
unsigned long convertformat(char *infile, int intype, char *outfile, int outtype, int *numspec, int *numreac);
char *CRNamtostr(int **V, int n, int m, int entries, int open);
void CRNamtofile(int **V, int n, int m, const char *outfname, int open);
char *cmplxCRN(char *di6, int n, int m, int minlayers, bool sourceonly, int *totcmplx);
matrix Vmatfromimat(int **imat, int n1, int m1, int *maxv);
matrix Vmatfromiimat(int ***imat, int n1, int m1, int *maxv);
matrix Vmatfrompatmat(int **imat, int n1, int m1, int *maxv);
matrix VmatMAfromlstoich(int **imat, int n1, int m1, int *maxv);
matrix MARHSfromlstoich(int **imatSi, int **imatSl, int n1, int m1);
matrix reacJac(int **AM, int n, int m, int *numv, bool minus);
matrix reacJac(char *di6, int n, int m, int *numv, bool minus);
bool CRNweakrev(int **V, int n, int m);
bool CRNweakrev2(int Rvec[][2], int *xc, int m);
bool connected(int **Sl, int **Sr, int n, int m);
unsigned long reacreportfile(char *fname, char *reactype, char *fitler, int n, int m, int maxpppdeg, int debug);
unsigned long filterCRNs(const char *fname, const char *outfname, int n, int m, const char filt[], int iconst, int debug);
int **DSR(char *di6, int n, int m);
void DSR(int **AM, int n, int m);
int **DSRfromPN(int **AM, int n, int m);
int **DSR(int **S, int **V, int n, int m, int minus);
int DSRCondStar(int **DSRAM, int n, int m, int *report, int debug);
int DSRCondStarPN(int **PNAM, int n, int m, int *report, int rev, int debug);
int DSRCondStar(char *di6, int n, int m, int *report, int rev, int debug);
int DSRCondStar(int **S, int **V, int n, int m, int minus, int *report, int rev, int debug);
int DSR2CondStar(int **S, int **V, int n, int m, int minus, int *report, int rev, int debug);
int DSR2CondStar(int **DSRAM, int n, int m, int *report, int rev, int debug);
int DSR2CondStar(char *di6, int n, int m, int *report, int rev, int debug);

matrix reacJMAeq(int **Si, int **Sil, int n, int m, int ***Q, matrix *QX, int *numv);
matrix reacJMAeq(int **AM, int n, int m, int ***Q, matrix *QX, int *numv);
matrix reacJMAeq(char *di6, int n, int m, int ***Q, matrix *QX, int *numv);
int degen(matrix QX, int n, int rk, int level, int maxpppdeg, int *pppdegused, int debug);
int degen(int **QX, int n, int rk, int level, int debug);
int degen2(matrix QX, int n, int rk, int maxpppdeg, int *pppdegused, int debug);
int degen2(int **QX, int n, int rk, int debug);
int idegen(matrix QX, int n, int rk, int level, int debug);
int idegen(int **QX, int n, int rk, int level, int debug);
int reacQMAposdef(int **Q, matrix QX, int n, bool minus, bool strict, int maxpppdeg, int debug);
int reacQMAposdef(int **AM, int n, int m, bool minus, int maxpppdeg, bool strict);
int reacQMAposdef(char *di6, int n, int m, bool minus, int maxpppdeg, bool strict);

int **poskerbasis(int **Si, int n, int m, int *tot, int *rk);
int reacJMAdegenerate(int **S, int **Sl, int n, int m, int debug);
int reacJMAdegenerate(int **AM, int n, int m, int debug);
int reacJMAdegenerate(char *di6, int n, int m, int debug);
int genS2from2(int **imatS, int **imatV, int n, int m, int ***mat1, int ***mat2);
int genS2from2(int **imatS, matrix imatV, int n, int m, int ***mat1, matrix *mat2);

void CRNPN2(int **SR, int **RS, int n, int m, bool *V1, int Vlen);
bool *CRNPN3(int **SR, int **RS, int n, int m, int minlayers, int *layers);
int weakendo(int **S, int **Sl, int n, int m, int debug);
int endotactic_part(int **S, int **Sl, int n, int m, int debug);
int endotactic(int **S, int **Sl, int n, int m, int debug);
int endotactic(int **AM, int n, int m, int debug);
int endotactic(char *di6, int n, int m, int debug);
