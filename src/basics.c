/* Copyright (C) 2010-2024, Murad Banaji
 *
 * This file is part of CRNcode
 *
 * CRNcode is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * CRNcode is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CRNcode: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 

 */

#include "basics.h"
#include <limits.h>

//
//
// Numbers, lists and strings
//
//

int sgn(int i){
  if(i==0)
    return 0;
  if(i>0)
    return 1;
  return -1;
}

//parity
int par(int k){
  if (k%2==0)
    return 1;
  return -1;
}

//Euclid's algorithm
int gcd(int a, int b){
  int temp;
  if(a<=0 || b <=0){
    fprintf(stderr, "the routine gcd requires positive integers.\n");
    exit(0);
  }
  while(b){
    temp=a%b;a=b;b=temp;
  }

  return a;
}

//overloading: allow integers <=0, and a list of arbitrary length
int gcd(int *a, unsigned long num){
  unsigned long i=0;
  int g=1;//in case of only zeros
  if(num<=0){
    fprintf(stderr, "ERROR in gcd: requires at least one integer. EXITING.\n");
    exit(0);
  }

  while(i<num-1 && a[i]==0)//skip initial zeros
    i++;
  if(a[i])
    g=abs(a[i]);
  while(i<num-1){
    if(a[i+1])//nonzero
      g=gcd(g,abs(a[i+1]));
    i++;
  }

  return g;
}

//replace an integer list with integers divided by their gcd
//(assumed memory already allocated)
void gcdlist(int *lstin, unsigned long num, int *lstout){
  unsigned long s;
  int g=gcd(lstin,num);
  for(s=0;s<num;s++)
    lstout[s]=lstin[s]/g;
  return;
}

// a binary number of length l in reverse order to a base 10 integer
// base10(...) is the version for binary strings in forward order
int binltoint(bool *vec, int l){
  int i,k=0;
  for(i=0;i<l;i++)
    k+=(int)pow(2.0, (double)i)*vec[i];
  return k;
}

// again in reverse; bool must be of sufficient length and already allocated
void inttobinl(int k1, bool *vec, int l){
  int i,k=k1;
  for(i=0;i<l;i++){
    vec[i]=(bool)(k%2);
    k=(k-(int)vec[i])/2;
  }
  return;
}

//versions which use unsigned long
long lbinltoint(bool *vec, int l){
  int i;
  unsigned long k=0;
  for(i=0;i<l;i++)
    k+=(int)pow(2.0, (double)i)*vec[i];
  return k;
}

//version which uses unsigned long
void linttobinl(unsigned long k1, bool *vec, int l){
  int i;
  unsigned long k=k1;
  for(i=0;i<l;i++){
    vec[i]=(bool)(k%2);
    k=(k-(int)vec[i])/2;
  }
  return;
}


// integer vector vec of length n has only even entries
bool iseven(int *vec, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec[i]%2!=0)
      return 0;
  }
  return 1;
}

// terminal character
int isend(char c){
  if((c=='\n') || (c==13) || (c=='\0'))
    return 1;
  return 0;
}

int isonlyspace(char *s){
  int c, k=0;
  while((c=s[k++]) != '\0'){
    if(!(isspace(c)))
      return 0;
  }
  return 1;
}

// decide if string s is an unsigned integer
int ispureint(char *s){
  int k=0;
  while(s[k] && isspace(s[k])) //skip space
    k++;
  if(!isdigit(s[k])) // starts with a noninteger or empty (apart from spaces)
    return 0;

  while(s[k] && isdigit(s[k]))
    k++;
  if(s[k] && !isspace(s[k]))
    return 0;
  while(s[k] && isspace(s[k])) //skip space
    k++;
  if(s[k]) // something more
    return 0;

  return 1;
}

//Are two vectors or lists equal
int areequal(int *vec1, int *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[i]!=vec1[i])
      return 0;
  }
  return 1;
}

//Overloading (boolean vectors)
int areequal(bool *vec1, bool *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[i]!=vec1[i])
      return 0;
  }
  return 1;
}

//Overloading: Are two float vectors equal upto some specified tolerancy
//here hard-coded as 10^{-8}
int areequal(double *vec1, double *vec2, int n){
  int i,tol=1e-8;
  for(i=0;i<n;i++){
    if(fabs(vec2[i]-vec1[i])>tol)
      return 0;
  }
  return 1;
}

//Is boolean vector vec1 a subvector of vec2? (both have length n)
bool issubvec(bool *vec1, bool *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec1[i]>vec2[i])
      return 0;
  }
  return 1;
}


//Overloading: check if vec1 is a subvector of any boolean vector in 
//the array vec2 from vec2[strt] to vec2[finish-1] 
//(indices assumed to be nonnegative and strt<=finish)
//return first match, and -1 if not
int issubvec(bool *vec1, bool **vec2, int strt, int finish, int n){
  int i;
  for(i=strt;i<finish;i++){
    if(issubvec(vec1, vec2[i],n))
      return i;
  }
  return -1;
}

// Is the ordered int list A a subset of the ordered int list B?
bool AsubsB(int *A, int szA, int *B, int szB){
  int lastfree=0;
  int i,j;
  bool flg;

  if(szA > szB)
    return 0;

  //  printvec(A,szA);  printvec(B,szB);
  for(i=0;i<szA;i++){
    if(lastfree>=szB)// off the end of superset
      return 0;
    flg=0;
    for(j=lastfree;j<szB;j++){
      if(B[j]>A[i])
	return 0;
      if(B[j]==A[i]){
	lastfree=j+1;
	flg=1;
	break;
      }
    }
    if(!flg)
      return 0;
  }
  return 1;
}


bool supervec(int *A, int szA, int **B, int totB){
  int i;
  for(i=0;i<totB;i++){
    //fprintf(stderr, "checking:\n");printvec(A,szA);printvec(B[i]+1,B[i][0]);
    if(AsubsB(B[i]+1, B[i][0], A, szA)){//Does the support of A contain that of any previously gathered?
      //fprintf(stderr, "superset:\n");printvec(A,szA);printvec(B[i]+1,B[i][0]);
      return 1;
    }
  }
  return 0;
}

void merge(int *xc, int *tmpvec, int k, int *outvec, int m){
  int i;
  inittozero(outvec,m);
  for(i=0;i<k;i++)
    outvec[xc[i]]=tmpvec[i];
  return;
}



void inittozero(int *vec, int len){
  int k;
  for(k=0;k<len;k++)
    vec[k]=0;
} 

void inittozero(bool *vec, int len){
  int k;
  for(k=0;k<len;k++)
    vec[k]=0;
} 


void inittozero(double *vec, int len){
  int k;
  for(k=0;k<len;k++)
    vec[k]=0.0;
} 

void inittozero(int **mat, int n, int m){
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<m;k++){
      mat[j][k]=0;
    }
  }
} 

void inittozero(bool **mat, int n, int m){
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<m;k++){
      mat[j][k]=0;
    }
  }
} 

void inittoone(int *vec, int len){
  int k;
  for(k=0;k<len;k++)
    vec[k]=1;
} 

void inittoone(bool *vec, int len){
  int k;
  for(k=0;k<len;k++)
    vec[k]=1;
} 

//return union of v and v1 in v
void boolunion(bool *v, bool *v1, int n){
  int i;
  for(i=0;i<n;i++){
    if(v1[i])
      v[i]=1;
  }
}

char *getlinefromstr(long *pos, char *block){
  int c=0, i=0;
  long len;
  char *str, *p;
  if(*pos > (int)strlen(block) || (*pos < 0)){
    str=strdup("");
    return str;
  }
  p=strchr(block+(*pos), '\n');
  if(p){
    len=strlen(block+(*pos)) - strlen(p);
    str = (char*) malloc(sizeof(char) * (len + 2));
    i=0;
    while((c=block[(*pos)++])!='\n')
      str[i++]=c;
    str[i++] = c;
    str[i] = 0;
  }
  else{
    str=strdup(block+(*pos));
    (*pos)=strlen(block);
  }
  return str;
}

int getnumints(char *s){
  int i, k;
  i=0, k=0;
  while(s[k]){
    while(s[k] && !isdigit(s[k])&& (s[k]!='-')&& (s[k]!='+')){k++;}
    if(isdigit(s[k])|| (s[k]=='-')|| (s[k]=='+')){
      i++;
      while(isdigit(s[k])|| (s[k]=='-')|| (s[k]=='+')){k++;}
    }
  }
  return i;
}

int getnumul(char *s){
  int i, k;
  i=0, k=0;
  while(s[k]){
    while(s[k] && !isdigit(s[k])){k++;}
    if(isdigit(s[k])){
      i++;
      while(isdigit(s[k])){k++;}
    }
  }
  return i;
}

int getnumints1(char *s){
  int i, k;
  i=0, k=0;
  while(s[k]){
    while(s[k] && !isdigit(s[k])&& (s[k]!='-')&& (s[k]!='+')){k++;}
    if(isdigit(s[k])|| (s[k]=='-')|| (s[k]=='+')){
      i++;
      while((s[k]=='-')||(s[k]=='+')){k++;}
      while(isdigit(s[k])){k++;}
    }
  }
  return i;
}

// number of segments in a line, separated by space
int getnumsegs(char *s){
  int i=0,k=0;
  while(s[k]){
    while(s[k] && isspace(s[k])){k++;} // skip initial space
    if(s[k] && !isspace(s[k])){
      i++;
      while(s[k] && !isspace(s[k])){k++;}//skip nonspace
    }
  }
  return i;
}


// chops a string (allocates memory for the output)
char *strchop2(char *t, int n, int n1){
  char *s;
  int i=0;
  if((n<= (int)strlen(t)) && n1>0){
    s = (char*) malloc(sizeof(char) * (n1+1));
    while (t[n+i] && (i< n1)){s[i] = t[n+i];i++;}
    s[i] = '\0';
  }
  else
    s = strdup("");
  
  return s;
}

//nth "word" from a string (any characters allowed)
char *getnthwd(char *s, int n){
  int i, j, k=0;
  char *v=NULL;
  while(s[k]){
    while(s[k] && isspace(s[k])){k++;} // skip space
    for(i=0;i<n-1;i++){
      while(s[k] && !isspace(s[k])){k++;} // skip word
      while(s[k] && isspace(s[k])){k++;} // skip nonword
    }
    j=0;
    while(!isspace(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}

//nth signed integer in a string
char *getnthint(char *s, int n){
  int i, j, k;
  i=0, k=0;
  char *v=NULL;
  while(s[k] != 0){
    while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    for(i=0;i<n-1;i++){
      while(isdigit(s[k]) || (s[k]=='-')|| (s[k]=='+')){k++;} // skip word
      while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    }
    j=0;
    while((s[k]=='-')||(s[k]=='+')||isdigit(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}

char *getnthint1(char *s, int n){
  int i=0, j, k=0;
  char *v=NULL;
  while(s[k]){
    while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    for(i=0;i<n-1;i++){
      while((s[k]=='-')||(s[k]=='+')){j++;k++;} //skip sign
      while(isdigit(s[k])){j++;k++;} // skip int
      while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    }
    j=0;
    while((s[k]=='-')||(s[k]=='+')){j++;k++;} // get the word
    while(isdigit(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}

// get the first integer in a string. Return pointer to just after
int getint(char **s){
  char *v;
  int j, k=0, r;

  while((*s)[k] && !isdigit((*s)[k]) && ((*s)[k]!='-')&& ((*s)[k]!='+')){k++;} // skip nonwords
  j=0;
  while(((*s)[k]=='-')||((*s)[k]=='+')){j++;k++;} // get the word
  while(isdigit((*s)[k])){j++;k++;} // get the word
  v = strchop2((*s), k-j, j);
  r=atoi(v);
  //cout << "line = \"" << *s << "\"" << endl;
  //cout << "v = " << v << ", k = " << k << endl;
  free(v);
  (*s)=(*s)+k;
  //cout << "newline = \"" << *s << "\"" << endl;
  return r;
}

//get an integer vector from a string. First element is to hold length
int *getintvec(char *str){
  int *out;
  int j;
  char *tmp=str;
  int len=getnumints(str);

  out=(int*) malloc(sizeof(int)*(len+1));
  out[0]=len;
  for(j=0;j<len;j++)
    out[j+1]=getint(&tmp);
  return out;
}

// get the first unsigned long integer in a string. Return pointer to just after
int getul(char **s){
  char *v;
  int j, k=0;
  unsigned long r;

  while((*s)[k] && !isdigit((*s)[k])){k++;} // skip nonwords including signs and newlines
  j=0;
  while(isdigit((*s)[k])){j++;k++;} // get the word
  v = strchop2((*s), k-j, j);
  r=strtoul(v,NULL,10);
  free(v);
  (*s)+=k;
  return r;
}

unsigned long *getulvec(char *str){
  unsigned long *out;
  int j;
  char *tmp=str;
  int len=getnumul(str);
  out=(unsigned long*) malloc(sizeof(unsigned long)*(len+1));
  out[0]=len;
  for(j=0;j<len;j++)
    out[j+1]=getul(&tmp);
  return out;
}


// is the integer "i" in the list lst?
int isinlist(int i, int ilst[], int tot){
  int j;
  for(j=0;j<tot;j++){
    if(i==ilst[j])
      return 1;
  }
  return 0;
}

//overloading
int isinlist(long i, long ilst[], int tot){
  int j;
  for(j=0;j<tot;j++){
    if(i==ilst[j])
      return 1;
  }
  return 0;
}

// a comment line or separator line or empty line
int iscomline(char s[]){ 
  int i=0;
  while((isspace((int) s[i]))){i++;}
  if (!s[i] || ((s[i] == '/') && (s[i+1] == '/')) || ((s[i] == '/') && (s[i+1] == '*')) || (s[i] == '#')){
    return 1;
  }
  return 0;
}

// str is a divider; return the parts of str before and after 
// the divider in str1 and str2
void chop(char *str, char **str1, char **str2, const char *str3){
  char *line;
  long pos=0;
  int pass=0;
  (*str1)=strdup(str);
  (*str2)=strdup(str);
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0 && !strstr(line, str3)){
    if(pass==0){strcpy((*str1), line);pass=1;}else{strcat((*str1), line);}
    free(line);line = getlinefromstr(&pos, str);
  }
  free(line);
  line = getlinefromstr(&pos, str);
  pass=0;
  if(strcmp(line, "")==0)
    strcpy((*str2), line);
  else{
    while(strcmp(line, "")!=0){
      if(pass==0){strcpy((*str2), line);pass=1;}else{strcat((*str2), line);}
      free(line);line = getlinefromstr(&pos, str);
    }
  }
  free(line);
}




//trim space from left and right (allocates memory)
char *lrtrim(char s[]){
  char *p;
  int n;
  while((*s) && isspace(*s)){s++;}// left trim
  p=strdup(s);
  for(n=strlen(p)-1;n>=0;n--)
    if(!(isspace(p[n])))
      break;
  p[n+1] = '\0';
  return p;
}

int startswith(char *s, char *lstr){
  char *p=NULL;
  char *tmp;
  tmp = lrtrim(s);
  p=strstr(tmp, lstr);
  if(!p || strcmp(p, tmp)!=0){
    free(tmp);
    return 0;
  }
  free(tmp);
  return 1;
}

//reverse string in place
void reverse(char s[]){
  int c, i, j;
  for (i=0, j=strlen(s)-1;i<j;i++,j--){
    c=s[i];
    s[i] = s[j];
    s[j] = c;
  }
}

//reverse integer list in place
void reverse(int s[], int tot){
  int i, j;
  for (i=0, j=tot-1;i<j;i++,j--)
    iswap(s, i, j);
}

int endswith(char *s, char *rstr){
  char *tmp, *tmp1;
  int j;

  tmp=lrtrim(s);
  tmp1=lrtrim(rstr);
  reverse(tmp);
  reverse(tmp1);
  j=startswith(tmp, tmp1);
  free(tmp);free(tmp1);
  return j;
}




// Get elements of string s separated by character 'sep' into the array v.
// Allocates memory for v and its elements
int chemgts2(char *s, char ***v, char sep){
  // sep is the separator
  int i, j, k;
  int numgets=0;
  char *tmp;
  i=0, k=0;

  (*v)=NULL;
  while(s[k]){
    j=0;
    while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
    while((s[k] != sep) && !isend(s[k])){j++;k++;}
    if(j>0)
      numgets++;
  }
  if(numgets>0){
    (*v)=(char**) malloc(sizeof(char*) * numgets);
    i=0;k=0;
    while(s[k]){
      j=0;
      while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
      while((s[k] != sep) && !isend(s[k])){j++;k++;}
      if(j>0){
	tmp=strchop2(s, k-j, j);
	(*v)[i++] = lrtrim(tmp);
	free(tmp);
      }
    }
  }
  return numgets;
}

// Add integer k to the end of integer vector v
// No memory checking -- assumes that v has enough space
int addtovec(int *v, int *numv, int k){
  int i;
  for(i=0;i<(*numv);i++){
    if(v[i]==k)//already there
      return i;
  }
  v[*numv]=k;
  (*numv)++;
  return (*numv)-1;
}


// check if a string s already belongs to a list of strings t; and if 
// not to add it to the end of t. If the string is found put its position
// in ind. Returns new free position in the list.
int addv(int k, const char *s, char ***t, int *ind){
  int i=0;

  if(k==0){
    (*t) = (char**) malloc(sizeof(char*) * 1);
    (*t)[0] = strdup(s);
    (*ind)=0;
    return 1;
  }

  while(i<k){
    if (strcmp(s, (*t)[i])==0){(*ind)=i;return k;}
    i++;
  }

  (*t)=(char**) realloc((*t), sizeof(char*) *(k+1));
  (*t)[k] = strdup(s);
  (*ind)=k;
  return k+1;
}

// The version of addv where we don't need the position of s in the list
// even if s is found
int addv1(int k, const char *s, char ***t){
  int i=0;

  if(k==0){
    (*t) = (char**) malloc(sizeof(char*) * 1);
    (*t)[k] = strdup(s);
    return k+1;
  }

  while(i<k)
    if (strcmp(s, (*t)[i++]) == 0){return k;}

  (*t)=(char**) realloc((*t), sizeof(char*) *(k+1));
  (*t)[k] = strdup(s);
  return k+1;
}

// The version of addv where we assume s is absent from t (or don't care)
int addv2(int k, char *s, char ***t){
  if(k==0){
    (*t) = (char**) malloc(sizeof(char*) * 1);
    (*t)[k] = strdup(s);
    return k+1;
  }

  (*t)=(char**) realloc((*t), sizeof(char*) *(k+1));
  (*t)[k] = strdup(s);
  return k+1;
}

//Add to an array of integer vectors (even if already present)
int addnewtoarray(int ***outmat, int sz, int *xc, int len){
  int l;
  if(sz==0) // first entry
    (*outmat) = (int**) malloc(sizeof(int*) * 1);
  else
    (*outmat)=(int**) realloc((*outmat), sizeof(int*) *(sz+1));
  (*outmat)[sz]=(int *)malloc(sizeof(int) *(len));
  for(l=0;l<len;l++)
    (*outmat)[sz][l]=xc[l];
  return sz+1;
}

//overloading
int addnewtoarray(long ***outmat, int sz, long *xc, int len){
  int l;
  if(sz==0) // first entry
    (*outmat) = (long**) malloc(sizeof(long*) * 1);
  else
    (*outmat)=(long**) realloc((*outmat), sizeof(long*) *(sz+1));
  (*outmat)[sz]=(long *)malloc(sizeof(long) *(len));
  for(l=0;l<len;l++)
    (*outmat)[sz][l]=xc[l];
  return sz+1;
}

// Add to a list of integers
void addnewto1Darray(int **outvec, int sz, int newint){
  if(sz==0) // first entry
    (*outvec) = (int*) malloc(sizeof(int) * 1);
  else
    (*outvec)=(int*) realloc((*outvec), sizeof(int) *(sz+1));
  (*outvec)[sz]=newint;
  return;
}

// Add to a bool list
void addnewto1Darray(bool **outvec, int sz, bool newbool){
  if(sz==0) // first entry
    (*outvec) = (bool*) malloc(sizeof(bool) * 1);
  else
    (*outvec)=(bool*) realloc((*outvec), sizeof(bool) *(sz+1));
  (*outvec)[sz]=newbool;
  return;
}

//overloading (doubles)
void addnewto1Darray(double **outvec, int sz, double newint){
  if(sz==0) // first entry
    (*outvec) = (double*) malloc(sizeof(int) * 1);
  else
    (*outvec)=(double*) realloc((*outvec), sizeof(double) *(sz+1));
  (*outvec)[sz]=newint;
  return;
}

//overloading (long)
void addnewto1Dlarray(long **outvec, int sz, long newint){
  if(sz==0) // first entry
    (*outvec) = (long *) malloc(sizeof(long) * 1);
  else
    (*outvec)=(long *) realloc((*outvec), sizeof(long) *(sz+1));
  (*outvec)[sz]=newint;
  return;
}

/* routine to grow a list of integer strings. Returns new free position in the list. */
int addvec1(int lstlen, int *newvec, int veclen, int ***veclst, int check, int *pos){
  int i=0;
  (*pos)=-1;

  if(lstlen==0)
    (*veclst) = (int**) malloc(sizeof(int*) * 1);
  //either allow replacement (check=0) or not already present
  else if(!check || ((*pos)=isinarray((*veclst), lstlen, newvec, veclen))<0)
    (*veclst)=(int**) realloc((*veclst), sizeof(int*) *(lstlen+1));
  else
    return lstlen;

  (*veclst)[lstlen] = (int*) malloc(sizeof(int) * (veclen));
  for(i=0;i<veclen;i++)
    (*veclst)[lstlen][i] = newvec[i];
  return lstlen+1;
}

//overloading for a boolean vector
int addvec1(int lstlen, bool *newvec, int veclen, bool ***veclst, int check, int *pos){
  int i=0;
  (*pos)=-1;

  if(lstlen==0)
    (*veclst) = (bool**) malloc(sizeof(bool*) * 1);
  //either allow replacement (check=0) or not already present
  else if(!check || ((*pos)=isinarray((*veclst), lstlen, newvec, veclen))<0)
    (*veclst)=(bool**) realloc((*veclst), sizeof(bool*) *(lstlen+1));
  else
    return lstlen;

  (*veclst)[lstlen] = (bool*) malloc(sizeof(bool) * (veclen));
  for(i=0;i<veclen;i++)
    (*veclst)[lstlen][i] = newvec[i];
  return lstlen+1;
}

//overloading for a double vector
int addvec1(int lstlen, double *newvec, int veclen, double ***veclst, int check, int *pos){
  int i=0;
  (*pos)=-1;

  if(lstlen==0)
    (*veclst) = (double**) malloc(sizeof(double*) * 1);
  //either allow replacement (check=0) or not already present
  else if(!check || ((*pos)=isinarray((*veclst), lstlen, newvec, veclen))<0)
    (*veclst)=(double**) realloc((*veclst), sizeof(double*) *(lstlen+1));
  else
    return lstlen;

  (*veclst)[lstlen] = (double*) malloc(sizeof(double) * (veclen));
  for(i=0;i<veclen;i++)
    (*veclst)[lstlen][i] = newvec[i];
  return lstlen+1;

}

// is string s in string-array v? If so return its index, if not return -1
int isinarray(char *v[], int numv, const char *s){
  int i;
  for(i=0;i<numv;i++)
    if(strcmp(s, v[i]) == 0)
      return i;
  return -1;
}

// overloading: ul array: since there is no negative return, 
// **we offset the return value**
unsigned long isinarray(char *v[], unsigned long numv, const char *s){
  unsigned long i;
  for(i=0;i<numv;i++)
    if(!strcmp(s,v[i]))
      return i+1;
  return 0;
}

// overloading: is the integer vector s a member of the array v? 
long isinarray(int **v, long numv, int *s, int n){
  long i;
  for(i=0;i<numv;i++){
    if(areequal(v[i], s, n))
      return i;
  }
  return -1;
}

// Overloading
int isinarray(int **v, int numv, int *s, int n){
  int i;
  for(i=0;i<numv;i++){
    if(areequal(v[i], s, n))
      return i;
  }
  return -1;
}

// Overloading (bool)
int isinarray(bool **v, int numv, bool *s, int n){
  int i;
  for(i=0;i<numv;i++){
    if(areequal(v[i], s, n))
      return i;
  }
  return -1;
}

//overloading (double)
int isinarray(double **v, int numv, double *s, int n){
  int i;
  for(i=0;i<numv;i++){
    if(areequal(v[i], s, n))
      return i;
  }
  return -1;
}



//Is the integer vector s a member of the array v? 
//If so return its index, if not return -1
// assume first two entries are other data
int isinarrayoffset2(int **v, int numv, int *s, int n){
  int i;
  for(i=0;i<numv;i++){
    if(areequal(v[i]+2, s, n))
      return i;
  }
  return -1;
}



void veccp(int *cflst,long r,int *cflst1){
  long s;
  for(s=0;s<r;s++)
    cflst1[s]=cflst[s];
  return;
}

void veccp(int *vecin,int r,int *vecout){
  int s;
  for(s=0;s<r;s++)
    vecout[s]=vecin[s];
  return;
}

void veccp(bool *vecin,int r,bool *vecout){
  int s;
  for(s=0;s<r;s++)
    vecout[s]=vecin[s];
  return;
}

void veccp(double *vecin,int r,double *vecout){
  int s;
  for(s=0;s<r;s++)
    vecout[s]=vecin[s];
  return;
}

//memory allocating version of veccp
int *veccp(int *vecin,int r){
  int s;
  int *vec1=(int *)malloc((size_t) ((r)*sizeof(int)));
  for(s=0;s<r;s++)
    vec1[s]=vecin[s];
  return vec1;
}


// add b to a and return in a
void vecadd(int *a, int *b, int n){
  int i;
  for(i=0;i<n;i++)
    a[i]+=b[i];
  return;
}

//overloading
void vecadd(double *a, double *b, int n){
  int i;
  for(i=0;i<n;i++)
    a[i]+=b[i];
  return;
}

// add nXm matrix b to nXm matrix a and return in a
void matadd(int **a, int **b, int n, int m){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      a[i][j]+=b[i][j];
  return;
}

// minus an nXm matrix in place
void minusmat(int **a, int n, int m){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      a[i][j]=-a[i][j];
  return;
}

//nXn identity matrix
int **idmat(int n, bool minus){
  int i,j;
  int **a=imatrix(0,n-1,0,n-1);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i==j && !minus)
	a[i][j]=1;
      else if(i==j && minus)
	a[i][j]=-1;
      else
	a[i][j]=0;
    }
  }
  return a;
}

// minus an nXm symbolic matrix and return new
matrix minusmat(matrix Jin, int n, int m){
  int i,j;
  matrix Jout(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      Jout(i,j)=-Jin(i,j);
  }
  return Jout;
}



void vecsum2(int *out, int *a, int *b, int numv){
  inittozero(out,numv);
  vecadd(out,a,numv);
  vecadd(out,b,numv);
  return;
}

void vecsum3(int *out, int *a, int *b, int *c, int numv){
  inittozero(out,numv);
  vecadd(out,a,numv);
  vecadd(out,b,numv);
  vecadd(out,c,numv);
  return;
}


int *subtract(int *a, int *b, int numv){
  int i;
  int *out=(int *)malloc((size_t) (numv*sizeof(int)));
  for(i=0;i<numv;i++)
    out[i]=a[i]-b[i];
  return out;
}

//overloading
int **subtract(int **a, int **b, int n, int m){
  int i,j;
  int **out=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      out[i][j]=a[i][j]-b[i][j];

  return out;
}

int *halve(int *a, int n){
  int i;
  int *out=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++)
    out[i]=a[i]/2;
  return out;
}




// is vec2 > vec1 (backwards)
int isgrtr(int *vec1, int *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[n-1-i]>vec1[n-1-i])
      return 1;
    else if(vec2[n-1-i]<vec1[n-1-i])
      return -1;
  }
  return 0;
}

// is vec2 > vec1 (forwards)
int isgrtrf(int *vec1, int *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[i]>vec1[i])
      return 1;
    else if(vec2[i]<vec1[i])
      return -1;
  }
  return 0;
}


//gcd of entries in a vector

int gcd_vec(int *v, int n){
  int i,j=0,g;
  if(n==1){
    if(v[0])
      return abs(v[0]);
    else
      return 1;
  }
  while(!(v[j]))
    j++;
  if(j>=n)//zero vector
    return 1;
  g=abs(v[j]);
  for(i=j+1;i<n;i++)
    if(v[i])
      g=gcd(g,abs(v[i]));
  return g;
}

//replace an integer vector with the smallest collinear integer vector
int reduce_vec(int *v, int n){
  int i;
  int g=gcd_vec(v, n);
  for(i=0;i<n;i++)
    v[i]/=g;
  return g;
}

// are columns indexed by xc in the matrix A pairwise disjoint? 
// (i.e. no pair have nonzero entries in the same place)
bool disjointcols(int **A, int n, int m, int *xc, int len){
  int i,j,flag;
  for(i=0;i<n;i++){
    j=0;flag=0;
    while(j<len && flag<2){
      if(A[i][xc[j]])
	flag++;
      j++;
    }
    if(flag>=2)
      return 0;
  }
  return 1;
}

//overloading
bool disjointcols(int **A, int n, int m, int i, int j){
  int xc[2];
  xc[0]=i;xc[1]=j;
  return disjointcols(A, n, m, xc, 2);
}

//are columns indexed by xc sign compatible
//minus means no two must have the same nonzero sign
bool signcompatcols(int **A, int n, int m, int *xc, int len, bool minus){
  int i,j,first;
  for(i=0;i<n;i++){
    j=0;first=0;
    while(j<len){
      if(!first && A[i][xc[j]])//first nonzero
	first=A[i][xc[j]];
      else{
	if(!minus && first*A[i][xc[j]]<0)
	  return 0;
	if(minus && first*A[i][xc[j]]>0)
	  return 0;
      }
      j++;
    }
  }
  return 1;
}

//overloading
bool signcompatcols(int **A, int n, int m, int i, int j, bool minus){
  int xc[2];
  xc[0]=i;xc[1]=j;
  return signcompatcols(A, n, m, xc, 2, minus);
}


bool isreversecol(int **S, int n, int j, int k){
  int i;
  for(i=0;i<n;i++){
    if(S[i][j]!=-S[i][k])
      return 0;
  }
  return 1;
}

bool isequalcol(int **S, int n, int j, int k){
  int i;
  for(i=0;i<n;i++){
    if(S[i][j]!=S[i][k])
      return 0;
  }
  return 1;
}

//replace a matrix with a 0,1,-1,matrix
void sgn_mat(int **imat, int n, int m){
  int i, j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat[i][j])
	imat[i][j]/=abs(imat[i][j]);
    }
  }
}

// replace an integer matrix with the matrix each of whose rows is the
// smallest integer vector collinear with the original
// transpose means do columns
bool reduce_mat(int **imat, int n, int m, bool transpose){
  int i, **imatT,g;
  bool flag=0;
  if(transpose){
    imatT=transposemat(imat, n, m);
    for(i=0;i<m;i++){
      g=reduce_vec(imatT[i],n);
      if(g>1)
	flag=1;
    }
    transposemat(imatT, m, n, imat);
    free_imatrix(imatT,0,m-1,0,n-1);
  }
  else{
    for(i=0;i<n;i++){
      g=reduce_vec(imat[i],m);
      if(g>1)
	flag=1;
    }
  }

  return flag;
}



void growlist(int **lst, int *r){
  (*r)++;
  if((*r)==1)
    (*lst)=(int*)malloc((size_t) (sizeof(int)));
  else
    (*lst)=(int*) realloc((*lst), sizeof(int)*((*r)));
}

//boolean verson
void growlist(bool **lst, int *r){
  (*r)++;
  if((*r)==1)
    (*lst)=(bool*)malloc((size_t) (sizeof(bool)));
  else
    (*lst)=(bool*) realloc((*lst), sizeof(bool)*((*r)));
}

// modified from the original
// generate an array from a data file. Memory for each entry 
// in the datafile is allocated. 
char **genarraydat(const char *datfile, int *num){
  int linelen, i=0;
  char *line;
  char **vararray;
  unsigned long numl=0;
  int maxl=0;
  FILE *fd;

  *num=0; // initialise
  numl = numflines(datfile, &maxl); // array size and maximum line length
  if(numl <=0){return NULL;}

  fd = fopen(datfile, "r");
  if(!fd){
    fprintf(stderr, "ERROR in genarraydat - File: %s could not be opened...\n", datfile);
    return NULL;
  }
  vararray=charar(0, numl-1);
  line = (char*) malloc(sizeof(char) * (maxl));

  while((linelen = gtline(fd, line, maxl)) > 0){
      vararray[i++] = getnthwd(line, 1);
  }
  fclose(fd);
  (*num)=i;
  free(line);
  return vararray;
}

//overloading: number of lines is ul, and skip comment lines
char **genarraydat(const char *datfile, unsigned long *num){
  int linelen;
  char *line;
  char **vararray;
  unsigned long numl=0;
  int maxl=0;
  FILE *fd;

  (*num)=0; // initialise
  numl = numflines(datfile, &maxl); // array size and maximum line length
  if(numl <=0){return NULL;}

  fd = fopen(datfile, "r");
  if(!fd){
    fprintf(stderr, "ERROR in genarraydat - File: %s could not be opened...\n", datfile);
    return NULL;
  }
  vararray=charar(0, numl-1);
  line = (char*) malloc(sizeof(char) * (maxl));

  while((linelen = gtline(fd, line, maxl)) > 0 && !iscomline(line)){
      vararray[(*num)++] = getnthwd(line, 1);
  }
  fclose(fd);
  free(line);
  return vararray;
}



//
//
// Files
//
//

//
// get the number of non-empty lines and 
// maximum linelength in a file
//

unsigned long numflines(const char *fname, int *maxline){
  FILE *fd;
  int i=0, c=0;
  unsigned long numl=0;
  (*maxline)=0;

  if(!(fd= fopen(fname, "r"))){
    fprintf(stderr, "WARNING in numflines: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }

  while((c=getc(fd)) != EOF){
    i++;
    if (c=='\n'){
      numl++;
      (*maxline)=max(i, (*maxline));
      i=0;
    }
  }
  if(i!=0){ // final line does not end in newline
    numl++;
    (*maxline)=max(i, (*maxline));
  }
  (*maxline)+=2; // to hold the terminal character plus possible EOF

  fclose(fd);
  return numl;
}

int gtline(FILE *fp, char s[], int lim)
{
  /* store a line as a string, including the terminal newline character */
  int c=0, i;

  for(i=0; i<lim-1 && (c=getc(fp))!=EOF && c!='\n';++i){s[i] = c;}
  if (c == '\n') {s[i++] = c;}
  s[i] = '\0';
  if(i==lim-1 && c!='\n' && c!=0)
    fprintf(stderr, "WARNING in gtline: line length exceeded. *%s*\n", s);
  return i;
}


bool filehasstr(char *fname, char *str){
  FILE *fd;
  char oneline[1000];

  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "Error in filehasstr: file %s couldn't be opened. EXITING\n", fname);
    exit(0);
  }

  while(getline0(fd, oneline, 1000) > 0){
    if(strstr(oneline, str)){
      fclose(fd);
      return 1;
    }
  }
  fclose(fd);
  return 0;
}


// get the filename minus any directory names, and any extension

char *getfnameraw(char *fname){
  int i, len;
  char *fnametrue, *tmp;
  len=strlen(fname);
  i=len;
  while(i>=0 && (fname[i] != '/'))
    i--;
  i++;
  tmp=strchop2(fname, i, len-i);
  i=0;
  while(tmp[i] && tmp[i] != '.')
    i++;
  fnametrue = strchop2(tmp, 0, i);
  free(tmp);
  return fnametrue;
}

char *getdname(char *fname){
  int i, len;
  char *dname;
  len=strlen(fname);
  i=len;
  while(i>=0 && (fname[i] != '/'))
    i--;
  i++;
  dname=strchop2(fname, 0, i);
  return dname;
}

char *readfileintostr(const char fname[]){
  FILE *fd;
  unsigned long i=0;
  int c=0;
  char *str;

  // first parse - find length
  if(!(fd=fopen(fname, "r"))){
    fprintf(stderr, "ERROR in readfileintostr: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }

  while((c=getc(fd)) != EOF)
    i++; 
  i++;
  str = (char*) malloc(sizeof(char) * (i));

  // second parse - get number of non-comment lines
  rewind(fd);i=0;
  while((c=getc(fd)) != EOF)
    str[i++]=c;
  str[i++]=0;
 
  fclose(fd);
  return str;
}

//no terminal character
int getline0(FILE *fp, char s[], int lim){
  int c=0, i;

  for(i=0; i<lim-1 && (c=getc(fp))!=EOF && c!=13 && c!='\n';++i)
    s[i] = c;
  s[i] = '\0';
  if(i==lim-1 && c!='\n' && c!=EOF)
    fprintf(stderr, "WARNING in getline0: line length exceeded. %s\n", s);
  if (c == '\n' || c == 13) { // don't do anything
    return i+1;
  }
  else
    return i;
  // we return the length it would have if it had a newline character at the end
}

//
//
// permutations and combinations
//
//

unsigned long factorial(int x){
  int i;
  unsigned long factx = 1;
  for(i=1; i<=x ; i++ )
    factx *= i;
  return factx;
}

int choose(int n, int k){
  if(n<k)
    return 0;
  if(!k)
    return 1;
  return (n*choose(n-1, k-1))/k;
}

//nCk. not pretty, but works to reasonably large values
unsigned long comb(int n, int k){
  int i,j,ind1=2;
  unsigned long combx=1;
  int inds[n-k];
  if(k>n)
    return 0;
  for(i=0;i<n-k;i++)
    inds[i]=1;
  j=ind1;
  for(i=k+1; i<=n ; i++ ){//fprintf(stderr, "%d*%d\t", combx, i);
    combx *= i;
    j=ind1;
    while(j<=n-k){ // try to remove factors
      if(combx%j==0 && inds[j-1]==1 ){
	combx/=j; inds[j-1]=0; if(j==ind1){ind1++;}
      }
      j++;
    }
  }
  return combx;
}

//vector vec of length n1 holds initial combination on n integers, 
//n obviously assumed to be >= n1
void firstcomb(int *vec, int n, int n1){
  int i;
  if(n1>n)
    inittozero(vec,n1);
  else{
    for(i=0;i<n1;i++)
      vec[i]=i;
  }
  return;
}

//increment the combination of length n1 on n integers
int nextcomb(int *vec, int n, int n1){
  int i, j;
  for(i=0;i<n1;i++){
    if(vec[n1-1-i]< n-1-i){
      vec[n1-1-i]++;
      for(j=n1-i;j<n1;j++)
	vec[j]=vec[j-1]+1;
      return 1;
    }
  }
  // only get here if we fail
  firstcomb(vec, n, n1);
  return 0;
}

//skip forward k combinations
void nextcombk(int *vec, int n, int n1, long k){
  long i;
  for(i=0;i<k;i++)
    nextcomb(vec, n, n1);
  return;
}

//Take a combination and translate to stars and bars
//We interpret entries in vec as positions of the "bars"
//We output a vector of the number of "stars" between bars
//vecout of length n1+1
int starbar(int *vec, int n, int n1, int *vecout){
  int i,tot=0;
  vecout[0]=vec[0];tot+=vecout[0];
  for(i=1;i<n1;i++){
    vecout[i] = vec[i]-vec[i-1]-1;
    tot+=vecout[i];
  }
  vecout[n1]=n-n1-tot;
  return 0;
}

// return all permutations of [0...n-1]
// should be okay for n upto about 10
// No error checking
int **allperms(int n){
  int **tmp;
  int **tmpmat;
  int i,j,k;
  unsigned long r;
  tmp=imatrix(0, factorial(n)-1, 0, n);
  if(n==1){
    tmp[0][0]=0;
    tmp[0][1]=1; // parity
  }
  else{
    tmpmat=allperms(n-1);
    i=0;
    for(k=0;k<n;k++){
      for(r=0;r<factorial(n-1);r++){
	i=n*r+k;
	for(j=0;j<n-k-1;j++)
	  tmp[i][j]=tmpmat[r][j];
	tmp[i][n-k-1]=n-1;
	for(j=n-k;j<n;j++)
	  tmp[i][j]=tmpmat[r][j-1];
	tmp[i][j]=tmpmat[r][j-1]*par(k); // parity
      }
    }
    free_imatrix(tmpmat, 0,factorial(n-1)-1,0,n-1);
  }
  return tmp;
}


/* recursive routine to return all permutations of an arbitrary  */
/* n-vector [vec1...vecn]. A final digit represents the parity  */
/* of the permutation: 1 for even, and -1 for odd */
int **allperms1(int *vec, int n){
  int **tmp;
  int **tmpmat;
  int i,j,k;
  unsigned long r;
  tmp=imatrix(0, factorial(n)-1, 0, n);
  if(n==1){
    tmp[0][0]=vec[0];
    tmp[0][1]=1; // parity
  }
  else{
    tmpmat=allperms1(vec, n-1);
    i=0;
    for(k=0;k<n;k++){
      for(r=0;r<factorial(n-1);r++){
	i=n*r+k;
	for(j=0;j<n-k-1;j++)
	  tmp[i][j]=tmpmat[r][j];
	tmp[i][n-k-1]=vec[n-1];
	for(j=n-k;j<n;j++)
	  tmp[i][j]=tmpmat[r][j-1];
	tmp[i][j]=tmpmat[r][j-1]*par(k);  //parity
      }
    }
    free_imatrix(tmpmat, 0,factorial(n-1)-1,0,n-1);
  }
  return tmp;
}


// generates all combinations of [0...n-1] of length n1, ordered (recursive)
int **allcombsgen(int n, int n1){
  int **tmp;
  int **tmpmat;
  int i,j,k;
  long cn, r;
  tmp=imatrix(0, comb(n, n1)-1, 0, n1-1);
  if(n1==1){
    for(k=0;k<n;k++)
      tmp[k][0]=k;
  }
  else{
    tmpmat=allcombsgen(n, n1-1);
    i=0;
    cn=comb(n,n1-1);
    for(r=0;r<cn-1;r++){
      if(tmpmat[r][n1-2]<n-1){
	for(k=tmpmat[r][n1-2]+1;k<n;k++){
	  for(j=0;j<n1-1;j++)
	    tmp[i][j]=tmpmat[r][j];
	  tmp[i][j]=k;
	  i++;
	}
      }
    }
    free_imatrix(tmpmat, 0,cn-1,0,n1-2);
  }
  return tmp;
}

// Why int ** (change to int *)
// My implementation of Johnson-Trotter algorithm
int nextperm(int **vec, int **veclr, int *par, int n){
  int tmp=-1, tmpi=-1, tmplr;
  int i, flag=0;
  for(i=0;i<n;i++){
    if((*vec)[i]>tmp){
      if((*veclr)[i]==1 && i<n-1 && ((*vec)[i] > (*vec)[i+1])){ // right
	tmp = (*vec)[i];
	tmpi=i;
	tmplr=(*veclr)[i];
	flag=1;
      }
      else if ((*veclr)[i]==-1 && i>0 && ((*vec)[i] > (*vec)[i-1])){ // left
	tmp = (*vec)[i];
	tmpi=i;
	tmplr=(*veclr)[i];
	flag=2;
      }
    }
  }
  if(flag==1){
    (*vec)[tmpi]=(*vec)[tmpi+1];
    (*vec)[tmpi+1]=tmp;
    (*veclr)[tmpi]=(*veclr)[tmpi+1];
    (*veclr)[tmpi+1]=tmplr;
    *par = -(*par);
    for(i=0;i<n;i++){
      if((*vec)[i]>tmp)
	(*veclr)[i]=-(*veclr)[i];
    }
    return 1;
  }
  if(flag==2){
    (*vec)[tmpi]=(*vec)[tmpi-1];
    (*vec)[tmpi-1]=tmp;
    (*veclr)[tmpi]=(*veclr)[tmpi-1];
    (*veclr)[tmpi-1]=tmplr;
    *par = -(*par);
    for(i=0;i<n;i++){
      if((*vec)[i]>tmp)
	(*veclr)[i]=-(*veclr)[i];
    }
    return 1;
  }
  return 0;

}



/* Just for testing purposes: prints out the */
/* permutations on [1, ..., n] in transposition order */

int tr(){
  int *vec;
  int *veclr;
  int i, par=1, flag=1,t=0;
  int n=4;
  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<n;i++){
    vec[i]=i;
    veclr[i]=-1;
  }

  while(flag){
    for(i=0;i<n;i++){
      if(veclr[i]==-1)
	fprintf(stderr, "<%d  ", vec[i]+1);
      else
	fprintf(stderr, "%d>  ", vec[i]+1);
    }
    fprintf(stderr, "(%d)\n", par);
    flag=nextperm(&vec, &veclr, &par, n);
    t++;
  }
  return 0;

}


//a sequence i_init, ..., i_init+n1
int firstcombfrom(int *vec, int n, int n1, int i_init){
  int i;
  //fprintf(stderr, "n=%d, n1=%d, i_init = %d\n", n, n1, i_init);
  if(n<i_init+n1){
    for(i=0;i<n1;i++)
      vec[i]=0;
    //fprintf(stderr, "exiting here\n");
    return 0;
  }
  else{
    for(i=i_init;i<i_init+n1;i++)
      vec[i-i_init]=i;
  }
  //fprintf(stderr, "or here\n");
  return 1;
}


//The next number in a given base
int nextnum(int *vec, int n, int base){
  int i;
  for(i=0;i<n;i++){
    if(vec[n-1-i]< base-1){
      vec[n-1-i]++;
      return 1;
    }
    vec[n-1-i]=0;
  }
  return 0;
}

// each digit has a different base
int nextnumb(int *vec, int n, int *base){
  int i;
  for(i=0;i<n;i++){
    if(vec[n-1-i]< base[n-1-i]-1){
      vec[n-1-i]++;
      return 1;
    }
    vec[n-1-i]=0;
  }
  //couldn't increment
  return 0;
}


// a number of length n in some base to base 10
long base10(int *vec, int n, int base){
  long j=0;
  int i=0;
  for(i=0;i<n;i++)
    j+=vec[i]*((int)(pow((double)base,(double)(n-i-1))));
  return j;
}

// a number in base 10 to a number of length n in some other base
void basek(int *vec, int n, int base, long num){
  int i;
  for(i=0;i<n;i++){
    vec[i]=num/((int)(pow((double)base,(double)(n-i-1))));
    num-=vec[i]*((int)(pow((double)base,(double)(n-i-1))));
    if(i==0)
      vec[i]=vec[i]%base;
  }
  return;
}

// quick version of adding in base "base", k is in base 10.
void nextnumk(int *vec, int n, int base, long k){
  long tmp;
  tmp=base10(vec, n, base);
  tmp+=k;
  basek(vec, n, base, tmp);
}



//
//
// print output
//
//

//output to stderr/cerr
void printvec(int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%d ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printvec(double *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%.5f ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}


void printvec(bool *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%d ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printvec(long *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%ld ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printvec(unsigned long *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%ld ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

//output to stdout/cout
void printvec1(int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    printf("%d ", ivec[i]);
  printf("\n");
  return;
}

void printvec1(double *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    printf("%.5f ", ivec[i]);
  printf("\n");
  return;
}


void printvec1(bool *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    printf("%d ", ivec[i]);
  printf("\n");
  return;
}

void printvec1(long *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    printf("%ld ", ivec[i]);
  printf("\n");
  return;
}


//output to stderr/cerr
void printmat(int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%2d ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

//output to stderr/cerr
void printmat(double **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%.4f ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printmat(bool **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%2d ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printmat(matrix imat, int n, int m){
  int i,j;
  //cerr << "symbolic matrix:\n";
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      cerr << imat(i, j) << " ";
    cerr << endl;
  }
  cerr << endl;
  return;
}

void printmat(char **cmat, int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    fprintf(stderr, "%s\t", cmat[i]);
    for(j=0;j<m;j++){
      fprintf(stderr, "%2d  ", imat[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

// Output to stdout/cout
void printmat1(int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      printf("%2d ", imat[i][j]);
    printf("\n");
  }
  printf("\n");
  return;
}

void printmat1(bool **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      printf("%2d ", imat[i][j]);
    printf("\n");
  }
  printf("\n");
  return;
}

void printmat1(matrix imat, int n, int m){
  int i,j;
  //cout << "symbolic matrix:\n";
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      cout << imat(i, j) << " ";
    cout << endl;
  }
  cout << endl;
  return;
}

void printmat1(char **cmat, int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    printf("%s\t", cmat[i]);
    for(j=0;j<m;j++){
      printf("%2d  ", imat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  return;
}

//print the matrix in a form required for maxima input
void printmaximamat(int **imat, int n, int m){
  int i,j;
  fprintf(stderr, "SS:matrix(");
  for(i=0;i<n;i++){
    fprintf(stderr, "[");
    for(j=0;j<m;j++){
      fprintf(stderr, "%d",imat[i][j]);
      if(j<(m-1))
	fprintf(stderr, ",");
    }
   if(i<n-1){fprintf(stderr, "],");}else{fprintf(stderr, "]");}
  }
  fprintf(stderr, ");\n");
  return;
}

//overloading: take a symbolic matrix as input
void printmaximamat(matrix mat, int n, int m){
  int i,j;
  cerr << "SS:matrix(";
  for(i=0;i<n;i++){
    cerr<< "[";
    for(j=0;j<m;j++){
      cerr<< mat[i][j];
      if(j<(m-1))
	cerr << ",";
    }
   if(i<n-1){cerr << "],";}else{cerr<< "]";}
  }
  cerr<< ");\n";
  return;
}

//print the matrix in a form required for maxima input
//This time to stdout
void printmaximamat1(int **imat, int n, int m){
  int i,j;
  printf("SS:matrix(");
  for(i=0;i<n;i++){
    printf("[");
    for(j=0;j<m;j++){
      printf("%d",imat[i][j]);
      if(j<(m-1))
	printf(",");
    }
   if(i<n-1){printf("],");}else{printf("]");}
  }
  printf(");\n");
  return;
}

//overloading: take a symbolic matrix as input
//print to cout
void printmaximamat1(matrix mat, int n, int m){
  int i,j;
  cout << "SS:matrix(";
  for(i=0;i<n;i++){
    cout<< "[";
    for(j=0;j<m;j++){
      cout<< mat[i][j];
      if(j<(m-1))
	cout << ",";
    }
   if(i<n-1){cout << "],";}else{cout<< "]";}
  }
  cout<< ");\n";
  return;
}

void printtens(int ***tens, int n, int m){
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(!tens[i][j][0])
	fprintf(stderr, "%2d ", 0);
      else{
	fprintf(stderr, " [");
	for(k=1;k<=tens[i][j][0];k++)
	  fprintf(stderr, "%2d ", tens[i][j][k]);
	fprintf(stderr, "] ");
      }
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}


void fprintvec(FILE *fd, int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(fd, "%d ", ivec[i]);
  fprintf(fd, "\n");
  return;
}

//no space
void fprintveca(FILE *fd, int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(fd, "%d", ivec[i]);
  fprintf(fd, "\n");
  return;
}




void printsubmat(int **imat, int *vec1, int *vec2, int k1, int k2){
  int i,j,i1,j1;
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(vec1)
	i1=vec1[i];
      else
	i1=i;
      if(vec2)
	j1=vec2[j];
      else
	j1=j;
      fprintf(stderr, "%2d ",imat[i1][j1]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printsubmat(matrix imat, int *vec1, int *vec2, int n, int m){
  int i,j;
  cerr << "symbolic matrix:\n";
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      cerr << imat(vec1[i], vec2[j]) << "\t";
    cerr << endl;
  }
  cerr << endl;
  return;
}

void printindsubmat(int **imat, int *vec1, int *vec2, int k1, int k2){
  int i,j;
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0)
	fprintf(stderr, "0  ");
      else if(imat[vec1[i]][vec2[j]]<0)
	fprintf(stderr, "-v%d ",-imat[vec1[i]][vec2[j]]);
      else
	fprintf(stderr, "v%d ",imat[vec1[i]][vec2[j]]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}



void printmatpat(char **cmat, int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    fprintf(stderr, "%s\t", cmat[i]);
    for(j=0;j<m;j++){
      if(imat[i][j]==2)
	fprintf(stderr, " x  ");
      else
	fprintf(stderr, "%2d  ", imat[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void listprint(int *cfs, int **mons, int dim, long len){
  long k;
  for(k=0;k<len;k++){
    cout << cfs[k] << ", ";
    printvec1(mons[k],dim);
  }
}



//
//
// allocation and freeing
//
//

char **charar(long nl, long nh)
     /* This function allocates a set of character pointers */
{
  long nrow = nh-nl+1;
  char **m;
  m=(char **) malloc((size_t)((nrow)*sizeof(char*)));
  if (!m) fprintf(stderr, "allocation failure in charar()\n");
  m-=nl;

  return m;
}


void ifree(int **imat, int n){
  int j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}

void ifree(int **imat, long n){
  long j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
  //fprintf(stderr, "exiting ifree_l\n");fflush(stdout);fflush(stderr);
}

void ifree(double **imat, long n){
  long j;
  //fprintf(stderr, "entering ifree_l\n");fflush(stdout);fflush(stderr);
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
  //fprintf(stderr, "exiting ifree_l\n");fflush(stdout);fflush(stderr);
}

void ifree(double **imat, int n){
  int j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}


void lfree_i(long **imat, int n){
  int j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}

void ulfree_ul(unsigned long **imat, unsigned long n){
  unsigned long j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}


int freearraydat(char **array, int lim){
  int i=lim;
  if(!array){return 0;}
  while(i-- > 0)
    free(array[i]);
  free((FREE_ARG) (array));
  //free((char*)array);
  return 0;
}

//overloading number of lines is ul
int freearraydat(char **array, unsigned long lim){
  unsigned long i=lim;
  if(!array){return 0;}
  while(i-- > 0)
    free(array[i]);
  free((FREE_ARG) (array));
  return 0;
}



double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


bool **bmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a bool matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  bool **m;

  /* allocate pointers to rows */
  m=(bool **) malloc((size_t)((nrow+NR_END)*sizeof(bool*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(bool *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(bool)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

//simple version
void free_imat(int **A, int numA){
  int i;
  for(i=0;i<numA;i++)
    free ((char *)(A[i]));
  if(numA && A)
    free((char *) A);
}

void free_bmat(bool **A, int numA){
  int i;
  for(i=0;i<numA;i++)
    free ((char *)(A[i]));
  if(numA && A)
    free((char *) A);
}


void free_bmatrix(bool **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

ex **exmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a symb matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  ex **m;

  /* allocate pointers to rows */
  m=(ex **) malloc((size_t)((nrow+NR_END)*sizeof(ex*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(ex *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(ex)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_exmatrix(ex **m, long nrl, long nrh, long ncl, long nch)
/* free an ex matrix allocated by exmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}



//
//
// sorting and searching
//
//

// v is an ordered list of integer vectors each of length vlen
// there are n vectors in v. The algorithm perfoms a 
// binary search on v to find vector x
long binisearch(int **v, long n, int *x, int vlen){
  long low, high, mid;
  int cmp;
  low=0;
  high=n-1;
  while(low<=high){
    mid=(low+high)/2;
    if((cmp=isgrtrf(x,v[mid],vlen))>0){ // v[mid] > x
      /* printvec1(x,vlen); */
      /* printvec1(v[mid],vlen); */
      /* cout << "greater" << endl; */
      high=mid-1;
    }
    else if(cmp<0) // x > v[mid]
      low=mid+1;
    else // equal
      return mid;
  }
  return -1;
}


void iswap(int *v, int i, int j){
  int temp = v[i];
  v[i]=v[j];
  v[j]=temp;
}

//overloading
void iswap(double *v, int i, int j){
  double temp = v[i];
  v[i]=v[j];
  v[j]=temp;
}


void iswapi(int **v, int i, int j){
  int *temp = v[i];
  v[i]=v[j];
  v[j]=temp;
}

void ipswap(int **p, int i, int j){
  int *temp=p[i];
  p[i]=p[j];
  p[j]=temp;
}

void lswap(long *p, int i, int j){
  long temp=p[i];
  p[i]=p[j];
  p[j]=temp;
}

void lpswap(long **p, int i, int j){
  long *temp=p[i];
  p[i]=p[j];
  p[j]=temp;
}

//sort a list of integers (ascending order)
void qsortt(int *v, int left, int right){
  int i, last;
  if (left >= right)
    return;
  iswap(v,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++)
    if (v[i]<v[left])
      iswap(v,++last,i);
  iswap(v,left,last);
  qsortt(v, left, last-1);
  qsortt(v, last+1, right);
  return;
}

//sort two lists of integers based on first (ascending order)
void qsort2(int *v, int *v1, int left, int right){
  int i, last;
  if (left >= right)
    return;
  iswap(v,left,(left + right)/2);
  iswap(v1,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++)
    if (v[i]<v[left]){
      iswap(v,++last,i);
      iswap(v1,last,i);
    }
  iswap(v,left,last);
  iswap(v1,left,last);
  qsort2(v, v1, left, last-1);
  qsort2(v, v1, last+1, right);
  return;
}


// sort a list of integer vectors (monomials)

void qsorti(int **v, int dim, long left, long right)
{
  long i, last;

  if (left >= right)
    return;
  iswapi(v,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++)
    if (isgrtrf(v[i],v[left],dim)>0)
      iswapi(v,++last,i);
  iswapi(v,left,last);
  qsorti(v, dim, left, last-1);
  qsorti(v, dim, last+1, right);
}

// sort a list of integer vectors (monomials) and a list of
// coefficients

void qsorti2(int **v, int *cfs, int dim, long left, long right){
  long i, last;

  if (left >= right)
    return;
  /* cout << "swapping v[" << left << "] and v[" << (left + right)/2 << endl; */
  /* printvec1(v[left], dim); */
  /* printvec1(v[(left + right)/2], dim); */
  iswapi(v,left,(left + right)/2);
  /* cout << "now: " << endl; */
  /* printvec1(v[left], dim); */
  /* printvec1(v[(left + right)/2], dim); */
 

  iswap(cfs,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++){
    if (isgrtrf(v[i],v[left],dim)>0){
      ++last;
      iswapi(v,last,i);
      iswap(cfs,last,i);
    }
  }
  iswapi(v,left,last);
  iswap(cfs,left,last);
  qsorti2(v, cfs, dim, left, last-1);
  qsorti2(v, cfs, dim, last+1, right);
}

//overload

void qsorti2(int **v, double *cfs, int dim, long left, long right){
  long i, last;

  if (left >= right)
    return;
  iswapi(v,left,(left + right)/2);
  iswap(cfs,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++){
    if (isgrtrf(v[i],v[left],dim)>0){
      ++last;
      iswapi(v,last,i);
      iswap(cfs,last,i);
    }
  }
  iswapi(v,left,last);
  iswap(cfs,left,last);
  qsorti2(v, cfs, dim, left, last-1);
  qsorti2(v, cfs, dim, last+1, right);
}


int unionlvec(long *a, int na, long *b, int nb, long **outvec){
  int i;
  (*outvec) = (long*) malloc(sizeof(long) * (na+nb));
  for(i=0;i<na;i++)
    (*outvec)[i]=a[i];
  for(i=0;i<nb;i++)
    (*outvec)[na+i]=b[i];
  return na+nb;
}

//two integer lists are equal upto permutation
bool unordlistareeq(int *a, int na, int *b, int nb){
  int i;
  if(na!=nb)
    return 0;
  for(i=0;i<na;i++){
    if(!isinlist(a[i],b,nb))
       return 0;
  }
  return 1;
}

//overload
bool unordlistareeq(long *a, int na, long *b, int nb){
  int i;
  if(na!=nb)
    return 0;
  for(i=0;i<na;i++){
    if(!isinlist(a[i],b,nb))
       return 0;
  }
  return 1;
}



//
//
// Matrices
//
//

/* copy an integer matrix */

int **cpmat(int **mat, int n, int m){
  int i,j;
  int **tmp;
  tmp=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){tmp[i][j]=mat[i][j];}
  }
  return tmp;
}

//overloading the nonallocating version
void cpmat(int **mat, int **out, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){out[i][j]=mat[i][j];}
  }
  return;
}

// is column k of matrix mat (with n rows) equal or plus/minus an earlier column?
bool eqorminus(int **mat, int n, int k){
  int i=0,j;
  int flgminus=1;
  int flg=1;
  if(k==0)
    return 0;

  for(j=0;j<k;j++){ //each previous col
    i=0;flgminus=1;flg=1;
    while(i<n && mat[i][j]==0 && mat[i][k]==0)//initial zeros
      i++;
    if(mat[i][j]==mat[i][k])
      i++;
    else if (mat[i][j]==-mat[i][k]){
      flgminus=-1;i++;
    }
    else
      flg=0;
    while(flg && i<n){
      if(mat[i][j]!=flgminus*mat[i][k])
	flg=0;
      i++;
    }
    if(flg)
      return 1;
  }
  return 0;
}


// list repeated columns (possibly after a sign change) 
// of an n X m matrix mat. Store the information in a binary vector
// of what to keep and what to discard
int reps(int **mat, int n, int m, bool **keeps){
  int k, tot=1;
  (*keeps)[0]=1;
  for(k=1;k<m;k++){// keep the first column
    if(eqorminus(mat, n, k))
      (*keeps)[k]=0;
    else{
      (*keeps)[k]=1;tot++;
    }
  }
  return tot;
}


// does column k of matrix mat (with n rows) appear earlier?
//return the first instance or -1 if first 
int repcol(int **mat, int n, int k){
  int i=0,j;
  int flg=1;
  if(k==0)
    return -1;

  for(j=0;j<k;j++){ //each previous col
    i=0;flg=1;
    while(i<n){
      if(mat[i][j]!=mat[i][k]){
	flg=0;
	break;
      }
      i++;
    }
    if(flg)
      return j;
  }
  return -1;
}

//return total unique columns
int repcolsall(int **mat, int n, int m, int *pat){
  int k, tot=1;
  pat[0]=-1;
  for(k=1;k<m;k++){// keep the first column
    if((pat[k]=repcol(mat, n, k))==-1)
      tot++;
  }
  return tot;
}


// remove columns which, upto sign change, appear previously 
int **redmat(int **mat, int n, int m, int *m1){
  bool *keeps=(bool *)malloc((size_t) ((m)*sizeof(bool)));
  int **tmp;
  int i,j,tot=0;

  (*m1)=reps(mat,n,m,&keeps);
  tmp=imatrix(0, n-1, 0, (*m1)-1);
  for(i=0;i<m;i++){//column
    if(keeps[i]){
      for(j=0;j<n;j++)//row
	tmp[j][tot]=mat[j][i];
      tot++;
    }
  }
  free((char*)keeps);
  return tmp;
}


// remove repeated columns from a matrix: keep only first instance
int **redmat1(int **mat, int n, int m, int *m1, int *pat){
  int **tmp;
  int i,j,tot=0;

  (*m1)=repcolsall(mat,n,m,pat);

  tmp=imatrix(0, n-1, 0, (*m1)-1);
  for(i=0;i<m;i++){//column
    if(pat[i]==-1){//keep
      for(j=0;j<n;j++)//row
	tmp[j][tot]=mat[j][i];
      tot++;
    }
  }
  return tmp;
}



//Are all diagonal entries of a square matrix positive?
bool hasposdiag(int **imat, int n){
  int i;
  for(i=0;i<n;i++){
    if(imat[i][i]<=0)
      return 0;
  }
  return 1;
}

//row total
int rowsum(int **imat, int n, int m, int rowk){
  int i,tot=0;
  if(rowk>=n){
    fprintf(stderr, "ERROR in rowsum: j out of range. Exiting.\n");
    exit(0);
  }
  for(i=0;i<m;i++)
    tot+=imat[rowk][i];
  return tot;
}

//column total
int colsum(int **imat1, int nlen, int mlen, int j){
  int i, tot=0;
  if(j>=mlen){
    fprintf(stderr, "ERROR in colsum: j out of range. Exiting.\n");
    exit(0);
  }
  for(i=0;i<nlen;i++)
    tot+=imat1[i][j];
  return tot;
}



//row maximum: for nonnegative matrices really!
int rowmax(int **imat, int n, int m, int rowk){
  int i,tot=-1000000;
  for(i=0;i<m;i++)
    tot=max(tot,imat[rowk][i]);
  return tot;
}

//col maximum: for nonnegative matrices really!
int colmax(int **imat, int n, int m, int colk){
  int i,tot=-1000000;
  for(i=0;i<n;i++)
    tot=max(tot,imat[i][colk]);
  return tot;
}


// compute imat * diag(vec) and return in imat
// WARNING: modifies imat. No error checking
// integer matrix and vector
int **scalecols(int **imat, int n, int m, int *vec){
  int i,j;
  int **imatout=cpmat(imat,n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      imatout[i][j]*=vec[j];
  }
  return imatout;
}

//overloading (symbolic matrix and symbolic column vector)
//slightly confusing that here we use a column vector
//but in the scalar case we use (effectively) a row vector
matrix scalecols(matrix imat, int n, int m, matrix vec){
  int i,j;
  matrix J(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      J(i,j)=imat(i,j)*vec(j,0);
  }
  return J;
}

//overloading (integer matrix and symbolic column vector)
matrix scalecols(int **imat, int n, int m, matrix vec){
  int i,j;
  matrix J(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      J(i,j)=imat[i][j]*vec(j,0);
  }
  return J;
}


void scalerows(int **imat, int n, int m, int *vec){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      imat[i][j]*=vec[i];
  }
  return;
}

//overloading (same issue as above)
matrix scalerows(matrix imat, int n, int m, matrix vec){
  int i,j;
  matrix J(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      J(i,j)=imat(i,j)*vec(i,0);
  }
  return J;
}

matrix scalerows(int **imat, int n, int m, matrix vec){
  int i,j;
  matrix J(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      J(i,j)=imat[i][j]*vec(i,0);
  }
  return J;
}

int **transposemat(int **imat, int n, int m){
  int i,j;
  int **tmp=imatrix(0, m-1, 0, n-1);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      tmp[j][i]=imat[i][j];
  return tmp;
}

matrix transposemat(matrix imat, int n, int m){
  int i,j;
  matrix imatT(m,n);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      imatT(j,i)=imat(i,j);
  return imatT;
}

//overloading: version with memory already correctly allocated
void transposemat(int **imatin, int n, int m, int **imatout){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      imatout[j][i]=imatin[i][j];
  return;
}

//subtract a row from all rows
void translatematbyrow(int **imatin, int n, int m, int rownum, int **imatout){
  int i,j;
  if(rownum>=n){
    fprintf(stderr, "ERROR in translatematbyrow: row %d out of range (0-%d). Exiting.\n", rownum, n-1);
    exit(0);
  }
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      imatout[i][j]=imatin[i][j]-imatin[rownum][j];
  return;
}


//subtract a column from all columns
void translatematbycol(int **imatin, int n, int m, int colnum, int **imatout){
  int i,j;
  if(colnum>=m){
    fprintf(stderr, "ERROR in translatematbyrow: row %d out of range (0-%d). Exiting.\n", colnum, n-1);
    exit(0);
  }
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      imatout[i][j]=imatin[i][j]-imatin[i][colnum];
  return;
}

/* replace each column C with the column pair C|-C */
int **doublemat(int **imat, int n, int m){
  int **tmp;
  int i,j;
  tmp=imatrix(0, n-1, 0, 2*m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      tmp[i][2*j]=imat[i][j];
      tmp[i][2*j+1]=-imat[i][j];
    }
  }
  return tmp;
}



// convert an integer matrix to a symbolic matrix
matrix imattoexmat(int **A, int n, int m){
  matrix C(n,m);
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      C(i,j)=A[i][j];
  }
  return C;
}

matrix isubmattoexmat(int **A, int *vec1, int n1, int *vec2, int m1){
  matrix C(n1,m1);
  int i,j;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++)
      C(i,j)=A[vec1[i]][vec2[j]];
  }
  return C;
}

//overloading: boolean vector inputs
matrix isubmattoexmat(int **A, int n, int m, bool *vecx, bool *vecy){
  int n1=nonzentries(vecx, n);
  int m1=nonzentries(vecy, m);
  matrix C(n1,m1);
  int i,j, i1=0, j1=0;
  for(i=0;i<n;i++){
    if(vecx[i]){
      j1=0;
      for(j=0;j<m;j++){
	if(vecy[j]){
	  C(i1,j1)=A[i][j];
	  j1++;
	}
      }
      i1++;
    }
  }
  return C;
}

//Get a submatrix from boolean vectors
int **submatfrombool(int **A, int n, int m, bool *vecx, bool *vecy, int *n1, int *m1){
  int i,j, i1=0, j1=0;
  (*n1)=nonzentries(vecx, n);
  (*m1)=nonzentries(vecy, m);
  int **C=imatrix(0, (*n1)-1, 0, (*m1)-1);
  for(i=0;i<n;i++){
    if(vecx[i]){
      j1=0;
      for(j=0;j<m;j++){
	if(vecy[j]){
	  C[i1][j1]=A[i][j];
	  j1++;
	}
      }
      i1++;
    }
  }
  return C;
}

//We've identified a set of points (columns), stored in a
//bool vector "vecy"; now extract the submatrix
int **submatfromcolbool(int **A, int n, int m, bool *vecy, int *m1){
  int i,j, j1=0;
  (*m1)=nonzentries(vecy, m);
  int **C=imatrix(0, n-1, 0, (*m1)-1);

  for(j=0;j<m;j++){//each col
    if(vecy[j]){
      for(i=0;i<n;i++)//each row
	C[i][j1]=A[i][j];
      j1++;
    }
  }
  return C;
}


//We've identified a set of points (rows), stored in a
//bool vector "vecx"; now extract the submatrix
int **submatfromrowbool(int **A, int n, int m, bool *vecx, int *n1){
  int i,j, i1=0;
  (*n1)=nonzentries(vecx, n);
  int **C=imatrix(0, (*n1)-1, 0, m-1);

  for(i=0;i<n;i++){//each row
    if(vecx[i]){
      for(j=0;j<m;j++)//each col
	C[i1][j]=A[i][j];
      i1++;
    }
  }
  return C;
}


// row degree of a matrix
int rowdegree(int **imat, int n, int m){
  int i,j,tmp,deg=0;
  for(i=0;i<n;i++){//each row
    tmp=0;
    for(j=0;j<m;j++){
      if(imat[i][j])
	tmp++;
    }
    deg=max(deg, tmp);
  }
  return deg;
}

// column degree of a matrix
int coldegree(int **imat, int n, int m){
  int i,j,tmp,deg=0;
  for(j=0;j<m;j++){//each column
    tmp=0;
    for(i=0;i<n;i++){
      if(imat[i][j])
	tmp++;
    }
    deg=max(deg, tmp);
  }
  return deg;
}

// generates a submatrix of co-dimension 1 from matrix imat
// removing row i1 and column i2		      
int **submat1(int **imat, int n, int i1, int i2){
  int **tmp;
  int i,j,k,l;
  tmp=imatrix(0, n-2, 0, n-2);i=0;k=0;
  while(k<n){
    j=0;l=0;
    if(k!=i1){
      while(l<n){
	if(l!=i2){
	  tmp[i][j]=imat[k][l];
	  j++;
	}
	l++;
      }
      i++;
    }
    k++;
  }
  return tmp;
}

//submatrix with rows indexed by xc and cols by yc
int **submat(int **imat, int n, int m, int *xc, int n1, int *yc, int m1){
  int **tmp;
  int i,j;
  tmp=imatrix(0, n1-1, 0, m1-1);
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      tmp[i][j]=imat[xc[i]][yc[j]];
    }
  }
  return tmp;
}

//submatrix with rows indexed by xc and cols by yc
matrix submat(matrix imat, int *xc, int n1, int *yc, int m1){
  matrix tmp(n1,m1);
  int i,j;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      tmp(i,j)=imat(xc[i],yc[j]);
    }
  }
  return tmp;
}

// Extract the negative part of a matrix
int **submatminus(int **mat, int n, int m){
  int i,j;
  int **tmp;
  tmp=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if (mat[i][j]>0)
	tmp[i][j]=0;
      else
	tmp[i][j]=mat[i][j];
    }
  }
  return tmp;
}

// Extract a square submatrix
int **submatgen(int **imat, int n, int m, int *i1, int *i2, int dim){
  int **tmp;
  int i,j;
  tmp=imatrix(0, dim-1, 0, dim-1);
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      tmp[i][j]=imat[i1[i]][i2[j]];
    }
  }
  return tmp;
}

//count the nonzero entries
int nonzentries(int *ivec, int n){
  int i,tot=0;
  for(i=0;i<n;i++){
    if(ivec[i])
      tot++;
  }
  return tot;
}

//overloading: boolean vector
int nonzentries(bool *bvec, int n){
  int i,tot=0;
  for(i=0;i<n;i++){
    if(bvec[i])
      tot++;
  }
  return tot;
}

//overloading, matrix
int nonzentries(int **imat, int n, int m){
  int i,j,tot=0;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat[i][j])
	tot++;
    }
  }
  return tot;
}

//first nonzero entry (int vector)
int firstnonz(int *vec, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec[i])
      return i;
  }
  return -1;
}

//first nonzero entry (bool vector)
int firstnonz(bool *vec, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec[i])
      return i;
  }
  return -1;
}

//is nonzero entry (int vector)
bool nonz(int *vec, int n, int i){
  if(vec[i])
    return 1;
  return 0;
}

//is nonzero entry (int vector)
bool nonz(bool *vec, int n, int i){
  if(vec[i])
    return 1;
  return 0;
}



//determinant
int det(int **imat, int n){
  int i;
  int tot=0;
  int **tmpmat;
  if(n==1)
    return imat[0][0];
  for (i=0;i<n;i++){
    tmpmat=submat1(imat, n, 0, i);
    tot+=par(i)*imat[0][i]*det(tmpmat, n-1);
    free_imatrix(tmpmat, 0, n-1, 0, n-1);
  }
  return tot;
}

ex det(matrix imat, int n){
  int xc[n];
  ex detex;
  firstcomb(xc, n, n);
  detex=symbdetsubmat(imat,n,n,xc,xc,n);
  return detex;
}


int detsubmat(int **imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  int tot=0;
  int i2a[dim-1];
  if(dim==1)
    return imat[i1[0]][i2[0]];

  //  i2a=(int*) malloc(sizeof(int) * (dim-1));
  for (i=0;i<dim;i++){
    if(imat[i1[0]][i2[i]]!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot+=par(i)*imat[i1[0]][i2[i]]*detsubmat(imat, n, m, i1+1, i2a, dim-1);
    }

  }
  //  free((FREE_ARG) (i2a));
  return tot;

}

//Here x and y are boolean vectors of lengths n and m respectively
//assumed to have the same number of nonzero entries
int detsubmat(int **imat, int n, int m, bool *x, bool *y){
  int i,dimx,dimy,tot;
  if(!x){dimx=n;}else{dimx=nonzentries(x,n);}
  if(!y){dimy=m;}else{dimy=nonzentries(y,m);}
  if(dimx!=dimy){fprintf(stderr, "ERROR in detsubmat (%d X %d): need a nonempty square matrix. EXITING.\n",dimx,dimy);exit(0);}
  int i1[dimx];
  int i2[dimx];
  tot=0;
  if(!x){for(i=0;i<n;i++){i1[tot++]=i;}}
  else{for(i=0;i<n;i++){if(nonz(x,n,i)){i1[tot++]=i;}}}
  tot=0;
  if(!y){for(i=0;i<m;i++){i2[tot++]=i;}}
  else{for(i=0;i<m;i++){if(nonz(y,m,i)){i2[tot++]=i;}}}
  
  return detsubmat(imat, n, m, i1, i2, dimx);

}



// multiply two matrices A and B of dimensions nXr and rXm
// to get the integer matrix C of dimension nXm
void multAB(int **A, int **B, int **C, int n, int r, int m){
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      C[i][j]=0;
      for(k=0;k<r;k++)
	C[i][j]+=A[i][k]*B[k][j];
    }
  }
  return;
}

matrix multAB(matrix A, int **B, int n, int r, int m){
  int i,j,k;
  matrix C(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      C(i,j)=0;
      for(k=0;k<r;k++)
	C(i,j)+=A(i,k)*B[k][j];
    }
  }
  return C;
}

matrix multAB(int **A, matrix B, int n, int r, int m){
  int i,j,k;
  matrix C(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      C(i,j)=0;
      for(k=0;k<r;k++)
	C(i,j)+=A[i][k]*B(k,j);
    }
  }
  return C;
}

matrix multAB(matrix A, matrix B, int n, int r, int m){
  int i,j,k;
  matrix C(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      C(i,j)=0;
      for(k=0;k<r;k++)
	C(i,j)+=A(i,k)*B(k,j);
    }
  }
  return C;
}


// SV where S is nXm and V is mXn
matrix multAB(matrix S, matrix V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S(i,k)*V(k,j);
    }
  }
  return J;
}


//overloading
matrix multAB(int **S, matrix V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S[i][k]*V(k,j);
    }
  }
  return J;
}

//overloading
int **multAB(int **S, int **V, int n, int m){
  int **J=imatrix(0, n-1, 0, m-1);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J[i][j]=0;
      for(k=0;k<m;k++)
	J[i][j]+=S[i][k]*V[k][j];
    }
  }
  return J;
}



//output AB^t where A,B are nXm matrices
matrix multABT(matrix S, matrix V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S(i,k)*V(j,k);
    }
  }
  return J;
}

//overloading
matrix multABT(matrix S, int **V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S(i,k)*V[j][k];
    }
  }
  return J;
}

//overloading
matrix multABT(int **S, matrix V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S[i][k]*V(j,k);
    }
  }
  return J;
}

//overloaded
int **multABT(int **S, int **V, int n, int m){
  int **J=imatrix(0, n-1, 0, n-1);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J[i][j]=0;
      for(k=0;k<m;k++)
	J[i][j]+=S[i][k]*V[j][k];
    }
  }
  return J;
}


// the symbolic determinant of a symbolic matrix
// returns in expanded form (very inefficient!
// Leibniz formula: only use for small matrices)

ex symbdetsubmat(matrix imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  ex tot=0;
  int i2a[dim-1];
  if(dim==1)
    return imat(i1[0], i2[0]);
  for (i=0;i<dim;i++){
    if(imat(i1[0], i2[i])!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1);

    }
  }

  return expand(tot);

}




//
// Get the rank of an integer matrix
// Using GiNAC rank function
//

int matrank(int **A, int n, int m){
  matrix J = imattoexmat(A, n, m);
  return J.rank();
}

int submatrank(int **A, int *vec1, int n, int *vec2, int m){
  matrix J = isubmattoexmat(A, vec1, n, vec2, m);
  return J.rank();
}

//overloading: boolean vectors
int submatrank(int **A, int n, int m, bool *vecx, bool *vecy){
  matrix J = isubmattoexmat(A, n, m, vecx, vecy);
  return J.rank();
}



// Construct the second additive compound matrix
matrix AdComp2(matrix J, int ktot){
  int i,j, k, m, cnk;
  int **r;
  cnk=ktot*(ktot-1)/2;

  matrix J2(cnk,cnk);
  r=imatrix(0, cnk-1, 1, 2);

  j=0;
  for(k=0;k<ktot-1;k++){
    for(m=k+1;m<ktot;m++){
      r[j][1]=k;
      r[j][2]=m;
      j++;
    }
  }


  for(i=0;i<cnk;i++){
    for(k=0;k<cnk;k++){
      if((r[i][1] == r[k][1]) && (r[i][2] == r[k][2]))
	J2(i,k) = J(r[i][1],r[i][1])+J(r[i][2],r[i][2]);
      else if (r[i][1] == r[k][1])
	J2(i,k) = J(r[i][2],r[k][2]);
      else if (r[i][1] == r[k][2])
	J2(i,k) = -J(r[i][2],r[k][1]);
      else if (r[i][2] == r[k][1])
	J2(i,k) = -J(r[i][1],r[k][2]);
      else if (r[i][2] == r[k][2])
	J2(i,k) = J(r[i][1],r[k][1]);
      else
	J2(i,k) = 0;

    }
  }

  free_imatrix(r, 0, cnk-1, 1, 2);
  return J2;
}

// return the determinant of the second additive compound in expanded form
ex AdComp2Det(matrix J, int n, int debug){
  ex detex;
  int m=(n*(n-1))/2;
  int xc[m];
  matrix J2 = AdComp2(J,n);

  if(debug){
    fprintf(stderr, "\n###Entering AdComp2Det: Here are the matrix and its second additive compound.\n");
    printmat(J,n,n);
    printmat(J2,m,m);
  }

  firstcomb(xc,m,m);
  detex = symbdetsubmat(J2, m, m, xc, xc, m);
  if(debug){cerr << "det(J^[2]) = " << detex << endl;}
  return expand(detex);
}


// get the sum of principal k X k minors
// of the symbolic n X n matrix M
ex getminorsum0(matrix M, int n, int k){
  int xc[k];
  int flag1;
  ex out=0;
  if(!k){
    out=1;
    return out;
  }
  firstcomb(xc,n,k);

  flag1=1;
  while(flag1==1){
    out+=expand(symbdetsubmat(M,n,n,xc,xc,k));
    flag1=nextcomb(xc,n,k);
  }
  return out;      
}



//
//
// Graph theory
//
//

void strongconnect(int i, int **mat, int n, int *inds, int *lowlink, int *ind, int ***S, int *totS, int *Sc, int *Scpos){
  int j,k;

  /* fprintf(stderr, "current: "); */
  /* printvec1(Sc, *Scpos); */
  // Set the depth index for v to the smallest unused index
  inds[i]=*ind;lowlink[i]=*ind;(*ind)++;
  Sc[(*Scpos)++]=i;//add to stack

  // Consider successors of v
  for (j=0;j<n;j++){
    if(mat[i][j]){// each out-edge
      if(inds[j]<0){
	strongconnect(j, mat, n, inds, lowlink, ind, S, totS, Sc, Scpos);
	lowlink[i]=min(lowlink[i],lowlink[j]);
      }
      else if(isinlist(j, Sc, *Scpos))
	lowlink[i]=min(lowlink[i],inds[j]);
    }
  }
  // If node i is a root node, pop the stack and generate an SCC
  if (lowlink[i]==inds[i]){
    //Create SCC
    if((*totS)==0)
      (*S) = (int**) malloc(sizeof(int*) * 1);
    else
      (*S) =(int**) realloc((*S), sizeof(int*) *((*totS)+1));

    k=(*Scpos)-1;
    while(Sc[k]!=i)
      k--;

    (*S)[(*totS)] = (int*) malloc(sizeof(int) * ((*Scpos)-k+1));

    (*S)[(*totS)][0]=(*Scpos)-k;// first element is size
    for(j=k;j<(*Scpos);j++)
      (*S)[(*totS)][j-k+1]=Sc[j];
    (*totS)++;(*Scpos)=k;//reset stack
  }

}

// My implementation of Tarjan's algorithm
// Based on pseudocode on Wikipedia
// mat is the n X n adjacency matrix of the digraph
// output is the list of SCCs and totS, the number of SCCs
int **Tarjan(int **mat, int n, int *totS){
  //  output: set of strongly connected components (set of vertex-sets)
  int i, ind=0; //ind increments as we encounter vertices
  int *inds=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *lowlink=(int *)malloc((size_t) ((n)*sizeof(int)));
  int **S=NULL;// to store the SCCs
  int *Sc= (int *)malloc((size_t) ((n)*sizeof(int)));//stack
  int Scpos=0;//stack-counter

  for(i=0;i<n;i++)
    inds[i]=-1;

  for(i=0;i<n;i++){
    if(inds[i]<0)
      strongconnect(i,mat,n,inds,lowlink,&ind, &S,totS, Sc, &Scpos);
  }

  free((char *)(Sc));
  free((char *)inds);
  free((char *)lowlink);

  /* fprintf(stderr, "Total number of SCCs in complex graph = %d:\n", *totS); */
  /* for(i=0;i<(*totS);i++){ */
  /*   fprintf(stderr, "i=%d, sz=%d\n", i,S[i][0]); */
  /*   for(ind=0;ind<S[i][0];ind++) */
  /*     fprintf(stderr, "%d ", S[i][1+ind]); */
  /*   fprintf(stderr, "\n"); */
  /* } */
  return S;
}



//like strongconnect, but just count the SCCs, don't store them

void strongconnect_light(int i, int **mat, int n, int *inds, int *lowlink, int *ind, int *totS, int *Sc, int *Scpos){
  int j,k;

  /* fprintf(stderr, "current: "); */
  /* printvec1(Sc, *Scpos); */
  // Set the depth index for v to the smallest unused index
  inds[i]=*ind;lowlink[i]=*ind;(*ind)++;
  Sc[(*Scpos)++]=i;//add to stack

  // Consider successors of v
  for (j=0;j<n;j++){
    if(mat[i][j]){// each out-edge
      if(inds[j]<0){
	strongconnect_light(j, mat, n, inds, lowlink, ind, totS, Sc, Scpos);
	lowlink[i]=min(lowlink[i],lowlink[j]);
      }
      else if(isinlist(j, Sc, *Scpos))
	lowlink[i]=min(lowlink[i],inds[j]);
    }
  }
  // If node i is a root node, pop the stack
  if (lowlink[i]==inds[i]){
    k=(*Scpos)-1;
    while(Sc[k]!=i)
      k--;
    (*totS)++;(*Scpos)=k;//reset stack
  }

}


// Like Tarjan, but don't store the SCCs

int Tarjan_light(int **mat, int n){
  int i, ind=0; //ind increments as we encounter vertices
  int *inds=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *lowlink=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *Sc= (int *)malloc((size_t) ((n)*sizeof(int)));
  int Scpos=0;//stack-counter
  int totS=0;

  for(i=0;i<n;i++)
    inds[i]=-1;

  for(i=0;i<n;i++){
    if(inds[i]<0)
      strongconnect_light(i,mat,n,inds,lowlink,&ind, &totS, Sc, &Scpos);
  }

  free((char *)Sc);
  free((char *)inds);
  free((char *)lowlink);

  return totS;
}

//sum of positive elements in a vector
int sumpos(int *vec, int m){
  int i,tot=0;
  for(i=0;i<m;i++){
    if(vec[i]>0)
      tot+=vec[i];
  }
  return tot;
}
