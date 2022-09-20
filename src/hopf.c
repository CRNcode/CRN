/* Copyright (C) 2010-2022, Murad Banaji
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
#include "hopf.h"
#include "pos.h"
#include "symbolic.h"
#include "convex.h"
#include "semidef.h"
#include "analysereacs.h"
#include <limits.h>

//
// change all coefficients to have magnitude 1 at the cost of repetition
//

void elongate(int ***explst, int **cflst, long *r, int numv){
  long s,s1=0;
  int tmp[(*r)][numv], cftmp[(*r)];
  int i,j,k,t;

  for(s=0;s<(*r);s++){
    for(i=0;i<numv;i++)
      tmp[s][i]=(*explst)[s][i];
    cftmp[s]=(*cflst)[s];
  }
  ifree((*explst),(*r));
  free((char*)(*cflst));

  for(s=0;s<(*r);s++)
    s1+=abs(cftmp[s]);

  (*explst)= (int**) malloc(sizeof(int*)*s1);
  (*cflst)= (int*) malloc(sizeof(int)*s1);
  for(s=0;s<s1;s++)
    (*explst)[s]= (int*) malloc(sizeof(int)*numv);

  s1=0;
  for(s=0;s<(*r);s++){
    t=abs(cftmp[s]);
    for(j=0;j<t;j++){
      for(k=0;k<numv;k++){
	(*explst)[s1][k]=tmp[s][k];
	(*cflst)[s1]=cftmp[s]/t;
      }
      s1++;
    }
  }
  (*r)=s1;
  return;
}

// make homogeneous using last variable, based on maximum degree of two polynomials
// use with care!
void homify(int **lst1, long r1, int **lst2, long r2, int numv){
  long i;
  int j,tot1[r1],tot2[r2],totmax=0;
  // get highest degree
  for(i=0;i<r1;i++){
    tot1[i]=0;
    for(j=0;j<numv;j++)
      (tot1[i])+=lst1[i][j];
    if(tot1[i]>totmax)
      totmax=tot1[i];
  }
  for(i=0;i<r2;i++){
    tot2[i]=0;
    for(j=0;j<numv;j++)
      (tot2[i])+=lst2[i][j];
    if(tot2[i]>totmax)
      totmax=tot2[i];
  }
  // make homogeneous
  for(i=0;i<r1;i++)
    lst1[i][numv] = totmax-tot1[i];
  for(i=0;i<r2;i++)
    lst2[i][numv] = totmax-tot2[i];
  return;

}



// Get the Real part of det(J+wI)
// Variables of J assumed to be in canonical form v%d
// otherwise algorithm will fail with unpredictable consequences.
ex ReJpiI(matrix J, int n, int numv, bool hom){
  ex v[numv+2];
  char str1[10];
  int k;
  ex tmp0,tmp=0,tmp1=0;
  bool fl;

  for(k=0;k<=numv+1;k++){
    sprintf(str1, "v%d", k);
    v[k]=get_possymbol(str1);
  }

  fl=1;
  for(k=n;k>=0;k-=2){// starts with det(J) and alternates sign
    tmp0=getminorsum0(J,n,k);
    if(fl){tmp+=tmp0;fl=0;}else{tmp-=tmp0;fl=1;}
  }
  if(hom)
    tmp1=polyhom(tmp,v[numv+1],0);
  else
    tmp1=tmp;
  return expand(tmp1);
}
//if(!hom){tmp=getminorsum0(J,n,k);}else{tmp=expand(pow(v[numv+1],n-k)*getminorsum0(J,n,k));}

// Get the imaginary part of det(J+wI)
// Variables of J assumed to be in canonical form v%d
// otherwise algorithm will fail with unpredictable consequences.
// hom=2 means homogenise minimally
ex ImJpiI(matrix J, int n, int numv, int hom){
  ex v[numv+2];
  char str1[10];
  int k;
  ex tmp0,tmp=0,tmp1=0;
  bool fl;

  for(k=0;k<=numv+1;k++){
    sprintf(str1, "v%d", k);
    v[k]=get_possymbol(str1);
  }

  fl=1;
  for(k=n-1;k>=0;k-=2){// starts with J_{n-1} and alternates sign
    tmp0=getminorsum0(J,n,k);
    if(fl){tmp+=tmp0;fl=0;}else{tmp-=tmp0;fl=1;}
  }
  if(hom==1)
    tmp1=polyhom(tmp,v[numv+1],1);
  else if(hom==2)
    tmp1=polyhom(tmp,v[numv+1],0);
  else
    tmp1=tmp;

  return expand(tmp1);
}
//if(!hom){tmp=getminorsum0(J,n,k);}else if(hom==1){tmp=expand(pow(v[numv+1],n-k)*getminorsum0(J,n,k));}else if(hom==2){tmp=expand(pow(v[numv+1],n-k-1)*getminorsum0(J,n,k));}

//special n=3 at bifurcation form
//i.e. we evaluate along det(J^[2])=0
ex ReJpiI3(matrix J, int numv, bool hom){
  ex v[numv+2];
  char str1[10];
  int k;
  ex J1=0,J2=0,tmpout;

  for(k=0;k<=numv+1;k++){
    sprintf(str1, "v%d", k);
    v[k]=get_possymbol(str1);
  }
  J1=getminorsum0(J,3,1);
  J2=getminorsum0(J,3,2);

  if(!hom)
    tmpout=expand(J1*J2-J1);
  else
    tmpout=expand(J1*J2-pow(v[numv+1],2)*J1);

  return tmpout;
}

// 3X3 matrix J, admits pair of purely imaginary eigenvalues; assume that
// external routine has confirmed that J is indeed a 3 X 3 matrix.
// det(J^[2]) can change sign and J2 > 0 or J1*J3>0 (so that when
// J3 = J1*J2, J2 is positive). Note that these are sufficient, but 
// by no means necessary conditions for an imaginary pair: it is quite 
// possible for J2 and/or J1*J3 to be unsigned, but there exist solutions 
// to J2>0, J1*J2=J3. It is also not necessary for J1*J2-J3 to take
// all signs for it to be zero.
// J assumed to be homogeneous and in standard form
int admitsIpair3(matrix J, int numv, int maxpppdegree, int *pppdegused, int debug){
  int allflg,ret=0,retJ2,ret2,ret13=-10;
  ex J1=0,J2=0,J3=0,dtJ2;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering admitsIpair3.\n");}

  J1=getminorsum0(J,3,1);
  J2=getminorsum0(J,3,2);
  J3=getminorsum0(J,3,3);
  dtJ2=expand(J1*J2-J3);
  if((retJ2=ispospoly(expand(dtJ2), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==-3 && ((ret2=ispospoly(expand(J2), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2 || (ret13=ispospoly(expand(J3*J1), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2)){
/* if(ret2==-3) */
/* exit(0); */
    ret=1;
  }

  if(debug){fprintf(stderr, "Signs: J3-J1*J2: %d, J2: %d, J1*J3: %d.\n", retJ2, ret2, (ret13==-10)?ispospoly(expand(J1*J3), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull):ret13);}
  if(debug){fprintf(stderr, "Exiting admitsIpair3: returning %d.\n", ret);}
  return ret;
}

//The 4D version
int admitsIpair4(matrix J, int numv, int maxpppdegree, int *pppdegused, int debug){
  int allflg,ret=0,retJ2,retalt=-10,ret13;
  ex J1=0,J2=0,J3=0,J4=0,dtJ2,Jalt;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering admitsIpair4.\n");}

  J1=getminorsum0(J,4,1);
  J2=getminorsum0(J,4,2);
  J3=getminorsum0(J,4,3);
  J4=getminorsum0(J,4,4);
  dtJ2=expand(-J1*J1*J4-J3*J3+J1*J2*J3);
  Jalt=expand(J2*J2*(J1*J1*J4+J3*J3));
  if((retJ2=ispospoly(expand(dtJ2), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==-3 && ((ret13=ispospoly(expand(J3*J1), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2 || (retalt=ispospoly(expand(Jalt), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2)){
/* if(ret2==-3) */
/* exit(0); */
    ret=1;
  }

  if(debug){fprintf(stderr, "Signs: -J1*J1*J4-J3*J3+J1*J2*J3: %d, J1*J3: %d, J2*J2*(J1*J1*J4+J3*J3): %d.\n", retJ2, ret13, (retalt==-10)?ispospoly(expand(Jalt), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull):retalt);}
  if(debug){fprintf(stderr, "Exiting admitsIpair4: returning %d.\n", ret);}
  return ret;
}



//The 2D version
int admitsIpair2(matrix J, int numv, int maxpppdegree, int *pppdegused, int debug){
  int allflg,ret=0,ret1, ret2;
  ex J1=0,J2=0;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering admitsIpair2.\n");}

  J1=getminorsum0(J,2,1);//trace
  J2=getminorsum0(J,2,2);//determinant
  if((ret1=ispospoly(expand(J1), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==-3 && (ret2=ispospoly(expand(J2), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2){
/* if(ret2==-3) */
/* exit(0); */
    ret=1;
  }

  if(debug){
    fprintf(stderr, "Signs: J1: %d, J2: %d.\n", ret1, ret2);
    fprintf(stderr, "Exiting admitsIpair2: returning %d.\n", ret);
  }
  return ret;
}



// Get the Real and imaginary part of det(J+wI)
// Variables of J assumed to be in canonical form v%d
// otherwise algorithm will fail with unpredictable consequences.
void ReImJpiI(ex tmp1, ex tmp2, int numv, int ***explst1, int **cflst1, long *r1, int ***explst2, int **cflst2, long *r2, bool hom){
  ex v[numv+2];
  char str1[10];
  int i;
  ex tmp;
  int **lst;
  int nv;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  if(!hom){nv=numv;}else{nv=numv+1;}

  //real part
  polytointlist0(tmp1, &lst, cflst1, r1);
  (*explst1)=moncflisttoexplist0(lst,*cflst1,*r1,nv,1);// redefines cflst
  ifree(lst, *r1);
  elongate(explst1,cflst1,r1,nv);

  //imaginary part
  polytointlist0(tmp2, &lst, cflst2, r2);
  (*explst2)=moncflisttoexplist0(lst,*cflst2,*r2,nv,1);// redefines cflst
  ifree(lst, *r2);
  elongate(explst2,cflst2,r2,nv);

  //cout << "sanity check = " <<  expand(tmp1-expltopolysimp(*explst1,*cflst1,*r1,numv)) << endl;
  //cout << "sanity check = " <<  expand(tmp2-expltopolysimp(*explst2,*cflst2,*r2,numv)) << endl;

  // Added to lower powers and to deal with MA
  if(hom)
    homify(*explst1,*r1,*explst2,*r2,numv);

  return;
}

//Overloading: version with the matrix as input
void ReImJpiI(matrix J, int n, int numv, int ***explst1, int **cflst1, long *r1, int ***explst2, int **cflst2, long *r2, bool hom){
  ex tmp1=ReJpiI(J, n, numv, hom);
  ex tmp2=ImJpiI(J, n, numv, hom);

  ReImJpiI(tmp1, tmp2, numv, explst1, cflst1, r1, explst2, cflst2, r2, hom);
  return;
}

//@@ This can be modified to take the approach in degen (in pos.c), rather
// than directly examining the minor-sums 
// Check if the matrix has the potential for B-T bifurcations
// The J_rk and J_{rk-1} be simultaneously zero? (rk is the rank of 
// the stoichiometric matrix).
int BTpotential(matrix J, int n, int rk, int numv, int maxpppdegree, int *pppdegused, int debug){
  int allflg,ret=0,retf,retp;
  ex Jf=0,Jp=0, JJ;
  int debugfull=(debug<=0)?0:debug-1;
  int numvinit, numvfinal, deg;

  if(debug){fprintf(stderr, "\n###Entering BTpotential.\n");}

  Jf=getminorsum0(J,n,rk);
  Jp=getminorsum0(J,n,rk-1);

  if((retf=ispospoly(expand(Jf), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2 || retf==-2 || (retp=ispospoly(expand(Jp), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2 || retp==-2){
    if(debug){cerr << "One polynomial is positive so SOS definitely positive. Exiting BTpotential.\n";}
    return 0;
  }

  JJ=polyhomsimp(expand(Jf*Jf+Jp*Jp), &numvinit, &numvfinal, &deg, 1, 0, debugfull);
  ret=ispospolynonneg(JJ, numvfinal, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull);
  if(ret==2){
    if(debug){fprintf(stderr, "A double degeneracy is ruled out, Exiting BTpotential.\n");}
    return 0;
  }

  if(debug){fprintf(stderr, "A double degeneracy could not be ruled out. Exiting BTpotential.\n");}
  return 1;
}



// squares of two lists as a polynomial

ex lsts_to_poly(int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int numv){
  ex v[numv+1];
  char str1[10];
  int i;
  ex tmp=0, tmp1=0;
  long r;

  for(i=0;i<numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  for(r=0;r<r1;r++)
    tmp1+=expltomono(cflst1[r],explst1[r],numv);
  tmp+=expand(pow(tmp1,2));
  tmp1=0;

  for(r=0;r<r2;r++)
    tmp1+=expltomono(cflst2[r],explst2[r],numv);
  tmp+=expand(pow(tmp1,2));

  return tmp;
}

// The input is the polynomial tmp, the output is the polynomial tmp1 which is tmp 
// with a number (clq1+clq2) of perfect squares subtracted. These perfect squares 
// ("cliques") come respectively entirely from Re and Im and are given as input in 
// the form of the vectors used1 and used2 (set in BigSq1). If used1[i]=m, then this 
// tells us that the ith vector from Re figures in the mth clique extracted from Re. 
// Similarly, for used2. 

// Current dangers: it seems that we use the actual coefficients in constructing the 
// squares to remove. If these are larger than 1, we might be unnecessarily "using up" 
// our positive monomials which may serve to dominate something else?

ex remsqs(ex tmp, ex *sqs, int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int numv, int *used1, int *used2, int clq1, int clq2, bool q){
  ex tmp2=tmp;
  ex tmp1,extot=0;
  bool swtch=0;
  long s, m;
  if(!q){cout << endl << "current state of the squares:\n" << endl;}

  // The following "incompatible square" blocks should never execute:
  // in BigSq1 used1 and used2 entries are >=0. They are for legacy. Remove?
  for(s=0;s<r1;s++){
    if(used1[s]==-1){
      if(!swtch){
	extot+=expltomono(1,explst1[s],numv);swtch=1;
      }
      else
	extot-=expltomono(1,explst1[s],numv);
      cerr << "ERROR in remsqs: This should not execute.\n";exit(0);
    }
  }
  tmp1=expand(pow(extot,2));tmp2-=tmp1;(*sqs)+=pow(extot,2);
  if(extot!=0 && !q){cout << "(incompatible square):  " << pow(extot,2) << endl;}

  extot=0;swtch=0;
  for(s=0;s<r2;s++){
    if(used2[s]==-2){
      if(!swtch){
	extot+=expltomono(1,explst2[s],numv);swtch=1;
      }
      else
	extot-=expltomono(1,explst2[s],numv);
      cerr << "ERROR in remsqs: this should not execute.\n";exit(0);
    }
  }
  tmp1=expand(pow(extot,2));tmp2-=tmp1;(*sqs)+=pow(extot,2);
  if(extot!=0 && !q){cout << "(incompatible square):  " << pow(extot,2) << endl;}

  
  for(m=1;m<=clq1;m++){// each clique from Re
    extot=0;
    for(s=0;s<r1;s++){// or should we use coefficients 1?
      if(used1[s]==m)
	extot+=expltomono(cflst1[s],explst1[s],numv);
    }
    tmp1=expand(pow(extot,2));tmp2-=tmp1;(*sqs)+=pow(extot,2);
    if(!q){cout << pow(extot,2) << endl;}
  }

  for(m=1;m<=clq2;m++){// each clique from Im
    extot=0;
    for(s=0;s<r2;s++){// or should we use coefficients 1?
      if(used2[s]==m)
	extot+=expltomono(cflst2[s],explst2[s],numv);
    }
    tmp1=expand(pow(extot,2));tmp2-=tmp1;(*sqs)+=pow(extot,2);
    if(!q){cout << pow(extot,2) << endl;}
  }
  if(!q){cout << endl;}
  return tmp2;
}

// sequential version of J2pI_recur. Rather than being recursive, 
// it assumes a wrapper function which calls it on loop.
// negative returns are failures, 0 means continue, -2 means a fundamental
// failure. -1 means a failure: too many interations, 
// or no cliques found, but we can continue: hope to recover.
// 1 means success. 

int J2pI_seq(ex tmp, ex *tmpout, int n, int numv, int *totiter, float *bestq, ex *best, int *lastsplitfail, int maxlevel, int maxcliques, int maxiter, bool q){
 
  int dim;
  int **lst,**explst,*cflst;
  long r,s;
  float ns;
  int *needcount,*needexcl;
  int **needs, **neededby;
  int count=0,negcount=0;
  int maxsz=5000;
  int tight=0;//loose
  long **left=NULL,**right=NULL;
  int *numleft=NULL,*numright=NULL,*absnegmax=NULL;
  int totcliques=0,j, retval;
  float totalq=0;
  ex tmp3;

  if(*totiter>=maxiter)
    return -1;

  (*totiter)++;
  dim=polytointlist(tmp, &lst, &cflst, &r);//dim=order of monomials
  for(s=0;s<r;s++){
    count++;
    if(cflst[s]<0)
      negcount++;
  }

  // Nothing to remove?
  if(negcount==0){
    if(tmp!=0){ // We have genuinely removed the negative monomials
      fprintf(stderr, "neg count: 0, success!\n\n");
      if(!q){cout << "final polynomial (J2pI_seq) = " << tmp << endl;
	cout << "neg count: 0, success!\n\n";}
      ifree(lst, r);free((char *)cflst);
      return 1;
    }
    else{ // We have probably recovered the original squares
      if(!q){fprintf(stderr, "neg count: 0, tmp = 0, failure!\n\n");
	cout << "neg count: 0, tmp = 0, failure!\n\n";}
      ifree(lst, r);free((char *)cflst);
      return -2; 
    }
  }


  explst=moncflisttoexplist(lst,dim,cflst,r,numv+1,1);// redefines cflst
  ifree(lst, r);

  // Gather data on how remaining negative terms can be dominated: 
  // "who needs who?", "who is needed by who?"
  // This information will be used to construct cliques

  if(*totiter==1 && !q)// first time only
    listprint(cflst, explst, numv+1, r);
  needcount = (int *)malloc((size_t) (r*sizeof(int))); 
  needexcl = (int *)malloc((size_t) (r*sizeof(int)));
  inittozero(needcount,r);inittozero(needexcl,r);
    
  needs=(int **)malloc((size_t)(r*sizeof(int*)));
  neededby=(int **)malloc((size_t)(r*sizeof(int*)));
  for(s=0;s<r;s++){
    neededby[s]=(int *)malloc((size_t)(sizeof(int)));neededby[s][0]=0;
    needs[s] = (int*)malloc(sizeof(int));needs[s][0]=0; //number needed
  }

  for(s=0;s<r;s++){
    if(cflst[s]<0){
      if(*totiter==1)
	ns=numsplits(s,numv+1,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);// changed the final flag to make less verbose
      else
	ns=numsplits(s,numv+1,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);
      if(ns<0.0001){// not able to split any more
	if(!q){cout << "failure to split:\n   " << cflst[s] << ", "; printvec1(explst[s],numv+1);}
	veccp(explst[s],numv+1,lastsplitfail); // store the monomial which cannot be dominated
	ifree(needs, r);ifree(neededby, r);ifree(explst, r);
	free((char *)cflst); free((char *)needcount);free((char *)needexcl);
	return -1; 
      }
      totalq+=ns;
      //	cout << "numsplits[" <<  s << "]: " << ns << endl;
    }
  }
  totalq/=((float)negcount);

  // Low quality means there won't be much choice on how to 
  // dominate negative monomials
  if(!q){fprintf(stderr, "nc: %d/%d, quality: %f\n", negcount, count, totalq);}
  if(!q){cout << "quality:  " << totalq << ",   negcount: " << negcount << endl;}
  if(*totiter>1 && totalq>(*bestq)){
    (*bestq)=totalq;(*best)=tmp;
    if(!q){cout << "best q: " << *bestq << endl;}
  }

  totcliques=allcliques1(explst,cflst,r,numv+1,neededby,needs,needexcl,&left,&right,&numleft,&numright,&absnegmax,maxsz,maxcliques,tight,0,1,0);// penultimate flag keeps cliques small (see "excflag" in comments on "allcliques1"). quiet
  for(j=0;j<totcliques;j++){// list all good cliques
    tmp3=makesquare(explst, cflst, left[j],numleft[j],right[j],numright[j],numv+1,absnegmax[j],1);
    if(!q){cout << "clique found (J2pI_seq): " << tmp3 << endl;}
  }

  if(totcliques==0){ // Allow a looser approach to clique building
    totcliques=allcliques1(explst,cflst,r,numv+1,neededby,needs,needexcl,&left,&right,&numleft,&numright,&absnegmax,maxsz,maxcliques,tight,0,0,0);//quet
    for(j=0;j<totcliques;j++){// list all good cliques
      tmp3=makesquare(explst, cflst, left[j],numleft[j],right[j],numright[j],numv+1,absnegmax[j],1);
      if(!q){cout << "not tight clique found (J2pI_seq): " << tmp3 << endl;}
    }
  }
  ifree(neededby, r);ifree(needs, r);
  free((char *)needcount);free((char *)needexcl);

  // If cliques are found (whether tight or not), store the first one, 
  // assumed to be the "best" clique by some measure.
  // (It is actually removed in the wrapper function J2pI_seqWrap)
  if(totcliques){
    j=0;
    tmp3=makesquare(explst, cflst, left[j],numleft[j],right[j],numright[j],numv+1,absnegmax[j],1);
    if(!q){cout << endl << "removing: " << tmp3 << endl << endl;}
    *tmpout=tmp3;retval=0;
  }
  else
    retval=-1;

  if(left){lfree_i(left,totcliques);}if(right){lfree_i(right,totcliques);}
  ifree(explst, r);free((char *)cflst);
  if(numleft){free((char*)numleft);}if(numright){free((char*)numright);}if(absnegmax){free((char*)absnegmax);}

  //cout << tmp << endl;

  return retval; 
}



int J2pI_seqWrap(ex tmp, ex *sqs, ex *finalpos, int n, int numv, int *totiter, float *bestq, ex *best, int *lastsplitfail, int maxlevel, int maxcliques, int maxiter, bool q){
  ex tmpout, tmp1=tmp;
  int retval=0;
  while(!retval){
    tmpout=0;
    retval=J2pI_seq(tmp1, &tmpout, n, numv, totiter, bestq, best, lastsplitfail, maxlevel, maxcliques, maxiter, q);
    (*sqs)+=tmpout;
    tmp1=tmp1-expand(tmpout);
    if(retval==1)
      (*finalpos)=tmp1;
  }
  return retval;
}

// The idea is a product of polynomials (a sum of the powers)

int *sum(int *a, int *b, int n){
  int i;
  int *out=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++)
    out[i]=a[i]+b[i];
  return out;

}

 
// count the honest splits (i.e. which respect Re and Im) of a negative monomial

void numhonspl(int ** explst, int *cflst, long r, int numv, int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int *out){
  long s,s1,s2;
  int *vc;

  inittozero(out, r);

  for(s1=0;s1<r1;s1++){
    for(s2=s1+1;s2<r1;s2++){
      if(cflst1[s1]*cflst1[s2]<0){
	vc=sum(explst1[s1],explst1[s2],numv);
	if((s=binisearch(explst,r,vc,numv))>=0 && cflst[s]<0){
	  out[s]++;
	}
      }
    }
  }
  for(s1=0;s1<r2;s1++){
    for(s2=s1+1;s2<r2;s2++){
      if(cflst2[s1]*cflst2[s2]<0){
	vc=sum(explst2[s1],explst2[s2],numv);
	if((s=binisearch(explst,r,vc,numv))>=0 && cflst[s]<0){
	  out[s]++;
	}
      }
    }
  }

  return;
}

// honest splits for (negative) monomial allt[indx] from (Re^2 + Im^2) are computed by 
// looking at the lists in Re and Im separately (and unsquared).
// Within Re (say) each pair of monomials of opposite sign are multiplied, and 
// we check if the result is of the form of allt[indx]
// We set the needs1 and neededby1 flags accordingly. 

int honestsplits(long indx, int n, int **allt, int *cfs, long r, int **lst1, int *cfs1, int r1, int **lst2, int *cfs2, int r2, int *used1, int *used2, int **needcount1, int **needs1, int **neededby1, int **needcount2, int **needs2, int **neededby2, int **needexcl, bool q){
  long s,s1;
  int *vc;
  int j,j1,j2;
  int spl=0;
  int powsplit=0;
  int cf=cfs[indx];
  //fprintf(stderr, "entering honestsplits.\n");
  if(!q){
    cout << "---------------------" << endl;
    /* cout << cf << ", "; printvec1(allt[indx],n); */expltomonoprint(cfs[indx],allt[indx],n);
  }
  for(s=0;s<r1;s++){
    for(s1=s+1;s1<r1;s1++){
      if(cfs1[s]*cfs1[s1]<0 && !used1[s] && !used1[s1]){
	vc=sum(lst1[s],lst1[s1],n);
	if(areequal(vc, allt[indx], n)){
	  (needs1[indx][0])+=2;j=(needs1[indx][0])+1;
	  needs1[indx]=(int*) realloc((needs1[indx]), sizeof(int)*(j)); 
      
	  (neededby1[s][0])+=2;(neededby1[s1][0])+=2;
	  j1=(neededby1[s][0])+1;j2=(neededby1[s1][0])+1;
	  neededby1[s]=(int*) realloc((neededby1[s]), sizeof(int)*(j1));
	  neededby1[s1]=(int*) realloc((neededby1[s1]), sizeof(int)*(j2));
	  powsplit+=2*abs(cfs1[s]*cfs1[s1]);;
	  needs1[indx][j-2]=s;needs1[indx][j-1]=s1;
	  neededby1[s][j1-2]=indx;neededby1[s][j1-1]=s1;
	  neededby1[s1][j2-2]=indx;neededby1[s1][j2-1]=s;
	  ((*needcount1)[s])++;((*needcount1)[s1])++;
	  if(!q){
	    /* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(1,lst1[s],n);
	    /* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(1,lst1[s1],n);
	  }
	  spl++;

	}
	if(vc)
	  free((char*)vc);
      }
    }
  }

  for(s=0;s<r2;s++){
    for(s1=s+1;s1<r2;s1++){
      if(cfs2[s]*cfs2[s1]<0 && !used2[s] && !used2[s1]){
	vc=sum(lst2[s],lst2[s1],n);
	if(areequal(vc, allt[indx], n)){
	  (needs2[indx][0])+=2;j=(needs2[indx][0])+1;
	  needs2[indx]=(int*) realloc((needs2[indx]), sizeof(int)*(j)); 
      
	  (neededby2[s][0])+=2;(neededby2[s1][0])+=2;
	  j1=(neededby2[s][0])+1;j2=(neededby2[s1][0])+1;
	  neededby2[s]=(int*) realloc((neededby2[s]), sizeof(int)*(j1));
	  neededby2[s1]=(int*) realloc((neededby2[s1]), sizeof(int)*(j2));
	  powsplit+=2*abs(cfs2[s]*cfs2[s1]);
	  needs2[indx][j-2]=s;needs2[indx][j-1]=s1;
	  neededby2[s][j1-2]=indx;neededby2[s][j1-1]=s1;
	  neededby2[s1][j2-2]=indx;neededby2[s1][j2-1]=s;
	  ((*needcount2)[s])++;((*needcount2)[s1])++;
	  if(!q){
	    /* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(1,lst2[s],n);
	    /* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(1,lst2[s1],n);
	  }
	  spl++;
	}
	if(vc)
	  free((char*)vc);
      }
    }
  }
  if(!q){cout << "powsplit = " << powsplit << endl;}
  if(powsplit==abs(cf)){ //no choice
    if(!q){cout << "All splits needed (honestsplits).\n";}
    (*needexcl)[indx]=1;
  }

  if(powsplit<abs(cf)){
    if(!q){
      cout << "Warning: powers don't add up.\n";
      fprintf(stderr, "failure to split (honestsplits).\n");
    }
    return 0;
    //exit(0);
  }

  //fprintf(stderr, "exiting honestsplits. spl=%d\n", spl);
  return spl;
}

// modcliques takes a pair of (positive) monomials with indices (ul and ur) 
// assumed to be needed in a square and assumed both to come from Re or both 
// from Im. "k" tells us whether they are both from Re or Im. Assume that 
// both are not already used in the same clique (otherwise we have immediate 
// failure). If neither is used, then create a new clique from them. If one 
// is used in an existing clique, then add the other. If both are already 
// used, but in different cliques, then merge the two cliques. 

// Question: shouldn't ul and ur be "long" rather than "int"?

void modcliques(int k, int *used1, long r1, int *used2, long r2, int ul, int ur, int *clq1, int *clq2, bool q){
  int clqtmp;
  long s;
  if(k==1){// Comes from Re
    if(used1[ul] && used1[ur] && used1[ul]==used1[ur]){// both already used in the same clique
      fprintf(stderr, "Something has gone wrong (modcliques): both monomials appear already to belong to the same clique from Re.\n");
      exit(0);
    }
    if(used1[ul]==0 && used1[ur]==0){// both new -- make new clique
      (*clq1)++;
      if(!q){cout << "creating new cLique(1): " << (*clq1) << endl;}
      used1[ul]=(*clq1);used1[ur]=(*clq1);
    }
    else if(used1[ul]==0){//only left is new -- add to clique
      if(!q){cout << "adding to cLique(1): " << used1[ur] << endl;}
      used1[ul]=used1[ur];
    }
    else if(used1[ur]==0){//only right is new
      if(!q){cout << "adding to cLique(1): " << used1[ul] << endl;}
      used1[ur]=used1[ul];
    }
    else{//both in use, but in different cliques -- merge cliques
      if(!q){cout << "merging cLiques(1): " << used1[ul] << ", " << used1[ur] << endl;}
      clqtmp=used1[ur];
      for(s=0;s<r1;s++){
	if(used1[s]==clqtmp)
	  used1[s]=used1[ul];
      }
    }
  }
  else if(k==2){ // comes from Im
    if(used2[ul] && used2[ur] && used2[ul]==used2[ur]){// both already used in the same clique
      fprintf(stderr, "Something has gone wrong (modcliques): both monomials appear already to belong to the same clique from Im.\n");
      exit(0);
    }
    if(used2[ul]==0 && used2[ur]==0){// both new -- make new clique
      (*clq2)++;
      if(!q){cout << "creating new cLique(2): " << (*clq2) << endl;}
      used2[ul]=(*clq2);used2[ur]=(*clq2);
    }
    else if(used2[ul]==0){//only left is new -- add to right clique
      if(!q){cout << "adding to cLique(2): " << used2[ur] << endl;}
      used2[ul]=used2[ur];
    }
    else if(used2[ur]==0){//only right is new -- add to left clique
      if(!q){cout << "adding to cLique(2): " << used2[ul] << endl;}
      used2[ur]=used2[ul];
    }
    else{//both in use, but in different cliques -- merge cliques
      if(!q){cout << "merging cLiques(2): " << used2[ul] << ", " << used2[ur] << endl;}
      clqtmp=used2[ur];
      for(s=0;s<r2;s++){
	if(used2[s]==clqtmp)
	  used2[s]=used2[ul];
      }
    }
  }
  return;
}

// count the negative terms in an expression

int negcnt(ex tmp){
  int **lst,*cflst;
  long r, s;
  int negcount=0;
  polytointlist0(tmp, &lst, &cflst, &r);
  //extolst(tmp, numv, &explst, &cflst, &r);
  for(s=0;s<r;s++){
    if(cflst[s]<0)
      negcount++;
  }
  ifree(lst,r); free((char*)cflst); 
  return negcount;
}

// Is the product of the first two monomials, the third?

bool isprod(int *a, int *b, int *prod, int numv){
  int i;
  for(i=0;i<numv;i++){
    if((a[i]+b[i])!=prod[i])
      return 0;
  }
  return 1;
}


// [Assuming called with the switch "best"]. Findsplitpair1 starts with 
// the (negative) monomial explst[t1]. Checks if it is the product of two 
// monomials from either Re or Im, appearing with opposite signs, and 
// either at least one is unused, or both are used but appear in different 
// cliques. In each case does a dummy run of modcliques to see how this 
// affect the negative count. Takes the best result and stored the endpoints 
// as ul and ur for later use in modcliques. 

// [Without the switch "best"] Doesn't do any dummy runs of modcliques or 
// check if the negative count diminishes - just takes the first valid 
// splitting and stores the endpoints ul and ur. 

int findpairsplit1(ex tmp, int **explst, long t1, long *u1, long *u2, int **explst1, int **explst2, int *cflst1, int *cflst2, long r1, long r2, int *used1, int *used2, int clq1, int clq2, int numv, bool best, bool q){
  long s, s1, s2, u1tmp, u2tmp;
  int *used1tmp;
  int *used2tmp;
  int clq1tmp=clq1,clq2tmp=clq2;
  ex tmp2;
  // ncold keeps track of the old count of negative monomials
  int nc,ncold,retval=0;
  ex sqs=0;

  ncold=100*negcnt(tmp); // hopefully not all worse than this!
  /* //building cliques */
  /* for(s1=0;s1<r1;s1++){ */
  /*   for(s2=s1+1;s2<r1;s2++){ */
  /*     if(isprod(explst1[s1],explst1[s2],explst[t1],numv) && cflst1[s1]*cflst1[s2]<0 && used1[s1]*used1[s2]==0 && used1[s1]>=0 && used1[s2]>=0){ */
  /* 	if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;} */
  /* 	return 1; */
  /*     } */
  /*   } */
  /* } */
  /* for(s1=0;s1<r2;s1++){ */
  /*   for(s2=s1+1;s2<r2;s2++){ */
  /*     if(isprod(explst2[s1],explst2[s2],explst[t1],numv) && cflst2[s1]*cflst2[s2]<0 && used2[s1]*used2[s2]==0 && used2[s1]>=0 && used2[s2]>=0){ */
  /* 	if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;} */
  /* 	return 2; */
  /*     } */
  /*   } */
  /* } */

  //merging

  if(!best){ // Just take the first valid splitting
    for(s1=0;s1<r1;s1++){
      for(s2=s1+1;s2<r1;s2++){
	if(isprod(explst1[s1],explst1[s2],explst[t1],numv) && cflst1[s1]*cflst1[s2]<0 && (used1[s1]*used1[s2]==0 || used1[s1]!=used1[s2]) && used1[s1]>=0 && used1[s2]>=0){
	  if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	  return 1;
	}
      }
    }
    for(s1=0;s1<r2;s1++){
      for(s2=s1+1;s2<r2;s2++){
	if(isprod(explst2[s1],explst2[s2],explst[t1],numv) && cflst2[s1]*cflst2[s2]<0 && (used2[s1]*used2[s2]==0 || used2[s1]!=used2[s2]) && used2[s1]>=0 && used2[s2]>=0){
	  if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	  return 2;
	}
      }
    }
  }
  else{ // Look for the "best" splitting (maximum lowering of the neg count)

    // Does lots of dummy runs of modcliques and remsqs to see which clique 
    // modification (i.e., new clique, built clique or merged cliques) 
    // leads to the best improvement of the neg count. Stores the 
    // indices of the best modification. These will be used in the 
    // actual run of modcliques


    used1tmp=(int *)malloc((size_t) (r1*sizeof(int)));
    used2tmp=(int *)malloc((size_t) (r2*sizeof(int)));

    // Question: seems a rather wasteful way of checking how explst[t1] 
    // can be constructed?

    for(s1=0;s1<r1;s1++){
      for(s2=s1+1;s2<r1;s2++){
	if(isprod(explst1[s1],explst1[s2],explst[t1],numv) && cflst1[s1]*cflst1[s2]<0 && (used1[s1]*used1[s2]==0 || used1[s1]!=used1[s2]) && used1[s1]>=0 && used1[s2]>=0){
	  for(s=0;s<r1;s++){used1tmp[s]=used1[s];}for(s=0;s<r2;s++){used2tmp[s]=used2[s];}
	  clq1tmp=clq1;clq2tmp=clq2;
	  if(cflst1[s1]>0){u1tmp=s1;u2tmp=s2;}else{u1tmp=s2;u2tmp=s1;}
	  modcliques(1, used1tmp, r1, used2tmp, r2, u1tmp, u2tmp, &clq1tmp, &clq2tmp, q);
	  tmp2=remsqs(tmp, &sqs, explst1, cflst1, r1, explst2, cflst2, r2, numv, used1tmp, used2tmp, clq1tmp, clq2tmp, q);
	  nc=negcnt(tmp2); if(!q){cout << "nc = " << nc << endl;}
	  if(nc<ncold){//lower negative count
	    if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    retval=1;ncold=nc;
	  }
	}
      }
    }
    for(s1=0;s1<r2;s1++){
      for(s2=s1+1;s2<r2;s2++){
	if(isprod(explst2[s1],explst2[s2],explst[t1],numv) && cflst2[s1]*cflst2[s2]<0 && (used2[s1]*used2[s2]==0 || used2[s1]!=used2[s2]) && used2[s1]>=0 && used2[s2]>=0){
	  for(s=0;s<r1;s++){used1tmp[s]=used1[s];}for(s=0;s<r2;s++){used2tmp[s]=used2[s];}
	  clq1tmp=clq1;clq2tmp=clq2;
	  if(cflst2[s1]>0){u1tmp=s1;u2tmp=s2;}else{u1tmp=s2;u2tmp=s1;}
	  modcliques(2, used1tmp, r1, used2tmp, r2, u1tmp, u2tmp, &clq1tmp, &clq2tmp, q);
	  tmp2=remsqs(tmp, &sqs, explst1, cflst1, r1, explst2, cflst2, r2, numv, used1tmp, used2tmp, clq1tmp, clq2tmp, q);
	  nc=negcnt(tmp2);if(!q){cout << "nc = " << nc << endl;}
	  if(nc<ncold){
	    if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    retval=2;ncold=nc;
	  }
	}
      }
    }
    free((char*)used1tmp);free((char*)used2tmp);
  }

  return retval;
}

bool ueq2v(int *u, int *v, int n){
  int i;
  for(i=0;i<n;i++){
    if(2*v[i]!=u[i])
      return 0;
  }
  return 1;
}

bool flgset(int newf, int a, int b){
  if(newf==2){
    if(!a && !b)
      return 1;
  }
  else if(newf==1){
    if((!a || !b) && a>=0 && b>=0)
      return 1;
  }
  else{
    if(a>=0 && b>=0 && (a!=b || (a==0 && b==0)))
      return 1;
  }
  return 0;
}


// new=2 is a new pair; new=1 has at least one new member; new=0 allows both members to be in use

int validpair1(int **explst, long t1, long t2, long *u1, long *u2, int **explst1, int **explst2, int *cflst1, int *cflst2, long r1, long r2, int *used1, int *used2, int numv, int newf){
  long s1, s2;

  //try list 1
  for(s1=0;s1<r1;s1++){
    if(ueq2v(explst[t1],explst1[s1], numv)){
      for(s2=s1+1;s2<r1;s2++){
	//cout << used1[s1]<< ", " << used1[s2] << endl;
	//cout << "flgset = " << flgset(newf, used1[s1], used1[s2]) << endl;
	if(ueq2v(explst[t2],explst1[s2], numv) &&  flgset(newf, used1[s1], used1[s2])){
	  if(cflst1[s1]*cflst1[s2]<0){
	    if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 1;
	  }
	}
      }
    }
    else if(ueq2v(explst[t2],explst1[s1], numv)){
     for(s2=s1+1;s2<r1;s2++){
       //cout << used1[s1]<< ", " << used1[s2] << endl;
       //cout << "flgset = " << flgset(newf, used1[s1], used1[s2]) << endl;
	if(ueq2v(explst[t1],explst1[s2], numv) && flgset(newf, used1[s1], used1[s2])){
	  if(cflst1[s1]*cflst1[s2]<0){
	    if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 1;
	  }
	}
      }
    }
  }
  // try list 2
  for(s1=0;s1<r2;s1++){
    if(ueq2v(explst[t1],explst2[s1], numv)){
      for(s2=s1+1;s2<r2;s2++){
	//cout << used2[s1]<< ", " << used2[s2] << endl;
	//cout << "flgset = " << flgset(newf, used2[s1], used2[s2]) << endl;
	if(ueq2v(explst[t2],explst2[s2], numv) && flgset(newf, used2[s1], used2[s2])){
	  if(cflst2[s1]*cflst2[s2]<0){
	    if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 2;
	  }
	}
      }
    }
    else if(ueq2v(explst[t2],explst2[s1], numv)){
      for(s2=s1+1;s2<r2;s2++){
	//cout << used2[s1]<< ", " << used2[s2] << endl;
	//cout << "flgset = " << flgset(newf, used2[s1], used2[s2]) << endl;
	if(ueq2v(explst[t1],explst2[s2], numv) && flgset(newf, used2[s1], used2[s2])){
	  if(cflst2[s1]*cflst2[s2]<0){
	    if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 2;
	  }
	}
      }
    }
  }
  return 0;
}


int randseedpair(int **explst, int **explst1, int **explst2, int *cflst, int *cflst1, int *cflst2, long r, long r1, long r2, int *used1, int *used2, int numv, int **needs1, int **needs2, long *ul, long *ur, int flg){
  long s,s1,s2;
  int i,j;
  long u1t,u2t;
  long *p1,*p2;
  int *kk;
  int nump=0;
  int timeint;
  time_t timepoint;
  int retval;
  timeint = time(&timepoint); /* set and convert time to an integer */
  srand(timeint); /*use this integer as a random seed */
  //cout << "randnum = " << rand() << endl;
  //fprintf(stderr, "entering randseedpair\n");
  for(s=0;s<r;s++){ // new pair
    if(cflst[s]<0){
      if(needs1[s][0]){
	for(j=1;j<needs1[s][0];j+=2){
	  s1=needs1[s][j];s2=needs1[s][j+1];
	  if(flgset(flg,used1[s1],used1[s2])){
	    if(cflst1[s1]>0){u1t=s1;u2t=s2;}else{u1t=s2;u2t=s1;}
	    addnewto1Dlarray(&p1, nump, u1t);
	    addnewto1Dlarray(&p2, nump, u2t);
	    addnewto1Darray(&kk, nump, 1);
	    nump++;
	  }
	}
      }
      if(needs2[s][0]){
	for(j=1;j<needs2[s][0];j+=2){
	  s1=needs2[s][j];s2=needs2[s][j+1];
	  if(flgset(flg,used2[s1],used2[s2])){
	    if(cflst2[s1]>0){u1t=s1;u2t=s2;}else{u1t=s2;u2t=s1;}
	    addnewto1Dlarray(&p1, nump, u1t);
	    addnewto1Dlarray(&p2, nump, u2t);
	    addnewto1Darray(&kk, nump, 2);
	    nump++;
	  }
	}
      }
    }
  }

  //fprintf(stderr, "nump = %d\n",nump);
  if(!nump)
    retval=0;
  else{
    i=rand()%nump;
    //fprintf(stderr, "choosing %d\n", i);
    *ul=p1[i];*ur=p2[i];retval=kk[i];
  }
  //fprintf(stderr, "exiting randseedpair\n");
  //  free((char*)p1);free((char*)p2);free((char*)kk);
  return retval;
}


int randseedpair_old(int **explst, int **explst1, int **explst2, int *cflst, int *cflst1, int *cflst2, long r, long r1, long r2, int *used1, int *used2, int numv, int **needs, int *needexcl, long *ul, long *ur){
  long s;
  int i,j,k=0;
  long u1t,u2t;
  long *p1,*p2;
  int *kk;
  int nump=0;
  int timeint;
  time_t timepoint;
  int retval;
  timeint = time(&timepoint); /* set and convert time to an integer */
  srand(timeint); /*use this integer as a random seed */
  //cout << "randnum = " << rand() << endl;

  for(s=0;s<r;s++){ // new pair
    if(cflst[s]<0){
      if(needs[s][0]>=2){
	for(j=1;j<needs[s][0];j+=2){
	  if((k=validpair1(explst, needs[s][j], needs[s][j+1], &u1t, &u2t, explst1, explst2, cflst1, cflst2, r1, r2, used1, used2, numv,0))){ // new pair
	    addnewto1Dlarray(&p1, nump, u1t);
	    addnewto1Dlarray(&p2, nump, u2t);
	    addnewto1Darray(&kk, nump, k);
	    nump++;
	  }
	}
      }
    }
  }
  fprintf(stderr, "nump = %d\n",nump);
  if(!nump)
    retval=0;
  else{
    i=rand()%nump;
    *ul=p1[i];*ur=p2[i];retval=kk[i];
  }
  //  free((char*)p1);free((char*)p2);free((char*)kk);
  return retval;
}

void forcedclq(long indx, long r, int *needexcl, int **needs, int **neededby, int *used, int *tot, bool q){
  long m;
  int j,k,flg;
  if(!needexcl[indx]) // not forced
    return;
  if(used[indx]==-1)// already dealt with
    return;
  // flg=1 means odds on left; flg=2 means odds on right
  used[indx]=-1;
  flg=0;
  for(j=1;j<=needs[indx][0];j++){// step through pairs needed by indx
    if(used[needs[indx][j]]){
      if(j%2==used[needs[indx][j]]%2){
	flg=1;break;
      }
      else{
	flg=2;break;
      }
    }
  }
  // set used flags
  if(!flg || flg==1){ //odd,even...
    for(j=1;j<=needs[indx][0];j++){
      if(used[needs[indx][j]] && used[needs[indx][j]]%2!=j%2){
	if(!q){fprintf(stderr, "impossible situation(1).\n");}/* exit(0); */used[indx]=-1;(*tot)=0;return;
      }
      else{
	if(j%2==0 && !used[needs[indx][j]]){
	  used[needs[indx][j]]=2;(*tot)++;
	}
	else if(!used[needs[indx][j]]){
	  used[needs[indx][j]]=1;(*tot)++;
	}
      }
    }
  }
  else if(flg==2){ //even,odd...
    for(j=1;j<=needs[indx][0];j++){
      if(used[needs[indx][j]] && used[needs[indx][j]]!=j%2+1){
	if(!q){fprintf(stderr, "impossible situation(2).\n");}/* exit(0); */used[indx]=-1;(*tot)=0;return;
      }
      else{
	if(j%2==0 && !used[needs[indx][j]]){
	  used[needs[indx][j]]=1;(*tot)++;
	}
	else if(!used[needs[indx][j]]){
	  used[needs[indx][j]]=2;(*tot)++;
	}
      }
    }
  }

  for(j=1;j<=needs[indx][0];j++){//step through positive terms surrounding indx
    for(k=1;k<=neededby[needs[indx][j]][0];k+=2){// step through negative terms dominated by each positive term
      m=neededby[needs[indx][j]][k];
      forcedclq(m,r,needexcl,needs,neededby,used,tot,q);
    }
  }
  return;
}

int isfreehonest(int **explst, int *cflst, int numv, long r, int **explst1, int *cflst1, int r1, int **explst2, int *cflst2, int r2, int *used, int *used1, int *used2, int *used1tmp, int *used2tmp){
  long s,s1;
  int *vc;
  int flg=0, sgn=0;
  //fprintf(stderr, "entering isfreehonest\n");
  for(s=0;s<r;s++){
    if(used[s]>0 && !iseven(explst[s],numv)){
      //fprintf(stderr, "couldn't be halved.\n");
      return 0;
    }
  }
  for(s=0;s<r;s++){
    if(used[s]>0){// left or right in clique
      vc=halve(explst[s], numv);
      if((s1=binisearch(explst1,r1,vc,numv))>=0){
	if(flg==2){
	  //fprintf(stderr, "mixed(1).\n");
	  return 0;
	}
	if(used1[s1]){// not available?
	  //fprintf(stderr, "failing s1 = %ld\n", s1);
	  return 0;
	}
	else{
	  //fprintf(stderr, "s1 = %ld\n", s1);
	  used1tmp[s1]=1;
	}
	flg=1;
	if(!sgn){
	  if(used[s]==1)
	    sgn=cflst1[s1];
	  else
	    sgn=-cflst1[s1];
	}
	if((used[s]==1 && sgn*cflst1[s1]<0) || (used[s]==2 && sgn*cflst1[s1]>0)){
	  //fprintf(stderr, "sign change(1).\n");
	  return 0;
	}
      }
      else if((s1=binisearch(explst2,r2,vc,numv))>=0){
	if(flg==1){
	  //fprintf(stderr, "mixed(2).\n");
	  return 0;
	}
	if(used2[s1]){//available?
	  //fprintf(stderr, "failing s1 = %ld\n", s1);
	  return 0;
	}
	else{
	  //fprintf(stderr, "s1 = %ld\n", s1);
	  used2tmp[s1]=1;
	}
	flg=2;
	if(!sgn){
	  if(used[s]==1)
	    sgn=cflst2[s1];
	  else
	    sgn=-cflst2[s1];
	}
	if((used[s]==1 && sgn*cflst2[s1]<0) || (used[s]==2 && sgn*cflst2[s1]>0)){
	  //fprintf(stderr, "sign change(2).\n");
	  return 0;
	}

      }
      else // not found
	return 0;
    }
  }
  //fprintf(stderr, "returning %d\n", flg);
  return flg;
}


// find the first forced clique, or the largest forced clique
// firstbigall=0: first; firstbigall=1: biggest; firstbigall=-1: biggest; firstbigall=2: all.

ex honestforced(int **explst, int *cflst, int numv, long r, int *needexcl, int **needs, int **neededby,int **explst1, int *cflst1, int r1, int **explst2, int *cflst2, int r2, int *used1, int *used2, int *used1tmp, int *used2tmp, bool hon, int firstbigall, bool q){
  long s,s1;
  bool flg=0;
  int *used=(int *)malloc((size_t) (r*sizeof(int)));
  int *usedtmp=(int *)malloc((size_t) (r*sizeof(int)));
  int *used1tmp1=(int *)malloc((size_t) (r1*sizeof(int)));
  int *used2tmp1=(int *)malloc((size_t) (r2*sizeof(int)));
  ex exout;
  long *leftinds,*rightinds;
  int numleft=0, numright=0;
  int tot=0,mrk=0,totold=0,totold1=10000;
  int clqtyp;
  //fprintf(stderr, "entering honestforced\n");

  inittozero(used,r);
  inittozero(used1tmp,r1);
  inittozero(used1tmp1,r1);
  inittozero(used2tmp,r2);
  inittozero(used2tmp1,r2);

  for(s=0;s<r;s++){
    if(cflst[s]<0 && needexcl[s]){
      for(s1=0;s1<r;s1++){usedtmp[s1]=0;}
      tot=0;
      forcedclq(s, r, needexcl, needs, neededby, usedtmp, &tot,q);
      inittozero(used1tmp,r1);inittozero(used2tmp,r2);
      clqtyp=isfreehonest(explst, cflst, numv, r, explst1, cflst1, r1, explst2, cflst2, r2, usedtmp, used1, used2, used1tmp, used2tmp);
      if(!clqtyp && !q){cout << "found a dishonest forced clique\n";}
      if(!hon || (hon && clqtyp)){
	if(!firstbigall){// first
	  flg=1;for(s1=0;s1<r;s1++){used[s1]=usedtmp[s1];}
	  for(s1=0;s1<r1;s1++){used1tmp1[s1]=used1tmp[s1];}for(s1=0;s1<r2;s1++){used2tmp1[s1]=used2tmp[s1];}
	  break;
	}
	else if(firstbigall==1){//biggest
	  if(tot>totold){
	    if(!q){fprintf(stderr, "new bigger clique found: tot = %d.\n", tot);}
	    flg=1;totold=tot;
	    for(s1=0;s1<r;s1++){used[s1]=usedtmp[s1];}
	    for(s1=0;s1<r1;s1++){used1tmp1[s1]=used1tmp[s1];}for(s1=0;s1<r2;s1++){used2tmp1[s1]=used2tmp[s1];}
	  }
	}
	else if(firstbigall==-1){//smallest
	  if(tot<totold1){
	    fprintf(stderr, "new smaller clique found: tot = %d.\n", tot);
	    flg=1;totold1=tot;
	    for(s1=0;s1<r;s1++){used[s1]=usedtmp[s1];}
	    for(s1=0;s1<r1;s1++){used1tmp1[s1]=used1tmp[s1];}for(s1=0;s1<r2;s1++){used2tmp1[s1]=used2tmp[s1];}
	  }
	}
	else{//all
	  flg=1;mrk++;
	  if(tot>totold){
	    for(s1=0;s1<r;s1++){totold=tot;used[s1]=usedtmp[s1];}
	  }
	  for(s1=0;s1<r1;s1++){if(used1tmp[s1]){used1tmp1[s1]=mrk;}}for(s1=0;s1<r2;s1++){if(used2tmp[s1]){used2tmp1[s1]=mrk;}}
	}
      }
    }
  }
  if(flg){	  
    for(s1=0;s1<r1;s1++){used1tmp[s1]=used1tmp1[s1];}for(s1=0;s1<r2;s1++){used2tmp[s1]=used2tmp1[s1];}
    for(s=0;s<r;s++){
      if(used[s]==1){numleft++;}else if(used[s]==2){numright++;}
    }
    if(!q){fprintf(stderr, "numleft=%d, numright=%d\n", numleft, numright);}
    leftinds=(long *)malloc((size_t) (numleft*sizeof(long)));
    rightinds=(long *)malloc((size_t) (numright*sizeof(long)));

    numleft=0;numright=0;
    for(s=0;s<r;s++){
      if(used[s]==1)
	leftinds[numleft++]=s;
      else if(used[s]==2)
	rightinds[numright++]=s;
    }
    exout=makesquare(explst, cflst, leftinds,numleft,rightinds,numright,numv,10,1);
    free((char*)leftinds);free((char*)rightinds);
  }
  free((char*)used);free((char*)usedtmp);free((char*)used1tmp1);free((char*)used2tmp1);
  //fprintf(stderr, "exiting honestforced\n");
  return exout;
}


// spf =0 seems to often lead to recovery of the entire original squares and failure. 

int BigSq1(ex tmp, ex tmpfull, ex *sqs, ex *finalpos, int n, int numv, int **explst1, int *cflst1, int *used1, long r1, int **explst2, int *cflst2, int *used2, long r2, int *clq1, int *clq2, int *level, bool q){
  int **explst,*cflst;
  long r,s,s1,ul,ur;
  int k,totiter=0;
  //int level1=0;
  ex ex1,best,ff;
  int *needexcl,*needexcl0;
  int *needcount,**needs,**neededby;
  int *needcount1,**needs1,**neededby1;
  int *needcount2,**needs2,**neededby2;
  float ns0,bestq=0;
  int ns,maxlevel=800,maxiter=100;
  bool splitfail=0,simpsq;
  long failed=0,failed1=0;
  bool spf=1;// use the last failure at each stage
  int count, negcount,maxcliques=1000;//maxcliques=0 means no maximum
  int retstat,flg1=0,flg2=0;
  int *used1tmp,*used2tmp;
  int *out;
  int *lastsplitfail=(int *)malloc((size_t) (numv*sizeof(int)));

  // First we try our basic routine (J2pI_seqWrap, or J2pI_recurr). 
  // We may be successful, or we may find that we arrive at a failure 
  // to dominate some negative monomial. If our basic routine fails, 
  // then we move on to more sophisticated things. We use the last 
  // failure (i.e. the troublesome negative monomial) as the "seed". 

  // The recursive version (J2pI_recurr) can get stuck deep in the recursion, 
  // exploring the various branches. The sequential version doesn't 
  // face this problem. A failure instead provides a monomial which 
  // is the basis of another attempt. 

  inittozero(lastsplitfail, numv);
  (*level)++;
  if(!q){fprintf(stderr, "------------------------------\nBigSq level = %d\n", *level);}

  if(*level > 5){maxiter=10000;}
  //if((retstat=J2pI_recurr(tmp, numv-1, &level1, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, 0))==1){// success
  if((retstat=J2pI_seqWrap(tmp, sqs, finalpos, n, numv-1, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, q))==1){// success
    fprintf(stderr, "Seems to be successful.\n");
    return 1;
  }
  else if(retstat==-2){
    if(!q){fprintf(stderr, "Seems to have failed.\n");}
    return -1;
  }
  if(!q){cout << "failed on: "; printvec1(lastsplitfail, numv);}

  // We tried J2pI_seq until we were unable to split monomial lastsplitfail 
  // We will now have to try something else. 

  negcount=0;count=0;
  extolst(tmp, numv, &explst, &cflst, &r);
  for(s=0;s<r;s++){
    count++;
    if(cflst[s]<0)
      negcount++;
  }
  // Question: what is this about? How can this ever occur? Legacy?

  if(!negcount){
    fprintf(stderr, "no more negative terms! success! (BigSq1)\n");
    cout << "no more negative terms! (BigSq1)" << endl;
    cout << "final polynomial (BigSq1): " << tmp << endl;
    exit(0);
  }
  if(!q){cout << "negative terms: " << negcount << "/" << count << endl;}
  //cout << tmp << endl;

  // Just for output, and just at level 1. How many honest splits 
  // (respecting Re and Im) are there of each negative negative 
  // monomial? (Store response in "out"). 

  if(*level==1 && !q){
    cout << "ways of splitting:\n";
    out=(int *)malloc((size_t) (r*sizeof(int)));
    inittozero(out, r);
    numhonspl(explst, cflst, r, numv, explst1, cflst1, r1, explst2, cflst2, r2, out);
    for(s=0;s<r;s++){
      if(cflst[s]<0){
	cout << out[s] << ", " << cflst[s] << ", "; printvec1(explst[s], numv);
      }
    }
    free((char*)out);
    //    exit(0);
  }


  //fprintf(stderr, "got to this point in BigSq1\n");

  // Gather information about splits and honest splits of the full polynomial (Re^2 + Im^2)

  needexcl = (int *)malloc((size_t) (r*sizeof(int))); 
  needexcl0 = (int *)malloc((size_t) (r*sizeof(int))); 
  needs=(int **) malloc((size_t)(r*sizeof(int*)));
  needs1=(int **) malloc((size_t)(r*sizeof(int*)));
  needs2=(int **) malloc((size_t)(r*sizeof(int*)));
  for(s=0;s<r;s++){
    needexcl[s]=0;needexcl0[s]=0;
    needs[s]=(int*) malloc(sizeof(int));needs[s][0]=0;
    needs1[s]=(int*) malloc(sizeof(int));needs1[s][0]=0;
    needs2[s]=(int*) malloc(sizeof(int));needs2[s][0]=0;
  }

  needcount=(int *)malloc((size_t) (r*sizeof(int)));
  needcount1=(int *)malloc((size_t) (r1*sizeof(int)));
  needcount2=(int *)malloc((size_t) (r2*sizeof(int)));     
  neededby=(int **) malloc((size_t)(r*sizeof(int*)));
  neededby1=(int **) malloc((size_t)(r1*sizeof(int*)));
  neededby2=(int **) malloc((size_t)(r2*sizeof(int*)));
  for(s=0;s<r;s++){
    needcount[s]=0;
    neededby[s]=(int *) malloc((size_t)(sizeof(int)));neededby[s][0]=0;
  }
  for(s=0;s<r1;s++){
    needcount1[s]=0;
    neededby1[s]=(int *) malloc((size_t)(sizeof(int)));neededby1[s][0]=0;
  }
  for(s=0;s<r2;s++){
    needcount2[s]=0;
    neededby2[s]=(int *) malloc((size_t)(sizeof(int)));neededby2[s][0]=0;
  }

  for(s=0;s<r;s++){
    if(cflst[s]<0){
      if(*level==2)
	ns0=numsplits(s,numv,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);
      else
	ns0=numsplits(s,numv,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);
      if(ns0<0.0001){// not able to split any more
	if(!q){cout << "fundamental failure to split (Bigsq1):\n   " << cflst[s] << ", "; printvec1(explst[s],numv);}
	failed=s;splitfail=1;break;
      }
    }
  }
 

  if(!splitfail){// there are splits: compute honest splits
    for(s=0;s<r;s++){
      if(cflst[s]<0){
	if(*level==1)
	  ns=honestsplits(s,numv,explst,cflst,r,explst1,cflst1,r1,explst2,cflst2,r2,used1,used2,&needcount1,needs1,neededby1,&needcount2,needs2,neededby2,&needexcl0,q);
	else
	  ns=honestsplits(s,numv,explst,cflst,r,explst1,cflst1,r1,explst2,cflst2,r2,used1,used2,&needcount1,needs1,neededby1,&needcount2,needs2,neededby2,&needexcl0,1);
	if(ns<=0 && !q){// no honest splits
	  cout << "no honest splits but still splittable (Bigsq1):\n   " << cflst[s] << ", "; printvec1(explst[s],numv);
	}
      }
    }
  }


  // first create simple squares...?

  simpsq=0;

  // Only on the first pass, assuming that we didn't have a fundamental 
  // split failure (could this even happen?), check if there is an honest 
  // forced clique to remove and if so, choose the largest such (final flag 
  // in "honestforced") and store it for removal. We then continue as 
  // normal, finding further cliques.

  if(*level<=1 && !splitfail){
    used1tmp=(int *)malloc((size_t) (r1*sizeof(int)));
    used2tmp=(int *)malloc((size_t) (r2*sizeof(int)));

    ff=honestforced(explst, cflst, numv, r, needexcl, needs, neededby, explst1,  cflst1, r1, explst2, cflst2, r2, used1, used2, used1tmp, used2tmp, 1, 1, q);
    //cout << "ff = " << ff << endl;
 
    /* if(clqsremleavessplits(tmp, explst1, cflst1, r1, used1tmp, explst2, cflst2, r2, used2tmp, numv)) */
    /*   fprintf(stderr, "removal allowed.\n"); */
    /* else */
    /*   fprintf(stderr, "removal not allowed.\n"); */
    /* exit(0); */

    if(ff!=0){ // If there was an honest forced clique. 
      flg1=1;
      while(flg1){
  	flg1=0;
  	for(s=0;s<r1;s++){
  	  if(used1tmp[s]){
  	    if(!q){fprintf(stderr, "clique of type 1\n");}
  	    flg1=used1tmp[s];break;
  	  }
  	}
  	if(flg1){
  	  (*clq1)++;
  	  for(s=0;s<r1;s++){
  	    if(used1tmp[s]==flg1){
  	      used1[s]=(*clq1);
  	      used1tmp[s]=0;
  	    }
  	  }
  	}
      }
      flg2=1;
      while(flg2){
  	flg2=0;
  	for(s=0;s<r2;s++){
  	  if(used2tmp[s]){
  	    if(!q){fprintf(stderr, "clique of type 2\n");}
  	    flg2=used2tmp[s];break;
  	  }
  	}
  	if(flg2){
  	  (*clq2)++;
  	  for(s=0;s<r2;s++){
  	    if(used2tmp[s]==flg2){
  	      used2[s]=(*clq2);
  	      used2tmp[s]=0;
  	    }
  	  }
  	}
      }
      if(!q){cout << "largest forced...\n" << ff << endl;}
      simpsq=1;
    }
    free((char*)used1tmp);free((char*)used2tmp);
  }

  // Question: How can we ever get here (*level) should always be >=1. 
  // Is this redundant legacy code?

  if(!simpsq && *level<=0 && !splitfail){//if(*level==1){

    for(s=0;s<r1;s++){
      if(needcount1[s]==1 && needcount1[neededby1[s][2]]==1){
	s1=neededby1[s][2];
	if(!q){cout << "pair used only once.\n";}
	if(cflst1[s]*cflst1[s1]<0 && !used1[s] && !used1[s1]){
	  (*clq1)++;
	  if(!q){cout << "new basic cLique(1): " << (*clq1) << endl;}
	  used1[s]=(*clq1);used1[s1]=(*clq1);simpsq=1;
	}
	/* else if(!used1[s] && !used1[s1]){ */
	/*   cout << "incompatible square (1) found.\n"; */
	/*   used1[s]=-1;used1[s1]=-1; */
	/*   ifree(explst,r); free((char*)cflst);  */
	/*   ifree(needs1, r);ifree(neededby1, r1); */
	/*   free((char *)needcount1);free((char *)needcount2); */
	/*   ifree(needs2, r);ifree(neededby2, r2); */
	/*   free((char *)needexcl); */
	/*   return 0; */
	/* } */
      }
    }
    for(s=0;s<r2;s++){
      if(needcount2[s]==1 && needcount2[neededby2[s][2]]==1){
	s1=neededby2[s][2];
	if(!q){cout << "pair used only once.\n";}
	if(cflst2[s]*cflst2[s1]<0 && !used2[s] && !used2[s1]){
	  (*clq2)++;
	  if(!q){cout << "new basic cLique(2): " << (*clq2) << endl;}
	  used2[s]=(*clq2);used2[s1]=(*clq2);simpsq=1;
	}
	/* else if(!used2[s] && !used2[s1]){ */
	/*   cout << "incompatible square (2) found.\n"; */
	/*   used2[s]=-1;used2[s1]=-1; */
	/*   ifree(explst,r); free((char*)cflst);  */
	/*   ifree(needs1, r);ifree(neededby1, r1); */
	/*   free((char *)needcount1);free((char *)needcount2); */
	/*   ifree(needs2, r);ifree(neededby2, r2); */
	/*   free((char *)needexcl); */
	/*   return 0; */
	/* } */
      }

    }
  }

  if(simpsq){//success
    ifree(explst,r); free((char*)cflst); 
    ifree(needs, r);ifree(neededby, r);
    ifree(needs1, r);ifree(neededby1, r1);
    ifree(needs2, r);ifree(neededby2, r2);
    free((char *)needcount);free((char *)needcount1);free((char *)needcount2);
    free((char *)needexcl);free((char *)needexcl0);
    return 0; 
  }
  
  k=0;

  failed1=binisearch(explst,r,lastsplitfail,numv);

  // it is possible to fail because no more proper cliques can be found
  // even though nothing actually fails to split

  // findpairsplit1 does dummy runs of modcliques (create new clique or add 
  // to clique or merge cliques) to decide how best to deal with the negative 
  // monomial "failed1". used1 and used2 storing clique information aren't 
  // actually modified, but the monomial pair (ul and ur) is stored and 
  // modcliques is called with this information. 

  if(spf && failed1>=0){ // spf tells us to seed using the last split failure
    if(!q){cout << "Seeding (BigSq1): using last fundamental split failure to seed.\n";}
    k=findpairsplit1(tmpfull, explst, failed1, &ul, &ur, explst1, explst2, cflst1, cflst2, r1, r2, used1, used2, *clq1, *clq2, numv, 1, q);
  }
  else if(splitfail){
    if(!q){cout << "Seeding (BigSq1): split failure, finding seed.\n";}
    k=findpairsplit1(tmpfull, explst, failed, &ul, &ur, explst1, explst2, cflst1, cflst2, r1, r2, used1, used2, *clq1, *clq2, numv, 1, q);
  }
  else{// first look for totally new pairs (new cliques), then half new pairs (clique building), then old pairs (clique merging)
    if(!q){cout << "Seeding: no split failure, random allowed pair.\n";}
    k=randseedpair(explst, explst1, explst2, cflst, cflst1, cflst2, r, r1, r2, used1, used2, numv, needs1, needs2, &ul, &ur,2);
    if(!k)
      k=randseedpair(explst, explst1, explst2, cflst, cflst1, cflst2, r, r1, r2, used1, used2, numv, needs1, needs2, &ul, &ur,0);
    if(!k)
      k=randseedpair_old(explst, explst1, explst2, cflst, cflst1, cflst2, r, r1, r2, used1, used2, numv, needs, needexcl, &ul, &ur);
    /* k=firstseedpair(explst, explst1, explst2, cflst, cflst1, cflst2, r, r1, r2, used1, used2, numv, needs1, needs2, needexcl0, &ul, &ur); */
    /* if(!k) */
    /*   k=firstseedpair_old(explst, explst1, explst2, cflst, cflst1, cflst2, r, r1, r2, used1, used2, numv, needs, needexcl, &ul, &ur); */
  }

  if(!k){ 
    // This occurs if findpairsplit1 fails. If it is run with the flag 
    // "best", then this means that no splitting of the bad monomial 
    // improved the negative count from whatever baseline was chosen. 

    fprintf(stderr, "Something has gone wrong (BigSq1): k=%d.\n", k);
    if(!q){listprint(cflst, explst, numv, r);}
    for(s=0;s<r;s++){
      if(cflst[s]<0){
    	if(numsplits(s,numv,explst,cflst,r,&needcount,&needexcl,needs,neededby,0)<0.0001){// not able to split any more
    	  if(!q){cout << "failure to split(Bigsq1):\n   " << cflst[s] << ", "; printvec1(explst[s],numv);}
    	}
      }
    }
    /* k=delrandclique(used1,r1,*clq1,used2,r2,*clq2,&clqtmp); */
    /* if(k) */
    /*   cout << "deleting random clique(" << k << "): " << clqtmp << endl; */
    fprintf(stderr, "Seems to have failed.\n");
    ifree(explst,r); free((char*)cflst); 
    ifree(needs, r);ifree(neededby, r);  
    ifree(needs1, r);ifree(neededby1, r1);  ifree(needs2, r);ifree(neededby2, r2);
    free((char *)needcount);free((char *)needcount1);free((char *)needcount2);
    free((char *)needexcl);free((char *)needexcl0);
    return -1;
    //exit(0);
  }
  // This is where we actually update the cliques using ul and ur
  else 
    modcliques(k, used1, r1, used2, r2, ul, ur, clq1, clq2, q);

  ifree(explst,r); free((char*)cflst); 
  ifree(needs, r);ifree(neededby, r);  
  ifree(needs1, r);ifree(neededby1, r1);  ifree(needs2, r);ifree(neededby2, r2);
  free((char *)needcount);free((char *)needcount1);free((char *)needcount2);
  free((char *)needexcl);free((char *)needexcl0);

  return 0; // continue
}



// We remove squares (remsqs) and find squares (BigSq1) to remove. 
// Only the final remsqs and the final BigSq1 actually do the successful 
// removal, each prior failure serves to build up the list of squares to 
// use, held via used1 and used2. Remsqs does not alter used1 and used2. 
// BigSq1 on the other hand works both with the polynomial output by remsqs 
// and the full original one. It first tries a direct approach (J2pI_seq) 
// on the output polynomial. Assuming this fails, then the last failed 
// monomial is stored. The ways of dealing with this (new square, build 
// an existing square, merge squares) are examined by looking at the 
// original polynomial, and the one resulting in the "best" outcome from 
// the point of view of removing negative monomials is chosen. This means 
// we have to do dummy runs on constructing, building or merging squares 
// to see what the total effect is on the polynomial. With the data stored 
// in used1 and used2, we return to remsqs, remove all the squares so far, 
// and pass the new poly back to BigSq1.  

int BigSqWrap1(matrix J, int n, int numv, bool q){
  int **explst1,*cflst1,**explst2,*cflst2;
  int *used1,*used2;
  long r1,r2;
  ex tmp=0,tmp1,tmp2,ex1, extot,extot2,best;
  int level=0;
  int maxtries=20000,tries=0,iters;
  int clq1=0,clq2=0;
  ex sqs=0, finalpos;
  int suc=0;

  if(!q){cerr << "Entering BigSqWrap1" << endl;}fflush(stderr);
  ReImJpiI(J, n, numv, &explst1, &cflst1, &r1, &explst2, &cflst2, &r2,1);
  if(!q){cerr << "Completed ReImJpiI" << endl;}fflush(stderr);
  tmp=lsts_to_poly(explst1, cflst1, r1, explst2, cflst2, r2, numv+1);

  //  tmp=ex_J2pI_simphom(J, n, numv);
  if(!q){cerr << "tmp =" << tmp << endl;}
  sqs=0;

  used1=(int *)malloc((size_t) (r1*sizeof(int)));
  used2=(int *)malloc((size_t) (r2*sizeof(int)));

  if(!q){
    //Real and imaginary parts of det(J+iI)
    cout << "The two lists: \n";
    listprint(cflst1, explst1, numv+1, r1);
    cout << endl;
    listprint(cflst2, explst2, numv+1, r2);
  }
  iters=0;

  // Question. Is there anything to be gained by running this whole loop more than once, namely resetting the level

  while(iters<1){
    iters++;level=0;
    inittozero(used1, r1);
    inittozero(used2, r2);

    tries=0;clq1=0;clq2=0;

    while(tries<maxtries){
      tries++;
      if(!q){cout << "tries = " << tries << endl;}
      sqs=0;
      tmp2=remsqs(tmp, &sqs, explst1, cflst1, r1, explst2, cflst2, r2, numv+1, used1, used2, clq1, clq2, q);
      if((suc=BigSq1(tmp2, tmp, &sqs, &finalpos, n, numv+1, explst1, cflst1, used1, r1, explst2, cflst2, used2, r2, &clq1,&clq2,&level,q)))
	//if(BigSq1_old(tmp2, n, numv+1, explst1, cflst1, used1, r1, explst2, cflst2, used2, r2, &clq1,&clq2,&level))
	break;
    }
    if(!q){fprintf(stderr, "Final BigSq1 level = %d\n", level);}
    if(suc==1 && !q)
      cout << "expand(" << tmp << "-(" << sqs << ")-("<< finalpos << "));" << endl;
    //cout << "expand(" << tmp << "-("<< finalpos << "));" << endl;
    //cout << "final = " << finalpos << endl;

    if(tries==maxtries)
      fprintf(stderr, "Maximum allowed tries reached.\n");
  }

  ifree(explst1,r1); ifree(explst2,r2);  
  free((char*)cflst1); free((char*)cflst2);
  free((char*)used1); free ((char*)used2);
  return suc;
}



//Instead of BigSqWrap1
//Is det(J^2+I) positive on the positive orthant? (it is always definitely
//negative, so we are trying to confirm strict positivity)
//return codes zero (couldn't confirm positive) and 1 (definitely positive)
int J2pIwrap(matrix J, int nlen, int numv, int maxpppdegree, int *pppdegused, int debug){
  int **explst1,*cflst1,**explst2,*cflst2;
  long r1,r2;
  ex tmp,tmp1;
  ex tmpr,tmpi;//real and imaginary parts
  int numvinit, numvfinal, deg,allflg,ret1,retr,reti;
  int debugfull=(debug<=0)?0:debug-1;
  (*pppdegused)=-1;

  if(debug){cerr << "Entering J2pIwrap.\n";}

  tmpr=ReJpiI(J,nlen,numv,1);
  tmpi=ImJpiI(J,nlen,numv,1);

  //Check polys individually
  if((reti=ispospoly(expand(tmpi), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2 || reti==-2 || (retr=ispospoly(expand(tmpr), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull))==2 || retr==-2){
    if(debug){cerr << "One polynomial is positive so SOS definitely positive. Exiting J2pIwrap.\n";}
    return 1;
  }

  //Do SOS
  ReImJpiI(tmpr, tmpi, numv, &explst1, &cflst1, &r1, &explst2, &cflst2, &r2,1);
  tmp=lsts_to_poly(explst1, cflst1, r1, explst2, cflst2, r2, numv+1);
  ifree(explst1,r1); ifree(explst2,r2);free((char*)cflst1); free((char*)cflst2);

  tmp1=polyhomsimp(expand(tmp), &numvinit, &numvfinal, &deg, 1, 0, debugfull);
  if(debug){cerr << "Examining the polynomial: " << tmp1 << endl;}

  if(canbezero(tmp1, numvfinal, debugfull)){
    if(debug){cerr << "The polynomial can definitely take the value zero.\n";}
    return 0;
  }

  ret1=ispospolynonneg(tmp1, numvfinal, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull);

  if(ret1==2){  
    if(debug){cerr << "The polynomial is definitely positive. Exiting J2pIwrap.\n";}
    return 1;
  }

  if(ret1==0 && debug)
    cerr << "Definitely can be zero. Exiting J2pIwrap.\n";
  else if(ret1==1 && debug)
    cerr << "Couldn't determine if the polynomial is strictly positive on the positive orthant. Exiting J2pIwrap.\n";
  else if(ret1==-4 && debug)
    cerr << "This polynomial appears not to be nonnegative on the positive orthant. This algorithm should only be used for polynomials known to be nonnegative on the positive orthant. Exiting J2pIwrap.\n";
  return 0;

}



bool check3mat(int **J){
  if(!hasposdiag(J,3))
    return 0;
  if(J[0][1]*J[1][0]>0 && (J[0][0]*J[1][1]-J[0][1]*J[1][0]!=0))
    return 0;
  if(J[0][2]*J[2][0]>0 && (J[0][0]*J[2][2]-J[0][2]*J[2][0]!=0))
    return 0;
  if(J[1][2]*J[2][1]>0 && (J[1][1]*J[2][2]-J[1][2]*J[2][1]!=0))
    return 0;
  if(J[0][2]*J[2][1]*J[1][0]>0 && (J[0][0]*J[1][1]*J[2][2]-J[0][2]*J[2][1]*J[1][0]!=0))
    return 0;
  if(J[0][1]*J[1][2]*J[2][0]>0 && (J[0][0]*J[1][1]*J[2][2]-J[0][1]*J[1][2]*J[2][0]!=0))
    return 0;
  return 1;
}

//
// Try to rule out the possibility of Hopf bifurcation
// Return codes: 1 means Hopf bifurcation ruled out
// 0 means no such verdict
// The di6 file hold *irreversible* reactions
//
// filters: NULL means all checks
// ANULL means all non-MA checks
// JsquaredP0, i.e., is J^2 a P0 matrix?
// Jcomp2P0, i.e., is J^{2] a P0 matrix
// Jcomp2nonsing, i.e., is J^[2] nonsingular
// J2pInonsing, i.e., is J^2+I nonsingular
// MAJ2pInonsing, i.e., is J^2+I nonsingular for the mass action Jacobian
// Return codes: 3=dynamically trivial; 2=no Hopf with GK; 1=no Hopf with MA; 0: no result.

int Hopfforbid(int **imatSi, int **imatSl, int nlen, int mlen, const char *filter, int effort, int *pppdegused, int debug){
  int numv=0;
  //int numvMA;
  matrix J, JMA, JMAeq, J2, exV;
  int ret;
  //int ret1;
  int pos2flag=0;
  matrix rhs;
  int **Q=NULL;
  matrix QX;
  int maxpppdeg=effort-1; //maxpppdeg<0 means don't use csdp
  int debugfull=(debug<=0)?0:debug-1;
  (*pppdegused)=-1;

  if(debug){fprintf(stderr, "\n###Entering Hopfforbid.\n");}

  if(!filter)
    filter="ALL";//set pointer to default

  if(strcmp(filter,"ALL") && strcmp(filter,"MAonly") && strcmp(filter,"GKonly") && strcmp(filter,"JsquaredP0") && strcmp(filter,"MAJsquaredP0") && strcmp(filter,"Jcomp2P0") && strcmp(filter,"MAQnegsemidef") && strcmp(filter,"MAJcomp2P0") && strcmp(filter,"Jcomp2nonsing") && strcmp(filter,"MAJcomp2nonsing") && strcmp(filter,"J2pInonsing") && strcmp(filter,"MAJ2pInonsing")){
    fprintf(stderr, "ERROR in HopfForbid: \"%s\" is not an allowed filter. EXITING.\n", filter);
    exit(0);
  }

  if(!hasposrkervec(imatSi,nlen,mlen,1)){
    if(debug){fprintf(stderr, "Hopf bifurcation automatically ruled out, since the system is dynamically trivial.\n");}
    return 3;
  }


  //
  //General kinetics calculations
  //
  if(debug && (!strcmp(filter, "ALL") || !strcmp(filter, "GKonly"))){fprintf(stderr, "Doing general kinetics tests.\n");}
  exV=Vmatfrompatmat(imatSl, nlen, mlen, &numv);
  J=multABT(imatSi, exV, nlen, mlen);//minus Jacobian
  if(debugfull){fprintf(stderr, "Jacobian matrix:\n");printmat(J, nlen, nlen);}

  if(nlen <=4 && (!strcmp(filter, "ALL") || !strcmp(filter, "GKonly") || !strcmp(filter, "Jcomp2P0"))){ // no filter or only check if J^[2] is P0
    if(AdComp2isP0(J, nlen, maxpppdeg, debugfull)==1){
      if(debug){cerr << "Jcomp2P0: -J^[2] is a P0 matrix (so Hopf bifurcation ruled out).\n";}
      return 2;
    }
    else{
      if(debugfull || (debug && !strcmp(filter, "Jcomp2P0"))){cerr << "Failed to confirm if -J^[2] is a P0 matrix.\n";}
    }
  }
  if(!strcmp(filter, "Jcomp2P0"))
    return 0;

  if(nlen<=4 && (!strcmp(filter, "ALL") || !strcmp(filter, "GKonly") || !strcmp(filter, "Jcomp2nonsing"))){ // no filter or only check if J^[2] is nonsingular
    if((pos2flag=AdComp2DetPos(J, nlen, numv, maxpppdeg, pppdegused, debugfull))==2 || pos2flag==-2 || (nlen<=3 && pos2flag==0)){
      if(debug){cerr << "J^[2] has strictly signed determinant, ruling out Hopf bifurcation. Result: " << pos2flag << endl;}
      if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
      return 2;
    }
    else if(debugfull || (debug && !strcmp(filter, "Jcomp2nonsing"))){
      if(pos2flag==-3)
	cerr << "det(J^[2]) can definitely switch signs.\n";
      else if(pos2flag==0)
	cerr << "det(J^[2]) is identically zero.\n";
      else
	cerr << "failed to confirm if either det(J^[2]) has strictly signed determinant or can take all signs. Return code: " << pos2flag << endl;
    }
  }
  if(!strcmp(filter, "Jcomp2nonsing"))
    return 0;

  if(!strcmp(filter, "ALL") || !strcmp(filter, "GKonly") || !strcmp(filter, "JsquaredP0")){ // no filter or only check if J^2 is P0
    //here, checking sign symmetry seems to be faster than direct confirmation

    /* J2=multAB(J, J, nlen, nlen);//Jacobian squared */
    /* if(debug){printmat1(J2, nlen, nlen);} */
    /* if(isP0matorth(J2, nlen, 1, maxpppdeg, pppdegused, debug)==1){//third argument means skip the top dimension */
    if(signsym(J, nlen, maxpppdeg,debugfull)){
      if(debug){cerr << "JsquaredP0: J^2 is a P0 matrix, ruling out Hopf bifurcation.\n";}
      return 2;
    }
    if(debugfull || (debug && !strcmp(filter, "JsquaredP0"))){cerr << "Failed to confirm if J^2 is a P0 matrix.\n";}
  }
  if(!strcmp(filter, "JsquaredP0"))
    return 0;


  if(!strcmp(filter, "ALL") || !strcmp(filter, "GKonly") || !strcmp(filter, "J2pInonsing")){ // no filter or only check if J^2+I is nonsingular
    if((ret=J2pIwrap(J,nlen,numv,maxpppdeg,pppdegused,debugfull))==1){ // 1 for success, 0 for failure. 
      if(debug){cerr << "J2pInonsing: Using J2pIwrap, det(J^2+I) is positive on the positive orthant, ruling out Hopf bifurcation.\n\n";}
      return 2;
    }
    if(debugfull || (debug && (!strcmp(filter, "GKonly") || !strcmp(filter, "J2pInonsing")))){cerr << "J2pIwrap could not determine if det(J^2+I) is positive on the positive orthant.\n\n";}
  }
  if(!strcmp(filter, "GKonly") || !strcmp(filter, "J2pInonsing"))
    return 0;

  if(debug && !strcmp(filter, "GKonly")){fprintf(stderr, "Failed to rule out Hopf bifurcations in this system for general kinetics.\n");}

  //This is commented out because if we are going to do MA, then why not focus on equilibrium. But there may be good reasons to do MA more generally
  /* //Mass action Jacobian (not necessarily at equilibrium) */
  /* //rhs vector field (mass action) */
  /* rhs=MARHSfromlstoich(imatSi, imatSl, nlen, mlen); */
  /* if(debug){fprintf(stderr, "Mass action vector field:\n");cerr << rhs << endl;} */
  /* exV=VmatMAfromlstoich(imatSl, nlen, mlen, &numvMA); */
  /* //if(debug){printmat1(exV, nlen, mlen);} */
  /* JMA=multABT(imatSi, exV, nlen, mlen); */
  /* if(debug){fprintf(stderr, "Mass action Jacobian matrix:\n");printmat1(JMA, nlen, nlen);} */

  /* if(!filter || !strcmp(filter, "ALL") || !strcmp(filter, "GKonly") || !strcmp(filter, "J2pIMAnonsing")){ // no filter or only check if J^2+I is nonsingular */
  /*   if(!rapid && debug){cerr << "Trying Mass Action.\n";} */

  /*   if(!rapid){ */
  /* 	if((ret1=J2pIwrap(JMA,nlen,numvMA,maxpppdeg,pppdegused,debug))==1){ // 1 for success, -1 for failure.  */
  /* 	  cerr << "J2pIMAnonsing: Using J2pIwrap, MA: det(J^2+I) is positive on the positive orthant.\n\n"; */
  /* 	  return 1; */
  /* 	} */
  /* 	if(debug){cerr << "J2pIwrap (MA): Could not determine if det(J^2+I) is positive on the positive orthant.\n\n";} */
  /*   } */
  /* } */
  /* if(!strcmp(filter, "J2pIMAnonsing")) */
  /*   return 0; */


  //
  // MA calculations at equilibrium
  //
  if(debug){fprintf(stderr, "Doing mass action tests.\n");}

  JMAeq=reacJMAeq(imatSi, imatSl, nlen, mlen, &Q, &QX, &numv);//Q is the reduced form of the Jacobian matrix when rankS==mlen-1
  if(debugfull){fprintf(stderr, "MA Jacobian matrix:\n");printmat(JMAeq, nlen, nlen);}

  if(!strcmp(filter, "ALL") || !strcmp(filter, "MAonly") || !strcmp(filter, "MAQnegsemidef")){
    if(debugfull || (debug && !strcmp(filter, "MAQnegsemidef"))){fprintf(stderr, "Checking if JMAeq is similar to a negative semidefinite matrix.\n");}
    
    if(reacQMAposdef(Q, QX, nlen, 0, 0, maxpppdeg, debugfull)){
      if(debug){cerr << "MAQnegsemidef: JMAeq is similar to a neg semidef matrix.\n";}
      if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
      return 1;
    }
    if(debugfull || (debug && !strcmp(filter, "MAQnegsemidef"))){cerr << "Failed to confirm if JMAeq is similar to a neg semidef matrix.\n";}
  }
  if(!strcmp(filter, "MAQnegsemidef")){
    if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
    return 0;
  }

  if(nlen<=4 && (!strcmp(filter, "ALL") || !strcmp(filter, "MAonly") || !strcmp(filter, "MAJcomp2P0"))){
    if(debugfull || (debug && !strcmp(filter, "MAJcomp2P0"))){fprintf(stderr, "Checking if -JMAeq^[2] is a P0 matrix.\n");}
    if(Q && nlen==3 && check3mat(Q)){
      if(debug){cerr << "Dimension 3 and e-o-s condition on Q satisfied (so JMAeq^[2] is a P0 matrix).\n";printmat(Q, nlen, nlen);}
      free_imatrix(Q, 0, nlen-1, 0, nlen-1);
      return 1;
    }

    if(AdComp2isP0(JMAeq, nlen, maxpppdeg, debugfull)==1){
      if(debug){cerr << "MAJcomp2P0: -JMAeq^[2] is a P0 matrix, ruling out Hopf bifurcation.\n";}
      if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
      return 1;
    }
    if(debugfull || (debug && !strcmp(filter, "MAJcomp2P0"))){cerr << "Failed to confirm if -JMAeq^[2] is a P0 matrix.\n";}
  }
  if(!strcmp(filter, "MAJcomp2P0")){
    if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
    return 0;
  }

  if(nlen <=4 && (!strcmp(filter, "ALL") || !strcmp(filter, "MAonly") || !strcmp(filter, "MAJcomp2nonsing"))){ // no filter or only check if J^[2] is nonsingular
    //fourth argument = maxpppdegree negative means don't use csdp
    if((pos2flag=AdComp2DetPos(JMAeq, nlen, numv, maxpppdeg, pppdegused, debugfull))==2 || pos2flag==-2 || (nlen<=3 && pos2flag==0)){
      if(debug){cerr << "JMAeq^[2] has strictly signed determinant. Result: " << pos2flag << endl;}
      if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
      return 1;
    }
    else if(debugfull || (debug && !strcmp(filter, "MAJcomp2nonsing"))){
      if((pos2flag==-3 || pos2flag==0)){
	if(pos2flag==-3)
	  cerr << "det(J^[2]) can definitely switch signs.\n";
	else
	  cerr << "det(J^[2]) is identically zero.\n";
      }
      else{
	cerr << "failed to confirm if either det(J^[2]) has strictly signed determinant or can take all signs. Return code: " << pos2flag << endl;
      }
    }
  }
  if(!strcmp(filter, "MAJcomp2nonsing")){
    if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
    return 0;
  }


  J2=multAB(JMAeq, JMAeq, nlen, nlen);//Jacobian squared
  if((!strcmp(filter, "ALL") || !strcmp(filter, "MAonly") || !strcmp(filter, "MAJsquaredP0"))){
    if(debugfull){fprintf(stderr, "Checking if JMAeq^2 is a P0 matrix.\n");printmat(J2, nlen, nlen);}
    if(Q && signsym(Q,nlen)){
      if(debug){cerr << "Q is sign symmetric (so JMAeq^2 is a P0 matrix).\n";printmat(Q, nlen, nlen);}
      free_imatrix(Q, 0, nlen-1, 0, nlen-1);
      return 1;
    }
    else if(!Q && signsym(QX,nlen,maxpppdeg,debug)){
      if(debug){cerr << "QX is sign symmetric (so JMAeq^2 is a P0 matrix).\n";printmat(QX, nlen, nlen);}
    }
    if(debugfull || (debug && !strcmp(filter, "MAJsquaredP0"))){cerr << "Failed to confirm if JMAeq^2 is a P0 matrix.\n";}
  }
  if(!strcmp(filter, "MAJsquaredP0")){
    if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
    return 0;
  }

  if(!strcmp(filter, "ALL") || !strcmp(filter, "MAonly") || !strcmp(filter, "MAJ2pInonsing")){ // no filter or only check if J^2+I is nonsingular
    if(debugfull){fprintf(stderr, "Checking if det(JMAeq^2+I) is positive on the positive orthant.\n");}
    if((ret=J2pIwrap(JMAeq,nlen,numv,maxpppdeg,pppdegused,debugfull))==1){ //num variables (pre-homogenisation) = nlen - same as dimension
      if(debug){cerr << "MAJ2pInonsing: Using J2pIwrap, det(JMAeq^2+I) is positive on the positive orthant.\n\n";}
      if(Q){free_imatrix(Q, 0, nlen-1, 0, nlen-1);}
      return 1;
    }
    if(debugfull || (debug && !strcmp(filter, "MAJ2pInonsing"))){cerr << "J2pIwrap could not determine if det(JMAeq^2+I) is positive on the positive orthant.\n\n";}
  }

  if(Q)
    free_imatrix(Q, 0, nlen-1, 0, nlen-1);

  if(debug){
    if(!strcmp(filter, "MAonly"))
      fprintf(stderr, "Failed to rule out Hopf bifurcations in this system for mass action.\n");
    else
      fprintf(stderr, "Failed to rule out Hopf bifurcations in this system for general kinetics or mass action.\n");
  }
  return 0;
}

//Overloading. Input is reactions as a (human readable) string
int HopfForbid(char *str, const char *filter, int effort, int *pppdegused, int debug){
  int **imatS=NULL, **imatV=NULL, **imatSi, **imatSl, **stoichl, **stoichr;
  int ret, mlen=0, nlen=0, numSi, allrev, allgood;
  char **chems;
  bool haszerocomplex=0;
  

  if(filter && strcmp(filter,"ALL") && strcmp(filter,"MAonly") && strcmp(filter,"GKonly") && strcmp(filter,"JsquaredP0") && strcmp(filter,"MAJsquaredP0") && strcmp(filter,"Jcomp2P0") && strcmp(filter,"MAQnegsemidef") && strcmp(filter,"MAJcomp2P0") && strcmp(filter,"Jcomp2nonsing") && strcmp(filter,"MAJcomp2nonsing") && strcmp(filter,"J2pInonsing") && strcmp(filter,"MAJ2pInonsing")){
    fprintf(stderr, "ERROR in HopfForbid: %s is not an allowed filter. EXITING.\n", filter);
    exit(0);
  }

  //The flag at the end causes the reactions to be printed out
  if(!getallreacs(str, &imatS, &imatV, &imatSi, &imatSl, &stoichl, &stoichr, &chems, &haszerocomplex, &nlen, &mlen, &numSi, &allrev, &allgood, debug)){
    fprintf(stderr, "ERROR in HopfForbid: couldn't read reactions from string %s. EXITING.\n", str);
    exit(0);
  }

  //important that imatSl is *minus* the left stoich matrix

  ret = Hopfforbid(imatSi, imatSl, nlen, mlen, filter, effort, pppdegused, debug);

  free_imatrix(imatS, 0, nlen-1, 0, mlen-1);
  free_imatrix(imatV, 0, nlen-1, 0, mlen-1);
  free_imatrix(imatSi, 0, nlen-1, 0, numSi-1);
  free_imatrix(imatSl, 0, nlen-1, 0, numSi-1);
  free_imatrix(stoichl, 0, nlen-1, 0, numSi-1);
  free_imatrix(stoichr, 0, nlen-1, 0, numSi-1);
  freearraydat(chems, nlen);

  return ret;
}

//Overloading: input is PN AM
int HopfForbid(int **AM, int n, int m, const char *filter, int effort, int *pppdegused, int debug){
  int ret;
  int **S, **Sl;
  bool minus=1;
  AMtoSSl(AM,n,m,minus,&S,&Sl);
  ret=Hopfforbid(S, Sl, n, m, filter, effort, pppdegused, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//Overloading: input is di6 string
int HopfForbid(char *di6, int n, int m, const char *filter, int effort, int *pppdegused, int debug){
  int ret;
  int **S, **Sl;
  bool minus=1;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  ret=Hopfforbid(S, Sl, n, m, filter, effort, pppdegused, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}


//overloading (set filter to "ALL")
int HopfForbid(char *di6, int n, int m, int effort, int *pppdegused, int debug){
  return HopfForbid(di6, n, m, "ALL", effort, pppdegused, debug);
}



// Simplified processing
// check for hopf bifurcation capability
// in a list of crns stored in di6 format

// arguments: 
// input file with one bimolecular crn in di6 format per line
// output file for crns with possible Hopf bifurcation (or NULL)
// output file for crns forbidding Hopf bifurcation (or NULL)
// number of species
// number of reactions
// filter (see HopfForbid above)
// verbose? 1 for verbose output


// ruleouthopf returns 1 if we can be confident that Hopf bifurcation
// is forbidden, and returns 0 otherwise

void ruleouthopf(const char *fname, const char *outfhopf, const char *outfnohopf, int n, int m, const char *filter, int effort, int debug){
  unsigned long i=0;
  int ret;
  FILE *fd, *fd1, *fdin;
  int maxl=0;
  unsigned long numl=numflines(fname, &maxl);
  char oneline[maxl];
  int pppdegused=-1;

  if(outfhopf && !(fd=fopen(outfhopf, "w"))){
    fprintf(stderr, "ERROR in ruleouthopf: could not open file %s for writing. EXITING.\n", outfhopf);
    exit(0);
  }
  if(outfnohopf && !(fd1=fopen(outfnohopf, "w"))){
    fprintf(stderr, "ERROR in ruleouthopf: could not open file %s for writing. EXITING.\n", outfnohopf);
    exit(0);
  }
  
  if(!(fdin = fopen(fname, "r"))){
    fprintf(stderr, "ERROR in ruleouthopf: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }

  fprintf(stderr, "checking %ld CRNs.\n", numl);
  //cout << "here\n";
  while(getline0(fdin, oneline, maxl) > 0 && !(iscomline(oneline))){
    //array2=genarraydat(fname, &tot);
    cout << i+1 << "/" << numl << endl;
    if(debug){cerr << i+1 << "/" << numl << endl;}
    ret=HopfForbid(oneline,n,m,filter,effort,&pppdegused,debug);
    if(outfhopf && !ret)
      fprintf(fd,"%s\n",oneline);
    if(outfnohopf && ret)
      fprintf(fd1,"%s\n",oneline);
    if(debug){
      cerr << "ret=" << ret;
      if(ret && pppdegused>=0)
	cerr << "\nUsed SDP to degree " << pppdegused << "\n";
      cerr << "\n******\n";
    }
    i++;
  }

  fclose(fdin);

  if(outfhopf)
    fclose(fd);
  if(outfnohopf)
    fclose(fd1);

  return;
}

//overloading (no filters - run all checks)
void ruleouthopf(const char *fname, const char *outfhopf, const char *outfnohopf, int n, int m, int effort, int debug){
  ruleouthopf(fname, outfhopf, outfnohopf, n, m, "NULL", effort, debug);
  return;
}

