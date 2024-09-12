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
#include "symbolic.h"
#include "convex.h"
#include "analysereacs.h"
#include "pos.h"
#include "hopf.h"
#include "isomorph.h"

//#include <iostream>

/* note that these are not the most efficient */
/* algorithms. Also there is very little error checking  */
/* (to maximise speed), and if the routines are used carelessly,  */
/* it is possible to leak memory.  */

/* Note that in the case of large sparse matrices, the algorithms */
/* presented here serve to highlight the need for graph-theoretic */
/* algorithms */

/* A permutation of [1, ..., n] is just an n-vector */

/* calculate the unsigned term in an n X n matrix corresponding */
/* to permutation tm */

int unsterm(int **imat, int *tm, int n){
  int i,tot=1;
  for(i=0;i<n;i++)
    tot*=imat[i][tm[i]];
  return tot;
}

/* the unsigned term in the determinant of a k X k submatrix  */
/* of n X m matrix imat, corresponding to permutation tm */

int unsterm1(int **imat, int n, int m, int *vec1, int *tm, int k){
  int i,j,tot=1;
  for(i=0;i<k;i++){
    if((j=imat[vec1[i]][tm[i]])!=0)
      tot*=j;
    else
      return 0;
  }
  return tot;
}

/* same as unsterm1, except that all positive entries in the  */
/* submatrix are replaced with zeros */

int unsterm2(int **imat, int n, int m, int *vec1, int *tm, int k){
  int i,tot=1;
  for(i=0;i<k;i++)
    tot*=min(imat[vec1[i]][tm[i]], 0);
  return tot;
}

/* check if the k X k minor of n X m matrix imat indexed by the pair */
/* (vec1, vec2) is sign singular */
/* pms is a matrix storing all permutations of vec2 */
/* generated with allperms1. The final digit of each */
/* vector in pms is the parity of the permutation */
int minorisSS(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  unsigned long i;
  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;
  for(i=0;i<fk;i++){
    if(unsterm1(imat,n,m,vec1,pms[i],k)) /* nonzero term */
      return 0;
  }
  return 1;
}

/* like minorisSS except doesn't first compute */
/* all the permutations, so should avoid memory failure */

int minorisSS2(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk){
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;

  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  vec=(int *)malloc((size_t) ((k)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((k)*sizeof(int)));

  for(i=0;i<k;i++){vec[i]=vec2[i];veclr[i]=-1;}

  while(flag){
    if(unsterm1(imat, n, m, vec1, vec,k)) /* nonzero term */
      return 0;
    flag=nextperm(&vec, &veclr, &par, k);
    if(t%1000000==0){
      fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100));
    }
    t++;
  }
  free((char *) vec);free((char *) veclr);
  return 1;
}

/* check if the k X k minor of n X m matrix imat indexed by the pair */
/* (vec1, vec2) is sign nonsingular or sign singular */
/* pms is a matrix storing all permutations of vec2 */
/* generated with allperms1. The final digit of each */
/* vector in pms is the parity of the permutation */
int minorisSNSSS(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  int tmp, tm1=0;
  unsigned long i;

  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  for(i=0;i<fk;i++){
    if((tmp=unsterm1(imat, n, m, vec1, pms[i],k))){ /* nonzero term */
      if(!tm1){tm1=pms[i][k]*tmp;} /* store first nonzero term */
      else if(tm1*pms[i][k]*tmp<0) /* oppositely signed terms */
	return 0;
    }
  }
  return 1;
}

int minorisSNSSS_rec(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;
  if(qualdetsubmat(imat, n, m, vec1, vec2, k)!=2)
    return 1;
  return 0;

}

/* Exactly like minorisSNSSS, but */
/* returns 1/-1 for SNS, 2 for SS and 0 for neither. */

int minorisSNSSS1(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  int tmp, tm1=0;
  unsigned long i;

  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 2;

  for(i=0;i<fk;i++){
    if((tmp=unsterm1(imat, n, m, vec1, pms[i],k))){ /* nonzero term */
      if(!tm1){tm1=pms[i][k]*tmp;} /* store first nonzero term */
      else if(tm1*pms[i][k]*tmp<0) /* oppositely signed terms */
	return 0;
    }
  }
  if(!tm1)
    return 2;
  if(tm1<0)
    return -1;
  return 1;
}



/* like minorisSNSSS except doesn't first compute */
/* all the permutations, so should avoid memory failure */
int minorisSNSSS2(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;
  //unsigned long fk=factorial(k);

  /* worth checking for a row or column of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<k;i++){
    vec[i]=vec2[i];
    veclr[i]=-1;
  }

  while(flag){

    /* for(i=0;i<n;i++){ */
    /*   if(veclr[i]==-1) */
    /* 	fprintf(stderr, "<%d  ", vec[i]+1); */
    /*   else */
    /* 	fprintf(stderr, "%d>  ", vec[i]+1); */
    /* } */

    if((tmp=unsterm1(imat, n, m, vec1, vec,k))){ /* nonzero term */
      if(!tm1){tm1=par*tmp;} /* store first nonzero term */
      else if(tm1*par*tmp<0){ /* oppositely signed terms */
	free((char *) vec);free((char *) veclr);
	return 0;
      }
    }
    flag=nextperm(&vec, &veclr, &par, k);
    /* if(t%1000000==0){ */
    /*   fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100)); */
    /* } */
    t++;
  }
  free((char *) vec);free((char *) veclr);
  return 1;
}

/* Exactly like the previous routine, but */
/* returns 1/-1 for SNS, 2 for SS and 0 for neither. */
int minorisSNSSS3(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;
  unsigned long fk=factorial(k);

  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 2;

  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<k;i++){
    vec[i]=vec2[i];
    veclr[i]=-1;
  }

  while(flag){

    /* for(i=0;i<k;i++){ */
    /*   if(veclr[i]==-1) */
    /* 	fprintf(stderr, "<%d  ", vec[i]+1); */
    /*   else */
    /* 	fprintf(stderr, "%d>  ", vec[i]+1); */
    /* } */
    /* fprintf(stderr, "\n"); */

    if((tmp=unsterm1(imat, n, m, vec1, vec,k))){ /* nonzero term */
      if(!tm1){tm1=par*tmp;} /* store first nonzero term */
      else if(tm1*par*tmp<0){ /* oppositely signed terms */
	free((char *) vec);free((char *) veclr);
	return 0;
      }
    }
    flag=nextperm(&vec, &veclr, &par, k);
    if(t%1000000==0){
      fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100));
    }
    t++;
  }
  free((char *) vec);free((char *) veclr);

  if(!tm1)
    return 2;
  if(tm1<0)
    return -1;
  return 1;

}

/* Check if the n times n matrix imat is sign nonsingular */
/* or sign singular */

int matrixisSNSSS(int **imat, int n){
  int vec1[n];
  unsigned long fn;
  int **pms;
  int i, retval=0;

  for(i=0;i<n;i++)
    vec1[i]=i;

  if(n>9){
    retval = minorisSNSSS3(imat, n, n, vec1, vec1, n);
  }
  else{
    fn=factorial(n);
    pms=allperms(n);

    retval = minorisSNSSS1(imat, n, n, vec1, vec1, n, fn, pms);
    free_imatrix(pms,0, fn-1, 0, n);
  }
  return retval;

}




//Check for reversible pairs, and create new reversible stoichiometric and pattern matrices
//mout is the number of reactions after reversification
int reversify(int **S, int **V, int n, int m, int ***Sout, int ***Vout, int *mout){
  int j,k;
  int pairs[m][2];
  int totpairs=0,curcol,curpair;
  bool used[m];
  inittozero(used,m);

  //Find cols with reverse
  for(j=0;j<m-1;j++){
    if(used[j])
      continue;
    for(k=j+1;k<m;k++){
      if(!used[k] && isreversecol(S, n, j, k) && signcompatcols(V, n, m, j, k, 1)){//or disjointcols? oppositely signed columns signcompatcols... 1
	pairs[totpairs][0]=j;pairs[totpairs][1]=k;
	totpairs++;
	used[k]=1;//only stores second element in each pair
	break;
      }
    }
  }

  (*mout)=m-totpairs;
  //Now update the matrices
  (*Sout)=imatrix(0,n-1,0,(*mout)-1);
  (*Vout)=imatrix(0,n-1,0,(*mout)-1);
  curcol=0;curpair=0;
  for(j=0;j<m;j++){//each column
    if(curpair<totpairs && pairs[curpair][0]==j){//first in a pair
      for(k=0;k<n;k++){//column subtract column only for V
	(*Sout)[k][curcol]=S[k][j];
	(*Vout)[k][curcol]=V[k][j]-V[k][pairs[curpair][1]];
      }
      curpair++;curcol++;
    }
    else if(!used[j]){//not second element in a pair
      for(k=0;k<n;k++){//copy columns
	(*Sout)[k][curcol]=S[k][j];
	(*Vout)[k][curcol]=V[k][j];
      }
      curcol++;
    }
  }

  if(curpair!=totpairs){
    fprintf(stderr, "ERROR in reversify. There should be %d pairs, but %d found. EXITING.\n", totpairs, curpair);
    exit(0);
  }

    
  return totpairs;
}


//Check for reactions with equal vectors and compatible patterns
//mout is the number of reactions after merging
int mergecols(int **S, int **V, int n, int m, int ***Sout, int ***Vout, int *mout){
  int j,k;
  int pairs[m][2];
  int totpairs=0,curcol,curpair;
  bool used[m];
  inittozero(used,m);

  //Find cols with reverse
  for(j=0;j<m-1;j++){
    if(used[j])
      continue;
    for(k=j+1;k<m;k++){
      if(!used[k] && isequalcol(S, n, j, k) && signcompatcols(V, n, m, j, k, 0)){//compatibly signed columns
	pairs[totpairs][0]=j;pairs[totpairs][1]=k;
	totpairs++;
	used[k]=1;//only stores second element in each pair
	break;
      }
    }
  }

  (*mout)=m-totpairs;
  //Now update the matrices
  (*Sout)=imatrix(0,n-1,0,(*mout)-1);
  (*Vout)=imatrix(0,n-1,0,(*mout)-1);
  curcol=0;curpair=0;
  for(j=0;j<m;j++){//each column
    if(curpair<totpairs && pairs[curpair][0]==j){//first in a pair
      for(k=0;k<n;k++){//column subtract column only for V
	(*Sout)[k][curcol]=S[k][j];
	(*Vout)[k][curcol]=V[k][j]+V[k][pairs[curpair][1]];
      }
      curpair++;curcol++;
    }
    else if(!used[j]){//not second element in a pair
      for(k=0;k<n;k++){//copy columns
	(*Sout)[k][curcol]=S[k][j];
	(*Vout)[k][curcol]=V[k][j];
      }
      curcol++;
    }
  }

  if(curpair!=totpairs){
    fprintf(stderr, "ERROR in mergecols. There should be %d pairs, but %d found. EXITING.\n", totpairs, curpair);
    exit(0);
  }

    
  return totpairs;
}





/* Does the k X k minor of n X m matrix imat indexed by */
/* vec1 and vec2 have a row or column of zeros? */
/* Note that a failure does *not* imply that the minor is */
/* not sign singular - to really check this you'd need to check  */
/* if all terms in the determinant expansion are zero, a much  */
/* lengthier task */

int minorhas0rc(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int flag, i=0, j=0;  

  /* fprintf(stderr, "\n###Entering minorhas0rc, k=%d\n", k); */
  /* printmat(imat,n,m); */
  /* printvec(vec1,k);printvec(vec2,k); */

  while(j<k){
    flag=0;i=0;
    while(i<k && !flag){
      if(imat[vec1[i]][vec2[j]])
	flag=1; /* not a column of zeros */
      i++;
    }
    if(!flag){
      /* fprintf(stderr, "here col...\n"); */
      /* printsubmat(imat, vec1, vec2, k, k); */
      return 1; /* found a column of zeros */
    }
    j++;
  }

  i=0;j=0;

  while(i<k){
    flag=0;j=0;
    while(j<k && !flag){
      //fprintf(stderr, "i,j=%d,%d: %d\n", i,j,imat[vec1[i]][vec2[j]]);
      if(imat[vec1[i]][vec2[j]])
	flag=1; /* not a row of zeros */
      j++;
    }
    if(!flag){
      /* fprintf(stderr, "here row...\n"); */
      /* printsubmat(imat, vec1, vec2, k, k); */
      return 1; /* found a row of zeros */
    }
    i++;
  }

  /* fprintf(stderr, "here no rc...\n"); */
  /* printsubmat(imat, vec1, vec2, k, k); */
  return 0;
}

/* exactly as the previous routine, except checking imat-  */
/* (imat with all positive entries replaced with zeros) */

int minorhas0rc1(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int flag, i=0, j=0;  

  while(j<k){
    flag=0;i=0;
    while(i<k && !flag){
      if(imat[vec1[i]][vec2[j]] <0)
	flag=1; /* not a column of zeros */
      i++;
    }
    if(!flag)
      return 1; /* found a column of zeros */
    j++;
  }

  i=0;j=0;

  while(i<k){
    flag=0;j=0;
    while(j<k && !flag){
      if(imat[vec1[i]][vec2[j]] <0)
	flag=1; /* not a row of zeros */
      j++;
    }
    if(!flag)
      return 1; /* found a row of zeros */
    i++;
  }

  return 0;
}


/* check if the k X k minor of n X m matrix imat indexed by  */
/* (vec1, vec2) is sign nonsingular or singular */
/* 1) calculate all permutations of vec2 */
/* 2) calculate the signed terms corresponding to vec1 and */
/*    each permutation of vec2 */
/* 3) Sum as we calculate to get the determinant */
/* 4) If determinant is zero then we are done */
/* 5) Otherwise multiply each term by the determinant to check if any */
/*    product is less than zero */

int minorisSNSsing(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  int dt=0, tmp, tm1=0;
  unsigned long i;
  int flag=0;

  /* worth checking for a column or row of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  for(i=0;i<fk;i++){ // only store nonzero terms
    tmp=unsterm1(imat, n, m, vec1, pms[i],k);
    if (tmp){/* nonzero term */
      if(!tm1) /* first nonzero term */
	tm1=pms[i][k]*tmp;
      else if(pms[i][k]*tmp*tm1 <0)
	flag=1; /* not SNS */

      dt+=pms[i][k]*tmp;
    }
  }

  if(!flag || !dt) /* flag=0 means SNS; dt=0 means singular */
    return 1;

  return 0;
}

/* like the previous routine, except  */
/* does not require large memory storage */
/* for all permutations */

int minorisSNSsing1(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int dt=0, i, par=1, flag=0,nxt=1;
  unsigned long t=0;
  /* unsigned long fk=factorial(k); */

  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k)) /* sign singular */
    return 2;

  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<k;i++){
    vec[i]=vec2[i];
    veclr[i]=-1;
  }

  while(nxt){

    /* for(i=0;i<k;i++){ */
    /*   if(veclr[i]==-1) */
    /* 	fprintf(stderr, "<%d  ", vec[i]+1); */
    /*   else */
    /* 	fprintf(stderr, "%d>  ", vec[i]+1); */
    /* } */
    /* fprintf(stderr, "\n"); */

    if ((tmp=unsterm1(imat, n, m, vec1, vec,k))){/* nonzero term */
      if(!tm1){tm1=par*tmp;} /* first nonzero term */
      else if(par*tmp*tm1 <0){
	flag=1; /* not SNS */
	//fprintf(stderr, "here\n");
      }

      dt+=par*tmp;
    }

    nxt=nextperm(&vec, &veclr, &par, k);
    /* if(t%1000000==0){ */
    /*   fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100)); */
    /* } */
    t++;
  }
  free((char *) vec);free((char *) veclr);

  if(!flag || !dt) /* flag=0 means SNS; dt=0 means singular */
    return 1;

  if(!tm1) /* SS */
    return 2;
  if(!flag){ /* SNS */
    if(tm1<0)
      return -1;
    else
      return 1;
  }
  if(!dt) /* singular but not SS */
    return 3;

  return 0; /* neither SNS nor singular */

}


int minorisSNSsing2(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;
  if(detsubmat(imat, n, m, vec1, vec2, k)==0)
    return 1;
  if(qualdetsubmat(imat, n, m, vec1, vec2, k)!=2)
    return 1;
  return 0;

}

// return the determinant

int minornotsing(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 0;
  return detsubmat(imat, n, m, vec1, vec2, k);

}

// return 1 for not SNS

int minornotSNS(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */
  if(qualdetsubmat(imat, n, m, vec1, vec2, k)==2)
    return 1;
  return 0;

}


/* compatible minors */
/* take two matrices imat1, and imat2 of equal dimensions */
/* assume that imat1 is a matrix, while imat2 is a sign pattern */
/* 1) check if either of imat1[vec1|vec2] or imat2[vec1|vec2] */
/*    is SS. If so, then we are done. */
/* 2) check if imat1[vec1|vec2]=0. If so done. Else: */
/* 3) check if imat2(vec1|vec2) is SNS. If not, done. Else */
/* 4) check imat1[vec1|vec2]imat2[vec1|vec2] > 0 */

int submat_signpat_compat(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k){
  int dt=0,qd=0;
 
  /* worth checking for rows and columns of zeros first */
  if(minorhas0rc(imat2, n, m, vec1, vec2, k) || minorhas0rc(imat1, n, m, vec1, vec2, k))
    return 1;

  if((dt=detsubmat(imat1, n, m, vec1, vec2, k))==0)
    return 1;
  if((qd=qualdetsubmat1(imat2, n, m, vec1, vec2, k))==0)
    return 1;
  if(qd==2)
    return 0;
  if(dt*qd<0) // opp signs: keep track of this info
    return -1;
  if(dt*qd>0) // a positive product: keep track of this info
    return 2;

  fprintf(stderr, "Error in \"submat_signpat_compat\": reached the end without returning. EXITING.");
  exit(0);
}



/* imat1 is an n X m matrix, while imat2 is an n X m sign-pattern */
/* check if they are compatible. */
/* flag = 3 means the matrices are compatible and r-strongly compatible */
/* flag = 2 means the matrices are compatible but not r-strongly compatible */
/* flag = 1 means the matrices are r-strongly compatible but not compatible */
/* flag = 0 means the matrices are none of the above */

int mat_signpat_compat0(int **imat1, int imatrank, int **imat2, int n, int m, int debug){
  int k;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=3, flg;
  //  int imatrank=matrank(imat1,n,m);
  bool posprod=0;
  for(k=imatrank;k>=1;k--){
    if(k==imatrank-1 && !posprod) // not r-strongly compatible
      flag=2;
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    if(n+m>15){
      fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    }
    for(r1=0;r1<cnk;r1++){
      if(r1%100==0 && (n+m>15))
	fprintf(stderr, ".");
      for(r2=0;r2<cmk;r2++){
	flg=submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k);
	if(flg<=0){ // fail to be compatible: return 1 or 0
	  if(debug){fprintf(stderr, "\nnumerical/pattern submatrices which fail to be compatible:\n\n");
	    printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	    fprintf(stderr, "*** and ***\n");
	    printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	  }
	  if(k==imatrank){flag=0;}else if(flag>=2){flag=flag-2;}
	  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	  return flag;
	}
	else if(flg==2 && k==imatrank) // a strictly positive product
	  posprod=1;

      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  return flag;
}

// are a matrix and sign pattern negative r-compatible

int mat_signpat_negr_compat(int **imat1, int imatrank, int **imat2, int n, int m, int debug){
  int k=imatrank;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1, flg;
  int posprod=0;
  if(debug){fprintf(stderr, "\nChecking negative compatibility...\n");}
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n, k);cmk=comb(m, k);
  if(debug){fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);}
  for(r1=0;r1<cnk;r1++){
    if(r1%100==0)
      fprintf(stderr, ".");
    for(r2=0;r2<cmk;r2++){
      flg=submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k);
      if(flg!=-1 && flg!=1){// fail to be negatively compatible
	if(debug){
	  fprintf(stderr, "\nnumerical/pattern submatrices which fail to be negatively compatible:\n\n");
	  printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	  fprintf(stderr, "*** and ***\n");
	  printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	}
	flag=0;
	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	return flag;
      }
      else if(flg==-1)
	posprod=1;
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  if(flag && posprod)
    return 1;
  return 0;
}

// wrapper for mat_signpat_compat0

int mat_signpat_compat(int **imat1, int imatrank, int **imat2, int n, int m, int debug){
  int flg;
  if(n+m>15)
    fprintf(stderr, "Depending on the speed of your computer, this could take a while.\n");

  flg=mat_signpat_compat0(imat1, imatrank, imat2, n, m, debug);
  if(flg)
    return flg;
  else if(mat_signpat_negr_compat(imat1, imatrank, imat2, n, m, debug))
    return -1;
  return 0;
}


// are two matrices compatible - only check largest minors

int mat_signpat_compat1(int **imat1, int **imat2, int n, int m, int rapid){
  int k;
  int r;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;

  r=min(n,m);
  k=r;
  
  xcombs=allcombsgen(n,k);
  fprintf(stderr, "here\n");
  ycombs=allcombsgen(m,k);
  fprintf(stderr, "here1\n");
  cnk=comb(n, k);
  cmk=comb(m, k);
  fprintf(stderr, "cnk= %.0f, cmk = %.0f\n", (double) cnk, (double)cmk);
  for(r1=0;r1<cnk;r1++){
    fprintf(stderr, "r1/cnk = %.0f/%.0f\n", (double)r1, (double)cnk);
    for(r2=0;r2<cmk;r2++){
      if(submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k)<=0){
	fprintf(stderr, "\npair of submatrices which fail to be compatible:\n\n");
	printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	fprintf(stderr, "*** and ***\n");
	printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	flag=0; // exit here for speed
	if(rapid){
	  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	  return flag;
	}
      }

    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  return flag;
}




/* Check if n X m matrix imat1 and n X m sign pattern imat2 */
/* are compatible. Remove empty rows/columns first. This is */
/* done twice because in theory each removal can create a new */
/* row/column of zeros in the other matrix */
// imat2 **must** be a pattern matrix (no entries other than 
//-1,0,1,2 where 2 means unsigned)
int arecompat(int **imat1, int **imat2, int n, int m, int debug){
  int **imat1a, **imat2a, **imat1b, **imat2b;
  int n1, m1, n2, m2;
  int flag=0;
  int imatrank=matrank(imat1, n, m);
  // simplify first
  if(debug){fprintf(stderr, "Checking compatibility...\n");
    fprintf(stderr, "Removing rows and columns of zeros...\n");
  }
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1);
  simppair(imat1a, imat2a, n1, m1, &imat1b, &imat2b, &n2, &m2);
  if(debug){printmat(imat1b, n2, m2);
    printmat(imat2b, n2, m2);
    fprintf(stderr, "Checking if matrix and sign-pattern are compatible, %d-strongly compatible or %d-strongly negatively compatible...\n", imatrank, imatrank);
  }
  flag=mat_signpat_compat(imat1b, imatrank, imat2b, n2, m2, debug);

  if(debug){
    if(flag==3)
      fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are compatible and %d-strongly compatible.\n---------------------------------------\n",imatrank);
    else if(flag==2)
      fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are compatible but not %d-strongly compatible.\n---------------------------------------\n", imatrank);
    else if(flag==1)
      fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are %d-strongly compatible but not compatible.\n---------------------------------------\n", imatrank);
    else if(flag==-1)
      fprintf(stderr, "Finished checking negative compatibility: matrix and sign-pattern are %d-strongly negatively compatible.\n---------------------------------------\n", imatrank);
    else
      fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are not compatible or %d-strongly compatible or %d-strongly negatively compatible.\n---------------------------------------\n", imatrank, imatrank);
  }
  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
  free_imatrix(imat1b, 0, n2-1, 0, m2-1);
  free_imatrix(imat2b, 0, n2-1, 0, m2-1);
  return flag;
}


// Simple version of arecompat
// imat2 **must** be a pattern matrix (no entries other than 
// -1,0,1,2 where 2 means unsigned)
// the flag strict checks for strong compatibility of at least one pair of
// minors of dimension rank(imat1)
// expect means sign of products to expect
int arecompatsimp(int **imat1, int **imat2, bool maxonly, int n, int m, int *strict, int debug){
  int **imat1a, **imat2a, **imat1b, **imat2b;
  int n1, m1, n2, m2;
  int imatrank=matrank(imat1, n, m);
  int k;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flg;
  int kmin=maxonly?imatrank:1;
  int expect=0;
  (*strict)=0;
  // simplify first
  if(debug){
    fprintf(stderr, "\n###Entering arecompatsimp...\n");
    fprintf(stderr, "Removing rows and columns of zeros...\n");
  }
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1);
  simppair(imat1a, imat2a, n1, m1, &imat1b, &imat2b, &n2, &m2);
  if(debug){
    printmat(imat1b, n2, m2);
    printmat(imat2b, n2, m2);
    fprintf(stderr, "Checking if matrix and sign-pattern are compatible...\n");
  }
  if(!m2 || !n2)//nothing left
    return 1;


  for(k=imatrank;k>=kmin;k--){
    xcombs=allcombsgen(n2,k);
    ycombs=allcombsgen(m2,k);
    cnk=comb(n2, k);cmk=comb(m2, k);
    if(debug && n2+m2>15){
      fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    }
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	flg=submat_signpat_compat(imat1b, imat2b, n2, m2, xcombs[r1], ycombs[r2], k);

	if(!expect && flg==2)//expect positive
	  expect=1;
	else if(!expect && flg==-1)//expect negative
	  expect=-1;
	if(flg==0 || (expect==1 && flg==-1) || (expect==-1 && flg==2)){ // fail to be compatible
	  if(debug){fprintf(stderr, "\nnumerical/pattern submatrices which fail to be compatible:\n\n");
	    printsubmat(imat1b, xcombs[r1], ycombs[r2], k, k);
	    fprintf(stderr, "*** and ***\n");
	    printsubmat(imat2b, xcombs[r1], ycombs[r2], k, k);
	  }
	  //exit(0);
	  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
	  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
	  free_imatrix(imat1b, 0, n2-1, 0, m2-1);
	  free_imatrix(imat2b, 0, n2-1, 0, m2-1);

	  return 0;
	}
	else if(flg==2||flg==-1)
	  (*strict)=1;
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
  free_imatrix(imat1b, 0, n2-1, 0, m2-1);
  free_imatrix(imat2b, 0, n2-1, 0, m2-1);
  if(debug){
    if(expect==1 && maxonly && strict)
      fprintf(stderr, "The matrix and sign pattern are %d-strongly compatible.\n\n",imatrank);
    else if(expect==1 && !maxonly)
      fprintf(stderr, "The matrix and sign pattern are compatible.\n\n");
    else if(expect==-1 && maxonly && strict)
      fprintf(stderr, "The matrix and sign pattern are %d-strongly negatively compatible.\n\n",imatrank);
    else if(expect==-1 && !maxonly)
      fprintf(stderr, "The matrix and sign pattern are negatively compatible.\n\n");
  }

  return expect;
}

//no third argument means all minors, not only maximal ones
int arecompatsimp(int **imat1, int **imat2, int n, int m, int *strict, int debug){
  return arecompatsimp(imat1, imat2, 0, n, m, strict, debug);
}

// imat1 is an nXm integer matrix
// imat2 **must** be an nXm pattern matrix (no entries other than 
// -1,0,1,2 where 2 means unsigned)
// Compatibility of A^(r) and B^(r) checked
// return codes of ispospoly: -2 for negative, 0 for zero, 2 for positive, 
// and -3 for unsigned (these are the only possibilities)
int arecompatr(int **imat1, int **imat2, int r, int n, int m, int debug){
  int **imat1a, **imat2a, **imat1b, **imat2b;
  int n1, m1, n2, m2;
  int k;
  long r1, r2, cnk, cmk;
  int **xcombs, **ycombs;
  int flg, expect=0;

  // simplify first
  if(debug){
    fprintf(stderr, "\n###Entering arecompatr...\n");
    fprintf(stderr, "Removing rows and columns of zeros...\n");
  }
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1);
  simppair(imat1a, imat2a, n1, m1, &imat1b, &imat2b, &n2, &m2);
  if(debug){
    printmat(imat1b, n2, m2);
    printmat(imat2b, n2, m2);
    fprintf(stderr, "Checking if matrix and sign-pattern are compatible...\n");
  }
  if(!m2 || !n2){//nothing left
    if(debug){fprintf(stderr, "All products are zero.\n");}
    return 0;
  }

  k=r;
  xcombs=allcombsgen(n2,k);
  ycombs=allcombsgen(m2,k);
  cnk=comb(n2, k);cmk=comb(m2, k);
  if(debug && n2+m2>15){
    fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
  }
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      flg=submat_signpat_compat(imat1b, imat2b, n2, m2, xcombs[r1], ycombs[r2], k);

      if(!expect && flg==2)//expect positive
	expect=1;
      else if(!expect && flg==-1)//expect negative
	expect=-1;
      if(flg==0 || (expect==1 && flg==-1) || (expect==-1 && flg==2)){ // fail to be compatible
	if(debug){fprintf(stderr, "Found either an unsigned determinant, or oppositely signed pairs\n");}
	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	free_imatrix(imat1a, 0, n1-1, 0, m1-1);
	free_imatrix(imat2a, 0, n1-1, 0, m1-1);
	free_imatrix(imat1b, 0, n2-1, 0, m2-1);
	free_imatrix(imat2b, 0, n2-1, 0, m2-1);
	return -3;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
  free_imatrix(imat1b, 0, n2-1, 0, m2-1);
  free_imatrix(imat2b, 0, n2-1, 0, m2-1);
  if(expect==1){
    if(debug){fprintf(stderr, "The matrix and sign pattern are %d-strongly compatible.\n\n",k);}
    return 2;
  }
  else if(expect==-1){
    if(debug){fprintf(stderr, "The matrix and sign pattern are %d-strongly negatively compatible.\n\n",k);}
    return -2;
  }

  if(debug){fprintf(stderr, "All products are zero.\n");}
  return 0;
}



/* Check if n X m matrix mat is SSD. The routine */
/* is a wrapper for isSSD1 adding a step removing */
/* columns/rows of zeros. */

int isSSD(int **mat, int n, int m, int q){
  int flag=0;
  int **imat, n1, m1;
  int imatrank=matrank(mat,n,m);
  int rkbad=0; // rank of bad minors found
  fprintf(stderr, "Checking if the matrix is SSD or %d-SSD...\n",imatrank);
  if(n+m > 15)
    fprintf(stderr, "...please be patient. This could take time.\n");
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isSSD1(imat, n1, m1, 1, &rkbad, q);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  if(flag){
    fprintf(stderr, "Finished checking SSD: the matrix is SSD.\n---------------------------------------\n");
    return 2;
  }
  else if(rkbad==imatrank){
    fprintf(stderr, "Finished checking SSD: the matrix is not SSD or %d-SSD.\n---------------------------------------\n",imatrank);
    return 0;
  }
  flag=isSSD1(mat, n, m, 0, &rkbad, q);

  if(flag){
    fprintf(stderr, "Finished checking SSD: the matrix is %d-SSD but not SSD.\n---------------------------------------\n",imatrank);
    return 1;
  }
  else{
    fprintf(stderr, "Finished checking SSD: the matrix is not SSD or %d-SSD.\n---------------------------------------\n",imatrank);
    return 0;
  }
}



/* Check if n X m matrix mat is CSD (All square */
/* submatrices are sign nonsingular or sign singular) */

/* flag = 2 means CSD */
/* flag = 1 means r-CSD but not CSD */
/* flag = 0 means neither */


int isCSD1(int **imat, int n, int m, bool allm, int *rkbad, int q){
  int k,r;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  unsigned long fk;
  int **pms;
  int ret;
  int rtrig=4;
  int rmin=2;
  int imatrank=matrank(imat,n,m);

  if(allm){
    r=min(n,m);rmin=2;
  }
  else{
    r=imatrank;rmin=r;
  }

  for(k=r;k>=rmin;k--){
    //fprintf(stderr, "k=%d\n", k);
    xcombs=allcombsgen(n,k);ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    fk=factorial(k);
    if(k<=rtrig){// smallish matrix
      for(r2=0;r2<cmk;r2++){
	pms=allperms1(ycombs[r2], k);
	for(r1=0;r1<cnk;r1++){
	  if(k>imatrank)// only check for sign singularity
	    ret=minorisSS(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);
	  else
	    ret=minorisSNSSS(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);

	  if(!ret){
	    fprintf(stderr, "submatrix which fails to be sign nonsingular or sign singular:\n");
	    printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	    flag=0;if(!(*rkbad)){(*rkbad)=k;}// rank of bad submatrix
	    if(q){
	      free_imatrix(pms,0, fk-1, 0, k);
	      free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	      free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	      return flag;
	    }
	  }
	}
      }
    }
    else{// large matrix
      for(r2=0;r2<cmk;r2++){
	for(r1=0;r1<cnk;r1++){
	  if(k>imatrank)// only check for sign singularity
	    ret=minorisSS2(imat, n, m, xcombs[r1], ycombs[r2], k, fk);
	  else
	    ret=minorisSNSSS_rec(imat, n, m, xcombs[r1], ycombs[r2], k);//recursive version
	  //	    ret=minorisSNSSS2(imat, n, m, xcombs[r1], ycombs[r2], k,fk);

	  if(!ret){
	    fprintf(stderr, "submatrix which fails to be sign nonsingular or sign singular:\n");
	    printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	    flag=0;if(!(*rkbad)){(*rkbad)=k;}// rank of bad submatrix
	    if(q){
	      free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	      free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	      return flag;
	    }
	  }

	}
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  return flag;
}

/* Check if n X m matrix imat is SSD */

int isSSD1(int **imat, int n, int m, bool allm, int *rkbad, int q){
  int k, ret=0;
  int r=100;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  unsigned long fk,t=0;
  int **pms;
  int rmin=2,szswtch=0;
  int imatrank=matrank(imat,n,m);


  if(allm){
    r=imatrank;rmin=2;
  }
  else{
    r=imatrank;rmin=r;
  }

  // all minors of size greater than imatrank are definitely singular, so nothing to check
  
  for(k=r;k>=rmin;k--){
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    fk=factorial(k);

    t=0;
    for(r2=0;r2<cmk;r2++){
      /* if(t%100==0) */
      /* 	fprintf(stderr, "t=%.0f\n", (double) t); */
      /* t++; */

      if(k<=szswtch)
	pms=allperms1(ycombs[r2], k);
      for(r1=0;r1<cnk;r1++){
	if(k>szswtch) // large submatrix
	  ret=minorisSNSsing2(imat, n, m, xcombs[r1], ycombs[r2], k);
	else
	  ret=minorisSNSsing(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);

	if(!ret){
	  fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular (t=%.0f):\n", (double) t);
	  //fprintf(stderr, "det = %d\n", detsubmat(imat, n, m, xcombs[r1], ycombs[r2], k));
	  //fprintf(stderr, "qualdet = %d\n", qualdetsubmat(imat, n, m, xcombs[r1], ycombs[r2], k-1));
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0;if(!(*rkbad)){(*rkbad)=k;}// rank of bad submatrix
	  if(q){ // exit here for speed
	    if(k<=szswtch)
	      free_imatrix(pms,0, fk-1, 0, k);
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    return flag;
	  }
	}
      }
      if(k<=szswtch)
	free_imatrix(pms,0, fk-1, 0, k);
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  return flag; // SSD (and hence automatically r-SSD)
}

/* Check all minors of size minorsize of an n X m matrix imat */
/* to see if they are SNS or singular */
/* get rid of pms and fk - slow and unnecessary */

int checkminors(char *fname, int **imat, int **vmat, int n, int m, int minorsize){
  int k, ret=0, cutsz;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1, flag1;
  unsigned long fk;
  int **pms;
  int rtrig=10;
  FILE *fd, *fd1;
  long ct=1, ct1=1, ct2=1;
  int xc[minorsize], yc[minorsize];
  int failflag=0;
  long lab=1;
  int *cut1,*cut2;
  int totgraphs=0, totgold=0;
  int **allgraphs=NULL;

  fd = fopen(fname, "w");
  if(!fd){
    fprintf(stderr, "ERROR in checkminors: \"%s\" could not be opened for reading.\n", fname);
    return 0;
  }
  fd1=fopen("graphs.dot", "w");
  if(!fd1){
    fprintf(stderr, "ERROR in checkminors: \"graphs.dot\" could not be opened for reading.\n");
    return 0;
  }

  k=minorsize;
  fprintf(fd, "ratmx:true;\n");
  fprintf(stderr, "Dimensions: matrices are %d X %d\n\n", n,m);
  fprintf(stderr, "Expected rank = %d. Checking only minors of size %d\n\n", k, k);

  if(n>20 || m > 20){
    fprintf(stderr, "dimensions %d X %d are large. Routine may not terminate.\n ", n, m);

    firstcomb(xc, n, k);
    firstcomb(yc, m, k);
    flag=1;
    flag1=1;


    //    printvec(xc,k);
    //    printvec(yc,k);

    while(flag==1){
      flag1=1;
      while(flag1==1){
	ct1++;
	if(ct1%10000==0){
	  fprintf(stderr, "%ld e04\n", ct2);ct2++;
	  //	  printvec(xc,k);
	  //	  printvec(yc,k);
	}
	ret=minorisSNSsing2(imat, n, m, xc, yc, k);

	if(!ret){
	  //	  fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular:\n");
	  //	  fprintf(stderr, "rows:    ");printvec(xc, k);fprintf(stderr, "columns: ");printvec(yc, k);
	  //	  printsubmat(imat, xc, yc, k, k);
	  //	  printmaximaindsubmat(fd, vmat, xc, yc, k, k);
	  //	  fprintf(stderr, "qualitative determinant = %d\n\n", qualdetsubmat(imat, n, m, xc, yc, k));

	  cutsz=cuttails(vmat, k, k, xc, yc, &cut1, &cut2);
	  totgold=totgraphs;
	  totgraphs=addtwointstr(totgraphs, cut1, cut2, k-cutsz, k-cutsz, &allgraphs);
	  if(totgraphs>totgold){
	    printdotindsubmat(fd1, vmat, cut1, cut2, k-cutsz, k-cutsz, lab);lab++;
	  }
	  free((char *) cut1);free((char *) cut2);

	//	  printdotindsubmat(fd1, vmat, xc, yc, k, k, lab);lab++;
	  failflag=1;
	  /* if(q){ */
	  /*   free((char *) xc);free((char *) yc); */
	  /*   return 0; */
	  /* } */
	}
      
	flag1=nextcomb(yc, m, k);
      }
      flag=nextcomb(xc, n, k);
    }
    
    if(!failflag)
      fprintf(stderr, "The second additive compound of the Jacobian is well behaved.\n");
    else
      fprintf(stderr, "The second additive compound of the Jacobian fails to be well behaved.\n");

    fclose(fd);fclose(fd1);
    return 0;
  }

  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n, k);cmk=comb(m, k);
  fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
  fk=factorial(k);

  for(r2=0;r2<cmk;r2++){
    if(k<rtrig)
      pms=allperms1(ycombs[r2], k);
    for(r1=0;r1<cnk;r1++){
      ct1++;
      if(ct1%10000==0){
	fprintf(stderr, "%ld e04\n", ct2);ct2++;
	//	  printvec(xc,k);
	//	  printvec(yc,k);
      }

      if(k>=rtrig) // large submatrix
	ret=minorisSNSsing2(imat, n, m, xcombs[r1], ycombs[r2], k);
      else
	ret=minorisSNSsing(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);
      
      if(!ret){
	fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular:\n");
	for(ct=0;ct<k;ct++)
	  fprintf(stderr, "%d,",ycombs[r2][ct]+1);
	fprintf(stderr, "\n");

	printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	printmaximaindsubmat(fd, vmat, xcombs[r1], ycombs[r2], k, k);

	cutsz=cuttails(vmat, k, k, xcombs[r1], ycombs[r2], &cut1, &cut2);
	totgold=totgraphs;
	totgraphs=addtwointstr(totgraphs, cut1, cut2, k-cutsz, k-cutsz, &allgraphs);
	if(totgraphs>totgold){
	  printdotindsubmat(fd1, vmat, cut1, cut2, k-cutsz, k-cutsz, lab);lab++;
	}
	//	printdotindsubmat(fd1, vmat, xcombs[r1], ycombs[r2], n1, m1, lab);lab++;
       	free((char *) cut1);free((char *) cut2);
	flag=0; // exit here for speed
      }

    }
    if(k<rtrig)
      free_imatrix(pms,0, fk-1, 0, k);
  }

  for(k=0;k<totgraphs;k++)
    free ((char *)(allgraphs[k]));
  free((char *) allgraphs);

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  fclose(fd);fclose(fd1);
  fprintf(stderr, "totgraphs = %d\n", totgraphs);
  return flag;
}



// get the position of a vector [i1, ...,ir] in a list
// of lexicographically ordered r-vectors chosen from {1,...,n}
// starts at position 0


long getpos(int *vec, int n, int r){
  long tot=0;
  int i,j, k1=n-1, k2=r-1, k3=0;
  for(j=0;j<r;j++){
    for(i=0;i<vec[j]-k3;i++){
      tot+=comb(k1, k2);
      k1--;
    }
    k1--;k2--;k3=vec[j]+1;
  }
  return tot;
}



int S2(int **S, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1){
  int **xcombs;
  int **tmp, **tmp1;
  long tot, r1, r2;
  int j;
  xcombs=allcombsgen(n,2);
  tot=comb(n,2);
  tmp=imatrix(0, tot-1, 0, 2*m*tot-1);
  tmp1=imatrix(0, tot-1, 0, 2*m*tot-1);

  for(r1=0;r1<tot;r1++){ //each row
    for(r2=0;r2<tot;r2++){ // each block
      //block one
      if(xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=S[xcombs[r1][0]][j];
	}
      }
      else if(xcombs[r1][0]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=-S[xcombs[r1][1]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=0;
	}
      }

      // block 2
      if(xcombs[r1][0]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=S[xcombs[r1][1]][j];
	}
      }
      else if(xcombs[r1][1]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=-S[xcombs[r1][0]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=0;
	}
      }


      //imat2

      if(xcombs[r1][0]==xcombs[r2][0] && xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp1[r1][2*m*r2+j]=sgn(S[xcombs[r1][0]][j]);
	  tmp1[r1][2*m*r2+m+j]=sgn(S[xcombs[r1][1]][j]);
	}
      }

    }
  }

  free_imatrix(xcombs, 0, comb(n, 2)-1, 0, 1);


  simppair(tmp, tmp1, comb(n,2), 2*m*comb(n,2), imat1, imat2, n1, m1); 
    free_imatrix(tmp, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);
    free_imatrix(tmp1, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);


  return 1;
}

// like S2, but with V explicitly

int S2a(int **S, int **V, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1){
  int **xcombs;
  int **tmp, **tmp1;
  long tot, r1, r2;
  int j;
  xcombs=allcombsgen(n,2);
  tot=comb(n,2);
  tmp=imatrix(0, tot-1, 0, 2*m*tot-1);
  tmp1=imatrix(0, tot-1, 0, 2*m*tot-1);

  for(r1=0;r1<tot;r1++){ //each row
    for(r2=0;r2<tot;r2++){ // each block
      //block one
      if(xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=S[xcombs[r1][0]][j];
	}
      }
      else if(xcombs[r1][0]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=-S[xcombs[r1][1]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=0;
	}
      }

      // block 2
      if(xcombs[r1][0]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=S[xcombs[r1][1]][j];
	}
      }
      else if(xcombs[r1][1]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=-S[xcombs[r1][0]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=0;
	}
      }


      //imat2

      if(xcombs[r1][0]==xcombs[r2][0] && xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp1[r1][2*m*r2+j]=sgn(V[xcombs[r1][0]][j]);
	  tmp1[r1][2*m*r2+m+j]=sgn(V[xcombs[r1][1]][j]);
	}
      }

    }
  }

  free_imatrix(xcombs, 0, comb(n, 2)-1, 0, 1);


  simppair(tmp, tmp1, comb(n,2), 2*m*comb(n,2), imat1, imat2, n1, m1); 
  free_imatrix(tmp, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);
  free_imatrix(tmp1, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);


  return 1;
}

int S2b(int **S, int **V, int n, int m, int ***imat1, int **valsvec, int **indsvec, int *base, int *basetot, int *n1, int *m1){
  int **xcombs;
  long tot, r1, r2;
  int j, totlen;
  int *tmpvals, *tmpinds;
  int **tmpmat;
  xcombs=allcombsgen(n,2);
  tot=comb(n,2);

  tmpmat=imatrix(0, tot-1, 0, 2*m*tot-1);
  tmpvals=(int *)malloc(sizeof(int*) *(2*m*tot));
  tmpinds=(int *)malloc(sizeof(int*) *(2*m*tot));

  for(r1=0;r1<tot;r1++){ //each row
    for(r2=0;r2<tot;r2++){ // each block
      //block one
      if(xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+j]=S[xcombs[r1][0]][j];
	}
      }
      else if(xcombs[r1][0]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+j]=-S[xcombs[r1][1]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+j]=0;
	}
      }

      // block 2
      if(xcombs[r1][0]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+m+j]=S[xcombs[r1][1]][j];
	}
      }
      else if(xcombs[r1][1]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+m+j]=-S[xcombs[r1][0]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+m+j]=0;
	}
      }


      //valsvec and indsvec // change - only do once

      if(xcombs[r1][0]==xcombs[r2][0] && xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmpvals[2*m*r2+j]=sgn(V[xcombs[r1][0]][j]);
	  tmpvals[2*m*r2+m+j]=sgn(V[xcombs[r1][1]][j]);

	  tmpinds[2*m*r2+j]=xcombs[r1][0]*m+j;
	  tmpinds[2*m*r2+m+j]=xcombs[r1][1]*m+j;
	}
      }

    }
  }

  totlen=0;
  // get rid of zeros
  for(r2=0;r2<tot;r2++){
    for(j=0;j<2*m;j++){
      if(tmpvals[2*m*r2+j]!=0)
	totlen++;
    }
  }



  (*valsvec)=(int *)malloc(sizeof(int*) *(totlen));
  (*indsvec)=(int *)malloc(sizeof(int*) *(totlen));
  (*imat1)=imatrix(0, tot-1, 0, totlen-1);
  (*n1)=tot;(*m1)=totlen;

  totlen=0;
  for(r2=0;r2<tot;r2++){
    base[r2]=0;
    for(j=0;j<2*m;j++){
      if(tmpvals[2*m*r2+j]!=0){
	(*valsvec)[totlen]=tmpvals[2*m*r2+j];
	(*indsvec)[totlen++]=tmpinds[2*m*r2+j];
	(base[r2])++;
      }
    }
  }
  
  basetot[0]=0; // running total
  for(r2=1;r2<tot;r2++){
    basetot[r2]=basetot[r2-1]+base[r2-1];
  }

  for(r1=0;r1<tot;r1++){ //each row
    totlen=0;
    for(r2=0;r2<tot;r2++){ // each block
      for(j=0;j<2*m;j++){
	if(tmpvals[2*m*r2+j]!=0){
	  (*imat1)[r1][totlen++]=tmpmat[r1][2*m*r2+j];
	}

      }
    }
  }

  printvec(tmpvals, 2*m*tot);
  printvec((*valsvec), (*m1));
  printvec(base, tot);
  printvec(basetot, tot);


  free((char *) tmpvals);
  free((char *) tmpinds);
  free_imatrix(tmpmat, 0, tot-1, 0, 2*m*tot-1);
  free_imatrix(xcombs, 0, comb(n, 2)-1, 0, 1);

  return 1;
}





//detsk has dimension (0, comb(n, k-1), 0, comb(m, k-1))

int **detsk(int **imat, int n, int m, int k){
  int **xcombs;
  int **ycombs;
  int **tmp;
  unsigned long r1, r2;
  int **submat;

  tmp=imatrix(0, comb(n, k-1), 0, comb(m, k-1));
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  for(r1=0;r1<comb(n, k);r1++){
    for(r2=0;r2<comb(m, k);r2++){
      submat = submatgen(imat, n, m, xcombs[r1], ycombs[r2],k);
      tmp[r1][r2]=det(submat,k);
      free_imatrix(submat, 0, k-1, 0, k-1);
    }
  }

  free_imatrix(xcombs, 0, comb(n, k)-1, 0, k-1);
  free_imatrix(ycombs, 0, comb(m, k)-1, 0, k-1);
  return tmp;

}

/* Check if n X m matrix mat is CSD. The routine */
/* is a wrapper for isCSD1 adding a step removing */
/* columns/rows of zeros. */

int isCSD(int **mat, int n, int m, int q){
  int flag;
  int **imat, n1, m1;
  int imatrank=matrank(mat,n,m);
  int rkbad=0; // rank of bad minors found
  fprintf(stderr, "Checking if the matrix is CSD or %d-CSD...\n", imatrank);
  printmat(mat, n, m);
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isCSD1(imat, n1, m1, 1, &rkbad, q);
  free_imatrix(imat,0, n1-1, 0, m1-1);

  if(flag){
    fprintf(stderr, "Finished checking CSD: the matrix is CSD.\n---------------------------------------\n");
    return 2;
  }
  else if(rkbad==imatrank){
    fprintf(stderr, "Finished checking CSD: the matrix is not CSD or %d-CSD.\n---------------------------------------\n", imatrank);
    return 0;
  }

  flag=isCSD1(mat, n, m, 0, &rkbad, q);

  if(flag){
    fprintf(stderr, "Finished checking CSD: the matrix is %d-CSD but not CSD.\n---------------------------------------\n", imatrank);
    return 1;
  }
  else{
    fprintf(stderr, "Finished checking CSD: the matrix is not CSD or %d-CSD.\n---------------------------------------\n", imatrank);
    return 0;
  }
}


int doubleisWSD2(int **imat, int n, int m){
  int **tmp;int flag=0;
  tmp=doublemat(imat, n, m);
  printmat(tmp, n, 2*m);
  if (isWSD2(tmp, n, 2*m))
    flag=1;
  free_imatrix(tmp, 0, n-1, 0, 2*m-1);
  return flag;

}

// WSD1 is faster than WSD2

int doubleisWSD(int **mat, int n, int m, int q){
  int **imat, **tmp;
  int n1,m1,flag1=0,flag2=0;
  int rkbad=0;
  int imatrank=matrank(mat,n,m);
  imat=simpmat(mat, n, m, &n1, &m1); // first simplify and then double
  tmp=doublemat(imat, n1, m1);
  //fprintf(stderr, "Checking if this matrix is WSD...\n");
  //printmat(tmp, n1, 2*m1);
  flag1=isWSD1(tmp, n1, 2*m1, 1, &rkbad, q); //all
  free_imatrix(tmp, 0, n1-1, 0, 2*m1-1);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  if(flag1==0 && rkbad==imatrank){// no need to check rWSD
    fprintf(stderr, "Finished checking WSD: the matrix is neither r-strongly WSD nor WSD.\n---------------------------------------\n");
    return 0;
  }
  // for r-strong WSD, don't simplify
  tmp=doublemat(mat, n, m);
  fprintf(stderr, "Checking if this matrix is WSD or r-strongly WSD...\n");
  printmat(tmp, n, 2*m);
  flag2=isWSD1(tmp, n, 2*m, 0, &rkbad, q); // only maximal
  free_imatrix(tmp, 0, n-1, 0, 2*m-1);

  if(flag1 && flag2){
    fprintf(stderr, "Finished checking WSD: the matrix is WSD and r-strongly WSD.\n---------------------------------------\n");
    return 3;
  }
  else if (flag1){
    fprintf(stderr, "Finished checking WSD: the matrix is WSD but not r-strongly WSD.\n---------------------------------------\n");
    return 2;
  }
  else if(flag2){
    fprintf(stderr, "Finished checking WSD: the matrix is r-strongly WSD but not WSD.\n---------------------------------------\n");
    return 1;
  }
  else{
    fprintf(stderr, "Finished checking WSD: the matrix is neither r-strongly WSD nor WSD.\n---------------------------------------\n");
    return 0;
  }


}



/* A wrapper for isWSD1, first removing rows/columns of zeros */

int isWSD(int **mat, int n, int m, int q){
  int **imat;
  int n1,m1,flag=0;
  int rkbad=0;
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isWSD1(imat, n1, 2*m1, 1, &rkbad, q);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  return flag;

}


int isWSD2(int **imat, int n, int m){
  int k,r;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int **submat;
  int flag=1;
  r=min(n,m);
  for(k=2;k<=r;k++){
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n,k);cmk=comb(m,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	submat = submatgen(imat, n, m, xcombs[r1], ycombs[r2],k);
	if(!isWS(submat, k)){
	  fprintf(stderr, "submatrix which fails to det(S)det(S-)>=0:\n");
	  printmat(submat, k,k);
	  flag=0; // exit here for speed
	}
	free_imatrix(submat, 0, k-1, 0, k-1);
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }
  return flag;
}

/* Direct computation on n X m matrix imat to see */
/* if it is WSD. The first step (WSpair) is to find the */
/* matrix of exponents, and to remove rows/columns of zeros. */
/* Then compatibility is checked. */

/* flag = 3 means the matrices are compatible and r-strongly compatible */
/* flag = 2 means the matrices are compatible but not r-strongly compatible */
/* flag = 1 means the matrices are r-strongly compatible but not compatible */
/* flag = 0 means the matrices are none of the above */

int isWSD1(int **imat, int n, int m, bool allm, int *rkbad, int quick){
  int k,rmin,rmax,tt;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  int **imat2, **imat3, n1, m1;
  int mr1, mr2;
  bool posprod=0;

  WSpair(imat, n, m, &imat2, &imat3, &n1, &m1); // imat3 is sparser

  if(n1!=n || m1!=m){
    fprintf(stderr, "The simplified matrix pair...\n");
    printmat(imat2, n1, m1);
    printmat(imat3, n1, m1);
  }

  mr1=matrank(imat2,n1,m1);mr2=matrank(imat3,n1,m1);
  rmax=min(mr1,mr2);

  if(!allm && rmax<mr1){ // can't be r-strongly compatible
    free_imatrix(imat2, 0, n1-1, 0, m1-1);
    free_imatrix(imat3, 0, n1-1, 0, m1-1);
    return 0;
  }

  if(allm)
    rmin=2;
  else
    rmin=rmax;

  for(k=rmax;k>=rmin;k--){
    xcombs=allcombsgen(n1,k);ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	if((tt=fixedminorcompat(imat2, imat3, n1, m1, xcombs[r1], ycombs[r2],k))<0){
	  fprintf(stderr, "submatrix which fails det(S)det(S-) >=0:\n");
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0;if((*rkbad)==0){(*rkbad)=k;}
	  if(quick){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    free_imatrix(imat2, 0, n1-1, 0, m1-1);
	    free_imatrix(imat3, 0, n1-1, 0, m1-1);
	    return flag;
	  }
	}
	else if(tt==1)// a positive product
	  posprod=1;
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }
  free_imatrix(imat2, 0, n1-1, 0, m1-1);
  free_imatrix(imat3, 0, n1-1, 0, m1-1);

  if(!allm && !posprod)// not r-strongly WSD
    return 0;

  return flag;
}

//The brief version, simply checking for normality
//S is the stoich matrix, Sl is the exponent matrix
int normal(int **S, int **Sl, int n, int m, int rk){
  int flag, flag1;
  int xc[rk], yc[rk];

  firstcomb(xc,n,rk);
  flag=1;
  while(flag==1){
    firstcomb(yc,m,rk);
    flag1=1;
    while(flag1==1){
      if(detsubmat(S,n,m,xc,yc,rk)!=0 && detsubmat(Sl,n,m,xc,yc,rk)!=0)
	return 1;
      flag1=nextcomb(yc,m,rk);
    }
    flag=nextcomb(xc,n,rk);
  }
  return 0;
}

//overloading
int normal(int **AM, int n, int m){
  int **S, **Sl;
  int ret,Srank,minus=1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  Srank=matrank(S,n,m);
  ret=normal(S, Sl, n, m, Srank);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading di6 input
int normal(char *di6, int n, int m){
  int **S, **Sl;
  int ret,Srank,minus=1;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  Srank=matrank(S,n,m);
  ret=normal(S, Sl, n, m, Srank);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//Check if a numeric-symbolic submatrix-product is identically zero
int submat_signpat_zero(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k){ 
  if(minorhas0rc(imat2, n, m, vec1, vec2, k) || minorhas0rc(imat1, n, m, vec1, vec2, k) || detsubmat(imat1, n, m, vec1, vec2, k)==0 || qualdetsubmat1(imat2, n, m, vec1, vec2, k)==0)
    return 1;
  return 0;
}

//Normality for positive general kinetics
//Jacobian matrix not identically singular on stoich subspace
//S is the stoich matrix, V is the pattern matrix
int PGKnormal(int **S, int **V, int n, int m, int rk){
  int flag, flag1;
  int xc[rk], yc[rk];

  firstcomb(xc,n,rk);
  flag=1;
  while(flag==1){
    firstcomb(yc,m,rk);
    flag1=1;
    while(flag1==1){
      if(!submat_signpat_zero(S, V, n, m, xc, yc, rk))
	 return 1;
      flag1=nextcomb(yc,m,rk);
    }
    flag=nextcomb(xc,n,rk);
  }
  return 0;
}

//overloading
int PGKnormal(int **AM, int n, int m){
  int **S, **Sl;
  int ret,Srank,minus=1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  Srank=matrank(S,n,m);
  ret=PGKnormal(S, Sl, n, m, Srank);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading di6 input
int PGKnormal(char *di6, int n, int m){
  int **S, **Sl;
  int ret,Srank,minus=1;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  Srank=matrank(S,n,m);
  ret=PGKnormal(S, Sl, n, m, Srank);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}



/* Check if n X m matrices imat1 and imat2 are compatible */
/* flag = 3 means the matrices are compatible and r-strongly compatible */
/* flag = 2 means the matrices are compatible but not r-strongly compatible */
/* flag = 1 means the matrices are r-strongly compatible but not compatible */
/* flag = 0 means the matrices are none of the above */

int mats_compat0(int **imat1a, int imatrank, int **imat2a, int n1, int m1, int quick, int debug){
  int k;
  int r,rmin,s1,s2;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int cmpt=1; // compatible
  int rscmpt=1; // r-strongly compatible
  int tt;
  bool posprod=0;

  r=imatrank;
  rmin=1;

  for(k=r;k>=rmin;k--){
    xcombs=allcombsgen(n1,k);
    ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	//	fprintf(stderr, "in here\n");
	if((tt=fixedminorcompat(imat1a, imat2a, n1, m1, xcombs[r1], ycombs[r2],k)) < 0){
	  if(debug){
	    fprintf(stderr, "\npair of submatrices with determinants of opposite sign:\n\n");
	    printsubmat(imat1a, xcombs[r1], ycombs[r2], k, k);
	    fprintf(stderr, "*** and ***\n");
	    printsubmat(imat2a, xcombs[r1], ycombs[r2], k, k);

	    for(s1=0;s1<k;s1++){
	      fprintf(stderr, "[ ");
	      for(s2=0;s2<k;s2++)
		fprintf(stderr, "%d,%d  ",xcombs[r1][s1],ycombs[r2][s2]);
	      fprintf(stderr, "]\n");
	    }
	    fprintf(stderr, "_____________________\n\n");
	  }
	  if(k==imatrank){rscmpt=0;}else{cmpt=0;} // could return here for speed
	  if(quick){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);

	    if(rscmpt==0) // not r-strongly
	      return 0;
	    else
	      return 1;
	  }
	}
	else if(tt>0 && k==imatrank){
	  posprod=1;
	}
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
    if(k==imatrank && !posprod)
      rscmpt=0;
  }

  if(rscmpt && cmpt)
    return 3;
  else if(rscmpt && !cmpt)
    return 1;
  else if(!rscmpt && cmpt)
    return 2;
  return 0;

}


// are two matrices r-strongly negatively compatible?
int mats_negr_compat(int **imat1a, int imatrank, int **imat2a, int n1, int m1, int debug){
  int k, s1,s2;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int tt;
  bool nonzprod=0;

  k=imatrank;
  xcombs=allcombsgen(n1,k);
  ycombs=allcombsgen(m1,k);
  cnk=comb(n1,k);cmk=comb(m1,k);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      //	fprintf(stderr, "in here\n");
      if((tt=fixedminorcompat(imat1a, imat2a, n1, m1, xcombs[r1], ycombs[r2],k)) > 0){ // always fail immediately
	if(debug){
	  fprintf(stderr, "\npair of submatrices with determinants of same sign:\n\n");
	  printsubmat(imat1a, xcombs[r1], ycombs[r2], k, k);
	  fprintf(stderr, "*** and ***\n");
	  printsubmat(imat2a, xcombs[r1], ycombs[r2], k, k);
	  for(s1=0;s1<k;s1++){
	    fprintf(stderr, "[ ");
	    for(s2=0;s2<k;s2++)
	      fprintf(stderr, "%d,%d  ",xcombs[r1][s1],ycombs[r2][s2]);
	    fprintf(stderr, "]\n");
	  }

	  fprintf(stderr, "_____________________\n\n");
	}
	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	return 0;
      }
      else if(tt<0)
	nonzprod=1;
    }
  }

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  if(nonzprod)
    return 1;

  return 0;

}

// Wrapper for mats_compat0 and mats_negr_compat

int mats_compat(int **imat1, int **imat2, int n, int m, int quick, int debug){
  int flg, imatrank;
  int **imat1a, **imat2a;
  int n1,m1;

  imatrank=matrank(imat1,n,m);
  if(debug){fprintf(stderr, "Checking if (stoich and exp) matrices are compatible, %d-strongly compatible or %d-strongly negatively compatible.\n\n",imatrank,imatrank);}
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1); // imat2a is sparser
  if(debug){printmat(imat1a, n1, m1);printmat(imat2a, n1, m1);}

  flg=mats_compat0(imat1a, imatrank, imat2a, n1, m1, quick, debug);
  if(flg){
    free_imatrix(imat1a, 0, n1-1, 0, m1-1);
    free_imatrix(imat2a, 0, n1-1, 0, m1-1);
    if(flg==3){ // both succeeded
      if(debug){fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are compatible and %d-strongly compatible.\n---------------------------------------\n", imatrank);}
      return 3;
    }
    else if(flg==1){
      if(debug){fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are %d-strongly compatible, but not compatible.\n---------------------------------------\n", imatrank);}
      return 1;
    }
    else if(flg==2){
      if(debug){fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are compatible but not %d-strongly compatible.\n---------------------------------------\n", imatrank);}
      return 2;
    }
    return flg;
  }


  flg=mats_negr_compat(imat1a, imatrank, imat2a, n1, m1, debug);
  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
  if(flg){
    if(debug){fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are %d-strongly negatively compatible.\n---------------------------------------\n", imatrank);}
    return -1;
  }
  if(debug){fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are not compatible or %d-strongly compatible or %d-strongly negatively compatible.\n---------------------------------------\n", imatrank, imatrank);}

  return 0;
}


// Semiconcordance (the MA (reduced) Jacobian matrix everywhere has the same
// rank as the stoichiometric matrix)
// Sf=stoichiometric matrix; Slf=minus left stoichiometric matrix
int semiconcord(int **Sf, int **Slf, int n, int m, int debug){
  long r1,r2;
  int **S, **Sl;
  int n1,m1;
  int k=matrank(Sf,n,m);
  long cnk, cmk;
  int x1[k],y1[k];//storage
  int **xcombs, **ycombs;
  int tt;
  bool nonzprod=0;
  int prodsofar=0;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering semiconcord.\n");}
  if(debugfull){fprintf(stderr,"Matrices:\n");printmat(Sf,n,m);printmat(Slf,n,m);}

  simppair(Sf, Slf, n, m, &S, &Sl, &n1, &m1); 

  if(debugfull){fprintf(stderr, "Simplified pair:\n");printmat(S, n1, m1);printmat(Sl, n1, m1);}

  cnk=comb(n1,k);cmk=comb(m1,k);xcombs=allcombsgen(n1,k);ycombs=allcombsgen(m1,k);

  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      tt=fixedminorcompat(S, Sl, n1, m1, xcombs[r1], ycombs[r2],k);
      if(!prodsofar && tt){
	nonzprod=1;
	veccp(xcombs[r1],k,x1);veccp(ycombs[r2],k,y1);//for later debugging
	if(tt>0){prodsofar=1;}
	else if(tt<0){prodsofar=-1;}
      }
      else if((prodsofar>0 && tt<0) || (prodsofar<0 && tt>0)){//fail immediately
	if(debugfull){
	  fprintf(stderr, "Not semiconcordant. Earlier pair of %s\n", prodsofar>1?"same sign:":"opposite signs:");
	  printsubmat(S, x1, y1, k, k);
	  fprintf(stderr, "*** and ***\n");
	  printsubmat(Sl, x1, y1, k, k);
	  fprintf(stderr, "\nCurrent pair of %s\n", tt>0?"same sign:":"opposite signs:");
	  printsubmat(S, xcombs[r1], ycombs[r2], k, k);
	  fprintf(stderr, "*** and ***\n");
	  printsubmat(Sl, xcombs[r1], ycombs[r2], k, k);
	}
	else if(debug)
	  fprintf(stderr, "The network is not semiconcordant.\n");
	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	free_imatrix(S, 0, n1-1, 0, m1-1);
	free_imatrix(Sl, 0, n1-1, 0, m1-1);
	return 0;
      }
    }
  }

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  free_imatrix(S, 0, n1-1, 0, m1-1);
  free_imatrix(Sl, 0, n1-1, 0, m1-1);

  if(nonzprod){
    if(debug){fprintf(stderr, "The network is semiconcordant\n");}
    return 1;
  }

  if(debug){fprintf(stderr, "The network is not semiconcordant because it fails to be normal.\n");}
  return 0;
}

//overloading: PN AM input
int semiconcord(int **AM, int n, int m, int debug){
  int **S, **Sl;
  int ret,minus=1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  ret=semiconcord(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading di6 input
bool semiconcord(char *di6, int n, int m, int debug){
  int **S, **Sl;
  int ret,minus=1;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  ret=semiconcord(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//semiaccordance (iThe MA Jacobian matrix is everywhere P0-)
// Sf=stoichiometric matrix; Slf=minus left stoichiometric matrix
int semiaccord(int **Sf, int **Slf, int n, int m, int debug){
  long r1,r2;
  int **S, **Sl;
  int n1,m1;
  int Srank=matrank(Sf,n,m);
  long cnk, cmk;
  int **xcombs, **ycombs;
  int tt;
  int k;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering semiaccord.\n");}
  if(debugfull){fprintf(stderr,"Stoich and left stoich matrices:\n");printmat(Sf,n,m);printmat(Slf,n,m);}

  simppair(Sf, Slf, n, m, &S, &Sl, &n1, &m1); 

  if(debugfull){fprintf(stderr, "Simplified pair:\n");printmat(S, n1, m1);printmat(Sl, n1, m1);}

  for(k=Srank;k>=1;k--){
    xcombs=allcombsgen(n1,k);
    ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	if((tt=fixedminorcompat(S, Sl, n1, m1, xcombs[r1], ycombs[r2],k)) < 0){
	  if(debugfull){
	    fprintf(stderr, "\npair of submatrices with determinants of opposite sign:\n\n");
	    printsubmat(S, xcombs[r1], ycombs[r2], k, k);
	    fprintf(stderr, "*** and ***\n");
	    printsubmat(Sl, xcombs[r1], ycombs[r2], k, k);
	  }
	  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	  free_imatrix(S, 0, n1-1, 0, m1-1);
	  free_imatrix(Sl, 0, n1-1, 0, m1-1);
	  if(debug){fprintf(stderr, "The network fails to be semiaccordant\n");}
	  return 0;
	}
      }
    }
  }

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  free_imatrix(S, 0, n1-1, 0, m1-1);
  free_imatrix(Sl, 0, n1-1, 0, m1-1);
  if(debug){fprintf(stderr, "The network is semiaccordant\n");}
  return 1;
}

int semiaccord(int **AM, int n, int m, int debug){
  int **S, **Sl;
  int ret,minus=1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  ret=semiaccord(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading di6 input
bool semiaccord(char *di6, int n, int m, int debug){
  int **S, **Sl;
  int ret,minus=1;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  ret=semiaccord(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}



// the mass action Jacobian is sign nonsingular


int allminorsigns(int **imat1, int **imat2, int n, int m, int q){
  int k, r, flag;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int val;
  int sgn=0;
  long countp, countn;

  fprintf(stderr, "Checking all minors for mass action Jacobian.\n\n");

  r=min(n,m);

  for(k=1;k<=r;k++){
    sgn=0;
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n,k);cmk=comb(m,k);
    flag=0;countp=0;countn=0;
    //    for(r1=0;(r1<cnk  && flag==0);r1++){
    //      for(r2=0;(r2<cmk && flag==0);r2++){
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	val=fixedminorcompat(imat1, imat2, n, m, xcombs[r1], ycombs[r2],k);
	if(val>0)
	  countp++;
	else if(val<0)
	  countn++;
	if(sgn==0 && val > 0){
	  sgn=1;
	}
	else if(sgn==0 && val < 0){
	  sgn=-1;
	}
	else if((sgn >0 && val < 0) || (sgn <0 && val > 0)){
	  flag=1;
	}
   
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);

    if(flag==1){
      fprintf(stderr, "Minors of dimension %d: %.2f %% (%.0f:%.0f) of all nonzero minors are positive.\n", k, 100.0*((double) countp)/((double)countp + (double)countn), (double) countp, (double) countn);
    }
    else if(sgn==0){
      fprintf(stderr, "Minors of dimension %d are all zero.\n", k);
    }
    else if(sgn==1){
      fprintf(stderr, "Nonzero minors of dimension %d are all positive (%.0f).\n", k, (double) countp);
    }
    else if(sgn==-1){
      fprintf(stderr, "Nonzero minors of dimension %d are all negative (%.0f).\n", k, (double) countn);
    }

  }


  return 1;
}



int MAisSNS(int **imat1, int **imat2, int n, int m, int q){
  int k;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int val;
  int sgn=0;

  fprintf(stderr, "Checking if mass action Jacobian is SNS.\n\n");

  if(m<n)
    return 0;

  k=n;
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n,k);cmk=comb(m,k);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      //      fprintf(stderr, "%.2f percent\n", ((double)r2)/((double) cmk));
      val=fixedminorcompat(imat1, imat2, n, m, xcombs[r1], ycombs[r2],k);
      //      fprintf(stderr, "here\n");
      if(sgn==0 && val > 0)
	sgn=1;
      else if(sgn==0 && val < 0)
	sgn=-1;
      else if((sgn >0 && val < 0) || (sgn <0 && val > 0)){
	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	fprintf(stderr, "Finished checking if mass action Jacobian is SNS: it has terms of opposite sign.\n---------------------------------------\n");
	return 0;
      }
   
    }
  }

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);


  if(sgn==0){
    fprintf(stderr, "Finished checking if mass action Jacobian is SNS: it is sign singular.\n---------------------------------------\n");
    return 0;
  }
  else
    fprintf(stderr, "Finished checking if mass action Jacobian is SNS: it is.\n---------------------------------------\n");

  return 1;
}






int isWS(int **mat, int n){
  int **tmp;
  tmp=submatminus(mat, n, n);


/*   printmat(mat, n, n); */
/*   fprintf(stderr, "\n"); */
/*   printmat(tmp, n, n); */
/*   fprintf(stderr, "det1 = %d,det2 = %d\n", det(mat, n),det(tmp, n)); */
/*   fprintf(stderr, "\n  ****  \n"); */


  if(det(mat, n)*det(tmp, n) >=0){
    free_imatrix(tmp, 0, n-1, 0, n-1);
    return 1;
  }
  free_imatrix(tmp, 0, n-1, 0, n-1);
  return 0;
}

// get rid of empty rows and columns

void simppair(int **imat, int **imat1, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1){
  int i, j, s1, s2, flag;
  int *emptycols;
  int *emptyrows;
  int debug=0;
  (*n1)=0;(*m1)=0;

  emptycols=(int *)malloc((size_t) ((m)*sizeof(int)));
  emptyrows=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<n;i++)
    emptyrows[i]=1;

  for(i=0;i<m;i++)
    emptycols[i]=1;


  i=0;j=0;

  while(j<m){
    i=0;flag=0;
    while(i<n && flag==0){
      if(imat[i][j] !=0){
	flag=1; // not a column of zeros
      }
      i++;
    }
    if(flag==1){
      i=0;
      while(i<n && emptycols[j] == 1){
	if(imat1[i][j] !=0){
	  emptycols[j]=0; // not a column of zeros
	  (*m1)++;
	}
	i++;
      }
    }
    j++;
  }
  if(debug && m!=(*m1))
    fprintf(stderr, "numcols = %d --> %d\n", m, (*m1));

  i=0;j=0;

  while(i<n){
    j=0;flag=0;
    while(j<m && flag==0){
      if(imat[i][j] !=0)
	flag=1; // not a row of zeros
      j++;
    }
    if(flag==1){
      j=0;
      while(j<m && emptyrows[i] == 1){
	if(imat1[i][j] !=0){
	  emptyrows[i]=0; // not a row of zeros
	  (*n1)++;
	}
	j++;
      }
    }
    i++;
  }
  if(debug && n!=(*n1))
    fprintf(stderr, "numrows = %d --> %d\n", n, (*n1));

  if((*n1)==0 || (*m1)==0){ // one matrix is a matrix of zeros
    (*imat2)=NULL;(*imat3)=NULL;
    free((FREE_ARG) (emptyrows));free((FREE_ARG) (emptycols));
    return;
  }


  (*imat2)=imatrix(0, (*n1)-1, 0, (*m1)-1);
  (*imat3)=imatrix(0, (*n1)-1, 0, (*m1)-1);

  s1=0;
  for(i=0;i<n;i++){
    s2=0;
    if(emptyrows[i]==0){ // not empty row
      for(j=0;j<m;j++){
	if(emptycols[j]==0){ // not empty column
	  (*imat2)[s1][s2]=imat[i][j];
	  (*imat3)[s1][s2]=imat1[i][j];
	  s2++;
	}
      }
      s1++;
    }
  }


  free((FREE_ARG) (emptyrows));
  free((FREE_ARG) (emptycols));

  return;
}

void WSpair(int **imat, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1){
  int i, j;
  int **tmpmat;
  (*n1)=0;(*m1)=0;

  tmpmat=imatrix(0, n-1, 0, m-1);

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      tmpmat[i][j]=min(imat[i][j], 0);

  simppair(tmpmat, imat, n, m, imat3, imat2, n1, m1); // tmpmat --> imat3 (sparser matrix)

  free_imatrix(tmpmat, 0, n-1, 0, m-1);
  return;

}



int minorisWS(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k, unsigned long fk){
  int **pms;
  int dt=0, dt2=0;
  unsigned long i;

  // imat2 is sparser
  if(minorhas0rc(imat2, n, m, vec1, vec2, k))
    return 1;

  pms=allperms1(vec2, k);

  for(i=0;i<fk;i++)
    dt2+=unsterm1(imat2, n, m, vec1, pms[i],k)*pms[i][k];

  if(dt2!=0){
    for(i=0;i<fk;i++)
      dt+=unsterm1(imat1, n, m, vec1, pms[i],k)*pms[i][k];

    if(dt*dt2<0){
      //      fprintf(stderr, "dt = %d, dt2=%d\n", dt, dt2);
      free_imatrix(pms,0, fk-1, 0, k);
      return 0;

    }
  }
  free_imatrix(pms,0, fk-1, 0, k);
  return 1;
}

/* This function returns the sign of the product of the  */
/* two k X k submatrices imat1(vec1|vec2) and imat2(vec1|vec2)  */
/* of n X m matrices imat1 and imat2 */

int fixedminorcompat(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k){
  int dt=0, dt2=0, prod;

  if(minorhas0rc(imat2, n, m, vec1, vec2, k) || minorhas0rc(imat1, n, m, vec1, vec2, k))
    return 0;

  dt2=detsubmat(imat2, n, m, vec1, vec2, k);
  if(dt2!=0){
    dt=detsubmat(imat1, n, m, vec1, vec2, k);
    if((prod=dt*dt2)<0)
      return -1;
    else if(prod>0)
      return 1;
    else
      return 0;
  }
  return 0;
}

int minorisnonzero(int **imat1, int n, int m, int *vec1, int *vec2, int k){
  int dt=0;
  if(minorhas0rc(imat1, n, m, vec1, vec2, k) || ((dt=detsubmat(imat1, n, m, vec1, vec2, k))==0))
    return 0;
  if(dt>0)
    return 1;
  return -1;

}


/* change an integer matrix into a subscript matrix */
/* columns first and then rows */

int **mattoind(int **mat, int n, int m){
  int i,j;
  int **tmp;
  int ind=0;
  tmp=imatrix(0, n-1, 0, m-1);
  for(j=0;j<m;j++){
    for(i=0;i<n;i++){
      if(mat[i][j] > 0){
	ind++;tmp[i][j] = ind;
      }
      else if(mat[i][j] < 0){
	ind++;tmp[i][j] = -ind;
      }
      else
	tmp[i][j]=0;
    }
  }
  return tmp;
}

int pairsareequal(int *p1, int *p2, int s1, int *q1, int *q2, int s2){
  int i;
  for (i=0;i<s1;i++){
    if (p1[i]!=p2[i])
      return 0;
  }
  for (i=0;i<s2;i++){
    if (q1[i]!=q2[i])
      return 0;
  }
  return 1;
}



int addtwointstr(int k, int *s, int *s1, int len1, int len2, int ***t)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. */
{
  int i=0,j=0, flag;

  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+2));
    (*t)[k][0] = len1;
    (*t)[k][1] = len2;
    for(i=2;i<2+len1;i++)
      (*t)[k][i] = s[i-2];
    for(i=2+len1;i<2+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-2];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    if(len1!=(*t)[j][0] || len2!=(*t)[j][1]){
      j++;flag=0;continue;
    }
    for(i=2;i<2+len1;i++){
      if((*t)[j][i]!=s[i-2]){
	flag=0;
	break;
      }
    }
    if(flag){
      for(i=2+len1;i<2+len1+len2;i++){
	if((*t)[j][i]!=s1[i-len1-2]){
	  flag=0;
	  break;
	}
      }
    }
    if(flag) // found
      break;
    j++;
  }



  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+2));
    (*t)[k][0] = len1;
    (*t)[k][1] = len2;
    for(i=2;i<2+len1;i++)
      (*t)[k][i] = s[i-2];
    for(i=2+len1;i<2+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-2];
    return k+1;
  }
  return k;
}

int addtwosignedintstr(int k, int *s, int *s1, int len1, int len2, int sgn, int ***t)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. */
{
  int i=0,j=0, flag;

  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+3));
    (*t)[k][0] = sgn;
    (*t)[k][1] = len1;
    (*t)[k][2] = len2;
    for(i=3;i<3+len1;i++)
      (*t)[k][i] = s[i-3];
    for(i=3+len1;i<3+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-3];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    if(len1!=(*t)[j][0] || len2!=(*t)[j][1]){
      j++;flag=0;continue;
    }
    for(i=3;i<3+len1;i++){
      if((*t)[j][i]!=s[i-3]){
	flag=0;
	break;
      }
    }
    if(flag){
      for(i=3+len1;i<3+len1+len2;i++){
	if((*t)[j][i]!=s1[i-len1-3]){
	  flag=0;
	  break;
	}
      }
    }
    if(flag) // found
      break;
    j++;
  }



  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+3));
    (*t)[k][0] = sgn;
    (*t)[k][1] = len1;
    (*t)[k][2] = len2;
    for(i=3;i<3+len1;i++)
      (*t)[k][i] = s[i-3];
    for(i=3+len1;i<3+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-3];
    return k+1;
  }
  return k;
}

int addintstr(int k, int *s, int len1, int ***t, int *ind)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. all strings have length len1 */
{
  int i=0,j=0, flag;
  (*ind)=0;
  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1));
    for(i=0;i<len1;i++)
      (*t)[k][i] = s[i];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    for(i=0;i<len1;i++){
      if((*t)[j][i]!=s[i]){
	flag=0;
	break;
      }
    }
    if(flag){ // found
      (*ind)=j;
      break;
    }
    j++;
  }

  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1));
    for(i=0;i<len1;i++)
      (*t)[k][i] = s[i];
    (*ind)=k;
    return k+1;
  }
  return k;
}

int addoneintstr(int k, int *s, int len1, int ***t)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. The first entry of each member of the list is the length of the string*/
{
  int i=0,j=0, flag;

  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+1));
    (*t)[k][0] = len1;
    for(i=1;i<1+len1;i++)
      (*t)[k][i] = s[i-1];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    if(len1!=(*t)[j][0]){
      j++;flag=0;continue;
    }
    for(i=1;i<1+len1;i++){
      if((*t)[j][i]!=s[i-1]){
	flag=0;
	break;
      }
    }
    if(flag) // found
      break;
    j++;
  }

  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+1));
    (*t)[k][0] = len1;
    for(i=1;i<1+len1;i++)
      (*t)[k][i] = s[i-1];
    return k+1;
  }
  return k;
}



// Check if the set subs of size n1 is a siphon
// by examining the irreversible stoichiometric matrix
// and the matrix of powers

bool issiphon(int **S,int **M,int n,int m,int *subs,int n1){
  int i,j,k;
  int **tmp=cpmat(S,n,m);
  for(i=0;i<n1;i++){// each  member of the siphon
    for(j=0;j<m;j++){
      if(M[subs[i]][j]!=0){//substrate affects reac rate
	//set column of S to zero
	for(k=0;k<n;k++)
	  tmp[k][j]=0;
      }
    }
  }
  for(i=0;i<n1;i++){
    for(j=0;j<m;j++){
      if(tmp[subs[i]][j]!=0){// row not wiped out
	free_imatrix(tmp, 0, n, 0, m-1);
	return 0;
      }
    }
  }

  free_imatrix(tmp, 0, n, 0, m-1);
  return 1;
}

int checksiphons(int **mat1, int **mat2, int n, int m, int ***allsiphons, int *totminsiphons, int ***allminsiphons){

  int j,k;
  bool flg=0;
  int **xcombs;
  long i,cnk;
  int totsiphons=0;
  (*totminsiphons)=0;

  for(k=1;k<=n;k++){// each dimension
    xcombs=allcombsgen(n,k);
    cnk=comb(n,k);
    for(i=0;i<cnk;i++){//each subset

      if(issiphon(mat1,mat2,n,m,xcombs[i],k)){
	// add to minimal list if poss.
	flg=1;
	for(j=0;j<totsiphons;j++){//each existing siphon
	  if(AsubsB((*allsiphons)[j]+1,(*allsiphons)[j][0],xcombs[i],k)){
	    flg=0;break;//not minimal
	  }
	}
	if(flg)
	  (*totminsiphons)=addoneintstr((*totminsiphons),xcombs[i],k,allminsiphons);
	// add to siphon list
	totsiphons=addoneintstr(totsiphons,xcombs[i],k,allsiphons);
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  }


  /* fprintf(stderr, "totsiphons=%d\n", totsiphons); */
  /* for(j=0;j<totsiphons;j++){//each siphon */
  /*   for(k=1;k<=(*allsiphons)[j][0];k++) */
  /*     fprintf(stderr, "%d ", (*allsiphons)[j][k]); */
  /*   fprintf(stderr, "\n"); */
  /* } */

  /* fprintf(stderr, "totminsiphons=%d\n", *totminsiphons); */
  /* for(j=0;j<(*totminsiphons);j++){//each siphon */
  /*   for(k=1;k<=(*allminsiphons)[j][0];k++) */
  /*     fprintf(stderr, "%d ", (*allminsiphons)[j][k]); */
  /*   fprintf(stderr, "\n"); */
  /* } */
  return totsiphons;


}



// remove vertices of degree 1 and their unique neighbour.
// do this recursively until there is no change
// The matrix is assumed to be a submatrix of the original, described
// by row-set rw and column set col, and the output is a reduced
// row set rwout and a reduced column set colout

int cuttails(int **mat, int n, int m, int *rw, int *col, int **rwout, int **colout){
  int i, j, i1, j1;
  int rwcount, colcount, val;
  int ch=1;
  int rwtmp[n];
  int coltmp[m];
  int del=0;

  for(i=0;i<n;i++)
    rwtmp[i]=rw[i];
  for(j=0;j<m;j++)
    coltmp[j]=col[j];


  while(ch){
    ch=0;
    for(i=0;i<n;i++){
      if(rwtmp[i]>=0){
	rwcount=0;
	for(j=0;j<m;j++){ // count entries in row
	  if(coltmp[j]>=0){
	    if(mat[rw[i]][col[j]]!=0){rwcount++;val=j;}
	  }
	}
	if(rwcount==1){rwtmp[i]=-1;coltmp[val]=-1;ch=1;del++;/* fprintf(stderr, "cut: %d, %d\n", i, val); */} //singleton
      }
    }
    for(j=0;j<m;j++){
      if(coltmp[j]>=0){
	colcount=0;
	for(i=0;i<n;i++){ // count entries in column
	  if(rwtmp[i]>=0){
	    if(mat[rw[i]][col[j]]!=0){colcount++;val=i;}
	  }
	}
	if(colcount==1){rwtmp[val]=-1;coltmp[j]=-1;ch=1;del++;/* fprintf(stderr, "cut1: %d, %d\n", val, j); */}//singleton
      }
    }
  }
  //  fprintf(stderr, "del = %d\n", del);

  (*rwout)=(int *)malloc((size_t) ((n-del)*sizeof(int)));
  (*colout)=(int *)malloc((size_t) ((m-del)*sizeof(int)));
  i1=0;j1=0;
  for(i=0;i<n;i++){
    if(rwtmp[i]>=0){
      (*rwout)[i1]=rw[i];i1++;
    }
  }
  //  fprintf(stderr, "i1 = %d\n", i1);
  for(j=0;j<m;j++){
    if(coltmp[j]>=0){
      (*colout)[j1]=col[j];j1++;
    }
  }
  //  fprintf(stderr, "j1 = %d\n", j1);
  return del;

}

// remove rows/columns with one or fewer nonzero entries

int **simpmat(int **mat, int n, int m, int *n1, int *m1){
  int i,j,k,l,flag;
  int **tmp;
  int x[n], y[m];
  *n1=0;*m1=0;
  for(i=0;i<n;i++){
    flag=0;
    for(j=0;j<m;j++){ // count entries in row
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=2){x[i]=1;(*n1)++;}
    else{x[i]=0;}
  }

  for(j=0;j<m;j++){
    flag=0;
    for(i=0;i<n;i++){ // count entries in column
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=2){y[j]=1;(*m1)++;}
    else{y[j]=0;}
  }
  tmp=imatrix(0, (*n1)-1, 0, (*m1)-1);

  k=0;
  for(i=0;i<n;i++){
    if(x[i]>0){
      l=0;
      for(j=0;j<m;j++){
	if(y[j] > 0){
	  tmp[k][l] = mat[i][j];
	  l++;
	}
      }
      k++;
    }
  }

  return tmp;
}


// remove empty rows/columns with one or fewer nonzero entries

int **vsimpmat(int **mat, int n, int m, int *n1, int *m1){
  int i,j,k,l,flag;
  int **tmp;
  int x[n], y[m];
  *n1=0;*m1=0;
  for(i=0;i<n;i++){
    flag=0;
    for(j=0;j<m;j++){ // count entries in row
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=1){x[i]=1;(*n1)++;}
    else{x[i]=0;}
  }

  for(j=0;j<m;j++){
    flag=0;
    for(i=0;i<n;i++){ // count entries in column
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=1){y[j]=1;(*m1)++;}
    else{y[j]=0;}
  }
  tmp=imatrix(0, (*n1)-1, 0, (*m1)-1);

  k=0;
  for(i=0;i<n;i++){
    if(x[i]>0){
      l=0;
      for(j=0;j<m;j++){
	if(y[j] > 0){
	  tmp[k][l] = mat[i][j];
	  l++;
	}
      }
      k++;
    }
  }

  return tmp;
}




// the symbolic determinant of a symbolic matrix
// without expansion of the answer

ex symbdetsubmat_noexp(matrix imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  ex tot=0;
  //ex extmp;
  int i2a[dim-1];
  if(dim==1)
    return imat(i1[0], i2[0]);
  /* if(dim>4) */
  /*   cout << "dim= " << dim << endl; */
  for (i=0;i<dim;i++){

    if(imat(i1[0], i2[i])!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1);
    /*   extmp=symbdetsubmat(imat, n, m, i1+1, i2a, dim-1); */
    /*   if(!extmp.info(info_flags::positive)) */
    /* 	tot+=par(i)*imat(i1[0],i2[i])*expand(extmp); */
    /*   else{ */
    /* 	tot+=par(i)*imat(i1[0],i2[i])*extmp; */
    /* 	//cout << "signed expr:\n" << extmp << endl; */
    /*   } */
    /* 	//      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1); */
    /*   //      cout << tot << endl; */
    /* } */

    }
  }

  return tot;

}

//
// swap ith and jth columns of integer matrix *v n is the column length
void colswap(int ***v, int n, int i, int j){
  int k;
  int temp;
  for(k=0;k<n;k++){
    temp=(*v)[k][i];
    (*v)[k][i]=(*v)[k][j];
    (*v)[k][j]=temp;
  }
  return;
}

//overloading (boolean matrix)
void colswap(bool ***v, int n, int i, int j){
  int k;
  bool temp;
  for(k=0;k<n;k++){
    temp=(*v)[k][i];
    (*v)[k][i]=(*v)[k][j];
    (*v)[k][j]=temp;
  }
  return;
}

void colswapmat(matrix *v, int n, int i, int j){
  int k;
  ex temp;
  for(k=0;k<n;k++){
    temp=(*v)(k,i);
    (*v)(k,i)=(*v)(k,j);
    (*v)(k,j)=temp;
  }
  return;
}


//lex ordering: return 1 if a>b, -1 if a>b, 0 otherwise
int cmpare(int *a, int *b, int n){
  int i;
  for(i=0;i<n;i++){
    if(a[i]<b[i])
      return -1;
    else if(a[i]>b[i])
      return 1;
  }
  return 0;
}



//
// reverse lexicographic ordering
//

int cmparerev(int *a, int *b, int n){
  int i;
  for(i=n-1;i>=0;i--){
    if(a[i]<b[i])
      return -1;
    else if(a[i]>b[i])
      return 1;
  }
  return 0;
}

//overloading
int cmparerev(bool *a, bool *b, int n){
  int i;
  for(i=n-1;i>=0;i--){
    if(a[i]<b[i])
      return -1;
    else if(a[i]>b[i])
      return 1;
  }
  return 0;
}

// n is the column length

void lexsortcol(int ***v, int n, int left, int right)
{
  int i, last;

  if (left >= right)
    return;
  colswap(v,n,left,(left + right)/2);
  last = left;
  for (i=left+1; i<=right;i++)
    if (cmparerev((*v)[i],(*v)[left],n)<0)
      colswap(v,n,++last,i);
  colswap(v,n,left,last);
  lexsortcol(v,n,left,last-1);
  lexsortcol(v,n,last+1,right);
  return;
}

void lexsortcol4(bool ***v, bool ***v1, int ***imat, matrix *mmat, int n, int left, int right)
{
  int i,k,last;
  bool vleft[n];
  bool vi[n];

  if (left >= right)
    return;
  colswap(v,n,left,(left + right)/2);
  colswap(v1,n,left,(left + right)/2);
  colswap(imat,n,left,(left + right)/2);
  colswapmat(mmat,n,left,(left + right)/2);
  last = left;


  for (i=left+1; i<=right;i++){
    for(k=0;k<n;k++)
      vi[k]=(*v)[k][i];
    for(k=0;k<n;k++)
      vleft[k]=(*v)[k][left];
    if (cmparerev(vi,vleft,n)<0){
      last++;
      colswap(v,n,last,i);
      colswap(v1,n,last,i);
      colswap(imat,n,last,i);
      colswapmat(mmat,n,last,i);
    }
  }

  colswap(v,n,left,last);
  colswap(v1,n,left,last);
  colswap(imat,n,left,last);
  colswapmat(mmat,n,left,last);

  lexsortcol4(v,v1,imat,mmat,n,left,last-1);
  lexsortcol4(v,v1,imat,mmat,n,last+1,right);
  return;
}

// Reorder the columns of a bool matrix, an int matrix and 
// an ex matrix putting them in upwards lexicographic order

void multireordercollex(bool **bmatin, bool ***bmatout, bool **bmat2in, bool ***bmat2out, int **imatin, int ***imatout, matrix mmatin, matrix *mmatout, int n, int m){
  int i,j;
  (*bmatout)=bmatrix(0, n-1, 0, m-1);
  (*bmat2out)=bmatrix(0, n-1, 0, m-1);
  (*imatout)=imatrix(0, n-1, 0, m-1);
  (*mmatout)=matrix(n,m);
  //copy
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      (*bmatout)[i][j]=bmatin[i][j];
      (*bmat2out)[i][j]=bmat2in[i][j];
      (*imatout)[i][j]=imatin[i][j];
      (*mmatout)(i,j)=mmatin(i,j);
    }
  }
  lexsortcol4(bmatout,bmat2out,imatout,mmatout,n,0,m-1);
  return;
}



//
// extract the coefficient of a term from the determinant of a 
// matrix. The routine is very inefficient: 
// for a dense matrix, can be almost as bad as computing a determinant
int coeff(matrix J2, int *i1, int *i2, int dim, lst trm){
  int tot=0, j, l, m,r;
  int i2a[dim-1];
  lst trm1;
  int sgn;
  if(dim>12)
  cout<<"calling: dim = " << dim << endl;

  if(dim==1){
    //fprintf(stderr, "dim=1\n");
    if(is_a<add>(J2(i1[0],i2[0]))){ // a sum of terms
      //fprintf(stderr, "multiterm\n");
      for (size_t k=0; k!=J2(i1[0],i2[0]).nops(); ++k){
	if(trm[0]==J2(i1[0],i2[0]).op(k)){
	  //cout << "found" << endl;
	  return 1;
	}
	else if(trm[0]==-J2(i1[0],i2[0]).op(k))
	  return -1;
      }
    }
    else{ // a single term
      //fprintf(stderr, "single termm\n");
      if(trm[0]==J2(i1[0],i2[0]))
	return 1;
      else if(trm[0]==-J2(i1[0],i2[0]))
	return -1;
    }
    return 0;
  }

  //fprintf(stderr, "dim=%d\n", dim);
  for(j=0;j<dim;j++){
    if(J2(i1[0],i2[j])==0){}
    else if(is_a<add>(J2(i1[0],i2[j]))){
      //fprintf(stderr, "multiterm\n");
      for (size_t k = 0; k != J2(i1[0],i2[j]).nops(); ++k){
      	//cout << "k=" <<k<< endl;
      	l=0;sgn=0;
      	while(l<dim && sgn==0){
      	  if(trm[l]==J2(i1[0],i2[j]).op(k))
      	    sgn=1;
      	  else if(trm[l]==-J2(i1[0],i2[j]).op(k))
      	    sgn=-1;
	  if(sgn){
	    r=0;
	    for(m=0;m<dim;m++){
	      if(m!=j)
		i2a[r++]=i2[m];
	    }
	    trm1=lst();
	    for(m=0;m<dim;m++){
	      if(m!=l)
		trm1.append(trm[m]);
	    }
	    //cout << "trm1 = " << trm1 << endl;
	    tot+=par(j)*sgn*coeff(J2, i1+1, i2a, dim-1, trm1);
	    break; // found in list
	  }
	  l++;
	}
      }
    }
    else{
      //fprintf(stderr, "single term\n");
      l=0;sgn=0;
      while(l<dim && sgn==0){
	//cout << "trm[l] = " << trm[l] << endl;
	//cout << "J2(i1[0],i2[j]) = " << J2(i1[0],i2[j]) << endl;
	if(trm[l]==J2(i1[0],i2[j])){
	  //cout << "found_single\n";
	  sgn=1;
	}
	else if(trm[l]==-J2(i1[0],i2[j]))
	  sgn=-1;
	if(sgn){// lth member found
	  r=0;
	  for(m=0;m<dim;m++){
	    if(m!=j)
	      i2a[r++]=i2[m];
	  }
	  trm1=lst();
	  for(m=0;m<dim;m++){
	    if(m!=l)
	      trm1.append(trm[m]);
	  }
	  //printvec(i2a,dim-1);
	  //cout << "trm1 = " << trm1 << endl;
	  tot+=par(j)*sgn*coeff(J2, i1+1, i2a, dim-1, trm1);
	  break; // found in list
	}
	l++;
      }

    }

  }
  return tot;
}

// recursive function to:
// extract the coefficient of a term from the determinant of a matrix
// each entry in the matrix is a list of integers
// the first entry in each list, is the length
// subsequent entries can be seen as symbols

int coeff_int(int ****J2, int *i1, int *i2, int dim, int *trm){
  int tot=0, k, j, l, m,r;
  int i2a[dim-1];
  int trm1[dim-1];
  int sgn;
  //cout<<"calling\n";

  if(dim==1){
    //cout << "in here\n";
    for (k=1; k<=(*J2)[i1[0]][i2[0]][0]; k++){
      if(trm[0]==(*J2)[i1[0]][i2[0]][k]){
	//cout << "found" << endl;
	return 1;
      }
      else if(trm[0]==-(*J2)[i1[0]][i2[0]][k])
	return -1;
    }
    return 0;
  }


  //fprintf(stderr, "dim=%d\n", dim);
  for(j=0;j<dim;j++){
    //cout << "in here1 dim=" << dim << endl;
    if((*J2)[i1[0]][i2[j]][0]==0){}
    else{
      //cout << "i,j " << i1[0] << "," << i2[j] << endl;
      for (k=1; k<=(*J2)[i1[0]][i2[j]][0]; k++){
	//cout << trm[0] << endl;
	//cout << "i,j,k= " << i1[0] << "," << i2[j] << "," << k << endl;
	//cout << "in here2 (*J2)[i1[0]][i2[j]][k]=" << (*J2)[i1[0]][i2[j]][k] << endl;
      	//cout << "k=" <<k<< endl;
      	l=0;sgn=0;
      	while(l<dim && sgn==0){
      	  if(trm[l]==(*J2)[i1[0]][i2[j]][k])
      	    sgn=1;
      	  else if(trm[l]==-(*J2)[i1[0]][i2[j]][k])
      	    sgn=-1;
	  if(sgn){
	    r=0;
	    for(m=0;m<dim;m++){
	      if(m!=j)
		i2a[r++]=i2[m];
	    }
	    r=0;
	    for(m=0;m<dim;m++){
	      if(m!=l)
		trm1[r++]=trm[m];
	    }

	    //cout << "(*J2)[i1[0]][i2[j]][k] =" << (*J2)[i1[0]][i2[j]][k] << endl;
	    //cout << "trm[l] =" << trm[l] << endl;
	    //cout << "trm1 = " << trm1 << endl;
	    tot+=par(j)*sgn*coeff_int(J2, i1+1, i2a, dim-1, trm1);
	    break; // found in list
	  }
	  l++;
	}
      }
    }
  }
  return tot;
}

// routine to check the coeff function on a large matrix 
// whose determinant has been evaluated with GINAC

void check_coeff_fn(){
  int k=15;
  matrix trial(k,k);
  int xc[k], yc[k];
  lst trm;
  possymbol v1("v1"), v2("v2"), v3("v3"), v4("v4"), v5("v5"), v6("v6"), v7("v7"), v8("v8"), v9("v9"), v10("v10"), v11("v11");
  trial(0,0)=v2+v1+v3; trial(0,1)=-v4; trial(0,2)=v6; trial(0,5)=v6;
  trial(1,0)=-v3;trial(1,1)= v4+v5+v1;trial(1,2)=-v7;trial(1,5)=-v2;trial(1,8)=v6;
  trial(2,0)=v2;trial(2,1)=-v5;trial(2,2)=v7+v8+v1+v6;trial(2,3)=-v9;trial(2,6)=-v2;
  trial(3,2)=-v8;trial(3,3)=v10+v1+v9;trial(3,4)=-v11;trial(3,7)=-v2;trial(3,12)=-v6;
  trial(4,3)=-v10;trial(4,4)=v11+v1;trial(4,8)=-v2;trial(4,13)= -v6;
  trial(5,1)=-v1;trial(5,5)=v4+v2+v5+v3;trial(5,6)=-v7;trial(5,9)=-v6;
  trial(6,0)=v1;trial(6,2)=-v1;trial(6,5)=-v5;trial(6,6)=v7+v2+v8+v6+v3;trial(6,7)=-v9;trial(6,9)=-v4;
  trial(7,3)=-v1;trial(7,6)=-v8;trial(7,7)=v10+v2+v9+v3;trial(7,8)=-v11;trial(7,10)=-v4;trial(7,12)=v6;
  trial(8,4)=-v1;trial(8,7)=-v10;trial(8,8)=v2+v11+v3;trial(8,11)=-v4;trial(8,13)=v6;
  trial(9,1)=v1;trial(9,5)=-v2;trial(9,6)=-v3;trial(9,9)=v7+v4+v8+v5+v6;trial(9,10)=-v9;
  trial(10,7)=-v3;trial(10,9)=-v8;trial(10,10)=v10+v4+v5+v9;trial(10,11)=-v11;trial(10,12)=-v7;
  trial(11,8)=-v3;trial(11,10)=-v10;trial(11,11)=v4+v11+v5;trial(11,13)=-v7;
  trial(12,3)=-v1;trial(12,7)=v2;trial(12,10)=-v5;trial(12,12)=v10+v7+v8+v9+v6;trial(12,13)=-v11;
  trial(13,4)=-v1;trial(13,8)=v2;trial(13,11)=-v5;trial(13,12)=-v10;trial(13,13)=v7+v11+v8+v6;trial(13,14)=-v9; 
  trial(14,13)=-v8;trial(14,14)=v10+v11+v9;
  printmat1(trial, k, k);
  firstcomb(xc, k, k);
  firstcomb(yc, k, k);
  printvec(xc,k);
  printvec(yc,k);
  trm = monotolist(pow(v10,2)*pow(v7,3)*pow(v4,2)*pow(v2,3)*v11*pow(v8,2)*v1*v9);
  cout << trm << endl;
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  trm.remove_all();
  //120*v10^2*v7*v11^2*v8^2*v5^3*v1^2*v3^3
  trm = monotolist(-v10*v10*v7*v11*v11*pow(v8,2)*pow(v5,3)*pow(v1,2)*pow(v3,3));
  //trm = v10,v10,v7,v11,v11,v8,v8,v5,v5,v5,v1,v1,v3,v3,v3; // expect 120
  cout << trm << endl;
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  //48*v7*v4*v2*v11^2*v5^4*v9^2*v6^2*v3^2
  trm.remove_all();
  trm = {v7,v4,v2,v11,v11,v5,v5,v5,v5,v9,v9,v6,v6,v3,v3}; // expect 48
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  //v8^2*v5^2*v1^3*v9*v3*v10*v7^2*v4*v2*v11
  trm.remove_all();
  trm = {v8,v8,v5,v5,v1,v1,v1,v9,v3,v10,v7,v7,v4,v2,v11}; // expect 267
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;

  return;
}

void check_coeff_int(){
  int k=15;
  int i, j, m;
  int ***trial;
  trial =(int ***) malloc((size_t)(k*sizeof(int**)));
  for(i=0;i<k;i++){
    trial[i]=(int **) malloc((size_t)(k*sizeof(int*)));
    for(j=0;j<k;j++){
      trial[i][j]=(int *) malloc((size_t)(6*sizeof(int)));
    }
  }
 
  int *xc, *yc;
  int trm[15]={10,10,7,11,11,8,8,5,5,5,1,1,3,3,3};
  for(i=0;i<k;i++){
    for(j=0;j<k;j++){
      for(m=0;m<6;m++){
	trial[i][j][m]=0;
      }
    }
  }
 
  trial[0][0][0]=3;trial[0][0][1]=1; trial[0][0][2]=2; trial[0][0][3]=3; 
  trial[0][1][0]=1;trial[0][1][1]=-4; 
  trial[0][2][0]=1;trial[0][2][1]=6; 
  trial[0][5][0]=1;trial[0][5][1]=6;

  trial[1][0][0]=1;trial[1][0][1]=-3;
  trial[1][1][0]= 3;trial[1][1][1]= 1;trial[1][1][2]= 4;trial[1][1][3]= 5;
  trial[1][2][0]=1;trial[1][2][1]=-7;
  trial[1][5][0]=1;trial[1][5][1]=-2;
  trial[1][8][0]=1;trial[1][8][1]=6;

  trial[2][0][0]=1;trial[2][0][1]=2;
  trial[2][1][0]=1;trial[2][1][1]=-5;
  trial[2][2][0]=4;trial[2][2][1]=7;trial[2][2][2]=8;trial[2][2][3]=1;trial[2][2][4]=6;
  trial[2][3][0]=1;trial[2][3][1]=-9;
  trial[2][6][0]=1;trial[2][6][1]=-2;

  trial[3][2][0]=1;trial[3][2][1]=-8;
  trial[3][3][0]=3;trial[3][3][1]=10;trial[3][3][2]=1;trial[3][3][3]=9;
  trial[3][4][0]=1;trial[3][4][1]=-11;
  trial[3][7][0]=1;trial[3][7][1]=-2;
  trial[3][12][0]=1;trial[3][12][1]=-6;

  trial[4][3][0]=1;trial[4][3][1]=-10;
  trial[4][4][0]=2;trial[4][4][1]=11;trial[4][4][2]=1;
  trial[4][8][0]=1;trial[4][8][1]=-2;
  trial[4][13][0]=1;trial[4][13][1]= -6;

  trial[5][1][0]=1;trial[5][1][1]=-1;
  trial[5][5][0]=4;trial[5][5][1]=4;trial[5][5][2]=2;trial[5][5][3]=5;trial[5][5][4]=3;
  trial[5][6][0]=1;trial[5][6][1]=-7;
  trial[5][9][0]=1;trial[5][9][1]=-6;

  trial[6][0][0]=1;trial[6][0][1]=1;
  trial[6][2][0]=1;trial[6][2][1]=-1;
  trial[6][5][0]=1;trial[6][5][1]=-5;
  trial[6][6][0]=5;trial[6][6][1]=7;trial[6][6][2]=2;trial[6][6][3]=8;trial[6][6][4]=6;trial[6][6][5]=3;
  trial[6][7][0]=1;trial[6][7][1]=-9;
  trial[6][9][0]=1;trial[6][9][1]=-4;

  trial[7][3][0]=1;trial[7][3][1]=-1;
  trial[7][6][0]=1;trial[7][6][1]=-8;
  trial[7][7][0]=4;trial[7][7][1]=10;trial[7][7][2]=2;trial[7][7][3]=9;trial[7][7][4]=3;
  trial[7][8][0]=1;trial[7][8][1]=-11;
  trial[7][10][0]=1;trial[7][10][1]=-4;
  trial[7][12][0]=1;trial[7][12][1]=6;

  trial[8][4][0]=1;trial[8][4][1]=-1;
  trial[8][7][0]=1;trial[8][7][1]=-10;
  trial[8][8][0]=3;trial[8][8][1]=2;trial[8][8][2]=11;trial[8][8][3]=3;
  trial[8][11][0]=1;trial[8][11][1]=-4;
  trial[8][13][0]=1;trial[8][13][1]=6;

  trial[9][1][0]=1;trial[9][1][1]=1;
  trial[9][5][0]=1;trial[9][5][1]=-2;
  trial[9][6][0]=1;trial[9][6][1]=-3;
  trial[9][9][0]=5;trial[9][9][1]=7;trial[9][9][2]=4;trial[9][9][3]=8;trial[9][9][4]=5;trial[9][9][5]=6;
  trial[9][10][0]=1;trial[9][10][1]=-9;

  trial[10][7][0]=1;trial[10][7][1]=-3;
  trial[10][9][0]=1;trial[10][9][1]=-8;
  trial[10][10][0]=4;trial[10][10][1]=10;trial[10][10][2]=4;trial[10][10][3]=5;trial[10][10][4]=9;
  trial[10][11][0]=1;trial[10][11][1]=-11;
  trial[10][12][0]=1;trial[10][12][1]=-7;

  trial[11][8][0]=1;trial[11][8][1]=-3;
  trial[11][10][0]=1;trial[11][10][1]=-10;
  trial[11][11][0]=3;trial[11][11][1]=4;trial[11][11][2]=11;trial[11][11][3]=5;
  trial[11][13][0]=1;trial[11][13][1]=-7;

  trial[12][3][0]=1;trial[12][3][1]=-1;
  trial[12][7][0]=1;trial[12][7][1]=2;
  trial[12][10][0]=1;trial[12][10][1]=-5;
  trial[12][12][0]=5;trial[12][12][1]=10;trial[12][12][2]=7;trial[12][12][3]=8;trial[12][12][4]=9;trial[12][12][5]=6;
  trial[12][13][0]=1;trial[12][13][1]=-11;

  trial[13][4][0]=1;trial[13][4][1]=-1;
  trial[13][8][0]=1;trial[13][8][1]=2;
  trial[13][11][0]=1;trial[13][11][1]=-5;
  trial[13][12][0]=1;trial[13][12][1]=-10;
  trial[13][13][0]=4;trial[13][13][1]=7;trial[13][13][2]=11;trial[13][13][3]=8;trial[13][13][4]=6;
  trial[13][14][0]=1;trial[13][14][1]=-9;

  trial[14][13][0]=1;trial[14][13][1]=-8;
  trial[14][14][0]=3;trial[14][14][1]=10;trial[14][14][2]=11;trial[14][14][3]=9;

  xc=(int *)malloc((size_t) (k*sizeof(int)));
  yc=(int *)malloc((size_t) (k*sizeof(int)));
  firstcomb(xc, k, k);
  firstcomb(yc, k, k);
  printvec(xc,k);
  printvec(yc,k);
  //120*v10^2*v7*v11^2*v8^2*v5^3*v1^2*v3^3
  //trm={10,10,7,11,11,8,8,5,5,5,1,1,3,3,3};
  //trm = v10,v10,v7,v11,v11,v8,v8,v5,v5,v5,v1,v1,v3,v3,v3; // expect 120
  //cout << *trm << endl;
  cout << "coeff = " << coeff_int(&trial, xc, yc, k, trm) << endl;
  /* //48*v7*v4*v2*v11^2*v5^4*v9^2*v6^2*v3^2 */
  /* trm.remove_all(); */
  /* trm = v7,v4,v2,v11,v11,v5,v5,v5,v5,v9,v9,v6,v6,v3,v3; // expect 48 */
  /* cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl; */
  /* //v8^2*v5^2*v1^3*v9*v3*v10*v7^2*v4*v2*v11 */
  /* trm.remove_all(); */
  /* trm = v8,v8,v5,v5,v1,v1,v1,v9,v3,v10,v7,v7,v4,v2,v11; // expect 267 */
  /* cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl; */
  free((char *) xc);free((char *) yc);

  return;
}



// qualitative + (inputs/outputs are 0,1,-1 and 2 where 2 means unsigned)
// No checking of input
int qualp(int a, int b){
  if(a==0 || b==0)
    return a+b;
  if(a==1 && b==1)
    return 1;
  if(a==-1 && b==-1)
    return -1;
  if(a==2 || b==2 || (a*b == -1))
    return 2;
  return 0;//never reaches here

}


// qualitative X (inputs/outputs are 0,1,-1 and 2 where 2 means unsigned)
int qualt(int a, int b){
  if(a==0 || b==0) // takes priority
    return 0;
  if(a==2 || b==2)
    return 2;
  return a*b;
}



// the qualitative determinant (outputs are 
// 0,1,-1 and 2 where 2 means unsigned)
// applies only to a real matrix which defines a sign pattern
// namely, not a matrix pattern
// namely, no unsigned entries allowed
// probably really inefficient

int qualdetsubmat(int **imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  int tot=0;
  int i2a[dim-1];
  /* assume this check has already been done. */
  /* if(minorhas0rc(imat, n, m, i1, i2, dim)) */
  /*   return 0; */
  if(dim==1)
    return sgn(imat[i1[0]][i2[0]]);
  for (i=0;i<dim;i++){
    if(imat[i1[0]][i2[i]]!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot=qualp(tot,qualt(par(i),qualt(sgn(imat[i1[0]][i2[i]]),qualdetsubmat(imat, n, m, i1+1, i2a, dim-1))));
    }

  }
  /* fprintf(stderr, "qualdet=%d\n", tot); */
  /* printsubmat(imat, i1, i2, dim, dim); */
  return tot;

}


// the qualitative determinant (outputs are 
// 0,1,-1 and 2 where 2 means unsigned)
// applies to a matrix pattern, namely entries must be
// only 0,1,-1, and 2, where 2 means unsigned 

int qualdetsubmat1(int **imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  int tot=0;
  int i2a[dim-1];
  /* assume this check has already been done. */
  /* if(minorhas0rc(imat, n, m, i1, i2, dim)) */
  /*   return 0; */
  if(dim==1)
    return imat[i1[0]][i2[0]];
  for (i=0;i<dim;i++){
    if(imat[i1[0]][i2[i]]!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot=qualp(tot,qualt(par(i),qualt(imat[i1[0]][i2[i]],qualdetsubmat1(imat, n, m, i1+1, i2a, dim-1))));
    }

  }
  /* fprintf(stderr, "qualdet=%d\n", tot); */
  /* printsubmat(imat, i1, i2, dim,dim); */
  return tot;

}



// here vec must be one larger than r

long getpos1(int *vec, int omit, int n, int r){
  long tot=0;
  int i,j, k1=n-1, k2=r-2, k3=0;
  for(j=0;j<r;j++){
    if(j!=omit){
      for(i=0;i<vec[j]-k3;i++){
	tot+=comb(k1, k2);
	k1--;
      }
      k1--;k2--;k3=vec[j]+1;
    }
  }
  return tot;
}

int ***allminors(int **imat, int n, int m){
  int i, s;//largest minors
  int ***tmp;
  s=min(n,m);


  tmp=(int ***) malloc((size_t)((s)*sizeof(int**)));
  if (!tmp) {fprintf(stderr, "allocation failure.\n");return NULL;}

  tmp[0]=cpmat(imat, n, m);
  for (i=1;i<s;i++){
    fprintf(stderr, "s = %d, i=%d***\n", s, i);
    tmp[i]=detsk1(imat, n, m, i+1, tmp[i-1]);
  }
  return tmp;

}

void free_allminors(int ***t, int n, int m)
{
  int i, s;
  s=min(n,m);
  for (i=0;i<s;i++){
    free_imatrix(t[i], 0, comb(n, i), 0, comb(m, i));
  }
  free((FREE_ARG) (t));
}


// calculates the minor of n times m matrix imat 
// indexed by xcombs and ycombs (which are of sign k)
// assumes that there is a matrix of minors (dets) of
// submatrices of imat of size k-1

int minor1(int **imat, int *xcombs, int *ycombs, int n, int m, int k, int **dets){
  int i;
  long i1,i2;
  int tot=0;
  fprintf(stderr, "xcombs = ");
  for(i=0;i<k;i++){
    fprintf(stderr, "%d ", xcombs[i]);
  }
  fprintf(stderr, "ycombs = ");
  for(i=0;i<k;i++){
    fprintf(stderr, "%d ", ycombs[i]);
  }
  fprintf(stderr, "\n");

  i1=getpos1(xcombs, 0, n, k);
  fprintf(stderr, "i1 = %li\n", i1);
  for (i=0;i<k;i++){
    i2=getpos1(ycombs, i, m, k);
    fprintf(stderr, "i1, i2 = %li, %li, entry=%d, dets = %d\n", i1, i2, imat[xcombs[0]][ycombs[i]], dets[i1][i2]);
    tot+=par(i)*imat[xcombs[0]][ycombs[i]]*dets[i1][i2];
  }
  fprintf(stderr, "tot = %d\n", tot);
  return tot;

}



/* generate both factors in the clever factorisation of J[2] */
/* The second factor is stored as a set of signed integers */
/* which can be converted to strings when needed */

int genS2(int **imat, int n, int m, int ***mat1, int ***mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;
  int **tmp;

  tmp = mattoind(imat, n, m);

  cnk=n*(n-1)/2;cmk=n*m;
  (*mat1)=imatrix(0, cnk-1, 0, cmk-1);
  (*mat2)=imatrix(0, cnk-1, 0, cmk-1);
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      /* fprintf(stderr, "%d, %d, %d\n", r2%n,xcombs[r1][0],xcombs[r1][1]); */
      if(r2%n==xcombs[r1][0]){
	(*mat1)[r1][r2] = -imat[xcombs[r1][1]][r2/n];
	(*mat2)[r1][r2] = -tmp[xcombs[r1][1]][r2/n];
      }
      else if (r2%n==xcombs[r1][1]){
	(*mat1)[r1][r2] = imat[xcombs[r1][0]][r2/n];
	(*mat2)[r1][r2] = tmp[xcombs[r1][0]][r2/n];
      }
      else{
	(*mat1)[r1][r2] = 0;
	(*mat2)[r1][r2] = 0;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  free_imatrix(tmp, 0, n-1, 0, m-1);
  return 1;
}


int genS2symb(matrix mat_in, int n, int m, matrix *mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;

  cnk=n*(n-1)/2;cmk=n*m;
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      if(r2%n==xcombs[r1][0])
	(*mat2)(r1,r2) = -mat_in(xcombs[r1][1],r2/n);
      else if (r2%n==xcombs[r1][1])
	(*mat2)(r1,r2) = mat_in(xcombs[r1][0],r2/n);
      else
	(*mat2)(r1,r2) = 0;
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  return 1;
}

/* Given a factorisation of J into (imatS)(imatV)^t, generate 
// both factors in the natural factorisation of J[2] */
int genS2from2(int **imatS, int **imatV, int n, int m, int ***mat1, int ***mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;

  cnk=n*(n-1)/2;cmk=n*m;
  (*mat1)=imatrix(0, cnk-1, 0, cmk-1);
  (*mat2)=imatrix(0, cnk-1, 0, cmk-1);
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      /* fprintf(stderr, "%d, %d, %d\n", r2%n,xcombs[r1][0],xcombs[r1][1]); */
      if(r2%n==xcombs[r1][0]){
	(*mat1)[r1][r2] = -imatS[xcombs[r1][1]][r2/n];
	(*mat2)[r1][r2] = -imatV[xcombs[r1][1]][r2/n];
      }
      else if (r2%n==xcombs[r1][1]){
	(*mat1)[r1][r2] = imatS[xcombs[r1][0]][r2/n];
	(*mat2)[r1][r2] = imatV[xcombs[r1][0]][r2/n];
      }
      else{
	(*mat1)[r1][r2] = 0;
	(*mat2)[r1][r2] = 0;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  return 1;
}

//overloading - where V is a symbolic matrix
int genS2from2(int **imatS, matrix imatV, int n, int m, int ***mat1, matrix *mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;

  cnk=n*(n-1)/2;cmk=n*m;
  (*mat1)=imatrix(0, cnk-1, 0, cmk-1);
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      /* fprintf(stderr, "%d, %d, %d\n", r2%n,xcombs[r1][0],xcombs[r1][1]); */
      if(r2%n==xcombs[r1][0]){
	(*mat1)[r1][r2] = -imatS[xcombs[r1][1]][r2/n];
	(*mat2)(r1,r2) = -imatV(xcombs[r1][1],r2/n);
      }
      else if (r2%n==xcombs[r1][1]){
	(*mat1)[r1][r2] = imatS[xcombs[r1][0]][r2/n];
	(*mat2)(r1,r2) = imatV(xcombs[r1][0],r2/n);
      }
      else{
	(*mat1)[r1][r2] = 0;
	(*mat2)(r1,r2) = 0;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  return 1;
}



//kth multiplicative compound of imat
int **detsk1(int **imat, int n, int m, int k, int **dets){
  int **xcombs;
  int **ycombs;
  int **tmp;
  long r1, r2,cnk,cmk;


  tmp=imatrix(0, comb(n, k)-1, 0, comb(m, k)-1);
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n,k);cmk=comb(m,k);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      tmp[r1][r2]=minor1(imat, xcombs[r1], ycombs[r2], n, m, k, dets);
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  return tmp;

}



void printmaximaindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2){
  int i,j;
  fprintf(fd, "MM:matrix(");
  for(i=0;i<k1;i++){
    fprintf(fd, "[");
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0)
	fprintf(fd, "0");
      else if(imat[vec1[i]][vec2[j]]<0)
	fprintf(fd, "-v%d",-imat[vec1[i]][vec2[j]]);
      else
	fprintf(fd, "v%d",imat[vec1[i]][vec2[j]]);
      if(j<(k2-1))
	fprintf(fd, ",");
    }
    if(i<k1-1)
      fprintf(fd, "],");
    else
      fprintf(fd, "]");
  }
  fprintf(fd, ");\n");
  fprintf(fd, "expand(determinant(MM));\n");
  return;
}

// Print a submatrix as an (undirected) graph in the "dot" language

void printdotindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2, long lab){
  int i,j;
  fprintf(fd, "graph mat%ld {\n", lab);
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0){}
      else if(imat[vec1[i]][vec2[j]]<0)
	fprintf(fd, "\tR%d -- C%d [label=v%d] [color=red]\n",vec1[i],vec2[j],-imat[vec1[i]][vec2[j]]);
      else
	fprintf(fd, "\tR%d -- C%d [label=v%d]\n",vec1[i],vec2[j],imat[vec1[i]][vec2[j]]);
    }
  }
  fprintf(fd, "}\n\n");
  return;
}

// print a pair of submatrices as a directed graph in the dot language

void printdotpairsubmat(FILE *fd, int **imat, int **imat1, int *vec1, int *vec2, int k1, int k2, long lab){
  int i,j;
  fprintf(fd, "digraph mat%ld {\n", lab);
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0 && imat1[vec1[i]][vec2[j]]==0){}
      else if(imat[vec1[i]][vec2[j]]<0 && imat1[vec1[i]][vec2[j]]<0)
	fprintf(fd, "\tR%d -> C%d [label=v%d] [color=red] [dir=none]\n",vec1[i],vec2[j],-imat1[vec1[i]][vec2[j]]);
      else if(imat[vec1[i]][vec2[j]]>0 && imat1[vec1[i]][vec2[j]]>0)
	fprintf(fd, "\tR%d -> C%d [label=v%d] [dir=none]\n",vec1[i],vec2[j],imat1[vec1[i]][vec2[j]]);
      else{
	if(imat[vec1[i]][vec2[j]]<0)
	  fprintf(fd, "\tC%d -> R%d [label=v%d] [color=red]\n",vec1[i],vec2[j],-imat1[vec1[i]][vec2[j]]);
	if(imat1[vec1[i]][vec2[j]]<0)
	  fprintf(fd, "\tR%d -> C%d [label=v%d] [color=red]\n",vec1[i],vec2[j],-imat1[vec1[i]][vec2[j]]);
	if(imat[vec1[i]][vec2[j]]>0)
	  fprintf(fd, "\tC%d -> R%d [label=v%d]\n",vec1[i],vec2[j],imat1[vec1[i]][vec2[j]]);
	if(imat1[vec1[i]][vec2[j]]>0)
	  fprintf(fd, "\tR%d -> C%d [label=v%d]\n",vec1[i],vec2[j],imat1[vec1[i]][vec2[j]]);
      }
    }
  }
  fprintf(fd, "}\n\n");
  return;
}

void printmaximaindmat(int **imat, int k1, int k2){
  int i,j;
  fprintf(stderr, "MM:matrix(");
  for(i=0;i<k1;i++){
    fprintf(stderr, "[");
    for(j=0;j<k2;j++){
      if(imat[i][j]==0)
	fprintf(stderr, "0");
      else if(imat[i][j]<0)
	fprintf(stderr, "-v%d",-imat[i][j]);
      else
	fprintf(stderr, "v%d",imat[i][j]);
      if(j<(k2-1))
	fprintf(stderr, ",");
    }
    if(i<k1-1)
      fprintf(stderr, "],");
    else
      fprintf(stderr, "]");
  }
  fprintf(stderr, ");\n");
  return;
}


/* read a pair of matrices from a string */
/* matrices are separated by ***** */
int readmatpairfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ***mat2){
  char *str1, *str2;
  int nl1, ml1, nl2, ml2;

  chop(str, &str1, &str2, "*****");

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);
    return 0;
  }
  (*mat2)=readmatrixfromstr(str2, &nl2, &ml2);
  if(!(*mat2) || (nl1!=nl2) || (ml1!=ml2)){
    fprintf(stderr, "ERROR: matrices could not be read.\n");
    free_imatrix((*mat1), 0, nl1, 0, ml1);
    (*mat1)=NULL;
    if((*mat2)){
      free_imatrix((*mat2), 0, nl2, 0, ml2);
      (*mat2)=NULL;
    }

    free(str1);free(str2);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  free(str1);free(str2);
  return 1;

}

/* read a matrix and its rank from a string */
/* matrices are separated by ***** */

int readmatrankfromstr(char *str, int *nlen, int *mlen, int ***mat1, int *rk){
  char *str1, *str2;
  int nl1, ml1, nl2=0, ml2=0;
  int **mat2;

  chop(str, &str1, &str2, "*****");

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);
    fprintf(stderr, "ERROR: matrix could not be read.\n");
    return 0;
  }
  mat2=readmatrixfromstr(str2, &nl2, &ml2);
  if(!mat2 || nl2!=1 || ml2!=1){
    fprintf(stderr, "ERROR: rank could not be read. (It must be an integer)\n");
    free_imatrix((*mat1), 0, nl1-1, 0, ml1-1);
    (*mat1)=NULL;
    if(mat2){
      free_imatrix(mat2, 0, nl2-1, 0, ml2-1);
      mat2=NULL;
    }
    free(str1);free(str2);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  (*rk) = mat2[0][0];
  /*  fprintf(stderr, "rk = %d\n", *rk); */
  free_imatrix(mat2, 0, nl2-1, 0, ml2-1);
  free(str1);free(str2);
  return 1;

}

int readmatpairrankfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ***mat2, int *rk){
  char *str1, *strtmp, *str2, *str3;
  int nl1, ml1, nl2, ml2, nl3, ml3;
  int **mat3;

  chop(str, &str1, &strtmp, "*****");
  fprintf(stderr, "str1 = %s\n", str1);
  fprintf(stderr, "strtmp = %s\n", strtmp);
  chop(strtmp, &str2, &str3, "*****");
  fprintf(stderr, "str2 = %s\n", str2);
  fprintf(stderr, "str3 = %s\n", str3);

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*mat2)=readmatrixfromstr(str2, &nl2, &ml2);
  if(!(*mat2) || (nl1!=nl2) || (ml1!=ml2)){
    fprintf(stderr, "ERROR: matrices could not be read.\n");
    free_imatrix((*mat1), 0, nl1, 0, ml1);
    (*mat1)=NULL;
    if((*mat2)){
      free_imatrix((*mat2), 0, nl2, 0, ml2);
      (*mat2)=NULL;
    }

    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  mat3=readmatrixfromstr(str3, &nl3, &ml3);
  if(!mat3 || nl3!=1 || ml3!=1){
    //    fprintf(stderr, "ttnl3 = %d, ml3 = %d\n", nl3, ml3);
    fprintf(stderr, "ERROR: rank could not be read. (It must be an integer)\n");
    free_imatrix((*mat1), 0, nl1-1, 0, ml1-1);
    free_imatrix((*mat2), 0, nl2-1, 0, ml2-1);
    (*mat1)=NULL;(*mat2)=NULL;
    if(mat3){
      free_imatrix(mat3, 0, nl3-1, 0, ml3-1);
      mat3=NULL;
    }
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  (*rk) = mat3[0][0];
  /*  fprintf(stderr, "rk = %d\n", *rk); */
  free_imatrix(mat3, 0, nl3-1, 0, ml3-1);


  (*nlen)=nl1;(*mlen)=ml1;
  free(str1);free(str2);free(strtmp);free(str3);
  return 1;

}

int ***readiimatrixfromstr(char *str, int *nlen, int *mlen){
  int numlines=0, numtmp, i,j,nt;
  char *line;
  long pos=0;
  int ***imat;
  char *wd,*wd0;
  line = getlinefromstr(&pos, str);
  (*mlen)=0;(*nlen)=0;
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;
      numtmp=getnumsegs(line);
      if(numlines!=1 && numtmp!=(*mlen)){fprintf(stderr, "WARNING: matrix does not seem well defined.\n");}
      if(numtmp>(*mlen)){(*mlen)=numtmp;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }

  free(line);
  (*nlen)=numlines;

  // no matrix found

  if((*nlen)==0 || (*mlen)==0){
    fprintf(stderr, "No matrix found in string %s.\n", str);
    return NULL;
  }

  imat=(int***) malloc(sizeof(int**)*(*nlen));
  for(i=0;i<(*nlen);i++)
    imat[i]=(int**) malloc(sizeof(int*)*(*mlen));

  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;
      for(i=1;i<(*mlen)+1;i++){
	wd0=getnthwd(line, i);
	nt=getnumints1(wd0);
	imat[numlines-1][i-1]=(int*) malloc(sizeof(int)*(nt+1));
	imat[numlines-1][i-1][0]=nt;
	for(j=1;j<=nt;j++){
	  wd=getnthint1(wd0, j);
	  imat[numlines-1][i-1][j]=atoi(wd);
	  free(wd);
	}
	free(wd0);
      }
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return imat;

}

void freeiim(int ***mat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      free((char *)mat[i][j]);
    free((char *)mat[i]);
  }
  free((char *)mat);
  return;
}


int readiimatpairrankfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ****mat2, int *rk){
  char *str1, *strtmp, *str2, *str3;
  int nl1, ml1, nl2, ml2, nl3, ml3;
  int **mat3;

  chop(str, &str1, &strtmp, "*****");
  fprintf(stderr, "str1 = %s\n", str1);
  fprintf(stderr, "strtmp = %s\n", strtmp);
  chop(strtmp, &str2, &str3, "*****");
  fprintf(stderr, "str2 = %s\n", str2);
  fprintf(stderr, "str3 = %s\n", str3);

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*mat2)=readiimatrixfromstr(str2, &nl2, &ml2);
  if(!(*mat2) || (nl1!=nl2) || (ml1!=ml2)){
    fprintf(stderr, "ERROR: matrices could not be read.\n");
    free_imatrix((*mat1), 0, nl1, 0, ml1);
    (*mat1)=NULL;
    if((*mat2)){
      freeiim((*mat2), nl2, ml2);
      (*mat2)=NULL;
    }

    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  mat3=readmatrixfromstr(str3, &nl3, &ml3);
  if(!mat3 || nl3!=1 || ml3!=1){
    //    fprintf(stderr, "ttnl3 = %d, ml3 = %d\n", nl3, ml3);
    fprintf(stderr, "ERROR: rank could not be read. (It must be an integer)\n");
    free_imatrix((*mat1), 0, nl1-1, 0, ml1-1);
    freeiim((*mat2), nl2, ml2);
    (*mat1)=NULL;(*mat2)=NULL;
    if(mat3){
      free_imatrix(mat3, 0, nl3-1, 0, ml3-1);
      mat3=NULL;
    }
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  (*rk) = mat3[0][0];
  /*  fprintf(stderr, "rk = %d\n", *rk); */
  free_imatrix(mat3, 0, nl3-1, 0, ml3-1);


  (*nlen)=nl1;(*mlen)=ml1;
  free(str1);free(str2);free(strtmp);free(str3);
  return 1;

}



int **readmatrixfromstr(char *str, int *nlen, int *mlen){
  int numlines=0, numtmp, i;
  char *line;
  long pos=0;
  int **imat;
  char *wd;
  line = getlinefromstr(&pos, str);
  (*mlen)=0;(*nlen)=0;
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;

      numtmp=getnumints(line);
      if(numlines!=1 && numtmp!=(*mlen)){fprintf(stderr, "WARNING: matrix does not seem well defined.\n");}
      if(numtmp>(*mlen)){(*mlen)=numtmp;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }

  free(line);
  (*nlen)=numlines;

  // no matrix found

  if((*nlen)==0 || (*mlen)==0){
    fprintf(stderr, "No matrix found in string %s.\n", str);
    return NULL;
  }

  // otherwise create the matrix

  imat=imatrix(0, numlines-1, 0, (*mlen)-1);
  
  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;
      for(i=1;i<(*mlen)+1;i++){
	wd=getnthint(line, i);
	imat[numlines-1][i-1]=atoi(wd);
	free(wd);
      }
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return imat;

}

int isreac(char *str){
  if(strstr(str, "<-->") || strstr(str, "<==>") || strstr(str, "<->")|| strstr(str, "<>") || strstr(str, "<=>")|| strstr(str, "-->")|| strstr(str, "==>")|| strstr(str, "->")|| strstr(str, ">")|| strstr(str, "=>")|| strstr(str, "<--")|| strstr(str, "<==")|| strstr(str, "<-")|| strstr(str, "<")|| strstr(str, "<="))
    return 1;
  return 0;
}

// S is the stoichiometric matrix
// Vpat is the pattern matrix of V
// Si is the irreversible stoichiometric matrix
// Sil is *minus* the corresponding matrix of powers (assuming mass action)
// stoichl is the left stoichiometric matrix
// stoichr is the right stoichiometric matrix

int getallreacs(char *str, int ***S, int ***Vpat, int ***Si, int ***Sil, int ***stoichl, int ***stoichr, char ***chems, bool *haszerocomplex, int *n, int *m, int *cols3, int *allrev, int *allgood, int debug){

  long pos=0;
  int pos1;
  char *line;
  int totchems=0;
  int i, ind;
  int irrcols=0, numcols=0;
  int lstoic, rstoic;

  char **leftchems, **rightchems;
  int *leftstoics,*rightstoics;
  int numlines=0, numleft, numright, rev=0;
  (*allrev)=1;
  (*allgood)=1;
  (*chems)=NULL;

  if(debug){cerr << str << endl;}
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      if(!getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev)){
	freearraydat(leftchems, numleft);
	freearraydat(rightchems, numright);
	freearraydat((*chems), totchems);
	free((FREE_ARG) (leftstoics));
	free((FREE_ARG) (rightstoics));
	return 0;
      }
      if(rev==0){
	irrcols++;(*allrev)=0;
      }
      else
	irrcols+=2;

      for(i=0;i<numleft;i++)
	totchems= addv1(totchems, leftchems[i], chems);
      for(i=0;i<numright;i++)
	totchems= addv1(totchems, rightchems[i], chems);

      freearraydat(leftchems, numleft);
      freearraydat(rightchems, numright);
      free((FREE_ARG) (leftstoics));
      free((FREE_ARG) (rightstoics));
      numlines++;
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  (*m)=numlines;
  (*n)=totchems;
  (*cols3)=irrcols;

  (*S)=imatrix(0, (*n)-1, 0, (*m)-1);
  (*Vpat)=imatrix(0, (*n)-1, 0, (*m)-1);
  (*Si)=imatrix(0, (*n)-1, 0, irrcols-1);
  (*Sil)=imatrix(0, (*n)-1, 0, irrcols-1);
  (*stoichl)=imatrix(0, (*n)-1, 0, irrcols-1);
  (*stoichr)=imatrix(0, (*n)-1, 0, irrcols-1);


  //pass two

  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev);
      if(numleft==0||numright==0)
	(*haszerocomplex)=1;

      for(i=0;i<(*n);i++){//each substrate
	lstoic=0;pos1=0;
	while(pos1<numleft){
	  if((ind=isinarray(leftchems+pos1, numleft-pos1, (*chems)[i]))>=0){
	    pos1+=ind;
	    lstoic+=leftstoics[pos1];
	    pos1++;
	  }
	  else
	    pos1=numleft;
	}
	rstoic=0;pos1=0;
	while(pos1<numright){
	  if((ind=isinarray(rightchems+pos1, numright-pos1, (*chems)[i]))>=0){
	    pos1+=ind;
	    rstoic+=rightstoics[pos1];
	    pos1++;
	    //	    fprintf(stderr, "pos1 = %d, chem = %s, stoic = %d\n", ind, rightchems[ind], rstoic);
	  }
	  else
	    pos1=numright;
	}

	(*stoichl)[i][numcols]=lstoic;
	(*stoichr)[i][numcols]=rstoic;
	if(rev){
	  (*stoichl)[i][numcols+1]=rstoic;
	  (*stoichr)[i][numcols+1]=lstoic;
	}

	if(lstoic>0 && rstoic > 0){ // on both sides of reac
	  (*allgood)=0;
	  (*S)[i][numlines]=rstoic-lstoic;
	  if(rev==1){
	    (*Vpat)[i][numlines]=2;
	    (*Si)[i][numcols]=rstoic-lstoic;
	    (*Si)[i][numcols+1]=-rstoic+lstoic;
	    (*Sil)[i][numcols]=-lstoic;
	    (*Sil)[i][numcols+1]=-rstoic;

	    //	    (*Sil)[i][numcols+1]=-rightstoics[ind];
	  }
	  else{
	    (*Vpat)[i][numlines]=-1;
	    (*Si)[i][numcols]=rstoic-lstoic;
	    (*Sil)[i][numcols]=-lstoic;
	  }

	}
	else if(lstoic>0){ // only on left
	  (*S)[i][numlines]=-lstoic;
	  (*Vpat)[i][numlines]=-1;
	  if(rev==1){
	    (*Si)[i][numcols]=-lstoic;
	    (*Sil)[i][numcols]=-lstoic;
	    (*Si)[i][numcols+1]=lstoic;
	    (*Sil)[i][numcols+1]=0;
	  }
	  else{
	    (*Si)[i][numcols]=-lstoic;
	    (*Sil)[i][numcols]=-lstoic;
	  }
	}
	else if(rstoic>0){ // only on right
	  (*S)[i][numlines]=rstoic;
	  if(rev==1){
	    (*Vpat)[i][numlines]=1;
	    (*Si)[i][numcols]=rstoic;
	    (*Sil)[i][numcols]=0;
	    (*Si)[i][numcols+1]=-rstoic;
	    (*Sil)[i][numcols+1]=-rstoic;
	  }
	  else{
	    (*Vpat)[i][numlines]=0;
	    (*Si)[i][numcols]=rstoic;
	    (*Sil)[i][numcols]=0;
	  }
	}

	else{ // doesn't figure in the reaction
	  (*S)[i][numlines]=0;
	  (*Vpat)[i][numlines]=0;

	  (*Si)[i][numcols]=0;
	  (*Sil)[i][numcols]=0;
	  if(rev==1){
	    (*Si)[i][numcols+1]=0;
	    (*Sil)[i][numcols+1]=0;
	  }

	}

      }

      freearraydat(leftchems, numleft);
      freearraydat(rightchems, numright);
      free((FREE_ARG) (leftstoics));
      free((FREE_ARG) (rightstoics));
      numlines++;if(rev==1){numcols+=2;}else{numcols++;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return 1;

}

//The light version: just return irreversible stoichiometric matrix and the 
//left and right stoichiometric matrices; species names; total species and reactions
//Note that Sil is the left stoichiometric matrix, and not its minus.
int getallreacs(char *str, int ***Si, int ***Sil, int ***Sir, char ***chems, int *n, int *m){
  long pos=0;
  int pos1;
  char *line;
  int totchems=0;
  int i, ind;
  int irrcols=0, numcols=0;
  int lstoic, rstoic;

  char **leftchems, **rightchems;
  int *leftstoics,*rightstoics;
  int numlines=0, numleft, numright, rev=0;
  (*chems)=NULL;

  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      if(!getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev)){
	freearraydat(leftchems, numleft);
	freearraydat(rightchems, numright);
	freearraydat((*chems), totchems);
	free((FREE_ARG) (leftstoics));
	free((FREE_ARG) (rightstoics));
	return 0;
      }
      if(rev==0)
	irrcols++;
      else
	irrcols+=2;

      for(i=0;i<numleft;i++)
	totchems= addv1(totchems, leftchems[i], chems);
      for(i=0;i<numright;i++)
	totchems= addv1(totchems, rightchems[i], chems);

      freearraydat(leftchems, numleft);
      freearraydat(rightchems, numright);
      free((FREE_ARG) (leftstoics));
      free((FREE_ARG) (rightstoics));
      numlines++;
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  (*m)=irrcols;
  //fprintf(stderr, "found %d reactions.\n", *m);fflush(stderr);
  (*n)=totchems;

  (*Si)=imatrix(0, (*n)-1, 0, (*m)-1);
  (*Sil)=imatrix(0, (*n)-1, 0, (*m)-1);
  (*Sir)=imatrix(0, (*n)-1, 0, (*m)-1);


  //pass two

  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev);

      for(i=0;i<(*n);i++){//each substrate
	lstoic=0;pos1=0;
	while(pos1<numleft){
	  if((ind=isinarray(leftchems+pos1, numleft-pos1, (*chems)[i]))>=0){
	    pos1+=ind;
	    lstoic+=leftstoics[pos1];
	    pos1++;
	  }
	  else
	    pos1=numleft;
	}
	rstoic=0;pos1=0;
	while(pos1<numright){
	  if((ind=isinarray(rightchems+pos1, numright-pos1, (*chems)[i]))>=0){
	    pos1+=ind;
	    rstoic+=rightstoics[pos1];
	    pos1++;
	    //	    fprintf(stderr, "pos1 = %d, chem = %s, stoic = %d\n", ind, rightchems[ind], rstoic);
	  }
	  else
	    pos1=numright;
	}

	(*Si)[i][numcols]=rstoic-lstoic;
	(*Sil)[i][numcols]=lstoic;
	(*Sir)[i][numcols]=rstoic;
	if(rev==1){
	  (*Si)[i][numcols+1]=-rstoic+lstoic;
	  (*Sil)[i][numcols+1]=rstoic;
	  (*Sir)[i][numcols+1]=lstoic;
	}

      }

      freearraydat(leftchems, numleft);
      freearraydat(rightchems, numright);
      free((FREE_ARG) (leftstoics));
      free((FREE_ARG) (rightstoics));
      numlines++;if(rev==1){numcols+=2;}else{numcols++;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return 1;

}

//from human readable to di6
char *reacstrtodi6(char *str, int minlayers, int *layers, int *n, int *m){
  int **S, **Sl, **Sr;
  char **chems;
  int r1,r2,clen,Vlen;
  bool *V;
  char *out;
  getallreacs(str, &S, &Sl, &Sr, &chems, n, m);
  V=CRNPN3(Sl,Sr,(*n),(*m),minlayers,layers);
  if((*n)+(*m)>62/(*layers)){
    fprintf(stderr, "ERROR in reacstrtodi6: Total CRN size too large. Exiting.\n");
    exit(0);
  }

  r1=(*n)+(*m);r2=(*layers)*(*layers)*r1*r1;
  Vlen=(r2%6==0)?r2:r2+6-r2%6;
  clen=Vlen/6;

  out=(char *)malloc((size_t) ((clen+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+(*layers)*r1;
  out[2+clen]=0;

  amtodig(V,clen,out);
  //fprintf(stderr, "out = %s\n", out);

  free((char *)V);
  free_imatrix(S, 0, (*n)-1, 0, (*m)-1);
  free_imatrix(Sl, 0, (*n)-1, 0, (*m)-1);
  free_imatrix(Sr, 0, (*n)-1, 0, (*m)-1);
  freearraydat(chems, (*n));
  return out;

}

int genMAXMAreacs(char *fname, int **imat1, int **imat2, int n, int m){
  int i, j;
  FILE *fd;

 fd = fopen(fname, "w");
  if(!fd){
    fprintf(stderr, "ERROR in getMAXMAreacs: \"%s\" could not be opened for writing.\n", fname);
    return 0;
  }

  for (j=0;j<m;j++){
    fprintf(fd, "f%d:k%d", j, j);
    for(i=0;i<n;i++){
      if(imat2[i][j]!=0){
	if(imat2[i][j]==-1){
	  fprintf(fd, "*S%d", i);
	}
	else{
	  fprintf(fd, "*(S%d**%d)", i, -imat2[i][j]);
	}
      }
    }
    fprintf(fd, ";\n");
  }
    fprintf(fd, "\n\n");

  for(i=0;i<n;i++){
    fprintf(fd, "F%d: ", i);
    for(j=0;j<m;j++){
      if(imat1[i][j]!=0){
	if(imat1[i][j]==1){
	  fprintf(fd, "+f%d", j);
     
	}
	else if(imat1[i][j]==-1){
	  fprintf(fd, "-f%d", j);
	}
	else if (imat1[i][j]<0){
	  fprintf(fd, "-%d*f%d", -imat1[i][j], j);

	}
	else{
	  fprintf(fd, "+%d*f%d", imat1[i][j], j);

	}

      }
    }
    fprintf(fd, ";\n");

  }
    fprintf(fd, "\n\n");

  fprintf(fd, "jj:matrix(");
  for(i=0;i<n;i++){
    fprintf(fd, "[");
    for(j=0;j<n;j++){
      fprintf(fd, "diff(F%d, S%d)", i, j);
      if(j<n-1){
	fprintf(fd, ", ");
      }
    }
    fprintf(fd, "]");
    if(i<n-1){
      fprintf(fd, ", ");
    }

  }
  fprintf(fd, ");\n\n");

  fclose(fd);
  return 1;

}


//The number in an expression such as 2A0. 1 if empty
int split1(char *s){
  char *p;
  int j=0, ret=1;
  // get the numeric part
  p=strdup(s);
  while(s[j] && isdigit(s[j])){j++;}
  p[j]=0;
  if(!isonlyspace(p))
    ret=atoi(p);
  free(p);
  return ret;
}

// get the character part in an expressino such as 2A0. No error checking
char *split2(char *s){
  char *p;
  int j=0;
  while(s[j] && isdigit(s[j])){j++;}
  p=strdup(s+j);
  return p;
}


int getreac(char *str, char ***leftchems, int **leftstoics, int *numleft, char ***rightchems, int **rightstoics, int *numright, int *rev){
  char *p, *left, *right;
  int j,sp;
  char **v1;
  char **tmp;
  int num1;
  int flag=1;
  if((p=strstr(str, "<-->")) || (p=strstr(str, "<==>")) || (p=strstr(str, "<->"))|| (p=strstr(str, "<=>")) || (p=strstr(str, "<>"))){
    if(strstr(p,"<-->") || strstr(p,"<==>"))
      sp=4;
    else if(strstr(p,"<->") || strstr(p,"<=>"))
      sp=3;
    else
      sp=2;
    j=strlen(str)-strlen(p);
    right=strdup(str+j+sp);
    left=strdup(str);
    left[j]=0;
    (*rev)=1;
  }
  else if((p=strstr(str, "-->")) || (p=strstr(str, "==>")) || (p=strstr(str, "->")) || (p=strstr(str, "=>")) || (p=strstr(str, ">"))){
    if(strstr(p,"-->") || strstr(p,"==>"))
      sp=3;
    else if(strstr(p,"->") || strstr(p,"=>"))
      sp=2;
    else
      sp=1;

    j=strlen(str)-strlen(p);
    right=strdup(str+j+sp);
    left=strdup(str);
    left[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, "<--")) || (p=strstr(str, "<==")) || (p=strstr(str, "<-")) || (p=strstr(str, "<=")) || (p=strstr(str, "<"))){
    if(strstr(p,"<--") || strstr(p,"<=="))
      sp=3;
    else if(strstr(p,"<-") || strstr(p,"<="))
      sp=2;
    else
      sp=1;

    j=strlen(str)-strlen(p);
    left=strdup(str+j+sp);
    right=strdup(str);
    right[j]=0;
    (*rev)=0;
  }
  else{
    fprintf(stderr, "ERROR in reaction %s\n", str);
    flag=0;
    left=strdup("");right=strdup("");
  }

  //left of reaction
  (*numleft)=chemgts2a(left, &v1, '+');
  (*leftstoics)=(int*) malloc(sizeof(int) * (*numleft));
  (*leftchems)=(char**) malloc(sizeof(char*) * (*numleft));
  for(j=0;j<(*numleft);j++){//allow a space between stoich and species name
    num1=chemgts2a(v1[j], &tmp, ' ');

    if(num1==1){
      (*leftstoics)[j]=split1(v1[j]);(*leftchems)[j]=split2(v1[j]);
    }
    else if(num1==2){
      if(!ispureint(tmp[0])){
	fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
      }
      (*leftstoics)[j]=atoi(tmp[0]);(*leftchems)[j]=strdup(tmp[1]);
    }
    else{
      (*leftstoics)[j]=0;(*leftchems)[j]=strdup("");
      fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
    }
    freearraydat(tmp, num1);
  }

  freearraydat(v1, (*numleft));

  //right of reaction

  (*numright)=chemgts2a(right, &v1, '+');
  (*rightstoics)=(int*) malloc(sizeof(int) * (*numright));
  (*rightchems)=(char**) malloc(sizeof(char*) * (*numright));
  for(j=0;j<(*numright);j++){
    num1=chemgts2a(v1[j], &tmp, ' ');  
  
    if(num1==1){
      (*rightstoics)[j]=split1(v1[j]);(*rightchems)[j]=split2(v1[j]);
      //fprintf(stderr, "%s: rstoic = %d, chem = %s\n", v1[j], (*rightstoics)[j], (*rightchems)[j]);
    }
    else if(num1==2){
     if(!ispureint(tmp[0])){
	fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
      }
      (*rightstoics)[j]=atoi(tmp[0]);(*rightchems)[j]=strdup(tmp[1]);
    }
    else{
      (*rightstoics)[j]=0;(*rightchems)[j]=strdup("");
      fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
    }
    freearraydat(tmp, num1);
  }

  freearraydat(v1, (*numright));


  free(left);
  free(right);
  return flag;

}



// For splitting a complex, or a stoichiometry-species pair
int chemgts2a(char *s, char ***v, char sep){
  // sep is the separator
  int i, j, k;
  int numgets=0;
  char *tmp, *tmp1;
  i=0, k=0;

  (*v)=NULL;
  while(s[k]){//first pass just to count
    j=0;
    while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
    while((s[k] != sep) && !isend(s[k])){j++;k++;}
    if(j>0){
      tmp1=strchop2(s, k-j, j);
      tmp=lrtrim(tmp1);
      if(!ispureint(tmp) || atoi(tmp)!=0)
	numgets++;
      free(tmp);free(tmp1);
    }
  }
  if(numgets>0){// found blocks
    (*v)=(char**) malloc(sizeof(char*) * numgets);
    i=0;k=0;
    while(s[k]){
      j=0;
      while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
      while((s[k] != sep) && !isend(s[k])){j++;k++;}
      if(j>0){
	tmp1=strchop2(s, k-j, j);
	tmp=lrtrim(tmp1);
	if(!ispureint(tmp) || atoi(tmp)!=0)
	  (*v)[i++] = lrtrim(strchop2(s, k-j, j));
	free(tmp);free(tmp1);
      }
    }
  }
  return numgets;
}



/* check if the matrix stored in file fname */
/* is sign nonsingular or sign singular */

int strmatisSNSSS(char *fname){

  int **imat1;
  char *str;
  int mlen=0, nlen=0;
  int retval;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
    fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }

  imat1=readmatrixfromstr(str, &nlen, &mlen);
  if(!(*imat1)){
    fprintf(stderr, "ERROR: Couldn't read a matrix from the data in file \"%s\". EXITING. \n", fname);
    free(str);
    return -1;
  }
  fprintf(stderr, "matrix found: \n");
  printmat(imat1, nlen, mlen);

  if(nlen!=mlen){
    fprintf(stderr, "The matrix is not square. Exiting.\n\n");
    exit(-1);
  }

  if((retval=matrixisSNSSS(imat1, nlen))==1){
    fprintf(stderr, "The matrix is sign nonsingular with positive determinant\n");
  }
  else if(retval==-1)
    fprintf(stderr, "The matrix is sign nonsingular with negative determinant\n");
  else if(retval==2)
    fprintf(stderr, "The matrix is sign singular\n");
  else{
    fprintf(stderr, "The matrix is neither sign nonsingular nor sign singular\n");
  }
  free(str);
  free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  return 0;

}

/* check if the matrix stored in file fname */
/* is SSD */


int strmatisSSD(char *fname, int q){

  int **imat1;
  char *str;
  int mlen=0, nlen=0;
  int retval;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
    fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }

  imat1=readmatrixfromstr(str, &nlen, &mlen);
  if(!(*imat1)){
    fprintf(stderr, "ERROR: Couldn't read a matrix from the data in file \"%s\". EXITING. \n", fname);
    free(str);
    return -1;
  }
  fprintf(stderr, "matrix: \n");
  printmat(imat1, nlen, mlen);


  if((retval=isSSD(imat1, nlen, mlen, q))){
    fprintf(stderr, "The matrix is SSD.\n");
  }
  else{
    fprintf(stderr, "The matrix is not SSD\n");
  }
  free(str);
  free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  return 0;

}




/* check if the matrix stored in file fname */
/* is CSD */


int strmatisCSD(char *fname, int q){

  int **imat1;
  char *str;
  int mlen=0, nlen=0;
  int retval;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
    fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }

  imat1=readmatrixfromstr(str, &nlen, &mlen);
  if(!(*imat1)){
    fprintf(stderr, "ERROR: Couldn't read a matrix from the data in file \"%s\". EXITING. \n", fname);
    free(str);
    return -1;
  }
  fprintf(stderr, "matrix: \n");
  printmat(imat1, nlen, mlen);


  if((retval=isCSD(imat1, nlen, mlen, q))){
    fprintf(stderr, "The matrix is CSD.\n");
  }
  else{
    fprintf(stderr, "The matrix is not CSD\n");
  }
  free(str);
  free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  return 0;

}

// For each siphon:
// does there exist a nonnegative vector orthogonal to the siphon face 
// and in ker(Gamma^t)?
// In the language of Angeli et. al., Petri nets paper, we check
// whether each siphon contains the support of a P-semiflow.
// Uses linear programming

bool structpersist(int **Si, int nlen, int mlen, int **allminsiphons, int totminsiphons){
  int k,j,l,m,siphlen;
  int *ia, *ja;
  double *ar,z;
  int r;
  int **siphmat;
  bool goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(nlen+1)*(nlen+mlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(nlen+1)*(nlen+mlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(nlen+1)*(nlen+mlen+1))*sizeof(double)));
  for(k=0;k<totminsiphons;k++){

    siphlen=nlen-allminsiphons[k][0];
    //fprintf(stderr, "siphlen=%d\n", siphlen);
    //make a siphon matrix
    siphmat=imatrix(0,nlen-1,0,siphlen-1);
    inittozero(siphmat,nlen,siphlen);//initialise

    r=0;
    for(j=0;j<nlen;j++){
      if(!isinlist(j,allminsiphons[k]+1,allminsiphons[k][0])){
	siphmat[j][r++]=1;
      }
    }
    /* printmat(siphmat, nlen, siphlen); */

    lp= glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    //number of constraints (rows) = col dimension of Gamma + boundedness constraint + siphlen
    glp_add_rows(lp, 1+mlen+siphlen);
    glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 10.0);
    for(j=2;j<mlen+2;j++)
      glp_set_row_bnds(lp, j, GLP_FX, 0.0, 0.0);
    for(j=mlen+2;j<mlen+2+siphlen;j++)
      glp_set_row_bnds(lp, j, GLP_FX, 0.0, 0.0);

    //number of variables (cols) = row dimension of Gamma
    glp_add_cols(lp, nlen);
    for(j=1;j<nlen+1;j++){
      glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
      glp_set_obj_coef(lp, j, 1.0); // objective function
    }

    m=1;
    for(j=1;j<nlen+1;j++){//col indices
      ia[m]=1;ja[m]=j;ar[m]=1.0;//boundedness constraint
      m++;
      for(l=2;l<2+mlen;l++){//row indices
	ia[m]=l;ja[m]=j;ar[m]=Si[j-1][l-2];
	m++;
      }
      for(l=2+mlen;l<2+mlen+siphlen;l++){//row indices
	ia[m]=l;ja[m]=j;ar[m]=siphmat[j-1][l-2-mlen];// face constraints
	m++;
      }
    }
    //fprintf(stderr, "m=%d\n", m);
    glp_load_matrix(lp, m-1, ia, ja, ar);
    glp_simplex(lp, &parm);//can be NULL
    z = glp_get_obj_val(lp);
    /* fprintf(stderr, "z = %.2f\n", z); */
 
    glp_delete_prob(lp);
    free_imatrix(siphmat, 0,nlen-1,0,siphlen-1);
    if(z<0.01){// assume fail
      goodflag=0;
      break;
    }
  }
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));
  return goodflag;
}

//overloading to take Si and Sil as input
bool structpersist(int **Si, int **Sil, int n, int m, int debug){
  int totsiphons, totminsiphons, ret;
  int **allsiphons=NULL, **allminsiphons=NULL;
  int debugfull=(debug<=0)?0:debug-1;
  //count the siphons
  if(debug){fprintf(stderr, "\n###Entering structpersist (checking if the network has critical siphons).\n");}

  totsiphons=checksiphons(Si,Sil,n,m,&allsiphons,&totminsiphons,&allminsiphons);
  if(debugfull){
    if(totsiphons){
      fprintf(stderr, "All siphons (minimal ones have *):\n");
      printsiphons(allsiphons, totsiphons, allminsiphons, totminsiphons);
    }
    else
      fprintf(stderr, "The network has no siphons.\n");
  }
  ret=structpersist(Si, n, m, allminsiphons, totminsiphons);
  free_imat(allsiphons,totsiphons);
  free_imat(allminsiphons,totminsiphons);
  if(debug)
    fprintf(stderr, "%s\n",ret?"The network is structurally persistent (no critical siphons)":"The network is not structurally persistent (it has critical siphons)");
  return ret;
}


//overloading: PN AM input
bool structpersist(int **AM, int n, int m, int debug){
  int **S, **Sl;
  int ret,minus=1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  ret=structpersist(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading di6 input
bool structpersist(char *di6, int n, int m, int debug){
  int **S, **Sl;
  int ret,minus=0;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  ret=structpersist(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}



// is the integer vector s a member of the array v? If so return its index, 
// if not return -1
// assume first entry is length of vector
int isinarray4(int **v, int numv, int *s){
  int i;
  for(i=0;i<numv;i++){
    if(unordlistareeq(v[i]+1,v[i][0],s+1,s[0]))
      return i;
  }
  return -1;
}

// Confirm that each member of A is also in B
// E.g. the CCs are also SCCs
bool is_sublist(int **A, int numA, int **B, int numB){
  int i;
  for(i=0;i<numA;i++){
    if(isinarray4(B,numB,A[i])<0)
      return 0;
  }
  return 1;
}

// Which member of the partition A is integer i a member of?
// Return -1 if it isn't in the partition
// Assume that the first entry of each integer string in A
// is the length of the remainder

int partmember(int i, int **A, int totA){
  int j,k;
  for(j=0;j<totA;j++){
    for(k=1;k<A[j][0]+1;k++){
      if(i==A[j][k])
	return j;
    }
  }
  return -1;
}

// return the index of the zero complex in a list of complexes, 
// or -1 if there is none
int zerocmplx(int **cmplxs, int totcmplx, int n){
  int i,j;
  bool flg;
  for(i=0;i<totcmplx;i++){
    flg=1;
    for(j=0;j<n;j++){
      if(cmplxs[i][j]){// not zero complex
	flg=0;break;
      }
    }
    if(flg)
      return i;
  }
  return -1;
}


//is b a sublist of a. (First entry is vector length)
bool unordsublist(int *a, int *b){
  int i;
  if(b[0]>a[0])
    return 0;
  for(i=1;i<b[0]+1;i++){
    if(!isinlist(b[i],a+1,a[0]))
       return 0;
  }
  return 1;
}

// is a given SCC b in a given CC a? - just check first member as
// each SCC is in a CC iff any element of it is
bool SCCinCC(int *SCCj,int *CCj){
  if(SCCj[0]>CCj[0] || CCj[0]<1 || SCCj[0]<1)
    return 0;
  if(!isinlist(SCCj[1],CCj+1,CCj[0]))
    return 0;
  return 1;
}

//Is an SCC terminal in an CC
//Assume we've checked the SCC is in the CC already

int isterm(int *SCCj,int *CCj,int **cmpmat,int totcmplx){
  int i,j;
  if(SCCj[0]==CCj[0]) // same size: sole terminal SCC in CC
    return 2;
  for(i=1;i<SCCj[0]+1;i++){//each complex in SCC
    for(j=0;j<totcmplx;j++){//outedges from SCCj[i]
      if(cmpmat[SCCj[i]][j] && !isinlist(j,SCCj+1,SCCj[0]))//outedge elsewhere
	return 0;//...so not terminal
    }
  }
  return 1;
}

// Check if the system is weakly reversible
// We do this by computing the complex graph, 
// computing all its SCCs and all its CCs
// and checking that each CC is an SCC
// It could be more efficient to check each
// CC is an SCC as it is computed. 

bool weak_rev(int **imatir, int Srank, int **stoichl, int **stoichr, int n, int m, int *numcomp, int *numlink, bool haszero, bool *zeronotterm, bool *zeroinitial, bool *def1flg, int *deficiency, char **chems, bool statswitch, bool htmlswitch, int debug){
  //stoichl is the left stoichiometric matrix; stoichr is the right stoichiometric matrix
  int i,j,jj,k,m1,totcmplx=0,zeroint;
  int **stoichlt, **stoichrt; //for transposed matrices
  int **cmplxs=NULL;
  int ind1,ind2;
  int **cmpmat; //for SCCs (the digraph)
  int **cmpmat1; //for CCs (the graph)
  int edg[m][2];
  int xc[n];//row list
  int yc[m];//col list
  int **SCC=NULL, **CC=NULL;
  int totSCC=0,totCC=0;
  int tcnt,tflg,CCrnk,totdef,CCdef;
  bool flag=0,onetflg;
  (*zeronotterm)=0;
  (*zeroinitial)=0;

  stoichlt=transposemat(stoichl, n, m);
  stoichrt=transposemat(stoichr, n, m);
  for(i=0;i<m;i++){// get the complexes and the edges of the incidence graph
    totcmplx=addintstr(totcmplx, stoichlt[i], n, &cmplxs, &ind1);
    totcmplx=addintstr(totcmplx, stoichrt[i], n, &cmplxs, &ind2);
    edg[i][0]=ind1;edg[i][1]=ind2;
  }
  (*numcomp)=totcmplx;

  free_imatrix(stoichlt,0,m-1,0,n-1);
  free_imatrix(stoichrt,0,m-1,0,n-1);

  if(debug){
    fprintf(stderr, "Complexes:\n");
    for(i=0;i<totcmplx;i++){
      fprintf(stderr, "{ ");
      for(j=0;j<n;j++){
	if(cmplxs[i][j]==1)
	  fprintf(stderr, "%s ", chems[j]);
	else if(cmplxs[i][j])
	  fprintf(stderr, "%d%s ", cmplxs[i][j],chems[j]);
      }
      fprintf(stderr, "}\n");
    }
    fprintf(stderr, "\n");
  }

  //Create the adjacency matrix of the complex digraph (cmpmat) 
  //and its symmetrised version (cmpmat1)
  cmpmat=imatrix(0,totcmplx-1,0,totcmplx-1);
  cmpmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  for(i=0;i<totcmplx;i++){
    for(j=0;j<totcmplx;j++){
      cmpmat[i][j]=0;
      cmpmat1[i][j]=0;
    }
  }
  for(i=0;i<m;i++){
    cmpmat[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][1]][edg[i][0]]=1;
  }

  // Extract the SCCs and CCs using Tarjan's algorithm
  // First element of a CC or SCC is its size
  // Then the list of complexes
  SCC=Tarjan(cmpmat, totcmplx, &totSCC);
  CC=Tarjan(cmpmat1, totcmplx, &totCC);
  (*numlink)=totCC;

  if(haszero){//if the zero complex is there...
    zeroint=zerocmplx(cmplxs, totcmplx, n); //...find its index
    k=partmember(zeroint,SCC,totSCC); //it is in the kth SCC
    for(jj=0;jj<totcmplx;jj++){// for each complex...
      if(partmember(jj,SCC,totSCC)==k){// ...in the same SCC as zero
	for(j=0;j<totcmplx;j++){ // for each complex...
	  if(cmpmat[jj][j] && (partmember(j,SCC,totSCC)!=k)){// outedge kth SCC
	    (*zeronotterm)=1;break;// so not terminal
	  }
	}
      }
      if((*zeronotterm))
	break;
    }
    //check if zero is initial (so we know whether {0} is an equilibrium)
    for(j=0;j<totcmplx;j++){
      if(j!=zeroint && cmpmat[zeroint][j]){// outedge from zero
	(*zeroinitial)=1;break;//some source reactions
      }
    }

  }

  /* //Which SCCs are terminal? */
  /* term=(int *)malloc((size_t) ((totSCC)*sizeof(int))); */
  /* for(i=0;i<totSCC;i++) */
  /*   term[i]=1; */

  /* for(i=0;i<totcmplx;i++){ */
  /*   k=partmember(i,SCC,totSCC); */
  /*   for(j=0;j<totcmplx;j++){ */
  /*     if(cmpmat[i][j] && (partmember(j,SCC,totSCC)!=k)){// outedge */
  /* 	term[i]=0;break; */
  /*     } */
  /*   } */
  /* } */

  if(debug){
    fprintf(stderr, "Complex incidence matrix:\n");
    printmat(cmpmat,totcmplx,totcmplx);
  }

  /* fprintf(stderr, "terminal?\n"); */
  /* printvec(term,totSCC); */

  if(debug){
    fprintf(stderr, "CCs of complex graph:\n");
    for(i=0;i<totCC;i++){
      fprintf(stderr, "    { ");
      for(j=1;j<CC[i][0]+1;j++)
	fprintf(stderr, "%d ", CC[i][j]);
      fprintf(stderr, "}\n");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "SCCs of complex graph:\n");
    for(i=0;i<totSCC;i++){
      fprintf(stderr, "    { ");
      for(j=1;j<SCC[i][0]+1;j++)
	fprintf(stderr, "%d ", SCC[i][j]);
      fprintf(stderr, "}\n");
    }
    fprintf(stderr, "\n");
  }

  if(is_sublist(CC,totCC,SCC,totSCC)){
    if(statswitch)
      fprintf(stdout, "WRflag\n");//for batch checking
    flag=1; // weakly reversible
  }

  (*deficiency)=totcmplx-totCC-Srank;//rank of CRN
  //fprintf(stderr, "CRNdef=%d\n", (*deficiency));


  if(statswitch &&!(*deficiency)){
    fprintf(stdout, "DEF0\n");//for batch checking
  }

  if((*deficiency)){//nonzero then does CRN satisfy the def. one thm?
    //Does each CC contain exactly one terminal SCC?
    onetflg=1;//assume all good initially
    for(i=0;i<totCC;i++){//each CC
      tcnt=0; //no of terminal SCCs in CC
      for(j=0;j<totSCC;j++){//each SCC
	if(SCCinCC(SCC[j],CC[i]) && (tflg=isterm(SCC[j],CC[i],cmpmat,totcmplx))){//is in there and is terminal
	  //fprintf(stderr, "%d is terminal in %d, flg=%d\n", j, i,tflg);
	  if(tflg==2){//sole terminal CC
	    tcnt=1;
	    break;
	  }
	  else
	    tcnt++;
	  if(tcnt>=2){
	    onetflg=0;
	    break;
	  }
	}
      }
      if(onetflg==0)
	break;
    }
    //does each subnetwork have deficiency <=1?
    (*def1flg)=0;totdef=0;
    if(onetflg){//Each CC contains exactly one terminal SCC
      (*def1flg)=1;
      for(i=0;i<n;i++){xc[i]=i;}//all rows
      //compute the deficiency of each CC = numcomp-1-Srank;
      for(i=0;i<totCC;i++){//each CC
	m1=0;
	//make the reduced matrix from imatir
	for(jj=0;jj<m;jj++){//each reaction
	  if(isinlist(edg[jj][0],CC[i]+1,CC[i][0]))//only check left complex
	    yc[m1++]=jj;
	}
	CCrnk=submatrank(imatir,xc,n,yc,m1);//rank of subnetwork
	CCdef=CC[i][0]-1-CCrnk;//deficiency of subnetwork
	if(CCdef>=2){
	  (*def1flg)=0;
	  break;//
	}
	totdef+=CCdef;//total deficiency
	//fprintf(stderr, "%d subndef=%d\n", i, CCdef);
      }
    }
    //Is the sum of subnetwork deficiencies equal to the network deficiency?
    if((*def1flg) && totdef!=(*deficiency))
      (*def1flg)=0;
  }
  if(statswitch && (*def1flg)){
    fprintf(stdout, "DEF1\n");//for batch checking
  }

  //http://reaction-networks.net/wiki/Deficiency_theory#Deficiency_one_theorem

  free_imat(SCC,totSCC);
  free_imat(CC,totCC);
  free_imat(cmplxs,totcmplx);

  free_imatrix(cmpmat,0,totcmplx-1,0,totcmplx-1);
  free_imatrix(cmpmat1,0,totcmplx-1,0,totcmplx-1);

  return flag;
}

// Is the Petri Net graph of left:imat1 (nXm) and right:imat2 strongly connected?

bool DSRSC(int **imat1, int **imat2, int n, int m){
  int i,j,**SCC,**imat=imatrix(0,n+m-1,0,n+m-1);
  int totSCC=0;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++)
      imat[i][j]=0;
    for(j=n;j<n+m;j++)
      imat[i][j]=imat1[i][j-n];
  }

  for(i=n;i<n+m;i++){
    for(j=0;j<n;j++)
      imat[i][j]=imat2[j][i-n];
    for(j=n;j<n+m;j++)
      imat[i][j]=0;
  }

  /* printmat(imat,n+m,n+m); */

  SCC=Tarjan(imat, n+m, &totSCC);
  free_imatrix(imat,0,n+m-1,0,n+m-1);
  free_imat(SCC,totSCC);
  if(totSCC==1)
    return 1;
  return 0;
}

// given a matrix with no more than two nonzero entries in each column
// check if its DSR graph has o-cycles: rows can be re-signed to 
// ensure no column has more than one positive entry or more than one 
// negative entry
// The basic idea follows M. Banaji, Cycle structure in SR and DSR graphs: implications for multiple equilibria and stable oscillation in chemical reaction networks, in K. Jensen, S. Donatelli and J. Kleijn (eds): Transactions On Petri Nets and Other models of Concurrency (ToPNoC), volume V, series: LNCS, volume 6900 (2012) 
// An SR-graph with S-degree <= 2 can be R-sorted if and only if it contains no o-cycles. (If it contains an o-cycle it can't be R-sorted is trivial since this would change the o-cycle to an e-cycle; the other direction for a connected graph is Lemma 2) in the ref above

int hasoloops(int **mat, int *resgn, int n, int m){

  int i, j, k, p, q,flg=1;
  resgn[0]=1;
  for(i=1;i<n;i++)
    resgn[i]=0;

  while(flg){
    flg=0;
    for(j=0;j<m;j++){// each column
      p=0;// no nonzero entries found yet
      for(i=0;i<n;i++){ // each entry
	if(!p && mat[i][j]){// first nonzero entry in col j
	  p=mat[i][j];q=i;
	}
	else if(mat[i][j]*p > 0){ // second entry of same sign in col j
	  if(!(resgn[q]) && !(resgn[i])){
	    resgn[i]=-1;
	    resgn[q]=1;
	    for(k=0;k<m;k++) // resign row i
	      mat[i][k]=-mat[i][k];
	    flg=1;// something resigned
	    break;
	  }
	  else if(!(resgn[i])){
	    resgn[i]=-1;
	    for(k=0;k<m;k++) // resign row i
	      mat[i][k]=-mat[i][k];
	    flg=1; // something resigned
	    break;
	  }
	  else if(!(resgn[q])){
	    resgn[q]=-1;
	    for(k=0;k<m;k++) // resign row q
	      mat[q][k]=-mat[q][k];
	    flg=1;// something resigned
	    break;
	  }
	  else // can't resign: failed
	    return 1;
	}
      }
      if(flg) // some row got resigned, start again. 
	break;
    }
  }// exit loop if no changes

  for(i=0;i<n;i++)//non-resigned rows
    if(!(resgn[i]))
      resgn[i]=1;
  return 0;
}

int hasoloopst(int **mat, int *resgn, int n, int m){
  int **tmp=transposemat(mat,n,m);
  int i,j,flg;
  flg=hasoloops(tmp,resgn,m,n);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      mat[i][j]=tmp[j][i];

  free_imatrix(tmp,0,m-1,0,n-1);
  return flg;
}





// symmetrise a pattern matrix

int **symmetrise(int **A, int n){
  int i,j;
  int **B=imatrix(0, n-1, 0, n-1);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      B[i][j]=0;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(A[i][j]){
	B[i][j]=1;B[j][i]=1;
      }
    }
  }
  return B;
}

bool connected_prod(int **A, int **B, int n, int m){
  // Is the i-graph of the product of AB^t
  // connected (not necessarily strongly)?
  // Assumes that the I-graph is well-defined
  int **C,**C1,**CC=NULL;
  int totCC=0;
  bool flg=0;
  C=multABT(A,B,n,m);
  C1=symmetrise(C,n);
  CC=Tarjan(C1,n,&totCC);
  if(totCC==1)
    flg=1;
  free_imatrix(C,0,n-1,0,n-1);
  free_imatrix(C1,0,n-1,0,n-1);
  free_imat(CC,totCC);
  return flg;
}

// From Graph-theoretic characterizations of monotonicity
// of chemical networks in reaction coordinates
// David Angeli, Patrick De Leenheer and Eduardo Sontag
// JMB (2009)
// Assumes it has been checked that no chemicals appear on both 
// sides of a reaction. 
// Assumes we have checked that stoichiometry classes are bounded
// Assumes structural persistence


int ALS_JMB_2009(int **im1, int **im2, int n, int m, int allgood, int bdclass, int persistflag){
  int *resgn=(int *)malloc((size_t) ((m)*sizeof(int)));
  //local copies
  int **imat1=cpmat(im1, n, m);
  int **imat2=cpmat(im2, n, m);
  int i,j;
  int flg=0;

  if(!allgood || bdclass!=1 || !persistflag)
    return 0;

  //  fprintf(stderr, "rowd = %d\n", rowdegree(imat1, n, m));
  if(rowdegree(imat1, n, m)<=2 && !hasoloopst(imat1, resgn, n, m)){//orthant monotone (gets resigned to preserve the nonnegative orthant)
    //fprintf(stderr, "here\n");
    // resign imat2
    for(i=0;i<n;i++)
      for(j=0;j<m;j++)
	imat2[i][j]=resgn[i]*imat2[i][j];

    // is the I-graph of the reverse product connected?
    if(connected_prod(imat2,imat1,m,n)){// strong monotonicity
      //if(hasposrkervec(imat1,n,m,1)==1)//Condition (9) in ref.
      if(!hasposlimvec(imat1,n,m))//Condition (9) in ref.
	flg=1;
      else if(!hasposrkervec(imat1,n,m,0))//Condition (10) in ref.
	flg=2;
    }
  }
  free_imatrix(imat1,0,n-1,0,m-1);
  free_imatrix(imat2,0,n-1,0,m-1);
  free((char*)resgn);
  return flg;
}


int simple_CRN(int **im, int n, int m){
  int rd,cd,flg=0;
  int **tmp,*resgn;
  int **imat=cpmat(im,n,m);//local copy
  // first for the matrix
  if((rd=rowdegree(imat, n, m))<=2 || (cd=coldegree(imat, n, m))<=2){
    if(rd){
      resgn=(int *)malloc((size_t) ((n)*sizeof(int)));
      if(!hasoloops(imat, resgn, n, m))
	flg=1;
      free((char*)resgn);
    }
    else if(cd){
      tmp=transposemat(imat, n, m);
      resgn=(int *)malloc((size_t) ((m)*sizeof(int)));
      if(!hasoloops(tmp, resgn, m, n))
	flg=2;
      free((char*)resgn);
      free_imatrix(tmp,0,m-1,0,n-1);
    }
    
  }
  free_imatrix(imat,0,n-1,0,m-1);
  return flg;
}

void printsiphons(int **allsiphons, int totsiphons, char **chems){
  int i,j;
  //  fprintf(stderr, "tot=%d\n",totsiphons);
  for(i=0;i<totsiphons;i++){
    fprintf(stderr, "    {");
    for(j=1;j<allsiphons[i][0];j++)
      fprintf(stderr, "%s ", chems[allsiphons[i][j]]);
    fprintf(stderr, "%s}\n", chems[allsiphons[i][j]]);
  }
  fprintf(stderr, "\n");
}



//overloading: add ** to the minimal ones
void printsiphons(int **allsiphons, int totsiphons, int **allminsiphons, int totminsiphons, char **chems){
  int i,j;
  //  fprintf(stderr, "tot=%d\n",totsiphons);
  for(i=0;i<totsiphons;i++){
    fprintf(stderr, "    {");
    for(j=1;j<allsiphons[i][0];j++)
      fprintf(stderr, "%s ", chems[allsiphons[i][j]]);
    fprintf(stderr, "%s}", chems[allsiphons[i][j]]);
    if(isinarray4(allminsiphons, totminsiphons, allsiphons[i])!=-1)
      fprintf(stderr, "*\n");
    else
      fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}

//overloading (default chemical names)
void printsiphons(int **allsiphons, int totsiphons){
  int i,j;
  //  fprintf(stderr, "tot=%d\n",totsiphons);
  for(i=0;i<totsiphons;i++){
    fprintf(stderr, "    {");
    for(j=1;j<allsiphons[i][0];j++)
      fprintf(stderr, "A%d ", allsiphons[i][j]);
    fprintf(stderr, "A%d}\n", allsiphons[i][j]);
  }
  fprintf(stderr, "\n");
}

//overloading: default chemical names, and ** to minimal ones
void printsiphons(int **allsiphons, int totsiphons, int **allminsiphons, int totminsiphons){
  int i,j;
  //  fprintf(stderr, "tot=%d\n",totsiphons);
  for(i=0;i<totsiphons;i++){
    fprintf(stderr, "    {");
    for(j=1;j<allsiphons[i][0];j++)
      fprintf(stderr, "A%d ", allsiphons[i][j]);
    fprintf(stderr, "A%d}", allsiphons[i][j]);
    if(isinarray4(allminsiphons, totminsiphons, allsiphons[i])!=-1)
      fprintf(stderr, "*\n");
    else
      fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}

//Overloading to take Si and Sil as input and print to stdout
void printsiphons(int **Si, int **Sil, int n, int m){
  int i,j,totsiphons, totminsiphons;
  int **allsiphons=NULL, **allminsiphons=NULL;

  totsiphons=checksiphons(Si,Sil,n,m,&allsiphons,&totminsiphons,&allminsiphons);
  if(totsiphons){
    printf("All siphons (minimal ones have *):\n");
    for(i=0;i<totsiphons;i++){
      printf("    {");
      for(j=1;j<allsiphons[i][0];j++)
	printf("A%d ", allsiphons[i][j]);
      printf("A%d}", allsiphons[i][j]);
      if(isinarray4(allminsiphons, totminsiphons, allsiphons[i])!=-1)
	printf("*\n");
      else
	printf("\n");
    }
    printf("\n");
  }
  else
    printf("The network has no siphons.\n");

  free_imat(allsiphons,totsiphons);
  free_imat(allminsiphons,totminsiphons);

  return;
}

//Overloading PN AM input
void printsiphons(int **AM, int n, int m){
  int **S, **Sl;
  int minus=1;
  int entries=nonzentries(AM,n+m,n+m);
  char *str=CRNamtostr(AM, n, m, entries, 0);
  printf("Network:\n%s", str);
  free(str);
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  printsiphons(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return;
}

//overloading di6 input
void printsiphons(char *di6, int n, int m){
  int **S, **Sl;
  int minus=0;
  char *str=di6toreacstr((char*)di6, n, m, 0);
  printf("Network:\n%s", str);
  free(str);
  di6toSSl(di6, n, m, minus, &S, &Sl);
  printsiphons(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return;
}

//molecularity of CRN
int molecularity(int **Si, int **Sil, int n, int m){
  int j;
  int tmpl,molleft=0;
  int tmpr,molright=0;
  for(j=0;j<m;j++){
    tmpl=colsum(Sil,n,m,j);
    tmpr=colsum(Si,n,m,j)+tmpl;
    molleft=max(molleft,tmpl);
    molright=max(molright,tmpr);
  }
  return max(molleft,molright);
}

//Overloading PN AM input
int molecularity(int **AM, int n, int m){
  int **S, **Sl;
  int minus=1,ret;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  ret=molecularity(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading di6 input
int molecularity(char *di6, int n, int m){
  int **S, **Sl;
  int minus=0,ret;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  ret=molecularity(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}


//Print a report on the molecularity of the network
void printmolecularity(int **Si, int **Sil, int n, int m){
  int j,totmol=0, totmolleft=0;
  int molleft[m];
  for(j=0;j<m;j++){
    molleft[j]=colsum(Sil,n,m,j);
    totmolleft+=molleft[j];
    totmol+=colsum(Si,n,m,j)+2*colsum(Sil,n,m,j);
  }
  qsortt(molleft,0,m-1);
  reverse(molleft,m);

  printf("total molecularity: %d\n", totmol);
  printf("left molecularities: ");printvec1(molleft,m);
  printf("total left molecularity: %d\n", totmolleft);
  return;
}

//Overloading PN AM input
void printmolecularity(int **AM, int n, int m){
  int **S, **Sl;
  int minus=1;
  int entries=nonzentries(AM,n+m,n+m);
  char *str=CRNamtostr(AM, n, m, entries, 0);
  printf("Network:\n%s", str);
  free(str);
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  printmolecularity(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return;
}

//overloading di6 input
void printmolecularity(char *di6, int n, int m){
  int **S, **Sl;
  int minus=0;
  char *str=di6toreacstr((char*)di6, n, m, 0);
  printf("Network:\n%s", str);
  free(str);
  di6toSSl(di6, n, m, minus, &S, &Sl);
  printmolecularity(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return;
}



// Returns -2 if unresignable
// or there are signed columns of different signs

int hasnegpos(int **imat, int n, int m){
  int i,j,p,q=0;
  for(j=0;j<m;j++){// each column
    p=0;// nothing found yet
    for(i=0;i<n;i++){
      if(imat[i][j]){
	if(!p && imat[i][j]<0)
	  p=-1;
	else if(!p && imat[i][j]>0)
	  p=1;
	else if(p*imat[i][j]<0){
	  p=0;break;//mixed column
	}
      }
    }
    if(p){//not a mixed column
      if(!q && p<0)
	q=-1;
      else if(!q && p>0)
	q=1;
      else if(q*p<0) // has signed columns of different sign
	return -2;
    }
  }
  return q;

}


// is a row v2 a nonzero multiple of row v1

int get_mult(int *v1, int *v2, int n){
  int i,t=0;
  for(i=0;i<n;i++){
    if(!(v1[i] && v2[i]) && (v1[i] || v2[i]))
      return 0;
    else if(!t && v1[i] && v2[i])
      t=v2[i]/v1[i];
    if(v2[i]!=t*v1[i])
      return 0;
  }
  return t;
}

void addrow(int ***mat, int n, int m){
  if(!n || !(*mat))
    (*mat) = (int**) malloc(sizeof(int*) * 1);
  else
    (*mat) =(int**) realloc((*mat), sizeof(int*)*(n+1));
  (*mat)[n]=(int*) malloc(sizeof(int) * m);
}

// Carry out the factorisation
// partition rows into pairwise linearly dependent sets
// for each such set create a single row in the second factor
// and a column in the first factor
// output the new "middle" dimension and the two matrices

int special_fact(int **mat, int n, int m, int ***F1, int ***F2o){
  int **F2=NULL,**F1t=NULL;
  int i,j,g,t,r=0;
  int dlt[n];
  /* cout << "entering special_fact.\n"; */
  for(i=0;i<n;i++)//initially all rows unused
    dlt[i]=0;

  for(i=0;i<n;i++){//each row
  if(!(dlt[i])){//not yet used: new row and column
    addrow(&F2,r,m);addrow(&F1t,r,n);r++;dlt[i]=1;
      for(j=0;j<m;j++)
	F2[r-1][j]=mat[i][j];
      g=reduce_vec(F2[r-1], m);
      for(j=0;j<i;j++)// set initial entries to zero (definitely all dealt with)
	F1t[r-1][j]=0;
      F1t[r-1][i]=g; // set entry to g
      for(j=i+1;j<n;j++){
	if(dlt[j])//used
	  F1t[r-1][j]=0;
	else if((t=get_mult(F2[r-1],mat[j],m))){//multiple of reduced?
	  F1t[r-1][j]=t;dlt[j]=1;//row now dealt with
	}
	else
	  F1t[r-1][j]=0;
      }
    }
  }
  /* printmat(mat,n,m); */
  /* printmat(F2,r,m); */
  /* printmat(F1t,r,n); */

  (*F2o)=cpmat(F2,r,m);
  (*F1)=transposemat(F1t,r,n);
  free_imat(F2,r);free_imat(F1t,r);
  /* cout << "exiting special_fact.\n"; */
  return r;
}

bool matseq(int **A, int **B, int n, int m){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      if(A[i][j]!=B[i][j])
	return 0;
  return 1;
}


// Routine to check if a matrix has a PM1 (Pete-Murad-1) 
// factorisation or not. The first factor must be a matrix 
// with exactly one nonzero entry in each row. The second 
// must have left nullspace of dimension 1 and including 
// a strictly positive vector

int hasPM1factors(int **Gamma, int **Vpat, int n, int m, int allgood, int debug){
  int i,j,h;
  int **Lambda=NULL, **Theta=NULL;
  int cd,rdim=0;
  int *resgn;
  int **Gamma_new;
  bool genflg=0,scflg=0,Kflg=0;
  if(!allgood)
    return 0;

  rdim=special_fact(Gamma,n,m,&Lambda,&Theta);
  resgn =(int *)malloc((size_t) ((rdim)*sizeof(int)));

  /* printmat(Lambda,n,rdim); */
  /* printmat(Theta,rdim,m); */

  if((cd=coldegree(Theta, rdim, m))<=2 && !hasoloops(Theta, resgn, rdim, m) && matrank(Theta, rdim, m)==rdim-1){//also resigns the second factor
    genflg=1;
    if(debug){
      fprintf(stderr, "This seems to be a monotone system.\n\n");
      //printvec(resgn,rdim);
    }
    //resign the first factor
    for(j=0;j<rdim;j++){//each column
      for(i=0;i<n;i++)//each entry
	Lambda[i][j]=resgn[j]*Lambda[i][j];
    }
    if(debug){
      printmat(Lambda,n,rdim);
      printmat(Theta,rdim,m);
    }
    // if it helps, reverse signs of Lambda and Theta

    if((h=hasnegpos(Lambda,n,rdim))==-1){
      for(j=0;j<rdim;j++){
	for(i=0;i<n;i++)
	  Lambda[i][j]=-Lambda[i][j];
	for(i=0;i<m;i++)
	  Theta[j][i]=-Theta[j][i];
      }
    }
    if(h!=-2)
      Kflg=1;

    //sanity check
    Gamma_new=imatrix(0,n-1,0,m-1);
    multAB(Lambda,Theta,Gamma_new,n,rdim,m);
    if(!matseq(Gamma,Gamma_new,n,m)){
      fprintf(stderr, "something went wrong: matrix product doesn't work out.\n");
      printmat(Gamma, n, m);
      fprintf(stderr, "doesn't factor as:\n\n");
      printmat(Lambda, n, rdim);printmat(Theta, rdim, m);
      exit(0);
    }
    free_imatrix(Gamma_new,0,n-1,0,m-1);

    if(debug){
      fprintf(stderr, "The stoichiometric matrix:\n\n");
      printmat(Gamma, n, m);
      fprintf(stderr, "factors into:\n\n");
      printmat(Lambda, n, rdim);
      printmat(Theta, rdim, m);
    }

    if(debug){
      if(h!=-2)
	fprintf(stderr, "K(Lambda) does not include negative vectors.\n\n");
      else
	fprintf(stderr, "However K(Lambda) includes negative vectors.\n\n");
    }
    scflg=DSRSC(Gamma, Vpat, n, m);
    if(debug){
      if(scflg)      
	fprintf(stderr, "The DSR graph is strongly connected.\n\n");
      else
	fprintf(stderr, "The DSR graph appears not to be strongly connected.\n\n");
    }
  }


  free((char *)resgn);
  free_imatrix(Lambda, 0,n-1,0,rdim-1);
  free_imatrix(Theta,0,rdim-1,0,m-1);

  if(genflg && scflg && Kflg)
    return 1;
  if(genflg && scflg)//not necessarily bounded SCs
    return 2;

  return 0;
}

//guess the format of CRNs in a file
//ret=1: di6, ret=2: reacstr
int guessreacformat(char *fname){
  FILE *fd;
  char line[1000];
  int linelen;
  int ret=0;
  int *reac;
  if(!(fd=fopen(fname, "r"))){
    fprintf(stderr, "ERROR in guessreacformat: could not open \"%s\" for reading.\n", fname);exit(0);
  }
  linelen = gtline(fd, line, 1000);
  while(linelen>0 && iscomline(line))
    linelen = gtline(fd, line, 1000);
  fclose(fd);


  //first non-comment line
  if(linelen>0 && line[0]==38 && (endswith(fname, (char*)".di6") || endswith(fname, (char*)".d6")))//& character, must be di6
    ret=1;
  else if(linelen>0 && isreac(line))
    ret=2;
  else{
    reac=getintvec(line);
    if(reac[0]>=2 && reac[0]%2==0){//could be Sauro or simpstr
      if(endswith(fname, (char*)".Sauro") || endswith(fname, (char*)".sauro"))
	ret=3;
      else if(endswith(fname, (char*)".simp"))
	ret=4;
      else//some integer vector format
	ret=5;
    }
    free((char*)reac);
  }

  return ret;
}

int guessoutputformat(char *fname){
  int ret=0;

  //first non-comment line
  if(endswith(fname, (char*)".di6") || endswith(fname, (char*)".d6"))
    ret=1;
  else if(endswith(fname, (char*)".reac") || endswith(fname, (char*)".reacs"))
    ret=2;
  else if(endswith(fname, (char*)".Sauro") || endswith(fname, (char*)".sauro"))
    ret=3;
  else if(endswith(fname, (char*)".simp"))
    ret=4;
  else if(endswith(fname, (char*)".dot"))
    ret=5;

  return ret;
}



// S is the stoichiometric matrix
// Vpat is the pattern matrix of V
// Si is the irreversible stoichiometric matrix
// Sil is the corresponding matrix of powers (assuming mass action)
// stoichl is the left stoichiometric matrix
// stoichr is the right stoichiometric matrix
// input can be reactions as a string, or a file

int analysereacs(const char inputstr[], int q, bool htmlswitch, bool statswitch){
  int k,**S, **Vpat, **Si, **Sil, **Sa, **stoichl, **stoichr;
  char *str;
  int mlen=0, nlen=0, mlena=0;
  char **chems;
  long pos=0;
  char *line;
  int type=-1;
  int allrev, allgood;
  int mirr;
  int csdflag=0;
  int ssdflag=0;
  int wsdflag=0;
  int compatflag=0;
  int MAcompatflag=0;
  char IC1maxpp[250];
  char IC1maxp[250];
  char IC1max[150];
  char IC1maxstr[500];
  char IC1maxppstr[1000];
  char IC1maxpstr[1000];
  char IC2[150];
  char IC2max[150];
  char IC2maxstr[800];
  char IC1maxIC2maxstr[800];
  char IC1maxppIC2maxstr[800];
  char IC1maxpIC2maxstr[800];
  char failsIC1[500];
  char IC1[200];
  char IC1p[500];
  char IC1pp[500];
  char MAIC2max[500];
  char MAIC1[500];
  char MAIC1p[800];
  char MAIC1pp[800];
  char MAIC1IC2max[800];
  char MAIC1pIC2max[1000];
  char MAIC1ppIC2max[1000];
  char notSSD[500];
  char notWSD[500];
  char notrWSD[600];
  char notrcmpt[600];
  char notrcmpt1[700];
  char feinbergdef0[500];
  char panteadef0[500];
  char feinbergdef1[500];
  char ALSstr[500];
  char PMglob[500];
  char PMloc[500];
  char BanajiPanteastr[300];
  char condstarref[300];
  char deftheor[200];
  char weakrstr[200];
  char defstr[200];
  char siphonstr[200];
  char structpers[200];
  char bddclass[300];
  char unbddclass[300];
  char max1_eq_per_class[200];
  char max1_pos_eq_per_class[200];
  char max1_eq_per_nontriv_class[200];
  int totsiphons=1,totminsiphons=0,**allsiphons=NULL,**allminsiphons=NULL;
  bool persistflag=1,def1flg=0;
  int bdclass, dynnontriv,PC1=0;
  int Srank,numlink,numcomp,deficiency,simp_flg,ALS_flg,PM1_flg;
  bool weakr, haszerocomplex=0,zeronotterm,zeroinitial,SSPO=1,posSSPO=1,zerounique=0;
  int condstar;
  char *tmpstr;
  int report=0;
  int nohopf=0,hopfeffort=0;//minimum effort=0; -1 means don't try
  int normalCRN;
  int degenCRN=0;
  int debug=0;
  int pppdegused=-1;

  //Version x=year 2012+x, .y=month number, .z = revision number
  fprintf(stdout, "Analysereacs version 10.3.0. (Please note that this is work in progress.)\n\n");

  pos=0;
  line = getlinefromstr(&pos, (char*)inputstr);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){//reactions
      type=0;
      str=strdup(inputstr);
    }
    free(line);
    line = getlinefromstr(&pos, (char*)inputstr);
  }
  free(line);

  pos=0;
  if(type==-1 && access(inputstr, F_OK) != -1){//file found
    str=readfileintostr(inputstr);
    if(isonlyspace(str)){
      fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", inputstr);
      exit(0);
    }
    line = getlinefromstr(&pos, str);
    while(strcmp(line, "")!=0){//reactions found
      if(!iscomline(line) && isreac(line))
	type=0; 
      free(line);
      line = getlinefromstr(&pos, str);
    }
    free(line);
  }

  if(type==-1){
    fprintf(stderr, "The format could not be understood. EXITING.\n");exit(0);
  }

  // 
  // str now definitely contains contains reactions
  //

  if(type==0){
    if(!getallreacs(str, &S, &Vpat, &Si, &Sil, &stoichl, &stoichr, &chems, &haszerocomplex, &nlen, &mlen, &mirr, &allrev, &allgood,debug)){
      free(str);
      fprintf(stderr, "ERROR: Couldn't read the reactions in file \"%s\". EXITING. \n", inputstr);
      return -1;
    }

    if(nlen==0 || mlen==0){
      fprintf(stderr, "ERROR: Couldn't read the reactions in file \"%s\". EXITING. \n", inputstr);
      free(str);
      free_imatrix(S, 0, nlen-1, 0, mlen-1);
      free_imatrix(Vpat, 0, nlen-1, 0, mlen-1);
      free_imatrix(Si, 0, nlen-1, 0, mirr-1);
      free_imatrix(Sil, 0, nlen-1, 0, mirr-1);
      free_imatrix(stoichl, 0, nlen-1, 0, mirr-1);
      free_imatrix(stoichr, 0, nlen-1, 0, mirr-1);
      freearraydat(chems, nlen);
      return -1;
    }
    free(str);

 
    Srank=matrank(S, nlen, mlen);
    fprintf(stderr, "The stoichiometric matrix (rank = %d):\n\n", Srank);
    printmat(chems, S, nlen, mlen);
    fprintf(stderr, "The pattern matrix for -V^T:\n\n");
    printmatpat(chems, Vpat, nlen, mlen);
    fprintf(stderr, "The irreversible stoichiometric matrix:\n\n");
    printmat(chems, Si, nlen, mirr);
    fprintf(stderr, "The matrix of powers (mass action):\n\n");
    printmat(chems, Sil, nlen, mirr);
    fprintf(stderr, "The left stoichiometric matrix:\n\n");
    printmat(chems, stoichl, nlen, mirr);
    fprintf(stderr, "The right stoichiometric matrix:\n\n");
    printmat(chems, stoichr, nlen, mirr);

    fprintf(stderr, "\n_________________________________\n\n");


    if(htmlswitch){
      strcpy(bddclass, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Stoichiometric_subspace_and_stoichiometry_classes\" target=\"_blank\">Stoichiometry classes are bounded</a>");
      strcpy(unbddclass, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Stoichiometric_subspace_and_stoichiometry_classes\" target=\"_blank\">Stoichiometry classes are unbounded</a> (the stoichiometric subspace includes a nonnegative vector)");
      strcpy(structpers, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Persistence_condition_2_.28PC2.29\" target=\"_blank\">structurally persistent</a>");
      strcpy(siphonstr, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Siphon\" target=\"_blank\">siphons</a>");
    }
    else{
      strcpy(bddclass, "Stoichiometry classes are bounded");
      strcpy(unbddclass, "Stoichiometry classes are unbounded (the stoichiometric subspace includes a nonnegative vector)");
      strcpy(structpers, "structurally persistent");
      strcpy(siphonstr, "siphons");
    }


    // Is the system weakly reversible?
    weakr=weak_rev(Si, Srank, stoichl, stoichr, nlen, mirr, &numcomp, &numlink, haszerocomplex, &zeronotterm, &zeroinitial, &def1flg, &deficiency, chems, statswitch, htmlswitch, debug);

    //count the siphons
    totsiphons=checksiphons(Si,Sil,nlen,mirr,&allsiphons,&totminsiphons,&allminsiphons);

    if(totsiphons){
      fprintf(stderr, "%s of the system:\n", siphonstr);
      printsiphons(allsiphons, totsiphons, chems);
      fprintf(stderr, "Minimal %s of the system:\n", siphonstr);
      printsiphons(allminsiphons, totminsiphons, chems);
      persistflag=structpersist(Si, nlen, mirr, allminsiphons, totminsiphons);
      if(!persistflag){
	//fprintf(stdout, "The system is not %s (it has critical siphons).\n\n", structpers);
	if(totsiphons==1 && allsiphons[0][0]==nlen)
	  zerounique=1; // only siphon is the one with all species, but this is critical
      }

    }
    else{ //there are no siphons at all (as an indirect consequence stoichiometry classes must be unbounded, for otherwise there would at least be the trivial siphon.)
      PC1=1;
      //	fprintf(stdout, "The system has no %s, and is %s. There are no boundary equilibria.\n\n", siphonstr, structpers);
    }

   // test if the irreversible stoichiometric matrix has no positive vectors in its kernel
    if(allgood && allrev) //simply reversible
      dynnontriv=1;
    else if(weakr) //weakly reversible
      dynnontriv=1;
    else
      dynnontriv=1-hasposlimvec(Si, nlen, mirr);


    /* if(dynnontriv){ */
    /*   int **Q; */
    /*   matrix QX; */
    /*   int numv; */
    /*   matrix JMA=reacJMAeq(Si, Sil, nlen, mirr, &Q, &QX, &numv); */
    /*   cerr << QX << endl; exit(0); */
    /* } */


    bdclass=1;
    if(hasposimvec(S, nlen, mlen))
      bdclass=0;//definitely unbounded

    condstar=DSRCondStar(S,Vpat, nlen, mlen,0,&report,2,1);
    normalCRN=normal(Si, Sil, nlen, mirr, Srank);
    if(dynnontriv)
      degenCRN=reacJMAdegenerate(Si, Sil, nlen, mirr,debug);

    if(dynnontriv && hopfeffort>=0 && nlen<=6 && nonzentries(Sil,nlen,mirr)<=10)
      nohopf=Hopfforbid(Si, Sil, nlen, mirr, "ALL", hopfeffort, &pppdegused, 0);//penultimate argument: effort;

    fprintf(stdout, "Quick report. The network:\n\thas %d species,\n\thas %d (irreversible) reactions,\n\thas rank %d\n\t%s\n\t%s\n\t%s%s\n\thas deficiency %d\n\tis %sweakly reversible\n\thas %d siphons\n\tis %sstructurally persistent\n\t%s Condition *%s.\n", nlen, mirr, Srank,bdclass?"has bounded stoichiometric classes":"has unbounded stoichiometric classes", dynnontriv?"is dynamically nontrivial":"is dynamically trivial", normalCRN?"is normal":"is not normal", (dynnontriv && degenCRN)?"\n\tis degenerate (with MA kinetics, all equilibria are degenerate)":"", deficiency, weakr?"":"not ", totsiphons, persistflag?"":"not ", condstar?"satisfies":"fails", condstar?" perhaps after some manipulation (forbids multiple positive nondegenerate equilibria)":"");
    if(nohopf==2)
      fprintf(stdout, "\tForbids Andronov-Hopf bifurcation on the positive orthant (general kinetics)\n");
    else if(nohopf==1)
      fprintf(stdout, "\tForbids Andronov-Hopf bifurcation on the positive orthant (mass action kinetics)\n");
    fprintf(stdout, "\n");

    /* if(!bdclass) */
    /*   fprintf(stdout, "%s.\n\n",unbddclass); */
    /* else */
    /*   fprintf(stdout, "%s (each stoichiometry class contains at least one equilibrium).\n\n", bddclass); */

    if(totsiphons==1 && bdclass==1 && persistflag){ // full set is a siphon, but corresponds to trivial stoichiometry class
      PC1=2;
      //fprintf(stdout, "The only boundary equilibrium of the system is the trivial one. The system is %s. (The only %s is the set of all species and it is not critical.)\n\n", structpers, siphonstr);
    }

    if((deficiency) && !(def1flg)){
      if(htmlswitch)
	fprintf(stdout, "The network is not deficiency zero and fails the conditions of the deficiency one theorem (see <a href=\"http://reaction-networks.net/wiki/Deficiency_theory\" target=\"_blank\">the wiki</a> for more details.)\n\n");
      else
	fprintf(stdout, "The network is not deficiency zero and fails the conditions of the deficiency one theorem.\n\n");
    }


    //count the siphons
    if(statswitch && !totsiphons)
      fprintf(stdout, "SIPH0\n");
    if(SSPO && posSSPO){ //deficiency zero tests didn't already answer these questions
      if(!dynnontriv && !totsiphons){
	if(!zeroinitial)//zero is an equilibrium
	  fprintf(stdout, "The only equilibrium of the system is a trivial one (the system has no siphons, and no positive equilibria).\n\n");
	else
	  fprintf(stdout, "The system has no equilibria (the system has no siphons, no positive equilibria, and zero is not an equilibrium).\n\n");
      }
      /* else if(dynnontriv==-1 && !totsiphons){//shouldn't execute */
      /* 	if(!zeroinitial)//zero is an equilibrium */
      /* 	  fprintf(stdout, "The only equilibrium of the system appears to be the trivial one.\n\n"); */
      /* 	else */
      /* 	  fprintf(stdout, "The system appears to have no equilibria.\n\n"); */
      /* } */
      else if(!dynnontriv)
	fprintf(stdout, "The system is dynamically trivial: it admits no positive limit sets, as the kernel of the (irreversible) stoichiometric matrix includes no positive vector.\n\n");
      /* else if(dynnontriv==-1)//shouldn't execute */
      /* 	fprintf(stdout, "All equilibria appear to be boundary equilibria.\n\n"); */
      else if(!totsiphons)
	fprintf(stdout, "This system has no %s, and is %s (the boundary includes no equilibria).\n\n", siphonstr, structpers);
      else if(PC1==2)
	fprintf(stdout, "The only siphon is the set of all species and it is not critical: The system is %s, and the only boundary equilibrium of the system is the trivial one.\n\n", structpers);
      else if(persistflag)
	fprintf(stdout, "The system has %s, but is %s: no nontrivial stoichiometry class includes boundary equilibria.\n\n", siphonstr, structpers);
      /* else */
      /* 	fprintf(stdout, "The system is not %s (it has critical siphons).\n\n", structpers); */
    }




    if(htmlswitch){
      strcpy(BanajiPanteastr, "<a href=\"https://epubs.siam.org/doi/10.1137/15M1034441\" target=\"_blank\">Banaji and Pantea, SIAM J. Appl. Dynamical Systems, 15(2), pp807-869, 2016 (arxiv:1309.6771)</a>");
      strcpy(condstarref, "<a href=\"https://arxiv.org/abs/0903.1190\" target=\"_blank\">Banaji and Craciun, Commun. Math. Sci., 7(4), pp867-900, 2010 (arxiv:0903.1190)</a>");
      strcpy(defstr, "<a href=\"http://reaction-networks.net/wiki/Deficiency#Deficiency\" target=\"_blank\">deficiency</a>");
      strcpy(weakrstr, "<a href=\"http://reaction-networks.net/wiki/Weakly_reversible#Weak_reversibility\" target=\"_blank\">weakly reversible</a>");
      strcpy(deftheor, "<a href=\"http://reaction-networks.net/wiki/Deficiency_theory\" target=\"_blank\">deficiency theory</a>");
      strcpy(feinbergdef0, "Theorem 6.1.1 in Feinberg (<a href=\"http://www.sciencedirect.com/science/article/pii/0009250987800994\" target=\"_blank\">Chem. Eng. Sci. 42(10), 1987</a>)");
      strcpy(panteadef0, "Theorem 6.3 and preceding remarks in Pantea (<a href=\"http://epubs.siam.org/doi/abs/10.1137/110840509\" target=\"_blank\">SIAM J. Math. Anal. 44(3), 2012</a>)");
      strcpy(feinbergdef1, "Theorem 6.2.1 in Feinberg (<a href=\"http://www.sciencedirect.com/science/article/pii/0009250987800994\" target=\"_blank\">Chem. Eng. Sci. 42(10), 1987</a>)");
      strcpy(ALSstr, "Theorem 2 in Angeli, De Leenheer and Sontag (<a href=\"http://link.springer.com/article/10.1007/s00285-009-0309-0\" target=\"_blank\">J. Math. Biol. 61(4), 2010</a>)");
      strcpy(PMglob, "Theorem 2.2 in Donnell and Banaji (<a href=\"http://epubs.siam.org/doi/abs/10.1137/120898486\" target=\"_blank\">SIADS, 12(2), 2013</a>)");
      strcpy(PMloc, "Theorem 2.1 in Donnell and Banaji (<a href=\"http://epubs.siam.org/doi/abs/10.1137/120898486\" target=\"_blank\">SIADS, 12(2), 2013</a>)");
      strcpy(IC1maxpp, "<a title=\"IC1''\"href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.27.27_.28IC1.27.27.29\">IC1''</a> and <a title=\"PC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Persistence_condition_1_.28PC1.29\">PC1</a>");
      strcpy(IC1maxp, "<a title=\"IC1''\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.27.27_.28IC1.27.27.29\">IC1''</a> and <a title=\"PC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Persistence_condition_2_.28PC2.29\">PC2</a>");
      strcpy(IC1max, "<a title=\"IC1''\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.27.27_.28IC1.27.27.29\">IC1''</a>");
      strcpy(IC2max, "<a title=\"IC2''\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2.27.27_.28IC2.27.27.29\">IC2''</a>");
      strcpy(IC2, "<a title=\"IC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2_.28IC2.29\">IC2</a>");
      strcpy(IC1pp, "<a title=\"IC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1_.28IC1.29\">IC1</a> and <a title=\"PC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Persistence_condition_1_.28PC1.29\">PC1</a>");
      strcpy(IC1p, "<a title=\"IC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1_.28IC1.29\">IC1</a> and <a title=\"PC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Persistence_condition_2_.28PC2.29\">PC2</a>");
      strcpy(IC1, "<a title=\"IC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1_.28IC1.29\">IC1</a>");
    }
    else{
      strcpy(BanajiPanteastr, "Banaji and Pantea, SIAM J. Appl. Dynamical Systems, 15(2), pp807-869, 2016 (arxiv:1309.6771)");
      strcpy(condstarref, "Banaji and Craciun, Commun. Math. Sci., 7(4), pp867-900, 2010 (arxiv:0903.1190)");
      strcpy(defstr, "deficiency");
      strcpy(weakrstr, "weakly reversible");
      strcpy(deftheor, "deficiency theory");
      strcpy(feinbergdef0, "Theorem 6.1.1 in Feinberg (Chem. Eng. Sci. 42(10), 1987)");
      strcpy(panteadef0, "Theorem 6.3 and preceding remarks in Pantea (SIAM J. Math. Anal. 44(3), 2012)");
      strcpy(feinbergdef1, "Theorem 6.2.1 in Feinberg (Chem. Eng. Sci. 42(10), 1987)");
      strcpy(ALSstr, "Theorem 2 in Angeli, De Leenheer and Sontag (J. Math. Biol. 61(4), 2010)");
      strcpy(PMglob, "Theorem 2.2 in Donnell and Banaji (SIADS, 12(2), 2013)");
      strcpy(PMloc, "Theorem 2.1 in Donnell and Banaji (SIADS, 12(2), 2013)");
      strcpy(IC1maxpp, "IC1'' and PC1");
      strcpy(IC1maxp, "IC1'' and PC2");
      strcpy(IC1max, "IC1''");
      strcpy(IC2max, "IC2''");
      strcpy(IC2, "IC2");
      strcpy(IC1pp, "IC1 and PC1");
      strcpy(IC1p, "IC1 and PC2");
      strcpy(IC1, "IC1");
    }

    // Set some strings associated with IC1*
    if(dynnontriv){ //IC1*, and positive equilibria possible
      if(bdclass)
	strcpy(max1_eq_per_class, "each stoichiometry class other than {0} is nontrivial and includes exactly one equilibrium. This equilibrium is positive"); //bounded classes and no nontrivial boundary equilibria (pp)
      else
	strcpy(max1_eq_per_class, "no stoichiometry class includes more than one equilibrium"); //unbounded classes and no nontrivial boundary equilibria (pp)

      if(bdclass)
	strcpy(max1_eq_per_nontriv_class, "each nontrivial stoichiometry class includes exactly one equilibrium, and this equilibrium is positive"); //bounded classes and structural persistence (p)
      else
	strcpy(max1_eq_per_nontriv_class, "no nontrivial stoichiometry class includes more than one equilibrium, and any equilibrium on a nontrivial stoichiometry class must be positive"); //unbounded classes and structural persistence (p)

      strcpy(max1_pos_eq_per_class, "no stoichiometry class includes more than one positive equilibrium");
      sprintf(failsIC1, "fails condition %s", IC1);
    }
    else{// no positive equilibria
      if(zeroinitial)
	strcpy(max1_eq_per_class, "the system has no equilibria");// (pp and 0 not an equilibrium)
      else
	strcpy(max1_eq_per_class, "the system has no nontrivial equilibria");// (pp)
      strcpy(max1_eq_per_nontriv_class, "no nontrivial stoichiometry class includes any equilibria");// (p)
      // non-structurally persistent systems
      if(bdclass)
	strcpy(max1_pos_eq_per_class, "the system has no positive equilibria (though it does have boundary equilibria)");
      else if(zerounique)
	strcpy(max1_pos_eq_per_class, "the system has no equilibria other than {0}");
      else
	strcpy(max1_pos_eq_per_class, "the system has no positive equilibria");
      sprintf(failsIC1, "fails condition %s (but has no positive equilibria, so this does not imply multiple positive equilibria on a stoichiometry class)", IC1);
    }

    //The general strings


    sprintf(IC1maxpstr, "General kinetics: %s. The system satisfies conditions %s", max1_eq_per_nontriv_class, IC1maxp);
    sprintf(IC1maxpIC2maxstr, "General kinetics: %s. The fully open system is injective. The system satisfies conditions %s and %s", max1_eq_per_nontriv_class, IC1maxp, IC2max);
    sprintf(MAIC1p, "Mass action kinetics: %s. The system satisfies conditions %s", max1_eq_per_nontriv_class, IC1p);
    sprintf(MAIC1pIC2max, "Mass action kinetics: %s. The fully open system is injective. The system satisfies conditions %s and %s", max1_eq_per_nontriv_class, IC1p, IC2max);
    sprintf(IC1maxppIC2maxstr, "General kinetics: %s. The fully open system is injective. The system satisfies conditions %s and %s", max1_eq_per_class, IC1maxpp, IC2max);
    sprintf(IC1maxppstr, "General kinetics: %s. The system satisfies conditions %s", max1_eq_per_class, IC1maxpp);
    sprintf(IC1maxstr, "General kinetics: %s. The system satisfies condition %s", max1_pos_eq_per_class, IC1max);
    sprintf(IC2maxstr, "General kinetics: the system %s. It satisfies condition %s: the fully open system is injective", failsIC1, IC2max);
    sprintf(IC1maxIC2maxstr, "General kinetics: %s, and the fully open system is injective. The system satisfies conditions %s and %s", max1_pos_eq_per_class, IC1max, IC2max);
    sprintf(MAIC1, "Mass action kinetics: %s. The system satisfies condition %s", max1_pos_eq_per_class, IC1);
    sprintf(MAIC1pp, "Mass action kinetics: %s. The system satisfies conditions %s", max1_eq_per_class, IC1pp);
    sprintf(MAIC2max, "Mass action kinetics: the fully open system is injective. The system satisfies condition %s", IC2max);
    sprintf(MAIC1IC2max, "Mass action kinetics: %s, and the fully open system is injective. The system satisfies conditions %s and %s", max1_pos_eq_per_class, IC1, IC2max);
    sprintf(MAIC1ppIC2max, "Mass action kinetics: %s. The fully open system is injective. The system satisfies conditions %s and %s", max1_eq_per_class, IC1pp, IC2max);
    sprintf(notSSD, "The system fails condition %s. There exists a choice of physical power-law kinetics and inflows and outflows such that the fully open system has multiple positive equilibria", IC2);
    sprintf(notrWSD, "Mass action kinetics: the system %s", failsIC1);
    sprintf(notWSD, "Mass action kinetics: there exists a choice of rate constants and inflows and outflows such that the system fails condition %s (the fully open system is noninjective)", IC2);
    sprintf(notrcmpt, "The system %s", failsIC1);
    sprintf(notrcmpt1, "The system %s. There exists a choice of physical power-law kinetics such that the system has multiple positive equilibria on some stoichiometry class", failsIC1);

 

    if(deficiency==0){
      if(weakr){//weakly reversible
	if(totsiphons==0){
	  fprintf(stdout, "This is a %s %s zero network. According to %s, with mass action kinetics: each nontrivial stoichiometry class admits exactly one positive equilibrium, and this equilibrium is locally asymptotically stable relative to its stoichiometry class. In fact, as the system has no siphons, this equilibrium is globally asymptotically stable relative to its stoichiometry class (including the boundary).\n\n", weakrstr, defstr, feinbergdef0);
	}
	else if(persistflag){
	  fprintf(stdout, "This is a %s %s zero network. According to %s, with mass action kinetics: each nontrivial stoichiometry class admits exactly one positive equilibrium, and this equilibrium is locally asymptotically stable relative to its stoichiometry class. In fact, as the system is structurally persistent, this equilibrium is globally asymptotically stable relative to its stoichiometry class (including the boundary).\n\n", weakrstr, defstr, feinbergdef0);
	}
	else if(Srank>3){
	  fprintf(stdout, "This is a %s %s zero network. According to %s, with mass action kinetics: each nontrivial stoichiometry class admits exactly one positive equilibrium, and this equilibrium is locally asymptotically stable relative to its stoichiometry class. There are no positive, nontrivial periodic orbits.\n\n", weakrstr, defstr, feinbergdef0);
	}
	else{
	  fprintf(stdout, "This is a %s %s zero network with stoichiometric subspace of dimension %d. With mass action kinetics: (i) According to %s, each nontrivial stoichiometry class admits exactly one positive equilibrium which attracts the relative interior of its stoichiometry class; (ii) By %s, this equilibrium is also locally asymptotically stable relative to its stoichiometry class. \n\n", weakrstr, defstr, Srank, panteadef0, feinbergdef0);
	}
	strcat(notrWSD, ". By ");strcat(notrWSD, deftheor);
	strcat(notrWSD, " however, with mass action kinetics, the system has precisely one positive equilibrium on each nontrivial stoichiometry class");
      }
      else if(zeronotterm){
	SSPO=0;//no steady states or POs
	fprintf(stdout, "This is a %s zero network and the zero complex does not lie in a terminal strong linkage class. According to %s, for general kinetics there are no equilibria at all and no nontrivial periodic orbits, including on the boundary.\n\n", defstr, feinbergdef0);
	strcat(notrcmpt, ". However, there are no equilibria at all and no nontrivial periodic orbits, including on the boundary (for general kinetics)");
	strcat(notrWSD, ". However, there are no equilibria at all and no nontrivial periodic orbits, including on the boundary (for general kinetics)");
      }
      //https://reaction-networks.net/wiki/Reaction_graph#Deficiency
      else{
	posSSPO=0;// no positive steady states or POs
	fprintf(stdout, "The network has %s zero, but is not %s. According to %s, for general kinetics there are no positive equilibria or positive nontrivial periodic orbits.\n\n", defstr, weakrstr, feinbergdef0);
	if(dynnontriv){//shouldn't execute!
	  strcat(notrcmpt, ". However, there are no positive equilibria at all (for general kinetics) - see ");
	  strcat(notrcmpt, feinbergdef0);
	  strcat(notrWSD, ". However, there are no positive equilibria at all (for general kinetics) - see ");
	  strcat(notrWSD, feinbergdef0);
	}
      }
    }
    else if(def1flg){
      if(weakr){
	fprintf(stdout, "This network is %s and satisfies the conditions of the deficiency one theorem. According to %s, for mass action kinetics, the system has precisely one positive equilibrium on each nontrivial stoichiometry class.\n\n", weakrstr, feinbergdef1);
	strcat(notrWSD, ". By ");strcat(notrWSD, deftheor);
	if(PC1 && bdclass)
	  strcat(notrWSD, " however, with mass action kinetics, each stoichiometry class other than {0} includes a unique equilibrium, and this equilibrium is positive");
	else if(PC1)
	  strcat(notrWSD, " however, with mass action kinetics, each stoichiometry class includes a unique equilibrium, and this equilibrium is positive");
	else if(persistflag)
	  strcat(notrWSD, " however, with mass action kinetics, each nontrivial stoichiometry class includes a unique equilibrium, and this equilibrium is positive");
	else
	  strcat(notrWSD, " however, with mass action kinetics, the system has precisely one positive equilibrium on each nontrivial stoichiometry class");

      }
      else{
	fprintf(stdout, "This network satisfies the conditions of the deficiency one theorem (but is not %s). According to %s, for mass action kinetics, the system has no more than one positive equilibrium on each nontrivial stoichiometry class.\n\n", weakrstr, feinbergdef1);
	strcat(notrWSD, ". By deficiency theory however, with mass action kinetics, the system has no more than one positive equilibrium on each nontrivial stoichiometry class - see ");
	strcat(notrWSD, feinbergdef1);
      }
    }
    else{
      fprintf(stdout, "The network is %s%s.\n\n", weakr?"":"not ", weakrstr);

      fprintf(stdout, "The network has %s %d.\n\n", defstr, deficiency);
    }
    // end of deficiency theory calculations

    if((simp_flg=simple_CRN(S,nlen,mlen))){
      if(simp_flg==1)
	fprintf(stderr, "The system has no o-cycles and S-degree <=2\n\n");
      else
	fprintf(stderr, "The system has no o-cycles and R-degree <=2\n\n");
    }

    /* fprintf(stderr, "numcomp = %d\n", numcomp); */
    /* fprintf(stderr, "numlink = %d\n", numlink); */
    /* fprintf(stderr, "Srank = %d\n", Srank); */

    // check if the stoichiometric matrix has positive vectors in its left-kernel


    if((ALS_flg=ALS_JMB_2009(S,Vpat,nlen,mlen,allgood,bdclass,persistflag))){
      if(ALS_flg==1){
	if(statswitch)
	  fprintf(stdout, "ALS1 (global convergence)\n");
	fprintf(stdout, "According to %s, the system (with general kinetics) is globally convergent in the following sense: all positive initial conditions converge to an equilibrium which is the unique equilibrium on its stoichiometry class.\n\n", ALSstr);
      }
      else if(ALS_flg==2){
	if(statswitch)
	  fprintf(stdout, "ALS2 (generic QC)\n");
	fprintf(stdout, "According to %s, the system (with general kinetics) is generically quasiconvergent: almost all positive initial conditions converge to the set of equilibria (the measure of the set of possibly non-convergent initial conditions is zero).\n\n", ALSstr);
      }
    }

    if((PM1_flg=hasPM1factors(S, Vpat, nlen, mlen, allgood, 0))){
      if(PM1_flg==1){
	if(!totsiphons){
	  if(statswitch)
	    fprintf(stdout, "PMflag1 (global no siph)\n");
	  fprintf(stdout, "According to %s, the system (with general kinetics) is globally convergent in the following sense: all initial conditions converge to an equilibrium which is the unique equilibrium on its stoichiometry class, and is (for stoichiometry classes other than {0}) positive.\n\n", PMglob);
	}
	else if(persistflag){
	  if(statswitch)
	    fprintf(stdout, "PMflag2 (global with siph)\n");
	  fprintf(stdout, "According to %s, the system (with general kinetics) is globally convergent in the following sense: all initial conditions on any nontrivial stoichiometry class converge to an equilibrium which is positive and is the unique equilibrium on its stoichiometry class.\n\n", PMglob);
	}
      }
      else if (PM1_flg==2){
	if(statswitch)
	  fprintf(stdout, "PMflag3 (local)\n");
	fprintf(stdout, "According to %s, the system (with general kinetics) is locally convergent in the following sense: every positive initial condition is locally asymptotically stable.\n\n", PMloc);
      }
    }

    if(condstar==1)//5th arg means not minus
      fprintf(stdout, "Condition * in %s holds%s. The network (with general kinetics) forbids multiple positive nondegenerate equilibria.\n\n", condstarref, report==1?": no e-cycles found":"");
    else if(condstar==2)//Reversification made a difference
      fprintf(stdout, "After reversifying the network, Condition * in %s holds%s. The network (with general kinetics) forbids multiple positive nondegenerate equilibria.\n\n", condstarref, report==1?": no e-cycles found":"");
    else
      fprintf(stdout, "Condition * in %s fails%s\n\n", condstarref, report==-1?": es-cycles found":": there are e-cycles with odd intersection");


    fprintf(stdout, "Claims based on theory in %s:\n\n", BanajiPanteastr);
    //free siphons
    for(k=0;k<totsiphons;k++)
      free ((char *)(allsiphons[k]));
    if(allsiphons)
      free((char *) allsiphons);
    for(k=0;k<totminsiphons;k++)
      free ((char *)(allminsiphons[k]));
    if(allminsiphons)
      free((char *) allminsiphons);

    if(allgood && allrev){ // simply reversible CRN - no reactants on both sides and all reversible
      Sa=redmat(S, nlen, mlen, &mlena);//remove redundant cols
      csdflag=isCSD(Sa, nlen, mlena, q);
      if(csdflag!=2){
	ssdflag=isSSD(Sa, nlen, mlena, q);
	if(ssdflag!=2)
	  wsdflag=doubleisWSD(Sa, nlen, mlena, q);
      }
      if(csdflag==2){//CSD
	if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	  fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1'' and IC2'' with *any stoichiometries*.)\n", IC1maxppIC2maxstr);
	else if(persistflag)
	  fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1'' and IC2'' with *any stoichiometries*.)\n", IC1maxpIC2maxstr);
	else
	  fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1'' and IC2'' with *any stoichiometries*.)\n", IC1maxIC2maxstr);
      }
      else if(ssdflag){
	if(ssdflag==2 && csdflag==1){//r-CSD but not CSD
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1'' with *any stoichiometries*.)\n", IC1maxppIC2maxstr);
	  else if(persistflag)
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1'' with *any stoichiometries*.)\n", IC1maxpIC2maxstr);
	  else
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1'' with *any stoichiometries*.)\n", IC1maxIC2maxstr);
	}
	else if(ssdflag==2){//SSD but not CSD
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n", IC1maxppIC2maxstr);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n", IC1maxpIC2maxstr);
	  else
	    fprintf(stdout, "%s.\n", IC1maxIC2maxstr);
	}
	else if(csdflag==1){//r-CSD (and r-SSD) but not CSD or SSD
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1'' with *any stoichiometries*.) %s.\n", IC1maxppstr, notSSD);
	  else if(persistflag)
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1'' with *any stoichiometries*.) %s.\n", IC1maxpstr, notSSD);
	  else
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1'' with *any stoichiometries*.) %s.\n", IC1maxstr, notSSD);
	}
	else{//r-SSD (but not SSD)
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s. %s.\n", IC1maxppstr, notSSD);
	  else if(persistflag)
	    fprintf(stdout, "%s. %s.\n", IC1maxpstr, notSSD);
	  else
	    fprintf(stdout, "%s. %s.\n", IC1maxstr, notSSD);
	}
      }
      else{// neither SSD nor r-SSD (as reversible, nonzero vector in ker(Gamma_ir) is automatic)
	fprintf(stdout, "General kinetics: %s. %s.\n", notrcmpt1,notSSD);
      }

      if(csdflag!=2 && ssdflag!=2){ //Not strongest conclusions: check mass action
	if(wsdflag==3){ //WSD and r-strongly WSD
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n", MAIC1ppIC2max);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n", MAIC1pIC2max);
	  else
	    fprintf(stdout, "%s.\n", MAIC1IC2max);
	}
	else if(wsdflag==2) //WSD but not r-strongly WSD
	  fprintf(stdout, "%s.\n%s.\n", MAIC2max, notrWSD);
	else if(wsdflag==1){ //r-strongly WSD but not WSD
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n%s.\n", MAIC1pp, notWSD);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n", MAIC1p, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n", MAIC1, notWSD);
	}
      }
      free_imatrix(Sa, 0, nlen-1, 0, mlena-1);
      if(csdflag==0 && ssdflag==0 && wsdflag==0)
	fprintf(stdout, "%s.\n%s.\n", notrWSD,notWSD);
    }
    //finished with the simply reversible case
    //Beginning the general case: some species on both sides, 
    //and/or some irreversible reactions.
    else{
      //if(!allgood || (!allrev && csdflag==0 && ssdflag==0 && wsdflag==0)){
  
      /* flag = 3 means the matrices are compatible and r-strongly compatible */
      /* flag = 2 means the matrices are compatible but not r-strongly compatible */
      /* flag = 1 means the matrices are r-strongly compatible but not compatible */
      /* flag = 0 means the matrices are none of the above */

      compatflag=arecompat(S, Vpat, nlen, mlen, debug);
      if(compatflag==3){//stoichiometric matrix and rate-pattern are compatible and r-strongly compatible -- accordance and concordance
	if(statswitch)
	  fprintf(stdout, "NR***\nGKIC1\nGKIC2max\nMAIC1\nMAIC2max\n");//for batch checking
	if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	  fprintf(stdout, "%s.\n", IC1maxppIC2maxstr);
	else if(persistflag)
	  fprintf(stdout, "%s.\n", IC1maxpIC2maxstr);
	else
	  fprintf(stdout, "%s.\n", IC1maxIC2maxstr);
      }
      else if(compatflag==2){ //stoichiometric matrix and rate-pattern are compatible but not r-strongly compatible -- accordance but not concordance (implying structural discordance)
	if(statswitch)
	  fprintf(stdout, "NR***\nGKIC2max\nMAIC2max\n");//for batch checking
	if(SSPO && posSSPO && dynnontriv)//possibly interior equilibria
	  tmpstr=notrcmpt1;
	else
	  tmpstr=notrcmpt;

	MAcompatflag=mats_compat(Si, Sil, nlen, mirr, q, debug);
	if(MAcompatflag==3 || MAcompatflag==1 || MAcompatflag==-1){// (stoich and exp) matrices are r-strongly compatible or r-strongly negatively compatible
	  if(statswitch)
	    fprintf(stdout, "MAIC1\n");//for batch checking
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC2maxstr, tmpstr, MAIC1pp);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC2maxstr, tmpstr, MAIC1p);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC2maxstr, tmpstr, MAIC1);
	}
	else{
	  if(SSPO && posSSPO && dynnontriv)//possibly interior equilibria
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC2maxstr,notrWSD,notrcmpt1);
	  else
	    fprintf(stdout, "%s.\n%s.\n", IC2maxstr,notrWSD);
	}
      }
      else if(compatflag==1 || compatflag==-1){//  matrix and sign-pattern are r-strongly compatible or r-strongly negatively compatible but are not compatible -- concordance but not accordance
	if(statswitch)
	  fprintf(stdout, "NR***\nGKIC1\nMAIC1\n");//for batch checking
	MAcompatflag=mats_compat(Si, Sil, nlen, mirr, q, debug);
	if(MAcompatflag==3 || MAcompatflag==2){// (stoich and exp) matrices are compatible
	  if(statswitch)
	    fprintf(stdout, "MAIC2max\n");//for batch checking
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1maxppstr, notSSD, MAIC2max);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1maxpstr, notSSD, MAIC2max);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1maxstr, notSSD, MAIC2max);
	}
	else//not compatible
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1maxppstr, notSSD, notWSD);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1maxpstr, notSSD, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1maxstr, notSSD, notWSD);
      }
      else{//neither concordance nor accordance. check MA

	if(SSPO && posSSPO && dynnontriv)//possibly interior equilibria
	  tmpstr=notrcmpt1;
	else
	  tmpstr=notrcmpt;

	MAcompatflag=mats_compat(Si, Sil, nlen, mirr, q, debug);
	if(MAcompatflag==3){// (stoich and exp) matrices are compatible and r-strongly compatible: semiconcordance and semiaccordance
	  if(statswitch)
	    fprintf(stdout, "NR***\nMAIC1\nMAIC2max\n");//for batch checking
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC1ppIC2max);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC1pIC2max);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC1IC2max);
	}
	else if(MAcompatflag==2){ //(stoich and exp) matrices are compatible but not r-strongly compatible or r-strongly negatively compatible: semiaccordance, but not semiconcordance (network is not normal)
	  if(statswitch)
	    fprintf(stdout, "NR***\nMAIC2max\n");//for batch checking
	  fprintf(stdout, "%s.\n%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC2max, notrWSD);
	}
	else if(MAcompatflag==1 || MAcompatflag==-1){ //(stoich and exp) matrices are r-strongly compatible or r-strongly negatively compatible, but not compatible: semiconcordance, but not semiaccordance
	  if(statswitch)
	    fprintf(stdout, "NR***\nMAIC1\n");//for batch checking
	  if(PC1)// no siphons, or bdd classes and only siphon face is {0}
	    fprintf(stdout, "%s.\n%s.\n%s.\n%s.\n", tmpstr, notSSD, MAIC1pp, notWSD);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n%s.\n", tmpstr, notSSD, MAIC1p, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n%s.\n", tmpstr, notSSD, MAIC1, notWSD);
	}
	else{
	  if(statswitch)
	    fprintf(stdout, "NR***\n");//for batch checking
	  //	allminorsigns(Si, Sil, nlen, mirr, q);
	  if(SSPO && posSSPO && dynnontriv)//possibly interior equilibria
	    fprintf(stdout, "%s.\n%s.\n%s.\n%s.\n", notrcmpt1, notSSD, notrWSD, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, notrWSD, notWSD);
	}
      }
    }
 

    free_imatrix(S, 0, nlen-1, 0, mlen-1);
    free_imatrix(Vpat, 0, nlen-1, 0, mlen-1);
    free_imatrix(Si, 0, nlen-1, 0, mirr-1);
    free_imatrix(Sil, 0, nlen-1, 0, mirr-1);
    free_imatrix(stoichl, 0, nlen-1, 0, mirr-1);
    free_imatrix(stoichr, 0, nlen-1, 0, mirr-1);
    freearraydat(chems, nlen);
  }

  fprintf(stdout, "\n");
  return 0;

}

//Overloading to allow di6 input
int analysereacs(const char *di6, int n, int m, int open, int q, bool htmlswitch, bool statswitch){
  int ret;
  char *str=di6toreacstr((char*)di6, n, m, open);
  ret=analysereacs(str, q, htmlswitch, statswitch);
  free(str);
  return ret;
}




void printindvec(char *s, int *v1, int *v2, int k){
  int j1, j2;

  for(j1=0;j1<k;j1++){
    for(j2=0;j2<k;j2++){
      fprintf(stderr, "%s[%d,%d] ", s, v1[j1], v2[j2]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}



int **joinmats(int **imat1, int **imat2, int n, int m){
  int **newmat=NULL;
  int i, j;
  newmat=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat2[i][j]>0 || imat1[i][j]>0)
	newmat[i][j]=1;
      else if(imat1[i][j]<0 || imat2[i][j]<0)
	newmat[i][j]=-1;
      else
	newmat[i][j]=0;
    }
  }
  return newmat;
}

bool insetint(int i, int *j, int sz, int k){
  int m;
  if(k==i)
    return 1;
  for(m=0;m<sz;m++){
    if(k==j[m])
      return 1;
  }
  return 0;
}
void addon(int *v1, int *v2, int n, int k){
  int i;
  for(i=0;i<n;i++)
    v2[i]=v1[i];
  v2[n]=k;
  return;
}

bool rowadjto(bool **p, int j, int *c, int dim){
  int r;
  for(r=0;r<dim;r++)
    if(p[j][c[r]])
      return 1;
  return 0;
}

bool coladjto(bool **p, int j, int *c, int dim){
  int r;
  for(r=0;r<dim;r++)
    if(p[c[r]][j])
      return 1;
  return 0;
}

bool **inttobool(int **imat, int n, int m){
  int i, j;
  bool **bmat=bmatrix(0,n-1,0,m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat[i][j])
	bmat[i][j]=1;
      else
	bmat[i][j]=0;
    }
  }
  return bmat;
}

bool **mattobool(matrix mmat, int n, int m){
  int i, j;
  bool **bmat=bmatrix(0,n-1,0,m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(mmat(i,j)!=0)
	bmat[i][j]=1;
      else
	bmat[i][j]=0;
    }
  }
  return bmat;
}

int **mattoint(matrix mmat, int n, int m){
  int i, j;
  int **imat=imatrix(0,n-1,0,m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(mmat(i,j)!=0)
	imat[i][j]=1;
      else
	imat[i][j]=0;
    }
  }
  return imat;
}



// take a list of integers and a 
// list of integer lists
// for each existing list check if each integer belongs to it
// and if not create a new list which is the concatenation
// of the existing list and the new integer
// The aim is to count matchings

int **addintstostrs(long *k, int *newint, int numnew, int **t, int numstrs, int len)  
{
  int i=0,j=0, m=0,flag;
  long r=0;
  int **out=NULL;

  if(numstrs==0){ // no strings
    out = (int**) malloc(sizeof(int*) * numnew); 
    for(i=0;i<numnew;i++){ // add each new int
      out[i] = (int*) malloc(sizeof(int) * 1);
      out[i][0] = newint[i];
    }
    (*k)=numnew;
    return out;
  }

  while(j<numstrs){ //existing strings
    for(m=0;m<numnew;m++){ // new integers
      flag=1;
      for(i=0;i<len;i++){ // is newint[m] in existing string j?
	if(t[j][i]==newint[m]){
	  flag=0; // yes it is
	  break;
	}
      }
      if(flag){// no it isn't so create concatenated str
	if(r==0){
	  out = (int**) malloc(sizeof(int*) * 1);
	  out[r] = (int*) malloc(sizeof(int) * (len+1));
	}
	else{
	  out=(int**) realloc(out, sizeof(int*) *(r+1));
	  out[r] = (int*) malloc(sizeof(int) * (len+1));
	}
	for(i=0;i<len;i++)
	  out[r][i]=t[j][i];
	out[r][len]=newint[m];
	r++;
      }
    }
    j++;
  }
  (*k)=r;
  return out;
}

// level is the size of the lists
// numlists is the number of lists
// k is the column number 
// r is the number of new lists

int **growlists(int **lists, int level, int numlists, int *rowset, int setlen, int k, bool **bmat, long *r){
  int j,p=0;
  int **newlists=NULL;
  (*r=0);
  int rs[setlen];

  for(j=0;j<setlen;j++){
    if(bmat[rowset[j]][k])
      rs[p++]=j;
  }
  /* cout << "p=" << p << endl; */
  /* cout << "numlists = " << numlists << endl; */
  
  newlists= addintstostrs(r, rs, p, lists, numlists, level);
  /* for(j=0;j<(*r);j++) */
  /*   printvec(newlists[j], level+1); */

  return newlists;
}


//
//for a square bool matrix
//

int **getalllists(int dim, bool **bmat, long *r){
  int **lists=NULL;
  int **newl=NULL;
  int i;
  long j, tot;
  int xc[dim];
  firstcomb(xc, dim, dim);
  *r=0;
  tot=0;
  for(i=0;i<dim;i++){
    cout << "i = " << i << endl;
    cout << "*r = " << *r << endl;
    newl=growlists(lists, i, tot, xc, dim, i, bmat, r);
    if(lists){
      for(j=0;j<tot;j++)
	free((char *) lists[j]);
      free((char *) lists);
    }
    lists=newl; // reset pointer
    tot=(*r);
  }
  return lists;
}

// parity of a permutation

int permpar(int *p, int n){
  int transcount=0;
  int pcpy[n];
  int i, j, p0;
  for(i=0;i<n;i++)
    pcpy[i]=p[i];
  for(i=0;i<n;i++){
    p0=pcpy[i];
    if(p0!=i){
      for(j=i+1;j<n;j++){
	if(pcpy[j]==i){
	  pcpy[i]=i;pcpy[j]=p0;
	  transcount++;
	  break;
	}
      }
    }
  }
  if(transcount%2==0)
    return 1;
  return -1;
}



struct inode
{
  int *iarr;
  int cff;
  struct inode *left;
  struct inode *right;
};

struct inode *talloc(void){
  return (struct inode *)malloc(sizeof(struct inode));
}


int *idup(int *i, int n){
  int j;
  int *k=(int *)malloc((size_t) (n*sizeof(int)));
  for(j=0;j<n;j++)
    k[j]=i[j];
  return k;
}

//
// extracts the coefficient only once
//

struct inode *addtree(struct inode *p, int *w, int n, int *rclist, int dim, int ****J2){
  int cond;
  if(p==NULL){
    p=talloc();
    p->iarr = idup(w,n);
    p->cff=coeff_int(J2, rclist, rclist, dim, w);
    p->left=p->right=NULL;
    //fprintf(stderr, "%d:  ", p->cff);printvec(w,n);
  }
  else if((cond=cmpare(w, p->iarr,n))==0){}// already there
  else if(cond<0)
    p->left=addtree(p->left,w,n,rclist,dim,J2);
  else
    p->right=addtree(p->right,w,n,rclist,dim,J2);
  return p;
}

//
// augments the coefficient
//

struct inode *addtree1(struct inode *p, int cf, int *w, int n){
  int cond;
  if(p==NULL){
    p=talloc();
    p->iarr = idup(w,n);
    p->cff=cf;
    p->left=p->right=NULL;
    //fprintf(stderr, "%d:  ", p->cff);printvec(w,n);
  }
  else if((cond=cmpare(w, p->iarr,n))==0){
    p->cff+=cf;
  }// already there
  else if(cond<0)
    p->left=addtree1(p->left,cf,w,n);
  else
    p->right=addtree1(p->right,cf,w,n);
  return p;
}



void treeprint(struct inode *p, int n){
  if(p!=NULL){
    treeprint(p->left, n);
    fprintf(stderr, "%d,  ", p->cff);
    printvec(p->iarr, n);
    treeprint(p->right, n);
  }
}

long treeread(char *fname, int ***list, int **cff, int *ord){
  FILE *fd;
  char oneline[1000];
  long k=0;
  int j;
  int len=0;
  char *tmp;
  (*ord)=0;
  *list=NULL;
  fd = fopen(fname, "r");
  if(!fd){
    fprintf(stderr, "ERROR in treeread: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  //cout << "here\n";
  while(getline0(fd, oneline, 1000) > 0){
    if(!(iscomline(oneline))){
      tmp=oneline;
      if(k==0){
	//cout << "k = " << k << endl;
	(*ord)=getnumints(oneline);
	if(!(*ord)){
	  fprintf(stderr, "ERROR in treeread: not a valid tree - no integers found on the first line:\n%s", oneline);
	  exit(0);
	}
	else if((*ord)==1){
	  fprintf(stderr, "ERROR in treeread: not a valid tree - only one integers found on the first line:\n%s", oneline);
	  exit(0);
	}

	(*list)=(int**) malloc(sizeof(int*));
	(*cff)=(int*) malloc(sizeof(int));
	(*list)[k]=(int*) malloc(sizeof(int)*((*ord)-1));
	//cout << "k = " << k << endl;
	(*cff)[k]=getint(&tmp);
	//cout << "ord = " << *ord << endl;
	for(j=0;j<(*ord)-1;j++){
	  //cout << "tmp = \"" << tmp << "\"" << endl;
	  (*list)[k][j]=getint(&tmp);

	}
	//cout << "k = " << k << endl;
      }
      else{
	len=getnumints(oneline);
	if(len!=(*ord)){
	  fprintf(stderr, "ERROR in treeread: not a valid tree - bad line:\n%s", oneline);
	  exit(0);
	}
	(*list)=(int**) realloc((*list), sizeof(int*)*(k+1));
	(*list)[k]=(int*) malloc(sizeof(int)*(len-1));
	(*cff)=(int*) realloc((*cff), sizeof(int)*(k+1));
	(*cff)[k]=getint(&tmp);
	for(j=0;j<len-1;j++)
	  (*list)[k][j]=getint(&tmp);
      }
      k++;
    }
  }
  fclose(fd);
  return k;
}




int idiff(int *a, int *b, int n){
  int k=0, k1=0;
  int tot=0;
  while(k<n && k1<n){
    if(a[k]<b[k1]){
      tot++;k++;
    }
    else if(a[k]>b[k1]){
      tot++;k1++;
    }
    else{k++;k1++;}
  }
  return tot;
}

bool idiffis1(int *a, int *b, int n){
  int k=0, k1=0;
  int tot=0;
  while(k<n && k1<n){
    if(a[k]<b[k1]){
      tot++;if(tot>1){return 0;}k++;
    }
    else if(a[k]>b[k1]){
      tot++;if(tot>1){return 0;}k1++;
    }
    else{k++;k1++;}
  }
if(!tot) // same string
    return 0;
  return 1;
}

// compute a determinant given the list of matchings

ex detterm(matrix mmat, int n, int m, int *xc, int *yc, int **lists, int numlists, int dim){
  ex out, tmp;
  int i,j;
  for(i=0;i<numlists;i++){
    tmp=1;
    for(j=0;j<dim;j++)
      tmp*=mmat(xc[lists[i][j]],yc[j]);
    out+=permpar(lists[i], dim)*tmp;
    //out+=tmp;
  }
  return out;
}


//
// identity matrix
//

int **makeidmat(int n){
  int **imat=imatrix(0,n-1,0,n-1);
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i==j)
	imat[i][j]=1;
      else
	imat[i][j]=0;
    }
  }
  return imat;
}


// Reorder the columns of a bool matrix putting the sparsest first
bool **reordercolsparse(bool **bmat, int n, int m){
  bool **outmat=bmatrix(0, n-1, 0, m-1);
  int i,j,k,c,r[m],s[m]; // r[m] are column totals, s[m] is the new order
  for(i=0;i<m;i++){ //each col
    r[i]=0;
    for(j=0;j<n;j++){ // each row
      if(bmat[j][i])
	(r[i])++;
    }
  }
  c=0;k=0;
  while(c<n+1){
    for(i=0;i<m;i++){ // each col
      if(r[i]==c) // total=c
	s[k++]=i;
    }
    c++;
  }
  printvec(s, m);
  for(j=0;j<n;j++)
    for(i=0;i<m;i++)
      outmat[j][i]=bmat[j][s[i]];

  return outmat;
}

// Reorder the columns of a bool matrix, an int matrix and 
// an ex matrix putting the sparsest first

void multireordercolsparse(bool **bmatin, bool ***bmatout, bool **bmat2in, bool ***bmat2out, int **imatin, int ***imatout, matrix mmatin, matrix *mmatout, int n, int m){
  int i,j,k,c,r[m],s[m]; // r[m] are column totals, s[m] is the new order
  (*bmatout)=bmatrix(0, n-1, 0, m-1);
  (*bmat2out)=bmatrix(0, n-1, 0, m-1);
  (*imatout)=imatrix(0, n-1, 0, m-1);
  for(i=0;i<m;i++){ //each col
    r[i]=0;
    for(j=0;j<n;j++){ // each row
      if(bmatin[j][i])
	(r[i])++;
    }
  }
  c=0;k=0;
  while(c<n+1){
    for(i=0;i<m;i++){ // each col
      if(r[i]==c) // total=c
	s[k++]=i;
    }
    c++;
  }
  fprintf(stderr, "reordering vector:\n");
  printvec(s, m);
  for(j=0;j<n;j++){
    for(i=0;i<m;i++){
      (*bmatout)[j][i]=bmatin[j][s[i]];
      (*bmat2out)[j][i]=bmat2in[j][s[i]];
      (*imatout)[j][i]=imatin[j][s[i]];
      (*mmatout)(j,i)=mmatin(j,s[i]);
    }
  }
  return;
}



int **imatfromiimat(int ***iimat, int n1, int m1){
  int **imat=imatrix(0, n1-1, 0, m1-1);
  int i,j;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      if(iimat[i][j][0]==0)
	imat[i][j]=0;
      else
	imat[i][j]=iimat[i][j][1];
    }
  }
  return imat;
}

// The routine converts an integer matrix to a symbolic matrix
// maxv stores the largest variable subscript

matrix Vmatfromimat(int **imat, int n1, int m1, int *maxv){
  matrix exmat(n1, m1);
  char str1[5];
  int i,j;
  (*maxv)=0;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      if(imat[i][j] !=0){
	if(abs(imat[i][j])>(*maxv)) //largest variable subscript
	  (*maxv)=abs(imat[i][j]);

  	if(imat[i][j] >0){
  	  sprintf(str1, "v%d", imat[i][j]);
  	  exmat(i,j)=get_possymbol(str1);
  	}
  	else if(imat[i][j] <0){
  	  sprintf(str1, "v%d", -imat[i][j]);
  	  exmat(i,j)=-get_possymbol(str1);
  	}
      }
    }
  }
  return exmat;
}

matrix Vmatfromiimat(int ***imat, int n1, int m1, int *maxv){
  matrix exmat(n1, m1);
  char str1[5];
  int i,j,k;
  (*maxv)=0;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      exmat(i,j)=0;
      for(k=1;k<=imat[i][j][0];k++){
	if(imat[i][j][k] !=0){
	  if(abs(imat[i][j][k])>(*maxv)) //largest variable subscript
	    (*maxv)=abs(imat[i][j][k]);

	  if(imat[i][j][k] >0){
	    sprintf(str1, "v%d", imat[i][j][k]);
	    exmat(i,j)+=get_possymbol(str1);
	  }
	  else if(imat[i][j][k] <0){
	    sprintf(str1, "v%d", -imat[i][j][k]);
	    exmat(i,j)-=get_possymbol(str1);
	  }
	}
      }
    }
  }
  return exmat;
}

matrix Vmatfrompatmat(int **imat, int n1, int m1, int *maxv){
  matrix exmat(n1, m1);
  char str1[5];
  int i,j;
  (*maxv)=1;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      if(imat[i][j] !=0){
	sprintf(str1, "v%d", (*maxv)++);
  	if(imat[i][j] >0)
  	  exmat(i,j)=get_possymbol(str1);
	else
  	  exmat(i,j)=-get_possymbol(str1);  	
      }
    }
  }
  return exmat;
}

// imat is presumed to be minus the left stoichiometry matrix. We return -V^T
// minus the transpose of the derivative of the reaction rates w.r.t. 
// species concentrations

matrix VmatMAfromlstoich(int **imat, int n1, int m1, int *maxv){
  matrix exmat(n1, m1);
  char str1[5];
  ex rates[m1];
  ex x[n1];
  ex k[m1];
  int i,j;
  // variables
  for(i=0;i<n1;i++){
    sprintf(str1, "v%d", i+1);
    x[i]=get_possymbol(str1);
  }
  // rate constants
  for(i=0;i<m1;i++){
    sprintf(str1, "k%d", i+1);
    k[i]=get_possymbol(str1);
  }
  // reaction rates
  for(j=0;j<m1;j++){
    rates[j]=k[j];
    for(i=0;i<n1;i++){
      if(imat[i][j])
	rates[j]*=pow(x[i],-imat[i][j]);
    }
    //cout << "rates:" << rates[j] << endl;
  }

  // V^T
  for(i=0;i<n1;i++)
    for(j=0;j<m1;j++)
      exmat(i,j)=imat[i][j]*rates[j]/x[i];

  (*maxv)=n1+m1;
  return exmat;
}

// RHS of a MA reaction system from the irreversible stoich matrix
// and the left irreversible stoichiometric matrix
// We can be free with variables as this is only used with polyhomsimp
// which will put it in canonical form

matrix MARHSfromlstoich(int **imatSi, int **imatSl, int n1, int m1){
  matrix rhs(n1,1);
  char str1[5];
  ex rates[m1];
  ex x[n1];
  ex k[m1];
  int i,j;
  // variables
  for(i=0;i<n1;i++){
    sprintf(str1, "v%d", i+1);
    x[i]=get_possymbol(str1);
  }
  // rate constants
  for(i=0;i<m1;i++){
    sprintf(str1, "k%d", i+1);
    k[i]=get_possymbol(str1);
  }
  // reaction rates
  for(j=0;j<m1;j++){
    rates[j]=k[j];
    for(i=0;i<n1;i++){
      if(imatSl[i][j])
	rates[j]*=pow(x[i],-imatSl[i][j]);
    }
    //cout << "rates[" << j << "]: " << rates[j] << endl;
  }

  //
  for(i=0;i<n1;i++){
    rhs(i,0)=0;
    for(j=0;j<m1;j++)
      rhs(i,0)+=imatSi[i][j]*rates[j];
    //cout << rhs(i,0) << endl;
  }

  return rhs;
}



bool hassignedpminors(matrix J, int n){
  int xc[n];
  int k;
  int flag;
  ex detex,ndetex;
  bool outf=1;
  for(k=1;k<=n;k++){
    firstcomb(xc, n, k);flag=1;
    while(flag==1){
      //printvec1(xc,k);
      detex=symbdetsubmat(J, n, n, xc, xc, k);
      ndetex=-detex;
      printvec1(xc,k);
      if(detex==0)
      	cout << "0" << endl;
      else if(detex.info(info_flags::positive))
      	cout << "+" << endl;
      else if(ndetex.info(info_flags::positive))
	cout << detex << endl;
      //cout << "-" << endl;
      /* if(detex!=0 && !(detex.info(info_flags::positive))&& !(ndetex.info(info_flags::positive))){ */
      /* 	cerr << detex << endl; */
      /* 	printsubmat(J,xc,xc,k,k); */
      /* 	cerr << "An apparently unsigned principal minor found.\n"; */
      /* 	outf=0; */
      /* } */


      /* else if(detex==0) */
      /* 	cout << "0" << endl; */
      /* else if(detex.info(info_flags::positive)) */
      /* 	cout << "+" << endl; */
      /* else if(detex.info(info_flags::negative)) */
      /* 	cout << "-" << endl; */

      flag=nextcomb(xc, n, k);
    }
  }
  return outf;
  
}

bool hassignedminors(matrix J, int n){
  int xc[n],yc[n];
  int k;
  int flag,flag1;
  ex detex,ndetex;
  bool outf=1;
  for(k=n-1;k<=n-1;k++){
    firstcomb(xc, n, k);flag=1;
    while(flag==1){
      firstcomb(yc, n, k);flag1=1;
      while(flag1==1){
	printvec1(xc,k);printvec1(yc,k);
	detex=symbdetsubmat(J, n, n, xc, xc, k);
	ndetex=-detex;
	if(detex!=0 && !(detex.info(info_flags::positive))&& !(ndetex.info(info_flags::positive))){
	  cerr << detex << endl;
	  printsubmat(J,xc,xc,k,k);
	  cerr << "An apparently unsigned principal minor found.\n";
	  outf=0;
	}
	else if(detex==0)
	  cout << "0" << endl;
	else if(detex.info(info_flags::positive))
	  cout << "+" << endl;
	else if(detex.info(info_flags::negative))
	  cout << "-" << endl;

	flag1=nextcomb(yc, n, k);
      }
      flag=nextcomb(xc, n, k);
    }
  }
  return outf;
  
}





ex getminorsum(matrix M, int n, int k){
  int xc[n], yc[n];
  int j,flag,flag1;
  ex out=0;
  if(k==0){
    out=1;
    return out;
  }
  firstcomb(xc,n,k);

  flag1=1;
  while(flag1==1){//diagonal terms
    out+=expand(pow(symbdetsubmat(M,n,n,xc,xc,k),2));
    flag1=nextcomb(xc,n,k);
  }

  firstcomb(xc,n,k);
  firstcomb(yc,n,k);
  flag1=1;
  while(flag1==1){//diagonal terms
    flag=1;
    for(j=0;j<k;j++)
      yc[j]=xc[j];
    flag=nextcomb(yc,n,k);
    while(flag==1){
      out+=expand(2*symbdetsubmat(M,n,n,xc,yc,k)*symbdetsubmat(M,n,n,yc,xc,k));
      flag=nextcomb(yc,n,k);
    }
    flag1=nextcomb(xc,n,k);
  }
  return out;      
}

//
// Formats and processing
//

// take two matrices (transpose of the left and right stoichiometric matrices)
// and return in V1 the bipartite labelled digraph in the nauty sense 
// Assumed maximum trimolecular.
// V1 must have size at least 4*(n+m)*(n+m)
void CRNPN2(int **SR, int **RS, int n, int m, bool *V1, int Vlen){
  int i,j;
  int r=n+m;
  int p;

  //initialise
  for(i=0;i<Vlen;i++)
    V1[i]=0;

  //diagonal elements (relationship between layers)
  for(i=0;i<r;i++){
    V1[2*r*(i+r)+i]=1;
    V1[2*r*i+i+r]=1;
  }

  //
  for(i=0;i<n;i++){
    for(j=n;j<r;j++){
      p=SR[i][j-n];
      if(p==1 || p==3)
	V1[2*r*i+j]=1;
      if(p==2 || p==3)
	V1[2*r*(i+r)+j+r]=1;
    }
  }
  for(i=n;i<r;i++){
    for(j=0;j<n;j++){
      p=RS[j][i-n];
      if(p==1 || p==3)
	V1[2*r*i+j]=1;
      if(p==2 || p==3)
	V1[2*r*(i+r)+j+r]=1;
    }
  }
  return;
}



// take two matrices and return the bipartite labelled digraph 
// in the nauty sense (replaced edge-labels with vertex labels)
// V1 must have size at least l^2*(n+m)*(n+m), where l is the number of layers
// Using total stoichiometry, not maximum stoichiometry to set layers
bool *CRNPN3(int **SR, int **RS, int n, int m, int minlayers, int *layers){
  int i,j,k,l;
  int r=n+m;
  int r2;
  int dim0=r*r;
  int V0[dim0];
  int *vec;
  int maxstoich=0;//largest edge label
  bool *V1;
  int Vlen;

  //first the adjacency matrix of the edge-labelled bipartite digraph
  inittozero(V0,dim0);

  for(i=0;i<n;i++){
    for(j=n;j<r;j++){
      V0[r*i+j]=SR[i][j-n];
      maxstoich=SR[i][j-n]>maxstoich?SR[i][j-n]:maxstoich;
    }
  }
  for(i=n;i<r;i++){
    for(j=0;j<n;j++){
      V0[r*i+j]=RS[j][i-n];
      maxstoich=RS[j][i-n]>maxstoich?RS[j][i-n]:maxstoich;
    }
  }
  //printvec(V0,dim0);

  l=0;
  while(maxstoich>0){
    l++;
    maxstoich/=2;
  }
  if(l<minlayers)
    l=minlayers;
  (*layers)=l;

  r2=l*l*r*r;
  Vlen=(r2%6==0)?r2:r2+6-r2%6;
  V1=(bool *)malloc((size_t) (Vlen*sizeof(bool)));
  vec=(int *)malloc((size_t) (l*sizeof(int)));

  // Now the nauty lift
  for(i=0;i<Vlen;i++)//initialise
    V1[i]=0;

  //diagonals
  for(i=0;i<l-1;i++){
    for(j=i+1;j<l;j++){
      for(k=0;k<r;k++){
	V1[l*r*(i*r+k)+j*r+k]=1;
	V1[l*r*(j*r+k)+i*r+k]=1;
      }
    }
  }

  //main
  for(i=0;i<r;i++){
    for(j=0;j<r;j++){
      basek(vec, l, 2, V0[r*i+j]);
      for(k=0;k<l;k++){
	if(vec[l-k-1])
	  V1[l*r*(k*r+i)+k*r+j]=1;
      }
    }
  }
  free((char*)vec);
  return V1;
}


//Sauro multigraph format description at https://arxiv.org/pdf/0901.3067.pdf
//E.g.: 3 reactions (vertices 0,1,2) 5 species (vertices 3,4,5,6,7)
//These are the first two digits; followed by edges (vertex-pairs, S-R and R-S)
//E.g.: 3 5 0 3 0 4 5 0 7 0 6 1 1 7 1 7 2 6 7 2 7 2

//A reaction Sauro format to one in layered di6 format
char *Saurotodi6(char *s, int minlayers, int *layers, int *n, int *m){
  int i, l;
  int dg1, dg2;//assume reac then species
  int **Sl, **Sr;
  int r1, r2, clen, Vlen;
  char *out=NULL;
  bool *V;
  int debug=0;
  int *reac=getintvec(s);

  if(reac[0]<=4 || !reac[1] || !reac[2] || reac[0]%2!=0){
    fprintf(stderr, "ERROR in Saurotodi6: \"%s\" is not a valid CRN.\n", s);
    exit(0);
  }

  (*n)=reac[2];(*m)=reac[1];
  //initialise matrices
  Sl=imatrix(0, (*n)-1, 0, (*m)-1);Sr=imatrix(0, (*n)-1, 0, (*m)-1);
  inittozero(Sl,(*n),(*m));inittozero(Sr,(*n),(*m));

  //create left and right stoichiometric matrices
  for(i=2;i<reac[0];i+=2){
    dg1=reac[i+1];dg2=reac[i+2];
    if(dg1||dg2){
      if(dg1<dg2)//R-to-S
	(Sr[dg2-(*m)][dg1])++;
      else//S-to-R
	(Sl[dg1-(*m)][dg2])++;
    }
  }
  free((char*)reac);

  if(debug){printmat(Sl,(*n),(*m));printmat(Sr,(*n),(*m));}

  V=CRNPN3(Sl,Sr,(*n),(*m),minlayers,&l);
  (*layers)=l;
  //fprintf(stderr, "layers=%d\n", l);

  if((*n)+(*m)>62/l){
    fprintf(stderr, "ERROR in Saurotodi6: Total CRN size too large. Exiting.\n");
    exit(0);
  }

  //assumes l layers
  r1=(*n)+(*m);r2=l*l*r1*r1;
  Vlen=(r2%6==0)?r2:r2+6-r2%6;
  clen=Vlen/6;
  out=(char *)malloc((size_t) ((clen+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+l*r1;
  out[2+clen]=0;

  amtodig(V,clen,out);
  if(debug){fprintf(stderr, "out = %s\n", out);}

  free((char *)V);
  free_imatrix(Sl, 0, (*n)-1, 0, (*m)-1);
  free_imatrix(Sr, 0, (*n)-1, 0, (*m)-1);
  return out;
}

//A reaction Sauro format to CRN PN AM
int **SaurotoAM(char *s, int *n, int *m){
  int i;
  int dg1, dg2;//assume reac then species
  int debug=0;
  int *reac=getintvec(s);
  int **AM;

  if(reac[0]<=4 || !reac[1] || !reac[2] || reac[0]%2!=0){
    fprintf(stderr, "ERROR in SaurotoAM: \"%s\" is not a valid CRN.\n", s);
    exit(0);
  }

  (*n)=reac[2];(*m)=reac[1];//dimensions
  //initialise matrices
  AM=imatrix(0, (*n)+(*m)-1, 0, (*n)+(*m)-1);
  inittozero(AM,(*n)+(*m),(*n)+(*m));

  //create AM
  for(i=2;i<reac[0];i+=2){
    dg1=reac[i+1];dg2=reac[i+2];
    if(dg1||dg2){
      if(dg1<dg2)//R-to-S
	(AM[dg1+(*n)][dg2-(*m)])++;
      else//S-to-R
	(AM[dg1-(*m)][dg2+(*n)])++;
    }
  }
  free((char*)reac);

  if(debug){printmat(AM,(*n)+(*m),(*n)+(*m));}

  return AM;
}



//A reaction simpstr format to one in layered di6 format
//simpstr format: we assume knowledge of the number of species and reactions
//a list of integers separated by spaces or commas: first column of Sl, then 
//first column of Sr, then second column of Sl, then second column of Sr..., 
//and so forth. If there are n species, then each block of 2*n entries
//corresponds to a reaction. 
char *simpstrtodi6(char *s, int minlayers, int *layers, int n, int m){
  int i, l;
  int nR,nS;
  int **Sl, **Sr;
  int r1, r2, clen, Vlen;
  char *out=NULL;
  bool *V;
  int debug=0;
  int *reac=getintvec(s);

  if(reac[0]!=2*n*m){
    fprintf(stderr, "ERROR in simpstrtodi6: \"%s\" is not a valid CRN.\n", s);
    exit(0);
  }

  nS=n;nR=m;
  //initialise matrices (not necessary)
  Sl=imatrix(0, nS-1, 0, nR-1);
  Sr=imatrix(0, nS-1, 0, nR-1);
  inittozero(Sl,nS,nR);
  inittozero(Sr,nS,nR);

  if(debug){printvec(reac+1,reac[0]);}
  //create left and right stoichiometric matrices
  for(i=0;i<reac[0];i++){
    if((i/n)%2==0)//Sl
      Sl[i-n*(i/n)][(i/n)/2]=reac[i+1];
    else//Sr
      Sr[i-n*(i/n)][(i/n-1)/2]=reac[i+1];
  }
  free((char*)reac);

  if(debug){printmat(Sl,nS,nR);printmat(Sr,nS,nR);}

  V=CRNPN3(Sl,Sr,nS,nR,minlayers,&l);
  (*layers)=l;
  //fprintf(stderr, "layers=%d\n", l);

  if(nS+nR>62/l){
    fprintf(stderr, "ERROR in Saurotodi6: Total CRN size too large. Exiting.\n");
    exit(0);
  }

  //assumes l layers
  r1=nS+nR;r2=l*l*r1*r1;
  Vlen=(r2%6==0)?r2:r2+6-r2%6;
  clen=Vlen/6;
  out=(char *)malloc((size_t) ((clen+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+l*r1;
  out[2+clen]=0;

  amtodig(V,clen,out);
  if(debug){fprintf(stderr, "out = %s\n", out);}

  free((char *)V);
  free_imatrix(Sl, 0, nS-1, 0, nR-1);
  free_imatrix(Sr, 0, nS-1, 0, nR-1);

  return out;
}


//CRN PN AM of a reaction simpstr format
int **simpstrtoAM(char *s, int n, int m){
  int i;
  int **AM;
  int debug=0;
  int *reac=getintvec(s);

  if(reac[0]!=2*n*m){
    fprintf(stderr, "ERROR in simpstrtoAM: \"%s\" is not a valid CRN.\n", s);
    exit(0);
  }

  AM=imatrix(0,n+m-1,0,n+m-1);
  inittozero(AM,n+m,n+m);
  if(debug){printvec(reac+1,reac[0]);}
  //create left and right stoichiometric matrices
  for(i=0;i<reac[0];i++){
    if((i/n)%2==0)//Sl
      AM[i-n*(i/n)][(i/n)/2+n]=reac[i+1];
    else//Sr^t
      AM[(i/n-1)/2+n][i-n*(i/n)]=reac[i+1];
  }
  free((char*)reac);

  if(debug){printmat(AM,n+m,n+m);}

  return AM;
}




// Take a reaction as a single line in Sauro format
// Print it as a human readable CRN to file fout. 
int Saurotoreacprint(char *s, char *fout){
  int i, j, jj, k;
  char *tmp;
  int nR,nS;
  int dg1, dg2;//assume reac then substrate
  int **matl, **matr;
  FILE *fd;
  bool flg;
  fd = fopen(fout, "w");
  if(!fd){
    fprintf(stderr, "ERROR in Saurotoreacprint: \"%s\" could not be opened for writing.\n", fout);
    exit(0);
  }
  //create a left array and a right array
  i=0;k=0;j=0;
  while(s[k] && !isdigit((int) s[k])){k++;}
  while(s[k] && isdigit((int) s[k])){j++;k++;}
  tmp=strchop2(s, k-j, j);nR=atoi(tmp);free(tmp);//number of reactions
  j=0;
  while(s[k] && !isdigit((int) s[k])){k++;}
  while(s[k] && isdigit((int) s[k])){j++;k++;}
  tmp=strchop2(s, k-j, j);nS=atoi(tmp);free(tmp);//number of species
  j=0;

  if(!nR || !nS){
    fprintf(stderr, "ERROR in Saurotoreacprint: \"%s\" is not a valid CRN as either no reactions or no species.\n", s);
    exit(0);
  }

  //initialise matrices
  matl=imatrix(0, nS-1, 0, nR-1);
  matr=imatrix(0, nS-1, 0, nR-1);
  for(i=0;i<nS;i++){
    for(jj=0;jj<nR;jj++){
      matl[i][jj]=0;matr[i][jj]=0;
    }
  }

  //create left and right stoichiometric matrices
  while(s[k]){
    while(s[k] && !isdigit((int) s[k])){k++;}
    while(s[k] && isdigit((int) s[k])){j++;k++;}
    tmp=strchop2(s, k-j, j);dg1=atoi(tmp);free(tmp);
    j=0;

    while(s[k] && !isdigit((int) s[k])){k++;}
    while(s[k] && isdigit((int) s[k])){j++;k++;}
    tmp=strchop2(s, k-j, j);dg2=atoi(tmp);free(tmp);
    j=0;

    if(dg1<dg2)//R-to-S
      (matr[dg2-nR][dg1])++;
    else//S-to-R
      (matl[dg1-nR][dg2])++;
  }

  //output to file
  for(i=0;i<nR;i++){
    flg=0;
    for(j=0;j<nS;j++){
      if(matl[j][i]){
	if(!flg){
	  if(matl[j][i]!=1)
	    fprintf(fd, "%dA%d",matl[j][i],j);
	  else
	    fprintf(fd, "A%d",j);
	  flg=1;
	}
	else{
	  if(matl[j][i]!=1)
	    fprintf(fd, "+%dA%d",matl[j][i],j);
	  else
	    fprintf(fd, "+A%d",j);
	}
      }
    }
    fprintf(fd, " --> ");
    flg=0;
    for(j=0;j<nS;j++){
      if(matr[j][i]){
	if(!flg){
	  if(matr[j][i]!=1)
	    fprintf(fd, "%dA%d",matr[j][i],j);
	  else
	    fprintf(fd, "A%d",j);
	  flg=1;
	}
	else{
	  if(matr[j][i]!=1)
	    fprintf(fd, "+%dA%d",matr[j][i],j);
	  else
	    fprintf(fd, "+A%d",j);
	}
      }
    }
    fprintf(fd, "\n");
  }

  free_imatrix(matl, 0, nS-1, 0, nR-1);
  free_imatrix(matr, 0, nS-1, 0, nR-1);

  fclose(fd);
  return 0;
}

// Read a list of reactions in Sauro Format 
// (http://128.208.17.26/NetworkEnumeration)
// and output these as individual files

int SauroSplit(const char *fname, const char *prefix){
  FILE *fd;
  int i=0;
  char oneline[1000];
  char fout[30];
  fd = fopen(fname, "r");
  if(!fd){
    fprintf(stderr, "ERROR in SauroSplit: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  while(getline0(fd, oneline, 1000) > 0){//get each line
    if(!iscomline(oneline)){
      i++;sprintf(fout, "%s%6.6d",prefix,i);
      fprintf(stderr, "fout = %s\n", fout);
      Saurotoreacprint(oneline, fout);
    }
  }
  fclose(fd);
  return 0;
}





// Take the adjacency matrix of a CRN petri net and print as a set of reactions to a file
void CRNamtofile(int **V, int n, int m, const char *outfname, int open){
  FILE *fd;
  int j,k,flg;

  if(!(fd=fopen(outfname, "w"))){
    fprintf(stderr, "ERROR in CRNamtofile: could not open file for writing.\n");exit(0);
  }

  for(k=0;k<m;k++){
    flg=0;
    for(j=0;j<n;j++){
      if(V[j][k+n]){
	if(!flg){
	  if(V[j][k+n]!=1)
	    fprintf(fd, "%dA%d",V[j][k+n],j);
	  else
	    fprintf(fd, "A%d",j);
	  flg=1;
	}
	else{
	  if(V[j][k+n]!=1)
	    fprintf(fd, "+%dA%d",V[j][k+n],j);
	  else
	    fprintf(fd, "+A%d",j);
	}
      }
    }
    fprintf(fd, "-->");
    flg=0;
    for(j=0;j<n;j++){
      if(V[k+n][j]){
	if(!flg){
	  if(V[k+n][j]!=1)
	    fprintf(fd, "%dA%d",V[k+n][j],j);
	  else
	    fprintf(fd, "A%d",j);
	  flg=1;
	}
	else{
	  if(V[k+n][j]!=1)
	    fprintf(fd, "+%dA%d",V[k+n][j],j);
	  else
	    fprintf(fd, "+A%d",j);
	}
      }
    }
    fprintf(fd, "\n");
  }
  if(open){//Add inflow and outflow reacs
    for(j=0;j<n;j++)
      fprintf(fd, "<-->A%d\n",j);
  }

  fclose(fd);
}

// Take the adjacency matrix of a CRN petri net and print as a set of reactions to a string. Allocates memory using number of entries
// Species named A0, A1, ...

char *CRNamtostr(int **V, int n, int m, int entries, int open){
  int j,k,flg;
  char *str;
  unsigned int maxlen=5*entries+5*m+7*n+10;
  if(n<=0||m<=0)
    return NULL;
  str=(char *)malloc((size_t) ((maxlen)*sizeof(char)));
  str[0]=0;

  for(k=0;k<m;k++){
    flg=0;
    for(j=0;j<n;j++){
      if(V[j][k+n]){
	if(!flg){//first species in reac
	  if(V[j][k+n]!=1)
	    sprintf(str+strlen(str), "%dA%d",V[j][k+n],j);
	  else
	    sprintf(str+strlen(str), "A%d",j);
	  flg=1;
	}
	else{
	  if(V[j][k+n]!=1)
	    sprintf(str+strlen(str), "+%dA%d",V[j][k+n],j);
	  else
	    sprintf(str+strlen(str), "+A%d",j);
	}
      }
    }
    sprintf(str+strlen(str), "-->");
    flg=0;
    for(j=0;j<n;j++){
      if(V[k+n][j]){
	if(!flg){//first species on RHS
	  if(V[k+n][j]!=1)
	    sprintf(str+strlen(str), "%dA%d",V[k+n][j],j);
	  else
	    sprintf(str+strlen(str), "A%d",j);
	  flg=1;
	}
	else{
	  if(V[k+n][j]!=1)
	    sprintf(str+strlen(str), "+%dA%d",V[k+n][j],j);
	  else
	    sprintf(str+strlen(str), "+A%d",j);
	}
      }
    }
    sprintf(str+strlen(str), "\n");
  }
  if(open){//Add inflow and outflow reacs
    for(j=0;j<n;j++)
      sprintf(str+strlen(str), "<-->A%d\n",j);
  }

  return str;
}

char *SSltostr(int **S, int **Sl, int n, int m){
  char *str;
  int entries;
  int **AM=SSltoAM(S, Sl, n, m, 0);
  entries=nonzentries(AM,n+m,n+m);
  str=CRNamtostr(AM, n, m, entries, 0);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}

// CRN Petri Net adjacency matrix to a simple string: 
// LHS then RHS of each reaction as row vectors
// To interpret it we need to know the number of species and reactions
char *CRNamtosimpstr(int **V, int n, int m){
  int j,k;
  char *str;
  unsigned int maxlen=6*n*m+1;
  if(n<=0||m<=0)
    return NULL;
  str=(char *)malloc((size_t) ((maxlen)*sizeof(char)));
  str[0]=0;

  for(k=0;k<m;k++){//each reaction
    for(j=0;j<n;j++)//LHS species
      sprintf(str+strlen(str), "%d,",V[j][k+n]);
    for(j=0;j<n;j++)//RHS species
      sprintf(str+strlen(str), "%d,",V[k+n][j]);
  }
  return str;
}

// di6 format to simple string format (see above)
char *di6tosimpstr(char *di6, int n, int m){
  char *str;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=CRNamtosimpstr(AM, n, m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}

// convert a file of di6 reaction digraphs (two-layer)
// to a file of CRNs in simple string format (see above)
unsigned long di6tosimpstrfile(const char *infile, const char *outfile, int n, int m){
  int linelen;
  char *line;
  unsigned long numl=0, tot=0;
  int maxl=0;
  FILE *fd, *fd1;
  char *str, *str1;

  if((numl=numflines(infile, &maxl))<=0){//no reactions
    if(!(fd1=fopen(outfile, "w"))){
	fprintf(stderr, "ERROR in di6tosimpstrfile - File: %s could not be opened...\n", outfile);
	exit(0);
    }
    fclose(fd1);
    return 0;
  }
  line = (char*) malloc(sizeof(char) * (maxl));

  if(!(fd=fopen(infile, "r"))){
    fprintf(stderr, "ERROR in di6tosimpstrfile - File: %s could not be opened...\n", infile);
    exit(0);
  }
  //outfile=NULL means print to stderr
  if(outfile && !(fd1=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in di6tosimpstrfile - File: %s could not be opened...\n", outfile);
    exit(0);
  }

  while((linelen = gtline(fd, line, maxl)) > 0){
    if(strcmp(line,"") && !iscomline(line)){
      if((tot+1)%100000==0)
	fprintf(stderr, "*");
      str=getnthwd(line, 1);
      str1=di6tosimpstr(str,n,m);
      if(outfile)
	fprintf(fd1, "%s\n", str1);
      else
	fprintf(stderr, "%s\n", str1);
      free(str);free(str1);
      tot++;
    }
  }
  fclose(fd);
  if(outfile)
    fclose(fd1);
  free(line);
  return numl;
}

//No error checking on memory. In using this routine, we assume that all 
//networks are genuine. If there are supposed to be n species in some networks, 
//but one is absent, the network will be returned as one with fewer species. 
//Also all reactions are treated as potentially trimolecular at the moment
unsigned long reactodi6file(const char *fname, const char *fout, const char *sepline, int *n, int *m){
  //fixed but surely enough for these purposes
  char line[1000];
  char block[5000];
  int linelen, pass=0;
  FILE *fd, *fd1;
  char *di6;
  int layers;
  unsigned long tot=0;
  int n0=0,m0=0;
  
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "ERROR in reactodi6file: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  if(!(fd1 = fopen(fout, "w"))){
    fprintf(stderr, "ERROR in reactodi6file: \"%s\" could not be opened for writing.\n", fout);
    exit(0);
  }

  linelen=gtline(fd,line,1000);
  while(1){
    if((tot+1)%100000==0)
      fprintf(stderr, "*");
    while(linelen > 0 && !strstr(line, sepline)){
      if(pass==0){strcpy(block, line);pass=1;}else{strcat(block, line);}
      linelen=gtline(fd,line,1000);
    }

    //minimum of two layers (treat all reactions as potentially trimolecular)
    di6=reacstrtodi6(block, 2, &layers, n, m);
    if(!n0){n0=*n;m0=*m;}
    else if(n0!=*n || m0!=*m){
      fprintf(stderr, "WARNING in reactodi6file. The networks do not seem to have consistent dimensions (or maybe not all are genuine). Expecting %d species and %d reactions, but last network had %d species and %d reactions:\n%s\n", n0, m0, *n, *m, block);
    }

    fprintf(fd1, "%s\n", di6);
    tot++;pass=0;
    free(di6);
    if(strstr(line, sepline))
      linelen=gtline(fd,line,1000);
    if(linelen==0)
      break;
  }
  return tot;
}

//Sauro to di6 format
unsigned long Saurotodi6file(const char *fname, const char *fout, int *n, int *m){
  //fixed but surely enough for these purposes
  char line[5000];
  int linelen;
  FILE *fd, *fd1;
  char *di6;
  int layers;
  unsigned long tot=0;
  int n0=0,m0=0;
  
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "ERROR in Saurotodi6file: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  if(!(fd1 = fopen(fout, "w"))){
    fprintf(stderr, "ERROR in Saurotodi6file: \"%s\" could not be opened for writing.\n", fout);
    exit(0);
  }

  linelen=gtline(fd,line,5000);
  while(1){
    if((tot+1)%100000==0)
      fprintf(stderr, "*");
    while(!strcmp(line, "") && iscomline(line))//skip comment/empty lines
      linelen=gtline(fd,line,5000);

    //minimum of two layers (treat all reactions as potentially trimolecular)
    di6=Saurotodi6(line, 2, &layers, n, m);
    if(!n0){n0=*n;m0=*m;}
    else if(n0!=*n || m0!=*m){
      fprintf(stderr, "WARNING in Saurotodi6file. The networks do not seem to have consistent dimensions (or maybe not all are genuine). Expecting %d species and %d reactions, but last network had %d species and %d reactions:\n%s\n", n0, m0, *n, *m, line);
    }

    fprintf(fd1, "%s\n", di6);
    free(di6);
    tot++;linelen=gtline(fd,line,5000);
    if(linelen==0)
      break;
  }
  return tot;
}


//simpstr to di6 format
unsigned long simpstrtodi6file(const char *fname, const char *fout, int n, int m){
  //fixed but surely enough for these purposes
  char line[5000];
  int linelen;
  FILE *fd, *fd1;
  char *di6;
  int layers;
  unsigned long tot=0;
  int n0=0,m0=0;
  
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "ERROR in simpstrtodi6file: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  if(!(fd1 = fopen(fout, "w"))){
    fprintf(stderr, "ERROR in simpstrtodi6file: \"%s\" could not be opened for writing.\n", fout);
    exit(0);
  }

  linelen=gtline(fd,line,5000);
  while(1){
    if((tot+1)%100000==0)
      fprintf(stderr, "*");
    while(!strcmp(line, "") && iscomline(line))//skip comment/empty lines
      linelen=gtline(fd,line,5000);

    //minimum of two layers (treat all reactions as potentially trimolecular)
    di6=simpstrtodi6(line, 2, &layers, n, m);
    if(!n0){n0=n;m0=m;}
    else if(n0!=n || m0!=m){
      fprintf(stderr, "WARNING in simpstrtodi6file. The networks do not seem to have consistent dimensions (or maybe not all are genuine). Expecting %d species and %d reactions, but last network had %d species and %d reactions:\n%s\n", n0, m0, n, m, line);
    }

    fprintf(fd1, "%s\n", di6);
    free(di6);
    tot++;linelen=gtline(fd,line,5000);
    if(linelen==0)
      break;
  }
  return tot;
}

//Wrapper for all the format conversions
unsigned long convertformat(char *infile, int intype, char *outfile, int outtype, int *numspec, int *numreac){
  unsigned long numout=0;
  if(intype==1 && outtype==4)
    numout=di6tosimpstrfile(infile,outfile,*numspec,*numreac);
  else if(intype==1 && outtype==2)
    numout=di6toreacfile(infile,outfile,*numspec,*numreac,0);
  else if(intype==1 && outtype==3)
    numout=di6toSaurofile(infile,outfile,*numspec,*numreac);
  else if(intype==1 && outtype==5)
    numout=di6todotfile(infile,outfile,*numspec,*numreac);
  else if(intype==2 && outtype==1)
    numout=reactodi6file(infile, outfile, "******", numspec, numreac);
  else if(intype==3 && outtype==1)
    numout=Saurotodi6file(infile, outfile, numspec, numreac);
  else if(intype==4 && outtype==1)
    numout=simpstrtodi6file(infile, outfile, *numspec, *numreac);
  else{
    fprintf(stderr, "ERROR in convertformat: couldn't carry out the conversion: perhaps this input-output combination is currently not supported.\n");
    exit(0);
  }
  return numout;
}




// Take the adjacency matrix of a CRN petri net and return 
// the stoichiometric matrix
int **CRNamtostoichmat(int **V, int n, int m){
  int j,k;
  int **stoichmat=imatrix(0, n-1, 0, m-1);

  for(k=0;k<m;k++){//each reac
    for(j=0;j<n;j++)//each species
      stoichmat[j][k]=V[k+n][j]-V[j][k+n];//RHS stoic-LHS stoic
  }
  return stoichmat;
}

// Take the adjacency matrix of a CRN petri net and return 
// the left stoichiometric matrix
int **CRNamtoleftstoichmat(int **V, int n, int m){
  int j,k;
  int **leftstoichmat=imatrix(0, n-1, 0, m-1);

  for(k=0;k<m;k++){//each reac
    for(j=0;j<n;j++)//each species
      leftstoichmat[j][k]=V[j][k+n];//LHS stoic
  }
  return leftstoichmat;
}

// Take the adjacency matrix of a CRN petri net and return 
// the left stoichiometric matrix augmented with a row of ones
int **CRNamtoSil1(int **V, int n, int m){
  int j,k;
  int **leftstoichmat1=imatrix(0, n, 0, m-1);

  for(k=0;k<m;k++){//each reac
    for(j=0;j<n;j++)//each species
      leftstoichmat1[j][k]=V[j][k+n];//LHS stoic
  }
  for(k=0;k<m;k++)
    leftstoichmat1[n][k]=1;

  return leftstoichmat1;
}


// di6 to adjacency matrix as a vector. 
// Assumes graph has fewer than 63 vertices
// Little error checking
int *digraph6toam(const char d6[], int *n){
  int l=strlen(d6);
  int k=0,p;
  //int v[6];
  int *vec=(int *)malloc((size_t) (((l-2)*6)*sizeof(int)));
  if(l<=2 || d6[0]!='&'){
    fprintf(stderr, "ERROR in digraph6toam: %s does not appear to be a well-formed string. EXITING.\n", d6);
    exit(0);
  }
  (*n)=(int)(d6[1])-63;//# vertices
  //fprintf(stderr, "n = %d\n", *n);
  while(k<l-2){
    p=int(d6[k+2])-63;
    //fprintf(stderr, "k=%d, char=%c, p = %d\n", k,d6[k+2],p);
    //basek(v,6,2,p);
    //printvec(v,6);
    basek(vec+6*k, 6, 2, p);
    k++;
  }

  //printvec(vec, (l-2)*6);
  return vec;
}


// Take a CRN in Nauty digraph6 format (two layers) and 
// convert to an adjacency matrix of a CRN petri net graph
// last argument tells us how many entries (sparsity)
int **di6toCRNam(char *di6, int n, int m, int *entries){
  int j,k;
  int *out2;
  int r1=n+m;
  int **V=imatrix(0, r1-1, 0, r1-1);
  int p,p1;
  int n1,partner[r1];
  (*entries)=0;

  out2=digraph6toam(di6, &n1);
  if(n1!=2*r1){
    fprintf(stderr, "Wrong dimensions in di6toCRNam. Exiting. n1=%d, 2*r1 = %d\n", n1, 2*r1);
    exit(0);
  }
  //find partner in upper level: only needed if labelg was used
  for(j=0;j<r1;j++){
    for(k=0;k<r1;k++){
      if(out2[2*r1*j+k+r1]){//partner in next level found
	partner[j]=k+r1-j;
	break;
      }
    }
  }

  for(j=0;j<r1;j++){
    for(k=0;k<r1;k++){
      p=out2[2*r1*j+k];p1=out2[2*r1*(j+partner[j])+k+partner[k]];//
      if(p && p1){V[j][k]=3;(*entries)++;}else if(p){V[j][k]=1;(*entries)++;}else if(p1){V[j][k]=2;(*entries)++;}else{V[j][k]=0;}
    }
  }
  free((char*)out2);
  //printmat(V,r1,r1);

  return V;
}

// Take a CRN in Nauty digraph6 format (l layers) and 
// convert to an adjacency matrix of a CRN petri net graph
// 

int **di6toCRNam1(char *di6, int n, int m, int *totV, int *entries){
  int i,j,k,l;
  int *out2;
  int r1=n+m;
  int **V=imatrix(0, r1-1, 0, r1-1);
  bool *p;
  int n1;
  int **partner;
  (*entries)=0;

  out2=digraph6toam(di6, &n1);
  if(n1%r1!=0){
    fprintf(stderr, "Wrong dimensions in di6toCRNam1. Exiting. n1=%d, r1 = %d\n", n1, r1);
    exit(0);
  }
  l=n1/r1;//number of layers

  partner=imatrix(0, l-1, 0, r1-1);
  //find partners in upper levels
  for(j=0;j<r1;j++){//each base vertex
    for(i=1;i<l;i++){//each level
      for(k=0;k<r1;k++){
	if(out2[l*r1*j+(i*r1+k)]){//partner in level i found
	  partner[i][j]=i*r1+k;
	  //fprintf(stderr, "partner[%d][%d]=%d\n",i,j,i*r1+k);
	  break;
	}
      }
    }
  }

  p=(bool *)malloc((size_t) (l*sizeof(bool)));

  (*totV)=0;
  for(j=0;j<r1;j++){
    for(k=0;k<r1;k++){
      p[0]=out2[l*r1*j+k];
      for(i=1;i<l;i++)
	p[i]=out2[l*r1*partner[i][j]+partner[i][k]];
      if((V[j][k]=binltoint(p,l)))
	(*entries)++;
      //fprintf(stderr, "V[%d][%d]=%d\n", j,k,V[j][k]);
      (*totV)+=V[j][k];
    }
  }
  free((char*)out2);
  free((char*)p);
  free_imatrix(partner, 0, l-1, 0, r1-1);
  //printmat(V,r1,r1);
  //fprintf(stderr, "totV=%d\n",*totV);

  return V;

}


// Get the number of layers in a CRN in Nauty digraph6 format 
int numlayers(char *di6, int n, int m){
  int r1=n+m;
  int n1;
  int *out2=digraph6toam(di6, &n1);
  if(n1%r1!=0){
    fprintf(stderr, "Wrong dimensions in di6toCRNam1. Exiting. n1=%d, r1 = %d\n", n1, r1);
    exit(0);
  }
  free((char*)out2);
  return n1/r1;
}


//di6 format to a human-readable reaction string
char *di6toreacstr1(char *di6, int n, int m, int open){
  char *str;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=CRNamtostr(AM, n, m, entries, open);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}

//l layers
char *di6toreacstr(char *di6, int n, int m, int open){
  char *str;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=CRNamtostr(AM, n, m, entries, open);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}


// Jacobian matrix of (irreversible) reaction system (general kinetics)
// from di6 string or PN adjacency matrix; 
// minus means we want minus the Jacobian matrix

matrix reacJac(int **AM, int n, int m, int *numv, bool minus){
  char str1[5];
  int k,j;
  matrix Dvt(n,m),J;
  int **S=imatrix(0, n-1, 0, m-1);

  (*numv)=0;

  for(k=0;k<m;k++){//k=reaction number
    for(j=0;j<n;j++){//j=species number
      S[j][k]=AM[k+n][j]-AM[j][k+n];//stoich matrix
      if(AM[j][k+n]){
	sprintf(str1, "v%d", (*numv)++ +1);
	Dvt(j,k)=get_possymbol(str1);//transpose of Dv
	if(minus)
	  Dvt(j,k)*=-1;
      }
    }
  }
  J=multABT(S, Dvt, n, m);//Jacobian
  free_imatrix(S, 0, n-1, 0, m-1);
  return J;
}

//overloading
matrix reacJac(char *di6, int n, int m, int *numv, bool minus){
  int entries;
  matrix J;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  J=reacJac(AM,n,m,numv,minus);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return J;
}

//generate stoich matrix, and left stoich matrix (possibly with sign change)
//from CRN Petri-Net AM
void AMtoSSl(int **AM, int n, int m, bool minus, int ***S, int ***Sl){
  int k,j;
  (*S)=imatrix(0, n-1, 0, m-1);
  (*Sl)=imatrix(0, n-1, 0, m-1);
  for(k=0;k<m;k++){//k=reaction number
    for(j=0;j<n;j++){//j=species number
      (*S)[j][k]=AM[k+n][j]-AM[j][k+n];//stoich matrix
      if(!minus)
	(*Sl)[j][k]=AM[j][k+n];//left stoich matrix
      else
	(*Sl)[j][k]=-AM[j][k+n];//minus left stoich matrix
    }
  }
  return;
}

// generate the CRN PN adjacency matrix from stoich and left stoich matrix 
// (possibly with sign change)
int **SSltoAM(int **S, int **Sl, int n, int m, bool minus){
  int k,j;
  int **AM=imatrix(0, n+m-1, 0, n+m-1);
  inittozero(AM,n+m,n+m);
  for(k=0;k<m;k++){//k=reaction number
    for(j=0;j<n;j++){//j=species number
      if(!minus)
	AM[j][k+n]=Sl[j][k];
      else
	AM[j][k+n]=-Sl[j][k];

      AM[k+n][j]=S[j][k]+AM[j][k+n];
    }
  }
  return AM;
}

//di6toSSl
void di6toSSl(char *di6, int n, int m, bool minus, int ***S, int ***Sl){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  AMtoSSl(AM, n, m, minus, S, Sl);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}

//Get the "complex CRN": each complex is the LHS of one reaction with empty RHS
//To check if CRNs are equivalent up to complexes
//sourceonly=0: get all complexes
//sourceonly=1: get only source complexes, one copy of each
//sourceonly=2: get source complexes with multiplicity
char *cmplxCRN(char *di6, int n, int m, int minlayers, int sourceonly, int *totcmplx){
  //stoichl is the left stoichiometric matrix; stoichr is the right stoichiometric matrix
  int i,j,k,totcmplx1=0;
  int **stoichlt, **stoichrt; //for transposed matrices
  int **cmplxs=NULL;
  int **Sl, **Sr;
  int ind1, ind2;
  int l, entries;
  int r1, r2, clen, Vlen;
  char *out=NULL;
  bool *V;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  //  int numl=numlayers(di6, n, m);//maintain the layers
  int numl=minlayers;//...or use minimum

  //complexes
  if(sourceonly==2){//sources only with multiplicity
    totcmplx1=m;
    Sl=imatrix(0, n-1, 0, totcmplx1-1);Sr=imatrix(0, n-1, 0, totcmplx1-1);
    inittozero(Sl,n,totcmplx1);inittozero(Sr,n,totcmplx1);
    for(k=0;k<m;k++){//k=reaction number
      for(j=0;j<n;j++){//j=species number
	Sl[j][k]=AM[j][k+n];
      }
    }
   
  }
  else{//either sources only, or all complexes, one copy each
    //CRN stoichiometric matrices
    stoichlt=imatrix(0, m-1, 0, n-1);stoichrt=imatrix(0, m-1, 0, n-1);
    inittozero(stoichlt,m,n);inittozero(stoichrt,m,n);
    for(k=0;k<m;k++){//k=reaction number
      for(j=0;j<n;j++){//j=species number
	stoichlt[k][j]=AM[j][k+n];
	stoichrt[k][j]=AM[k+n][j];
      }
    }
    for(i=0;i<m;i++){
      totcmplx1=addintstr(totcmplx1, stoichlt[i], n, &cmplxs, &ind1);
      if(!sourceonly)
	totcmplx1=addintstr(totcmplx1, stoichrt[i], n, &cmplxs, &ind2);
    }
    Sl=imatrix(0, n-1, 0, totcmplx1-1);Sr=imatrix(0, n-1, 0, totcmplx1-1);
    inittozero(Sl,n,totcmplx1);inittozero(Sr,n,totcmplx1);
    for(k=0;k<totcmplx1;k++){//k=complex number
      for(j=0;j<n;j++)//j=species number
	Sl[j][k]=cmplxs[k][j];
    }
    free_imatrix(stoichlt,0,m-1,0,n-1);free_imatrix(stoichrt,0,m-1,0,n-1);
    free_imat(cmplxs,totcmplx1);
  }

  V=CRNPN3(Sl,Sr,n,totcmplx1,numl,&l);
  if(n+totcmplx1>62/l){
    fprintf(stderr, "ERROR in cmplxCRN: Total CRN size too large. Exiting.\n");
    exit(0);
  }

  r1=n+totcmplx1;r2=l*l*r1*r1;
  Vlen=(r2%6==0)?r2:r2+6-r2%6;
  clen=Vlen/6;
  out=(char *)malloc((size_t) ((clen+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+l*r1;
  out[2+clen]=0;

  amtodig(V,clen,out);

  //fprintf(stderr, "here\n");
  free((char *)V);
  free_imatrix(Sl,0,n-1,0,totcmplx1-1);	
  free_imatrix(Sr,0,n-1,0,totcmplx1-1);	
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);


  (*totcmplx)=totcmplx1;
  return out;
}



// The main factor in the mass action Jacobian matrix in the special case:
// dynamically nontrivial, 1D kernel
int **reacQMA1(int **AM, int n, int m, bool minus){
  int **S, **Sl, **Sa;
  int **Q;
  int poskervec[m];

  AMtoSSl(AM, n, m, minus, &S, &Sl);

  if(matrank(S,n,m)!=m-1){
    fprintf(stderr, "ERROR in reacQMA1. The rank of the network is not one less than the number of reactions. EXITING.\n");
    exit(0); 
  }
//  if(hasposrkervec(S,n,m,1)!=1){
  if(hasposlimvec(S,n,m)){
    fprintf(stderr, "ERROR in reacQMA1. The network is dynamically trivial. EXITING.\n");
    exit(0); 
  }

  posintkervec(S, n, m, poskervec,1);
  Sa=scalecols(S,n,m,poskervec);
  //printmat(Sa, n, m);
  Q=multABT(Sa,Sl,n,m);

  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  free_imatrix(Sa, 0, n-1, 0, m-1);
  return Q;
}


// Overloading to return the matrix
matrix reacJacMA1(int **AM, int n, int m, bool minus, bool q){
  int **Q=reacQMA1(AM, n, m, minus);
  matrix J;
  matrix xm1(n,1);

  xm1=makexvec(n);
  //printmat(xm1, n, 1);
  J=scalecols(Q,n,n,xm1);
  if(!q){fprintf(stderr, "Jacobian matrix:\n");printmat(J, n, n);}

  free_imatrix(Q, 0, n-1, 0, n-1);

  return J;
}

//The mass action Jacobian matrix at an equilibrium, in terms of 
//the variables and additional positive parameters (not rate constants)
//To be used when we already know the system is dynamically nontrivial
//the case where dim(ker(Si))=1 is special, and in that case we lose the 
//additional variable: it factors out and we can ignore it
//Only in this case do we also return Q
//Also returns the rank of Si 
matrix reacJMAeq(int **Si, int **Sil, int n, int m, int ***Q, matrix *QX, int *numv, int *rk){
  int totbasis;
  matrix bas, Sa, Sb, J;
  int **Sia;
  //Setting this equal to one attempts to minimise the magnitudes of entries 
  //in J: can be important for numerical stability of SDP
  //However, reduction alters Q, and risks, for example, changing Q
  //from negative semidefinite to not so. So we only reduce
  int reduce=1;
  int **basis;
  int **basisT;
  (*Q)=NULL;

  basisT=poskerbasis(Si,n,m,&totbasis,rk);
  if(totbasis>1 && reduce){//Equivalent to replacing Q with QD
    reduce_mat(Sil, n, m, 0);//divide each row by its gcd
  }
  if(!basisT){
    fprintf(stderr, "ERROR in JMAeq: the system appears to be dynamically trivial. EXITING.\n");exit(0);
  }
  matrix xm1(n,1);
  matrix xm2(totbasis,1);
  basis=transposemat(basisT,totbasis,m);

  //fprintf(stderr, "rank of Sil = %d\n", matrank(Sil,n,m));//exit(0);

  //printmat(Si,n,m);
  //printmat(Sil,n,m);
  //printmat(basis, m, totbasis);
  //exit(0);
  //column vectors of unknowns: 1 to n are x_i^{-1} where x_i are the variables
  //n+1 to n+totbasis are the new positive variables
  xm1=makexvec(1,n);
  xm2=makexvec(n+1,totbasis);
  Sb=scalerows(Sil, n, m, xm1);

  if(totbasis==1){//special case
    Sia=scalecols(Si, n, m, basisT[0]);//we need to use basisT here (row vec)
    (*Q)=multABT(Sia,Sil,n,m);
    //printmat((*Q),n,n);
    //if(reduce){reduce_mat((*Q), n, n, 1);}//divide each column by its gcd
    //printmat((*Q),n,n);
    J=multABT(Sia,Sb,n,m);
    free_imatrix(Sia, 0, n-1, 0, m-1);
    (*numv)=n;
  }
  else{
    bas=multAB(basis,xm2,m,totbasis,1);
    Sa=scalecols(Si, n, m, bas);
    (*QX)=multABT(Sa,Sil,n,m);
    //printmat(Si,n,m);
    //cerr<<bas<<endl;
    //cerr<<Sa<<endl;exit(0);
    J=multABT(Sa,Sb,n,m);
    (*numv)=n+totbasis;
  }
  free_imat(basisT,totbasis);
  //free_imatrix(basisT, 0, totbasis-1, 0, m-1); 
  free_imatrix(basis, 0, m-1, 0, totbasis-1);
  return J;
}

//version where we don't need the rank
matrix reacJMAeq(int **Si, int **Sil, int n, int m, int ***Q, matrix *QX, int *numv){
  int rk;
  matrix J=reacJMAeq(Si, Sil, n, m, Q, QX, numv, &rk);
  return J;
}

//overloading: input is PN adjacency matrix
matrix reacJMAeq(int **AM, int n, int m, int ***Q, matrix *QX, int *numv){
  bool minus=0;
  int **S, **Sl;
  matrix J;
  AMtoSSl(AM, n, m, minus, &S, &Sl);

  J=reacJMAeq(S, Sl, n, m, Q, QX, numv);

  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return J;
}

//overloading: input is PN adjacency matrix, also return the rank
matrix reacJMAeq(int **AM, int n, int m, int ***Q, matrix *QX, int *numv, int *rk){
  bool minus=0;
  int **S, **Sl;
  matrix J;
  AMtoSSl(AM, n, m, minus, &S, &Sl);

  J=reacJMAeq(S, Sl, n, m, Q, QX, numv, rk);

  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return J;
}

//deprecated: remove
//version which also returns stoichiometric matrix
matrix reacJMAeqaa(int **AM, int n, int m, int ***S, int ***Q, matrix *QX, int *numv){
  int **Sl=imatrix(0,n-1,0,m-1);
  matrix J;
  bool minus=0;

  AMtoSSl(AM, n, m, minus, S, &Sl);

  J=reacJMAeq((*S), Sl, n, m, Q, QX, numv);
  free_imatrix(Sl,0,n-1,0,m-1);
  return J;
}

//overloading: di6 input
matrix reacJMAeq(char *di6, int n, int m, int ***Q, matrix *QX, int *numv){
  int entries;
  matrix J;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  J=reacJMAeq(AM,n,m,Q,QX,numv);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return J;
}

//overloading: di6 input, also return rank
matrix reacJMAeq(char *di6, int n, int m, int ***Q, matrix *QX, int *numv, int *rk){
  int entries;
  matrix J;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  J=reacJMAeq(AM,n,m,Q,QX,numv,rk);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return J;
}

//print Q or QX (first factor of JMAeq) to stdout
void printreacQ(char *di6, int n, int m){
  int **Q;
  matrix QX;
  matrix J;
  int numv;
  J=reacJMAeq(di6, n, m, &Q, &QX, &numv);

  if(Q){
    printmat1(Q,n,n);
    free_imatrix(Q, 0, n-1, 0, n-1);
  }
  else
    printmat1(QX,n,n);
  
  return;
}


//Print stoichiometric and exponent matrices
void printSSl(char *di6, int n, int m){
  int **S, **Sl;
  int minus=0;
  di6toSSl(di6, n, m, minus, &S, &Sl);
  printmat1(S,n,m);
  printmat1(Sl,n,m);
  return;
}

//print Q or QX (first factor of JMAeq) to stdout
//in format readable by maxima
void printreacQ_max(char *di6, int n, int m){
  int **Q;
  matrix QX;
  matrix J;
  int numv;
  J=reacJMAeq(di6, n, m, &Q, &QX, &numv);

  if(Q){
    printmaximamat1(Q,n,n);
    free_imatrix(Q, 0, n-1, 0, n-1);
  }
  else
    printmaximamat1(QX,n,n);
  
  return;
}

//print JMAeq to stdout
void printreacJMAeq(char *di6, int n, int m){
  int **Q;
  matrix QX;
  matrix J;
  int numv;
  J=reacJMAeq(di6, n, m, &Q, &QX, &numv);
  printmat1(J,n,n);
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return;
}

//print JMAeq to stdout in a form readable by maxima
void printreacJMAeq_max(char *di6, int n, int m){
  int **Q;
  matrix QX;
  matrix J;
  int numv;
  J=reacJMAeq(di6, n, m, &Q, &QX, &numv);
  printmaximamat1(J,n,n);
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return;
}

//print the positive general kinetics Jacobian matrix
void printreacJ(char *di6, int n, int m){
  int numv;
  matrix J=reacJac(di6,n,m,&numv,0);
  printmat1(J,n,n);
  return;
}

//print a surjective matrix W such that W[A|1]=0. 
void printW(int **Sil, int **S, int n, int m, int debug){
  int **W=NULL;
  int tot;
  int rk;
  int i,bez=1;

  W=minA1tkerbasis(Sil,n,m,&tot,&rk,&bez,0,debug);
  fprintf(stderr, "deg = ");
  for(i=0;i<tot;i++)
    fprintf(stderr, "%d ", sumpos(W[i],m));
  fprintf(stderr, "(%d)\n",bez);
  printmat(W, tot, m);
  free_imat(W,tot);
  
  return;
}

//overloading: input is PN adjacency matrix, also return the rank
void printW(int **AM, int n, int m, int debug){
  bool minus=0;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  printW(Sl,S,n,m,debug);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return;
}

//overloading: di6 input
void printW(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  char *str=di6toreacstr((char*)di6, n, m, 0);
  fprintf(stderr, "Network:\n%s", str);
  free(str);
  printW(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}



//CylinderScr: output for mathematica: the better version
//Searching for multiple positive equilibria
char *CylinderScr(int **Si, int **Sil, int n, int m){
  int i,i1,j,t,count=0;
  int numeqs=3;//number of equilibria
  int timeout=10;//maximum seconds
  int **basisT;
  int rk;
  int totbasis;
  char *str=(char *)malloc((size_t) (n*m*100*sizeof(char)));//should be safe
  //printmat(Si,n,m);
  //printmat(Sil,n,m);
  basisT=poskerbasis(Si,n,m,&totbasis,&rk);

  sprintf(str, "Print[TimeConstrained[CylindricalDecomposition[");
  //main part
  for(t=0;t<numeqs;t++){//eq number
    for(j=0;j<m;j++){//reaction
      if(count>0)
	sprintf(str+strlen(str), " && ");
      sprintf(str+strlen(str), "k%d",j+1);
      for(i=0;i<n;i++){
	if(Sil[i][j])
	  sprintf(str+strlen(str), "*x%d%d",i+1,t+1);
      }
      sprintf(str+strlen(str), " == ");
      for(i=0;i<totbasis;i++){
	if(i>0)
	  sprintf(str+strlen(str), "+%d*a%d%d",basisT[i][j],i+1,t+1);
	else
	  sprintf(str+strlen(str), "%d*a%d%d",basisT[i][j],i+1,t+1);
      }
      count++;
    }
  }

  for(j=0;j<m;j++)//reaction
    sprintf(str+strlen(str), " && k%d>0",j+1);

  for(j=0;j<n;j++){//species
    for(t=0;t<numeqs;t++){//eq no
      sprintf(str+strlen(str), " && x%d%d>0",j+1,t+1);
    }
  }

  for(j=0;j<totbasis;j++){//convex coordinates
    for(t=0;t<numeqs;t++){//eq no
      sprintf(str+strlen(str), " && a%d%d>0",j+1,t+1);
    }
  }

  for(t=0;t<numeqs;t++){//eq no
    for(i1=t+1;i1<numeqs;i1++){//eq no
      sprintf(str+strlen(str), " && (");
      for(j=0;j<n;j++){//species
	if(j>0)
	  sprintf(str+strlen(str), " || x%d%d!=x%d%d",j+1,t+1,j+1,i1+1);
	else
	  sprintf(str+strlen(str), "x%d%d!=x%d%d",j+1,t+1,j+1,i1+1);
      }
      sprintf(str+strlen(str), ")");
    }
  }

  sprintf(str+strlen(str), ", {");
  for(j=0;j<m;j++){//reaction
    if(j>0)
      sprintf(str+strlen(str), ", k%d",j+1);
    else
      sprintf(str+strlen(str), "k%d",j+1);
  }
  for(t=0;t<numeqs;t++){//eq no
    for(j=0;j<totbasis;j++){//convex coordinates
      sprintf(str+strlen(str), ", a%d%d",j+1,t+1);
    }
  }
  for(i1=0;i1<numeqs;i1++){//eq no
    for(j=0;j<n;j++){//species
      sprintf(str+strlen(str), ", x%d%d",j+1,i1+1);
    }
  }
  sprintf(str+strlen(str), "}],%d,999]]",timeout);


  fprintf(stderr, "%s\n\n",str);
  //exit(0);
  free(str);
  free_imat(basisT,totbasis);
  return NULL;
  //return str;
}

//overloading: input is PN adjacency matrix, also return the rank
char *CylinderScr(int **AM, int n, int m){
  bool minus=0;
  char *str;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  str=CylinderScr(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return str;
}

//overloading: di6 input
char *CylinderScr(char *di6, int n, int m){
  char *str;
  matrix J;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=CylinderScr(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}



//For nondegenerate (n,n+2,n) CRNs
//Searching for numeq positive equilibria
char *multinn2nscr(int **Si, int **Sil, int n, int m, int numeq){
  int i,j,count=0;
  int timeout=10;//maximum seconds
  int **basisT;
  int rk;
  int totbasis;
  char *str=(char *)malloc((size_t) (n*m*100*sizeof(char)));//should be safe
  int W[m];
  int **Sil1=imatrix(0, n, 0, m-1);
  int degree=0;

  if(m-n!=2){
    fprintf(stderr, "ERROR in \"multinn2nscr\": networks must have 2 more reactions than species.\n");
    exit(0);
  }

  //Augment with a row of ones
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      Sil1[i][j]=Sil[i][j];
    }
  }
  for(j=0;j<m;j++){
    Sil1[n][j]=1;
  }
  
  //printmat(Si,n,m);
  //printmat(Sil1,n+1,m);
  basisT=poskerbasis(Si,n,m,&totbasis,&rk);
  if(rk!=n){
    fprintf(stderr, "ERROR in \"multinn2nscr\": networks must have rank equal to the number of species.\n");
    exit(0);
  }
  if(totbasis!=2){
    fprintf(stderr, "ERROR in \"multinn2nscr\": network must be dynamically nontrivial.\n");
    exit(0);
  }
  intkervec(Sil1, n+1, m, W,1);

  for(j=0;j<m;j++){
    if(W[j]>0)
      degree+=W[j];
  }
  //printmat(basisT,totbasis,m);

  //printvec(W,m);
  //exit(0);


  sprintf(str, "Print[\"degree = %d\"]\n",degree);

  if(degree>=numeq){//only bother if polynomial is of sufficient degree
    sprintf(str+strlen(str), "f[x_]:=");
    count=0;
    for(i=0;i<m;i++){
      if(W[i]>0){
	if(count)
	  sprintf(str+strlen(str), "*");
	if(W[i]==1)
	  sprintf(str+strlen(str), "(%d*x + %d*(1-x))",basisT[0][i],basisT[1][i]);
	else
	  sprintf(str+strlen(str), "(%d*x + %d*(1-x))^%d",basisT[0][i],basisT[1][i],W[i]);
	count++;
      }
    }
    sprintf(str+strlen(str), "-T*");
    count=0;
    for(i=0;i<m;i++){
      if(W[i]<0){
	if(count)
	  sprintf(str+strlen(str), "*");
	if(W[i]==-1)
	  sprintf(str+strlen(str), "(%d*x + %d*(1-x))",basisT[0][i],basisT[1][i]);
	else
	  sprintf(str+strlen(str), "(%d*x + %d*(1-x))^%d",basisT[0][i],basisT[1][i],-W[i]);
	count++;
      }
    }


    //sprintf(str+strlen(str), "\nPrint[TimeConstrained[FindInstance[f[x1]==0&&f[x2]==0&&f[x3]==0&&T>0&&0<x1<x2<x3<1,{T,x1,x2,x3}],%d,999]]",timeout);

    sprintf(str+strlen(str), "\nPrint[TimeConstrained[FindInstance[");
    for(i=0;i<numeq;i++){
      if(i>0)
	sprintf(str+strlen(str), " && ");
      sprintf(str+strlen(str), " f[x%d]==0",i+1);
    }
    sprintf(str+strlen(str), " && T>0 && 0 < ");
    for(i=0;i<numeq;i++)
      sprintf(str+strlen(str), "x%d < ",i+1);

    sprintf(str+strlen(str), "1,{T");
    for(i=0;i<numeq;i++)
      sprintf(str+strlen(str), ", x%d",i+1);

    sprintf(str+strlen(str), "}],%d,999]]",timeout);

  }


  fprintf(stderr, "%s\n\n",str);
  //exit(0);
  free(str);
  free_imat(basisT,totbasis);
  free_imatrix(Sil1,0, n, 0, m-1);
  return NULL;
  //return str;
}


//overloading: input is PN adjacency matrix, also return the rank
char *multinn2nscr(int **AM, int n, int m, int numeq){
  bool minus=0;
  char *str;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  str=multinn2nscr(S, Sl, n, m, numeq);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return str;
}

//overloading: di6 input
char *multinn2nscr(char *di6, int n, int m, int numeq){
  char *str;
  matrix J;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=multinn2nscr(AM,n,m, numeq);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}


//For nondegenerate (n,n+3,n) CRNs
//Searching for numeq positive equilibria
//Assume dynamically nontrivial
char *multinn3nscr(int **Si, int **Sil, int n, int m, int numeq){
  int i,j,count=0;
  int timeout=10;//maximum seconds
  int **basisT;
  int rk, rk1;
  int totbasis;
  char *str=(char *)malloc((size_t) (n*m*100*sizeof(char)));//should be safe
  int **W;
  int totW;
  int **Sil1=imatrix(0, n, 0, m-1);
  int debug=0;
  int deg;

  if(m-n!=3){
    fprintf(stderr, "ERROR in \"multinn3nscr\": networks must have 3 more reactions than species.\n");
    exit(0);
  }

  //Augment with a row of ones
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      Sil1[i][j]=Sil[i][j];
    }
  }
  for(j=0;j<m;j++){
    Sil1[n][j]=1;
  }
  
  //printmat(Si,n,m);
  //printmat(Sil1,n+1,m);
  basisT=poskerbasis(Si,n,m,&totbasis,&rk);
  if(rk!=n){
    fprintf(stderr, "ERROR in \"multinn3nscr\": networks must have rank equal to the number of species.\n");
    exit(0);
  }
  if(totbasis<3){
    fprintf(stderr, "ERROR in \"multinn3nscr\": ker_+(Gamma) should have at least three nonnegative extremal vectors.\n");
    exit(0);
  }

  W=minA1tkerbasis(Sil, n, m, &totW, &rk1, &deg, 0, debug);

  //sprintf(str, "Print[\"degree = %d\"]\n",degree);
  str[0]=0;
  for(j=0;j<2;j++){
    sprintf(str+strlen(str), "f%d[x_,y_]:=",j+1);
    count=0;
    for(i=0;i<m;i++){
      if(W[j][i]>0){
	if(count)
	  sprintf(str+strlen(str), "*");
	if(W[j][i]==1)
	  sprintf(str+strlen(str), "(%d*x + %d*y + %d*(1-x-y))",basisT[0][i],basisT[1][i],basisT[2][i]);
	else
	  sprintf(str+strlen(str), "(%d*x + %d*y + %d*(1-x-y))^%d",basisT[0][i],basisT[1][i],basisT[2][i],W[j][i]);
	count++;
      }
    }
    sprintf(str+strlen(str), "-T%d*",j+1);
    count=0;
    for(i=0;i<m;i++){
      if(W[j][i]<0){
	if(count)
	  sprintf(str+strlen(str), "*");
	if(W[j][i]==-1)
	  sprintf(str+strlen(str), "(%d*x + %d*y + %d*(1-x-y))",basisT[0][i],basisT[1][i],basisT[2][i]);
	else
	  sprintf(str+strlen(str), "(%d*x + %d*y + %d*(1-x-y))^%d",basisT[0][i],basisT[1][i],basisT[2][i],-W[j][i]);
	count++;
      }
    }
    sprintf(str+strlen(str), "\n");
  }


    //sprintf(str+strlen(str), "\nPrint[TimeConstrained[FindInstance[f[x1]==0&&f[x2]==0&&f[x3]==0&&T>0&&0<x1<x2<x3<1,{T,x1,x2,x3}],%d,999]]",timeout);

    sprintf(str+strlen(str), "\nPrint[TimeConstrained[FindInstance[");
    for(i=0;i<numeq;i++){
      if(i>0)
	sprintf(str+strlen(str), " && ");
      sprintf(str+strlen(str), " f1[x%d,y%d]==0 && f2[x%d,y%d]==0",i+1,i+1,i+1,i+1);
    }
    sprintf(str+strlen(str), " && T1>0 && T2>0");
    for(i=0;i<m;i++){
      for(j=0;j<numeq;j++){
	sprintf(str+strlen(str), "&& %d*x%d + %d*y%d + %d*(1-x%d-y%d)>0",basisT[0][i],j+1,basisT[1][i],j+1,basisT[2][i],j+1,j+1);
      }
    }
    for(i=0;i<numeq-1;i++){
      for(j=i+1;j<numeq;j++){
	sprintf(str+strlen(str), "&& (x%d!=x%d || y%d!=y%d)",i+1,j+1,i+1,j+1);
      }
    }
    //add not equal
    sprintf(str+strlen(str), ",{T1, T2");
    for(i=0;i<numeq;i++)
      sprintf(str+strlen(str), ", x%d, y%d",i+1, i+1);

    sprintf(str+strlen(str), "}],%d,999]]",timeout);



  fprintf(stderr, "%s\n\n",str);
  //fprintf(stderr, "totbasis=%d, totW=%d\n", totbasis, totW);
  //exit(0);
  free(str);
  free_imat(basisT,totbasis);
  free_imat(W, totW);
  free_imatrix(Sil1,0, n, 0, m-1);
  return NULL;
  //return str;
}

//overloading: input is PN adjacency matrix, also return the rank
char *multinn3nscr(int **AM, int n, int m, int numeq){
  bool minus=0;
  char *str;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  str=multinn3nscr(S, Sl, n, m, numeq);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return str;
}

//overloading: di6 input
char *multinn3nscr(char *di6, int n, int m, int numeq){
  char *str;
  matrix J;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=multinn3nscr(AM,n,m, numeq);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}


//multi342scr: Just for potentially nondegenerate (3,4,2)
//trial
char *multi342scr(int **Si, int **Sil, int n, int m){
  int i,j;
  int **basisT;
  int rk;
  int totbasis;
  char *str=(char *)malloc((size_t) (n*m*100*sizeof(char)));//should be safe
  int W[n];
  int **Sil1=imatrix(0, n, 0, m-1);
  int **Sil1t;
  int **Sit;
  possymbol x("x");
  char str1[5];
  ex eq;
  int Tpos;

  if(n!=3||m!=4){
    fprintf(stderr, "ERROR in \"multi342scr\": networks must have 3 species and four reactions.\n");
    exit(0);
  }

  //Augment with a row of ones
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      Sil1[i][j]=Sil[i][j];
    }
  }
  for(j=0;j<m;j++){
    Sil1[n][j]=1;
  }
  Sil1t=transposemat(Sil1,n+1,m);

  
  //printmat(Si,n,m);
  //printmat(Sil1t,m,n+1);
  if(matrank(Sil1t,m,n+1)==4){
    matrix G=imattoexmat(Sil1t,m,n+1);
    matrix G1=inverse(G);
    cout << G1 << endl;
  
    basisT=poskerbasis(Si,n,m,&totbasis,&rk);
    if(rk!=2){
      fprintf(stderr, "ERROR in \"multi342scr\": networks must have rank 2.\n");
      exit(0);
    }
    //printmat(basisT,totbasis,m); 
    Sit=transposemat(Si,n,m);
    intkervec(Sit, m, n, W,1);
    //printvec(W,n);
    matrix h(1,m);
    matrix hh(1,n);
    for(j=0;j<m;j++)
      h(0,j) = expand(basisT[0][j]*x+basisT[1][j]*(1-x));
    //cout << h << endl;
    for(j=0;j<n;j++){
      hh(0,j) = 1;
      for(i=0;i<m;i++){
	hh(0,j)*=pow(h(0,i),G1(j,i));
	//cout << h(0,i) << " " << G1(i,j) << " " << hh(0,j) << endl;
      }
      //expand(hh(0,j));
    }
    //cout << hh << endl;
    eq=0;Tpos=1;
    for(i=0;i<n;i++){
      if(W[i]!=0){
	if(W[i]<0)
	  Tpos=0;
	sprintf(str1, "t%d", i);
	eq+=W[i]*get_possymbol(str1)*hh(0,i);
      }
    }
    //expand(eq);
    fprintf(stderr, "f[x_] := ");
    for(i=0;i<n;i++){
      if(W[i]!=0){
	if(W[i]<0)
	  Tpos=0;
	fprintf(stderr, "%s", W[i]>0?"+":"-");
	if(abs(W[i])!=1)
	  fprintf(stderr, "%d*", abs(W[i]));
	fprintf(stderr, "t%d", i);
	for(j=0;j<m;j++){
	  if(G1(i,j)!=0 && h(0,j)!=1){
	    if(G1(i,j)==1)
	      cerr << "*(" << h(0,j) << ")";
	    else
	      cerr << "*(" << h(0,j) << ")^(" << G1(i,j) << ")";
	  }
	}
	
      }
    }
    cerr << endl;
    
    //cerr << eq << endl;
    fprintf(stderr, "Print[TimeConstrained[FindInstance[f[x1]==T && f[x2]==T && f[x3]==T && f[x4]==T &&  0 < x1 < x2 < x3 < x4 < 1 &&t0>0 && t1>0 && t2>0 && t3>0 %s, {T,x1,x2,x3,x4,t0,t1,t2,t3}],10,999]]\n", Tpos>0?"&& T>0":"");

    free_imatrix(Sit,0, m-1, 0, n-1);
    free_imat(basisT,totbasis);
  }

  free(str);
  free_imatrix(Sil1,0, n, 0, m-1);
  free_imatrix(Sil1t,0, m-1, 0, n);
  return NULL;
  //return str;
}

//overloading: input is PN adjacency matrix, also return the rank
char *multi342scr(int **AM, int n, int m){
  bool minus=0;
  char *str;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  str=multi342scr(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return str;
}

//overloading: di6 input
char *multi342scr(char *di6, int n, int m){
  char *str;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=multi342scr(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}

//print vector in kernel of [A|1], which is assumed to be 1-D
void printA1kervec(int **Si, int **Sil, int n, int m){
  int i,j;
  int W[n+1];
  int **Sil1t=imatrix(0, m-1, 0, n);

  //Augment with a column of ones
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      Sil1t[j][i]=Sil[i][j];
    }
  }
  for(j=0;j<m;j++){
    Sil1t[j][n]=1;
  }
  if(matrank(Sil1t,m,n+1)!=n){
    fprintf(stderr, "ERROR in \"printA1kervec\": [A|1] must have rank n.\n");
  }
  intkervec(Sil1t, m, n+1, W,1);

  printvec(W,n+1);
  free_imatrix(Sil1t,0, m-1, 0, n);
  return;
}

//overloading: input is PN adjacency matrix, also return the rank
void printA1kervec(int **AM, int n, int m){
  bool minus=0;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  printA1kervec(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return;
}

//overloading: di6 input
void printA1kervec(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  char *str=di6toreacstr((char*)di6, n, m, 0);
  fprintf(stderr, "Network:\n%s", str);
  free(str);
  printA1kervec(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}



//solvability condition assuming a single polynomial: e.g. (n,n+2,n) CRNs
//not identically degenerate CRNs. 
int printsolvability(int **Si, int **Sil, int n, int m, int debug){
  int degree=0;
  int tot;
  int rk;
  int **W=minA1tkerbasis(Sil,n,m,&tot,&rk,&degree,0,debug);
  printmat(W,tot,m);
  fprintf(stderr, "degree = %d\n",degree);
  free_imat(W,tot);
  return degree;
}


//overloading: input is PN adjacency matrix, also return the rank
int printsolvability(int **AM, int n, int m, int debug){
  bool minus=0,degree;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  degree=printsolvability(S, Sl, n, m, debug);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return degree;
}

//overloading: di6 input
int printsolvability(char *di6, int n, int m, int debug){
  int entries,degree,totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  char *str=di6toreacstr((char*)di6, n, m, 0);
  fprintf(stderr, "Network:\n%s", str);
  free(str);
  degree=printsolvability(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return degree;
}


//Bezout degree of solvability condition (assuming square)
//not identically degenerate CRNs. 
int solvabilitydegree(int **Si, int **Sil, int n, int m){
  int degree=0;
  int debug=0;
  int tot;
  int rk;
  int **W=minA1tkerbasis(Sil,n,m,&tot,&rk,&degree,0,debug);
  free_imat(W,tot);
  return degree;
}

//overloading: input is PN adjacency matrix, also return the rank
int solvabilitydegree(int **AM, int n, int m){
  bool minus=0;
  int **S, **Sl;
  int degree;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  degree=solvabilitydegree(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return degree;
}

//overloading: di6 input
int solvabilitydegree(char *di6, int n, int m){
  int entries;
  int totV;
  int degree;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  //char *str=di6toreacstr((char*)di6, n, m, 0);
  //fprintf(stderr, "Network:\n%s", str);
  //free(str);
  degree=solvabilitydegree(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return degree;
}


//Mixed volume of the Q-polynomials, assuming r = m-1
void printQmixed(int **Si, int **Sil, int n, int m){
  int i,degree=0;
  int debug=0;
  int tot,rk,V;
  int **Q=minA1tkerbasis(Sil,n,m,&tot,&rk,&degree,1,debug);
  //replace last col with zero col (constant on RHS of each equation)
  for(i=0;i<tot;i++)
    Q[i][n]=0;
  fprintf(stderr, "Q^t (last col replaced with zeros):\n");
  printmat(Q,tot,n+1);
  V=PolytopeVol(NULL,Q,tot,n+1,debug);
  fprintf(stderr, "Mixed Vol = %d\n", V);
  free_imat(Q,tot);
  return;
}

//overloading: input is PN adjacency matrix, also return the rank
void printQmixed(int **AM, int n, int m){
  bool minus=0;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  printQmixed(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return;
}

//overloading: di6 input
void printQmixed(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  //char *str=di6toreacstr((char*)di6, n, m, 0);
  //fprintf(stderr, "Network:\n%s", str);
  //free(str);
  printQmixed(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}



//Mixed volume of the Q-polynomials, assuming r = m-1
int Qmixed(int **Si, int **Sil, int n, int m){
  int i,degree=0;
  int debug=0;
  int tot,rk,V;
  int **Q=minA1tkerbasis(Sil,n,m,&tot,&rk,&degree,1,debug);
  for(i=0;i<tot;i++)//replace last col with zeros (constant on RHS)
    Q[i][n]=0;
  V=PolytopeVol(NULL,Q,tot,n+1,debug);
  free_imat(Q,tot);
  return V;
}

//overloading: input is PN adjacency matrix, also return the rank
int Qmixed(int **AM, int n, int m){
  bool minus=0;
  int **S, **Sl;
  int V;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  V=Qmixed(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return V;
}

//overloading: di6 input
int Qmixed(char *di6, int n, int m){
  int entries;
  int totV;
  int V;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  //char *str=di6toreacstr((char*)di6, n, m, 0);
  //fprintf(stderr, "Network:\n%s", str);
  //free(str);
  V=Qmixed(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return V;
}



//print the inverse of [A|1] assuming square and invertible
void printA1inv(int **Si, int **Sil, int n, int m){
  int i,j;
  int **Sil1=imatrix(0, n, 0, m-1);
  int **Sil1t;

  if(n!=m-1){
    fprintf(stderr, "ERROR in \"printA1inv\": networks must have one more reaction than species.\n");
    exit(0);
  }

  //Augment with a row of ones
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      Sil1[i][j]=Sil[i][j];
    }
  }
  for(j=0;j<m;j++){
    Sil1[n][j]=1;
  }
  Sil1t=transposemat(Sil1,n+1,m);

  
  //printmat(Si,n,m);
  //printmat(Sil1t,m,n+1);
  if(matrank(Sil1t,m,n+1)!=m){
    fprintf(stderr, "ERROR in \"printA1inv\": [A|1] must have full rank.\n");
  }
    
  matrix G=imattoexmat(Sil1t,m,n+1);
  matrix G1=inverse(G);
  cerr << G1 << endl;
  
  free_imatrix(Sil1,0, n, 0, m-1);
  free_imatrix(Sil1t,0, m-1, 0, n);
  return;
  //return str;
}

//overloading: input is PN adjacency matrix, also return the rank
void printA1inv(int **AM, int n, int m){
  bool minus=0;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  printA1inv(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return;
}

//overloading: di6 input
void printA1inv(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  char *str=di6toreacstr((char*)di6, n, m, 0);
  fprintf(stderr, "Network:\n%s", str);
  free(str);
  printA1inv(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}


//CylinderScr: output for mathematica
//Searching for multiple positive equilibria
char *CylinderScr3(int **Si, int **Sil, int n, int m){
  int i,i1,j,t,count=0;
  int numeqs=3;//number of equilibria
  int timeout=10;//maximum seconds
  char *str=(char *)malloc((size_t) (n*m*100*sizeof(char)));//should be safe
  //printmat(Si,n,m);
  //printmat(Sil,n,m);

  sprintf(str, "Print[TimeConstrained[CylindricalDecomposition[");
  //main part
  for(t=0;t<numeqs;t++){//eq number
    for(i=0;i<n;i++){//species
      if(count>0)
	sprintf(str+strlen(str), " && ");
      for(j=0;j<m;j++){//reaction
	if(Si[i][j]!=0){
	  if(Si[i][j]>0&&j>0)
	    sprintf(str+strlen(str), "+");
	  if(Si[i][j]==1)
	    sprintf(str+strlen(str), "k%d",j+1);
	  else if(Si[i][j]==-1)
	    sprintf(str+strlen(str), "-k%d",j+1);
	  else
	    sprintf(str+strlen(str), "%d*k%d",Si[i][j],j+1);
	  for(i1=0;i1<n;i1++){
	    if(Sil[i1][j]>0){
	      if(Sil[i1][j]>1)
		sprintf(str+strlen(str), "*x%d%d^%d",i1+1,t+1,Sil[i1][j]);
	      else
		sprintf(str+strlen(str), "*x%d%d",i1+1,t+1);

	    }
	  }
	}
      }
      sprintf(str+strlen(str), " == 0");
      count++;
    }
  }
  for(j=0;j<m;j++)//reaction
    sprintf(str+strlen(str), " && k%d>0",j+1);

  for(j=0;j<n;j++){//species
    for(t=0;t<numeqs;t++){//eq no
      sprintf(str+strlen(str), " && x%d%d>0",j+1,t+1);
    }
  }

  for(t=0;t<numeqs;t++){//eq no
    for(i1=t+1;i1<numeqs;i1++){//eq no
      sprintf(str+strlen(str), " && (");
      for(j=0;j<n;j++){//species
	if(j>0)
	  sprintf(str+strlen(str), " || x%d%d!=x%d%d",j+1,t+1,j+1,i1+1);
	else
	  sprintf(str+strlen(str), "x%d%d!=x%d%d",j+1,t+1,j+1,i1+1);
      }
      sprintf(str+strlen(str), ")");
    }
  }

  sprintf(str+strlen(str), ", {");
  for(j=0;j<m;j++){//reaction
    if(j>0)
      sprintf(str+strlen(str), ", k%d",j+1);
    else
      sprintf(str+strlen(str), "k%d",j+1);
  }
  for(i1=0;i1<numeqs;i1++){//eq no
    for(j=0;j<n;j++){//species
      sprintf(str+strlen(str), ", x%d%d",j+1,i1+1);
    }
  }
  sprintf(str+strlen(str), "}],%d,999]]",timeout);


  fprintf(stderr, "%s\n",str);
  //exit(0);
  free(str);
  return NULL;
  //return str;
}

//overloading: input is PN adjacency matrix, also return the rank
char *CylinderScr3(int **AM, int n, int m){
  bool minus=0;
  char *str;
  int **S, **Sl;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  str=CylinderScr3(S, Sl, n, m);
  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(S,0,n-1,0,m-1);
  return str;
}

//overloading: di6 input
char *CylinderScr3(char *di6, int n, int m){
  char *str;
  matrix J;
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=CylinderScr3(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}


// Take a digraph Nauty digraph6 format (two layers) and 
// convert to a sign pattern 


int **di6tosignpat(char *di6, int n){
  int j,k;
  int *out2;
  int **V=imatrix(0, n-1, 0, n-1);
  int p,p1;
  int n1,partner[n];

  out2=digraph6toam(di6, &n1);
  if(n1!=2*n){
    fprintf(stderr, "Wrong dimensions in di6tosignpat. Exiting. n1=%d, n = %d\n", n1, n);
    exit(0);
  }
  //find partner in upper level: only needed if labelg was used
  for(j=0;j<n;j++){
    for(k=0;k<n;k++){
      if(out2[2*n*j+k+n]){//partner in next level found
	partner[j]=k+n-j;
	break;
      }
    }
  }

  for(j=0;j<n;j++){
    for(k=0;k<n;k++){
      p=out2[2*n*j+k];p1=out2[2*n*(j+partner[j])+k+partner[k]];//
      if(p){V[j][k]=1;}else if(p1){V[j][k]=-1;}else{V[j][k]=0;}
    }
  }
  free((char*)out2);
  //printmat(V,n,n);

  return V;

}

//di6 to Sauro format
// memory allocation assumes no integers more than 2 digits long
char *di6toSauro(char *di6, int n, int m){
  int i,j,k;
  int entries;
  int totV;
  int **V=di6toCRNam1(di6, n, m, &totV, &entries);
  char *str;

  str=(char *)malloc((size_t) ((8*entries+20)*sizeof(char)));
  str[0]=0;

  sprintf(str+strlen(str), "%d %d ", m, n);

  for(k=0;k<m;k++){//k=reaction number
    for(j=0;j<n;j++){//j=species number
      if(V[j][k+n]){//species to reaction
	for(i=0;i<V[j][k+n];i++){
	  sprintf(str+strlen(str), "%d %d ",j+m,k);
	}
      }
    }
    for(j=0;j<n;j++){
      if(V[k+n][j]){//reaction to species
	for(i=0;i<V[k+n][j];i++){
	  sprintf(str+strlen(str), "%d %d ",k,j+m);
	}
      }
    }
  }
  free_imatrix(V, 0, n+m-1, 0, n+m-1);
  return str;
}


// convert a file of di6 reaction digraphs (two-layer)
// to a file of CRNs (human readable)
// We need to know n=#species and m=#reacs as
// the di6 string only stores the total #vertices
// the switch "open" means add the reactions <-->Ai
unsigned long di6toSaurofile(const char *infile, const char *outfile, int n, int m){
  int linelen;
  char *line;
  unsigned long numl=0, tot=0;
  int maxl=0;
  FILE *fd, *fd1;
  char *str, *str1;


  if((numl=numflines(infile, &maxl))<=0){//no reactions
    if(!(fd1=fopen(outfile, "w"))){
	fprintf(stderr, "ERROR in di6toSaurofile - output file \"%s\" could not be opened...\n", outfile);
	exit(0);
    }
    fclose(fd1);
    return 0;
  }
  line = (char*) malloc(sizeof(char) * (maxl));

  if(!(fd=fopen(infile, "r"))){
    fprintf(stderr, "ERROR in di6toSaurofile - input file \"%s\" could not be opened...\n", infile);
    exit(0);
  }
  //outfile=NULL means print to stderr
  if(outfile && !(fd1=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in di6toSaurofile - output file: \"%s\" could not be opened...\n", outfile);
    exit(0);
  }

  while((linelen = gtline(fd, line, maxl)) > 0){
    if(strcmp(line,"") && !iscomline(line)){
      if((tot+1)%100000==0)
	fprintf(stderr, "*");
      str=getnthwd(line, 1);
      str1=di6toSauro(str,n,m);
      if(outfile)
	fprintf(fd1, "%s\n", str1);
      else
	fprintf(stderr, "%s\n", str1);
      free(str);free(str1);
      tot++;
    }
  }
  fclose(fd);
  if(outfile)
    fclose(fd1);
  free(line);
  return numl;
}



// no error checking on length of outstr
void di6toSauro2(char *di6, int n, int m, char *outstr){
  int i,j,k;
  int entries;
  int totV;
  int **V=di6toCRNam1(di6, n, m, &totV, &entries);
  char tmpstr[30];

  sprintf(tmpstr, "%d %d ", m, n);
  strcpy(outstr, tmpstr);

  for(k=0;k<m;k++){//k=reaction number
    for(j=0;j<n;j++){//j=species number
      if(V[j][k+n]){//species to reaction
	for(i=0;i<V[j][k+n];i++){
	  sprintf(tmpstr, "%d %d ",j+m,k);
	  strcat(outstr, tmpstr);
	}
      }
    }
    for(j=0;j<n;j++){
      if(V[k+n][j]){//reaction to species
	for(i=0;i<V[k+n][j];i++){
	  sprintf(tmpstr, "%d %d ",k,j+m);
	  strcat(outstr, tmpstr);
	}
      }
    }
  }
  free_imatrix(V, 0, n+m-1, 0, n+m-1);
  return;

}

void di6toSauro1(FILE *fd, char *di6, int n, int m){
  int i,j,k;
  int entries;
  int totV;
  int **V=di6toCRNam1(di6, n, m, &totV, &entries);

  fprintf(fd, "%d %d ", m, n);

  for(k=0;k<m;k++){//k=reaction number
    for(j=0;j<n;j++){//j=species number
      if(V[j][k+n]){//species to reaction
	for(i=0;i<V[j][k+n];i++){
	  fprintf(fd, "%d %d ",j+m,k);
	}
      }
    }
    for(j=0;j<n;j++){
      if(V[k+n][j]){//reaction to species
	for(i=0;i<V[k+n][j];i++){
	  fprintf(fd, "%d %d ",k,j+m);
	}
      }
    }
  }
  fprintf(fd, "\n");
  free_imatrix(V, 0, n+m-1, 0, n+m-1);
  return;

}

// convert a file of di6 reaction digraphs (two-layer)
// to a file of CRNs (human readable)
// We need to know n=#species and m=#reacs as
// the di6 string only stores the total #vertices
// the switch "open" means add the reactions <-->Ai

unsigned long di6toreacfile(const char *infile, const char *outfile, int n, int m, int open){
  int linelen;
  char *line;
  unsigned long numl=0;
  int maxl=0;
  FILE *fd, *fd1;
  char *str, *str1;
  unsigned int tot=0;

  if((numl=numflines(infile, &maxl))<=0){//no reactions
    if(!(fd1=fopen(outfile, "w"))){
	fprintf(stderr, "ERROR in di6toreacfile - File: %s could not be opened...\n", outfile);
	exit(0);
    }
    fclose(fd1);
    return 0;
  }
  line = (char*) malloc(sizeof(char) * (maxl));

  if(!(fd=fopen(infile, "r"))){
    fprintf(stderr, "ERROR in di6toreacfile - File: %s could not be opened...\n", infile);
    exit(0);
  }
  //outfile=NULL means print to stderr
  if(outfile && !(fd1=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in di6toreacfile - File: %s could not be opened...\n", outfile);
    exit(0);
  }

  while((linelen = gtline(fd, line, maxl)) > 0){
    if((tot+1)%100000==0)
      fprintf(stderr, "*");
    str=getnthwd(line, 1);
    str1=di6toreacstr(str,n,m,open);
    if(outfile)
      fprintf(fd1, "//%d/%ld\n%s******\n", tot+1, numl, str1);
    else
      fprintf(stderr, "%s******\n", str1);
    tot++;
    free(str);free(str1);
  }
  fclose(fd);
  if(outfile)
    fclose(fd1);
  free(line);
  return numl;
}




//Visualise Petri Net of CRN from the adjacency matrix
//In dot language
char *printPNgraph(int **AM, int n, int m, unsigned long labl){
  int i,j;
  char *str=(char*) malloc(sizeof(char) * (30*n*m+50*(n+m)+100));
  str[0]=0;

  sprintf(str+strlen(str), "digraph CRN%ld {\n", labl);
  sprintf(str+strlen(str), "\tnode[label=\"\"];\n");
  for(i=0;i<n;i++)//species
    sprintf(str+strlen(str), "\tS%d [width=0.2] [height=0.2] [color=blue]\n",i+1);
  for(i=0;i<m;i++)//reactions
    sprintf(str+strlen(str), "\tR%d [width=0.2] [height=0.2]\n",i+1);
  for(i=0;i<n;i++){//species
    for(j=0;j<m;j++){//reactions
	//sprintf(str+strlen(str), "\tS%d -> R%d [label=%d]\n",i+1,j+1,AM[i][n+j]);
      if(AM[i][n+j]>1)
	sprintf(str+strlen(str), "\tS%d -> R%d [color=red]\n",i+1,j+1);
      else if(AM[i][n+j])	
	sprintf(str+strlen(str), "\tS%d -> R%d\n",i+1,j+1);


	//sprintf(str+strlen(str), "\tR%d -> S%d [label=%d]\n",j+1,i+1,AM[n+j][i]);
      if(AM[n+j][i]>1)
	sprintf(str+strlen(str), "\tR%d -> S%d [color=red]\n",j+1,i+1);
      else if(AM[n+j][i])
	sprintf(str+strlen(str), "\tR%d -> S%d\n",j+1,i+1);
    }
  }
  sprintf(str+strlen(str), "}\n######\n");

  return str;
}

//overloading to take di6 input
char *printPNgraph(char *di6, int n, int m, unsigned long labl){
  int entries;
  char *str;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  str=printPNgraph(AM, n, m, labl);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return str;
}


// convert a file of di6 reaction digraphs (two-layer)
// to a file of PN graphs in dot language
unsigned long di6todotfile(const char *infile, const char *outfile, int n, int m){
  int linelen;
  char *line;
  unsigned long numl=0;
  int maxl=0;
  FILE *fd, *fd1;
  char *str, *str1;

  if((numl=numflines(infile, &maxl))<=0){//no reactions
    if(!(fd1=fopen(outfile, "w"))){
	fprintf(stderr, "ERROR in di6todotfile - File: %s could not be opened...\n", outfile);
	exit(0);
    }
    fclose(fd1);
    return 0;
  }
  line = (char*) malloc(sizeof(char) * (maxl));

  if(!(fd=fopen(infile, "r"))){
    fprintf(stderr, "ERROR in di6todotfile - File: %s could not be opened...\n", infile);
    exit(0);
  }
  //outfile=NULL means print to stderr
  if(outfile && !(fd1=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in di6todotfile - File: %s could not be opened...\n", outfile);
    exit(0);
  }

  numl=0;
  while((linelen = gtline(fd, line, maxl)) > 0){
    numl++;
    str=getnthwd(line, 1);
    str1=printPNgraph(str,n,m,numl);
    if(outfile)
      fprintf(fd1, "%s\n", str1);
    else
      fprintf(stderr, "%s\n", str1);
    free(str);free(str1);
  }
  fclose(fd);
  if(outfile)
    fclose(fd1);
  free(line);
  return numl;
}


//
// Filter crn lists
//

// return the rank of a CRN from its left and right stoichiometric matrices
int checkrankLR(int **L, int **R, int n, int m){
  int ret;
  int **S=subtract(R,L,n,m);
  matrix exS=imattoexmat(S,n,m);//in order to use the rank() function in GiNAC
  ret=exS.rank();
  free_imatrix(S, 0, n-1, 0, m-1);
  return ret;
}

// return the rank of a CRN from its stoichiometric matrix
int checkrankS(int **S, int n, int m){
  matrix exS=imattoexmat(S,n,m);//in order to use the rank() function in GiNAC
  return exS.rank();
}

// return the rank of a CRN from its adjacency matrix
int checkrank(int **AM, int n, int m){
  int ret;
  int **S=CRNamtostoichmat(AM,n,m);//stoichiometric matrix
  ret=checkrankS(S,n,m);
  free_imatrix(S, 0, n-1, 0, m-1);
  return ret;
}



//l layers. Overloading: return the rank of a CRN in di6 format
int checkrank(char *di6, int n, int m){
  int entries,ret;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=checkrank(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//2 layers, deprecated
int checkrankold(char *di6, int n, int m){
  int entries,ret;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  ret=checkrank(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


// return the rank of the reactant matrix from CRN adjacency matrix
int checkrankleft(int **AM, int n, int m){
  int ret;
  int **S=CRNamtoleftstoichmat(AM,n,m);//left stoichiometric matrix
  ret=checkrankS(S,n,m);
  free_imatrix(S, 0, n-1, 0, m-1);
  return ret;
}

// Overloading: return the rank the reactant matrix of a CRN in di6 format
int checkrankleft(char *di6, int n, int m){
  int entries,ret;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=checkrankleft(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


// return the rank of the reactant matrix from CRN adjacency matrix
int checkrankSil1(int **AM, int n, int m){
  int ret;
  int **S=CRNamtoSil1(AM,n,m);//left stoichiometric matrix with added row of ones
  ret=checkrankS(S,n+1,m);
  free_imatrix(S, 0, n, 0, m-1);
  return ret;
}

// Overloading: return the rank the reactant matrix of a CRN in di6 format
int checkrankSil1(char *di6, int n, int m){
  int entries,ret;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=checkrankSil1(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



//depth first seach on the Petri Net graph, stored as
//left and right stoich matrices. Recursive.
void dfs(int k, int **Sl, int **Sr, int n, int m, int *v, int *tot){
  int j;
  if(v[k])// already visited
    return;
  (*tot)++;
  v[k]=1;
  if(k>=n){// r-vertex: column k-n
    for(j=0;j<n;j++){
      if(Sr[j][k-n] || Sl[j][k-n]){// omit Sl for strongly connected
	if(!(v[j]))
	  dfs(j,Sl,Sr,n,m,v,tot);// follow S-vertex j
      }
    }
  }
  else{//s-vertex
    for(j=0;j<m;j++){
      if(Sl[k][j] || Sr[k][j]){ // omit Sr for strongly connected
	if(!(v[n+j]))
	  dfs(n+j,Sl,Sr,n,m,v,tot);// follow R-vertex j
      }
    }
  }
  return;
}


// Is a PetriNet graph connected?
// Sl and Sr are the left and right stoichiometric matrices
bool connected(int **Sl, int **Sr, int n, int m){
  int v[n+m];// visited
  int tot=0;
  inittozero(v,n+m);
  dfs(0,Sl,Sr,n,m,v, &tot);
  if(tot==n+m)
    return 1;
  return 0;
}

// Is a PetriNet graph connected?
int PNconnect(int **Sl, int **Sr, int n, int m){
  return connected(Sl, Sr, n, m);
}

//overloading: PN AM input
int PNconnect(int **AM, int n, int m){
  int i,j,flg=0;
  int **Sl=imatrix(0,n-1,0,m-1);
  int **Sr=imatrix(0,n-1,0,m-1);

  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      Sr[j][i]=AM[i+n][j];Sl[j][i]=AM[j][i+n];
    }
  }
  flg=PNconnect(Sl, Sr, n, m);

  free_imatrix(Sl,0,n-1,0,m-1);
  free_imatrix(Sr,0,n-1,0,m-1);
  return flg;
}


// Overloading: di6 input
int PNconnect(char *di6, int n, int m){
  int entries,flg=0;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  flg=PNconnect(AM, n, m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return flg;
}

// Is a PetriNet graph strongly connected? (Probably inefficient:
// No need to find all the SCCs)
int PNstrongconnect(int **AM, int n, int m){
  int comps=Tarjan_light(AM, n+m);
  if(comps==1)//1 SCC
    return 1;
  return 0;
}


int PNstrongconnect(char *di6, int n, int m){
  int entries,flg=0;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  flg=PNstrongconnect(AM, n, m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return flg;
}



// Is the CRN dynamically nontrivial? Input: Petri Net adjacency matrix
bool DN(int **AM, int n, int m){
  bool ret;
  int i,j;
  int **St=imatrix(0,m-1,0,n-1);//transposed stoich matrix
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      St[i][j]=AM[i+n][j]-AM[j][i+n];
    }
  }
  ret=hasposimvec(St,m,n);
  free_imatrix(St, 0, m-1, 0, n-1);
  return 1-ret;
}

//l layers
bool DN(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  bool ret=DN(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//2 layers, deprecated
bool DNold(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=DN(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// Bounded stoichiometry classes?
bool bdclass(int **AM, int n, int m, int debug){
  bool ret;
  int i,j;
  int **S=imatrix(0,n-1,0,m-1);//stoich matrix
  if(debug){fprintf(stderr, "\n###Entering bdclass: checking if stoichiometric classes are bounded.\n");}
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      S[j][i]=AM[i+n][j]-AM[j][i+n];
    }
  }
  ret=hasposimvec(S,n,m);
  free_imatrix(S, 0, n-1, 0, m-1);
  if(debug){fprintf(stderr, "Result: stoichiometric classes are %s.\n",ret?"unbounded":"bounded");}
  return 1-ret;
}

bool bdclass(char *di6, int n, int m, int debug){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=bdclass(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// Homogeneous: i.e., 1^t\Gamma = 0
bool homogen(int **AM, int n, int m, int debug){
  int i,j,tot;
  int **S=imatrix(0,n-1,0,m-1);//stoich matrix
  if(debug){fprintf(stderr, "\n###Entering homogen: checking if the stoichiometric subspace has a vector of ones in its left kernel.\n");}
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      S[j][i]=AM[i+n][j]-AM[j][i+n];
    }
  }
  for(i=0;i<m;i++){// reactions
    tot=0;
    for(j=0;j<n;j++)// species
      tot+=S[j][i];
    if(tot){//not homogeneous
      free_imatrix(S, 0, n-1, 0, m-1);
      if(debug){fprintf(stderr, "Result: not homogeneous.\n");}
      return 0;
    }
  }
  free_imatrix(S, 0, n-1, 0, m-1);
  if(debug){fprintf(stderr, "Result: homogeneous.\n");}
  return 1;
}

bool homogen(char *di6, int n, int m, int debug){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=homogen(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


// Check if the system is weakly reversible
// We do this by computing the complex graph, 
// computing all its SCCs and all its CCs
// and checking that each CC is an SCC
// additional info such as deficiency not computed

bool CRNweakrev(int **AM, int n, int m){
  int i,j,totcmplx=0;
  int **Slt, **Srt; //transposed left and right stoichiometry matrices
  int **cmplxs=NULL;
  int ind1,ind2;
  int **cmpmat; //for SCCs (the digraph)
  int **cmpmat1; //for CCs (the graph)
  int edg[m][2];
  int **SCC=NULL, **CC=NULL;
  int totSCC=0,totCC=0;

  Slt=imatrix(0,m-1,0,n-1);Srt=imatrix(0,m-1,0,n-1);
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      Srt[i][j]=AM[i+n][j];Slt[i][j]=AM[j][i+n];
    }
  }

  for(i=0;i<m;i++){// get the complexes and the edges of the incidence graph
    totcmplx=addintstr(totcmplx, Slt[i], n, &cmplxs, &ind1);
    totcmplx=addintstr(totcmplx, Srt[i], n, &cmplxs, &ind2);
    edg[i][0]=ind1;edg[i][1]=ind2;
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(Srt,0,m-1,0,n-1);

  //Create the adjacency matrix of the complex digraph (cmpmat) 
  //and its symmetrised version (cmpmat1)
  cmpmat=imatrix(0,totcmplx-1,0,totcmplx-1);
  cmpmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  for(i=0;i<totcmplx;i++){
    for(j=0;j<totcmplx;j++){
      cmpmat[i][j]=0;
      cmpmat1[i][j]=0;
    }
  }
  for(i=0;i<m;i++){
    cmpmat[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][1]][edg[i][0]]=1;
  }

  // Extract the SCCs and CCs using Tarjan's algorithm
  // First element of a CC or SCC is its size
  // Then the list of complexes
  SCC=Tarjan(cmpmat, totcmplx, &totSCC);
  CC=Tarjan(cmpmat1, totcmplx, &totCC);


  free_imat(SCC,totSCC);
  free_imat(CC,totCC);
  free_imat(cmplxs,totcmplx);

  free_imatrix(cmpmat,0,totcmplx-1,0,totcmplx-1);
  free_imatrix(cmpmat1,0,totcmplx-1,0,totcmplx-1);

  if(totSCC==totCC)
    return 1;
  return 0;

}

//overloading
bool CRNweakrev(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=CRNweakrev(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//The version of the weak reversibility check where reactions are stored
//as LHS,RHS pairs, and xc holds the indices of the reactions
//for use within "genbimol", for example
bool CRNweakrev2(int Rvec[][2], int *xc, int m){
  int i,j,totcmplx=0;
  int **cmat; //for SCCs (the digraph) 
  int **cmat1; //for CCs (the graph)
  int edg[m][2];
  int totSCC=0,totCC=0;
  int cmplx[2*m];

  //get complexes and edges of complex digraph
  for(i=0;i<m;i++){
    edg[i][0]=addtovec(cmplx,&totcmplx,Rvec[xc[i]][0]);
    edg[i][1]=addtovec(cmplx,&totcmplx,Rvec[xc[i]][1]);
  }

  //make the adjacency matrices of complex digraph and graph
  cmat=imatrix(0,totcmplx-1,0,totcmplx-1);
  cmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  for(i=0;i<totcmplx;i++){
    for(j=0;j<totcmplx;j++){
      cmat[i][j]=0;
      cmat1[i][j]=0;
    }
  }

  for(i=0;i<m;i++){
    cmat[edg[i][0]][edg[i][1]]=1;
    cmat1[edg[i][0]][edg[i][1]]=1;
    cmat1[edg[i][1]][edg[i][0]]=1;
  }

  // Extract the SCCs and CCs using Tarjan's algorithm
  // First element of a CC or SCC is its size
  // Then the list of complexes
  totSCC = Tarjan_light(cmat, totcmplx);
  totCC = Tarjan_light(cmat1, totcmplx);

  free_imatrix(cmat,0,totcmplx-1,0,totcmplx-1);
  free_imatrix(cmat1,0,totcmplx-1,0,totcmplx-1);

  if(totSCC==totCC)
    return 1;
  return 0;

}


//Deficiency of the CRN
int deficiency(int **AM, int n, int m){
  int i,j,totcmplx=0;
  int **Slt, **Srt; //transposed left and right stoichiometry matrices
  int **cmplxs=NULL;
  int ind1,ind2;
  int **cmpmat1; //for CCs (the graph)
  int edg[m][2];
  int totCC=0;
  int **S;//stoichiometric matrix (not transposed)
  matrix exS;
  int Srank;

  Slt=imatrix(0,m-1,0,n-1);Srt=imatrix(0,m-1,0,n-1);S=imatrix(0,m-1,0,n-1);
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      Srt[i][j]=AM[i+n][j];Slt[i][j]=AM[j][i+n];
      S[j][i]=Srt[i][j]-Slt[i][j];//RHS stoic-LHS stoic
    }
  }

  //rank
  exS=imattoexmat(S,n,m);//in order to use the rank() function in GiNAC
  Srank=exS.rank();
  free_imatrix(S, 0, n-1, 0, m-1);

  for(i=0;i<m;i++){// get the complexes and the edges of the incidence graph
    totcmplx=addintstr(totcmplx, Slt[i], n, &cmplxs, &ind1);
    totcmplx=addintstr(totcmplx, Srt[i], n, &cmplxs, &ind2);
    edg[i][0]=ind1;edg[i][1]=ind2;
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(Srt,0,m-1,0,n-1);

  //Create the symmetrised adjacency matrix of the complex digraph (cmpmat1)
  cmpmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  inittozero(cmpmat1,totcmplx,totcmplx);
  for(i=0;i<m;i++){
    cmpmat1[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][1]][edg[i][0]]=1;
  }

  totCC=Tarjan_light(cmpmat1, totcmplx);

  free_imat(cmplxs,totcmplx);

  free_imatrix(cmpmat1,0,totcmplx-1,0,totcmplx-1);

  return totcmplx-totCC-Srank;//deficiency

}

//overloading
int deficiency(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  int ret=deficiency(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Number of linkage classes of the CRN
int linkage_classes(int **AM, int n, int m){
  int i,j,totcmplx=0;
  int **Slt, **Srt; //transposed left and right stoichiometry matrices
  int **cmplxs=NULL;
  int ind1,ind2;
  int **cmpmat1; //for CCs (the graph)
  int edg[m][2];
  int totCC=0;

  Slt=imatrix(0,m-1,0,n-1);Srt=imatrix(0,m-1,0,n-1);
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      Srt[i][j]=AM[i+n][j];Slt[i][j]=AM[j][i+n];
    }
  }

  for(i=0;i<m;i++){// get the complexes and the edges of the incidence graph
    totcmplx=addintstr(totcmplx, Slt[i], n, &cmplxs, &ind1);
    totcmplx=addintstr(totcmplx, Srt[i], n, &cmplxs, &ind2);
    edg[i][0]=ind1;edg[i][1]=ind2;
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(Srt,0,m-1,0,n-1);

  //Create the symmetrised adjacency matrix of the complex digraph (cmpmat1) 
  //and its symmetrised version (cmpmat1)
  cmpmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  inittozero(cmpmat1,totcmplx,totcmplx);
  for(i=0;i<m;i++){
    cmpmat1[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][1]][edg[i][0]]=1;
  }

  //Connected components
  totCC=Tarjan_light(cmpmat1, totcmplx);

  free_imat(cmplxs,totcmplx);
  free_imatrix(cmpmat1,0,totcmplx-1,0,totcmplx-1);

  return totCC;
}

//Overloading di6 input
int linkage_classes(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  int ret=linkage_classes(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



// Weakly reversible and deficiency zero
// A lot of redundancy if we separately do the two checks
// Takes as input the adjacency matrix of the PN graph

bool WRdef0(int **AM, int n, int m){
  int i,j,totcmplx=0;
  int **Slt, **Srt; //transposed left and right stoichiometry matrices
  int **cmplxs=NULL;
  int ind1,ind2;
  int **cmpmat; //for SCCs (the digraph)
  int **cmpmat1; //for CCs (the graph)
  int edg[m][2];
  int **SCC=NULL, **CC=NULL;
  int totSCC=0,totCC=0;
  int **S;//stoichiometric matrix (not transposed)
  matrix exS;
  int Srank,def;

  Slt=imatrix(0,m-1,0,n-1);Srt=imatrix(0,m-1,0,n-1);S=imatrix(0,m-1,0,n-1);
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      Srt[i][j]=AM[i+n][j];Slt[i][j]=AM[j][i+n];
      S[j][i]=Srt[i][j]-Slt[i][j];//RHS stoic-LHS stoic
    }
  }
  //rank
  exS=imattoexmat(S,n,m);//in order to use the rank() function in GiNAC
  Srank=exS.rank();
  free_imatrix(S, 0, n-1, 0, m-1);

  for(i=0;i<m;i++){// get the complexes and the edges of the incidence graph
    totcmplx=addintstr(totcmplx, Slt[i], n, &cmplxs, &ind1);
    totcmplx=addintstr(totcmplx, Srt[i], n, &cmplxs, &ind2);
    edg[i][0]=ind1;edg[i][1]=ind2;
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(Srt,0,m-1,0,n-1);

  //Create the adjacency matrix of the complex digraph (cmpmat) 
  //and its symmetrised version (cmpmat1)
  cmpmat=imatrix(0,totcmplx-1,0,totcmplx-1);
  cmpmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  for(i=0;i<totcmplx;i++){
    for(j=0;j<totcmplx;j++){
      cmpmat[i][j]=0;
      cmpmat1[i][j]=0;
    }
  }
  for(i=0;i<m;i++){
    cmpmat[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][1]][edg[i][0]]=1;
  }

  // Extract the SCCs and CCs using Tarjan's algorithm
  // First element of a CC or SCC is its size
  // Then the list of complexes
  SCC=Tarjan(cmpmat, totcmplx, &totSCC);
  CC=Tarjan(cmpmat1, totcmplx, &totCC);

  free_imat(SCC,totSCC);
  free_imat(CC,totCC);
  free_imat(cmplxs,totcmplx);

  free_imatrix(cmpmat,0,totcmplx-1,0,totcmplx-1);
  free_imatrix(cmpmat1,0,totcmplx-1,0,totcmplx-1);

  def=totcmplx-totCC-Srank;//deficiency
  if(totSCC==totCC && def==0)
    return 1;
  return 0;

}

//overloading
int WRdef0(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=WRdef0(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// genuine if there are no unused species
// returns 1 if every species participates in some reaction
// otherwise returns 0

bool genuine(int **AM, int n, int m){
  int flag,i,j;
  for(j=0;j<n;j++){// each species
    flag=0;//potentially an unused species
    for(i=0;i<m;i++){// reactions
      if(AM[i+n][j] || AM[j][i+n]){//not unused
	flag=1;
	break;
      }
    }
    if(!flag)//unused species
      return 0;
  }
  return 1;//genuine
}

//l layers
bool genuine(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  bool ret=genuine(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//2 layers only
bool genuineold(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=genuine(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



// check if every species has some net production or consumption in some reac
// returns 0 if this is the case, otherwise returns 1


bool triv_spec(int **AM, int n, int m){
  int flag,i,j;
  for(j=0;j<n;j++){// each species
    flag=0;//potentially a trivial species
    for(i=0;i<m;i++){// reactions
      if(AM[i+n][j]!=AM[j][i+n]){//net production or consumption
	flag=1;
	break;
      }
    }
    if(!flag)//trivial species
      return 1;
  }
  return 0;//no trivial species
}

//overloading
bool triv_spec(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=triv_spec(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// zero-one CRN: no stoichiometry larger than one
bool zeroone(int **AM, int n, int m){
  int i,j;
  for(j=0;j<n;j++){// each species
    for(i=0;i<m;i++){// reactions
      if(AM[i+n][j]>1 || AM[j][i+n]>1){//not zero-one
	return 0;
      }
    }
  }
  return 1;//zero-one
}

//overloading
bool zeroone(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=zeroone(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


// nonautocatalytic if no species ever figures on both sides of any
// reaction
bool nonautocat(int **AM, int n, int m){
  int i,j;
  for(j=0;j<n;j++){// each species
    for(i=0;i<m;i++){// reactions
      if(AM[i+n][j] && AM[j][i+n]){//autocatalytic
	return 0;
      }
    }
  }
  return 1;//not autocatalytic
}

//overloading
bool nonautocat(char *di6, int n, int m){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  bool ret=nonautocat(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Is the Jacobian matrix P0?
//The version which uses symbolic algebra
//slower at least for smallish matrices
int reacJisP01(int **AM, int n, int m, bool minus, int maxpppdeg, int debug){
  int numv, ret;
  matrix J=reacJac(AM,n,m,&numv,minus);
  int pppdegused;
  ret=isP0matorth(J,n,0,maxpppdeg,&pppdegused,debug);
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}


//overloading
int reacJisP01(char *di6, int n, int m, bool minus, int maxpppdeg, int debug){
  int numv, ret;
  matrix J=reacJac(di6,n,m,&numv,minus);
  int pppdegused;
  ret=isP0matorth(J,n,0,maxpppdeg,&pppdegused,debug);
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}

//Is the Jacobian matrix P0/P0-? (P0- = Accordance)
//The version which uses compatibility checking
//faster, at least for small matrices
//Here AM is the Petri Net AM (not DSR AM)
int reacJisP0(int **AM, int n, int m, bool minus, int debug){
  int **Sl, **S;
  int ret,strict;
  int debugfull=(debug<=0)?0:debug-1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  //important: if Sl is not a pattern, then arecompatsimp fails
  //So must change Sl to a sign matrix
  sgn_mat(Sl,n,m);
  ret=arecompatsimp(S,Sl,n,m,&strict,debugfull);
  free_imatrix(S,0,n-1,0,m-1);
  free_imatrix(Sl,0,n-1,0,m-1);
  if(ret==1){//can also be -1, but this is not accordance
    if(debug){fprintf(stderr, "The network is accordant (The Jacobian matrix is a P0- matrix).\n\n");}
    return 1;
  }

  if(debug){fprintf(stderr, "The network fails to be accordant.\n");}
  return 0;
}

//overloading
int reacJisP0(char *di6, int n, int m, bool minus, int debug){
  int entries;
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  int ret=reacJisP0(AM,n,m,minus,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

int accord(int **AM, int n, int m, int debug){
  if(debug){fprintf(stderr, "\n###Entering accord: checking accordance.\n");}
  return reacJisP0(AM,n,m,1,debug);
}

int accord(char *di6, int n, int m, int debug){
  if(debug){fprintf(stderr, "\n###Entering accord: checking accordance.\n");}
  return reacJisP0(di6,n,m,1,debug);
}

//@@
//Direct approach, but using r-(strong)-compatibility might be better
//The signs of all the J_i (sums of principal minors) of J
int *allJsgns(matrix J, int nlen, int numv, int maxpppdegree, int *pppdegused, int debug){
  int *Jsgn=(int*) malloc(sizeof(int)*(nlen+1));
  int k,allflg;
  ex tmp;
  int debugfull=(debug<=0)?0:debug-1;
  if(debug){fprintf(stderr, "\n###Entering allJsgns: matrix.\n");printmat(J,nlen,nlen);}

  Jsgn[0]=2;//trivial (empty sum)

  for(k=1;k<=nlen;k++){// starts with tr(J) and alternates sign
    tmp=getminorsum0(J,nlen,k);
    Jsgn[k]=ispospoly(expand(tmp), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull);
  }
  if(debug){
    fprintf(stderr, "Signs: ");printvec(Jsgn,nlen+1);
    fprintf(stderr, "Exiting allJsgns.\n");
  }

  return Jsgn;
}

//Version with an integer matrix and pattern matrix as input
int *allJsgns(int **S, int **Slpat, int nlen, int mlen, int debug){
  int *Jsgn=(int*) malloc(sizeof(int)*(nlen+1));
  ex tmp;
  int k;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering allJsgns.\n");}

  Jsgn[0]=2;//trivial (empty sum)

  for(k=1;k<=nlen;k++){// starts with tr(J) and alternates sign
    Jsgn[k]=arecompatr(S, Slpat, k, nlen, mlen, debugfull);
  }
  if(debug){
    fprintf(stderr, "Signs: ");printvec(Jsgn,nlen+1);
    fprintf(stderr, "Exiting allJsgns.\n");
  }

  return Jsgn;
}

//Signs of all principal minor sums for the GK Jacobian 
//(Signs are as returned by ispospoly)
int *reacJsgns(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int numv;
  matrix J=reacJac(AM,n,m,&numv,0);
  return allJsgns(J,n,numv,maxpppdegree,pppdegused,debug);
}

//overloading
int *reacJsgns(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int *ret=reacJsgns(AM, n, m, maxpppdegree, pppdegused, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Version using compatibility: should be faster
int *reacJsgns(int **AM, int n, int m, int debug){
  int **Sl, **S;
  int *Jsgn;
  int debugfull=(debug<=0)?0:debug-1;
  if(debug){fprintf(stderr, "\n###Entering allJsgns: matrix.\n");}
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  //important: if Sl is not a pattern, then arecompatsimp fails
  //So must change Sl to a sign matrix
  sgn_mat(Sl,n,m);
  Jsgn=allJsgns(S,Sl,n,m,debugfull);
  free_imatrix(S,0,n-1,0,m-1);
  free_imatrix(Sl,0,n-1,0,m-1);

  if(debug){
    fprintf(stderr, "Signs: ");printvec(Jsgn,n+1);
    fprintf(stderr, "Exiting allJsgns.\n");
  }

  return Jsgn;
}

//overloading (version using compatibility)
int *reacJsgns(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int *ret=reacJsgns(AM, n, m, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// For use when using Descartes' rule of signs
// maximum number of possible sign changes in a sign vector length n+1, 
// where n is the degree of the polynomial. 
// The output uses the coding of ispospoly. 
// if the sign is -3 or -4, then we assume a sign change, but increase
// (*quest) to reflect the uncertainty.
int signchanges(int *vec, int n, int *quest){
  int i,c=1;//1 means positive, -1 means negative
  int tot=0;
  (*quest)=0;
  for(i=1;i<=n;i++){
    if(vec[i]==-3|| vec[i]==-4){
      c=-c;tot++;(*quest)++;
    }
    else if(c>=0 && (vec[i]==-1|| vec[i]==-2)){
      c=-1;tot++;
    }
    else if(c<=0 && (vec[i]==1|| vec[i]==2)){
      c=1;tot++;
    }
  }
  return tot;
}

//change the sign of every other element (coding of ispospoly)
//vec is assumed to be the coefficient vector of a monic poly of degree n
//So vec has length n+1, but vec[0]=1.
//E.g. if we change from examining spectrum of J to that of -J.
void semiflip(int *vec, int n){
  int i;
  for(i=1;i<=n;i+=2){
    if(vec[i] && vec[i]!=-3 && vec[i]!=-4)
      vec[i]=-vec[i];
  }
  return;
}

//maximum number of eigenvalues (+ve, -ve: using Descartes rule of signs; and also identically zero)
int realspec(matrix J, int n, int numv, int maxpppdegree, int *pppdegused, int *maxpos, int *maxneg, int *totzero, int *quest, int debug){
  int i;
  int *sgnch=allJsgns(J,n,numv,maxpppdegree,pppdegused,debug);

  (*totzero)=0;
  for(i=n;i>=0;i--){
    if(sgnch[i]==0)//identically zero
      (*totzero)++;
    else
      break;
  }
  //Note the length of the vector is n+1 (but argument must be n);
  (*maxneg)=signchanges(sgnch,n,quest);
  semiflip(sgnch,n+1);
  (*maxpos)=signchanges(sgnch,n,quest);
  free((char*)sgnch);
  return  (*maxneg)+(*maxpos)+(*totzero);
}

int reacJrealspec(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int *maxpos, int *maxneg, int *totzero, int *quest, int debug){
  int numv;
  matrix J=reacJac(AM,n,m,&numv,0);
  return realspec(J, n, numv, maxpppdegree, pppdegused, maxpos, maxneg, totzero, quest, debug);
}



//overloading
int reacJrealspec(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int *maxpos, int *maxneg, int *totzero, int *quest, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJrealspec(AM, n, m, maxpppdegree, pppdegused, maxpos, maxneg, totzero, quest, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

int printreacJrealspec(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int maxpos, maxneg, totzero, quest;
  int maxreal=reacJrealspec(di6, n, m, maxpppdegree, pppdegused, &maxpos, &maxneg, &totzero, &quest, debug);
  if(debug){
    fprintf(stderr, "dim=%d, uncertainty:%s, maxreal=%d, totzero=%d, maxpos=%d, maxneg=%d\n",n,quest?"yes":"no", maxreal,totzero,maxpos,maxneg);
    if(maxreal < n)
      fprintf(stderr, "Total real eigenvalues always less than the dimension.\n");
  }
  return maxreal;
}

//MA at equilibria version
int reacJMArealspec(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int *maxpos, int *maxneg, int *totzero, int *quest, int debug){
  int numv;
  int ret;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  ret=realspec(J, n, numv, maxpppdegree, pppdegused, maxpos, maxneg, totzero, quest, debug);
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}

//overloading
int reacJMArealspec(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int *maxpos, int *maxneg, int *totzero, int *quest, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMArealspec(AM, n, m, maxpppdegree, pppdegused, maxpos, maxneg, totzero, quest, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

int printreacJMArealspec(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int maxpos, maxneg, totzero, quest;
  int maxreal=reacJMArealspec(di6, n, m, maxpppdegree, pppdegused, &maxpos, &maxneg, &totzero, &quest, debug);
  if(debug){
    fprintf(stderr, "dim=%d, uncertainty: %s, maxreal=%d, totzero=%d, maxpos=%d, maxneg=%d\n",n,quest?"yes":"no", maxreal,totzero,maxpos,maxneg);
    if(maxreal < n)
      fprintf(stderr, "Total real eigenvalues always less than the dimension.\n");
  }
  return maxreal;
}



//The mass action Jacobian matrix at equilibria admits a pair of purely imaginary eigenvalues
int reacJadmitsIpair(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int numv;
  if(n<2 || n>4){
    fprintf(stderr, "ERROR in reacJadmitsIpair: only for 2D, 3D or 4D systems.\n");
    exit(0);
  }
  matrix J=reacJac(AM,n,m,&numv,0);
  if(n==2)
    return admitsIpair2(J,numv,maxpppdegree,pppdegused,debug);
  else if(n==3)
    return admitsIpair3(J,numv,maxpppdegree,pppdegused,debug);
  return admitsIpair4(J,numv,maxpppdegree,pppdegused,debug);
}

//Overloading
int reacJadmitsIpair(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int numv;
  if(n<2 || n>4){
    fprintf(stderr, "ERROR in reacJadmitsIpair: only for 2D, 3D or 4D systems.\n");
    exit(0);
  }
  matrix J=reacJac(di6,n,m,&numv,0);
  if(n==2)
    return admitsIpair2(J,numv,maxpppdegree,pppdegused,debug);
  else if(n==3)
    return admitsIpair3(J,numv,maxpppdegree,pppdegused,debug);
  return admitsIpair4(J,numv,maxpppdegree,pppdegused,debug);
}



//Better to use concord?
//Is J nonsingular?
int reacJnonsing(int **AM, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(AM,n,m,&numv,0);
  ret=DetSgn(J, n, numv, maxpppdegree, debug);
  if(ret==2 || ret==-2)
    return 1;
  return 0;
}

//Overloading
int reacJnonsing(char *di6, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  ret=DetSgn(J, n, numv, maxpppdegree, debug);
  if(ret==2 || ret==-2)
    return 1;
  return 0;
}

//Concordance: are S and +-Sl r-strongly compatible
int concord(int **AM, int n, int m, int debug){
  int **Sl, **S;
  int ret,strict;
  int debugfull=(debug<=0)?0:debug-1;
  if(debug){fprintf(stderr, "\n###Entering concord: checking concordance.\n");}
  AMtoSSl(AM, n, m, 1, &S, &Sl);
  //important: if Sl is not a pattern, then arecompatsimp fails
  //So must change Sl to a sign matrix
  sgn_mat(Sl,n,m);
  //third argument means only check rank(S) X rank(S) minors
  ret=arecompatsimp(S,Sl,1,n,m,&strict,debugfull);
  free_imatrix(S,0,n-1,0,m-1);
  free_imatrix(Sl,0,n-1,0,m-1);
  if(ret&&strict){//either ret==1 || ret==-1 will do
    if(debug){fprintf(stderr, "The network is concordant.\n");}
    return 1;
  }
  if(debug){fprintf(stderr, "The network is not concordant.\n");}
  return 0;
}

//overloading
int concord(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=concord(AM, n, m, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Concordance: directly using symbolic algebra
//seems to be generally slower
int concord1(int **AM, int n, int m, int maxpppdegree, int debug){
  int **Sl, **S;
  ex minorsum;
  int ret=0,numv,Srank,t, allflg,pppdegused;
  matrix J=reacJac(AM,n,m,&numv,0);
  bool minus=1;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  Srank=matrank(S,n,m);
  minorsum=getminorsum0(J, n, Srank);
  t=ispospoly(expand(minorsum), numv, &allflg, 0, 0, maxpppdegree, &pppdegused, debug);
  if(t==2 || t==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    if(debug){fprintf(stderr, "reqacJconcord1: the system is concordant.\n");} 
    ret=1;
  }
  free_imatrix(S,0,n-1,0,m-1);
  free_imatrix(Sl,0,n-1,0,m-1);
  return ret;
}

int concord1(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=concord1(AM, n, m, maxpppdegree, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Is the Jacobian matrix squared P0?
//direct calculation
int reacJsquaredisP0_direct(int **AM, int n, int m, int maxpppdeg, int debug){
  int numv, ret;
  matrix J=reacJac(AM,n,m,&numv,0);
  matrix J2=multAB(J, J, n, n);//Jacobian squared
  int pppdegused;
  ret=isP0matorth(J2,n,1,maxpppdeg,&pppdegused,debug);//skip top dimension
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}

//overloading
int reacJsquaredisP0_direct(char *di6, int n, int m, int maxpppdeg, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  matrix J2=multAB(J, J, n, n);//Jacobian squared
  int pppdegused;
  ret=isP0matorth(J2,n,1,maxpppdeg,&pppdegused,debug);//skip top dimension
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}

//Is J^2, the GK Jacobian matrix squared, P0? Equivalent to
//sign symmetry of J. This seems to be much more efficient to calculate
//at least for low dimensions
int reacJsquaredisP0(int **AM, int n, int m, int maxpppdeg, int debug){
  int numv,ret;
  matrix J=reacJac(AM,n,m,&numv,0);
  if(debug){fprintf(stderr, "\n###Entering reacJsquaredisP0.\n");} 
  ret=signsym(J,n,maxpppdeg, debug);
  if(debug){fprintf(stderr, "Exiting reacJsquaredisP0: returning %d.\n", ret);} 
  return ret;
}

//overloading
int reacJsquaredisP0(char *di6, int n, int m, int maxpppdeg, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  if(debug){fprintf(stderr, "\n###Entering reacJsquaredisP0.\n");} 
  ret=signsym(J,n,maxpppdeg, debug);  
  if(debug){fprintf(stderr, "Exiting reacJsquaredisP0: returning %d.\n", ret);} 
  return ret;
}


//Is J^[2] P0? "minus" means is -J^[2] P0?
int reacJcomp2isP0(int **AM, int n, int m, bool minus, int maxpppdeg, int debug){
  int numv;
  matrix J=reacJac(AM,n,m,&numv,minus);
  return AdComp2isP0(J, n, maxpppdeg, debug);
}

//overloading
int reacJcomp2isP0(char *di6, int n, int m, bool minus, int maxpppdeg, int debug){
  int numv;
  matrix J=reacJac(di6,n,m,&numv,minus);
  return AdComp2isP0(J, n, maxpppdeg, debug);
}

//Is J^[2] nonsingular?
int reacJcomp2nonsing(int **AM, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  int pppdegused;
  matrix J=reacJac(AM,n,m,&numv,0);
  ret=AdComp2DetPos(J, n, numv, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}


//overloading
int reacJcomp2nonsing(char *di6, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  int pppdegused;
  ret=AdComp2DetPos(J, n, numv, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}

//Is det(J^[2]) strictly signed?
int reacJcomp2detsigned(int **AM, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(AM,n,m,&numv,0);
  int pppdegused;
  ret=AdComp2DetPos(J, n, numv, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2 || ret==0){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}


//overloading
int reacJcomp2detsigned(char *di6, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  int pppdegused;
  ret=AdComp2DetPos(J, n, numv, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2 || ret==0){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}

//Is det(J^[2]) definitely unsigned? (mixed vertices)
int reacJcomp2detunsigned(int **AM, int n, int m, int debug){
  int numv,ret;
  matrix J=reacJac(AM,n,m,&numv,0);
  int pppdegused;
  ret=AdComp2DetPos(J, n, numv, -1, &pppdegused, debug);
  if(ret==-3)
    return 1;
  return 0;
}


//overloading
int reacJcomp2detunsigned(char *di6, int n, int m, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  int pppdegused;
  ret=AdComp2DetPos(J, n, numv, -1, &pppdegused, debug);
  if(ret==-3)
    return 1;
  return 0;
}



//Is det(J^2+I) strictly positive?
int reacJ2pIdetpos(int **AM, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(AM,n,m,&numv,0);
  int pppdegused;
  ret=J2pIwrap(J,n,numv,maxpppdegree,&pppdegused, debug);
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}

//overloading
int reacJ2pIdetpos(char *di6, int n, int m, int maxpppdegree, int debug){
  int numv,ret;
  matrix J=reacJac(di6,n,m,&numv,0);
  int pppdegused;
  ret=J2pIwrap(J,n,numv,maxpppdegree,&pppdegused,debug);
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}


//To test if the determinant at MA equilibria is identically singular
int reacJMAsingular(int **AM, int n, int m){
  int ret, numv;
  int **Q;
  matrix QX;
  ex dt;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  if(Q){
    if(det(Q,n)==0){ret=1;}else{ret=0;}
  }
  else{
    if((dt=det(QX, n))==0){ret=1;}else{ret=0;}
  }
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}



//overloading
int reacJMAsingular(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAsingular(AM, n, m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Is the mass action Jacobian determinant at equilibria always nonsingular?
int reacJMAnonsing(int **AM, int n, int m, int maxpppdegree, int debug){
  int ret, dt, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  if(Q){
    if(det(Q,n)==0){ret=0;}else{ret=1;}
  }
  else{
    dt=DetSgn(QX, n, numv, maxpppdegree, debug);
    if(dt==2 || dt==-2){ret=1;}else{ret=0;}
  }

  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}

//overloading
int reacJMAnonsing(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAnonsing(AM, n, m, maxpppdegree, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Does the mass action Jacobian determinant at equilibria always
// have positive determinant?
int reacJMAdetpos(int **AM, int n, int m, int maxpppdegree, int debug){
  int ret, dt, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  if(Q){
    if(det(Q,n)>0){ret=1;}else{ret=0;}
  }
  else{
    dt=DetSgn(QX, n, numv, maxpppdegree, debug);
    if(dt==2){ret=1;}else{ret=0;}
  }

  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}

//overloading
int reacJMAdetpos(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAdetpos(AM, n, m, maxpppdegree, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//@@
//Direct approach: but using QX might be quicker
//Signs of all the principal minor sums of the MA Jacobian matrix at equilibria
int *reacJMAsgns(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int *ret;
  int numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  ret=allJsgns(J,n,numv,maxpppdegree,pppdegused,debug);
  if(debug){
    fprintf(stderr, "Signs: ");printvec(ret,n+1);
  }

  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}

//overloading
int *reacJMAsgns(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int *ret=reacJMAsgns(AM, n, m, maxpppdegree, pppdegused, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// The mass action Jacobian matrix at equilibria definitely admits a pair 
// of purely imaginary eigenvalues. (Sufficient, but not necessary)
// Only for 2D-4D systems
int reacJMAadmitsIpair(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int ret;
  int numv;
  int **Q;
  matrix QX;
  if(n<2 || n>4){
    fprintf(stderr, "ERROR in reacJMAadmitsIpair: only for 2D, 3D or 4D systems.\n");
    exit(0);
  }
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  if(n==2)
    ret=admitsIpair2(J,numv,maxpppdegree,pppdegused,debug);
  else if(n==3)
    ret=admitsIpair3(J,numv,maxpppdegree,pppdegused,debug);
  else
    ret=admitsIpair4(J,numv,maxpppdegree,pppdegused,debug);

  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}

//overloading
int reacJMAadmitsIpair(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAadmitsIpair(AM, n, m, maxpppdegree, pppdegused, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



// QX is the first factor in J=QX D, where D is a diagonal matrix of unknowns
// which do not figure in QX. J (& QX) are assumed to be generically of rank rk 
// (first factor is a constant matrix of rank rk). We check for potential 
// vanishing of terms in the charpoly of J - since the unknowns in D do not 
// appear in QX, it is sufficient to find either an unsigned principal minor 
// of QX or a pair of principal minors of QX which are oppositely signed. 
// level=0 means we check rk X rk principal minors; level=1 means 
// (rk-1)X(rk-1) principal minors; and so forth. 
// return 0 means nondegenerate; 1 means can be degenerate; 2 means identically degenerate
int degen(matrix QX, int n, int rk, int level, int maxpppdeg, int *pppdegused, int debug){
  int xc[n];
  int flag, deg;
  ex detex, detex_simp;
  int numv, numvinit,allflag;
  int dcfs=0;//assuming integer coeffs (or use polyhasdcfs if not sure)
  int polysgn;
  int sgn=0;
  int debugfull=(debug<=0)?0:debug-1;
  int failflag=0;
  (*pppdegused)=-1;

  if(debug){fprintf(stderr, "\n###Entering degen. Examining %d X %d principal minors of matrix:\n", rk-level, rk-level);printmat(QX, n, n);}

  firstcomb(xc, n, rk-level);flag=1;
    
  while(flag==1){
    detex = symbdetsubmat(QX, n, n, xc, xc, rk-level);
    detex_simp=polyhomsimp(expand(detex), &numvinit, &numv, &deg, 1, dcfs, debug);
    polysgn=ispospoly(detex_simp, numv, &allflag, 0, dcfs, maxpppdeg, pppdegused, debugfull);
    if(polysgn==-3){//unsigned principal minor: return
      if(debug){
	fprintf(stderr, "Exiting degen. An unsigned %d X %d principal minor found.\n", rk-level,rk-level);
	if(debugfull){ 
	  cerr<< detex_simp << endl;
	  printsubmat(QX,xc,xc,rk-level,rk-level);
	}
      }
      return 1;
    }
    else if(polysgn==-4){//couldn't determine sign
      failflag=1;
    }
    else if(((polysgn==2 || polysgn==1) && sgn==-1) || ((polysgn==-2||polysgn==-1) && sgn==1)){//mixed minors
      if(debug){
	fprintf(stderr, "Exiting degen. Found both positive and negative %d X %d principal minors.\n", rk-level, rk-level);
      }
      return 1;
    }
    else if(polysgn==2 || polysgn==1)
      sgn=1;
    else if(polysgn==-2 || polysgn==-1)
      sgn=-1;
    flag=nextcomb(xc, n, rk-level);
  }

  if(sgn==0){
    if(debug){fprintf(stderr, "Identically degenerate at level %d.\n", level);}
    return 2;
  }

  if(failflag){
    if(debug){fprintf(stderr, "Exiting degen. There were minors whose signs could not be determined.\n");}
    return -1;
  }

  if(debug){fprintf(stderr, "Exiting degen. Not degenerate at level %d.\n", level);}
  return 0;  
}

//overloading for integer matrix
int degen(int **QX, int n, int rk, int level, int debug){
  int xc[n];
  int flag;
  int detex;
  int sgn=0;

  if(debug){fprintf(stderr, "\n###Entering degen. Examining %d X %d principal minors of matrix:\n", rk-level, rk-level);printmat(QX, n, n);}

  firstcomb(xc, n, rk-level);flag=1;
    
  while(flag==1){
    detex = detsubmat(QX, n, n, xc, xc, rk-level);
    if((detex>0 && sgn==-1) || (detex<0 && sgn==1)){//mixed minors
      if(debug){
	fprintf(stderr, "Exiting degen. Found both positive and negative %d X %d principal minors.\n", rk-level, rk-level);
      }
      return 1;
    }
    else if(detex>0)
      sgn=1;
    else if(detex<0)
      sgn=-1;
    flag=nextcomb(xc, n, rk-level);
  }

  if(sgn==0){
    if(debug){fprintf(stderr, "Exiting degen. Identically degenerate at level %d.\n", level);}
    return 2;
  }

  if(debug){fprintf(stderr, "Exiting degen. Not degenerate at level %d.\n", level);}
  return 0;  
}

//Just want to check if identically degenerate
int idegen(matrix QX, int n, int rk, int level, int debug){
  int xc[n];
  int flag;
  ex detex;

  if(debug){fprintf(stderr, "\n###Entering idegen. Examining %d X %d principal minors of matrix:\n", rk-level, rk-level);printmat(QX, n, n);}

  firstcomb(xc, n, rk-level);flag=1;
    
  while(flag==1){
    detex = symbdetsubmat(QX, n, n, xc, xc, rk-level);
    if(expand(detex)!=0){
      if(debug){
	fprintf(stderr, "Exiting idegen. Found a nonzero %d X %d principal minor.\n", rk-level, rk-level);
      }
      return 0;
    }
    flag=nextcomb(xc, n, rk-level);
  }

  if(debug){fprintf(stderr, "Exiting idegen. All %d X %d principal minors are zero.\n", rk-level, rk-level);}
  return 1; 
}


//Rapid version - just want to check if identically degenerate
//overloading for integer matrix
int idegen(int **QX, int n, int rk, int level, int debug){
  int xc[n];
  int flag;

  if(debug){fprintf(stderr, "\n###Entering idegen. Examining %d X %d principal minors of matrix:\n", rk-level, rk-level);printmat(QX, n, n);}

  firstcomb(xc, n, rk-level);flag=1;
  while(flag==1){
    if(detsubmat(QX, n, n, xc, xc, rk-level)!=0){
      if(debug){
	fprintf(stderr, "Exiting idegen. Found a nonzero %d X %d principal minor.\n", rk-level, rk-level);
      }
      return 0;
    }
    flag=nextcomb(xc, n, rk-level);
  }

  if(debug){fprintf(stderr, "Exiting idegen. All %d X %d principal minors are zero.\n", rk-level, rk-level);}
  return 1;  
}



// Checking for second order degeneracies 
// We assume it has been checked that the rank of QX (& J) can drop below rk. 
// Here we check if the rank of J can drop further: can next term in 
// its characteristic polynomial can also be zero? We find out by checking if
// QX has any (rk-1)X(rk-1) principal minors which are unsigned; or of differ
// in sign
int degen2(matrix QX, int n, int rk, int maxpppdeg, int *pppdegused, int debug){
  return degen(QX, n, rk, 1, maxpppdeg, pppdegused, debug);
}

//overloading: integer matrix
int degen2(int **QX, int n, int rk, int debug){
  return degen(QX, n, rk, 1, debug);
}


//Check if mass action equilibria can be second-order degenerate (already
//assuming they can be degenerate)
int reacJMAdegen2(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int ret, rk;
  int numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv, &rk);
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering reacJMAdegen2.\n");}

  //Alternative approach using J directly 
  /* int allflg; */
  /* ex Jp=getminorsum0(J,n,rk-1); */
  /* int retp=ispospoly(expand(Jp), numv, &allflg, 0, 0, maxpppdegree, pppdegused, debugfull); */
  /* if(retp==2 || retp==-2) */
  /*   return 0; */
  /* return 1; */

  if(Q){
    ret=degen2(Q, n, rk, debugfull);
    free_imatrix(Q, 0, n-1, 0, n-1);
    if(debug){fprintf(stderr, "Exiting reacJMAdegen2. %s degen2.\n", ret?"Satisfies":"Fails");}
    return ret;
  }

  ret=degen2(QX,n,rk,maxpppdegree,pppdegused,debugfull);
  if(debug){fprintf(stderr, "Exiting reacJMAdegen2. %s degen2.\n", ret?"Satisfies":"Fails");}
  return ret;
}

//overloading
int reacJMAdegen2(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAdegen2(AM,n,m,maxpppdegree,pppdegused,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



// Check if the network has the potential for Bogdanov-Takens bifurcations
// We are testing whether the highest order terms in the characteristic polynomial 
// of the Jacobian matrix at MA equilibria can simultaneously be zero. (We discount
// the terms which are automatically zero.) The algorithm returns 0 if B-T bifurcations
// are ruled out. It returns 1 if they could not be ruled out (this does not
// mean that they must occur.)
int reacJMABT(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int ret, rk;
  int numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv, &rk);
  int debugfull=(debug<=0)?0:debug-1;

  if(debug)
    fprintf(stderr, "\n###Entering reacJMABT.\n");

  ret=BTpotential(J, n, rk, numv, maxpppdegree, pppdegused, debugfull);

  if(debug)
    fprintf(stderr, "Exiting reacJMABT. %s.\n", ret?"B-T bifurcations could not be ruled out":"B-T bifurcations were ruled out");

  return ret;
}

//overloading
int reacJMABT(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMABT(AM,n,m,maxpppdegree,pppdegused,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}




//Are equilibria of the mass action network always nondegenerate?
//(Only use on DN networks: reacJMAeq has no meaning and fails otherwise)
int reacMAeqconcord(int **AM, int n, int m, int maxpppdegree, int debug){
  int **Q;
  int ret=0,Srank;
  int numv;
  matrix QX;
  int pppdegused=-1;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv, &Srank);

  if(debug){fprintf(stderr, "\n###Entering reacMAeqconcord (are mass action equilibria always nondegenerate?)\n");}

  //alternative which directly uses J
  /* ex Jf=getminorsum0(J,n,Srank); */
  /* int allflg; */
  /* int retf=ispospoly(expand(Jf), numv, &allflg, 0, 0, maxpppdegree, &pppdegused, debug); */
  /* if(retf==2 || retf==-2) */
  /*   return 1; */
  /* return 0; */


  if(Q){
    if(!degen(Q,n,Srank,0, debug)){ret=1;}else{ret=0;}
    free_imatrix(Q, 0, n-1, 0, n-1);
    if(debug){fprintf(stderr, "The MA network %s at equilibria.\n", ret?"is nondegenerate":"can be degenerate");}
  }
  else{
    if(!degen(QX,n,Srank,0,maxpppdegree,&pppdegused,debug)){ret=1;}else{ret=0;}
    if(debug){fprintf(stderr, "The MA network %s nondegenerate at equilibria.\n", ret?"is":"could not be confirmed to be");}
  }

  return ret;
}

//overloading
int reacMAeqconcord(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacMAeqconcord(AM,n,m,maxpppdegree,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



//Is JMA positive definite? or negative definite
//untested algorithm (probably rare for it to hold)
//(Only use on DN networks)
int reacJMAposdef(int **AM, int n, int m, bool minus, int maxpppdeg, bool strict, int debug){
  int ret, numv,i,j;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  matrix Js(n,n);
  int pppdegused;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      Js(i,j) = J(i,j)+J(j,i);
      if(minus)
	Js(i,j)=-Js(i,j);
    }
  }
  if(strict)
    ret=isposdef(Js,n,maxpppdeg,&pppdegused,1);
  else
    ret=isP0matorth(Js,n,0,maxpppdeg,&pppdegused,1);
  if(Q){free_imatrix(Q, 0, n-1, 0, n-1);}
  if(ret && debug && pppdegused>=0)
    fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  return ret;
}

//overloading
int reacJMAposdef(char *di6, int n, int m, bool minus, int maxpppdeg, bool strict, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAposdef(AM, n, m, minus, maxpppdeg, strict, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Is the first factor in JMAeq positive definite?
//In the case that the stoichiometric matrix has a rank-1 kernel
//Q is an integer matrix; otherwise it is NULL, and QX is a symbolic matrix
int reacQMAposdef(int **Q, matrix QX, int n, bool minus, bool strict, int maxpppdeg, int debug){
  int ret,i,j;
  matrix QXs(n,n);
  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering reacQMAposdef. Checking if the mass action Jacobian matrix is similar to a %s %sdefinite matrix.\n", minus?"negative":"positive", strict?"":"semi");}

  if(Q){//special case. Symmetrise
    int **Qs=transposemat(Q,n,n);
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	Qs[i][j]+=Q[i][j];
	if(minus)
	  Qs[i][j]=-Qs[i][j];
      }
    }
    if(debugfull){fprintf(stderr, "Symmetrised matrix Qs:\n");printmat(Qs,n,n);}
    if(strict)
      ret=isposdef(Qs,n);
    else{
      ret=isP0matorth(Qs,n,0,debugfull);
      /* if(!ret){//try reduction */
      /* 	reduce_mat(Q, n, n, 1); */
      /* 	int **Qs1=transposemat(Q,n,n); */
      /* 	for(i=0;i<n;i++){ */
      /* 	  for(j=0;j<n;j++){ */
      /* 	    Qs1[i][j]+=Q[i][j]; */
      /* 	    if(minus) */
      /* 	      Qs1[i][j]=-Qs1[i][j]; */
      /* 	  } */
      /* 	} */
      /* 	ret=isP0matorth(Qs1,n,0,debugfull); */
      /* 	if(ret){ */
      /* 	  printmat(Qs1,n,n); */
      /* 	  fprintf(stderr, "got here\n");exit(0); */
      /* 	} */
      /* } */
    }
    free_imatrix(Qs, 0, n-1, 0, n-1);
  }
  else{//Note: QX only defined if Q is NULL; symmetrise
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	QXs(i,j) = QX(i,j)+QX(j,i);
	if(minus)
	  QXs(i,j)=-QXs(i,j);
      }
    }
    if(strict)
      ret=isposdef(QXs,n,maxpppdeg,&pppdegused,debugfull);
    else
      ret=isP0matorth(QXs,n,0,maxpppdeg,&pppdegused,debugfull);
    if(ret && debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  }

  if(debug && ret)
    fprintf(stderr, "The mass action Jacobian matrix is similar to a %s %sdefinite matrix.\n", minus?"negative":"positive", strict?"":"semi");
  else if(debug && !ret)
    fprintf(stderr, "Could not determine that the mass action Jacobian matrix is similar to a %s %sdefinite matrix.\n", minus?"negative":"positive", strict?"":"semi");

  return ret;
}

//overloading (Only use on DN networks)
int reacQMAposdef(int **AM, int n, int m, bool minus, bool strict, int maxpppdeg, int debug){
  int ret, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  ret=reacQMAposdef(Q, QX, n, minus, strict, maxpppdeg, debug);
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}



//overloading
int reacQMAposdef(char *di6, int n, int m, bool minus, bool strict, int maxpppdeg, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacQMAposdef(AM, n, m, minus, strict, maxpppdeg, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//deprecated
//Is Q positive definite? or negative definite
//This doesn't imply the same for J=QD; but Q negative definite
//does imply that J is Hurwitz.
int reacQMA1posdef(int **AM, int n, int m, bool minus, bool strict){
  int ret;
  int **Q = reacQMA1(AM, n, m, 0);
  int **Qt = transposemat(Q, n, n);
  matadd(Qt,Q,n,n);
  if(minus)
    minusmat(Qt,n,n);
  if(strict)
    ret=isposdef(Qt,n);
  else
    ret=isP0matorth(Qt,n,0,1);
  free_imatrix(Q, 0, n-1, 0, n-1);
  free_imatrix(Qt, 0, n-1, 0, n-1);
  return ret;
}

//overloading
int reacQMA1posdef(char *di6, int n, int m, bool minus, bool strict){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacQMA1posdef(AM, n, m, minus, strict);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



//Is the special MA Jacobian matrix P0 or P0-? 
//Only need to check QX
//(Only use on DN networks)
int reacJMAisP0(int **AM, int n, int m, bool minus, int maxpppdeg, int debug){
  int ret, numv;
  int **Q;
  matrix QX,QX1;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;
  //if(!Q){cerr << QX << endl;exit(0);}
  if(Q){
    if(minus)
      minusmat(Q,n,n);
    ret=isP0matorth(Q,n,0,debugfull);
    free_imatrix(Q, 0, n-1, 0, n-1);
  }
  else{
    if(minus){
      QX1=minusmat(QX,n,n);
      ret=isP0matorth(QX1,n,0,maxpppdeg,&pppdegused,debugfull);
    }
    else
      ret=isP0matorth(QX,n,0,maxpppdeg,&pppdegused,debugfull);
    if(ret && debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  }
  if(debug){
    if(ret)
      fprintf(stderr, "The mass action Jacobian matrix at equilibria is P0%s\n",minus?"-":"");
    else
      fprintf(stderr, "The mass action Jacobian matrix at equilibria could not be confirmed to be P0%s\n",minus?"-":"");
  }
  return ret;
}

//overloading
//(Only use on DN networks)
int reacJMAisP0(char *di6, int n, int m, bool minus, int maxpppdeg, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAisP0(AM, n, m, minus, maxpppdeg, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Does the special MA Jacobian matrix have mixed trace?
//(Only use on DN networks)
int reacJMAmixedtrace(int **AM, int n, int m, bool minus, int maxpppdeg, int debug){
  int ret, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;
  //if(!Q){cerr << QX << endl;exit(0);}
  if(Q){
    ret=mixedtrace(Q,n,debugfull);
    free_imatrix(Q, 0, n-1, 0, n-1);
  }
  else{
    ret=mixedtrace(QX,n,maxpppdeg,&pppdegused,debugfull);
    if(ret && debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  }

  if(debug){
    if(ret==1)
      fprintf(stderr, "The mass action Jacobian matrix at equilibria has mixed trace\n");
    else if(ret==-1)
      fprintf(stderr, "It could not be confirmed whether the mass action Jacobian matrix at equilibria has mixed trace.\n");
    else
      fprintf(stderr, "The mass action Jacobian matrix at equilibria does not have mixed trace\n");
  }
  return ret;
}

//overloading
//(Only use on DN networks)
int reacJMAmixedtrace(char *di6, int n, int m, bool minus, int maxpppdeg, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAmixedtrace(AM, n, m, minus, maxpppdeg, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}





//J at MA equilibria is P0-
int MAeqaccord(int **AM, int n, int m, int maxpppdeg, int debug){
  if(debug){fprintf(stderr, "\n###Entering MAeqaccord (checking if the MA Jacobian matrix at equilibria is P0-)\n");}
  return reacJMAisP0(AM,n,m,1,maxpppdeg,debug);
}
//overloading
int MAeqaccord(char *di6, int n, int m, int maxpppdeg, int debug){
  if(debug){fprintf(stderr, "\n###Entering MAeqaccord (checking if the MA Jacobian matrix at equilibria is P0-)\n");}
  return reacJMAisP0(di6,n,m,1,maxpppdeg,debug);
}

  
//deprecated
//Is the special MA Jacobian matrix squared P0? (1D pos kernel)
int reacJMA1squaredisP0(int **AM, int n, int m){
  int **Q = reacQMA1(AM, n, m, 0);
  int ret;
  ret=signsym(Q,n);
  free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}
//overloading
int reacJMA1squaredisP0(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1squaredisP0(AM, n, m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Is JMAeq^2 P0. Necessary and sufficient to examine Q/QX for sign symmetry
int reacJMAsquaredisP0(int **AM, int n, int m, int maxpppdeg, int debug){
  int ret, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  if(debug){fprintf(stderr, "\n###Entering reacJMAsquaredisP0: checking if the J^2 is always P0 at MA equilibria.\n");}
  if(Q){
    ret=signsym(Q,n);
    free_imatrix(Q, 0, n-1, 0, n-1);
    if(debug){fprintf(stderr, "J^2 is %salways P0 at MA equilibria.\n",ret?"":"not ");}
  }
  else{
     ret=signsym(QX,n, maxpppdeg, debug);
     if(debug){fprintf(stderr, "%sJ^2 is always P0 at MA equilibria.\n",ret?"":"Could not confirm that ");}
  }

  return ret;
}

//overloading
int reacJMAsquaredisP0(char *di6, int n, int m, int maxpppdeg, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAsquaredisP0(AM, n, m, maxpppdeg, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//deprecated
//Is J^[2] for the special MA Jacobian matrix (1D pos kernel) P0 or P0-? 
int reacJMA1comp2isP0(int **AM, int n, int m, bool minus, int maxpppdeg, int debug){
  matrix J=reacJacMA1(AM, n, m, minus, 1);
  return AdComp2isP0(J, n, maxpppdeg, debug);
}
//overloading
int reacJMA1comp2isP0(char *di6, int n, int m, bool minus, int maxpppdeg, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1comp2isP0(AM, n, m, minus, maxpppdeg, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Not deprecated
//Is J^[2] for the special MA Jacobian matrix P0-?
//Done via factorisation of J^[2]; mainly for testing purposes
//for systems with 1D pos kernel
//misses some, but much faster than reacJMAcomp2isP0
int reacJMA1comp2simppair(int **AM, int n, int m, bool minus, int debug){
  int ret, numv, strict;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  int **S2mat, **V2imat;
  int **Di=idmat(n, minus);//minus
  long cnk=n*(n-1)/2;
  long cmk=n*n;

  if(!Q){
    fprintf(stderr, "ERROR in reacJMA1comp2simppair: only for systems with 1D positive kernel. EXITING.\n");
    exit(0);
  }

  genS2from2(Q, Di, n, n, &S2mat, &V2imat);
  if(debug){
    printmat(S2mat,cnk,cmk);
    printmat(V2imat,cnk,cmk);
  }

  ret=arecompatsimp(S2mat, V2imat, cnk, cmk, &strict, debug);

  free_imatrix(S2mat, 0, cnk-1, 0, cmk-1);
  free_imatrix(V2imat, 0, cnk-1, 0, cmk-1);
  free_imatrix(Di,0,n-1,0,n-1);
  free_imatrix(Q, 0, n-1, 0, n-1);
  if(ret==1)//can also be minus 1
    return 1;
  return 0;
}

//overloading
int reacJMA1comp2simppair(char *di6, int n, int m, bool minus, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1comp2simppair(AM, n, m, minus, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Is J^[2] for the special MA Jacobian matrix P0 or P0-? 
int reacJMAcomp2isP0(int **AM, int n, int m, bool minus, int maxpppdeg, int debug){
  int ret, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);

  if(minus){
    matrix Jm=minusmat(J,n,n);
    ret=AdComp2isP0(Jm, n, maxpppdeg, debug);
  }
  else
    ret=AdComp2isP0(J, n, maxpppdeg, debug);

  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}


//overloading
int reacJMAcomp2isP0(char *di6, int n, int m, bool minus, int maxpppdeg, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAcomp2isP0(AM, n, m, minus, maxpppdeg, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Does the DSR2 graph for the MA1 system satisfy Condition *?
//1D kernel only
int DSR2MA1CondStar(int **AM, int n, int m, int debug){
  int ret, report, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  int **Di=idmat(n, 1);//minus

  if(!Q){
    fprintf(stderr, "ERROR in DSR2MA1CondStar: only for systems with 1D positive kernel. EXITING.\n");
    exit(0);
  }
   
  if(debug){fprintf(stderr, "Q= ");printmat(Q,n,n);}
  ret=DSR2CondStar(Q, Di, n, n, 0, &report, 0, debug);

  free_imatrix(Di,0,n-1,0,n-1);
  free_imatrix(Q, 0, n-1, 0, n-1);
  return ret;
}

//overloading: di6 input
int DSR2MA1CondStar(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=DSR2MA1CondStar(AM, n, m, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;

}
 

//deprecated
//Is J^[2] for the special MA Jacobian matrix (1D pos kernel) nonsingular?
int reacJMA1comp2nonsing(int **AM, int n, int m, int maxpppdegree){
  int ret;
  matrix J=reacJacMA1(AM, n, m, 0, 1);
  int pppdegused;
  int debug=1;
  ret=AdComp2DetPos(J, n, n, maxpppdegree, &pppdegused, 0);
  if(ret==2 || ret==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}


//overloading
int reacJMA1comp2nonsing(char *di6, int n, int m, int maxpppdegree){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1comp2nonsing(AM, n, m, maxpppdegree);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Sign of the determinant of J^[2] where J is the mass action Jacobian
//matrix at equilibria
int reacJMAcomp2det(int **AM, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int ret, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  int debugfull=(debug<=0)?0:debug-1;
  if(debug){fprintf(stderr, "\n###Entering reacJMAcomp2det: examining the determinant of J^[2] (where J is the MA Jacobian)\n");}

  ret=AdComp2DetPos(J, n, numv, maxpppdegree, pppdegused, debugfull);
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);

  if(debug){fprintf(stderr, "reacJMAcomp2det returned %d (return codes: -4: unable to determine; -3: definitely unsigned; -2: strictly negative; -1: nonpositive; 0: the zero polynomial; 1: nonnegative; 2: strictly positive)\n", ret);}
  return ret;
}

//overloading
int reacJMAcomp2det(char *di6, int n, int m, int maxpppdegree, int *pppdegused, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAcomp2det(AM, n, m, maxpppdegree, pppdegused, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Is J^[2] for the MA Jacobian matrix at equilibria nonsingular
int reacJMAcomp2nonsing(int **AM, int n, int m, int maxpppdegree, int debug){
  int pppdegused;
  int ret=reacJMAcomp2det(AM, n, m, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}

//overloading
int reacJMAcomp2nonsing(char *di6, int n, int m, int maxpppdegree, int debug){
  int pppdegused;
  int ret=reacJMAcomp2det(di6, n, m, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused); 
    return 1;
  }
  return 0;
}

// does J^[2] for the MA Jacobian matrix
// have strictly signed determinant? (Good enough to preclude 
// Hopf bifurcations in dimension <=3)
int reacJMAcomp2detsigned(int **AM, int n, int m, int maxpppdegree, int debug){
  int pppdegused;
  int ret=reacJMAcomp2det(AM, n, m, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2 || ret==0){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}

//overloading
int reacJMAcomp2detsigned(char *di6, int n, int m, int maxpppdegree, int debug){
  int pppdegused;
  int ret=reacJMAcomp2det(di6, n, m, maxpppdegree, &pppdegused, debug);
  if(ret==2 || ret==-2 || ret==0){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    return 1;
  }
  return 0;
}

// does J^[2] for the MA Jacobian matrix
// definitely have unsigned determinant? (mixed vertices)
int reacJMAcomp2detunsigned(int **AM, int n, int m, int debug){
  int pppdegused;
  int ret=reacJMAcomp2det(AM, n, m, -1, &pppdegused, debug);
  if(ret==-3)
    return 1;
  return 0;
}

//overloading
int reacJMAcomp2detunsigned(char *di6, int n, int m, int debug){
  int pppdegused;
  int ret=reacJMAcomp2det(di6, n, m, -1, &pppdegused, debug);
  if(ret==-3)
    return 1;
  return 0;
}


//Are partial derivatives of det(J^[2]) in terms of x never simultaneously zero?
//If so, then we can automatically satisfy the transversality condition
//at any potential Hopf bifurcation points
int reacJMAcomp2detnonstationary(int **AM, int n, int m, int maxpppdegree, int debug){
  int numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);
  int debugfull=(debug<=0)?0:debug-1;
  ex detex=AdComp2Det(J, n, debugfull);
  char **pvars;
  ex v[numv];
  ex d, d1=0;
  int i,t;
  int allflg;
  int pppdegused;
  if(debug){fprintf(stderr, "\n###Entering reacJMAcomp2detnonstationary (are partial derivatives of det(J^[2]) in terms of x ever simultaneously zero?)\n");}

  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);
  polyvars(detex, &pvars);
  polyvarss(detex, v, pvars, numv);//get the variables themselves
  freearraydat(pvars, numv);
  for(i=0;i<n;i++){
    d=detex.diff(ex_to<symbol>(v[i]));
    t=heuristic_squares(d,numv,debugfull);
    if(debugfull){cerr << "Examining " << d << "; outcome: " << t << endl;}
    if(t==2 || t==-2){
      if(debug){fprintf(stderr, "reacJMAcomp2detnonstationary: some partial derivative of det(J^[2]) is always %s.\n",(t==2)?"positive":"negative");} 
      return 1;//nondegenerate
    }
    d1+=pow(d,2);
  }

  t=ispospoly(expand(d1), numv, &allflg, 0, 0, maxpppdegree, &pppdegused, debugfull);
  //cerr << t << endl;
  if(t==2 || t==-2){
    if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
    if(debug){fprintf(stderr, "reacJMAcomp2detnonstationary: partial derivatives of det(J^[2]) are never simultaneously zero.\n");} 
    return 1;
  }

  if(debug){fprintf(stderr, "reacJMAcomp2detnonstationary: can't confirm whether partial derivatives of det(J^[2]) are ever simultaneously zero.\n");} 
  return 0;
}

//overloading
int reacJMAcomp2detnonstationary(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAcomp2detnonstationary(AM, n, m, maxpppdegree, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



//deprecated
int reacJMA1comp2detsigned(int **AM, int n, int m, int maxpppdegree){
  int ret;
  matrix J=reacJacMA1(AM, n, m, 0, 1);
  int pppdegused;
  ret=AdComp2DetPos(J, n, n, maxpppdegree, &pppdegused, 0);
  if(ret==2 || ret==-2 || ret==0)
    return 1;
  return 0;
}

//deprecated
int reacJMA1comp2detsigned(char *di6, int n, int m, int maxpppdegree){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1comp2detsigned(AM, n, m, maxpppdegree);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//deprecated
int reacJMA1comp2detunsigned(int **AM, int n, int m){
  int ret;
  matrix J=reacJacMA1(AM, n, m, 0, 1);
  int pppdegused;
  ret=AdComp2DetPos(J, n, n, -1, &pppdegused, 1);
  if(ret==-3)
    return 1;
  return 0;
}

//deprecated
int reacJMA1comp2detunsigned(char *di6, int n, int m){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1comp2detunsigned(AM, n, m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

// Are equilibria of the MA network *all* degenerate? 
// (Assume the network is not dynamically trivial.) 
int reacJMAdegenerate(int **S, int **Sl, int n, int m, int debug){
  int ret;
  int **Q;
  matrix QX;
  matrix J;
  int numv,rk;

  if(debug){fprintf(stderr, "\n###Entering reacJMAdegenerate: checking if all MA equilibria are degenerate.\n");}

  if(hasposlimvec(S,n,m)){
    fprintf(stderr, "ERROR in reacJMAdegenerate. The network is dynamically trivial. EXITING.\n");
    exit(0); 
  }

  J=reacJMAeq(S, Sl, n, m, &Q, &QX, &numv, &rk);
  if(Q){
    if(idegen(Q,n,rk,0,debug)){
      if(debug){fprintf(stderr, "MA equilibria are always degenerate.\n");}
      ret=1;
    }
    else{
      if(debug){fprintf(stderr, "MA equilibria not always degenerate.\n");}
      ret=0;
    }
    free_imatrix(Q, 0, n-1, 0, n-1);
    return ret;
  }

  if(idegen(QX,n,rk,0,debug)){
    if(debug){fprintf(stderr, "MA equilibria are always degenerate.\n");}
    ret=1;
  }
  else{
    if(debug){fprintf(stderr, "MA equilibria not always degenerate.\n");}
    ret=0;
  }
  return ret;
}

//overloading
int reacJMAdegenerate(int **AM, int n, int m, int debug){
  int **Sl, **S;
  int ret;
  bool minus=0;
  AMtoSSl(AM, n, m, minus, &S, &Sl);
  ret=reacJMAdegenerate(S, Sl, n, m, debug);
  free_imatrix(S,0,n-1,0,m-1);
  free_imatrix(Sl,0,n,0,m-1);
  return ret;
}


//overloading, l layers
int reacJMAdegenerate(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  int ret=reacJMAdegenerate(AM, n, m, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//2 layers - deprecated
int reacJMAdegenerateold(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMAdegenerate(AM, n, m, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Is det(J^2+I) for the special MA Jacobian matrix (1D pos kernel) positive?
int reacJ2pIMAdetpos(int **AM, int n, int m, int maxpppdegree, int debug){
  int ret, numv;
  int **Q;
  matrix QX;
  matrix J=reacJMAeq(AM, n, m, &Q, &QX, &numv);

  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;
  if(debug){fprintf(stderr, "\n###Entering reacJ2pIMAdetpos: examining the sign of det(J^2+I) (where J is the MA Jacobian)\n");}
  ret=J2pIwrap(J,n,numv,maxpppdegree,&pppdegused,debugfull);
  if(debug && pppdegused>=0)
      fprintf(stderr, "Used SDP to degree %d\n", pppdegused);
  if(Q)
    free_imatrix(Q, 0, n-1, 0, n-1);

  if(debug){fprintf(stderr, "%sdet(J^2+I) is strictly positive\n", ret?"":"Unable to determine whether ");}
  return ret;
}

//overloading
int reacJ2pIMAdetpos(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJ2pIMAdetpos(AM, n, m, maxpppdegree, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//deprecated
//Is det(J^2+I) for the special MA Jacobian matrix (1D pos kernel) positive?
int reacJ2pIMA1detpos(int **AM, int n, int m, int maxpppdegree, int debug){
  matrix J=reacJacMA1(AM, n, m, 0, 1);
  int pppdegused;
  return J2pIwrap(J,n,n,maxpppdegree,&pppdegused,debug);
}


//deprecated
//overloading
int reacJ2pIMA1detpos(char *di6, int n, int m, int maxpppdegree, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJ2pIMA1detpos(AM, n, m, maxpppdegree, debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//deprecated
//Are partial derivatives of det(J^[2]) in terms of x never simultaneously zero?
//If so, then we can automatically satisfy the transversality condition
//at any potential Hopf bifurcation points
int reacJMA1comp2detnonstationary(int **AM, int n, int m, int maxpppdegree, bool q){
  matrix J=reacJacMA1(AM, n, m, 0, q);
  ex detex=AdComp2Det(J, n, 1-q);
  char **pvars;
  int numv=polyvars(detex, &pvars);
  ex v[numv];
  ex d, d1=0;
  int i,t;
  int allflg;
  int pppdegused;
  polyvarss(detex, v, pvars, numv);//get the variables themselves
  freearraydat(pvars, numv);
  for(i=0;i<n;i++){
    d=detex.diff(ex_to<symbol>(v[i]));
    t=heuristic_squares(d,numv,1-q);
    if(!q){cerr << "Examining " << d << "; outcome: " << t << endl;}
    if(t==2 || t==-2){
      if(!q){fprintf(stderr, "MA1Jcomp2nondegen: some derivative of det(J^[2]) is never zero.\n");} 
      return 1;//nondegenerate
    }
    d1+=pow(d,2);
  }
  t=ispospoly(expand(d1), numv, &allflg, 0, 0, maxpppdegree, &pppdegused, q);
  cerr << t << endl;
  if(t==2 || t==-2){
    if(!q){fprintf(stderr, "MA1Jcomp2nondegen: derivatives of det(J^[2]) are never simultaneously zero.\n");} 
    return 1;
  }

  if(!q){fprintf(stderr, "MA1Jcomp2nondegen: can't confirm that derivatives of det(J^[2]) are never simultaneously zero.\n");} 
  return 0;
}

//overloading (deprecated)
int reacJMA1comp2detnonstationary(char *di6, int n, int m, int maxpppdegree, bool q){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
//Petri Net adjacency matrix
  int ret=reacJMA1comp2detnonstationary(AM, n, m, maxpppdegree, q);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



//read a single reaction from an open file (either a single line, or a block)
int reacread(FILE *fd, bool isblock, char line[], int maxl, char block[], int maxblock){
  int linelen;
  int blocklen=0;
  int pass=0;

  if(!isblock){
    linelen=getline0(fd, line, maxl);
    while(linelen>0 && iscomline(line))
      linelen=getline0(fd, line, maxl);
    return linelen;
  }

  linelen=gtline(fd,line,1000);
  while(linelen > 0 && !strstr(line, "******")){
    if(!iscomline(line)){
      blocklen+=linelen;
      if(blocklen>maxblock-10){
	fprintf(stderr, "ERROR in reacread. Block length got too long. EXITING.\n");exit(0);
      }
      if(pass==0){strcpy(block, line);pass=1;}else{strcat(block, line);}
    }
    linelen=gtline(fd,line,1000);
  }
  return blocklen;
}


//Create report on a single reaction given via the PN AM, the 
//stoichiometric matrix, and the left stoichiometric matrix
//Note: Sl is *minus* the left stoichiometric matrix
int reacreport(int **AM, int **S, int **Sl, int n, int m, int Srank, int maxpppdeg, char ***outlist, int **outvec, char *filter, int debug){
  int totv=0;
  int gen, nonautocatalytic, nontriv, bdc;
  int connect, Spersist, endotact, WR, norm;
  int DSRreport, DSRCstar,acc,semiacc,MAeqacc;
  int conc,semiconc,MAeqconc;
  int nohopf,pppdegused,QMAnegsemidef;
  int ipairGK,ipairMA;
  int MAdegen;
  char *str;
  int entries;

  if(debug){
    fprintf(stderr, "\n########Entering reacreport.\n");
    entries=nonzentries(AM,n+m,n+m);
    str=CRNamtostr(AM, n, m, entries, 0);
    fprintf(stderr, "%s\n", str);
    free(str);
  }

  //Always calculate: used often
  nontriv= (int)(DN(AM,n,m));

  //basic tests (should all be quick and work for large systems
  if(!filter || strstr(filter, "basic")){
    //species
    addnewto1Darray(outvec, totv, n);
    totv=addv2(totv, (char*)"numspec", outlist);
    //reactions
    addnewto1Darray(outvec, totv, m);
    totv=addv2(totv, (char*)"numreac", outlist);
    //rank
    addnewto1Darray(outvec, totv, Srank);
    totv=addv2(totv, (char*)"rank", outlist);
    //genuine
    gen= (int)(genuine(AM,n,m));
    addnewto1Darray(outvec, totv, gen);
    totv=addv2(totv, (char*)"genuine", outlist);
    //nonautocatalytic
    nonautocatalytic= (int)(nonautocat(AM,n,m));
    addnewto1Darray(outvec, totv, nonautocatalytic);
    totv=addv2(totv, (char*)"non-autocatalytic", outlist);
    //DN
    addnewto1Darray(outvec, totv, nontriv);
    totv=addv2(totv, (char*)"DN", outlist);
    //bdclass
    if(Srank>=n){bdc=0;}else{bdc=(int)(bdclass(AM,n,m,debug));}
    addnewto1Darray(outvec, totv, bdc);
    totv=addv2(totv, (char*)"bdclass", outlist);
    //PNconnect (not decomposable)
    connect=PNconnect(AM,n,m);
    addnewto1Darray(outvec, totv, connect);
    totv=addv2(totv, (char*)"PNconnect", outlist);
    //structurally persistent (no critical siphons)
    Spersist=(int)(structpersist(AM,n,m,debug));
    addnewto1Darray(outvec, totv, Spersist);
    totv=addv2(totv, (char*)"structpersist", outlist);

    //normal?
    norm=(int)(normal(AM,n,m));
    addnewto1Darray(outvec, totv, norm);
    totv=addv2(totv, (char*)"normal", outlist);
    //weakly reversible (!DN ==> !WR (WR ==> 1 in kernel of Si), !normal --> !WR (Lemma E.1 in bp_SIADS))
    if(!nontriv || !norm){WR=0;}else{WR=(int)(CRNweakrev(AM,n,m));}
    addnewto1Darray(outvec, totv, WR);
    totv=addv2(totv, (char*)"WR", outlist);
    //DSR condition * (check both)
    DSRCstar=DSRCondStarPN(AM, n, m, &DSRreport, 2, debug);
    addnewto1Darray(outvec, totv, DSRCstar);
    totv=addv2(totv, (char*)"DSR Condition *", outlist);
    //MAdegen
    if(!nontriv){
      MAdegen=-1;
      addnewto1Darray(outvec, totv, MAdegen);
      totv=addv2(totv, (char*)"MAdegen not relevant (dynamically trivial network)", outlist);
    }
    else{
      MAdegen=reacJMAdegenerate(AM,n,m,debug);
      addnewto1Darray(outvec, totv, MAdegen);
      totv=addv2(totv, (char*)"MA equilibria always degenerate?", outlist);
    }
  }

  if(!filter || strstr(filter, "endo")){//2 means strongly endotactic
    endotact=(int)(endotactic(AM,n,m,debug));
    addnewto1Darray(outvec, totv, endotact);
    totv=addv2(totv, (char*)"endotactic", outlist);
  }


  //accordance, concordance, etc
  if(!filter || strstr(filter, "accord") || strstr(filter, "concord")){
    //accordant
    //  if(DSRCstar){acc=1;}else{acc=accord(AM,n,m,debug);}
    acc=accord(AM,n,m,debug);
    addnewto1Darray(outvec, totv, acc);
    totv=addv2(totv, (char*)"accordant", outlist);
    //semiaccordant
    if(acc){semiacc=1;}else{semiacc=semiaccord(AM,n,m,debug);}
    addnewto1Darray(outvec, totv, semiacc);
    totv=addv2(totv, (char*)"semiaccordant", outlist);
    //MA eq accordant: J is P0- at MA equilibria
    if(!nontriv){
      MAeqacc=-1;
      addnewto1Darray(outvec, totv, MAeqacc);
      totv=addv2(totv, (char*)"MAeq accordant test not relevant (dynamically trivial network)", outlist);
    }
    else{
      if(semiacc){MAeqacc=1;}else{MAeqacc=MAeqaccord(AM,n,m,maxpppdeg,debug);}
      addnewto1Darray(outvec, totv, MAeqacc);
      totv=addv2(totv, (char*)"MAeq accordant (MA Jacobian at equilibria is P0-)", outlist);
    }


    //concordant
    if(acc && norm){conc=1;}else{conc=concord(AM,n,m,debug);}
    addnewto1Darray(outvec, totv, conc);
    totv=addv2(totv, (char*)"concordant", outlist);
    //semiconcordant
    if(conc){semiconc=1;}else{semiconc=semiconcord(AM,n,m,debug);}
    addnewto1Darray(outvec, totv, semiconc);
    totv=addv2(totv, (char*)"semiconcordant", outlist); 
    //MAeq concordant
    if(!nontriv){
      MAeqconc=-1;
      addnewto1Darray(outvec, totv, MAeqconc);
      totv=addv2(totv, (char*)"MAeq concordant test not relevant (dynamically trivial network)", outlist);
    }
    else{
      if(semiconc){MAeqconc=1;}else{MAeqconc=reacMAeqconcord(AM,n,m,maxpppdeg,debug);}
      addnewto1Darray(outvec, totv, MAeqconc);
      totv=addv2(totv, (char*)"MAeqconcordant (MA Jacobian has maximal rank at equilibria?)", outlist);
    }
  }

  //Connected with the possibility of Hopf bifurcation, stability of equilibria, etc
  if(!filter || strstr(filter, "hopf") || strstr(filter, "Hopf") || strstr(filter, "stable")){
    //Are hopf bifurcations forbidden?
    if(!nontriv){nohopf=3;}else{nohopf=HopfForbid(AM,n,m,"ALL",maxpppdeg+1,&pppdegused,debug);}
    addnewto1Darray(outvec, totv, nohopf);
    totv=addv2(totv, (char*)"Hopf bifurcations forbidden?", outlist);
    //QMAnegsemidef
    if(!nontriv){
      QMAnegsemidef=-1;
      addnewto1Darray(outvec, totv, QMAnegsemidef);
      totv=addv2(totv, (char*)"QMAnegsemidef test not relevant (dynamically trivial network)", outlist);
    }
    else{
      QMAnegsemidef=reacQMAposdef(AM,n,m,1,0,maxpppdeg,debug);
      addnewto1Darray(outvec, totv, QMAnegsemidef);
      totv=addv2(totv, (char*)"QMAnegsemidef (MA Jacobian at equilibria in the closure of Hurwitz matrices)", outlist);
    }
    //Does the system with general kinetics definitely admit a nonzero imaginary eigenvalue?
    if(nohopf>=2 || n<=1 || n>=5){ipairGK=0;}else{ipairGK=reacJadmitsIpair(AM,n,m,maxpppdeg,&pppdegused,debug);}
    addnewto1Darray(outvec, totv, ipairGK);
    totv=addv2(totv, (char*)"Definitely admits a nonzero imaginary eigenvalue (general kinetics)?", outlist);
    //Does the system with MA kinetics definitely admit an imaginary eigenvalue at some equilibrium?
    if(nohopf>=1 || n<=1 || n>=5){ipairMA=0;}else{ipairMA=reacJMAadmitsIpair(AM,n,m,maxpppdeg,&pppdegused,debug);}
    addnewto1Darray(outvec, totv, ipairMA);
    totv=addv2(totv, (char*)"Definitely admits a nonzero imaginary eigenvalue at an equilibrium (MA kinetics)?", outlist);
  }


  if(debug){fprintf(stderr, "\n########Exiting reacreport: carried out %d tests.\n___________________________\n", totv);}
  return totv;
} 

//wrapper to produce report on a single reaction in any format
//type=1: di6 input; type=2: reacstr; type 3: Sauro; type 4: simpstr
//types 2 and 3 do not use n1 and m1: these are then dummy variables
int reacreport(char *str, int type, int n1, int m1, int maxpppdeg, char ***outlist, int **outvec, char *filter, int debug){
  int entries,ret,Srank;
  int **AM;
  int **S, **Sl, **Sr;
  bool minus=1;
  int n=n1,m=m1;//for reacstr and Sauro, n1,m1 are dummies
  char **chems;
  int totV;

  if(type==1){//di6
    AM=di6toCRNam1(str, n, m, &totV, &entries);
    AMtoSSl(AM, n, m, minus, &S, &Sl);
  }
  else if(type==2){//reacstr
    getallreacs(str, &S, &Sl, &Sr, &chems, &n, &m);
    free_imatrix(Sr, 0, n-1, 0, m-1);//not needed
    freearraydat(chems,n);
    AM=SSltoAM(S,Sl, n, m, 0);//not minus
  }
  else if(type==3){//Sauro
    AM=SaurotoAM(str, &n, &m);
    AMtoSSl(AM, n, m, minus, &S, &Sl);
  }
  else if(type==4){//simpstr
    AM=simpstrtoAM(str, n, m);
    AMtoSSl(AM, n, m, minus, &S, &Sl);
  }
  else{
    fprintf(stderr, "ERROR in \"reacreport\": input type not recognised. EXITING.\n");exit(0);
  }
  Srank=matrank(S,n,m);
  ret=reacreport(AM,S,Sl,n,m,Srank,maxpppdeg,outlist,outvec,filter,debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//This is the master wrapper program to read each reaction from a file
//and produce the report
unsigned long reacreportfile(char *fname, char *reactype, char *filter, int n, int m, int maxpppdeg, int debug){
  FILE *fd;
  int i,maxl=0;
  numflines(fname, &maxl);//will warn and exit if can't be read
  unsigned long tot=0;
  char line[maxl];
  int maxblock=300*maxl;
  char block[maxblock];
  int type;
  bool isblock=0;
  int reacstrlen;
  char **outlist;
  int *outvec;
  int totreport;

  if(!reactype && ((type=guessreacformat(fname))<=0 || type >=5)){
    fprintf(stderr, "ERROR in \"reacfileread\". Please specify the reaction format. EXITING.\n");
    exit(0);
  }
  else if(!strcmp(reactype, "di6") || !strcmp(reactype, "d6"))
    type=1;
  else if(!strcmp(reactype, "reacstr") || !strcmp(reactype, "reacs"))
    type=2;
  else if(!strcmp(reactype, "Sauro") || !strcmp(reactype, "sauro"))
    type=3;
  else if(!strcmp(reactype, "simp") || !strcmp(reactype, "simpstr"))
    type=4;
  else{
    fprintf(stderr, "ERROR in \"reacfileread\". Please specify the reaction format. EXITING.\n");
    exit(0);
  }


  if(type==2)
    isblock=1;
  
  fd=fopen(fname, "r");

  while((reacstrlen=reacread(fd, isblock, line, maxl, block, maxblock))>0){
    tot++;
    printf("%ld: report\n", tot);
    if(isblock)
      totreport=reacreport(block, type, n, m, maxpppdeg, &outlist, &outvec, filter, debug);
    else
      totreport=reacreport(line, type, n, m, maxpppdeg, &outlist, &outvec, filter, debug);
    //do something
    for(i=0;i<totreport;i++)
      printf("%s: %d\n",outlist[i],outvec[i]);
    printf("******\n");

    freearraydat(outlist, totreport);
    free((char*)outvec);
  }

  fclose(fd);
  return tot;
}

//Facial structure of the reactant vertex polytope of CRN
void getfacestructCRN(int **AM, int n, int m, int debug){
  int **S, **Sl, **Sl1;
  int m1;
  bool **allfacets;
  int numfacets;
  int pat[m];
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  Sl1=redmat1(Sl,n,m,&m1, pat);//remove repeated cols
  bool verts[m1];
  numfacets=getfacestruct(Sl1, n, m1, verts, &allfacets, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  free_imatrix(Sl1, 0, n-1, 0, m1-1);
  free_bmat(allfacets,numfacets);
  return;
}

//overloading di6 input
void getfacestructCRN(char *di6, int n, int m, int debug){
  int entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  getfacestructCRN(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}

bool facetisparallel(int **Splus, int n, int m, int Srk, int **vertmat, int numverts, bool *facet){
  int i,j;
  int newrk;
  int sz;
  long *inds=booltoinds(facet,numverts,&sz);
  for(i=0;i<n;i++){
    for(j=0;j<sz;j++){
      Splus[i][m+j]=vertmat[i][inds[j]]-vertmat[i][inds[0]];
    }
  }

  newrk=matrank(Splus,n,m+sz);//rank of CRN

  if(newrk==Srk)//facet is parallel
    return 1;
  return 0;
}


//partition a CRN
//assume we have already checked that rank(S)< rank(N) where N=Newton polytope
int *partitionCRN(int **S, int **Sl, int n, int m, int *numpart){
  int j, ind=0, tot=0;
  int *part=(int *)malloc((size_t) ((m)*sizeof(int)));
  int **St=transposemat(S,n,m);
  int **Slt=transposemat(Sl,n,m);
  int **Sltshift=imatrix(0, m-1, 0, n-1);
  inittozero(part,m);

  while(tot<m){//counts total
    ind++;
    for(j=0;j<m;j++){//first remaining zero
      if(part[j]==0){
	part[j]=ind;
	tot++;
	break;
      }
    }
    translatematbyrow(Slt, m, n, j, Sltshift);
    for(j=0;j<m;j++){
      if(part[j]==0 && isincone(St,(long*)NULL,m,n,Sltshift[j])){
	part[j]=ind;tot++;
      }
    }
  }

  free_imatrix(St,0,m-1,0,n-1);
  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(Sltshift,0,m-1,0,n-1);
  (*numpart)=ind;
  return part;
}

//Extract a subCRN
int subCRN(int **S, int **Sl, int n, int m, int ***subS, int ***subSl, int *part, int ind){
  int i,j,m1=0;
  for(j=0;j<m;j++){
    if(part[j]==ind)
      m1++;
  }
  (*subS)=imatrix(0,n-1,0,m1-1);
  (*subSl)=imatrix(0,n-1,0,m1-1);

  m1=0;
  for(j=0;j<m;j++){
    if(part[j]==ind){
      for(i=0;i<n;i++){
	(*subS)[i][m1]=S[i][j];
	(*subSl)[i][m1]=Sl[i][j];
      }
      m1++;
    }
  }
  return m1;
}


//Check if a CRN is "weakly endotactic"
//i.e., no reaction vectors point out of N
//a necessary but not sufficient condition for endotacticity
int weakendo(int **S, int **Sl, int n, int m, int debug){
  int j, j0=-1, ret=1;
  int **St=transposemat(S,n,m);
  int **Slt=transposemat(Sl,n,m);
  int **Sltshift=imatrix(0, m-1, 0, n-1);
  if(debug){fprintf(stderr, "\n###Entering weakendo.\n");}

  for(j=0;j<m;j++){
    translatematbyrow(Slt, m, n, j, Sltshift);
    if(!isincone(Sltshift, (long*)NULL, m, n, St[j])){
      j0=j+1;ret=0;break;
    }
  }

  free_imatrix(St,0,m-1,0,n-1);
  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(Sltshift,0,m-1,0,n-1);

  if(debug){
    if(!ret)
      fprintf(stderr, "An outward pointing reaction (%d) found.\n", j0);
    else
      fprintf(stderr, "No outward pointing reactions found.\n");
  }
  return ret;
}



//Check if a (sub)CRN is endotactic/strongly endotactic;
//Recursive if there are inessential faces
//returns 2 for strong, 1 for endo but not strong, and 0 for not endo
//Assume we have broken into components if necessary using "partitionCRN" so that
//the dimension of the Newton polytope is <= the rank of the network
int endotactic_part(int **S, int **Sl, int n, int m, int debug){
  int i, i1, j, j1, j2, m1, numfacets;
  int numverts;
  //facets of dimension 2 to rk-1
  bool **allfacets=NULL;
  int rank, *rankvec;
  int **Sl1, **vertmat;
  int **St1=imatrix(0, m-1, 0, n);
  int **Slt1=imatrix(0, m-1, 0, n);
  int pat[m];
  int ret=2;//strong endotactic until proven otherwise
  int goodreac;
  Sl1=redmat1(Sl,n,m,&m1,pat);//remove repeated cols (store the pattern)
  bool verts[m1];
  int totverts;
  int debugfull=(debug<=0)?0:debug-1;
  int **redS, **redSl;
  bool iness[m];//to store inessential face
  int totiness;
  int subret;
  char *str;

  if(debug){fprintf(stderr, "\n###Entering endotactic_part(...)\n");}

  if(!weakendo(S,Sl, n, m, debug)){
    if(debug){fprintf(stderr, "The (sub)CRN fails to be weakly endotactic.\n");}
    return 0;
  }

  //alternative version (small networks): use getfacestruct instead of getfacestruct1; comment out "free((char*)rankvec);"
  //vectors in allfacets have length numverts, and refer to vertmat
  //numfacets=getfacestruct(Sl1, n, m1, verts, &allfacets, debugfull);
  numfacets=getfacestruct1(Sl1, n, m1, verts, &allfacets, &rank, &rankvec, &totverts, debugfull);

  //matrix of extremal vectors of the cone
  vertmat=conemat(Sl1, n, m1, verts, &numverts);//dimensions (n+1) X numverts

  if(debug){fprintf(stderr, "Vertex matrix:\n");printmat(vertmat,n+1,numverts);}
  if(debug){fprintf(stderr, "Found %d proper faces of dimension >=2\n", numfacets);}

  //transposed and augmented matrices
  for(j=0;j<m;j++){
    for(i=0;i<n;i++){
      St1[j][i]=S[i][j];Slt1[j][i]=Sl[i][j];
    }
    St1[j][n]=0;Slt1[j][n]=1;
  }
  //if(debug){fprintf(stderr, "Augmented transposed matrices: \n");printmat(St1,m,n+1);printmat(Slt1,m,n+1);}
  //check for inessential faces (we have already checked not weakly endotactic)
  for(i=0;i<numfacets;i++){//each face
    if(debug){fprintf(stderr, "Checking face %d: \n", i+1);printvec(allfacets[i],numverts);}
    goodreac=0;inittozero(iness,m);totiness=0;
    for(j=0;j<m;j++){//find each reactant complex in face
      if(isinrcone(vertmat, allfacets[i], n+1, numverts, Slt1[j])){
	iness[j]=1;totiness++;//keep track of reactions in case face is inessential
	if(debug){fprintf(stderr, "\tReaction %d originates in the face at: ", j+1);printvec(Slt1[j],n);}
	//check if reaction vector is in span of face
	if(!isinrspan(vertmat, allfacets[i], n+1, numverts, St1[j], 0)){
	  goodreac=1;//reaction pointing off face: not inessential
	  if(debug){fprintf(stderr, "\t\tReaction %d points off the face.\n",j+1);}
	  break;
	}
	else{
	  if(debug){fprintf(stderr, "\t\tReaction %d points along the face.\n", j+1);}
	}
      }
    }
    if(!goodreac){//inessential face with no good reactions
      if(debug){fprintf(stderr, "**The face has no good reactions. Removing and trying again.\n\n");}
      ret=1;//fails strong endo, but could still be endo
      //remove everything stored in iness to create redS, redSl
      redS=imatrix(0,n-1,0,m-totiness-1);redSl=imatrix(0,n-1,0,m-totiness-1);
      j2=0;
      for(j1=0;j1<m;j1++){
	if(!(iness[j1])){
	  for(i1=0;i1<n;i1++){
	    redS[i1][j2]=S[i1][j1];redSl[i1][j2]=Sl[i1][j1];
	  }
	  j2++;
	}
      }
      if(debug){str=SSltostr(redS, redSl, n, m-totiness);fprintf(stderr, "New subCRN:\n%s\n",str);free(str);printmat(redS,n,m-totiness);printmat(redSl,n,m-totiness);}
      subret=endotactic_part(redS, redSl, n, m-totiness, debug);
      free_imatrix(redS,0,n-1,0,m-totiness-1);free_imatrix(redSl,0,n-1,0,m-totiness-1);
      if(subret<ret){ret=subret;}
      if(ret==0)//failure
	break;
    }
    else if(debug){fprintf(stderr, "**The face has a good reaction.\n\n");}
  }

  free_imatrix(Sl1, 0, n-1, 0, m1-1);
  free_imatrix(St1, 0, m-1, 0, n);
  free_imatrix(Slt1, 0, m-1, 0, n);
  free_imatrix(vertmat, 0, n, 0, numverts-1);
  free_bmat(allfacets, numfacets);
  if(numfacets)
    free((char*)rankvec);

  /* if(debug) */
  /*     fprintf(stderr, "The (sub)network %s strongly endotactic.\n",ret?"is":"fails to be"); */

  return ret;
}



int endotactic(int **S, int **Sl, int n, int m, int debug){
  int i, m1, ret=2, subret;
  int Srk=matrank(S,n,m);//rank of CRN
  int Slrk;
  bool verts[m];
  int **SlC;
  int *part;
  int numpart;
  int sz;
  int **subS,**subSl;
  char *str;
  inittoone(verts,m);
  SlC=conemat(Sl,n,m,verts,&m1);
  Slrk=matrank(SlC,n+1,m1)-1;
  free_imatrix(SlC,0,n,0,m1-1);
  //Only decompose if Srk< Slrk
  if(Srk<Slrk){
    part=partitionCRN(S, Sl, n, m, &numpart);//partition
    if(debug){fprintf(stderr, "partitioning network: ");printvec(part,m);}
    for(i=0;i<numpart;i++){
      sz=subCRN(S, Sl, n, m, &subS, &subSl, part, i+1);
      if(debug){str=SSltostr(subS, subSl, n, sz);fprintf(stderr, "processing subCRN %d: \n%s\n", i+1, str);free(str);}

      //process each part
      subret=endotactic_part(subS, subSl, n, sz, debug);
      free_imatrix(subS,0,n-1,0,sz-1);
      free_imatrix(subSl,0,n-1,0,sz-1);
      if(subret<ret){ret=subret;}
      if(ret==0){
	free((char*)part);
	if(debug){fprintf(stderr, "The subCRN fails to be endotactic.\n#####\n");}
	return 0;
      }
      if(debug){fprintf(stderr, "The subCRN is %sendotactic.\n#####\n",ret==2?"strongly ":"");}

    }
    free((char*)part);
  }
  else//process whole
    ret=endotactic_part(S, Sl, n, m, debug);

  if(debug){
    if(ret==2){fprintf(stderr, "The CRN is strongly endotactic.\n");}
    else if(ret==1){fprintf(stderr, "The CRN is endotactic but not strongly endotactic.\n");}
    else{fprintf(stderr, "The CRN fails to be endotactic.\n");}
  }

  return ret;
}


//overloading: PN input
int endotactic(int **AM, int n, int m, int debug){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=endotactic(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input
int endotactic(char *di6, int n, int m, int debug){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=endotactic(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//overloading: PN input
int PolytopeVol(int **AM, int n, int m, int debug){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=PolytopeVol(S, Sl, n, m, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}


//overloading: di6 input (rank is ignored here)
int PolytopeVol(char *di6, int n, int m, int debug){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=PolytopeVol(AM,n,m,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  fprintf(stderr, "Vol=%d\n", ret);
  return ret;
}



//largely untested and very inefficient
//computes the mixed volume in the case with and without
//conservation laws. 
int MixedVolCRN(int **S, int **Sl, int n, int m, int sourceonly, int rank, int debug){
  int i,j,k;
  bool **b;
  int **Sl1, **Sl2;
  int nummin,V,Vf;
  int xc[n];
  int flag;
  int PV;
  int **St,**basisSt;
  int tot,rk,deg;
  int pat[n];
  int lins;
  
  if(sourceonly){//sources only
    if(rank==n){//full rank case
      V=PolytopeVol(S, Sl, n, m, debug);
      fprintf(stderr, "MixedVol=%d\n",V);
      return V;
    }
    else{
      //sources only: r copies of Newton polytope and n-r copies of
      //basic-simplex
      
      //add in basic sources 0, x_i

      Sl1=imatrix(0, n-1, 0, m+n);
      cpmat(Sl, Sl1, n, m);
      for(i=0;i<n;i++){
	for(j=m;j<m+n+1;j++){
	  if(j==m)
	    Sl1[i][j]=0;
	  else{
	    if(i==j-m-1)
	      Sl1[i][j]=1;
	    else
	      Sl1[i][j]=0;
	  }
	}
      }
      b=bmatrix(0, n-1, 0, m+n);
      inittozero(b,n,m+n+1);
      for(i=0;i<rank;i++){//rank copies of Newton polytope
	for(j=0;j<m;j++){b[i][j]=1;}
      }
      for(i=rank;i<n;i++){//n-rank copies of basic simplex
	for(j=m;j<m+n+1;j++){b[i][j]=1;}
      }
      
      //printmat(Sl1,n,m+n+1);//
      //printmat(b,n,m+n+1);//

      for(i=0;i<n;i++){xc[i]=i;}
      Sl2=MinkowskiN(Sl1, n, m+n+1, b, xc, n, &nummin);
      //printmat(Sl2,n,nummin);//
      V=PolytopeVol(NULL, Sl2, n, nummin, debug);
      free_imatrix(Sl2, 0, n-1, 0, nummin-1);
      //fprintf(stderr, "PV=%d\n",V);//
      if(debug){fprintf(stderr, "MixedVol (running) = %d/%d\n", V, (int)(factorial(n)));}
      for(k=1;k<n;k++){
	firstcomb(xc, n, k);
	flag=1;
	while(flag==1){
	  //fprintf(stderr, "checking: ");printvec(xc,k);//
	  Sl2=MinkowskiN(Sl1, n, m+n+1, b, xc, k, &nummin);
	  //printmat(Sl2,n,nummin);//
	  PV=0;if(nummin>=n){PV=PolytopeVol(NULL, Sl2, n, nummin, debug);}
	  //fprintf(stderr, "PV=%d\n",PV);//
	  V+=(int)(pow(-1.0,(double(n-k))))*PV;
	  if(debug){fprintf(stderr, "MixedVol (running) = %d/%d\n", V, (int)(factorial(n)));}//exit(0);
	  //printmat(Sl2,n,nummin);
	  free_imatrix(Sl2, 0, n-1, 0, nummin-1);
	  flag=nextcomb(xc, n, k);
	  //fprintf(stderr, "V=%d/%d\n",V,(int)factorial(n));
	}
      }
      free_imatrix(Sl1, 0, n-1, 0, m+n);
      if(V%((int)factorial(n))!=0){
	fprintf(stderr, "Oops: something went wrong with the Mixed Volume calculation: not an integer. EXITING.\n");exit(0);
      }
      Vf=V/((int)factorial(n));//should be an integer
      fprintf(stderr, "MixedVol=%d\n",Vf);
      free_bmatrix(b,0,n-1,0,m+n);
      return Vf;
    }//end of not-full-rank case
  }//end of sources only

  //
  // Full CRN (i.e., not sources only)
  //

  //conservation laws
  if(matrank(S,n,m)<n){

    St=transposemat(S,n,m);
    basisSt=minkerbasis(St,m,n,&tot,&rk,&deg,debug);
    inittoone(pat,n);
    for(i=0;i<tot;i++){//each basis vector
      j=0;
      //Choose which species to determine from conservation laws
      //printvec(basisSt[j],n);
      while(j<n && (!(basisSt[i][j]) || !pat[j])){j++;}
	pat[j]=0;
    }
    //fprintf(stderr, "pat: "); printvec(pat,n);


    //augment with basic simplex
    Sl1=imatrix(0, n-1, 0, m+n);
    cpmat(Sl, Sl1, n, m);
    for(i=0;i<n;i++){
      for(j=m;j<m+n+1;j++){
	if(j==m)
	  Sl1[i][j]=0;
	else{
	  if(i==j-m-1)
	    Sl1[i][j]=1;
	  else
	    Sl1[i][j]=0;
	}
      }
    }
    b=bmatrix(0, n-1, 0, m+n);
    inittozero(b,n,m+n);
    lins=0;
    for(i=0;i<n;i++){
      if(pat[i]){//reaction poly
	for(j=0;j<m;j++){
	  if(S[i][j])//net production or consumption in jth reaction
	    b[i][j]=1;
	}
      }
      else{//conservation poly
	for(j=m;j<n+m+1;j++){
	  //printvec(basisSt[lins],n);
	  if(j==m || basisSt[lins][j-m-1])
	    b[i][j]=1;
	}
	lins++;
      }
    }
    printmat(Sl1,n,n+m+1);
    printmat(b,n,n+m+1);

    //Full volume
    for(i=0;i<n;i++){xc[i]=i;}
    Sl2=MinkowskiN(Sl1, n, m+n+1, b, xc, n, &nummin);
    V=PolytopeVol(NULL, Sl2, n, nummin, debug);
    free_imatrix(Sl2, 0, n-1, 0, nummin-1);
    //fprintf(stderr, "V=%d/%d\n",V,(int)factorial(n));exit(0);

    for(k=1;k<n;k++){
      firstcomb(xc, n, k);
      flag=1;
      while(flag==1){
	Sl2=MinkowskiN(Sl1, n, m+n+1, b, xc, k, &nummin);
	V+=(int)(pow(-1.0,(double(n-k))))*PolytopeVol(NULL, Sl2, n, nummin, debug);
	//printmat(Sl2,n,nummin);
	free_imatrix(Sl2, 0, n-1, 0, nummin-1);
	flag=nextcomb(xc, n, k);
	//fprintf(stderr, "V=%d/%d\n",V,(int)factorial(n));
      }
    }
    if(V%((int)factorial(n))!=0){
      fprintf(stderr, "Oops: something went wrong with the Mixed Volume calculation: not an integer. EXITING.\n");exit(0);
    }
    Vf=V/((int)factorial(n));//should be an integer
    fprintf(stderr, "MixedVol=%d\n",Vf);
    free_bmatrix(b,0,n-1,0,m+n-1);
    free_imatrix(St,0,m-1,0,n-1);
    free_imat(basisSt,tot);
    free_imatrix(Sl1, 0, n-1, 0, m+n);
  }
  else{//no conservation laws

    b=bmatrix(0, n-1, 0, m-1);
    inittozero(b,n,m);
    for(i=0;i<n;i++){
      for(j=0;j<m;j++){
	if(S[i][j])//net production or consumption in jth reaction
	  b[i][j]=1;
      }
    }

    //Full volume
    for(i=0;i<n;i++){xc[i]=i;}
    Sl2=MinkowskiN(Sl, n, m, b, xc, n, &nummin);
    V=PolytopeVol(NULL, Sl2, n, nummin, debug);
    free_imatrix(Sl2, 0, n-1, 0, nummin-1);
    //fprintf(stderr, "V=%d/%d\n",V,(int)factorial(n));

    for(k=1;k<n;k++){
      firstcomb(xc, n, k);
      flag=1;
      while(flag==1){
	Sl2=MinkowskiN(Sl, n, m, b, xc, k, &nummin);
	V+=(int)(pow(-1.0,(double(n-k))))*PolytopeVol(NULL, Sl2, n, nummin, debug);
	//printmat(Sl2,n,nummin);
	free_imatrix(Sl2, 0, n-1, 0, nummin-1);
	flag=nextcomb(xc, n, k);
	//fprintf(stderr, "V=%d/%d\n",V,(int)factorial(n));
      }
    }
    if(V%((int)factorial(n))!=0){
      fprintf(stderr, "Oops: something went wrong with the Mixed Volume calculation: not an integer. EXITING.\n");exit(0);
    }
    Vf=V/((int)factorial(n));//should be an integer
    fprintf(stderr, "MixedVol=%d\n",Vf);
    free_bmatrix(b,0,n-1,0,m-1);

  }
  return Vf;

}

//overloading: PN input
int MixedVolCRN(int **AM, int n, int m, int sourceonly, int rank, int debug){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=MixedVolCRN(S, Sl, n, m, sourceonly, rank, debug);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input (rank is ignored here)
int MixedVolCRN(char *di6, int n, int m, int rank, int debug){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=MixedVolCRN(AM,n,m,0,rank,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//overloading: di6 input (rank essential here)
int MixedVolCRNsrc(char *di6, int n, int m, int rank, int debug){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=MixedVolCRN(AM,n,m,1,rank,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Solvability and Mixed Vol comparison (for (n,n+k,n) CRNs)
void printsolvabilityMixedVol(char *di6, int n, int m, int debug){
  int entries,degree,V,totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  degree=solvabilitydegree(AM,n,m);
  V=MixedVolCRN(AM,n,m,0,-1,debug);
  if(degree>V)
    fprintf(stderr, "****");
  else if(V>degree)
    fprintf(stderr, "####");
  fprintf(stderr, "solvability: %d, Mixed Vol: %d\n", degree, V);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}

//Solvability and Mixed Vol using sources only comparison (for (n,n+k,n) CRNs)
void printsolvabilityMixedVolsrc(char *di6, int n, int m, int rank, int debug){
  int entries,degree,V,totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  degree=solvabilitydegree(AM,n,m);
  V=MixedVolCRN(AM,n,m,1,rank,debug);
  if(degree>V)
    fprintf(stderr, "****");
  else if(V>degree)
    fprintf(stderr, "####");
  fprintf(stderr, "solvability: %d, Mixed Vol (sources only): %d\n", degree, V);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}

//Q-mixed volume and Mixed Vol using sources only comparison (for (n,m,m-1) CRNs)
void printQmixedVolsrc(char *di6, int n, int m, int rank, int debug){
  int entries,V,totV, QM;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  QM=Qmixed(AM,n,m);
  V=MixedVolCRN(AM,n,m,1,rank,debug);
  if(QM>V)
    fprintf(stderr, "****");
  else if(V>QM)
    fprintf(stderr, "####");
  printQmixed(AM,n,m);
  fprintf(stderr, "Q-mixed volume: %d, Mixed Vol (sources only): %d\n", QM, V);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return;
}



//Check if a CRN has any complex with three or more reactions and such 
//that one of these lies in the span of the others. In this case, 
//dynamical equivalence is not equivalent to simple equivalence. 
//Apply this test after removing CRNs with redundant reactions
int spanreacs(int **S, int **Sl, int n, int m){
  int i,j,i1;
  int **Slt=transposemat(Sl, n, m);
  int **St=transposemat(S, n, m);
  int part[m];
  bool partb[m];
  int ind=1;
  int nR;
  inittozero(part,m);
  //partition the reactant complexes (R-Cs)
  for(i=0;i<m;i++){
    if(!part[i]){
      part[i]=ind++;
      for(j=i+1;j<m;j++){
	if(!part[j] && areequal(Slt[i],Slt[j],n))
	  part[j]=part[i];
      }
    }
  }

  //Go through each R-C
  for(i1=1;i1<ind;i1++){
    nR=0;
    for(i=0;i<m;i++){//set boolean vector
      if(part[i]==i1){partb[i]=1;nR++;}else{partb[i]=0;}
    }
    if(nR>=3){//more than two reactions on this R-C
      for(i=0;i<m;i++){
	if(partb[i]){
	  partb[i]=0;
	  //printmat(St,m,n);printvec(St[i],n);
	  if(isinspan(St,partb,m,n,St[i],0)){//redundant reac.
	    free_imatrix(Slt,0,m-1,0,n-1);
	    free_imatrix(St,0,m-1,0,n-1);
	    return 1;
	  }
	  partb[i]=1;
	}
      }
    }
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(St,0,m-1,0,n-1);
  return 0;
}

//overloading: PN input
int spanreacs(int **AM, int n, int m){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=spanreacs(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input
int spanreacs(char *di6, int n, int m){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=spanreacs(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//affinely dependent sources
//implies equilibria are identically degenerate in (2,n,2) networks
int affdepsources(int **S, int **Sl, int n, int m){
  int j;
  int ret=0;
  int **Slp=imatrix(0, n, 0, m-1);
  cpmat(Sl, Slp, n, m);
  for(j=0;j<m;j++)
    Slp[n][j]=1;

  if(matrank(Slp,n+1,m)<n+1)
    ret=1;
  free_imatrix(Slp,0,n,0,m-1);

  return ret;
}

//overloading: PN input
int affdepsources(int **AM, int n, int m){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=affdepsources(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input
int affdepsources(char *di6, int n, int m){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=affdepsources(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Count the number of distinct source complexes in a CRN
int countsources(int **S, int **Sl, int n, int m){
  int i,j;
  int **Slt=transposemat(Sl, n, m);
  int part[m];
  int ind=1;
  inittozero(part,m);
  //partition the reactant complexes (R-Cs)
  for(i=0;i<m;i++){
    if(!part[i]){
      part[i]=ind++;
      for(j=i+1;j<m;j++){
	if(!part[j] && areequal(Slt[i],Slt[j],n))
	  part[j]=part[i];
      }
    }
  }
  free_imatrix(Slt,0,m-1,0,n-1);
  return ind-1;
}

//overloading: PN input
int countsources(int **AM, int n, int m){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=countsources(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input, l layers
int countsources(char *di6, int n, int m){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=countsources(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//2 layers only - deprecated
int countsourcesold(char *di6, int n, int m){
  int ret, entries;
  int **AM=di6toCRNam(di6, n, m, &entries);
  ret=countsources(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}



//Check if a CRN has any reactions which could be removed to leave
//a dynamically equivalent CRN
int redundantreacs(int **S, int **Sl, int n, int m){
  int i,j,i1;
  int **Slt=transposemat(Sl, n, m);
  int **St=transposemat(S, n, m);
  int part[m];
  bool partb[m];
  int ind=1;
  int nR;
  inittozero(part,m);
  //partition the reactant complexes (R-Cs)
  for(i=0;i<m;i++){
    if(!part[i]){
      part[i]=ind++;
      for(j=i+1;j<m;j++){
	if(!part[j] && areequal(Slt[i],Slt[j],n))
	  part[j]=part[i];
      }
    }
  }

  //Go through each R-C
  for(i1=1;i1<ind;i1++){
    nR=0;
    for(i=0;i<m;i++){//set boolean vector
      if(part[i]==i1){partb[i]=1;nR++;}else{partb[i]=0;}
    }
    if(nR>=2){//more than one reaction on this R-C
      for(i=0;i<m;i++){
	if(partb[i]){
	  partb[i]=0;
	  //printmat(St,m,n);printvec(St[i],n);
	  if(isincone(St,partb,m,n,St[i])){//redundant reac.
	    free_imatrix(Slt,0,m-1,0,n-1);
	    free_imatrix(St,0,m-1,0,n-1);
	    return 1;
	  }
	  partb[i]=1;
	}
      }
    }
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(St,0,m-1,0,n-1);
  return 0;
}



//overloading: PN input
int redundantreacs(int **AM, int n, int m){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=redundantreacs(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input
int redundantreacs(char *di6, int n, int m){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=redundantreacs(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}

//Does the CRN have repeated reactions?
int repeatedreacs(int **S, int **Sl, int n, int m){
  int i,j;
  int **Slt=transposemat(Sl, n, m);
  int **St=transposemat(S, n, m);
  for(i=0;i<m-1;i++){
    for(j=i+1;j<m;j++){
      if(areequal(Slt[i],Slt[j],n) && areequal(St[i],St[j],n)){
	free_imatrix(Slt,0,m-1,0,n-1);
	free_imatrix(St,0,m-1,0,n-1);
	return 1;
      }

    }
  }

  free_imatrix(Slt,0,m-1,0,n-1);
  free_imatrix(St,0,m-1,0,n-1);
  return 0;
}

//overloading: PN input
int repeatedreacs(int **AM, int n, int m){
  int ret, **S, **Sl;
  AMtoSSl(AM, n, m, 0, &S, &Sl);
  ret=repeatedreacs(S, Sl, n, m);
  free_imatrix(S, 0, n-1, 0, m-1);
  free_imatrix(Sl, 0, n-1, 0, m-1);
  return ret;
}

//overloading: di6 input
int repeatedreacs(char *di6, int n, int m){
  int ret, entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=repeatedreacs(AM,n,m);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}




// filter CRNs in NAUTY 2 layer format
// Runs a number of different possible tests on a CRN, determined by
// the filter
unsigned long filterCRNs(const char *fname, const char *outfname, int n, int m, const char filt[], int iconst, int debug){
  unsigned long num=0;
  int report;
  unsigned long i=0;
  FILE *fd, *fdin;
  int maxl=0;
  unsigned long numl=numflines(fname, &maxl);
  char oneline[maxl];
  char *str;
  int val;
  int pppdegused;

  //const char *allfilts[2]={"hello","my"};

  // Dynamically isomorphic CRNs (i.e., giving rise to the same set of
  // dynamical systems under any scaling invariant kinetics)
  // Uses NAUTY and writes the translation table to stderr
  if(!strcmp(filt, "isomorph"))
    return dynamicshortfile(fname, outfname, "nofilter", n, m, 0);
  else if(!strcmp(filt, "dynisomorphMA"))
    return dynamicshortfile(fname, outfname, "MA", n, m, 0);
  else if(!strcmp(filt, "dynisomorphGK"))
    return dynamicshortfile(fname, outfname, "GK", n, m, 0);
  else if(!strcmp(filt, "dynisomorphNN1"))
    return dynamicshortfile(fname, outfname, "NN1", n, m, 0);
  else if(!strcmp(filt, "isomorphcmplx"))//same complexes (up to isomorphism)
    return dynamicshortfile(fname, outfname, "allcomplexes", n, m, 0);
  else if(!strcmp(filt, "isomorphnewton"))//isomorphic newton polytopes
    return dynamicshortfile(fname, outfname, "sourcecomplexes", n, m, 0);
  else if(!strcmp(filt, "isomorphnewton1"))//isomorphic sources
    return dynamicshortfile(fname, outfname, "sourcecomplexesfull", n, m, 0);
   else if(!strcmp(filt, "orderbynetprod"))
    return dynamicshortfile(fname, outfname, "netprod", n, m, 0);

  if(!(fdin = fopen(fname, "r"))){
    fprintf(stderr, "ERROR in filterCRNs: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }

  if(!outfname){//pseudotest: print to stdout
    printf("checking %ld CRNs.\n", numl);
    if(!strcmp(filt, "multinn2nscr")||!strcmp(filt, "multinn3nscr")||!strcmp(filt, "multi342scr")||!strcmp(filt, "CylinderScr3")||!strcmp(filt, "CylinderScr")){
      fprintf(stderr, "#!/usr/bin/env wolframscript\n(* ::Package:: *)\n\n");
    }
    
    while(getline0(fdin, oneline, maxl) > 0){
      if(iscomline(oneline))//ignore
	continue;

      i++;
      if(!strcmp(filt, "multinn2nscr")||!strcmp(filt, "multinn3nscr")||!strcmp(filt, "multi342scr")||!strcmp(filt, "CylinderScr3")||!strcmp(filt, "CylinderScr")){
	 fprintf(stderr, "Print[\"%ld/%ld\"]\n",i,numl);
	 str=di6toreacstr(oneline,n,m,0);
	 fprintf(stderr, "Print[\"%s\"]\n",str);
	 free(str);
      }
      if(!strcmp(filt, "MixedVol")||!strcmp(filt, "MixedVolsrc")||!strcmp(filt, "printsolvabilityMixedVol")||!strcmp(filt, "printsolvabilityMixedVolsrc")||!strcmp(filt, "printQmixed")||!strcmp(filt, "printQmixedVolsrc")){
	fprintf(stderr,"//%ld/%ld\n",i,numl);//to stdout
	str=di6toreacstr(oneline,n,m,0);
	fprintf(stderr, "%s",str);
	free(str);
      }

      printf("//%ld/%ld\n",i,numl);//to stdout
      if(!strcmp(filt, "printreacs")){
	str=di6toreacstr(oneline,n,m,0);
	printf("%s", str);
	free(str);
	num++;
      }
      else if(!strcmp(filt, "linkageclasses")){
	printf("%d linkage classes\n", linkage_classes(oneline, n, m));
	num++;
      }
      else if(!strcmp(filt, "deficiency")){
	printf("deficiency: %d\n", deficiency(oneline, n, m));
	num++;
      }
      //real spectrum at MA equilibria
      else if(!strcmp(filt, "JMArealspec")){
	printreacJMArealspec(oneline, n, m, iconst-1, &pppdegused, debug);
	num++;
      }
      else if(!strcmp(filt, "printsiphons")){
	printsiphons(oneline, n, m);
	num++;
      }
      //print the MA Jacobian matrix at equilibria
      else if(!strcmp(filt, "printJ")){
	printreacJ(oneline, n, m);
	num++;
      }
      //print the MA Jacobian matrix at equilibria
      else if(!strcmp(filt, "printJMAeq")){
	printreacJMAeq(oneline, n, m);
	num++;
      }
      //print the MA Jacobian matrix at equilibria in maxima form
      else if(!strcmp(filt, "printJMAeqmax")){
	printreacJMAeq_max(oneline, n, m);
	num++;
      }
      //print the first factor in the MA Jacobian matrix at equilibria
      else if(!strcmp(filt, "printQ")){
	printreacQ(oneline, n, m);
	num++;
      }
      //print the matrix whose rows are a basis of ker[A|1]^t
      else if(!strcmp(filt, "printW")){
	printW(oneline, n, m, debug);
	num++;
      }
      //print the the stoichiometric and exponent matrices
      else if(!strcmp(filt, "printSSl")){
	printSSl(oneline, n, m);
	num++;
      }
      //print the first factor in the MA Jacobian matrix at equilibria
      //in maxima form
      else if(!strcmp(filt, "printQmax")){
	printreacQ_max(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "printmolecularity")){
	printmolecularity(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "CylinderScr3")){
	CylinderScr3(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "CylinderScr")){
	CylinderScr(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "multinn2nscr")){
	multinn2nscr(oneline, n, m, iconst);//searching for iconst PNE
	num++;
      }
      else if(!strcmp(filt, "multinn3nscr")){
	multinn3nscr(oneline, n, m, iconst);//searching for iconst PNE
	num++;
      }
      else if(!strcmp(filt, "multi342scr")){
	multi342scr(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "printA1inv")){
	printA1inv(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "printA1kervec")){
	printA1kervec(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "printsolvability")){
	printsolvability(oneline, n, m, debug);
	num++;
      }
      else if(!strcmp(filt, "printsolvabilityMixedVol")){
	printsolvabilityMixedVol(oneline, n, m, debug);
	num++;
      }
      else if(!strcmp(filt, "printsolvabilityMixedVolsrc")){
	printsolvabilityMixedVolsrc(oneline, n, m, iconst, debug);
	num++;
      }
      else if(!strcmp(filt, "printQmixed")){
	printQmixed(oneline, n, m);
	num++;
      }
      else if(!strcmp(filt, "printQmixedVolsrc")){
	printQmixedVolsrc(oneline, n, m, iconst, debug);
	num++;
      }
      else if(!strcmp(filt, "getfacestruct")){
	getfacestructCRN(oneline, n, m, debug);
	num++;
      }
      //Volume of Newton Polytope
      else if(!strcmp(filt, "PolytopeVol")){
	PolytopeVol(oneline, n, m, debug);
	num++;
      }
      //Mixed volume (4th arg: rank ignored)
      else if(!strcmp(filt, "MixedVol")){
	MixedVolCRN(oneline, n, m, -1, debug);
	num++;
      }
      //Mixed volume LHS only (4th arg: rank essential)
      else if(!strcmp(filt, "MixedVolsrc")){
	MixedVolCRNsrc(oneline, n, m, iconst, debug);
	num++;
      }
      else{
	fprintf(stderr, "\"%s\" is not a recognised pseudotest. EXITING.\n",filt);
	exit(0);
      }
      printf("******\n");
    }

    fclose(fdin);
    return num;
  }

  if(!(fd=fopen(outfname, "w"))){
    fprintf(stderr, "ERROR in filterCRNs: \"%s\" could not be opened for writing.\n", outfname);
    exit(0);
  }

  fprintf(stdout, "checking %ld CRNs.\n", numl);
  //cout << "here\n";
  while(getline0(fdin, oneline, maxl) > 0){
    if(iscomline(oneline))//ignore
      continue;

    i++;
    printf("%ld/%ld\n",i,numl);//to stdout
    if(debug){
      fprintf(stderr,"%ld/%ld\n",i,numl);
      str=di6toreacstr(oneline,n,m,0);
      fprintf(stderr, "%s", str);
      free(str);
    }

    //Start of main checks
    if(!strcmp(filt, "rank")){
      if((val=checkrank(oneline, n, m))==iconst){//rank
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "rank=%d\n", val);
    }
    else if(!strcmp(filt, "notrank")){
      if((val=checkrank(oneline, n, m))!=iconst){//rank
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "rank=%d\n", val);
    }
    //molecularity
    else if(!strcmp(filt, "molecularity")){
      if((val=molecularity(oneline, n, m))==iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "molecularity=%d\n", val);
    }
    else if(!strcmp(filt, "notmolecularity")){
      if((val=checkrank(oneline, n, m))!=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "molecularity=%d\n", val);
    }
    //number of sources
    else if(!strcmp(filt, "sources")){
      if((val=countsources(oneline, n, m))==iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "sources=%d\n", val);
    }
    else if(!strcmp(filt, "notsources")){
      if((val=countsources(oneline, n, m))!=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "sources=%d\n", val);
    }
    //at least this number of sources
    else if(!strcmp(filt, "sourcesplus")){
      if((val=countsources(oneline, n, m))>=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "sources=%d\n", val);
    }
    else if(!strcmp(filt, "notsourcesplus")){
      if((val=countsources(oneline, n, m))<iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "sources=%d\n", val);
    }
    //affine dependence among sources
    else if(!strcmp(filt, "affdepsources")){
      if(affdepsources(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notaffdepsources")){
      if(!affdepsources(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }

    //rank of reactant matrix
    else if(!strcmp(filt, "leftrank")){
      if((val=checkrankleft(oneline, n, m))==iconst){//rank
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "left rank=%d\n", val);
    }
    else if(!strcmp(filt, "notleftrank")){
      if((val=checkrankleft(oneline, n, m))!=iconst){//rank
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "left rank=%d\n", val);
    }

    //rank of reactant matrix augmented with row of ones
    else if(!strcmp(filt, "leftrank1")){
      if((val=checkrankSil1(oneline, n, m))==iconst){//rank
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "rank of [A|1]=%d\n", val);
    }
    else if(!strcmp(filt, "notleftrank1")){
      if((val=checkrankSil1(oneline, n, m))!=iconst){//rank
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "rank of [A|1]=%d\n", val);
    }

    //degree of solvability polynomial (use for (n,n+k,n))
    else if(!strcmp(filt, "solvabilitydegree")){
      if((val=solvabilitydegree(oneline, n, m))==iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "solvability degree=%d\n", val);
    }
    else if(!strcmp(filt, "notsolvabilitydegree")){
      if((val=solvabilitydegree(oneline, n, m))!=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "solvability degree=%d\n", val);
    }
    //Like solvability degree, but return CRN with solvability degree at least as given
    else if(!strcmp(filt, "solvabilitydegreeplus")){
      if((val=solvabilitydegree(oneline, n, m))>=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "solvability degree=%d\n", val);
    }
    else if(!strcmp(filt, "notsolvabilitydegreeplus")){
      if((val=solvabilitydegree(oneline, n, m))<iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "solvability degree=%d\n", val);
    }


    //mixed volume of Q
    else if(!strcmp(filt, "Qmixed")){
      if((val=Qmixed(oneline, n, m))==iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "Qmixed Vol=%d\n", val);
    }
    else if(!strcmp(filt, "notQmixed")){
      if((val=solvabilitydegree(oneline, n, m))!=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "Qmixed Vol=%d\n", val);
    }


    //dynamically nontrivial
    else if(!strcmp(filt, "DN") || !strcmp(filt, "dynnontriv")){
      if(DN(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDN") || !strcmp(filt, "notdynnontriv")){
      if(!DN(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //bounded stoichiometric classes
    else if(!strcmp(filt, "bdclass")){
      if(bdclass(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notbdclass")){
      if(!bdclass(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //homogeneous (i.e., 1^t\Gamma = 0)
    else if(!strcmp(filt, "homogen")){
      if(homogen(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "nothomogen")){
      if(!homogen(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //weakly reversible
    else if(!strcmp(filt, "WR")){
      if(CRNweakrev(oneline,n,m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notWR")){
      if(!CRNweakrev(oneline,n,m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //connected Petri Net graph
    else if(!strcmp(filt, "PNconnect")){
      if(PNconnect(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Disconnected Petri Net graph
    else if(!strcmp(filt, "notPNconnect")){
      if(!PNconnect(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Deficiency
    else if(!strcmp(filt, "deficiency")){
      if(deficiency(oneline, n, m)==iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }    
    else if(!strcmp(filt, "notdeficiency")){
      if(deficiency(oneline, n, m)!=iconst){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Weakly reversible and deficiency 0
    else if(!strcmp(filt, "WRdef0")){
      if(WRdef0(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }    
    else if(!strcmp(filt, "notWRdef0")){
      if(!WRdef0(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Genuine (no unused species)
    else if(!strcmp(filt, "genuine")){
      if(genuine(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }    
    else if(!strcmp(filt, "notgenuine")){
      if(!genuine(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //trivial species (a species with no net change in any reaction)
    else if(!strcmp(filt, "trivspec")){
      if(triv_spec(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }    
    else if(!strcmp(filt, "nottrivspec")){
      if(!triv_spec(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //non-autocatalytic (no species on both sides of any reaction)
    else if(!strcmp(filt, "nonautocat") || !strcmp(filt, "notautocat")){
      if(nonautocat(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }    
    else if(!strcmp(filt, "autocat")){
      if(!nonautocat(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //zero-one (all stoichiometries <=1)
    else if(!strcmp(filt, "zeroone")){
      if(zeroone(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    } 
    else if(!strcmp(filt, "notzeroone")){
      if(!zeroone(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //structurally persistent
    else if(!strcmp(filt, "structpersist")){
      if(structpersist(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notstructpersist")){
      if(!structpersist(oneline, n, m,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //endotactic
    else if(!strcmp(filt, "endotactic")){
      if(endotactic(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notendotactic")){
      if(!endotactic(oneline, n, m,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //strongly endotactic
    else if(!strcmp(filt, "strongendo")){
      if(endotactic(oneline, n, m, debug)==2){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notstrongendo")){
      if(endotactic(oneline, n, m,debug)!=2){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //repeated reactions?
    else if(!strcmp(filt, "repeatedreacs")){
      if(repeatedreacs(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notrepeatedreacs")){
      if(!repeatedreacs(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //redundant reactions?
    else if(!strcmp(filt, "redundantreacs")){
      if(redundantreacs(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notredundantreacs")){
      if(!redundantreacs(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //redundant reactions?
    else if(!strcmp(filt, "spanreacs")){
      if(spanreacs(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notspanreacs")){
      if(!spanreacs(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Condition * (without reversifying/merging)
    else if(!strcmp(filt, "DSRCondStar")){
      if(DSRCondStar(oneline, n, m, &report, 0, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSRCondStar")){
      if(!DSRCondStar(oneline, n, m, &report, 0, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Condition * (with reversifying/merging)
    else if(!strcmp(filt, "DSRCondStarRev")){
      if(DSRCondStar(oneline, n, m, &report, 1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSRCondStarRev")){
      if(!DSRCondStar(oneline, n, m, &report, 1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Condition * (without and then with reversifying/merging)
    else if(!strcmp(filt, "DSRCondStarBoth")){
      if(DSRCondStar(oneline, n, m, &report, 2, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSRCondStarBoth")){
      if(!DSRCondStar(oneline, n, m, &report, 2, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Condition * (without reversifying/merging)
    else if(!strcmp(filt, "DSR2CondStar")){
      if(DSR2CondStar(oneline, n, m, &report, 0, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSR2CondStar")){
      if(!DSR2CondStar(oneline, n, m, &report, 0, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Condition * (with reversifying/merging)
    else if(!strcmp(filt, "DSR2CondStarRev")){
      if(DSR2CondStar(oneline, n, m, &report, 1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSR2CondStarRev")){
      if(!DSR2CondStar(oneline, n, m, &report, 1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Condition * (without and then with reversifying/merging)
    else if(!strcmp(filt, "DSR2CondStarBoth")){
      if(DSR2CondStar(oneline, n, m, &report, 2, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSR2CondStarBoth")){
      if(!DSR2CondStar(oneline, n, m, &report, 2, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Jacobian matrix is nonsingular (general  kinetics)
    else if(!strcmp(filt, "Jnonsing")){
      if(reacJnonsing(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJnonsing")){
      if(!reacJnonsing(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Jacobian matrix is concordant (general  kinetics)
    //different from nonsing, because we don't make assumptions about 
    //the rank of the stoichiometric matrix
    //The first version is for debugging really - using symbolic algebra seems to be slower
    else if(!strcmp(filt, "concord1")){
      if(concord1(oneline, n, m, iconst-1, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notconcord1")){
      if(!concord1(oneline, n, m, iconst-1, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "concord") || !strcmp(filt, "concordant")){//do positively compatible first
      if(concord(oneline, n, m, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notconcord") || !strcmp(filt, "notconcordant")){
      if(!concord(oneline, n, m, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Semiconcordant (MA Jacobian matrix is a nonsingular transformation on the stoichiometric subspace)
    else if(!strcmp(filt, "semiconcord") || !strcmp(filt, "semiconcordant")){//do positively compatible first
      if(semiconcord(oneline, n, m, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notsemiconcord") || !strcmp(filt, "notsemiconcordant")){
      if(!semiconcord(oneline, n, m, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //accordance: minus Jacobian matrix is P0 (general  kinetics)
    else if(!strcmp(filt, "accord") || !strcmp(filt, "accordant") || !strcmp(filt, "JisP0minus")){
      if(reacJisP0(oneline, n, m, 1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notaccord") || !strcmp(filt, "notaccordant") || !strcmp(filt, "notJisP0minus")){
      if(!reacJisP0(oneline, n, m, 1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Semiaccordant (MA Jacobian matrix is P0- everywhere)
    else if(!strcmp(filt, "semiaccord") || !strcmp(filt, "semiaccordant")){//do positively compatible first
      if(semiaccord(oneline, n, m, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notsemiaccord") || !strcmp(filt, "notsemiaccordant")){
      if(!semiaccord(oneline, n, m, debug)){
    	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Jacobian matrix squared is P0 (general  kinetics)
    else if(!strcmp(filt, "JsquaredisP0")){
      if(reacJsquaredisP0(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJsquaredisP0")){
      if(!reacJsquaredisP0(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Second additive compound of Jacobian matrix is P0- (general  kinetics)
    else if(!strcmp(filt, "Jcomp2isP0minus")){
      if(reacJcomp2isP0(oneline, n, m, 1,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJcomp2isP0minus")){
      if(!reacJcomp2isP0(oneline, n, m, 1,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Second additive compound of Jacobian matrix is P0 (general  kinetics)
    else if(!strcmp(filt, "Jcomp2isP0")){
      if(reacJcomp2isP0(oneline, n, m, 0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJcomp2isP0")){
      if(!reacJcomp2isP0(oneline, n, m, 0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Second additive compound of Jacobian matrix is nonsingular (general  kinetics)
    else if(!strcmp(filt, "Jcomp2nonsing")){
      if(reacJcomp2nonsing(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJcomp2nonsing")){
      if(!reacJcomp2nonsing(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Second additive compound of Jacobian matrix has strictly signed determinant (general  kinetics) (In dimensions <=3 an identically zero determinant rules out Hopf bifurcation, but not more generally.)
    else if(!strcmp(filt, "Jcomp2detsigned")){
      if(reacJcomp2detsigned(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJcomp2detsigned")){
      if(!reacJcomp2detsigned(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //J^[2] has unsigned determinant (mixed vertices)
    //Note: this is *not* the negation of Jcomp2detsigned
    else if(!strcmp(filt, "Jcomp2detunsigned")){
      if(reacJcomp2detunsigned(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJcomp2detunsigned")){
      if(!reacJcomp2detunsigned(oneline, n, m, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //admits an imaginary pair (3D and 4D only, sufficient, but not necessary)
    else if(!strcmp(filt, "JadmitsIpair")){
      if(reacJadmitsIpair(oneline, n, m,iconst-1,&pppdegused,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJadmitsIpair")){
      if(!reacJadmitsIpair(oneline, n, m,iconst-1,&pppdegused,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Is det(J^2+I) for general kinetics strictly positive?
    else if(!strcmp(filt, "J2pIdetpos")){
      if(reacJ2pIdetpos(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJ2pIdetpos")){
      if(!reacJ2pIdetpos(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //
    // The rest are as above, but for the MA Jacobian matrix
    //
    //normal
    else if(!strcmp(filt, "normal")){
      if(normal(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notnormal")){
      if(!normal(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //not structurally discordant
    else if(!strcmp(filt, "PGKnormal") || !strcmp(filt, "notstructdiscord")){
      if(PGKnormal(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //structurally discordant
    else if(!strcmp(filt, "notPGKnormal") || !strcmp(filt, "structdiscord")){
      if(!PGKnormal(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }

    // Identically singular determinant at equilibria
    else if(!strcmp(filt, "JMAsingular")){
      if(reacJMAsingular(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //This means not identically singular (different from always nonsingular)
    else if(!strcmp(filt, "notJMAsingular")){
      if(!reacJMAsingular(oneline, n, m)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Jacobian matrix is *always* nonsingular (mass action, at equilibria)
    //rules out saddle node bifurcations
    else if(!strcmp(filt, "JMAnonsing")){
      if(reacJMAnonsing(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAnonsing")){
      if(!reacJMAnonsing(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Jacobian matrix (mass action, at +ve equilibria) has positive determinant
    else if(!strcmp(filt, "JMAdetpos")){
      if(reacJMAdetpos(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAdetpos")){
      if(!reacJMAdetpos(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Jacobian matrix at equilibria always has rank of stoichiometric matrix
    //better than nonsing, because we don't make assumptions about the rank
    else if(!strcmp(filt, "MAeqconcord")){
      if(reacMAeqconcord(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notMAeqconcord")){
      if(!reacMAeqconcord(oneline, n, m, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "QMAnegdef")){
      if(reacQMAposdef(oneline, n, m, 1,1, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notQMAnegdef")){
      if(!reacQMAposdef(oneline, n, m, 1,1, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "QMAnegsemidef")){
      if(reacQMAposdef(oneline, n, m, 1,0, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notQMAnegsemidef")){
      if(!reacQMAposdef(oneline, n, m, 1,0, iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "QMApossemidef")){
      if(reacQMAposdef(oneline, n, m, 0,0,iconst-1, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notQMApossemidef")){
      if(!reacQMAposdef(oneline, n, m, 0,0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "MAeqaccord") || !strcmp(filt, "JMAisP0minus")){
      if(reacJMAisP0(oneline, n, m, 1,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notMAeqaccord") || !strcmp(filt, "notJMAisP0minus")){
      if(!reacJMAisP0(oneline, n, m, 1,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAisP0")){
      if(reacJMAisP0(oneline, n, m, 0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAisP0")){
      if(!reacJMAisP0(oneline, n, m, 0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //JMA has mixed trace means the trace can definitely take all signs
    else if(!strcmp(filt, "mixedtrace")){
      if(reacJMAmixedtrace(oneline, n, m, 1,iconst-1,debug)==1){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //JMA notmixedtrace means the trace is definitely of constant sign
    else if(!strcmp(filt, "notmixedtrace")){
      if(reacJMAmixedtrace(oneline, n, m, 1,iconst-1,debug)==0){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAsquaredisP0")){
      if(reacJMAsquaredisP0(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAsquaredisP0")){
      if(!reacJMAsquaredisP0(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "DSR2MA1CondStar")){
      if(DSR2MA1CondStar(oneline, n, m,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notDSR2MA1CondStar")){
      if(DSR2MA1CondStar(oneline, n, m,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMA1comp2simppair")){
      if(reacJMA1comp2simppair(oneline, n, m,1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMA1comp2simppair")){
      if(!reacJMA1comp2simppair(oneline, n, m,1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAcomp2isP0minus")){
      if(reacJMAcomp2isP0(oneline, n, m,1,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAcomp2isP0minus")){
      if(!reacJMAcomp2isP0(oneline, n, m,1,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAcomp2isP0")){
      if(reacJMAcomp2isP0(oneline, n, m,0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAcomp2isP0")){
      if(!reacJMAcomp2isP0(oneline, n, m,0,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAcomp2nonsing")){
      if(reacJMAcomp2nonsing(oneline, n, m,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAcomp2nonsing")){
      if(!reacJMAcomp2nonsing(oneline, n, m,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAcomp2detsigned")){
      if(reacJMAcomp2detsigned(oneline, n, m,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAcomp2detsigned")){
      if(!reacJMAcomp2detsigned(oneline, n, m,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Note that this is not the negation of the last
    else if(!strcmp(filt, "JMAcomp2detunsigned")){
      if(reacJMAcomp2detunsigned(oneline, n, m,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAcomp2detunsigned")){
      if(!reacJMAcomp2detunsigned(oneline, n, m,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //admits an imaginary pair (3D only, sufficient, but not necessary)
    else if(!strcmp(filt, "JMAadmitsIpair")){
      if(reacJMAadmitsIpair(oneline, n, m,iconst-1,&pppdegused,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Note: this does *not* mean that JMA forbids an imaginary pair.
    else if(!strcmp(filt, "notJMAadmitsIpair")){
      if(!reacJMAadmitsIpair(oneline, n, m,iconst-1,&pppdegused,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    //Can mass action equilibria be second order degenerate (assuming they 
    //can be degenerate, checked notMAeqconcord)
    else if(!strcmp(filt, "MAeqdegen2")){
      if(reacJMAdegen2(oneline, n, m, iconst-1, &pppdegused, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notMAeqdegen2")){
      if(!reacJMAdegen2(oneline, n, m, iconst-1, &pppdegused, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "MABTbif")){
      if(reacJMABT(oneline, n, m, iconst-1, &pppdegused, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notMABTbif")){
      if(!reacJMABT(oneline, n, m, iconst-1, &pppdegused, debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "JMAcomp2detnonstationary")){
      if(reacJMAcomp2detnonstationary(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJMAcomp2detnonstationary")){
      if(!reacJMAcomp2detnonstationary(oneline, n, m, iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "J2pIMAdetpos")){
      if(reacJ2pIMAdetpos(oneline, n, m,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "notJ2pIMAdetpos")){
      if(!reacJ2pIMAdetpos(oneline, n, m,iconst-1,debug)){
	num++;fprintf(fd, "%s\n", oneline);
      }
    }
    else if(!strcmp(filt, "MAdegenerate")){
      if((val=reacJMAdegenerate(oneline, n, m,debug))){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "%s\n", val?"degenerate":"not degenerate");
    }
    else if(!strcmp(filt, "notMAdegenerate")){
      if(!(val=reacJMAdegenerate(oneline, n, m,debug))){
	num++;fprintf(fd, "%s\n", oneline);
      }
      if(debug)
	fprintf(stderr, "%s\n", val?"degenerate":"not degenerate");
    }
    else{
      fprintf(stderr, "\"%s\" is not a recognised filter. EXITING.\n",filt);
      exit(0);
    }
    if(debug){fprintf(stderr,"******\n");}
  }

  fclose(fd);fclose(fdin);
  return num;
}



int backtrack(int i, int **mat, int n, int *s, int ***C, int *totC, int *Sc, int *Scpos, int *marked, int *markedstack, int *markedstackpos){
  int f=0,j,k;
  Sc[(*Scpos)++]=i;//add to stack
  marked[i]=1;
  markedstack[(*markedstackpos)++]=i;//add to marked stack

  for (j=0;j<n;j++){
    if(mat[i][j]){// each out-edge
      if(j<*s){//do nothing or delete mat[i][j]?
      }
      else if(j==*s){//circuit
	//printvec(Sc, *Scpos);
	//store
	if((*totC)==0)
	  (*C) = (int**) malloc(sizeof(int*) * 1);
	else
	  (*C) =(int**) realloc((*C), sizeof(int*) *((*totC)+1));


	(*C)[(*totC)] = (int*) malloc(sizeof(int) * ((*Scpos)+1));

	(*C)[(*totC)][0]=(*Scpos);// first element is size
	for(k=0;k<(*Scpos);k++)
	  (*C)[(*totC)][k+1]=Sc[k];
	(*totC)++;

	f=1;
      }
      else if(marked[j]==0){
	if(backtrack(j, mat, n, s, C, totC, Sc, Scpos, marked, markedstack, markedstackpos)==1)
	  f=1;
      }
    }
  }

  if(f){//circuit found
    while((k=markedstack[(*markedstackpos)-1])!=i){
      marked[k]=0;
      (*markedstackpos)--;
    }
    marked[i]=0;
    (*markedstackpos)--;
  }
  (*Scpos)--;

  return f;

}


int **TarjanCircuits(int **mat, int n, int *totC){//adjacency matrix as input
  int s,k;
  int *marked=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *Sc= (int *)malloc((size_t) ((n)*sizeof(int)));//stack
  int Scpos=0;//stack-counter
  int *markedstack= (int *)malloc((size_t) ((n)*sizeof(int)));//marked stack
  int markedstackpos=0;//stack-counter
  int **C=NULL;

  inittozero(marked,n);

  for (s=0;s<n;s++){
    backtrack(s, mat, n, &s, &C, totC, Sc, &Scpos, marked, markedstack, &markedstackpos);
    for(k=markedstackpos-1;k>=0;k--)
      marked[markedstack[k]]=0;
    markedstackpos=0;
  }
  return C;

}

//Test program for TarjanCircuits
//https://stackoverflow.com/questions/25898100/ing-cycles-in-a-graph-using-tarjans-algorithm
void checkTarjanCircuits1(){
  int j, nlen,mlen;
  int **C;//to hold the circuits
  int totC=0;

  //int **mat=readmatrixfromstr((char *)"1 1 0 0 0 0 0 0 0 0\n0 0 0 0 1 0 1 1 0 0\n0 0 0 0 1 0 1 1 0 0\n0 0 0 0 1 0 1 1 0 0\n0 0 1 1 0 0 0 0 0 0\n0 0 1 1 0 0 0 0 0 0\n0 0 0 0 0 1 0 0 1 0\n0 0 0 0 0 1 0 0 1 0\n0 0 0 0 0 0 0 0 0 0\n0 0 0 0 0 0 0 0 0 0\n", &nlen, &mlen);
  int **mat=readmatrixfromstr((char *)"0  0  0 -1 -1  0  0 \n 0  0  0  0 -1 -1  0 \n 0  0  0 -1  0 -1  0 \n-1  0  1  0  0  0  0 \n-1 -1  0  0  0  0  0 \n 0  1 -1  0  0  0  0 \n 1  0  0  0  0  0  0 ", &nlen, &mlen);

  fprintf(stderr, "Adjacency matrix:\n");
  printmat(mat,nlen, mlen);
  C=TarjanCircuits(mat, nlen, &totC);
  fprintf(stderr, "total cycles: %d\n", totC);
  for(j=0;j<totC;j++)
    printvec(C[j]+1,C[j][0]); 
  free_imat(C,totC);
}



// Construct a DSR graph
//switch sign of the second factor
int **DSR(char *di6, int n, int m){
  int entries,i,j;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);

  //redefine AM to give AM of DSR
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      AM[i+n][j]-=AM[j][i+n];
      AM[j][i+n]=-AM[j][i+n];
    }
  }
  return AM;
}

//overloading: Modifies AM from Petri Net AM to DSR AM
void DSR(int **AM, int n, int m){
  int i,j;
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      AM[i+n][j]-=AM[j][i+n];
      AM[j][i+n]=-AM[j][i+n];
    }
  }
  return;
}

//version which creates a new matrix
int **DSRfromPN(int **AM, int n, int m){
  int i,j;
  int **DSRAM=imatrix(0, n+m-1, 0, n+m-1);
  inittozero(DSRAM,n+m,n+m);
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      DSRAM[i+n][j]=AM[i+n][j]-AM[j][i+n];
      DSRAM[j][i+n]=-AM[j][i+n];
    }
  }
  return DSRAM;
}

//overloading use S and V (or minus V) to construct DSR AM
//Here S is the stoichiometric matrix, and V is +- the pattern
//matrix for the reaction rate derivatives
int **DSR(int **S, int **V, int n, int m, int minus){
  int i,j;
  int **AM=imatrix(0, n+m-1, 0, n+m-1);
  for(i=0;i<n+m;i++){
    for(j=0;j<n+m;j++){
      if((i<n && j<n) || (i>=n && j>=n))
	AM[i][j]=0;
      else if (i<n && j>=n){
	if(minus)
	  AM[i][j]=-V[i][j-n];
	else
	  AM[i][j]=V[i][j-n];
      }
      else//i>=n && j<n)
	AM[i][j]=S[j][i-n];//St

    }
  }
  return AM;
}



void shiftforward(int *u, int len){
  int i,a=u[len-1];
  for(i=len-1;i>0;i--)
    u[i]=u[i-1];
  u[0]=a;
  return;
}

void shiftforward(int *u, int len, int n){
  int i;
  for(i=0;i<n;i++)
    shiftforward(u,len);
  return;
}

void shiftback(int *u, int len){
  int i,a=u[0];
  for(i=0;i<len-1;i++)
    u[i]=u[i+1];
  u[len-1]=a;
  return;
}

void shiftback(int *u, int len, int n){
  int i;
  for(i=0;i<n;i++)
    shiftback(u,len);
  return;
}

//return first index if edge found; otherwise -1
int checkedge(int v1, int v2, int *circ, int len){
  int i;
  if(v1==circ[len-1] && v2==circ[0])
    return len-1;
  for(i=0;i<len-1;i++){
    if(v1==circ[i] && v2==circ[i+1])
      return i;
  }
  return -1;
}

// align two cycles. Return 1 if some matching edges found
// and first match shifted to start; return 0 otherwise
int align(int *circ1, int len1, int *circ2, int len2,int debug){
  int t, q=0, flag=0;
  if(debug){fprintf(stderr, "\n###Entering align(...)\n");}

  while(q<len1 && checkedge(circ1[len1-1],circ1[0],circ2, len2)!=-1){
    shiftforward(circ1, len1);
    flag=1;
    q++;
  }
  if(q==len1){//never should occur: one cycle can't be a sub-cycle of the other
    fprintf(stderr, "ERROR in align: cycles are identical,q=%d. EXITING.\n",q);
    printvec(circ1,len1);
    printvec(circ2,len2);
    exit(0);
  }

  if(flag){//found edge and sent to start of circ 2; aligned
    t=checkedge(circ1[0],circ1[1],circ2, len2);
    shiftback(circ2, len2,t);
    if(debug){fprintf(stderr, "found edge and sent to start of circ 2; aligned\n");
      printvec(circ1,len1);
      printvec(circ2,len2);
    }
    return 1;
  }

  //search forwards
  q=0;
  while(q<len1-1 && (t=checkedge(circ1[0],circ1[1],circ2, len2))==-1){
    shiftback(circ1, len1);
    q++;
  }

  if(q==len1-1)//no matching edges
    return 0;

  //matching edge found - shift to start
  shiftback(circ2, len2,t);
  return 1;

}

// Check if two cycles have odd intersection (all components of the
// intersection have an odd number of edges)
// recursive algorithm. Call with oddsofar=0
int oddintersect(int *circ1, int len1, int *circ2, int len2, int *oddsofar, int debug){
  int i,q;
  int *circ1cp,*circ2cp, ret;
  if(debug){
    fprintf(stderr, "\n###Entering oddintersect. Checking paths\n");
    printvec(circ1,len1);
    printvec(circ2,len2);
  }
  if(len1<=1 || len2<=1)
    return (*oddsofar);

  //Beware: this alters the cycles
  if(!align(circ1,len1,circ2,len2,debug)){//no common edges, trivially satisfied 
    if(debug){fprintf(stderr, "No common edges, returning %d\n",*oddsofar);}
    return (*oddsofar);
  }

  if(debug){
    fprintf(stderr, "Common edges found. Aligned cycles:\n");
    printvec(circ1,len1);
    printvec(circ2,len2);
  }

  //cycles aligned and do have edge intersection. First component
  q=0;i=0;
  while (i<len1-1 && i<len2-1 && circ1[i]==circ2[i] && circ1[i+1]==circ2[i+1]){
    i++;q++;
  }

  if(q%2==0){//even component
    if(debug){fprintf(stderr, "intersection has even component: so not odd.\n");}
    return 0;
  }

  if(debug){fprintf(stderr, "intersection has odd component.\n");}
  (*oddsofar)=1;

  //So far odd
  circ1cp=veccp(circ1+q,len1-q);
  circ2cp=veccp(circ2+q,len2-q);

  ret=oddintersect(circ1cp, len1-q, circ2cp, len2-q, oddsofar,debug);
  free((char*)circ1cp);free((char*)circ2cp);
  return ret;
}

// Condition *
// report stores info on the outcome
// failure: -1=non-es cycles; -2=e-cycles with odd intersection
// success: 1=all odd
// Here DSRAM is the DSR adjacency matrix, *not* the PN adjacency matrix.
int DSRCondStar(int **DSRAM, int n, int m, int *report, int debug){
  int j,k,Csgn,p,q,Cpar;
  int lbltmp,Clbla, Clblb;//labels
  int totC=0;
  int **C=TarjanCircuits(DSRAM, n+m, &totC);
  int par[totC];
  int oddsofar;
  bool evenflag=0;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){
    fprintf(stderr, "There are %d cycles\n", totC);
    fprintf(stderr, "DSR adjacency matrix:\n");
    printmat(DSRAM,n+m,n+m);
    if(debugfull)
      fprintf(stderr, "\nCycles:\n");
  }
  for(j=0;j<totC;j++){
    //cycle sign
    Csgn=1;
    for(k=0;k<C[j][0];k++){//each edge
      p=C[j][k+1];
      if(k==C[j][0]-1){q=C[j][1];}else{q=C[j][k+2];}
      Csgn*=DSRAM[p][q];
    }
    //parity
    Cpar=(int)(pow((double)(-1),(double)(C[j][0]/2)))*Csgn;
    par[j]=Cpar;

    if(debug && debugfull){
      printvec(C[j]+1,C[j][0]); 
      fprintf(stderr, "%s\n", Cpar<0?"o-cycle":"e-cycle");
    }

    if(par[j]>0){//e-cycle: check if s-cycle
      evenflag=1;
      if(C[j][0]==2){//short e-cycle, hence non-es-cycle
	(*report)=-1;
	if(debug)
	  fprintf(stderr, "Found a non-es-cycle (short e-cycle). The network in this form fails Condition *.\n");
	free_imat(C,totC);
	return 0;
      }
      else{//not short
	Clbla=1;Clblb=1;
	for(k=0;k<C[j][0];k++){//each edge
	  p=C[j][k+1];
	  if(k==C[j][0]-1){q=C[j][1];}else{q=C[j][k+2];}
	  if(p>q){//edge from S (bottom-left block)
	    lbltmp=abs(DSRAM[p][q]);
	  }
	  else{//edge from V (top-right block)
	    if(DSRAM[q][p]*DSRAM[p][q]<=0){//No label to import
	      (*report)=-1;
	      if(debug)
		fprintf(stderr, "Found a non-es-cycle (absent label). The network in this form fails Condition *.\n");
	      free_imat(C,totC);
	      return 0;
	    }
	    lbltmp=abs(DSRAM[q][p]);
	  }
	  if(k%2==0)
	    Clbla*=lbltmp;
	  else
	    Clblb*=lbltmp;
	}
	if(!Clbla||!Clblb){//absent label: non-es-cycle
	  (*report)=-1;
	  if(debug)
	    fprintf(stderr, "Found a non-es-cycle (absent label). The network in this form fails Condition *.\n");
	  free_imat(C,totC);
	  return 0;
	}
	if(abs(Clbla-Clblb)!=0){//non-es-cycle
	  (*report)=-1;
	  if(debug){fprintf(stderr, "Found a non-es-cycle. The network in this form fails Condition *.\n");}
	  free_imat(C,totC);
	  return 0;
	}
      }
    }
  }

  if(evenflag==0){//no e-cycles
    (*report)=1;
    if(debug && debugfull){fprintf(stderr, "No e-cycles found\n");}
    free_imat(C,totC);
    return 0;
  }



  //So far so good. Check for odd intersections
  if(debug){fprintf(stderr, "No non-es-cycles found. Checking for odd intersections.\n");}
  for(j=0;j<totC-1;j++){
    if(par[j]>0){
      for(k=j+1;k<totC;k++){
	if(par[k]>0){//two even cycles
	  oddsofar=0;
	  if(debug && debugfull){
	    fprintf(stderr, "Checking:\n");
	    printvec(C[j]+1,C[j][0]);
	    printvec(C[k]+1,C[k][0]);
	  }
	  if(oddintersect(C[j]+1,C[j][0],C[k]+1,C[k][0],&oddsofar,debugfull)){//e-cycles with odd intersection
	    (*report)=-2;
	    if(debug){
	      fprintf(stderr, "Odd intersection between es-cycles. The network in this form fails Condition *.\n");
	      printvec(C[j]+1,C[j][0]); 
	      printvec(C[k]+1,C[k][0]); 
	    }
	    free_imat(C,totC);
	    return 0;
	  }
	  if(debug && debugfull){
	    fprintf(stderr, "Intersection not odd.\n\n");
	  }
	}
      }
    }
  }

  if(debug){fprintf(stderr, "Finished checking for odd intersections: none found. The network in this form satisfies Condition *\n\n");}

  free_imat(C,totC);
  (*report)=0;
  return 1;
}




// Overloading.
// rev=0 means don't reversify/merge cols, rev=1 means reversify & merge cols, rev=2 means 
// try without reversifying, and if that fails try with reversifying
// In case rev=2, return 2 if first fails and second succeeds
// Be careful with minus: for the DSR adjacency matrix we use S^t (bottom left) and -Dv (top right)
int DSRCondStar(int **S, int **V, int n, int m, int minus, int *report, int rev, int debug){
  int ret=0;
  int **Sout, **Vout, **S1, **Sout1, **Vout1, mout, mout1;
  int **DSRAM=DSR(S,V,n,m,minus);

  if(rev!=1)
    ret = DSRCondStar(DSRAM, n, m, report, debug);

  if(rev==1 || (rev==2 && !ret)){//reversify
    S1=cpmat(S,n,m);
    reduce_mat(S1, n, m, 1);//divide each column by its gcd
    reversify(S1, V, n, m, &Sout, &Vout, &mout);
    mergecols(Sout, Vout, n, mout, &Sout1, &Vout1, &mout1);
    //if(mout1!=mout || mout!=m){fprintf(stderr, "merged: %d, %d, %d\n", m, mout, mout1);exit(0);}
    ret=DSRCondStar(Sout1,Vout1, n, mout1,minus,report,0,debug);
    if(rev==2 && ret)//original failed and reversify succeeded
      ret++;
    free_imatrix(S1,0,n-1,0,m-1);
    free_imatrix(Sout,0,n-1,0,mout-1);
    free_imatrix(Vout,0,n-1,0,mout-1);
    free_imatrix(Sout1,0,n-1,0,mout1-1);
    free_imatrix(Vout1,0,n-1,0,mout1-1);
  }

  free_imatrix(DSRAM, 0, n+m-1, 0, n+m-1);
  if(ret==2 && debug){
    fprintf(stderr, "Reversifying made a difference!\n");
    //exit(0);
  }
  return ret;
}


//Overloading.
// rev=0 means don't reversify, rev=1 means reversify, rev=2 means 
// try without reversifying, and if that fails try with reversifying
// In case rev=2, return 2 if first fails and second succeeds
// Note: AM is the PN AM, not the DSR AM
int DSRCondStarPN(int **PNAM, int n, int m, int *report, int rev, int debug){
  int ret=0;
  int **Sout, **Vout, **Sout1, **Vout1, mout, mout1;
  int **V, **S;
  int **DSRAM=DSRfromPN(PNAM,n,m);

  if(debug){fprintf(stderr, "\n###Entering DSRCondStarPN. Checking condition *...\n");} 

  if(rev!=1)
    ret = DSRCondStar(DSRAM, n, m, report, debug);

  if(rev==1 || (rev==2 && !ret)){//reversify
    AMtoSSl(PNAM, n, m, 1, &S, &V);
    reduce_mat(S, n, m, 1);//divide each column by its gcd
    reversify(S, V, n, m, &Sout, &Vout, &mout);
    if(mout<m){
      if(debug){fprintf(stderr, "\nSimplifying the network and checking condition * again...\n");} 
      if(mout==m || !(ret=DSRCondStar(Sout,Vout, n, mout,0,report,0,debug))){
	mergecols(Sout, Vout, n, mout, &Sout1, &Vout1, &mout1);
	if(mout1<mout){
	  if(debug){fprintf(stderr, "\nSimplifying the network and checking condition * again...\n");} 
	  ret=DSRCondStar(Sout1,Vout1, n, mout1,0,report,0,debug);
	}
	free_imatrix(Sout1,0,n-1,0,mout1-1);
	free_imatrix(Vout1,0,n-1,0,mout1-1);
      }
    }

    if(rev==2 && ret)//original failed and reversify succeeded
      ret++;
    free_imatrix(V,0,n-1,0,m-1);
    free_imatrix(S,0,n-1,0,m-1);
    free_imatrix(Sout,0,n-1,0,mout-1);
    free_imatrix(Vout,0,n-1,0,mout-1);
  }

  free_imatrix(DSRAM,0,n+m-1,0,n+m-1);    
  if(ret==2 && debug){
    fprintf(stderr, "Network manipulation made a difference!\n\n");
    //exit(0);
  }
  return ret;
}

//Overloading. di6 input
// rev=0 means don't reversify, rev=1 means reversify, rev=2 means 
// try without reversifying, and if that fails try with reversifying
// In case rev=2, return 2 if first fails and second succeeds
int DSRCondStar(char *di6, int n, int m, int *report, int rev, int debug){
  int ret=0,entries;
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  ret=DSRCondStarPN(AM,n,m,report,rev,debug);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  return ret;
}


//Does Condition * hold for the DSR 2 graph?
int DSR2CondStar(int **S, int **V, int n, int m, int minus, int *report, int rev, int debug){
  int ret=0;
  int **S2f, **V2f, **S2, **V2, **Sout, **Vout, **S2a, **Sout1, **Vout1, mout, mout1;
  int cnk=n*(n-1)/2;
  int nm=n*m;
  int n0,m0;
  int **DSRAM;

  genS2from2(S, V, n, m, &S2f, &V2f);
  simppair(S2f,V2f,cnk,nm, &S2, &V2, &n0, &m0);
  if(debug){
    fprintf(stderr, "S, V:\n");printmat(S, n, m);printmat(V, n, m);
    fprintf(stderr, "S2f, V2f:\n");printmat(S2f, cnk, nm);printmat(V2f, cnk, nm);
    fprintf(stderr, "S2, V2:\n");printmat(S2, n0, m0);printmat(V2, n0, m0);
    //exit(0);
  }
  free_imatrix(S2f,0,cnk-1,0,mout-1);
  free_imatrix(V2f,0,cnk-1,0,mout-1);


  DSRAM=DSR(S2,V2,n0,m0,minus);

  if(rev!=1)
    ret = DSRCondStar(DSRAM, n0, m0, report, debug);

  if(rev==1 || (rev==2 && !ret)){//reversify and merge
    S2a=cpmat(S2,n0,m0);
    reduce_mat(S2, n0, m0, 1);//divide each column by its gcd
    reversify(S2a, V2, n0, m0, &Sout, &Vout, &mout);
    mergecols(Sout, Vout, n0, mout, &Sout1, &Vout1, &mout1);
    //if(mout1!=mout || mout!=m0){fprintf(stderr, "merged: %d, %d, %d\n", m0, mout, mout1);exit(0);}
    ret=DSRCondStar(Sout1,Vout1, n0, mout1,minus,report,0,debug);
    if(rev==2 && ret)//original failed and reversify succeeded
      ret++;
    free_imatrix(S2a,0,n0-1,0,m0-1);
    free_imatrix(Sout,0,n0-1,0,mout-1);
    free_imatrix(Vout,0,n0-1,0,mout-1);
    free_imatrix(Sout1,0,n0-1,0,mout1-1);
    free_imatrix(Vout1,0,n0-1,0,mout1-1);
  }

  free_imatrix(S2, 0, n0-1, 0, m0-1); 
  free_imatrix(V2, 0, n0-1, 0, m0-1); 

  free_imatrix(DSRAM, 0, n0+m0-1, 0, n0+m0-1);
  if(ret==2 && debug){
    fprintf(stderr, "Reversifying made a difference!\n");
    //exit(0);
  }
  return ret;
}

//Overloading: DSR (not DSR2) AM as input
int DSR2CondStar(int **DSRAM, int n, int m, int *report, int rev, int debug){
  int i,j,ret=0;
  int **S, **V;

  V=imatrix(0,n-1,0,m-1);//already assumed to be minus V by construction
  S=imatrix(0,n-1,0,m-1);
  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      S[j][i]=DSRAM[i+n][j];
      V[j][i]=DSRAM[j][i+n];
    }
  }
  //5th arg means V has correct sign already
  ret=DSR2CondStar(S, V, n, m, 0, report, rev, debug);
  free_imatrix(S, 0, n-1, 0, m-1); 
  free_imatrix(V, 0, n-1, 0, m-1); 
  return ret;
}


//Overloading. di6 input
int DSR2CondStar(char *di6, int n, int m, int *report, int rev, int debug){
  int ret=0;
  int **DSRAM=DSR(di6,n,m);
  ret=DSR2CondStar(DSRAM,n,m,report,rev,debug);
  free_imatrix(DSRAM, 0, n+m-1, 0, n+m-1);
  return ret;
}


void checkDSRCondStar(const char *fname, int debug){
  int n, m;
  int **Smat, **Vmat;
  int report=0;
  char *str=readfileintostr(fname);
  fprintf(stderr, "Checkng:\n%s\n", str);
  readmatpairfromstr(str,&n,&m,&Smat,&Vmat);
  if(debug){
    fprintf(stderr, "n=%d, m=%d\n", n,m);fflush(stderr);
    printmat(Smat,n,m);
    printmat(Vmat,n,m);
  }
  if(DSRCondStar(Smat,Vmat,n,m,0,&report,0,debug))
    fprintf(stderr, "Condition * holds%s\n", report==1?": no e-cycles found":"");
  else
    fprintf(stderr, "Condition * fails%s\n", report==-1?": es-cycles found":": e-cycles with odd intersection");
  free_imatrix(Smat, 0, n-1,0,m-1);
  free_imatrix(Vmat, 0, n-1,0,m-1);
  free(str);
  return;
}



// A basis of minimal vectors for ker[A|1]^t (or ker[A|1]).
// "minimal" means with as many zeros as possible, corresonding
// in the case ker[A|1]^t to dependencies amongst sources on faces
// of minimal dimension in the Newton Polytope
// "transpose" means ker[A|1]
int **minA1tkerbasis(int **Sil, int nt, int mt, int *tot, int *rk, int *deg, int transpose, int debug){
  int **basis=NULL;
  int i, j;
  int **Sil0=imatrix(0, nt, 0, mt-1);
  int **Sil1;
  int n,m;

  if(debug){fprintf(stderr, "\n###Entering minA1tkerbasis.\n");}
  
  (*tot)=0;

  //Augment with a row of ones
  for(i=0;i<nt;i++){
    for(j=0;j<mt;j++){
      Sil0[i][j]=Sil[i][j];
    }
  }
  for(j=0;j<mt;j++)
    Sil0[nt][j]=1;
  
  if(transpose){
    n=mt-1;m=nt+1;
    Sil1=transposemat(Sil0,nt+1,mt);
  }
  else{
    n=nt;m=mt;
    Sil1=cpmat(Sil0,nt+1,mt);
  }
  free_imatrix(Sil0,0, nt, 0, mt-1);
  //fprintf(stderr, "here: n = %d, m= %d\n",n,m);
  if(debug){printmat(Sil1,n+1,m);}
  
  basis=minkerbasis(Sil1,n+1,m,tot,rk,deg,debug);
  
  if(debug){fprintf(stderr, "###Exiting minA1tkerbasis.\n");}
  free_imatrix(Sil1,0, n, 0, m-1);
  return basis;
}


