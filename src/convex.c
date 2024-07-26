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
#include "convex.h"
#include "symbolic.h"
#include <limits.h>


//my implementation using GLPK
//Working to the left! each point is a *row* of imat1
//The desired vector *vec is a row-vector. 
//nlen is the number of rows (i.e., points)
//mlen is the length of each row (i.e., the dimension of the space)
//We can choose some or all rows (in the latter case, sets inds=NULL)
//x^t[1|imat1] = [1,vec] with x_i>=0

int isinconvhull(int **imat1, long *inds, int nlen, int mlen, int *vec){
  int j,l,m,lpstat;
  int *ia, *ja;
  double *ar;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  /* fprintf(stderr, "Entering isinconvhull...\n"); */
  /* fflush(stderr); */

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(nlen)*(mlen+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  //dimension of the space (number of variables in the polynomial)+boundedness constraint
  glp_add_rows(lp, 1+mlen);
  glp_set_row_bnds(lp, 1, GLP_FX, 1.0, 1.0);//sum x_i=1
  for(j=2;j<mlen+2;j++)
    glp_set_row_bnds(lp, j, GLP_FX, vec[j-2], vec[j-2]);


  //add the points
  glp_add_cols(lp, nlen);m=1;
  for(j=1;j<nlen+1;j++){//each point
    glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);// to ensure positivity
    glp_set_obj_coef(lp, j, 1.0); // objective function (not relevant as only feasibility)

    ia[m]=1;ja[m]=j;ar[m]=1.0;//boundedness constraint
    m++;
    for(l=2;l<2+mlen;l++){//the coordinates of the point
      ia[m]=l;ja[m]=j;
      if(inds)
	ar[m]=imat1[inds[j-1]][l-2];
      else
	ar[m]=imat1[j-1][l-2];
      m++;
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);
  //glp_write_lp(lp, NULL,"tmpfile.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);
  if(lpstat==GLP_NOFEAS) // no feasible solution
    goodflag=0;

  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  /* fprintf(stderr, "Exiting isinconvhull...\n"); */
  /* fflush(stderr); */

  return goodflag;
}

//Is a given (row) vector in the (positive) span of the rows
// of a given matrix. Simpler than isinconvhull.
//The desired vector is *vec
//nlen is the number of rows of the matrix
//mlen is the length of each row (the number of columns)
//We can choose some or all rows (in the latter case, sets inds=NULL)
//Solve x^t[imat1] = [vec] with x_i>=0

int isinspan(int **imat1, long *inds, int nlen, int mlen, int *vec, bool pos){
  int j,l,m,lpstat;
  int *ia, *ja;
  double *ar;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  /* fprintf(stderr, "Entering isincone...\n"); */
  /* fflush(stderr); */

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  //dimension of the space (number of variables)
  glp_add_rows(lp, mlen);
  for(j=1;j<mlen+1;j++)
    glp_set_row_bnds(lp, j, GLP_FX, vec[j-1], vec[j-1]);

  //add the points
  glp_add_cols(lp, nlen);m=1;
  for(j=1;j<nlen+1;j++){//each point
    if(pos)
      glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);// to ensure positivity
    else
      glp_set_col_bnds(lp, j, GLP_FR, 0.0, 0.0);// free
    glp_set_obj_coef(lp, j, 1.0); // objective function (not relevant as only feasibility)

    for(l=1;l<1+mlen;l++){//the coordinates of the point
      ia[m]=l;ja[m]=j;
      if(inds)
	ar[m]=imat1[inds[j-1]][l-1];
      else
	ar[m]=imat1[j-1][l-1];
      m++;
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);
  //glp_write_lp(lp, NULL,"tmpfile.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);
  if(lpstat==GLP_NOFEAS) // no feasible solution
    goodflag=0;

  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  /* fprintf(stderr, "Exiting isincone...\n"); */
  /* fflush(stderr); */

  return goodflag;
}

long *booltoinds(bool *vec, int len, int *n1){
  int i,i1=0;
  (*n1)=nonzentries(vec, len);
  long *out=(long*)malloc(sizeof(long)*(*n1));
  for(i=0;i<len;i++){
    if(vec[i])
      out[i1++]=i;
  }
  return out;
}

//overloading
int isinspan(int **imat1, bool *boolinds, int nlen, int mlen, int *vec, bool pos){
  int ret, n1;
  long *inds=booltoinds(boolinds, nlen, &n1);
  ret=isinspan(imat1, inds, n1, mlen, vec, pos);
  free((char*)inds);
  return ret;
}

//=is in the positive span
int isincone(int **imat1, long *inds, int nlen, int mlen, int *vec){
  return isinspan(imat1, inds, nlen, mlen, vec, 1);
}

//overloading
int isincone(int **imat1, bool *boolinds, int nlen, int mlen, int *vec){
  return isinspan(imat1, boolinds, nlen, mlen, vec, 1);
}


//The version where points are columns of the matrix
//"inds" is now a subset of the column indices, and 
//the target vector "vec" is of length "nlen"
int isinrspan(int **imat1, long *inds, int nlen, int mlen, int *vec, bool pos){
  int flg;
  int **tmp=transposemat(imat1,nlen,mlen);
  flg=isinspan(tmp, inds, mlen, nlen, vec, pos);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return flg;
}

int isinrspan(int **imat1, bool *boolinds, int nlen, int mlen, int *vec, bool pos){
  int flg;
  int **tmp=transposemat(imat1,nlen,mlen);
  flg=isinspan(tmp, boolinds, mlen, nlen, vec, pos);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return flg;
}


//=is in the positive span
int isinrcone(int **imat1, long *inds, int nlen, int mlen, int *vec){
  return isinrspan(imat1, inds, nlen, mlen, vec, 1);
}

//overloading
int isinrcone(int **imat1, bool *boolinds, int nlen, int mlen, int *vec){
  return isinrspan(imat1, boolinds, nlen, mlen, vec, 1);
}

//For debugging
void checkisincone(){
  int **imat=imatrix(0,2,0,2);
  int *vec=(int*) malloc(sizeof(int) * 2);
  int ret;

  imat[0][0]=-1;
  imat[0][1]=1;
  imat[1][0]=1;
  imat[1][1]=2;
  imat[2][0]=1;
  imat[2][1]=3;
  vec[0]=-10;
  vec[1]=11;

  printmat(imat,3,2);
  printvec(vec,2);

  ret=isincone(imat, (long*)NULL, 3, 2, vec);
  if(ret)
    fprintf(stderr, "yes, vector is in cone.\n");
  else
    fprintf(stderr, "no, vector is not in cone.\n");


  free_imatrix(imat,0,2,0,2);
  free((char*)vec);

}


//Is imat1[ind] in the convex hull of remaining (row) vectors?
//Beware repeated rows!
int isvert(int **imat1, int *inds, int nlen, int mlen, int ind){
  int i,tot=0,i1;
  long indsnew[nlen];

  for(i=0;i<nlen;i++){//Remove ind if necessary
    if(inds==NULL)
      i1=i;
    else
      i1=inds[i];
    if(i1!=ind)
      indsnew[tot++]=i1;
  }
  if(tot==0)//vacuously true
    return 1;

  if(isinconvhull(imat1, indsnew, tot, mlen, imat1[ind]))
    return 0;
  return 1;
}


//Is imat1[ind] in the positive span of remaining (row) vectors?
//Beware repeated rows!
int isextremeray(int **imat1, int *inds, int nlen, int mlen, int ind){
  int i,tot=0,i1;
  long indsnew[nlen];

  for(i=0;i<nlen;i++){//Remove ind if necessary
    if(inds==NULL)
      i1=i;
    else
      i1=inds[i];
    if(i1!=ind)
      indsnew[tot++]=i1;
  }
  if(tot==0)//vacuously true
    return 1;

  if(isincone(imat1, indsnew, tot, mlen, imat1[ind]))
    return 0;
  return 1;
}



//Assume no repeated poitns (rows): we check if each point is in the 
//convex hull of the remainder (which fails if there are repeated points)
//get the vertices and store in a boolean vector (length nlen, pre-allocated)
//points are *rows* of matrix. return is the number of vertices
int getvertbool(int **imat1, int nlen, int mlen, bool *verts){
  int i, numverts=0;
  for(i=0;i<nlen;i++){
    if(isvert(imat1, NULL, nlen, mlen, i)){
      verts[i]=1;
      numverts++;
    }
    else
      verts[i]=0;
  }

  return numverts;
}



//Assume no repeated points (columns)
//get the vertices and store in a boolean vector (length mlen, pre-allocated)
//points are *columns* of matrix. return is the number of vertices
int getvertboolT(int **imat1, int nlen, int mlen, bool *verts){
  int numverts=0;
  int **tmp=transposemat(imat1,nlen,mlen);
  numverts=getvertbool(tmp, mlen, nlen, verts);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return numverts;
}

//like getvertbool, except extreme rays of a cone. Vectors are rows of imat1
int getextremeraybool(int **imat1, int nlen, int mlen, bool *verts){
  int i, numverts=0;
  for(i=0;i<nlen;i++){
    if(isextremeray(imat1, NULL, nlen, mlen, i)){
      verts[i]=1;
      numverts++;
    }
    else
      verts[i]=0;
  }

  return numverts;
}


//like getvertboolT, except look for extreme rays of the cone
//spanned by the columns of imat1
int getextremerayboolT(int **imat1, int nlen, int mlen, bool *verts){
  int numverts=0;
  int **tmp=transposemat(imat1,nlen,mlen);
  numverts=getextremeraybool(tmp, mlen, nlen, verts);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return numverts;
}



bool has_nonz_signed_col(int **imat, int n, int m){
  int i,j,p;
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
    if(p)//not a mixed column
      return 1;
  }
  return 0;
}


// Does the transpose of imat1 have a positive vector in its kernel? 
// Checked only up to some tolerance
// "strict" means a strictly positive vector; otherwise nonnegative
int hasposkervec(int **imat1, int nlen, int mlen, bool strict){
  int j,l,m,lpstat;
  int *ia, *ja;
  double *ar,z;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;
  double tol=0.0;
  if(strict){
    if(has_nonz_signed_col(imat1, nlen, mlen))
      return 0;
    tol=0.01;
  }
  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(nlen)*(mlen+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);
  //number of constraints (rows) = col dimension of Gamma + boundedness constraint
  glp_add_rows(lp, 1+mlen);
  glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 10.0);
  for(j=2;j<mlen+2;j++)
    glp_set_row_bnds(lp, j, GLP_FX, 0.0, 0.0);

  //number of variables (cols) = row dimension of Gamma
  glp_add_cols(lp, nlen);m=1;
  for(j=1;j<nlen+1;j++){//col indices
    glp_set_col_bnds(lp, j, GLP_LO, tol, 0.0);// to ensure positivity: could introduce false negatives. 
    glp_set_obj_coef(lp, j, 1.0); // objective function

    ia[m]=1;ja[m]=j;ar[m]=1.0;//boundedness constraint
    m++;
    for(l=2;l<2+mlen;l++){//row indices
      ia[m]=l;ja[m]=j;ar[m]=imat1[j-1][l-2];
      m++;
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);

  //glp_write_lp(lp, NULL,"tmpfile.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);
  if(strict){// probably no positive vector
    if(lpstat!=GLP_OPT && lpstat!=GLP_FEAS)// infeasible
      goodflag=-1;//probably (unless a very marginally positive vector)!
  }
  else{// no nonnegative vector (must have z value 10 if it exists)
    if(lpstat==GLP_OPT && (z=glp_get_obj_val(lp))<0.01) // solution is zero
      goodflag=0;
  }
  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  return goodflag;
}

// LCM of expressions in a list. Uses GiNAC lcm
ex LCmult(ex *exlist, int n){
  int i;
  ex L=exlist[0];
  for(i=1;i<n;i++)
    L=lcm(L,exlist[i]);
  return L;
}

// Find a positive vector in the kernel of nlen X mlen matrix imat
// Return in the vector intvec assumed to be of length mlen
// Using GiNAC solver, and LCM
// Assumes that it has already been checked that:
// 1) ker(imat) includes a positive vector (e.g. via hasposkervec)
// 2) ker(imat) is 1D, e.g. via checking the rank of imat
void posintkervec(int **imat, int nlen, int mlen, int intvec[], bool q){
  int i;
  matrix J = imattoexmat(imat, nlen, mlen);
  matrix kern(nlen,1);
  matrix vars(mlen,1);
  ex x[mlen], y[mlen], z[mlen];
  char sstr[5];
  matrix sols;
  ex LCM;

  for(i=0;i<nlen;i++)//kernel vector
    kern(i,0)=0;

  if(!q){fprintf(stderr, "Checking matrix:\n");printmat(J,nlen,mlen);}

  for(i=0;i<mlen;i++){//vector of unknowns
    sprintf(sstr,"x%d",i);
    x[i]=get_possymbol(sstr);
    vars(i,0)=x[i];
  }

  sols=J.solve(vars,kern);
  //cerr<<sols<<endl;
  for(i=0;i<mlen;i++){
    x[i]=sols(i,0).subs(x[mlen-1]==1);
    y[i]=x[i].denom();
  }
  LCM=LCmult(y, mlen);
  //cerr<<LCM<<endl;
  for(i=0;i<mlen;i++){
    z[i]=LCM*x[i];
    if(!(ex_to<numeric>(z[i]).is_integer())){
      cerr<< "ERROR in posintkervec. Noninteger quantity: " << z[i] << endl;
      exit(0);
    }
    else
      intvec[i]=ex_to<numeric>(z[i]).to_int();
  }

  if(!q){fprintf(stderr, "positive vector in kernel: ");printvec(intvec,mlen);}
  return;
}


// Find a vector in the kernel of nlen X mlen matrix imat
// Return in the vector intvec assumed to be of length mlen
// Using GiNAC solver, and LCM
// Assumes that it has already been checked that:
// 2) ker(imat) is 1D, e.g. via checking the rank of imat
void intkervec(int **imat, int nlen, int mlen, int intvec[], bool q){
  int i;
  matrix J = imattoexmat(imat, nlen, mlen);
  matrix kern(nlen,1);
  matrix vars(mlen,1);
  ex x[mlen], y[mlen], z[mlen];
  char sstr[5];
  matrix sols;
  ex LCM;
  int lastnonz=mlen-1;

  for(i=0;i<nlen;i++)//kernel vector
    kern(i,0)=0;

  if(!q){fprintf(stderr, "Checking matrix:\n");printmat(J,nlen,mlen);}

  for(i=0;i<mlen;i++){//vector of unknowns
    sprintf(sstr,"x%d",i);
    x[i]=get_possymbol(sstr);
    vars(i,0)=x[i];
  }


  sols=J.solve(vars,kern);
  //cout<<sols<<endl;
  for(i=mlen-1;i>=0;i--){
    if(sols[i]!=0){
      lastnonz=i;
      //fprintf(stderr, "lastnonz=%d\n",lastnonz);
      break;
    }
  }
  //cerr<<sols<<endl;
  for(i=0;i<mlen;i++){
    x[i]=sols(i,0).subs(x[lastnonz]==1);
    y[i]=x[i].denom();
  }
  LCM=LCmult(y, mlen);
  //cerr<<LCM<<endl;
  for(i=0;i<mlen;i++){
    z[i]=LCM*x[i];
    if(!(ex_to<numeric>(z[i]).is_integer())){
      cerr<< "ERROR in intkervec. Noninteger quantity: " << z[i] << endl;
      exit(0);
    }
    else
      intvec[i]=ex_to<numeric>(z[i]).to_int();
  }

  if(!q){fprintf(stderr, "vector in kernel: ");printvec(intvec,mlen);}
  return;
}

// A basis of minimal vectors for ker M.
// "minimal" means with as many zeros as possible, corresonding
int **minkerbasis(int **M, int n, int m, int *tot, int *rk, int *deg, int debug){
  int kerdim;//dimension of kernel
  int **basis=NULL;
  int **basistmp;
  int tottmp=0;
  int **faces=NULL;
  int **redSi;
  int flag;
  int i, k0;
  int *alldegs;
  int xc[n],yc[m],tmpvec[m],tmpvec1[m],tmpvec2[m];

  if(debug){fprintf(stderr, "\n###Entering minkerbasis:\n");printmat(M,n,m);}
  
  (*tot)=0;

  (*rk)=matrank(M,n,m);//rank of augmented matrix=dim(affhull(.))+1
  kerdim=m-(*rk);
  

  //no vectors in kernel
  //  if(kerdim<=0 || hasposrkervec(Si,n,m,1)!=1)
  if(kerdim<=0){// deal with this possibility prior to using this routine
    (*deg)=0;
    fprintf(stderr, "ERROR in \"minkerbasis\": the kernel must not be trivial. EXITING.\n");
    return NULL;
  }

  //basis=imatrix(0, kerdim-1, 0, m-1);
  if(kerdim==1){
    intkervec(M, n, m, tmpvec1,1);
    (*tot)=addnewtoarray(&basis, *tot, tmpvec1, m);
    (*deg)=sumpos(tmpvec1,m);
    if(debug){fprintf(stderr, "deg = %d\n", (*deg));}
    if(debug){fprintf(stderr, "###Exiting minkerbasis.\n");}
    return basis;
  }


  firstcomb(xc, n, n);//all rows
  k0=1;
  while(k0<m){
    firstcomb(yc, m, k0);flag=1;
    while(flag==1){//each vector with support of size k0
      //this should avoid non-extreme vectors being found; intkervec fails otherwise
      if(!supervec(yc, k0, faces, tottmp)){
	//fprintf(stderr, "%d\n", *tot);
	redSi=submat(M,n,m,xc,n,yc,k0);
	if(debug){fprintf(stderr, "**Checking submatrix:\n");printmat(redSi,n,k0);}
	if(matrank(redSi, n, k0)<k0){// has nontrivial kernel
	  intkervec(redSi, n, k0, tmpvec,1);
	  merge(yc,tmpvec,k0,tmpvec1,m);
	  addnewtoarray(&basistmp, tottmp, tmpvec1, m);
	  //addnewtoarray(&basis, *tot, tmpvec1, m);
	  addnewto1Darray(&alldegs, tottmp, sumpos(tmpvec1,m));

	  if(debug){printvec(basistmp[tottmp],m);printvec(yc,k0);}
	  tmpvec2[0]=k0;
	  for(i=0;i<k0;i++)
	    tmpvec2[i+1]=yc[i];
	  addnewtoarray(&faces, tottmp, tmpvec2, k0+1);

	  //(*tot)++;
	  tottmp++;
	}
	free_imatrix(redSi,0,n-1,0,k0-1);
      }
      flag=nextcomb(yc, m, k0);
    }
    //if((*tot)<kerdim)
      k0++;
  }
  if(tottmp<kerdim){
    fprintf(stderr, "ERROR in \"minkerbasis\": could not find enough basis vectors. This shouldn't happen. EXITING.\n");
    exit(0);
  }

  //Extract the best "kerdim" in terms of Bezout degree
  int inds[tottmp];
  for(i=0;i<tottmp;i++)
    inds[i]=i;
  qsort2(alldegs,inds,0,tottmp-1);
  for(i=0;i<kerdim;i++){
    addnewtoarray(&basis, *tot, basistmp[inds[i]], m);
    (*tot)++;
  }
  (*deg)=1;
  for(i=0;i<kerdim;i++)
    (*deg)*=alldegs[i];

  
  free_imat(basistmp, tottmp);
  free_imat(faces, tottmp);
  if(debug){fprintf(stderr, "###Exiting minkerbasis.\n");}
  free((char*)alldegs);
  return basis;
}


//Extreme vectors of the positive kernel of Si as rows of output
int **poskerbasis(int **Si, int n, int m, int *tot, int *rk){
  (*rk)=matrank(Si,n,m);//rank
  int kerdim=m-(*rk);//dimension of kernel
  int xc[n],yc[m],tmpvec[m],tmpvec1[m],tmpvec2[m];
  int **basis=NULL;
  int **faces=NULL;
  int **redSi;
  int flag;
  int i, k0;
  int debug=0;
  (*tot)=0;
  if(debug){fprintf(stderr, "\n###Entering poskerbasis. Rank=%d\nSi:\n", (*rk));printmat(Si,n,m);}

  //no +ve vectors in kernel
  //  if(kerdim<=0 || hasposrkervec(Si,n,m,1)!=1)
  if(kerdim<=0 || hasposlimvec(Si,n,m))
    return NULL;

  //basis=imatrix(0, kerdim-1, 0, m-1);
  if(kerdim==1){
    posintkervec(Si, n, m, tmpvec1,1);
    (*tot)=addnewtoarray(&basis, *tot, tmpvec1, m);
    return basis;
  }


  firstcomb(xc, n, n);//all rows
  k0=2;//initial dimension set to two as we don't allow trivial reactions
//  while(k0<=m-kerdim+1 /* && (*tot)<kerdim */){
  while(k0<m){
    firstcomb(yc, m, k0);flag=1;
    while(flag==1){//each vector with support of size k0
      //this should avoid non-extreme vectors being found; posintkervec fails otherwise
      if(!supervec(yc, k0, faces, (*tot))){
	//fprintf(stderr, "%d\n", *tot);
	redSi=submat(Si,n,m,xc,n,yc,k0);
	if(debug){fprintf(stderr, "**Checking submatrix:\n");printmat(redSi,n,k0);}
	if(!hasposlimvec(redSi,n,k0)){
	  posintkervec(redSi, n, k0, tmpvec,1);
	  merge(yc,tmpvec,k0,tmpvec1,m);
	  addnewtoarray(&basis, *tot, tmpvec1, m);
	  if(debug){printvec(basis[(*tot)],m);printvec(yc,k0);}
	  tmpvec2[0]=k0;
	  for(i=0;i<k0;i++)
	    tmpvec2[i+1]=yc[i];
	  addnewtoarray(&faces, *tot, tmpvec2, k0+1);

	  (*tot)++;
	  if(kerdim<=2 && (*tot)==kerdim){//done
	    free_imatrix(redSi,0,n-1,0,k0-1);
	    free_imat(faces, *tot);
	    if(debug){fprintf(stderr, "Exiting poskerbasis.\n");}
	    return basis; 
	  }
	}
	free_imatrix(redSi,0,n-1,0,k0-1);
      }
      flag=nextcomb(yc, m, k0);
    }
    //if((*tot)<kerdim)
      k0++;
  }


  if((*tot)>kerdim){//only do this if we have a possibly non-simplicial cone
    bool verts[*tot];
    int numverts=getextremeraybool(basis, *tot, m, verts);
    if(numverts!=*tot){//not all vectors found are extreme vectors: can this happen?
      fprintf(stderr, "\tUnexpected behaviour in \"poskerbasis\".\n\tSome vectors found which were not extreme.\n\tExplore why. EXITING.\n");exit(0);
      int tot1;
      int **basisnew=submatfromrowbool(basis, *tot, m, verts, &tot1);
      free_imat(faces, *tot);
      free_imat(basis, *tot);
      (*tot)=tot1;
      return basisnew;
    }
    else if(debug)
      fprintf(stderr, "The positive part of the kernel forms a non-simplicial cone.\n");
  }
  free_imat(faces, *tot);
  if(debug){fprintf(stderr, "Exiting poskerbasis.\n");}
  return basis;
}




// Does imat1 have a nonnegative vector in its image?
int hasposimvec(int **imat1, int nlen, int mlen){
  int j,l,m,lpstat,csum;
  int *ia, *ja;
  double *ar;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(mlen+1)*(nlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(mlen+1)*(nlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(mlen)*(nlen+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);
  //number of constraints (rows) = row dimension of Gamma + boundedness constraint
  glp_add_rows(lp, 1+nlen);
  glp_set_row_bnds(lp, 1, GLP_DB, 1.0, 10.0);
  for(j=2;j<nlen+2;j++)
    glp_set_row_bnds(lp, j, GLP_LO, 0.0, 0.0);

  //number of variables (cols) = col dimension of Gamma
  glp_add_cols(lp, mlen);
  m=1;
  for(j=1;j<mlen+1;j++){
    glp_set_col_bnds(lp, j, GLP_FR, 0.0, 0.0);// free
    csum=colsum(imat1, nlen, mlen, j-1);
    glp_set_obj_coef(lp, j, csum); // objective function
    ia[m]=1;ja[m]=j;ar[m]=csum;//boundedness constraint
    m++;
    for(l=2;l<2+nlen;l++){//row indices
      ia[m]=l;ja[m]=j;ar[m]=imat1[l-2][j-1];
      m++;
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);

  //glp_write_lp(lp, NULL,"tmpfile1.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);

  if(lpstat==GLP_NOFEAS) // no feasible solution
    goodflag=0;

  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  return goodflag;
}

// Overloading: the submatrix version: rows indexed by inds
// Does imat1 have a nonnegative vector in its image?
int hasposimvec(int **imat1, int nlen, int mlen, bool *inds, int numinds){
  int j,l,m,lpstat,csum;
  int *ia, *ja;
  double *ar;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(mlen+1)*(numinds+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(mlen+1)*(numinds+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(mlen)*(numinds+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);
  //number of constraints (rows) = row dimension of Gamma + boundedness constraint
  glp_add_rows(lp, 1+numinds);
  glp_set_row_bnds(lp, 1, GLP_DB, 1.0, 10.0);
  for(j=2;j<nlen+2;j++){
    if(inds[j-2])
      glp_set_row_bnds(lp, j, GLP_LO, 0.0, 0.0);
  }

  //number of variables (cols) = col dimension of Gamma
  glp_add_cols(lp, mlen);
  m=1;
  for(j=1;j<mlen+1;j++){
    glp_set_col_bnds(lp, j, GLP_FR, 0.0, 0.0);// free
    csum=colsum(imat1, nlen, mlen, j-1);
    glp_set_obj_coef(lp, j, csum); // objective function
    ia[m]=1;ja[m]=j;ar[m]=csum;//boundedness constraint
    m++;
    for(l=2;l<2+nlen;l++){//row indices
      if(inds[j-2]){
	ia[m]=l;ja[m]=j;ar[m]=imat1[l-2][j-1];
	m++;
      }
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);

  //glp_write_lp(lp, NULL,"tmpfile1.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);

  if(lpstat==GLP_NOFEAS) // no feasible solution
    goodflag=0;

  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  return goodflag;
}


//Does the kernel of imat1 include a nonnegative (or strictly positive)
//vector
int hasposrkervec(int **imat1, int nlen, int mlen, bool strict){
  int flg;
  int **tmp=transposemat(imat1,nlen,mlen);
  flg=hasposkervec(tmp,mlen,nlen,strict);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return flg;
}

//nonnegative vector in the image of the transpose
//equivalent to no positive vector in the kernel
int hasposlimvec(int **imat1, int nlen, int mlen){
  int flg;
  int **tmp=transposemat(imat1,nlen,mlen);
  flg=hasposimvec(tmp,mlen,nlen);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return flg;
}

//Do the points indexed by "face" lie in a proper
//face of the (closed, convex) cone spanned by allpoints?
//assume we have already computed the extreme vectors of the cone
//and stored these in "verts"
//We check if a *positive* combination of points of face can equal
//a *positive* combination of points of verts
//This can be formulated as a questions of whether the kernel of
//an augmented matrix includes a positive (>> 0) vector; equivalently, whether
//the row space includes no vectors > 0.
//Details in CRNT_notes.pdf
int ispropsubface(int **allpoints, int n, int m, bool *verts, int numverts, bool *face, int numface){
  int i,j,k=0, ret;
  int **imatT=imatrix(0,numverts+numface-1,0,n-1);//construct transposed to save a transposition

  //vertices
  for(j=0;j<m;j++){//each column of allpoints
    if(verts[j]){
      for(i=0;i<n;i++)
	imatT[k][i]=allpoints[i][j];
      k++;
    }
  }

  //face
  for(j=0;j<m;j++){
    if(face[j]){
      for(i=0;i<n;i++)
	imatT[k][i]=-allpoints[i][j];
      k++;
    }
  }

  ret=hasposimvec(imatT, numverts+numface,n);
  free_imatrix(imatT,0,n-1,0,numverts+numface-1);

  return ret;
}


//"in" is a boolean vector of length "len". "comb" holds indices of a subset
//of nonzero entries of "comb". we convert "comb" to a boolean vector "out"
//of length "len" - a subvector of "in". No error checking.
void subboolvec(bool *in, int len, int *comb, int comblen, bool *out){
  int i, t, t1;
  inittozero(out, len);
  t=0;t1=0;
  for(i=0;i<len;i++){
    if(in[i]){
      if(comb[t1]==t){
	out[i]=1;
	t1++;
      }
      t++;
    }
  }
  return;
}

//Finds maximal set of vertices of same rank
//returns size of the set
int growsubface(int **allverts, int n, int numverts, bool *v, int sz, int rk){
  int i,tot=sz;
  bool fullvec[n];
  bool v1[numverts];
  inittoone(fullvec,n);
  veccp(v,numverts,v1);
  for(i=0;i<numverts;i++){
    if(!v[i]){
      v1[i]=1;
      if(submatrank(allverts, n, numverts, fullvec, v1)==rk){//dependent vertex
	v[i]=1;
	tot++;
      }
      else
	v1[i]=0;
    }
  }
  return tot;
}


//subset is a subset of the vertices
//verts is of length m with numverts nonzero entries
//it hold the vertices of the whole polytope
//"faceverts" is of length numverts with numfaceverts nonzero entries
//it holds the vertices of a face of rank "facerank"
//"current" is of length numverts with numcurrent nonzero entries
//it holds the vertices of the current subset of a face (also of rank "facerank") 
//"ranks" and "facets" are global
//initialise "ranks" to -1, initialise "facets" to 0;
//Note that this routine will fail if there are more than about 63 vertices
//as "index" cannot be set.

void visitfacets(int **allverts, int n, int numverts, bool *faceverts, int numfaceverts, int facerank, bool *current, int numcurrent, int *ranks, bool *facets, bool ***allfacets, int *totfacets, int debug){
  int k=numcurrent-1;
  int k1, pos;
  long r1, cnk;
  int **vcombs;
  bool v[numverts];
  unsigned long index;
  bool fullvec[n];

  if(facerank<=2)//We already know the vertices
    return;

  if(debug){fprintf(stderr, "\n###Entering visitfacets. Current set: ");printvec(current, numverts);}

  inittoone(fullvec, n);

  //generate combinations on numfaceverts-1 only - recursion can take deeper
  vcombs=allcombsgen(numcurrent,k);
  cnk=comb(numcurrent,k);
  for(r1=0;r1<cnk;r1++){//each subset of face of size numcurrent-1

    //convert vcombs[r1] to a vert bool
    subboolvec(current, numverts, vcombs[r1], k, v);
    index=lbinltoint(v, numverts);//get the index
    if(debug){
      fprintf(stderr, "\nExamining %ld: ", index);printvec(v, numverts);
      fprintf(stderr, "Subset of face: ");printvec(faceverts, numverts);
    }

    //already visited, then skip
    if(ranks[index]!=-1){//continue
      if(debug){fprintf(stderr, "Already visited this set; ignoring.\n");}
      continue;
    }

    ranks[index]=submatrank(allverts, n, numverts, fullvec, v);
    if(debug){fprintf(stderr, "rank=%d\n", ranks[index]);}

    //rank too low - shouldn't happen, but here for debugging
    if(ranks[index]<facerank-1){//should never occur: only going down one level at a time
      fprintf(stderr, "Somehow we reached a subset of rank which is too low. EXITING.\n");
      exit(0);
    }

    //recursion of two different kinds
    if(ranks[index]==facerank){
      visitfacets(allverts, n, numverts, faceverts, numfaceverts, facerank, v, k, ranks, facets, allfacets, totfacets, debug);
    }
    else if(ispropsubface(allverts, n, numverts, faceverts, numfaceverts, v, k)){
      if(debug){fprintf(stderr, "found a subface: ");printvec(v,numverts);}//may not be maximal
      k1=growsubface(allverts, n, numverts, v, k, facerank-1);

      index=lbinltoint(v, numverts);//get the index
      if(facets[index])//already dealt with by some other path
	continue;

      facets[index]=1;
      ranks[index]=facerank-1;
      if(debug){fprintf(stderr, "**found a new rank %d face: ", ranks[index]);printvec(v,numverts);}//may not be maximal
      (*totfacets)=addvec1((*totfacets), v, numverts, allfacets, 0, &pos);//penultimate means don't check if already there

      visitfacets(allverts, n, numverts, v, k1, ranks[index], v, k1, ranks, facets, allfacets, totfacets, debug);
    }

    //finished face and all subfaces: perhaps add to a "completed" list?

    //reaching here means lower rank, but not a subface, which means setting facets[index]=0 (default)
    // ignore
  }
  free_imatrix(vcombs, 0, cnk-1, 0, k-1);

  if(debug){fprintf(stderr, "Exiting visitfacets.\n");}
  return;
}

//The version where we don't store ranks, etc. 
void visitfacets2(int **allverts, int n, int numverts, bool *faceverts, int numfaceverts, int facerank, bool *current, int numcurrent, bool ***allfacets, int *totfacets, int debug){
  int k=numcurrent-1;
  int k1, pos;
  long r1, cnk;
  int **vcombs;
  bool v[numverts];
  //unsigned long index;
  bool fullvec[n];
  int currank;
  int newtotfacets;

  if(facerank<=2)//We already know the vertices
    return;

  if(debug){fprintf(stderr, "\n###Entering visitfacets2. Current set: ");printvec(current, numverts);}

  inittoone(fullvec, n);

  //generate combinations on numfaceverts-1 only - recursion can take deeper
  vcombs=allcombsgen(numcurrent,k);
  cnk=comb(numcurrent,k);
  for(r1=0;r1<cnk;r1++){//each subset of face of size numcurrent-1

    //convert vcombs[r1] to a vert bool
    subboolvec(current, numverts, vcombs[r1], k, v);
    //index=lbinltoint(v, numverts);//get the index
    if(debug){
      fprintf(stderr, "\nExamining: ");printvec(v, numverts);
      fprintf(stderr, "Subset of face: ");printvec(faceverts, numverts);
    }


    currank=submatrank(allverts, n, numverts, fullvec, v);
    if(debug){fprintf(stderr, "rank=%d\n", currank);}

    //rank too low - shouldn't happen, but here for debugging
    if(currank<facerank-1){//should never occur: only going down one level at a time
      fprintf(stderr, "Somehow we reached a subset of rank which is too low. EXITING.\n");
      exit(0);
    }

    //recursion of two different kinds
    if(currank==facerank){//not a new face
      visitfacets2(allverts, n, numverts, faceverts, numfaceverts, facerank, v, k, allfacets, totfacets, debug);
    }
    else if(ispropsubface(allverts, n, numverts, faceverts, numfaceverts, v, k)){
      if(debug){fprintf(stderr, "found a subface: ");printvec(v,numverts);}//may not be maximal
      k1=growsubface(allverts, n, numverts, v, k, facerank-1);

      
      newtotfacets=addvec1((*totfacets), v, numverts, allfacets, 1, &pos);//penultimate means check if already there
	
      if(newtotfacets>(*totfacets)){
	(*totfacets)=newtotfacets;
	if(debug){fprintf(stderr, "**found a new rank %d face: ", currank);printvec(v,numverts);}
	visitfacets2(allverts, n, numverts, v, k1, currank, v, k1, allfacets, totfacets, debug);
      }
    }
  }
  free_imatrix(vcombs, 0, cnk-1, 0, k-1);

  if(debug){fprintf(stderr, "Exiting visitfacets.\n");}
  return;
}




//We've identified a set of points (columns), stored in a
//bool vector "verts"; now extract the submatrix augmented with a
//row of ones (like "submatfromcolbool" in basics.c)
int **conemat(int **A, int n, int m, bool *verts, int *m1){
  int i,j, j1=0;
  (*m1)=nonzentries(verts, m);
  int **C=imatrix(0, n, 0, (*m1)-1);//extra row

  for(j=0;j<m;j++){//each col
    if(verts[j]){
      for(i=0;i<n;i++)//each row
	C[i][j1]=A[i][j];
      C[n][j1]=1;//extra row
      j1++;
    }
  }
  return C;
}

//dimension of affine hull of columns of A
int affrank(int **A, int n, int m){
  int i,ret;
  int **C=imatrix(0, n, 0, m);
  cpmat(A,C,n,m);
  for(i=0;i<m;i++)//row of ones
    C[n][i]=1;
  ret=matrank(C,n+1,m);
  free_imatrix(C,0,n,0,m-1);
  return ret-1;  
}

//Get the facial structure of the convex hull of a set of points which are
//the columns of the n X m matrix Sl. Assume no repeated columns: these should
//be removed in a preprocessing stage
int getfacestruct(int **Sl, int n, int m, bool *verts, bool ***allfacets, int debug){
  int m1;
  unsigned long k;
  int numverts=getvertboolT(Sl, n, m, verts);

  bool faceverts[numverts];
  inittoone(faceverts,numverts);//full
  int **vertmat=conemat(Sl,n,m,verts,&m1);
  int rk=matrank(vertmat,n+1,numverts);
  unsigned long numsubsets=(unsigned long)pow(2.0, (double)numverts);
  bool facets[numsubsets];
  int ranks[numsubsets];

  int numfacets=0;//the number of facets of rank 2 to rk-1 (i.e., not the vertices, and not the full polytope)


  if(debug){
    fprintf(stderr, "\n###Entering getfacestruct. Matrix:\n");
    printmat(Sl, n,m);
    fprintf(stderr, "Cone matrix:\n");
    printmat(vertmat, n+1, numverts);
  }

  //initialise
  for(k=0;k<numsubsets;k++){ranks[k]=-1;facets[k]=0;}

  visitfacets(vertmat,n+1, numverts, faceverts, numverts, rk, faceverts, numverts, ranks, facets, allfacets, &numfacets, debug);

  free_imatrix(vertmat,0,n,0,m1-1);
  return numfacets;
}


int getfacestruct2(int **Sl, int n, int m, bool *verts, bool ***allfacets, int debug){
  int m1;
  //unsigned long k;
  int numverts=getvertboolT(Sl, n, m, verts);

  bool faceverts[numverts];
  inittoone(faceverts,numverts);//full
  int **vertmat=conemat(Sl,n,m,verts,&m1);
  int rk=matrank(vertmat,n+1,numverts);
  //unsigned long numsubsets=(unsigned long)pow(2.0, (double)numverts);
  //bool facets[numsubsets];
  //int ranks[numsubsets];

  int numfacets=0;//the number of facets of rank 2 to rk-1 (i.e., not the vertices, and not the full polytope)


  if(debug){
    fprintf(stderr, "\n###Entering getfacestruct2. Matrix:\n");
    printmat(Sl, n,m);
    fprintf(stderr, "Cone matrix:\n");
    printmat(vertmat, n+1, numverts);
  }

  //initialise
  //for(k=0;k<numsubsets;k++){ranks[k]=-1;facets[k]=0;}

  visitfacets2(vertmat,n+1, numverts, faceverts, numverts, rk, faceverts, numverts, allfacets, &numfacets, debug);

  free_imatrix(vertmat,0,n,0,m1-1);
  return numfacets;
}


////////////////////////////////



int getfacestruct1(int **Sl, int n, int m, bool *verts, bool ***allfacets, int *rk, int **rankvec, int *totverts, int debug){
  int m1;
  int i, j, j1, k, flag, jflag;
  int numverts=getvertboolT(Sl, n, m, verts);//get the vertices
  bool v[numverts];
  bool greylist[numverts];
  bool faceverts[numverts];
  inittoone(faceverts,numverts);//full
  int **vertmat=conemat(Sl,n,m,verts,&m1);
  (*rk)=matrank(vertmat,n+1,numverts);

  int numfacets=0;//all bar the full polytope
  int newnumfacets;

  int sz, sz1;
  bool fullvec[n+1];
  bool fullvec1[numverts];
  int pos;

  bool **allfacetcands;//all candidates
  bool **tmpfacets;
  int totfacetcands=0;
  int tottmpfacets=0;
  int candtots[*rk+1];
  int curtot;
  bool *triagevec;
  int cmplt;

  (*totverts)=numverts;

  inittoone(fullvec, n+1);
  inittoone(fullvec1, numverts);

  if(debug){
    fprintf(stderr, "\n###Entering getfacestruct1. Matrix:\n");
    printmat(Sl, n,m);
    fprintf(stderr, "Cone matrix:\n");
    printmat(vertmat, n+1, numverts);
  }

  inittozero(v,numverts);

  //We first extract *candidates* for faces: maximal sets not intersecting the interior
  //In a later step we find which are faces

  candtots[0]=0;

  //add the vertices first
  for(k=0;k<numverts;k++){
    if(k>0){v[k-1]=0;}v[k]=1;
    addvec1(totfacetcands, v, numverts, &allfacetcands, 0, &pos);//no checking on vertices
    addnewto1Darray(&triagevec, totfacetcands, 1);
    totfacetcands++;
  }
  candtots[1]=numverts;

  for(k=1;k<(*rk)-1;k++){//Each rank
    // go through each existing face of rank k-1
    for(j=candtots[k-1];j<candtots[k];j++){

      jflag=1;//means base face is good to our knowledge so far
      inittozero(greylist,numverts);
      boolunion(greylist, allfacetcands[j], numverts);
      //augment with a single vertex in all possible ways
      for(i=0;i<numverts;i++){
	//The greylist is to avoid finding the same completions again and again
	//It must always at least contain the initial vertices: otherwise fundamental errors occur
	if(!greylist[i]){
	  veccp(allfacetcands[j], numverts, v);
	  sz=nonzentries(v, numverts);//current size
	  //fprintf(stderr, " building from:");printvec(v,numverts);
	  v[i]=1;//add a 1
	  //fprintf(stderr, "examining:");printvec(v,numverts);
	  //check if already a completion
	  if((cmplt=issubvec(v, allfacetcands, candtots[k], totfacetcands, numverts))==-1){//no existing completion
	    //doesn't intersect interior?
	    if(ispropsubface(vertmat, n+1, numverts, fullvec1, numverts, v, sz+1)){
	      //if(debug){fprintf(stderr, "found a subface: ");printvec(v,numverts);}//may not be maximal
	      if(k>1){
		sz1=growsubface(vertmat, n+1, numverts, v, sz+1, k+1);
		//if(debug){fprintf(stderr, "**Grown face, sz=%d, rank=%d\n",sz1,k+1);printvec(v,numverts);}

		//now check if allfacetcands[j] is a subface of new face
		if(!ispropsubface(vertmat, n+1, numverts, v, sz1, allfacetcands[j], sz)){//base face was invalid
		  // allfacetcands[j] was not a real face
		  triagevec[j]=0;
		  jflag=0;
		  //fprintf(stderr, "invalid base face\n");printvec(allfacetcands[j],numverts);
		  if(tottmpfacets){
		    free_bmat(tmpfacets,tottmpfacets);
		    tottmpfacets=0;
		  }
		  break;//finish with this base-face (not really a face)
		}

	      }
	      else
		sz1=sz+1;

	      //good so far: seems to be a face
	      tottmpfacets=addvec1(tottmpfacets, v, numverts, &tmpfacets, 1, &pos);
	      boolunion(greylist, v, numverts);//grow the greylist (only applies to current base-set)


	    }
	  }
	  //There is an existing completion, but it invalidates base-face
	  else if(!ispropsubface(vertmat, n+1, numverts, allfacetcands[cmplt], nonzentries(allfacetcands[cmplt], numverts), allfacetcands[j], sz)){
	    // allfacetcands[j] was not a real face
	    triagevec[j]=0;
	    jflag=0;
	    //fprintf(stderr, "invalid base face\n");printvec(allfacetcands[j],numverts);
	    if(tottmpfacets){
	      free_bmat(tmpfacets,tottmpfacets);
	      tottmpfacets=0;
	    }
	    break;//finish with this base-face (not really a face)
	  }
	}

      }

      //base face didn't intersect interior of any new ones (good news so far)
      if(jflag){//add the tempfacets to the main list
	for(i=0;i<tottmpfacets;i++){
	  newnumfacets=addvec1(totfacetcands, tmpfacets[i], numverts, &allfacetcands, 1, &pos);//penultimate means check if already there
	  //fprintf(stderr, "newnumfacets=%d\n",newnumfacets);
	  if(newnumfacets>totfacetcands){
//fprintf(stderr, "here: %d, %d, %d\n", i, tottmpfacets,totfacetcands);printvec(tmpfacets[i],numverts);exit(0);

	    addnewto1Darray(&triagevec, totfacetcands, 1);
	    totfacetcands=newnumfacets;
	    if(debug){fprintf(stderr, "**Potential faces so far (including this): %d. New candidate (rank %d, %d vertices):\n", totfacetcands, k+1, nonzentries(tmpfacets[i], numverts));printvec(tmpfacets[i],numverts);}
	  }
	}
	if(tottmpfacets){
	  free_bmat(tmpfacets,tottmpfacets);
	  tottmpfacets=0;
	}

      }
    }

    candtots[k+1]=totfacetcands;
    if(debug){fprintf(stderr, "candtots[%d] = %d\n", k+1, candtots[k+1]);}
  }

  //Add the whole polytope as the final facet
  totfacetcands=addvec1(totfacetcands, fullvec1, numverts, &allfacetcands, 0, &pos);
  candtots[k+1]=totfacetcands;
  //fprintf(stderr, "candtots[%d] = %d\n", k+1, candtots[k+1]);

  //triage: check if really facets

  //check and add in reverse order (don't add vertices or the whole polytope)
  for(k=(*rk)-1;k>=2;k--){
    curtot=0;
    //fprintf(stderr, "doing rank %d (%d to %d)\n", k, candtots[k-1], candtots[k]);
    //go through faces of this dimension
    for(j=candtots[k-1];j<candtots[k];j++){
      if(triagevec[j]){
	flag=1;
	//printvec(allfacetcands[j],numverts);
	//check all higher dimensional faces containing this one
	for(j1=candtots[k];j1<totfacetcands;j1++){
	  //(j1) is a face, contains (j) as a subset, but (j) isn't a proper subface of (j1)
	  if(triagevec[j1] && issubvec(allfacetcands[j], allfacetcands[j1], numverts) && !ispropsubface(vertmat, n+1, numverts, allfacetcands[j1], nonzentries(allfacetcands[j1], numverts), allfacetcands[j], nonzentries(allfacetcands[j], numverts))){
	    triagevec[j]=0;//not a face
	    flag=0;
	    break;//no need to check further
	  }
	}
	if(flag){//found a face
	  if(debug){fprintf(stderr, "**found a new proper rank %d face spanned by %d vertices:\n", k, nonzentries(allfacetcands[j], numverts));printvec(allfacetcands[j],numverts);}
	  numfacets=addvec1(numfacets, allfacetcands[j], numverts, allfacets, 0, &pos);//penultimate means don't check
	  curtot++;
	}
      }
    }
 
    addnewto1Darray(rankvec, ((*rk)-1-k), curtot);
  }

  if(debug)
    fprintf(stderr, "\n###Exiting getfacestruct1.\n\n");

  free((char*)triagevec);
  free_bmat(allfacetcands,totfacetcands);
  free_imatrix(vertmat,0,n,0,m1-1);
  return numfacets;
}

//returns number of triangles just added
int triangulate(bool *F, int m1, int rkF, bool **faces, int numfaces, int toprank, int *rankvec, bool ***Trs, int *totTrs, int *level, int debug){
  int pos;
  int v,i,i0,j;
  int newTrs;
  int totadded=0;
  (*level)++;
  if(debug){
    fprintf(stderr, "###\nEntering \"triangulate\" (level %d)\n", *level);
    fprintf(stderr, "now triangulating...\n");printvec(F,m1);
  }	
  if(rkF<=1 ||nonzentries(F,m1)==rkF+1){//1D or already simplicial
    (*totTrs)=addvec1((*totTrs), F, m1, Trs, 0, &pos);
    if(debug){fprintf(stderr, "(already triangulated)\n");}
    (*level)--;
    return 1;
  }
  v=firstnonz(F,m1);//first vertex

  if(rkF==toprank)
    i0=0;
  else{
    i0=0;
    for(j=0;j<toprank-rkF;j++)
      i0+=rankvec[j];
  }
  //fprintf(stderr, "i0=%d, i1=%d\n",i0,i0+rankvec[toprank-rkF]);exit(0);
  for(i=i0;i<i0+rankvec[toprank-rkF];i++){
    //fprintf(stderr, "(v=%d), i=%d: ",v,i);printvec(faces[i],m1);
    if(!nonz(faces[i],m1,v) && issubvec(faces[i],F,m1)){//subface doesn't include v
      if(debug){fprintf(stderr, "now checking (dim %d)...\n", rkF-1);printvec(faces[i],m1);}
      newTrs=triangulate(faces[i], m1, rkF-1, faces, numfaces, toprank, rankvec, Trs, totTrs, level, debug);
      if(debug){fprintf(stderr, "i=%d, Triangles so far: %d\n",i, *totTrs);}
      for(j=(*totTrs)-newTrs;j<(*totTrs);j++){//add in v
	(*Trs)[j][v]=1;
      }
      totadded+=newTrs;
    }
  }

  (*level)--;
  return totadded;
  
}


//Minkowski sum of two sets of points: return vertices of the convex hull
//in a matrix. Put the number of vertices in numverts
//Can set b1, b2 to be NULL: means use all columns
int **Minkowski2(int **Sl1, int n, int m1, bool *b1, int **Sl2, int m2, bool *b2, int *numverts){
  int **tmpmat=imatrix(0, n, 0, m1*m2);
  int **tmpmat2;
  int m0;
  int pat[m1*m2]; //not used
  int j1, j2, i, tot=0;
  int **vertmat;
  for(j1=0;j1<m1;j1++){
    for(j2=0;j2<m2;j2++){
      if((!b1 || b1[j1]) && (!b2 || b2[j2])){
	for(i=0;i<n;i++){
	  tmpmat[i][tot]=Sl1[i][j1]+Sl2[i][j2];
	}
	tot++;
      }
    }
  }
  //printmat(tmpmat,n,tot);
  tmpmat2=redmat1(tmpmat,n,tot,&m0,pat);//remove repeated cols
  //printmat(tmpmat2,n,m0);
  
  bool verts[m0];
  (*numverts)=getvertboolT(tmpmat2, n, m0, verts);//vertex indices
  vertmat=submatfromcolbool(tmpmat2, n, m0, verts, numverts);//vertices matrix
  //printmat(vertmat,n,*numverts);
  free_imatrix(tmpmat,0,n-1,0,m1*m2-1);
  free_imatrix(tmpmat2,0,n-1,0,m0-1);

  return vertmat;

}

//Minkowski sum of numR poltopes constructed from a common set of vertices
//vertices are columns of Sl
//Each row of b1 defines a polytope
int **MinkowskiN(int **Sl, int n, int m, bool **b, int *yc, int k, int *numverts){
  int i,m0;
  int mtmp[k-1];
  int **Stmp[k-1];
  int pat[m]; //not used
  int **Sl0, **Sl1;
  int **vertmat;
  int len;
  if(k<1){
    fprintf(stderr, "\"MinkowskiN\" requires at least 1 polytope. EXITING.\n");exit(0);
  }

  if(k==1){//Just get the vertices of convex hull
    Sl0=submatfromcolbool(Sl, n, m, b[yc[0]], &len);
    Sl1=redmat1(Sl0,n,len,&m0,pat);//remove repeated cols
    bool verts[m0];
    (*numverts)=getvertboolT(Sl1, n, m0, verts);//vertex indices
    vertmat=submatfromcolbool(Sl1, n, m0, verts, numverts);//vertices matrix
    free_imatrix(Sl0,0,n-1,0,len-1);
    free_imatrix(Sl1,0,n-1,0,m0-1);
    return vertmat;
  }
  //printmat(b,n,m);exit(0);
  Stmp[0]=Minkowski2(Sl, n, m, b[yc[0]], Sl, m, b[yc[1]], mtmp+0);
  for(i=1;i<k-1;i++){
    Stmp[i]=Minkowski2(Stmp[i-1], n, mtmp[i-1], NULL, Sl, m, b[yc[i+1]], mtmp+i);
  }

  for(i=0;i<k-2;i++)
    free_imatrix(Stmp[i],0,n,0,mtmp[i]);
  (*numverts)=mtmp[k-2];
  return Stmp[k-2];
}


// Get the volume of a polytope (n-dim volume where n is the
// dimension of the ambient space)
int PolytopeVol(int **S, int **Sl, int n, int m, int debug){
  int i, j, m0,m1, numfaces;
  //faces of dimension 2 to rk-1
  bool **allfaces=NULL;
  int rank, *rankvec;
  int pat[m];
  int **Sl0;
  int debugfull=(debug<=0)?0:debug-1;
  bool **Trs;
  int totTrs=0;
  int V=0;
  int totverts;
  int **Sl1;
  int level=0;
    
  if(debug){fprintf(stderr, "\n###Entering PolytopeVol\n");}
  //fprintf(stderr, "%d, %d, %d\n",n,m,affrank(Sl,n,m));
  if(affrank(Sl,n,m)<n){
    if(debug){fprintf(stderr, "Lower dimensional polytope: Volume = 0\n###Exiting PolytopeVol\n");}
    return 0;
  }

  Sl0=redmat1(Sl,n,m,&m0,pat);//remove repeated cols (store the pattern)
  bool verts[m0];
  totverts=getvertboolT(Sl0, n, m0, verts);//get the vertices
  Sl1=submatfromcolbool(Sl0, n, m0, verts, &m1);//keep only vertices

  
  //subtract first vertex from all vertices
  for(j=1;j<m1;j++){
    for(i=0;i<n;i++){
      Sl1[i][j]-=Sl1[i][0];
    }
  }
  for(i=0;i<n;i++)
    Sl1[i][0]=0;

  numfaces=getfacestruct1(Sl1, n, m1, verts, &allfaces, &rank, &rankvec, &totverts, debugfull);
    
  if(debug){fprintf(stderr, "Points:\n");printmat(Sl1,n,m1);}

  triangulate(verts, m1, rank-1, allfaces, numfaces, rank-1, rankvec, &Trs, &totTrs, &level, debugfull);

  for(i=0;i<totTrs;i++){
    if(debug){printvec(Trs[i],m1);}
    Trs[i][0]=0;//base vertex
    V+=abs(detsubmat(Sl1, n, m1, NULL, Trs[i]));
    if(debug){fprintf(stderr, "Current total Vol = %d/%d\n", V, (int)(factorial(rank-1)));}
  }

    
  //fprintf(stderr, "Total (%d)-Vol = %d/%d\n", rank-1, V, ((int)(factorial(rank-1))));
  free_imatrix(Sl0, 0, n-1, 0, m0-1);
  free_imatrix(Sl1, 0, n-1, 0, m1-1);
  free_bmat(allfaces, numfaces);
  free_bmat(Trs, totTrs);
  if(numfaces)
    free((char*)rankvec);

  if(debug){fprintf(stderr, "\n###Exiting PolytopeVol; Vol=%d/(%d!)\n",V,n);}

  return V;
}


