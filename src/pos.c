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

#include "pos.h"
#include "basics.h"
#include "symbolic.h"
#include "convex.h"
#include "semidef.h"
#include <limits.h>

int AMGMsimp(ex detex, int numv, int debug){
  size_t j,k,k1;
  ex prod=1;
  int hnv;
  int allflg;
  if(is_a<add>(detex) && detex.nops()==4){
    if((hnv=NewtonPolyVertSgn1(detex, numv, &allflg, debug))==-1 && !allflg){// all vertices negative
      for(k=0;k<4;k++){
	for(j=0;j<(detex.op(k)).nops();j++){
	  if(is_a<numeric>(detex.op(k).op(j)) && (ex_to<numeric>(detex.op(k).op(j))==3)){
	    for(k1=0;k1<4;k1++){
	      if(k1!=k)
		prod*=detex.op(k1);
	    }
	    prod+=pow((detex.op(k))/3, 3);
	    if(prod==0){
	      if(debug){fprintf(stderr, "Nonpositive by simple AM-GM.\n");}
	      return -1;
	    }
	  }
	}
      }
    }
    else if((hnv=NewtonPolyVertSgn1(detex, numv, &allflg, debug))==1 && !allflg){// all vertices positive
      for(k=0;k<4;k++){
	for(j=0;j<(detex.op(k)).nops();j++){
	  if(is_a<numeric>(detex.op(k).op(j)) && (ex_to<numeric>(detex.op(k).op(j))==3)){
	    for(k1=0;k1<4;k1++){
	      if(k1!=k)
		prod*=detex.op(k1);
	    }
	    prod+=pow((detex.op(k))/3, 3);
	    if(prod==0){
	      if(debug){fprintf(stderr, "Nonnegative by simple AM-GM.\n");}
	      return 1;
	    }
	  }
	}
      }
    }
  }
  return 0;
}

//
// Naive attempt to "dominate" negative monomials
// (Requires a homogeneous polynomial)
// Aimed at *strict* positivity
// tmp is the polynomial assumed to be in standard form
// (variables consist of a single letter followed by a number beginning at 1)
// modified so it outputs a decomposition if successful
//
// (numv+1) appears at several points in this algorithm which seems wasteful.
// TO DO: check why and remove if necessary

int J2pI_recurr(ex tmp, ex orig_poly, ex rems[], int numv, int *level, int *totiter, float *bestq, ex *best, int *lastsplitfail, int maxlevel, int maxcliques, int maxiter, int debug){
 
  int dim;
  int **lst,**explst,*cflst;
  long r,s;
  float ns;
  int *needcount,*needexcl;
  int **needs, **neededby;
  ex tmp3;
  int negcount=0;
  int maxsz=20000;
  int tight=0;//loose
  long **left=NULL,**right=NULL;
  int *numleft=NULL,*numright=NULL,*absnegmax=NULL;
  int totcliques=0,j;
  float totalq=0;
  int retstat;
  int debugfull=(debug<=0)?0:debug-1;

  if(*level>=maxlevel){
    if(debug){cerr << "maximum level reached: failure.\n";}
    return 0;
  }
  if(*totiter>=maxiter){
    if(debug){cerr << "maximum iterations reached: failure\n";}
    return -1;
  }

  (*level)++;(*totiter)++;
  if(debugfull){fprintf(stderr, "level = %d\n", *level);cerr << "tmp = " << tmp << endl;}
  dim=polytointlist(tmp, &lst, &cflst, &r);//dim=order of monomials

  //
  //cout << "dim = " << dim << endl;
  //printvec1(lst[0],numv);
  //treeprint3(cflst,lst,dim,r,numv,0);
  for(s=0;s<r;s++){
    if(cflst[s]<0)
      negcount++;
  }
  if(negcount==0){
    if(tmp!=0){
      if(debug){
	fprintf(stderr, "Success: positive. Certificate:\n");
	cerr << orig_poly << " = " << tmp;
	for(j=0;j<(*level)-1;j++)
	  cerr << "+(" << rems[j] << ")";
	cerr << endl << endl;
      }
      ifree(lst, r);free((char *)cflst);
      return 1;
    }
    else{
      if(debug){
	fprintf(stderr, "Partial success: nonnegative. Certificate:\n");
	cerr << orig_poly << " = " << tmp;
	for(j=0;j<(*level);j++)
	  cerr << "+(" << rems[j] << ")";
	cerr << endl << endl;
      }
      ifree(lst, r);free((char *)cflst);
      return -2;
    }
  }

  explst=moncflisttoexplist(lst,dim,cflst,r,numv+1,1);// redefines cflst
  ifree(lst, r);

  /* if(*totiter==1 && debug)// first time only */
  /*   listprint(cflst, explst, numv+1, r); */


  needcount = (int *)malloc((size_t) (r*sizeof(int))); 
  needexcl = (int *)malloc((size_t) (r*sizeof(int))); 
  for(s=0;s<r;s++){needcount[s]=0;needexcl[s]=0;}
    
  needs=(int **)malloc((size_t)(r*sizeof(int*)));
  neededby=(int **)malloc((size_t)(r*sizeof(int*)));
  for(s=0;s<r;s++){
    neededby[s]=(int *)malloc((size_t)(sizeof(int)));neededby[s][0]=0;
    needs[s] = (int*)malloc(sizeof(int));needs[s][0]=0; //number needed
  }
  for(s=0;s<r;s++){
    if(cflst[s]<0){
      if(*totiter==1 && debug)
	ns=numsplits(s,numv+1,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);// changed the final flag to make less verbose
      else
	ns=numsplits(s,numv+1,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);
      if(ns<0.0001){// not able to split any more
	if(debugfull){cerr << "failure to split:\n   " << cflst[s] << ", "; printvec(explst[s],numv+1);}
	veccp(explst[s],numv+1,lastsplitfail);
	(*level)--;
	ifree(needs, r);ifree(neededby, r);ifree(explst, r);
	free((char *)cflst); free((char *)needcount);free((char *)needexcl);
	//this is problematic: a failure to split need not be fatal
	//Consider changing
	if(*totiter==1)
	  return -1; // fatal
	else
	  return 0; // exit, but not fatal
      }
      totalq+=ns;
      //	cout << "numsplits[" <<  s << "]: " << ns << endl;
    }
  }
  totalq/=((float)negcount);
  if(debug){fprintf(stderr, "tot: %d, neg count: %d\n", (*totiter), negcount);}
  if(*totiter>1 && totalq>(*bestq)){
    (*bestq)=totalq; 
    (*best)=tmp;
    //if(debug){cerr << "best q: " << *bestq << endl;}
  }


  /* tmp3=anyforced(explst, cflst, numv+1, r, needexcl, needs, neededby,1); */
  /* cout << "forced clique:" << tmp3 << endl; */
  /* if(tmp3!=0 && (retstat=J2pI_recurr(tmp-expand(tmp3), orig_poly, rems, numv, level, totiter, bestq, best, lastsplitfail, maxlevel, maxcliques, maxiter, debug))){// success or fatal failure */
  /*   if(left){lfree_i(left,totcliques);}if(right){lfree_i(right,totcliques);} */
  /*   ifree(explst, r);free((char *)cflst);  */
  /*   if(numleft){free((char*)numleft);}if(numright){free((char*)numright);}if(absnegmax){free((char*)absnegmax);} */
  /*   return retstat; */
  /* } */

  totcliques=allcliques1(explst,cflst,r,numv+1,neededby,needs,needexcl,&left,&right,&numleft,&numright,&absnegmax,maxsz,maxcliques,tight,0,1,0);//quiet

  if(debug){// list all good cliques
    for(j=0;j<totcliques;j++){
      tmp3=makesquare(explst, cflst, left[j],numleft[j],right[j],numright[j],numv+1,absnegmax[j],0);
      cerr << "clique found: " << tmp3 << endl;
    }
  }


  if(totcliques==0){ // try without forcing
    maxcliques=20;
    totcliques=allcliques1(explst,cflst,r,numv+1,neededby,needs,needexcl,&left,&right,&numleft,&numright,&absnegmax,maxsz,maxcliques,tight,0,0,0);//quiet
    if(debug){// list all good cliques
      for(j=0;j<totcliques;j++){
	tmp3=makesquare(explst, cflst, left[j],numleft[j],right[j],numright[j],numv+1,absnegmax[j],0);
	cerr << "not tight clique found: " << tmp3 << endl;
      }
    }
  }
  ifree(neededby, r);ifree(needs, r);
  free((char *)needcount);free((char *)needexcl);

  for(j=0;j<totcliques;j++){
    tmp3=makesquare(explst, cflst, left[j],numleft[j],right[j],numright[j],numv+1,absnegmax[j],0);
    if(debug){cerr << "numleft = " << numleft[j] << ", numright = " << numright[j] << ", level = " << *level << endl;
      cerr << endl << "Removing: " << tmp3 << endl << "from\n" << tmp << endl << "leaving\n" << expand(tmp-tmp3) << endl;}
    rems[*level-1]=factor(tmp3);
    if((retstat=J2pI_recurr(tmp-expand(tmp3), orig_poly, rems, numv, level, totiter, bestq, best, lastsplitfail, maxlevel, maxcliques, maxiter, debug))){// success or fatal failure
      if(left){lfree_i(left,totcliques);}if(right){lfree_i(right,totcliques);}
      ifree(explst, r);free((char *)cflst); 
      if(numleft){free((char*)numleft);}if(numright){free((char*)numright);}if(absnegmax){free((char*)absnegmax);}
      return retstat;
    }
  }

  // freeing is very important - deeply recursive

  if(left){lfree_i(left,totcliques);}if(right){lfree_i(right,totcliques);}
  ifree(explst, r);free((char *)cflst);
  if(numleft){free((char*)numleft);}if(numright){free((char*)numright);}if(absnegmax){free((char*)absnegmax);}

  //cout << tmp << endl;

  (*level)--;

  return 0; // exit - go back up a level
}


//
// Wrapper for the heuristic attempt to write the polynomial
// as a sum of squares and positive monomials
// Requires a homogeneous expanded polynomial in standard form (v1,v2,...)
// Return codes match ispospoly
//

int heuristic_squares(ex detex, int numv, int debug){
  int level=0, totiter=0;
  int ret;
  ex best;
  float bestq;
  int lastsplitfail[numv];
  int maxcliques=1000, maxiter=1000, maxlevel=20;
  ex detx;
  int hnv0,allflg;
  ex rems[maxlevel];//squares removed
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering heuristic_squares\n");}

  if(detex==0){
    if(debug)
      fprintf(stderr, "The polynomial is identically zero.\n");
    return 0;
  }

  hnv0=NewtonPolyVertSgn1(detex, numv, &allflg, debugfull);
  if(hnv0==-1)
    detx=hnv0*detex;
  else if(hnv0==1)
    detx=detex;
  else
    return -3;//definitely mixed sign

  ret= J2pI_recurr(detx, detx, rems, numv, &level, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, debugfull);
  if(ret==1)//success (definite)
    return hnv0*2;
  else if(ret==-2)//partial success (semidefinite)
    return hnv0;

  if(debug){fprintf(stderr, "Exiting heuristic_squares\n");}
  return -4;//failure
}


// check if a polynomial is positive on the strictly positive orthant
// assumed to be expanded, homogeneous and in standard form (v1, v2, ...)
// 0 means this is the zero polynomial
// 1 means nonnegative, but not necessarily positive
// 2 means positive
// -1 means nonpositive, but not necessarily negative
// -2 means negative
// -3 means definitely can be both positive and negative
// -4 means failure
// allflg gets set to 1 if all terms of the same sign (otherwise zero)
// doublecfs means that the coeffs are not necessarily integers, so skip J2pI_recurr
int ispospoly(ex detx1, int numv, int *allflg, bool isfct, int doublecfs, int maxpppdeg, int *pppdegused, int debug){
  int ret1,hnv0,hnv;
  int level=0, totiter=0;
  int lastsplitfail[numv];
  int maxcliques=1000, maxiter=1000, maxlevel=20;
  ex best;
  float bestq;
  ex fct;
  ex prod=1;
  int tpn;
  int ldeg=2, ldeg1=2;
  //int amgm;
  int tryppp=1;
  int ppppreprocess=1;
  (*allflg)=0;
  ex detex, detx;
  int powval=0;
  ex rems[maxlevel];
  int debugfull=(debug<=0)?0:debug-1;


  (*pppdegused)=-1;

  if(maxpppdeg<0)
    tryppp=0;

  if(debug){fprintf(stderr, "\n###Entering ispospoly...\n");}

  detx=gcdpoly(detx1);
  if(debug && expand(detx-detx1)!=0)
    cerr << "Divided through by gcd of coefficients to get " << detx << endl;

  inittozero(lastsplitfail, numv);

  if(detx==0){
    if(debug)
      fprintf(stderr, "The polynomial is identically zero.\n");
    return 0;
  }
  //Do this at the outset and if necessary replace poly with -poly
  if((hnv0=NewtonPolyVertSgn1(detx, numv, allflg, debug))==2){//vertices of both signs
    if(debug){
      fprintf(stderr, "The Newton polytope of this polynomial has vertices with both positive and negative coefficients. It can be both positive and negative on the positive orthant.\n");
      cerr << detx << endl;
    }
    return -3;
  }

  if(!doublecfs){
    if(hnv0==-1)
      detex=getroot(-detx, &powval);
    else
      detex=getroot(detx, &powval);
  }
  else{
    if(hnv0==-1)
      detex=-detx;
    else
      detex=detx;
  }

  if(debug && powval){
    if(hnv0==1)
      cerr << "Examining the root: " << detex << " of " << detx << endl;
    else
      cerr << "Examining the resigned root: " << detex << " of " << detx  << endl;
  }
  else if(debug){
    if(hnv0==1)
      cerr << "Examining: " << detex << endl;
    else
      cerr << "Examining the resigned polynomial: " << detex << " (original = " << detx << ")" << endl;
  }

  if(powval){// nontrivial power
    if(powval%2==0){//positive power
      if(debug){
	if(hnv0==1)
	  fprintf(stderr, "The polynomial is an even power of another, so it is definitely nonnegative. Examining the root.\n");
	else
	  fprintf(stderr, "The polynomial is of the form -(P)^(2k), so it is definitely nonpositive. Examining the root.\n");
      }
    }
    else{
      if(debug)
	fprintf(stderr, "The polynomial is %san odd power of another: examining the root.\n",hnv0==1?"":"minus ");
    }
  }

  //Working with the original polynomial, or the resigned polynormial, or the root (if powval!=0), possibly after resigning.
  if((hnv=NewtonPolyVertSgn1(detex, numv, allflg, debug))==2){//vertices of both signs
    if(debug){
      fprintf(stderr, "The Newton polytope of the %spolynomial has vertices with both positive and negative coefficients. It can be both positive and negative on the positive orthant.\n", powval?"root ":"");
      //cerr << detex << endl;
    }
    if(powval && powval%2==0)
      return hnv0;
    return -3;
  }

  if(hnv==1 && (*allflg)){
    if(debug){
      fprintf(stderr, "The %spolynomial is strictly positive on the positive orthant (all terms are positive).\n", powval?"(resigned) root ":"");
      //cerr << detex << endl;
    }
    return hnv0*2;
  }

  if(hnv==-1 && (*allflg)){
    if(debug){
      fprintf(stderr, "The %spolynomial is strictly negative on the positive orthant (all terms are negative).\n", powval?"(resigned) root ":"");
      //cerr << detex << endl;
    }
    if(powval && powval%2==0)
      return hnv0*2;
    return -hnv0*2;
  }

  if(!doublecfs && hnv==1 && (ret1=J2pI_recurr(detex, detex, rems, numv, &level, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, debug))==1){// success
    if(debug){
      fprintf(stderr, "The %spolynomial is strictly positive on the positive orthant (but not all terms are positive).\n", powval?"(resigned) root ":"");
      //cerr << detex << endl;
    }
    return hnv0*2;
  }

  //nonnegative but not necessarily positive
  if(!doublecfs && hnv==1 && ret1==-2){
    if(tryppp && ppp(detex, ppppreprocess, maxpppdeg, 1, &ldeg,debugfull)==0){//success - definite
      if(debug)
	fprintf(stderr, "According to csdp the %spolynomial is positive on the positive orthant (left degree = %d).\n", powval?"(resigned) root ":"", ldeg);
      (*pppdegused)=ldeg;
      return hnv0*2;
    }
    if(debug)
      fprintf(stderr, "The %spolynomial is nonnegative on the positive orthant.%s\n", tryppp?"":" (Not trying csdp.)", powval?"(resigned) root ":"");
    return hnv0*1;
  }
  
  
  if(!doublecfs && hnv==-1 && (ret1=J2pI_recurr(-detex, -detex, rems, numv, &level, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, debug))==1){
    if(debug){
      fprintf(stderr, "The %spolynomial is strictly negative on the positive orthant (but not all terms are negative).\n", powval?"(resigned) root ":"");
      cerr << detex << endl;
    }
    if(powval && powval%2==0)
      return hnv0*2;
    return -hnv0*2;
  }

  if(!doublecfs && hnv==-1 && ret1==-2){
    if(tryppp && ppp(hnv*detex, ppppreprocess, maxpppdeg, 1, &ldeg,debugfull)==0){//success - definite
      if(debug)
	fprintf(stderr, "According to csdp the %spolynomial is negative on the positive orthant (left degree = %d).\n", powval?"(resigned) root ":"", ldeg);
      (*pppdegused)=ldeg;
      if(powval && powval%2==0)
	return hnv0*2;
      return -hnv0*2;
    }
    if(debug)
      fprintf(stderr, "The %spolynomial is nonpositive on the positive orthant.%s\n", tryppp?"":" (Not trying csdp.)", powval?"(resigned) root ":"");
    if(powval && powval%2==0)
      return hnv0*1;
    return -hnv0*1;
  }

  if(!(tpn=testposneg(detex, 0.001, 1.0, 10000, 1e-10))){//+-values found
    if(debug){fprintf(stderr, "All vertices of the Newton polytope of the %spolynomial are %s, but numerical evaluation gives both positive and negative values on the positive orthant.\n", powval?"(resigned) root ":"", hnv==1?"positive":"negative");
    cerr << detex << endl;}
    if(powval && powval%2==0)
      return hnv0*1;
    return -3;
    // testposneg: returns 0 if both positive and negative values are found
    // returns -1 if only nonpositive values found
    // returns 1 if only nonnegative values found

  }

  /* if(!doublecfs && (amgm=AMGMsimp(detex, numv, debug))){ */
  /*   if(debug) */
  /*     fprintf(stderr, "The polynomial is %s by a simple arithmetic mean-geometric mean argument.\n", amgm==1?"nonnegative":"nonpositive"); */
  /*   return amgm; */
  /* } */


  //failure so far:couldn't write as SOS+ or -SOS+ 
  if(debug){fprintf(stderr, "Numerics suggests that the %spolynomial is %s on the positive orthant.", powval?"(resigned) root ":"", tpn==1?"nonnegative":"nonpositive");if(tryppp){fprintf(stderr, " Trying to prove this with semidefinite programming upto degree %d.\n",maxpppdeg);}else{fprintf(stderr, " Not trying semidefinite programming.\n");}}

  if(tryppp && ppp(hnv*detex, ppppreprocess, maxpppdeg, 0, &ldeg1,debugfull)==0){//success - at least semidefinite
//exit(0);
    if(ppp(hnv*detex, ppppreprocess, maxpppdeg, 1, &ldeg,debugfull)==0){//success - definite
      if(debug)
	fprintf(stderr, "According to csdp the %spolynomial is %s on the positive orthant (left degree = %d).\n", powval?"(resigned) root ":"", hnv==1?"positive":"negative", ldeg);
      (*pppdegused)=ldeg;
      if(powval && powval%2==0)
	return hnv0*2;
      return hnv*hnv0*2;
    }
    else{
      if(debug){
	fprintf(stderr, "According to csdp the %spolynomial is %s on the positive orthant (left degree = %d).\n", powval?"(resigned) root ":"", hnv==1?"nonnegative":"nonpositive", ldeg1);
	//cerr<<detex<<endl;
      }
      (*pppdegused)=ldeg1;
      if(powval && powval%2==0)
	return hnv0*1;
      return hnv0*hnv;
    }
  }

  return -4;
  
}


// Version of ispospoly where we already know the polynomial is nonnegative
// E.g. For polynomial det(J^2+wI) [use with care]
// allflg gets set to 1 if all terms of the same sign (otherwise zero)
// doublecfs means that the coeffs are not necessarily integers, so skip J2pI_recurr
// return codes: positive (2), definitely can be zero (0), 
// failure - nonnegative but can't determine positivity (1), 
// failure, e.g. not satisfying premise of nonnegativity (-4)

int ispospolynonneg(ex detx, int numv, int *allflg, bool isfct, int doublecfs, int maxpppdeg, int *pppdegused, int debug){
  int ret1,hnv0,hnv;
  int level=0, totiter=0;
  int lastsplitfail[numv];
  int maxcliques=1000, maxiter=1000, maxlevel=20;
  ex best;
  float bestq;
  ex fct;
  ex prod=1;
  int ldeg=2;
  int tryppp=1;
  int ppppreprocess=1;
  (*allflg)=0;
  ex detex;
  int powval=0;
  ex rems[maxlevel];
  int debugfull=(debug<=0)?0:debug-1;

  (*pppdegused)=-1;

  if(debug){fprintf(stderr, "\n###Entering ispospolynonneg...\n\n");}

  inittozero(lastsplitfail, numv);

  if(detx==0){
    if(debug)
      fprintf(stderr, "The polynomial is identically zero.\n");
    return 0;
  }

  //Do this at the outset
  if((hnv0=NewtonPolyVertSgn1(detx, numv, allflg, debugfull))!=1){//vertices not all +ve
    if(debug)
      fprintf(stderr, "Not all vertices of the Newton polytope are positive. This algorithm is only for polynomials known to be nonnegative.\n");
    return -4;
  }
  if((*allflg)){
    if(debug)
      fprintf(stderr, "The polynomial is strictly positive on the positive orthant (all terms are positive).\n");
    return 2;
  }

  //Now assume hnv0=1

  if(!doublecfs)
    detex=getroot(detx, &powval);
  else
    detex=detx;

  if(debug && !doublecfs && powval)
    cerr << "Examining the root: " << detex << endl;


  //Working with the original polynomial, or the root (if powval!=0)
  if((hnv=NewtonPolyVertSgn1(detex, numv, allflg, debugfull))==2){//vertices of both signs
    if(powval && powval%2==0){//even power of indefinite poly
      if(debug)
	fprintf(stderr, "The Newton polytope of the root polynomial has vertices with both positive and negative coefficients. It can be both positive and negative on the positive orthant.\n");
      return 0; // definitely can be zero
    }
    else if(powval && powval%2!=0){//odd power of indefinite poly (shouldn't occur!)
      fprintf(stderr, "The Newton polytope of the root polynomial has vertices with both positive and negative coefficients. It can be both positive and negative on the positive orthant.\n");
      return -4; // fundamental failure
    }
    //no other possibilities should occur
  }

  // No need to check allflg - already done for the whole poly
  if(!doublecfs && (ret1=J2pI_recurr(detex, detex, rems, numv, &level, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, debugfull))==1){// success
    if(debug){
      fprintf(stderr, "The %spolynomial is strictly positive on the positive orthant (but not all terms are positive).\n", powval?"root ":"");
      cerr << detex << endl;
    }
    return 2;
  }

  if(tryppp && ppp(hnv*detex, ppppreprocess, maxpppdeg, 1, &ldeg,debugfull)==0){//success - definite
    if(debug)
      fprintf(stderr, "According to csdp the %spolynomial is positive on the positive orthant (left degree = %d).\n", powval?"(resigned) root ":"", ldeg);
    (*pppdegused)=ldeg;
    return 2;
  }

  return 1;
}

//overloading: default is maxpppdeg=0
int ispospolynonneg(ex detx, int numv, int *allflg, bool isfct, int doublecfs, int *pppdegused, int debug){
  return ispospolynonneg(detx, numv, allflg, isfct, doublecfs, 0, pppdegused, debug);
}

//deprecated

int ispospolynonneg_old(ex detex, int n, int numv, int *allflg, bool isfct, int doublecfs, int debug){
  int ret=-4,ret1,j,nng=0,hnv,sgn=1;
  int level=0, totiter=0;
  int lastsplitfail[numv];
  int maxcliques=1000, maxiter=1000, maxlevel=20;
  ex best;
  float bestq;
  ex fct;
  size_t k1;
  ex prod=1;
  int ldeg=2;
  int maxpppdeg=2;
  int tryppp=1;
  int ppppreprocess=1;
  (*allflg)=0;
  ex rems[maxlevel];

  inittozero(lastsplitfail, numv);

  if(detex==0){
    if(debug)
      fprintf(stderr, "The polynomial is identically zero.\n");
    return 0;
  }

  //if((hnv=hasnegvertex(detex,numv))==2){//vertices of both signs
  if((hnv=NewtonPolyVertSgn1(detex, numv, allflg, debug))==2){//vertices of both signs
    fprintf(stderr, ".");
    if(debug)
      fprintf(stderr, "The Newton polytope of this polynomial has vertices with both positive and negative coefficients. The polynomial can be both positive and negative on the positive orthant.\n");
    return -3;
  }
  else if(hnv==1 && (*allflg)){
    if(debug){
      fprintf(stderr, "The polynomial is strictly positive on the positive orthant (all terms are positive).\n");
      cerr << detex << endl;
    }
    return 2;
  }
  else if(hnv==-1 && (*allflg)){
    if(debug){
      fprintf(stderr, "The polynomial is strictly negative on the positive orthant (all terms are negative).\n");
      cerr << detex << endl;
    }
    return -2;
  }
  // made these quiet - otherwise a lot of output
  else if(!doublecfs && hnv==1 && (ret1=J2pI_recurr(detex, detex, rems, numv, &level, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, 0))==1){// success
    if(debug){
      fprintf(stderr, "The polynomial is strictly positive on the positive orthant (but not all terms are positive).\n");
      cerr << detex << endl;
    }
    if(isfct)
      cout << "J2pI_recurr successful on this factor\n";
    else
      cout << "J2pI_recurr successful\n";
    return 2;
  }
  else if(!doublecfs && hnv==-1 && (ret1=J2pI_recurr(-detex, -detex, rems, numv, &level, &totiter, &bestq, &best, lastsplitfail, maxlevel, maxcliques, maxiter, 0))==1){
    if(debug){
      fprintf(stderr, "The polynomial is strictly negative on the positive orthant (but not all terms are negative).\n");
      cerr << detex << endl;
    }
    if(isfct)
      cout << "J2pI_recurr successful on this factor\n";
    else
      cout << "J2pI_recurr successful\n";
    return -2;
  }
  //failure so far:couldn't write as SOS+ or -SOS+ 
  //and doesn't have mixed vertices. Factorise and try again.
  else{
    if(isfct)
      cout << "J2pI_recurr unsuccessful on this factor\n";
    else
      cout << "J2pI_recurr unsuccessful\n";

    if(tryppp && (isfct || is_a<add>(factor(detex))) && ppp(hnv*detex, ppppreprocess, maxpppdeg, 1, &ldeg,0)==0){//success
      if(debug){
	fprintf(stderr, "According to csdp this polynomial is %s on the positive orthant (left degree = %d).\n", hnv==1?"positive":"negative", ldeg);
	//cerr<<detex<<endl;
      }
      if(isfct)
	cout << "ppp successful (degree=" << ldeg << ") on this factor\n";
      else
	cout << "ppp successful (degree=" << ldeg << ")\n";

      return 2*hnv;
    }

    //To see what happened
    //if hnv==1 expect positive; if hnv==-1, expect negative
    if(!isfct){//not already a factor of something
      if(debug){cout << "trying factorisation.\n";}
      //fprintf(stdout, "something\n");
      fct=factor(detex);
      //weird behaviour: the next line causes ginac crashes. why?
      //cout << "some stuff " << fct << endl;
      //if(debug){cout << "factorisation: " << fct <<endl;}
      if(!is_a<add>(fct) && fct.nops()>1){//nontrivial factorisation
	if(debug){cout << "number of terms: " << fct.nops() << endl;}
	for (k1=0; k1!=fct.nops(); ++k1){
	  if(debug){cout << "factor " << k1 << ": " << fct.op(k1) << endl;}
	  //note: make this silent (confusing if success being reported for factors)
	  if((j=ispospolynonneg_old(expand(fct.op(k1)), n, numv, allflg, 1, doublecfs, 1))<=-3){//factor of undetermined sign
	    if(debug){fprintf(stdout, "factor of undetermined sign.\n");}
	    return j;
	  }
	  else if(j==1 || j==-1)
	    nng=1;//may be zero
	  sgn*=j;
	  /* cout << "k1 = " << k1 << ", sgn = " << sgn << endl; */
	  /* cout << "hnv = " << hnv << endl; */
	  //      cout << fct.op(k1) << endl;
	}
	if(sgn>0 && nng) // nonnegative
	  ret=1;
	else if(sgn>0) //strictly positive
	  ret=2;
	else if(sgn<0 && nng) // nonpositive
	  ret=-1;
	else if(sgn<0) //strictly negative
	  ret=-2;
	if(debug){fprintf(stdout, "Factorising worked!\n");}
      }
      else // no factorisation: definite failure
	ret=-4;
    }
    else //already a factor in a factorised expression: definite failure
      ret=-4;
  }

  return ret;
  
}

//return the factorised form
ex canfactor(ex detex, unsigned int *tot, int debug){
  unsigned int k1;
  ex fct=factor(detex);
  if(debug){fprintf(stderr, "\n###Entering canfactor.\n");}

  (*tot)=fct.nops();
  if(!is_a<add>(fct) && (*tot)>1){//nontrivial factorisation
    if(debug){
      cerr << "number of terms: " << fct.nops() << endl;
      for (k1=0; k1!=(*tot); ++k1)
	cerr << "factor " << k1 << ": " << fct.op(k1) << endl;
    }
  }
  else//important otherwise (*tot) are just terms
    (*tot)=1;

  if(debug){fprintf(stderr, "Leaving canfactor: was%s able to factor.\n", (*tot)==1?" not":"");}

  return fct;
}


// For use on polynomials known to be nonnegative (otherwise just use
// NewtonPolyVertSgn1).
// Confirm via factorisation if possible whether a polynomial 
// can definitely be zero. Returninging 0 means inconclusive
// 1 means definitely takes the value zero at some point

int canbezero(ex detex, int numv, bool debug){
  unsigned int k1;
  int allflg, powval;
  unsigned int tot;
  ex fct=canfactor(detex, &tot, debug);
  int debugfull=(debug<=0)?0:debug-1;

  if(tot<=1)//no need for debugging (given by canfactor)
    return 0;

  if(debug){fprintf(stderr, "\n###Entering canbezero: examining polynomial. total factors=%d\n", tot); cerr << detex << endl;}


  for (k1=0; k1<tot; ++k1){
    //cerr << "here: " << fct.op(k1) << endl;
    detex=getroot(fct.op(k1), &powval);

    if(NewtonPolyVertSgn1(expand(detex), numv, &allflg, debugfull)==2){//vertices of both signs
      if(debug){fprintf(stderr, "Found a factor with mixed vertices. Polynomial can be zero.\n");}
      return 1;
    }
  }

  if(debug){fprintf(stderr, "Leaving canbezero: unable to determine via factorisation if the polynomial can be zero.\n");}
  return 0;
}


	

//Optimised version of isfullpospoly
//with minimal checks
//Assume that all vertices are positive squares
//and that we are checking for non-negativity only
//The polynomial is homogeneous of known degree deg
//SDP is run with left degree upto degree maxldeg (even)
// 1 means nonnegative
// -3 means definitely has negative values
// -4 means failure


int ifppopt(int **lst, double *cfs, int polytot, int numv, int deg, int maxldeg, int *ldeg){
  int preprocess=1;
  int ret=-1;
  int ldegree=0;
  int debug=0;
  (*ldeg)=0;

  if(!testforpos(lst, cfs, polytot, numv, -1.0, 1.0, 10000, 1e-10))
    return -3;

  while(ret!=0 && ldegree<=maxldeg){
    ret=pppfull1(lst, cfs, (long)polytot, numv, deg, deg, preprocess, ldegree, ldegree, 0, debug);
    ldegree+=2;
  }

  if(!ret){
    *ldeg=ldegree-2;
    return 1;
  }

  return -4;
}




int ifppopt1(ex tmp, int numv, int maxldegree, int strict, int *ldeg){
  int **lst=NULL;
  double *cflst;
  long r;
  int polymindeg,polymaxdeg;
  int ret=-1;
  int debug=0;

  //cerr << "entering ppp\n";
  //cerr << "numv = " << numv << endl << tmp << endl << tmp1<< endl;

  if(debug)
    cout << "ifppopt1, checking polynomial:\n" << tmp << endl;

  extolst1(tmp, numv, &lst, &cflst, &r, &polymindeg, &polymaxdeg);

  if(polymindeg!=polymaxdeg){
    fprintf(stderr, "ERROR in pppfull: polynomial is not homogeneous.\n");
    exit(0);
  }

  ret=ifppopt(lst, cflst, r, numv, polymindeg, maxldegree, ldeg);


  ifree(lst, r);
  free((char*)cflst);

  return ret;
}



// Requires a homogeneous polynomial
// check if a polynomial is positive on all of R^n\{0}
// 0 means this is the zero polynomial
// 1 means nonnegative, but not necessarily positive
// 2 means positive
// -1 means nonpositive, but not necessarily negative
// -2 means negative
// -3 means definitely can be both positive and negative
// -4 means failure
// nonnegonly means only test for nonnegativity (not strict positivity)

int isfullpospoly(ex detex1, int numv, bool isfct, int maxpppdeg, int *pppdegused, int nonnegonly, int debug){
  int hnv;
  ex fct;
  ex prod=1;
  int tpn;
  int ldeg=2, ldeg1=2;
  int tryppp=1;
  int allpossq=0;
  ex detex;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering isfullpospoly.\n");}

  detex=gcdpoly(detex1);
  if(debug && expand(detex-detex1)!=0)
    cerr << "Divided through by gcd of coefficients to get " << detex << endl;


  (*pppdegused)=-1;
  if(maxpppdeg<0)
    tryppp=0;

  if(detex==0){
    if(debug)
      fprintf(stderr, "The polynomial is identically zero.\n");
    return 0;
  }
  else if(is_a<numeric>(detex)){
    if(ex_to<numeric>(detex)>0){
      if(debug)
	fprintf(stderr, "The polynomial is a positive constant.\n");
      return 2;
    }
    else if(ex_to<numeric>(detex)<0){
      if(debug)
	fprintf(stderr, "The polynomial is a negative constant.\n");
      return -2;
    }
  }
  hnv=NewtonPolyVertPosSq(detex, numv);

  if(hnv==3){
    fprintf(stderr, "The Newton polytope of this polynomial has nonsquare vertices. The polynomial can be both positive and negative on R^n\\{0}.\n");
    return -3;
  }
  else if(hnv==2){//vertices of both signs
    fprintf(stderr, ".");
    if(debug)
      fprintf(stderr, "The Newton polytope of this polynomial has vertices with both positive and negative coefficients. The polynomial can be both positive and negative on R^n\\{0}.\n");
    return -3;
  }
  else if(((hnv==1) || (hnv==-1)) && !(tpn=testposneg(detex, -1.0, 1.0, 10000, 1e-10))){//+-values found
    fprintf(stderr, "All vertices of the Newton polytope are %s, but numerical evaluation gives both positive and negative values for this polynomial on R^n.\n", hnv==1?"positive":"negative");
    cerr << detex << endl;
    return -3;
  }
  else if(hnv==10){//all monomials +ve squares: definitely nonnegative
    if(debug){
      fprintf(stderr, "All monomials are positive squares. The polynomial is nonnegative on R^n.\n");
      cerr << detex << endl;
    }
    if(nonnegonly || !tryppp)
      return 1;
    allpossq=1;hnv=1;
  }
  else if(hnv==-10){//all monomials -ve squares: definitely nonpositive
    if(debug)
      fprintf(stderr, "All monomials are negative squares. The polynomial is nonpositive on R^n.\n");
    if(nonnegonly || !tryppp)
      return 1;
    allpossq=-1;hnv=-1;
  }


  if(debug && !allpossq){//only evaluates if we did numerics
    fprintf(stderr, "Numerics suggests that the polynomial is %s on R^n\\{0}.", tpn==1?"nonnegative":"nonpositive");
    if(tryppp)
      fprintf(stderr, " Trying to prove this with factorisation and semidefinite programming upto degree %d.\n",maxpppdeg);
    else
      fprintf(stderr, " Not trying semidefinite programming.\n");
  }

  //Either already factorised or no nontrivial factorisation
  if(tryppp){
    if((allpossq || pppfull(hnv*detex, numv, 1, maxpppdeg, 0, &ldeg1,debugfull)==0)){//success - at least semidefinite
//exit(0);
      if(!nonnegonly && pppfull(hnv*detex, numv, 1, maxpppdeg, 1, &ldeg,debugfull)==0){//success - definite
	if(debug)
	  fprintf(stderr, "According to csdp this polynomial is %s on R^n\\{0} (left degree = %d).\n", hnv==1?"positive":"negative", ldeg);
	(*pppdegused)=ldeg;
	return 2*hnv;
      }
      else{
	if(debug){
	  fprintf(stderr, "According to csdp this polynomial is %s on R^n\\{0} (left degree = %d).\n", hnv==1?"nonnegative":"nonpositive", ldeg1);
	  //cerr<<detex<<endl;
	}
	(*pppdegused)=ldeg1;
	return hnv;
      }
    }
  }

  return -4;
  
}

//
// Know number of variables: don't recompute
// allflg if all terms (not just vertices)
// are of the same sign
// Doesn't require the polynomial to be homogeneous
// but does require variable names to begin with subscript 1
// (problem with all routines depending on monotointlist)

int NewtonPolyVertSgn1(ex tmp, int numv, int *allflg, int debug){
  long r,s;
  int **lst;
  double *cflst;
  int polymindeg,polymaxdeg;
  int flag=0;
  int allpos=0;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering NewtonPolyVertSgn1\n");}
  //fprintf(stderr, "got here\n");fflush(stdout);fflush(stderr);
  (*allflg)=0;
  extolst1(tmp, numv, &lst, &cflst, &r, &polymindeg, &polymaxdeg);

  //fprintf(stderr, "got here_a\n");fflush(stdout);fflush(stderr);
  if(debugfull){cerr << tmp << endl;
    for(s=0;s<r;s++){
      fprintf(stderr, "%.4f: ", cflst[s]);
      printvec(lst[s], numv);
    }
  }

  flag=vertsgns(lst, cflst, numv,r,&allpos,debugfull);
  if(allpos)// all monomials of the same sign
    (*allflg)=1;
  ifree(lst,r);free((char*)cflst);

  //fprintf(stderr, "got here2b\n");fflush(stdout);fflush(stderr);
  if(debug){
    if(allpos && flag==1)
      fprintf(stderr, "All terms are positive\n");
    else if(allpos && flag==-1)
      fprintf(stderr, "All terms are negative\n");
    else if(flag==1)
      fprintf(stderr, "All vertex terms are positive\n");
    else if(flag==-1)
      fprintf(stderr, "All vertex terms are negative\n");
    else
      fprintf(stderr, "The polynomial has mixed vertices\n");
    fprintf(stderr, "Exiting NewtonPolyVertSgn1\n\n");
  }

  return flag;
}


//
// Numerically test a polynomial for definiteness
// Return 0 if both positive and negative values are found
// Return -1 if only nonpositive values found
// Return 1 if only nonnegative values found
//

int testposneg(ex tmp, double min, double max, int tottries, double tol){
  int i,j;
  char **pvars;
  int numv=polyvars(tmp, &pvars);//get the variable names
  double tmpval1;
  ex tmpval;
  ex v[numv];
  exmap m;
  int flg=2;//initial state (also return state for the zero polynomial)

  //cerr << "entering testposneg" << endl;

  polyvarss(tmp, v, pvars, numv);//get the variables themselves
  freearraydat(pvars, numv);

  //Random seeding should be done externally

  for(j=0;j<tottries;j++){
    for(i=0;i<numv;i++)
      m[v[i]]=min+((double)rand())/((double)RAND_MAX)*(max-min);

    tmpval=evalf(tmp.subs(m));
    //cerr << "value = " << tmpval << endl;
    if (is_a<numeric>(tmpval)) {
      tmpval1=ex_to<numeric>(tmpval).to_double();
      //cout << tmpval1 << endl;
      if(flg==2){
	//nonzero evaluations
	if(tmpval1<-tol)
	  flg=-1;
	else if(tmpval>tol)
	  flg=1;
      }
      else if((flg==1 && tmpval<0.0) || (flg==-1 && tmpval>0.0))//unsigned
	return 0;
    }
    else{
      fprintf(stderr, "ERROR: something went wrong in testposneg: expression is not numeric.\n");
      exit(0);
    }
  }
  return flg;
}

// vecl --- vecmid --- vecr

int *sympt(int *vecmid, int *vecl, int cf, int n, bool pqswitch){
  int *vecr;
  int i;
  if(pqswitch){
    if(!iseven(vecl,n))
      return NULL;
  }
  if(cf<=0) // can't dominate negative with negative or zero
    return NULL;
  vecr=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++){
    vecr[i]=2*vecmid[i]-vecl[i];
    if(vecr[i]<0){ // only positive powers allowed
      free((char *)vecr);
      return NULL;
    }
  }
  return vecr;
}

//    2*sqrt(cfs[s])*sqrt(cfs[s+s1]);
int splt(int a, int b){
  if(a<0 || b<0){
    fprintf(stderr,"ERROR in splt a=%d,b=%d. Exiting.\n",a,b);
    exit(0);
  }
  if(a==0 || b==0)
    return 0;
  int i=sqrt(a);// integer part of sqrt
  int j=sqrt(b);// integer part of sqrt
  return 2*i*j*min(a/(i*i),b/(j*j));
}



// Find all the ways that a negative monomial allt[indx] may be dominated 
// symmetrically by pairs of vectors from the list allt
// assume that allt is ordered

// needs is a vector of pairs which can dominate vec.
// needs[i][0] is the total number of vectors stored (twice number of pairs). 
// needs grows dynamically as more pairs are found

// neededby is a vector of pairs with a similar format to needs
// neededby[i] tells us which vectors vector i (assumed positive) may be needed to dominate
// neededby[i][0] is the total number of vectors stored (twice number of pairs). 
// The first entry in each pair of neededby[i] refers to a negative monomial
// which can be dominated by vector i. The second entry refers to the other
// endpoint needed for this domination.  

// needexcl tells us that the situation is fragile: needexcl[i] is set to 1 if all 
// possible dominating pairs are needed to dominate the (negative) monomial i. 

// The return value of numsplits measure the fragility: 0.1 means that the 
// situation is fragile; larger numbers tell us that we have redundancy in 
// the ways of dominating this negative monomial. 

// numsplits doesn't care whether the monomials come from Re or Im. "honest" splits are those which respect Re and Im, namely we dominate a negative

float numsplits(long indx, int n, int **allt, int *cfs, long r, int **needcount, int **needexcl, int **needs, int **neededby, int debug){
  long s,s1;
  int *vc;
  int j,j1,j2;
  int spl=0;
  int powsplit=0;
  int cf=cfs[indx];
  //fprintf(stderr, "entering numsplits.\n");
  if(debug){
    cout << "---------------------" << endl;
    /* cout << cf << ", "; printvec1(allt[indx],n); */expltomonoprint(cfs[indx],allt[indx],n);
  }
  for(s=0;s<r;s++){
    vc=sympt(allt[indx],allt[s],cfs[s],n,0);
    if(vc && (s1=binisearch(allt+s,r-s,vc,n))>=0 && cfs[s+s1]>0){
      (needs[indx][0])+=2;
      j=(needs[indx][0])+1;
      needs[indx]=(int*) realloc((needs[indx]), sizeof(int)*(j)); 
      
      (neededby[s][0])+=2;
      (neededby[s+s1][0])+=2;
      j1=(neededby[s][0])+1;
      j2=(neededby[s+s1][0])+1;
      neededby[s]=(int*) realloc((neededby[s]), sizeof(int)*(j1));
      neededby[s+s1]=(int*) realloc((neededby[s+s1]), sizeof(int)*(j2));

      powsplit+=splt(cfs[s],cfs[s+s1]);
      needs[indx][j-2]=s;needs[indx][j-1]=s+s1;
      neededby[s][j1-2]=indx;neededby[s][j1-1]=s+s1;
      neededby[s+s1][j2-2]=indx;neededby[s+s1][j2-1]=s;
      ((*needcount)[s])++;
      ((*needcount)[s+s1])++;
      if(debug){
	/* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(cfs[s],allt[s],n);
	/* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(cfs[s+s1],allt[s+s1],n);
      }
      spl++;

    }
    if(vc)
      free((char*)vc);
  }
  if(powsplit==abs(cf)){ //no choice
    (*needexcl)[indx]=1;
  }
  if(debug){
    cout << "powsplit = " << powsplit << endl;
    if(powsplit<abs(cf)){
      cout << "Warning: powers don't add up.\n";
      fprintf(stderr, "failure to split.\n");
      //exit(0);
    }
    else if(powsplit==abs(cf))
      cout << "All splits needed (numsplits).\n";
  }
  if(powsplit<abs(cf)) // basically a failure to split in this way
    return -1.0;
  //fprintf(stderr, "exiting numsplits.\n");
  return sqrt(((float)powsplit)/((float)(abs(cf)))-0.99);
}



int testforpos(int **lst, int *cfs, int polytot, int numv, double min, double max, int tottries, double tol){
  int i,j;
  double m[numv];

  //Random seeding should be done externally
  for(j=0;j<tottries;j++){
    for(i=0;i<numv;i++)
      m[i]=min+((double)rand())/((double)RAND_MAX)*(max-min);

    if(devalpoly(lst, NULL, cfs, 1, polytot, numv, 0, m)<-tol)
      return 0;//possibly negative value found
  }
  return 1;//no negative values found (to numerical precision)
}

//overload

int testforpos(int **lst, double *cfs, int polytot, int numv, double min, double max, int tottries, double tol){
  int i,j;
  double m[numv];

  //Random seeding should be done externally
  for(j=0;j<tottries;j++){
    for(i=0;i<numv;i++)
      m[i]=min+((double)rand())/((double)RAND_MAX)*(max-min);

    if(devalpoly(lst, NULL, cfs, 1, polytot, numv, 0, m)<-tol)
      return 0;//possibly negative value found
  }
  return 1;//no negative values found (to numerical precision)
}

//
// Uses isinconvhull - my implementation using GLPK
//

int hasvertnonsq(int **lst, int numv, long r, int *allsq){
  long s;
  long totsq=0,totnonsq=0;
  int flagsq=1;
  int debug=0;
  long sqs[r], nonsqs[r];
  (*allsq)=0;

  if(debug){printf("entering hasvertnonsq\n");}
  for(s=0;s<r;s++){//list of nonsquares and squares
    if(notsq(lst[s], numv))
      nonsqs[totnonsq++]=s;
    else
      sqs[totsq++]=s;
  }

  if(totnonsq==0){
    (*allsq)=1;
    if(debug){printf("exiting hasvertnonsq (all squares)\n");}
    return 0;//well behaved: all monomials are squares
  }
  else if(totsq==0){
    if(debug){printf("exiting hasvertnonsq (no squares)\n");}
    return 1;//nonsquare vertex
  }  

  //Are any nonsquares outside the convex hull of the squares?
  for(s=0;s<totnonsq;s++){
    if(!isinconvhull(lst, sqs, totsq, numv, lst[nonsqs[s]])){
      flagsq=0;//nonsquare monomial outside convex hull of squares
      if(debug){
	fprintf(stderr, "nonsquare monomial: ");
	printvec(lst[nonsqs[s]],numv);
	fprintf(stderr, "is not in the convex hull of the squares.\n");
      }
      break;
    }
    else if(debug){
      fprintf(stderr, "nonsquare monomial: ");
      printvec(lst[nonsqs[s]],numv);
      fprintf(stderr, "is in the convex hull of the squares.\n");
    }
  }

  if(debug){printf("exiting hasvertnonsq\n");}
  if(!flagsq)
    return 1;//nonsquare vertex
  return 0;//well behaved: all vertices are squares

}


int NewtonPolyVertPosSq(ex tmp, int numv){
  long s,r;
  int **lst;
  double *cflst;
  int polymindeg,polymaxdeg;
  int ret=0;
  int debug=0;

  extolst1(tmp, numv, &lst, &cflst, &r, &polymindeg, &polymaxdeg);
  if(debug){
    for(s=0;s<r;s++){
      fprintf(stderr, "%.4f: ", cflst[s]);printvec(lst[s],numv);
    } 
  }

  ret=NewtonPolyVertPosSq_a(lst, cflst, numv, r);
  ifree(lst,r);free((char*)cflst);
  return ret;

}

//Are all vertices of the Newton polytope positive/negative squares?
//3 means vertices not all squares
//2 means vertices all squares, but of mixed sign
//1 means all vertices are squares and positive
//-1 means all vertices are squares and negative

int NewtonPolyVertPosSq_a(int **lst, int *cflst, int numv, long r){
  int sgnflag=0,sqflag=0;
  int allsq=0;
  int allpos=0;
  int debug=0;
  //fprintf(stderr, "got here\n");fflush(stdout);fflush(stderr);
  if(debug){fprintf(stderr, "Entering NewtonPolyVertPosSq_a.\n");} 
  if(r==0)//Return trivially true for the zero polynomial
    return 1;

  sqflag=hasvertnonsq(lst,numv,r,&allsq);
  if(!sqflag){//all vertices are squares
    sgnflag=vertsgns(lst, cflst, numv,r, &allpos,debug);
    if(allsq){
      if(allpos==1)// all monomials are positive squares: automatically positive and SOS
	sgnflag=10;
      else if (allpos==-1)// all monomials are negative squares: automatically negative and -SOS
	sgnflag=-10;
    }
  }

  //fprintf(stderr, "got here2b %d %d\n", sqflag, sgnflag);fflush(stdout);fflush(stderr);
  if(sqflag)//has nonsquare vertices 
    return 3;
  else{// 10 for all verts +ve squares, -10 for all verts -ve squares, 2 for mixed vertices, 1 for all verts positive, -1 for all verts negative
    return sgnflag;
  }

}


//overloaded to accept double cfs

int NewtonPolyVertPosSq_a(int **lst, double *cflst, int numv, long r){
  int sgnflag=0,sqflag=0;
  int allsq=0;
  int allpos=0;
  int debug=0;
  //fprintf(stderr, "got here\n");fflush(stdout);fflush(stderr);

  if(r==0)//Return trivially true for the zero polynomial
    return 1;

  sqflag=hasvertnonsq(lst,numv,r,&allsq);
  if(!sqflag){//all vertices are squares
    sgnflag=vertsgns(lst, cflst, numv,r, &allpos,debug);
    if(allsq){
      if(allpos==1)// all monomials are positive squares: automatically positive and SOS
	sgnflag=10;
      else if (allpos==-1)// all monomials are negative squares: automatically negative and -SOS
	sgnflag=-10;
    }
  }

  //fprintf(stderr, "got here2b %d %d\n", sqflag, sgnflag);fflush(stdout);fflush(stderr);
  if(sqflag)//has nonsquare vertices 
    return 3;
  else{// 10 for all verts +ve squares, -10 for all verts -ve squares, 2 for mixed vertices, 1 for all verts positive, -1 for all verts negative
    return sgnflag;
  }

}


//
// Uses isinconvhull - my implementation using GLPK
//

int vertsgns(int **lst, int *cflst, int numv, long r, int *allpos, int debug){
  long s;
  long neg=0,pos=0;
  int flagp=1,flagn=1;
  long indp[r], indn[r];
  if(debug){fprintf(stderr, "\n###Entering vertsgns\n");}

  (*allpos)=0;

  if(r==0)
    return 0;

  for(s=0;s<r;s++){
    if(cflst[s]<0)
      indn[neg++]=s;
    else if(cflst[s]>0)
      indp[pos++]=s;
  }

  if(neg==0){
    (*allpos)=1;
    return 1;
  }
  else if(pos==0){
    (*allpos)=-1;
    return -1;
  }

  for(s=0;s<neg;s++){
    if(!isinconvhull(lst, indp, pos, numv, lst[indn[s]])){
      flagp=0;//-ve monomial outside convex hull of +ve ones
      if(debug){
	fprintf(stderr, "negative monomial: ");
	printvec(lst[indn[s]],numv);
	fprintf(stderr, "is not in the convex hull of the positive monomials.\n");
      }
      break;
    }
    else if(debug){
      fprintf(stderr, "negative monomial: ");
      printvec(lst[indn[s]],numv);
      fprintf(stderr, "is in the convex hull of the positive monomials.\n");
    }
  }
  if(!flagp){
    for(s=0;s<pos;s++){
      if(!isinconvhull(lst, indn, neg, numv, lst[indp[s]])){
	flagn=0;//+ve vertex outside convex hull of -ve ones
	if(debug){
	  fprintf(stderr, "positive monomial: ");
	  printvec(lst[indp[s]],numv);
	  fprintf(stderr, "is not in the convex hull of the negative monomials.\n");
	}
	break;
      }
      else if(debug){
	fprintf(stderr, "positive monomial: ");
	printvec(lst[indp[s]],numv);
	fprintf(stderr, "is in the convex hull of the negative monomials.\n");
      }

    }
  }

  if(debug){fprintf(stderr, "Exiting vertsgns.\n");}
  if(flagp)
    return 1;
  else if(flagn)
    return -1;
  return 2;//code for mixed vertices

}

//overloaded to accept double cfs

//
// Uses isinconvhull - my implementation using GLPK
// Doesn't require the polynomial to be homogeneous
//

int vertsgns(int **lst, double *cflst, int numv, long r, int *allpos, int debug){
  long s;
  long neg=0,pos=0;
  int flagp=1,flagn=1;
  long indp[r], indn[r];
  if(debug){fprintf(stderr, "\n###Entering vertsgns\n");}

  (*allpos)=0;

  if(r==0)
    return 0;

  for(s=0;s<r;s++){
    if(cflst[s]<0.0)
      indn[neg++]=s;
    else if(cflst[s]>0.0)
      indp[pos++]=s;
  }

  if(neg==0){
    (*allpos)=1;
    return 1;
  }
  else if(pos==0){
    (*allpos)=-1;
    return -1;
  }

  for(s=0;s<neg;s++){
    if(!isinconvhull(lst, indp, pos, numv, lst[indn[s]])){
      flagp=0;//-ve monomial outside convex hull of +ve ones
      if(debug){
	fprintf(stderr, "negative monomial: ");
	printvec(lst[indn[s]],numv);
	fprintf(stderr, "is not in the convex hull of the positive monomials.\n");
      }
      break;
    }
    else if(debug){
      fprintf(stderr, "negative monomial: ");
      printvec(lst[indn[s]],numv);
      fprintf(stderr, "is in the convex hull of the positive monomials.\n");
    }
  }
  if(!flagp){
    for(s=0;s<pos;s++){
      if(!isinconvhull(lst, indn, neg, numv, lst[indp[s]])){
	flagn=0;//+ve vertex outside convex hull of -ve ones
	if(debug){
	  fprintf(stderr, "positive monomial: ");
	  printvec(lst[indp[s]],numv);
	  fprintf(stderr, "is not in the convex hull of the negative monomials.\n");
	}
	break;
      }
      else if(debug){
	fprintf(stderr, "positive monomial: ");
	printvec(lst[indp[s]],numv);
	fprintf(stderr, "is in the convex hull of the negative monomials.\n");
      }

    }
  }

  if(debug){fprintf(stderr, "Exiting vertsgns.\n");}
  if(flagp)
    return 1;
  else if(flagn)
    return -1;
  return 2;//code for mixed vertices

}

int *midpt(int *vecl, int *vecr, int cfl, int cfr, int n){
  int *vecmid;
  int tmp;
  int i;
  if(cfl<0 || cfr<0) // can't dominate negative with negative
    return NULL;
  vecmid=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++){
    tmp=vecl[i]+vecr[i];
    if(tmp%2!=0){// only even allowed
      free((char *)vecmid);
      return NULL;
    }
    else
      vecmid[i]=tmp/2;
  }
  return vecmid;
}


// A recursive routine. Starts with a positive monomial allt[m0] assumed to be needed (not necessarily exclusively) to dominate some negative monomial, and two lists of other (positive) monomials on the same side as allt[m0], and on the opposite side. 

// Find all the other negative monomials which require mpos to dominate them. Find the other ends of the dominating pairs. Continue recursively until the set grows no further. Check for all midpoints of left pairs and all midpoints of right pairs. 

// "neg" is the list of negative monomials in the group
// "all" is the total number of monomials in the group

// (*failflag)=1 is a fundamental failure (Question: 
// is (*failflag)=-1 ever used?). Can occur if we fail to 
// find midpoints; or the clique gets too large

// IMPORTANT NOTE: what is gathered
// doesn't have to include all the negative midpoints
// So it may not factorise perfectly

// excflag causes gather to fail fundamentally if some midpoint is not found. 
// It causes the recursion to be entered only if absolutely necessary, namely 
// some left or right term gathered is *exclusively* needed by some negative 
// monomial. 

void gather(int **allt, int *cfs, long r, int n, int **neededby,int **needs,int *needexcl,long **lt,long **rt,long **neg,long **mid,long **all,int *numlt,int *numrt,int *numneg,int *nummid,int m0,int startind,int maxsz, int *numsofar, bool onlft, int *failflag, bool tight, bool excflag){
  int j,j1;
  int *vc;
  long s1;
  bool flg1,added;
  //fprintf(stderr, "entering gather...\n");
  //  cout << "numsofar0=" << *numsofar << endl;
  if((*failflag)==1){// failed somewhere
    cout << "failed on entry\n";
    return;
  }
  //  cout << "here1\n";
  if((*numsofar)>maxsz){// too many terms
    (*failflag)=1;
    cout << "entered with too many terms\n";
    return;
  }

  // enter routine; check if a given term has all the necessary midterms. 
  // If so add these midterms and the term itself; if not then return
  // if the term is the first in a left or right list, the check is trivially successful

  if(onlft){ // on the left
    //first check for/add midpoints with existing left terms

    for(j1=0;j1<(*numlt);j1++){ // all must exist
      vc=midpt(allt[m0],allt[(*lt)[j1]],cfs[m0],cfs[(*lt)[j1]],n);
      if(!vc || (s1=binisearch(allt,r,vc,n))<0 || cfs[s1]<=0 || isinlist(s1, (*lt), (*numlt)) || isinlist(s1, (*rt), (*numrt)) || (abs(cfs[s1]) < splt(cfs[m0], cfs[(*lt)[j1]]))){// midpoint not found or negative or already a left or right term or has insufficient weight
    	if(vc)
    	  free((char*)vc);
    	if(tight || excflag){(*failflag)=1;/* cout << "major fail\n"; */}else{(*failflag)=-1;}
    	//cout << "failing to find midterm\n";printvec1(allt[m0],n);printvec1(allt[(*lt)[j1]],n);
    	return;
      }
      if(vc)
    	free((char*)vc);
    }

    // all midterms good

    for(j1=0;j1<(*numlt);j1++){
      vc=midpt(allt[m0],allt[(*lt)[j1]],cfs[m0],cfs[(*lt)[j1]],n);
      s1=binisearch(allt,r,vc,n);
      if(!isinlist(s1, (*mid), (*nummid))){
    	(*numsofar)++;
    	if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}
    	addnewto1Dlarray(mid,(*nummid),s1);(*nummid)++;
    	addnewto1Dlarray(all, (*numsofar)-1,s1);
    	//cout << "adding midterm\n   "; printvec1(allt[s1],n);
      }
      free((char*)vc);
    }

    // now add the actual (positive) monomial to leftterms and all

    (*numsofar)++;
    if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

    addnewto1Dlarray(lt, (*numlt),m0);(*numlt)++;
    addnewto1Dlarray(all, (*numsofar)-1,m0);
    //cout << "adding on left\n   "; printvec1(allt[m0],n);

  }
  else{ // on the right
    // first check for/add midpoints with existing rightterms

    // Question: If the only failure is that it already exists in the list,
    // is this a problem?
    for(j1=0;j1<(*numrt);j1++){
      vc=midpt(allt[m0],allt[(*rt)[j1]],cfs[m0],cfs[(*rt)[j1]],n);
      if(!vc || (s1=binisearch(allt,r,vc,n))<0 || cfs[s1]<=0 || isinlist(s1, (*lt), (*numlt)) || isinlist(s1, (*rt), (*numrt)) || (abs(cfs[s1]) < splt(cfs[m0], cfs[(*rt)[j1]]))){// midpoint not found or negative or already a left or right term or has insufficient weight
    	if(vc)
    	  free((char*)vc);
    	if(tight || excflag){(*failflag)=1;/* cout << "major fail\n"; */}else{(*failflag)=-1;}
	//cout << "failing to find midterm\n";printvec1(allt[m0],n);printvec1(allt[(*rt)[j1]],n);
    	return;
      }
      if(vc)
    	free((char*)vc);
    }

    // all midterms good

    for(j1=0;j1<(*numrt);j1++){
      vc=midpt(allt[m0],allt[(*rt)[j1]],cfs[m0],cfs[(*rt)[j1]],n);
      s1=binisearch(allt,r,vc,n);
      if(!isinlist(s1, (*mid), (*nummid))){
    	(*numsofar)++;
    	if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

    	addnewto1Dlarray(mid,(*nummid),s1);(*nummid)++;
    	addnewto1Dlarray(all, (*numsofar)-1,s1);
    	//cout << "adding midterm\n   "; printvec1(allt[s1],n);
      }
      free((char*)vc);
    }


    // now add the actual (positive) monomial to rightterms and all

    (*numsofar)++;
    if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

    addnewto1Dlarray(rt, (*numrt),m0);(*numrt)++;
    addnewto1Dlarray(all, (*numsofar)-1,m0);
    //cout << "adding on right\n   "; printvec1(allt[m0],n);
  }

  //
  // At this point we have added the positive monomial we started with. 
  // We have checked all the midterms with previous monomials on the same 
  // side, confirmed they exist and added them. We now take the added 
  // monomial as the starting point and look where it is needed: 
  // enter the recursion. If excflag, then we only 
  // enter the recursion provided there is an *exclusive* need for left 
  // or right monomial gathered. In other words, with excflag set, 
  // the default is don't enter the recursion unless we are forced to. 
  // (Keep squares small). 
  //
  // Question: do we not need to do a recursion on the midpoints added too? 
  // Might we be using up positive monomials which are needed for other 
  // dominations? Are we hoping these are automatically taken care of? 

  // For each opposite end-point associated with positive monomial m0...
  for(j=startind;j<neededby[m0][0]+1;j+=2){ 
    added=0;
    flg1=0;if(!excflag ||needexcl[neededby[m0][j]]){flg1=1;}
    if(onlft && !isinlist(neededby[m0][j+1], (*rt), (*numrt)) && flg1){// new right term
      /* if(needs[neededby[m0][j]][0]==2)//has to be tight - negative term allows only one domination */
      /* 	gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],1,maxsz,numsofar,0,failflag,1, excflag); */
      /* else */
      gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],1,maxsz,numsofar,0,failflag,tight, excflag);
      added=1;
      if((*failflag)==1)
	return;
      /* if((*failflag) && needexcl[neededby[m0][j]]){ */
      /* 	/\* printvec(allt[m0],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j]],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j+1]],n); *\/ */
      /* 	/\* printvec(needs[neededby[m0][j]],needs[neededby[m0][j]][0]+1); *\/ */
      /* 	/\* //exit(0); *\/ */
      /* 	(*failflag)==1; return; */
      /* } */
    }
    else if(!onlft && !isinlist(neededby[m0][j+1], (*lt), (*numlt)) && flg1){// add a left term
      /* if(needs[neededby[m0][j]][0]==2)//has to be tight - negative term allows only one domination */
      /* 	gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],maxsz,numsofar,1,failflag,1, excflag); */
      /* else */
      gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],1,maxsz,numsofar,1,failflag,tight, excflag);
      added=1;
      if((*failflag)==1)
	return;
      /* if((*failflag) && needexcl[neededby[m0][j]]){ */
      /* 	/\* printvec(allt[m0],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j]],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j+1]],n); *\/ */
      /* 	/\* printvec(needs[neededby[m0][j]],needs[neededby[m0][j]][0]+1); *\/ */
      /* 	/\* //exit(0); *\/ */
      /* 	(*failflag)==1; return; */
      /* } */
    }

    // should only do this if previous steps were fully successful:
    // add negative midpoints of any endpoints added to arrays "neg" and "all"


    if(added && !(*failflag) && !isinlist(neededby[m0][j], (*neg), (*numneg))){ 
      (*numsofar)++;
      if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

      addnewto1Dlarray(neg, (*numneg),neededby[m0][j]);(*numneg)++;
      addnewto1Dlarray(all, (*numsofar)-1,neededby[m0][j]);

      //cout << "adding negative term\n   "; printvec1(allt[neededby[m0][j]],n);
    }

  }
  //cout << "numsofar5=" << *numsofar << endl;
  //fprintf(stderr, "exiting gather...\n");
  return;

}

void sizesort4(int *cliquesz, long **lt, long **rt, int *numlt, int *numrt, int *absnegmax, int left, int right)
{
  int i,last;

  if (left >= right)
    return;
  iswap(cliquesz,left,(left+right)/2);
  iswap(absnegmax,left,(left+right)/2);
  lpswap(lt,left,(left+right)/2);
  lpswap(rt,left,(left+right)/2);
  iswap(numlt,left,(left+right)/2);
  iswap(numrt,left,(left+right)/2);
  last = left;

  for (i=left+1; i<=right;i++){
    if (cliquesz[i]>cliquesz[left]){// or <
      last++;
      iswap(cliquesz,last,i);
      iswap(absnegmax,last,i);
      lpswap(lt,last,i);
      lpswap(rt,last,i);
      iswap(numlt,last,i);
      iswap(numrt,last,i);
    }
  }

  iswap(cliquesz,left,last);
  iswap(absnegmax,left,last);
  lpswap(lt,left,last);
  lpswap(rt,left,last);
  iswap(numlt,left,last);
  iswap(numrt,left,last);

  sizesort4(cliquesz, lt, rt, numlt,numrt, absnegmax, left, last-1);
  sizesort4(cliquesz, lt, rt, numlt,numrt, absnegmax, last+1,right);
  return;
}

void lkeeponly(long **imat, int n, int m){
  int j;
  if(n>0){
    for(j=m;j<n;j++)
      free ((char *)(imat[j]));
  }
}

// see the routine "gather" to see the meaning of excflag. Roughly, 
// excflag causes cliques to remain small: they only grow if they absolutely 
// have to because some (positive) monomial gathered is exclusively required 
// for domination of some other negative monomial. 

// There seems to be some fragility here. Perhaps, permuting variable 
// names can lead to different cliques (of equal quality) coming first 
// in the natural order, and then divergence in how the algorithm proceeds. 
// Question. Can this actually cause the difference between failure and success?

int allcliques1(int **allt, int *cfs, long r, int n, int **neededby, int **needs, int *needexcl, long ***left, long ***right, int **numleft, int **numright, int **absnegmax, int maxsz, int maxgetcliques, bool tight, bool onlycmplt, bool excflag, int debug){
  long *lt=NULL,*rt=NULL,*neg=NULL,*mid=NULL,*all=NULL;
  int numlt=0,numrt=0,numneg=0,nummid=0;
  int *cliquequal=NULL;// "quality" of a clique
  int numsofar=0;
  int failflag=0;
  int *vc;
  long m0,s1;
  int k,k1;
  long *outvec;
  int numcliques=0;
  bool flg,cmplt;
  long **cliquearray=NULL;
  int *cliquesz=NULL;
  int tmpabsneg=0;
  int st,mx;
  int maxcliques=100, cliqueqtmp;

  // for each term in the full expansion of (Re^2 + Im^2)...

  for(m0=0;m0<r;m0++){ 
    if(excflag) // in this case only start gathering once on each monomial
      mx=2;
    else
      mx=neededby[m0][0];

    // m0 is a +ve monomial potentially needed by some -ve monomial:
    // use it to seed "gather" to form a clique
    if(cfs[m0]>0 && neededby[m0][0]>0){ 
      //cout << "checking: (mx="<<mx<<"): "; printvec1(allt[m0],n);
      for(st=1;st<mx;st+=2){
	numlt=0;numrt=0;numneg=0;nummid=0;
	numsofar=0;failflag=0;
	//fprintf(stderr, "entering gather.\n");
	gather(allt,cfs,r,n,neededby,needs,needexcl,&lt,&rt,&neg,&mid,&all,&numlt,&numrt,&numneg,&nummid,m0,st,maxsz,&numsofar,1,&failflag,tight,excflag);
	//fprintf(stderr, "exiting gather.\n");
	//cout << "gathered\n" << failflag <<" " << numsofar << endl;
	if(numsofar>=3 && numneg && (failflag!=1)){
	  //cout << "found clique\n";
	  //store the clique in outvec
	  unionlvec(lt,numlt,rt,numrt, &outvec);
	  flg=0;
	  for(k=0;k<numcliques;k++){
	    if(unordlistareeq(outvec,numlt+numrt,cliquearray[k],cliquesz[k])){// only check left and right halves
	      flg=1;
	      break;
	    }
	  }
	  if(!flg){// not already in array of cliques; add + build output info

	    if(onlycmplt){//only accept complete cliques
	      //is the clique complete?
	      cmplt=1;k=0;k1=0;
	      while(cmplt && k<numlt){
		while(cmplt && k1<numrt){
		  vc=midpt(allt[lt[k]],allt[rt[k1]],cfs[lt[k]],cfs[rt[k1]],n);
		  s1=binisearch(allt,r,vc,n);
		  if(s1<0 || !(isinlist(s1,neg,numneg)))
		    cmplt=0;
		  k1++;
		  if(vc){free((char*)vc);}
		}
		k++;
	      }
	    }

	    if(!onlycmplt || cmplt){
	      if(debug){
	      //fprintf(stderr, "found clique number %d (size: %d, neg: %d)\n", numcliques+1, numsofar, numneg);
		cout << "exiting with l,r,mid,n, tot " << numlt << ", " << numrt << ", " <<  nummid << ", " <<  numneg << ", " <<  numsofar << endl;
		cout << "found a clique of size " << numsofar << endl;
		for(k=0;k<numlt;k++){
		  cout << cfs[lt[k]] << ", "; printvec1(allt[lt[k]],n);
		}
		for(k=0;k<numrt;k++){
		  cout << cfs[rt[k]] << ", "; printvec1(allt[rt[k]],n);
		}
		for(k=0;k<numneg;k++){
		  cout << cfs[neg[k]] << ", "; printvec1(allt[neg[k]],n);
		}
		for(k=0;k<nummid;k++){
		  cout << cfs[mid[k]] << ", "; printvec1(allt[mid[k]],n);
		}
	      }
	      //	      fprintf(stderr, "found a clique\n");
	      addnewto1Darray(&cliquesz,numcliques,numlt+numrt);
	      addnewtoarray(&cliquearray,numcliques,outvec,numlt+numrt);
	      addnewtoarray(left,numcliques,lt,numlt);
	      addnewto1Darray(numleft,numcliques,numlt);
	      addnewtoarray(right,numcliques,rt,numrt);
	      addnewto1Darray(numright,numcliques,numrt);
	      cliqueqtmp=(10*numneg)/(numlt+numrt) + numneg/5;
	      //cliqueqtmp=(10*numneg)/(numlt+numrt) - numneg/10;
	      if(debug){cout << "clique quality (allcliques1): " << cliqueqtmp << endl;}
	      addnewto1Darray(&cliquequal,numcliques,cliqueqtmp);

	      for(k=0;k<numneg;k++)
		tmpabsneg=max(tmpabsneg,abs(cfs[neg[k]]));
	      addnewto1Darray(absnegmax,numcliques,tmpabsneg);

	      numcliques++;
	      if(maxgetcliques && numcliques>=maxgetcliques){
	      	if(debug){fprintf(stderr, "got %d cliques\n", maxgetcliques);}
		sizesort4(cliquequal, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
		//sizesort4(cliquesz, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
		if(numcliques>maxcliques){// free the ones which won't be returned
		  lkeeponly(*left, numcliques, maxcliques);
		  lkeeponly(*right, numcliques, maxcliques);
		}

	      	free((char*)outvec);
	      	if(lt){free((char*)lt);lt=NULL;}if(rt){free((char*)rt);rt=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(all){free((char*)all);all=NULL;}
	      	if(cliquearray){lfree_i(cliquearray, numcliques);}
	      	if(cliquesz){free((char *)cliquesz);}
	      	if(cliquequal){free((char *)cliquequal);}

		if(numcliques>maxcliques)
		  return maxcliques;
		return numcliques;
	      }
	    }
	  }
	  free((char*)outvec);
	}
	if(lt){free((char*)lt);lt=NULL;}if(rt){free((char*)rt);rt=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(all){free((char*)all);all=NULL;}
      }
    }
  }
  if(debug){fprintf(stderr, "got %d cliques\n", numcliques);
  cout << "total cliques found (allcliques1) = "<< numcliques << endl;}

  // Sort cliques by quality/size before returning?

  sizesort4(cliquequal, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
  //sizesort4(cliquesz, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
  
  if(numcliques>maxcliques){// free the ones which won't be returned
    lkeeponly(*left, numcliques, maxcliques);
    lkeeponly(*right, numcliques, maxcliques);
  }
  if(cliquearray){lfree_i(cliquearray, numcliques);}
  if(cliquesz){free((char *)cliquesz);}
  if(cliquequal){free((char *)cliquequal);}

  if(numcliques>maxcliques)
    return maxcliques;
  return numcliques;
}

// make a perfect square from a list
// takes the terms on the left and those on the right
// computes the midpoints
// assumes that consistency checks have been done earlier

ex makesquare(int **l, int *cfs, long *leftinds, int numleft, long *rightinds, int numright, int n, int maxnegcf, bool simp){
  ex extot=0,extot1=0;
  ex ex1=1;
  int i,j,k,k2;
  int minfac=1000;
  int *vc;
  int *kleft=(int *)malloc((size_t) (numleft*sizeof(int)));
  int *kright=(int *)malloc((size_t) (numright*sizeof(int)));

  if(maxnegcf==2){// simplest case
    minfac=1;
    inittoone(kleft,numleft);
    inittoone(kright,numright);
  }
  else if(simp){
    minfac=1;
    for(i=0;i<numleft;i++)
      kleft[i]=sqrt(cfs[leftinds[i]]);
    for(i=0;i<numright;i++)
      kright[i]=sqrt(cfs[rightinds[i]]);//integer square root
  }
  else{
    // get the smallest common factor
    for(i=0;i<numleft;i++){
      k=sqrt(cfs[leftinds[i]]);//integer square root
      kleft[i]=k;
      minfac=min(minfac,cfs[leftinds[i]]/(k*k));//integer division
    }
    for(i=0;i<numright;i++){
      k=sqrt(cfs[rightinds[i]]);//integer square root
      kright[i]=k;
      minfac=min(minfac,cfs[rightinds[i]]/(k*k));//integer division
    }
  }
  extot1=facsq(l, leftinds, kleft, numleft, rightinds, kright, numright, n, minfac);
  if(extot1!=-1){
    free((char*)kleft);free((char*)kright);
    return extot1;
  }


  for(i=0;i<numleft;i++){//each monomial on left
    ex1=expltomono(minfac*kleft[i]*kleft[i], l[leftinds[i]],n);
    extot+=ex1;
  }
  for(i=0;i<numright;i++){//each monomial on right
    ex1=expltomono(minfac*kright[i]*kright[i], l[rightinds[i]],n);
    extot+=ex1;
  }

  // left midpoints (positive)

  for(i=0;i<numleft;i++){
    for(j=i+1;j<numleft;j++){
      vc=midpt(l[leftinds[i]],l[leftinds[j]],cfs[leftinds[i]], cfs[leftinds[j]],n);
      k2=2*kleft[i]*kleft[j]*minfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot+=ex1;
    }
  }

  // right midpoints (positive)

  for(i=0;i<numright;i++){
    for(j=i+1;j<numright;j++){
      vc=midpt(l[rightinds[i]],l[rightinds[j]],cfs[rightinds[i]], cfs[rightinds[j]],n);
      k2=2*kright[i]*kright[j]*minfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot+=ex1;
    }
  }

 // mixed midpoints (negative)

  for(i=0;i<numleft;i++){
    for(j=0;j<numright;j++){
      vc=midpt(l[leftinds[i]],l[rightinds[j]],cfs[leftinds[i]], cfs[rightinds[j]],n);
      k2=2*kleft[i]*kright[j]*minfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot-=ex1;
    }
  }


  /* if(extot-expand(extot1)!=0){ */
  /*   fprintf(stderr, "serious problem!\n"); */
  /*   cout << extot << endl; */
  /*   cout << extot1 << endl; */
  /*   exit(0); */
  /* } */
  free((char*)kleft);free((char*)kright);
  /* if(extot1==-1) */
  /*   return extot; */
  return extot;
}



// Evaluate a polynomial at vals
// Fourth argument: if true, then the cfs are indexed by inds, otherwise, they are independent

double devalpoly(int **lst, int *inds, int *cflst, int cfsindexed, int tot, int numv, int offset, double *vals){
  double ret=0.0, ret1;
  int i,i1,j;
  /* fprintf(stderr, "Evaluating polynomial at:\n"); */
  /* for(i=0;i<numv;i++) */
  /*   fprintf(stderr, "%.2f ", vals[i]); */
  /* fprintf(stderr, "\n"); */
  for(i=0;i<tot;i++){
    if(inds)
      i1=inds[i];
    else
      i1=i;
    if(cfsindexed)
      ret1=cflst[i1];
    else
      ret1=cflst[i];
    for(j=0;j<numv;j++){
      if(lst[i1][j+offset])
	ret1*=pow(vals[j], (double)(lst[i1][j+offset]));
    }
    ret+=ret1;
  }
  /* fprintf(stderr, "Returning %.4f\n", ret); */
  return ret;
}

//overload

double devalpoly(int **lst, int *inds, double *cflst, int cfsindexed, int tot, int numv, int offset, double *vals){
  double ret=0.0, ret1;
  int i,i1,j;
  /* fprintf(stderr, "Evaluating polynomial at:\n"); */
  /* for(i=0;i<numv;i++) */
  /*   fprintf(stderr, "%.2f ", vals[i]); */
  /* fprintf(stderr, "\n"); */
  for(i=0;i<tot;i++){
    if(inds)
      i1=inds[i];
    else
      i1=i;
    if(cfsindexed)
      ret1=cflst[i1];
    else
      ret1=cflst[i];
    for(j=0;j<numv;j++){
      if(lst[i1][j+offset])
	ret1*=pow(vals[j], (double)(lst[i1][j+offset]));
    }
    ret+=ret1;
  }
  /* fprintf(stderr, "Returning %.4f\n", ret); */
  return ret;
}

int notsq(int *lst, int numv){
  int i;
  for(i=0;i<numv;i++){
    if(lst[i]%2==1)
      return 1;
  }
  return 0;
}


ex facsq(int **l, long *leftinds, int *kleft, int numleft, long *rightinds, int *kright, int numright, int numv, int minfac){
  int *exfaci;
  int i,*vc,*vc1;
  ex exout=0,exfac;
  exfaci=intersect(l, leftinds, numleft, rightinds, numright, numv);
  for(i=0;i<numleft;i++){
    vc=subtract(l[leftinds[i]],exfaci,numv);
    if(!iseven(vc,numv)){
      exout=-1;free((char*)vc);return -1;
    }
    vc1=halve(vc,numv);
    /* cout << "hh: " << expltomono(kleft[i],vc1,numv) << endl; */
    exout+=expltomono(kleft[i],vc1,numv);
    free((char*)vc1);free((char*)vc);
  }

  for(i=0;i<numright;i++){
    vc=subtract(l[rightinds[i]],exfaci,numv);
    if(!iseven(vc,numv)){
      exout=-1;free((char*)vc);return -1;
    }
    vc1=halve(vc,numv);
    /* cout << "hh: " << expltomono(kright[i],vc1,numv) << endl; */
    exout-=expltomono(kright[i],vc1,numv);
    free((char*)vc1);free((char*)vc);
  }
  exfac=expltomono(1,exfaci,numv);
  free((char*)exfaci);
  return minfac*exfac*pow(exout,2);
}


int *intersect(int **l, long *leftinds, int numleft, long *rightinds, int numright, int numv){// intersection of a set of monomials
  int *isect=(int *)malloc((size_t) (numv*sizeof(int)));
  int i,j;

  for(i=0;i<numv;i++)
    isect[i]=100;

  for(i=0;i<numv;i++){
    for(j=0;j<numleft;j++){
      isect[i]=min(isect[i],l[leftinds[j]][i]);
      if(!isect[i])
	break;
    }
    if(isect[i]){
      for(j=0;j<numright;j++){
	isect[i]=min(isect[i],l[rightinds[j]][i]);
	if(!isect[i])
	  break;
      }
    }
  }
  return isect;
}


// Check if a matrix is trivially a P matrix or a P0 matrix
// Return 1 means P matrix, return 2 means P0
// Just uses GiNaC flags on each minor (i.e., presumably all coeffs +ve)

int isPmatrix(matrix J, int n, int debug){
  int xc[n];
  int k, flag, ret=0, pflag=1;
  ex detex;

  if(debug){cerr << "Checking matrix" << endl << J;}

  for(k=1;k<=n;k++){
    firstcomb(xc, n, k);flag=1;
    
    while(flag==1){
      detex = symbdetsubmat(J, n, n, xc, xc, k);
      /* if(debug){printvec1(xc,k);} */
      if(detex!=0 && !(detex.info(info_flags::positive))){
	if(debug){
	  cerr << detex << endl;
	  printsubmat(J,xc,xc,k,k);
	  cerr << "A possibly negative principal minor found.\n";
	}
	return 0;
      }
      else if(detex==0)// may still be P0
	pflag=0;
      flag=nextcomb(xc, n, k);
    }
  }
  if(!pflag){
    if(debug){cerr << "The matrix appears to be a P0 matrix.\n";}
    ret=2;
  }
  else{
    if(debug){cerr << "The matrix appears to be a P matrix.\n";}
    ret=1;
  }
  return ret;
  
}

//Integer matrix is strictly P?
int isPmatrixstrict(int **J, int n){
  int xc[n];
  int k,flag;
  int detex;

  for(k=1;k<=n;k++){
    firstcomb(xc, n, k);flag=1;
    
    while(flag==1){
      detex = detsubmat(J, n, n, xc, xc, k);
      if(detex<=0)
	return 0;
      flag=nextcomb(xc, n, k);
    }
  }

  return 1;
}

//Integer matrix is positive definite? Sylvester's criterion
//Assumes symmetric matrix
//https://en.wikipedia.org/wiki/Sylvester%27s_criterion
int isposdef(int **J, int n){
  int xc[n];
  int k;
  int detex;

  for(k=1;k<=n;k++){
    firstcomb(xc, n, k);
    detex = detsubmat(J, n, n, xc, xc, k);
    if(detex<=0)
      return 0;
  }
  return 1;
}

//overloading (assumes symmetric matrix)
int isposdef(matrix J, int n, int maxpppdeg, int *pppdegused, int debug){
  int xc[n];
  int k, deg;
  int dcfs=0;//assuming integer coeffs (or use polyhasdcfs if not sure)
  ex detex, detex_simp;
  int numv, numvinit,allflag;
  (*pppdegused)=-1;

  for(k=1;k<=n;k++){
    firstcomb(xc, n, k);
    detex = symbdetsubmat(J, n, n, xc, xc, k);
    detex_simp=polyhomsimp(expand(detex), &numvinit, &numv, &deg, 1, dcfs, debug);
    if(ispospoly(detex_simp, numv, &allflag, 0, dcfs, maxpppdeg, pppdegused, debug)<0){
      if(debug){
	cerr << detex_simp << endl;
	printsubmat(J,xc,xc,k,k);
	cerr << "A possibly negative leading principal minor found. dim = " << k << endl;
      }
      return 0;
    }
  }
  return 1;
}

//
// J is an n X n matrix with entries which are polynomials over Z 
// Is J a P0 matrix on the positive orthant?
// The flag skiptop is to ignore the top dimension
// For example if we already know that J is the square of another matrix
// then the determinant is automatically nonnegative
int isP0matorth(matrix J, int n, bool skiptop, int maxpppdeg, int *pppdegused, int debug){
  int xc[n];
  int k, flag, deg;
  ex detex, detex_simp;
  int numv, numvinit,allflag;
  int dcfs=0;//assuming integer coeffs (or use polyhasdcfs if not sure)
  int nmax;
  (*pppdegused)=-1;
  if(skiptop)
    nmax=n-1;
  else
    nmax=n;

  if(debug){fprintf(stderr, "\n###Entering isP0matorth. Examining matrix:\n");printmat(J, n, n);}

  for(k=1;k<=nmax;k++){
    firstcomb(xc, n, k);flag=1;
    
    while(flag==1){
      detex = symbdetsubmat(J, n, n, xc, xc, k);
      //cerr<<detex<< endl;
      /* dcfs=polyhasdcfs(detex); *///If we're not sure
      detex_simp=polyhomsimp(expand(detex), &numvinit, &numv, &deg, 1, dcfs, debug);
      //cerr<<detex_simp<< endl;
      /* if(debug){printvec1(xc,k);} */
      //penultimate argument means that the cfs are assumed to be integers
      //previous argument is maximum degree of LHS multiplier
      if(ispospoly(detex_simp, numv, &allflag, 0, dcfs, maxpppdeg, pppdegused, debug)<0){
	if(debug){
	  cerr << detex_simp << endl;
	  printsubmat(J,xc,xc,k,k);
	  cerr << "Exiting isP0matorth. A possibly negative principal minor found. dim = " << k << endl;
	}
	return 0;
      }
      flag=nextcomb(xc, n, k);
    }
  }
  if(debug){cerr << "Exiting isP0matorth. The matrix appears to be a P0 matrix.\n";}
  return 1;  
}



//Overloading: integer matrix
int isP0matorth(int **J, int n, bool skiptop, int debug){
  int xc[n];
  int k, flag;
  int nmax;

  if(skiptop){nmax=n-1;}else{nmax=n;}

  for(k=1;k<=nmax;k++){
    firstcomb(xc, n, k);flag=1;    
    while(flag==1){
      if(detsubmat(J, n, n, xc, xc, k)<0){
	if(debug){cerr << "A negative principal minor found.\n";}
	return 0;
      }
      flag=nextcomb(xc, n, k);
    }
  }

  if(debug){cerr << "The matrix appears to be a P0 matrix.\n";}
  return 1;  
}

// overloading (no skipping the top dimension)
int isP0matorth(matrix J, int n, int maxpppdeg, int *pppdegused, int debug){
  return isP0matorth(J, n, 0, maxpppdeg, pppdegused, debug);
}

// overloading (integer matrix and no skipping the top dimension)
int isP0matorth(int **J, int n, int debug){
  return isP0matorth(J, n, 0, debug);
}

//Sum of principal minors positive
int isQmatrix(matrix J, int n){
  int xc[n];
  int k;
  int flag;
  ex detex;
  for(k=1;k<=n;k++){
    firstcomb(xc, n, k);flag=1;
    detex=0;
    while(flag==1){ // sum all minors of particular size
      //    printvec1(xc,k);
      detex+=symbdetsubmat(J, n, n, xc, xc, k);
      printvec1(xc,k);
      flag=nextcomb(xc, n, k);
    }
    if(detex!=0 && !(detex.info(info_flags::positive))){
      cout << detex << endl;
      cout << "An apparently negative coefficient found.\n";
      return 0;
    }
  }
  return 1;
}


//Is the square matrix sign-symmetric?
int signsym(int **J, int n){
  int xc[n], yc[n];
  int detex,detex1;
  int j,k,t,flag, flag1;

  for(k=1;k<=n-1;k++){//no need for top dimension
    firstcomb(xc, n, k);flag=1;
    while(flag==1){
      for(j=0;j<k;j++)
	yc[j]=xc[j];
      flag1=1;t=0;
      while(flag1==1){
	if(t){
	  detex = detsubmat(J, n, n, xc, yc, k);
	  detex1 = detsubmat(J, n, n, yc, xc, k);
	  if(detex*detex1<0){//not sign symmetric
	    return 0;
	  }
	}
	t++;
	flag1=nextcomb(yc, n, k);
      }
      flag=nextcomb(xc, n, k);
    }
  }
  return 1;
}

//overloading: symbolic matrix
int signsym(matrix J, int n, int maxpppdeg, int debug){
  int xc[n], yc[n];
  ex detex,detex1,detexprod,detex_simp;
  int j,k,t,flag, flag1, deg;
  int numv, numvinit,allflag;
  int dcfs=0;//assuming integer coeffs (or use polyhasdcfs if not sure)
  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering signsym.\n");}

  for(k=1;k<=n-1;k++){//no need for top dimension
    firstcomb(xc, n, k);flag=1;
    while(flag==1){
      for(j=0;j<k;j++)
	yc[j]=xc[j];
      flag1=1;t=0;
      while(flag1==1){
	if(t){
	  detex = symbdetsubmat(J, n, n, xc, yc, k);
	  detex1 = symbdetsubmat(J, n, n, yc, xc, k);
	  detexprod=expand(detex*detex1);
	  detex_simp=polyhomsimp(expand(detexprod), &numvinit, &numv, &deg, 1, dcfs, debugfull);
	  if(ispospoly(detex_simp, numv, &allflag, 0, dcfs, maxpppdeg, &pppdegused, debugfull)<0){
	    if(debug){fprintf(stderr, "Exiting signsym: the matrix is not sign symmetric.\n");}
	    return 0;
	  }
	}
	t++;
	flag1=nextcomb(yc, n, k);
      }
      flag=nextcomb(xc, n, k);
    }
  }

  if(debug){fprintf(stderr, "Exiting signsym: the matrix is sign symmetric.\n");}
  return 1;
}

//The sign of the determinant of J. Return codes are those of ispospoly
int DetSgn(matrix J, int n, int numv, int maxpppdeg, int debug){
  ex detex=det(J,n);
  int allflg,ret;
  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;
  if(debug){fprintf(stderr, "\n###Entering DetSgn.\nJ = ");printmat(J,n,n);cerr << "det: " << detex << endl;}
  ret=ispospoly(detex, numv, &allflg, 0, 0, maxpppdeg, &pppdegused, debugfull);
  if(debug){fprintf(stderr, "Exiting DetSgn.\n");}
  return ret;
}

// The coeff of lambda^{n-k} in the charpoly of J
// Usually where we want k to be the rank of Gamma, the first factor
int DetSgnk(matrix J, int k, int n, int numv, int maxpppdeg, int debug){
  ex detex=getminorsum0(J, n, k);
  int allflg;
  int pppdegused;
  return ispospoly(detex, numv, &allflg, 0, 0, maxpppdeg, &pppdegused, debug);
}



// check if the second additive compound matrix has positive determinant
// Returns 2 for positive, 1 for nonnegative, 0 for identically zero, 
// -1 for nonpositive, -2 for negative, -3 for mixed sign
// and -4 for unable to determine (fundamental failure)
// why the numv+1?
int AdComp2DetPos(matrix J, int n, int numv, int maxpppdeg, int *pppdegused, int debug){
  ex detex=AdComp2Det(J, n, debug);
  int allflg;
  (*pppdegused)=-1;
  return ispospoly(detex, numv, &allflg, 0, 0, maxpppdeg, pppdegused, debug);
}


// check if the second additive compound matrix is a P0 matrix

int AdComp2isP0(matrix J, int n, int maxpppdeg, int debug){
  matrix J2 = AdComp2(J,n);
  int pppdegused;
  if(debug){fprintf(stderr, "Additive compound:\n");
    printmat(J2,n*(n-1)/2,n*(n-1)/2);}
  return isP0matorth(J2, n*(n-1)/2, maxpppdeg, &pppdegused, debug);
}

// for use on a pre-processed polynomial
// i.e., (expanded, canonical variables and homogeneous)
// dcfs means coefficients are doubles
// deg is the degree
// The return values depend on the filter: positivity checks shadow
// ispospoly
// The filter all produces verbose output, even if debug is set to 0
int polytest(ex tmp, int numv, int deg, int dcfs, const char *filter, int maxpppdeg, int debug){
  ex tmp1, tmp2;
  int allflg=0, powval;
  int ret;
  char mesg[100];
  int pppdegused;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering polytest.\n");}

  cerr << "Examining the polynomial:\n\t" << tmp << endl;
  if(!strcmp(filter,"factor") || !strcmp(filter,"all")){
    if(!dcfs){
      tmp1=factor(tmp);
      if(tmp1!=tmp){
	cerr << "\nTrying to factorise (GiNaC): " << tmp << " = " << tmp1 << endl;
	if(isapower(tmp1)){
	  tmp2=getroot(tmp1,&powval);
	  cerr << "This polynomial is a power (" << powval << ") of another: " << tmp2 << endl;
	}
	if(strcmp(filter,"all"))
	  return 1;//successfully factorised
      }
      else{
	cerr << "Searched for a factorisation (GiNaC). None found.\n";
	if(strcmp(filter,"all"))
	  return 0;
      }
    }
    else{
      fprintf(stderr, "Will not try to factor a polynomial whose coefficients are not integers.\n");
      if(strcmp(filter,"all"))
	return 0;
    }
  }

  //Newton polytope
  if(!strcmp(filter,"newton") || !strcmp(filter,"all")){
    fprintf(stderr, "\n**Examining the Newton polytope.\n");
    //don't debug unless desired
    ret=NewtonPolyVertSgn1(tmp, numv, &allflg, debugfull);
    fprintf(stderr,"The result of NewtonPolyVertSgn1: %d\n", ret);
    fprintf(stderr, "This means that ");
    if(ret==0){
      fprintf(stderr, "this polynomial is zero.\n");
      if(strcmp(filter,"all"))
	return 0;
    }
    else if(ret==1 && allflg){
      fprintf(stderr, "all terms are positive.\n");
      if(strcmp(filter,"all"))
	return 2;
    }
    else if(ret==1 && !allflg){
      fprintf(stderr, "all vertex terms are positive (but not all terms are positive).\n");
      if(strcmp(filter,"all"))
	return 1;
    }
    else if(ret==-1 && allflg){
      fprintf(stderr, "all terms are negative.\n");
      if(strcmp(filter,"all"))
	return -2;
    }
    else if(ret==-1 && !allflg){
      fprintf(stderr, "all vertex terms are negative (but not all terms are negative).\n");
      if(strcmp(filter,"all"))
	return -1;
    }
    
    else if(ret==2)
      fprintf(stderr, "there are both positive and negative vertex monomials.\n");
    if(strcmp(filter,"all"))
      return -3;
  }

  //Heuristic attempt at SOS decomposition
  if(!strcmp(filter,"heuristic") || !strcmp(filter,"all")){
    fprintf(stderr, "\n**Trying a heuristic to find an SOS+(positive terms) decomposition\n");
    if(!dcfs){
      ret=heuristic_squares(tmp,numv,debugfull);
      fprintf(stderr,"The result of heuristic_squares (homogenised polynomial): %d\nThis means that on the positive orthant ", ret);
      if(ret==0){
	fprintf(stderr, "this polynomial is zero.\n");
	if(strcmp(filter,"all"))
	  return 0;
      }
      else if(ret==1){
	fprintf(stderr, "this polynomial is nonnegative.\n");
	if(strcmp(filter,"all"))
	  return 1;
      }
      else if(ret==-1){
	fprintf(stderr, "this polynomial is nonpositive.\n");
	if(strcmp(filter,"all"))
	  return -1;
      }
      else if(ret==2){
	fprintf(stderr, "this polynomial is (strictly) positive.\n");
	if(strcmp(filter,"all"))
	  return 2;
      }
      else if(ret==-2){
	fprintf(stderr, "this polynomial is (strictly) negative.\n");
	if(strcmp(filter,"all"))
	  return -2;
      }
      else if(ret==-3){
	fprintf(stderr, "this polynomial can take all signs.\n");
	if(strcmp(filter,"all"))
	  return -3;
      }
      else{
	fprintf(stderr, "the heuristic algorithm failed to find a simple SOS+(positive terms) decomposition.\n");
	if(strcmp(filter,"all"))
	  return -4;
      }
 
    }
    else{
      fprintf(stderr, "Will not try \"heuristic_squares\" on a polynomial whose coefficients are not integers.\n");
      return 0;
    }
  }


  //
  // Positive on the positive orthant?
  //
  if(!strcmp(filter,"posorth") || !strcmp(filter,"all")){
    fprintf(stderr, "\n**Examining the sign of the polynomial on the positive orthant using \"ispospoly\".\n");
    ret=ispospoly(tmp,numv,&allflg,0,dcfs,maxpppdeg,&pppdegused,debugfull);
    if(ret==2 || ret==-2)
      sprintf(mesg, " (with SDP degree %d)", pppdegused);
    else
      sprintf(mesg, " (tried SDP up to degree %d)", maxpppdeg);
    fprintf(stderr,"The result of ispospoly%s: %d\nThis means that on the positive orthant ", maxpppdeg<0?" (not using SDP)":mesg, ret);
    if(ret==0){
      fprintf(stderr, "this polynomial is zero.\n");
      if(strcmp(filter,"all"))
	return 0;
    }
    else if(ret==1){
      fprintf(stderr, "this polynomial is nonnegative.\n");
      if(strcmp(filter,"all"))
	return 1;
    }
    else if(ret==-1){
      fprintf(stderr, "this polynomial is nonpositive.\n");
      if(strcmp(filter,"all"))
	return -1;
    }
    else if(ret==2){
      fprintf(stderr, "this polynomial is (strictly) positive.\n");
      if(strcmp(filter,"all"))
	return 2;
    }
    else if(ret==-2){
      fprintf(stderr, "this polynomial is (strictly) negative.\n");
      if(strcmp(filter,"all"))
	return -2;
    }
    else if(ret==-3){
      fprintf(stderr, "this polynomial is indefinite.\n");
      if(strcmp(filter,"all"))
	return -3;
    }
    else{
      fprintf(stderr, "the algorithm failed to determine the sign of this polynomial.\n");
      if(strcmp(filter,"all"))
	return -4;
    }

  }

  //
  // Positive everywhere?
  //
  if(!strcmp(filter,"posall") || !strcmp(filter,"all")){
    fprintf(stderr, "\n**Examining the sign of this polynomial on all of R^n\\{0} using \"isfullpospoly\".\n");
    ret=isfullpospoly(tmp,numv,0,maxpppdeg,&pppdegused,0,debugfull);
    sprintf(mesg, " (trying SDP up to degree %d)", maxpppdeg);
    fprintf(stderr,"The result of isfullpospoly%s: %d\nThis means that on R^n\\{0} ", maxpppdeg<0?" (not using SDP)":mesg, ret);
    if(ret==0){
      fprintf(stderr, "this polynomial is zero.\n");
      if(strcmp(filter,"all"))
	return 0;
    }
    else if(ret==1){
      fprintf(stderr, "this polynomial is nonnegative.\n");
      if(strcmp(filter,"all"))
	return 1;
    }
    else if(ret==-1){
      fprintf(stderr, "this polynomial is nonpositive.\n");
      if(strcmp(filter,"all"))
	return -1;
    }
    else if(ret==2){
      fprintf(stderr, "this polynomial is (strictly) positive.\n");
      if(strcmp(filter,"all"))
	return 2;
    }
    else if(ret==-2){
      fprintf(stderr, "this polynomial is (strictly) negative.\n");
      if(strcmp(filter,"all"))
	return -2;
    }
    else if(ret==-3){
      fprintf(stderr, "this polynomial is indefinite.\n");
      if(strcmp(filter,"all"))
	return -3;
    }
    else if(ret==-4){
      fprintf(stderr, "the algorithm failed to determine the sign of this polynomial.\n");
      if(strcmp(filter,"all"))
	return -4;
    }
  }

  return 0;

}


int polytest(const char *file, const char *filter, int maxpppdeg, int debug){
  ex tmphom;
  int numvhom=0,numvin, deg, dcfs=0;

  //make homogeneous
  tmphom=polyfromfile(file, &numvin, &numvhom, &deg, 1, &dcfs, debug);

  return polytest(expand(tmphom), numvhom, deg, dcfs, filter, maxpppdeg, debug);

}
