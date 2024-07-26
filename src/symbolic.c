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
#include <limits.h>

const possymbol & get_possymbol(const string & s)
{
  static map<string, possymbol> directory;
  map<string, possymbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, possymbol(s))).first->second;
}

const symbol & get_symbol(const string & s)
{
  static map<string, symbol> directory;
  map<string, symbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, symbol(s))).first->second;
}


//
// Use after trying to factorise to check if something is a perfect power
//

int isapower(ex e){
  if(is_a<power>(e)){
    if(is_a<numeric>(e.op(1))){
      if((ex_to<numeric>(e.op(1))).is_integer())
	return ex_to<numeric>(e.op(1)).to_int();
    }
  }
  return 0;
}

// If a polynomial is a power of another, then replace with the root
// assumes the original polynomial is in expanded form

ex getroot(ex e, int *powval){
  ex tmp=factor(e);
  ex tmp1;
  (*powval)=0;//not a root
  if(is_a<power>(tmp)){
    if(is_a<numeric>(tmp.op(1)) && (ex_to<numeric>(tmp.op(1))).is_integer()){
      (*powval)=ex_to<numeric>(tmp.op(1)).to_int();
      tmp1=expand(tmp.op(0));
      return tmp1;
    }
  }

  return e;
}

//
// monomial to a list
//

lst monotolist(ex e){
  lst tmp;
  int j;
  //int coe=-1;
  //int r;

 if(is_a<add>(e)){
    cout << "failure in monotolist: not expecting a sum:\n";
    cout << e << endl;
    exit(0);
  }

  if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  if(is_a<symbol>(e.op(k).op(0)))
	    tmp.append(e.op(k).op(0));
	  else{
	    cout << "failure 1 in monotolist: expecting a symbol but got:\n";
	    cout << e.op(k).op(0) << endl;
	    exit(0);
	  }

	  //	   cout << e.op(k).op(0) << endl;
	}
      }
      else if(is_a<numeric>(e.op(k))){
	//coe=k;
      }
      else{
	if(is_a<symbol>(e.op(k)))
	  tmp.append(e.op(k));
	else{
	  cout << "failure 2 in monotolist: expecting a symbol but got:\n";
	  cout << e.op(k) << endl;
	  exit(0);
	}
	//	 cout << e.op(k) << endl;
      }
    }
  }
  /* if(coe >=0) */
  /*   tmp.append(coe); */
  /* else */
  /*   tmp.append(1); */

  return tmp;
}

// Converts a monomial of degree r from a list (e.g. x2^2*x4*x5 = (2,2,4,5))
// to exponential form (e.g. x2^2*x4*x5 = (0,2,0,1,1), if numv=5)
// assumes numbering starts at 1. Use with caution
void monovtomonoexp(int *k, int r, int *out, int numv){
  int j;
  inittozero(out,numv);
  for(j=0;j<r;j++)
    (out[k[j]-1])++;
  return;
}

// inhomogeneous version; assumes the first entry in k is the order
void monovtomonoexp0(int *k, int *out, int numv){
  int j;
  //printvec1(k,k[0]+1);
  inittozero(out,numv);
  for(j=1;j<k[0]+1;j++){
    (out[k[j]-1])++;
  }

  //printvec1(out,numv);

  return;
}


//
// convert an expression, assumed to be a monomial
// into a list of integers. Each symbol in the monomial is 
// assumed to be of the form vk, where "v" is a letter and "k"
// an integer. For example 3*x2^2*x4 becomes [2,2,4]
// The coefficient is returned in cf
// This function is only to be used with monomials put in standard form
//
void monotointlist(int *cf, ex e, int **lst, int *r){
  
  int j;
  (*cf)=1;
  *r=0;
  (*lst)=NULL; // NULL returned if expression is zero

  //cerr << "in int version" << endl;
  //cout << "e = " << e << endl;

  if(is_a<add>(e)){
    cout << "failure in monotointlist: not expecting a sum:\n" << e << endl;
    exit(0);
  }

  if(is_a<power>(e)){
    for(j=0;j<e.op(1);j++){
      if(is_a<symbol>(e.op(0))){
	growlist(lst, r);
	(*lst)[(*r)-1] = atoi(((ex_to<symbol>(e.op(0))).get_name()).c_str()+1);
      }
      else{
	cout << "failure 1 in monotointlist: expecting a symbol but got:\n" << e.op(0) << endl;
	exit(0);
      }
    }
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    growlist(lst, r);
    (*lst)[(*r)-1] = atoi(((ex_to<symbol>(e)).get_name()).c_str()+1);
  }
  else if(is_a<numeric>(e)){
    if((ex_to<numeric>(e)).is_integer())
      (*cf)=ex_to<numeric>(e).to_int();
    else{
      fprintf(stderr, "ERROR in monotointlist: expecting an integer but got ");
      cerr << e << endl;
      exit(0);
    }
  }
  else if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  //cout << e.op(k).op(0) << endl;
	  if(is_a<symbol>(e.op(k).op(0))){
	    //cout <<"symbol\n";
	    growlist(lst, r);
	    (*lst)[(*r)-1] = atoi(((ex_to<symbol>(e.op(k).op(0))).get_name()).c_str()+1);
	  }
	  else{
	    cout << "failure 1 in monotointlist: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	    exit(0);
	  }
	 
	}
      }
      else if(is_a<numeric>(e.op(k))){
	if((ex_to<numeric>(e.op(k))).is_integer())
	  (*cf)=ex_to<numeric>(e.op(k)).to_int();
	else{
	  fprintf(stderr, "ERROR in monotointlist: expecting an integer but got ");
	  cerr << e.op(k) << endl;
	  exit(0);
	}
      }
      else{
	//	   cout << e.op(k).op(0) << endl;
	if(is_a<symbol>(e.op(k))){
	  growlist(lst, r);
	  (*lst)[(*r)-1] = atoi(((ex_to<symbol>(e.op(k))).get_name()).c_str()+1);
	}
	else{
	  cout << "failure 2 in monotointlist: expecting a symbol but got:\n" << e.op(k) << endl;
	  exit(0);
	}

      }
    }
  }
  qsortt((*lst),0,(*r)-1);
}


// Overloading: version where coeffs are double
void monotointlist(double *cf, ex e, int **lst, int *r){
  
  int j;
  (*cf)=1;
  *r=0;
  (*lst)=NULL; // NULL returned if expression is zero

  //cerr << "e = " << e << endl;

  if(is_a<add>(e)){
    cout << "failure in monotointlist: not expecting a sum:\n" << e << endl;
    exit(0);
  }

  if(is_a<power>(e)){
    for(j=0;j<e.op(1);j++){
      if(is_a<symbol>(e.op(0))){
	growlist(lst, r);
	(*lst)[(*r)-1] = atoi(((ex_to<symbol>(e.op(0))).get_name()).c_str()+1);
      }
      else{
	cout << "failure 1 in monotointlist: expecting a symbol but got:\n" << e.op(0) << endl;
	exit(0);
      }
    }
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    growlist(lst, r);
    (*lst)[(*r)-1] = atoi(((ex_to<symbol>(e)).get_name()).c_str()+1);
  }
  else if(is_a<numeric>(e))
    (*cf)=ex_to<numeric>(e).to_double();
  else if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      //cerr << e.op(k) << endl;
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  //cout << e.op(k).op(0) << endl;
	  if(is_a<symbol>(e.op(k).op(0))){
	    //cout <<"symbol\n";
	    growlist(lst, r);
	    (*lst)[(*r)-1] = atoi(((ex_to<symbol>(e.op(k).op(0))).get_name()).c_str()+1);
	  }
	  else{
	    cout << "failure 1 in monotointlist: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	    exit(0);
	  }
	 
	}
      }
      else if(is_a<numeric>(e.op(k)))
	(*cf)=ex_to<numeric>(e.op(k)).to_double();
      else{
	//	   cout << e.op(k).op(0) << endl;
	if(is_a<symbol>(e.op(k))){
	  growlist(lst, r);
	  (*lst)[(*r)-1] = atoi(((ex_to<symbol>(e.op(k))).get_name()).c_str()+1);
	}
	else{
	  cout << "failure 2 in monotointlist: expecting a symbol but got:\n" << e.op(k) << endl;
	  exit(0);
	}

      }
    }
  }
  qsortt((*lst),0,(*r)-1);

}


//
// The version of monotointlist where the variables are known 
// and stored in some order (seems not to be used anywhere currently
// but possibly a better way to proceed than standardising)
//

int monotointlist1(ex e, char **vars, int numv, int **lst, int *r){
  
  int j;
  int ret=1;
  *r=0;
  (*lst)=NULL; // NULL returned if expression is zero

  //cout << "e = " << e << endl;

  if(is_a<add>(e)){
    cout << "failure in monotointlist1: not expecting a sum:\n" << e << endl;
    exit(0);
  }

  if(is_a<power>(e)){
    for(j=0;j<e.op(1);j++){
      if(is_a<symbol>(e.op(0))){
	growlist(lst, r);
	(*lst)[(*r)-1]=isinarray(vars, numv, ((ex_to<symbol>(e.op(0))).get_name()).c_str())+1;
      }
      else{
	cout << "failure 1 in monotointlist1: expecting a symbol but got:\n" << e.op(0) << endl;
	exit(0);
      }
	 
    }
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    growlist(lst, r);
    (*lst)[(*r)-1]=isinarray(vars, numv, ((ex_to<symbol>(e)).get_name()).c_str())+1;
  }
  else if(is_a<numeric>(e))
    ret=ex_to<numeric>(e).to_int();
  else if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  //cout << e.op(k).op(0) << endl;
	  if(is_a<symbol>(e.op(k).op(0))){
	    //cout <<"symbol\n";
	    growlist(lst, r);
	    (*lst)[(*r)-1]=isinarray(vars, numv, ((ex_to<symbol>(e.op(k).op(0))).get_name()).c_str())+1;
	  }
	  else{
	    cout << "failure 1 in monotointlist1: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	    exit(0);
	  }
	 
	}
      }
      else if(is_a<numeric>(e.op(k)))
	ret=ex_to<numeric>(e.op(k)).to_int();// discard numbers
      else{
	//	   cout << e.op(k).op(0) << endl;
	if(is_a<symbol>(e.op(k))){
	  growlist(lst, r);
	  (*lst)[(*r)-1]=isinarray(vars, numv, ((ex_to<symbol>(e.op(k))).get_name()).c_str())+1;
	}
	else{
	  cout << "failure 2 in monotointlist1: expecting a symbol but got:\n" << e.op(k) << endl;
	  exit(0);
	}

      }
    }
  }
  qsortt((*lst),0,(*r)-1);
  return ret;
}

// Does not assume prior knowledge of the order
// the first entry in lst is the order of the monomial

void monotointlist0(int *p, ex e, int **lst){
  int i,r,*q;
  monotointlist(p,e,&q,&r);
  (*lst)=(int *)malloc((size_t) (sizeof(int)*(r+1)));
  (*lst)[0]=r;
  for(i=1;i<=r;i++)
    (*lst)[i]=q[i-1];
  free((char*)q);
  return;
}

//overloading to accept coeffs which are doubles
void monotointlist0(double *p, ex e, int **lst){
  int i,r,*q;
  monotointlist(p,e,&q,&r);
  (*lst)=(int *)malloc((size_t) (sizeof(int)*(r+1)));
  (*lst)[0]=r;
  for(i=1;i<=r;i++)
    (*lst)[i]=q[i-1];
  free((char*)q);
  return;
}


// Does not assume prior knowledge of the order
// the first entry in lst is the order of the monomial



// r is the number of monomials
// the return is the degree (except the zero poly also returns 0). 
//(assumed homogeneous)

int polytointlist(ex e, int ***lst, int **cflst, long *r){
  int monord,ord1;
  if(e==0){
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*lst)[0]=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)[0]=0;
    (*lst)[0][0]=0;
    return 0;//probably better as -1
  }
  if(!(is_a<add>(e))){// a single term
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    monotointlist(&((*cflst)[0]),e,lst[0],&monord); //allocates
    return monord;
  }
  (*r)=e.nops(); // terms in sum
  //cout << "r= "<< *r << endl;
  (*lst)=(int **)malloc((size_t) (sizeof(int*)*(*r)));
  (*cflst)=(int *)malloc((size_t) (sizeof(int)*(*r)));
  for (size_t k=0; k!=e.nops(); ++k){
    //cout << "k= "<< k << endl;
    //cout << "first term= "<< e.op(k) << endl;
    monotointlist(&((*cflst)[k]),e.op(k),&((*lst)[k]),&ord1);
    //cout << "firstcoeff= "<< (*cflst)[k] << ", ";
    //cout << "ord1 = " << ord1 << endl;
    //printvec((*lst)[k], ord1);
    if(k==0)
      monord=ord1;
    else if(ord1!=monord){
      cout << "polynomial is not homogeneous" << endl;
      cout << "expected ord = " << monord << " term = " << e.op(k) <<" actual ord = " << ord1 << endl;
      cout << "polynomial = " << e << endl;
      exit(0);
    }
  }
  return monord;
}

//overloading

int polytointlist(ex e, int ***lst, double **cflst, long *r){
  int monord,ord1;
  if(e==0){
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*lst)[0]=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)=(double *)malloc((size_t) (sizeof(double)));
    (*cflst)[0]=0.;
    (*lst)[0][0]=0;
    return 0;
  }
  if(!(is_a<add>(e))){// a single term
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*cflst)=(double *)malloc((size_t) (sizeof(double)));
    monotointlist(&((*cflst)[0]),e,lst[0],&monord); //allocates
    return monord;
  }
  (*r)=e.nops(); // terms in sum
  //cout << "r= "<< *r << endl;
  (*lst)=(int **)malloc((size_t) (sizeof(int*)*(*r)));
  (*cflst)=(double *)malloc((size_t) (sizeof(double)*(*r)));
  for (size_t k=0; k!=e.nops(); ++k){
    //cout << "k= "<< k << endl;
    //cout << "first term= "<< e.op(k) << endl;
    monotointlist(&((*cflst)[k]),e.op(k),&((*lst)[k]),&ord1);
    //cout << "firstcoeff= "<< (*cflst)[k] << ", ";
    //cout << "ord1 = " << ord1 << endl;
    //printvec((*lst)[k], ord1);
    if(k==0)
      monord=ord1;
    else if(ord1!=monord){
      cout << "polynomial is not homogeneous" << endl;
      cout << "expected ord = " << monord << " term = " << e.op(k) <<" actual ord = " << ord1 << endl;
      exit(0);
    }
  }
  return monord;
}

//return the total degree of a monomial
int monodegree(ex e){
  int tot=0;

  if(e==0)
    return -1;

  if(is_a<add>(e)){
    cout << "failure in monodegree: not expecting a sum:\n" << e << endl;
    exit(0);
  }

  if(is_a<power>(e)){
    if(is_a<symbol>(e.op(0))){
      if(is_a<numeric>(e.op(1)))
	return ex_to<numeric>(e.op(1)).to_int();
      else{
	cout << "failure in monodegree: expecting an integer exponent but got " << e.op(1) << endl;
	exit(0);
      }
    }
    cout << "failure 1 in monodegree: expecting a symbol but got:\n" << e.op(0) << endl;
    exit(0);
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    return 1;
  }
  else if(is_a<numeric>(e)){// constant
    return 0;
  }
  else if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<power>(e.op(k))){
	if(is_a<symbol>(e.op(k).op(0))){
	  if(is_a<numeric>(e.op(k).op(1)))
	    tot+=ex_to<numeric>(e.op(k).op(1)).to_int();
	  else{
	    cout << "failure in monodegree: expecting a numerical exponent but got:" << e.op(k).op(1) << endl;
	    exit(0);
	  }
	}
	else{
	  cout << "failure 1 in monodegree: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	  exit(0);
	}
      }
      else if(is_a<symbol>(e.op(k)))
	tot+=1;
      else if(!is_a<numeric>(e.op(k))){
	cout << "failure 2 in monodegree: expecting a number or symbol but got:\n" << e.op(k) << endl;
	exit(0);
      }

    }
  }
  return tot;
}

//degree of a polynomial
int polydegree(ex e){
  int t,deg=-1;
  if(e==0)
    return -1;
  if(!(is_a<add>(e))){// a single term
    return monodegree(e);
  }
  for (size_t k=0; k!=e.nops(); ++k){
    if((t=monodegree(e.op(k)))>deg)
      deg=t;
  }
  return deg;
}

//homogenise a polynomial (assume expanded) using 
//variable var (assumed not to be one of the variables of the polynomial
ex polyhom(ex e, ex var, int offset){
  ex tot=0;
  int t,deg=polydegree(e);
  if(deg<=0 || !(is_a<add>(e)))// constant or a single term
    return e;

  for (size_t k=0; k!=e.nops(); ++k){
    t=monodegree(e.op(k));
    tot+=expand(pow(var,deg-t+offset)*e.op(k));
  }
  return tot;
}



// has non-integer coefficient?

int monohasdcf(ex e){
  if(is_a<add>(e)){
    cout << "failure in monohasdcf: not expecting a sum:\n" << e << endl;
    exit(0);
  }
  if(e==0 || is_a<power>(e) || is_a<symbol>(e))
    return 0;
  if(is_a<numeric>(e) && (ex_to<numeric>(e)).is_integer())
    return 0;

  for (size_t k=0; k!=e.nops(); ++k){
    if(is_a<numeric>(e.op(k))){
      if((ex_to<numeric>(e.op(k))).is_integer())
	return 0;
      else
	return 1;
    }
  }
  return 0;
}


// has noninteger coefficient?
// assume the polynomial is expanded
int polyhasdcfs(ex e){
  if(!(is_a<add>(e)))// a single term
    return monohasdcf(e);

  for (size_t k=0; k!=e.nops(); ++k){
    if(monohasdcf(e.op(k)))
      return 1;
  }
  return 0;
}

//Extract only the coefficient of a monomial (assumed to be integer)
int monocf(ex e){
  if(is_a<add>(e)){
    cout << "failure in monocf: not expecting a sum:\n" << e << endl;
    exit(0);
  }
  if(e==0)
    return 0;
  if(is_a<power>(e) || is_a<symbol>(e))
    return 1;
  if(is_a<numeric>(e)){
    if((ex_to<numeric>(e)).is_integer())
      return ex_to<numeric>(e).to_int();
    cout << "failure in monocf: noninteger coefficient:\n" << e << endl;
    exit(0);
  }

  for (size_t k=0; k!=e.nops(); ++k){
    if(is_a<numeric>(e.op(k))){
      if((ex_to<numeric>(e.op(k))).is_integer())
	return ex_to<numeric>(e.op(k)).to_int();
      cout << "failure in monocf: noninteger coefficient:\n" << e << endl;
      exit(0);
    }
  }
  return 1;
}


//Divide each of coeff. of a poly over Z by the gcd of the coeffs
ex gcdpoly(ex e){
  int g=1;
  int t;
  if(polyhasdcfs(e) || e==0)
    return e;

  if(!(is_a<add>(e)))// a single term
    return expand(e/abs(monocf(e)));

  g=0;
  for (size_t k=0; k!=e.nops(); ++k){//each term
    t=monocf(e.op(k));
    if(t){
      if(!g)//initial
	g=abs(t);
      else
	g=gcd(g,abs(t));
    }
  }

  if(!g)
    return e;

  return expand(e/g);
}



// polynomial not assumed homogeneous. 
// r is the number of monomials. 
// The first entry in each output vector is the degree of the monomial.
// variables assumed to be of the form a single letter followed by a number

void polytointlist0(ex e, int ***lst, int **cflst, long *r){
  //fprintf(stderr, "entering polytointlist0\n");fflush(stdout);fflush(stderr);
  if(e==0){// empty
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*lst)[0]=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)[0]=0;
    (*lst)[0][0]=0;
    //fprintf(stderr, "exiting polytointlist0_1\n");fflush(stdout);fflush(stderr);
    return;
  }
  if(!(is_a<add>(e))){// not a sum
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    monotointlist0(&((*cflst)[0]),e,&((*lst)[0])); //allocates
    //fprintf(stderr, "exiting polytointlist0_2\n");fflush(stdout);fflush(stderr);
    return;
  }
  (*r)=e.nops(); // terms in sum
  //cout << "ri= "<< *r << endl;
  (*lst)=(int **)malloc((size_t) (sizeof(int*)*(*r)));
  (*cflst)=(int *)malloc((size_t) (sizeof(int)*(*r)));
  for (size_t k=0; k!=e.nops(); ++k){
    monotointlist0(&((*cflst)[k]),e.op(k),&((*lst)[k]));
  }
  //fprintf(stderr, "exiting polytointlist0_3\n");fflush(stdout);fflush(stderr);
  return;
}

//overloading to accept coeffs which are doubles
void polytointlist0(ex e, int ***lst, double **cflst, long *r){
  //fprintf(stderr, "entering polytointlist0\n");fflush(stdout);fflush(stderr);
  if(e==0){// empty
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*lst)[0]=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)=(double *)malloc((size_t) (sizeof(double)));
    (*cflst)[0]=0.;
    (*lst)[0][0]=0;
    //fprintf(stderr, "exiting polytointlist0_1\n");fflush(stdout);fflush(stderr);
    return;
  }
  if(!(is_a<add>(e))){// not a sum
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*cflst)=(double *)malloc((size_t) (sizeof(double)));
    monotointlist0(&((*cflst)[0]),e,&((*lst)[0])); //allocates
    //fprintf(stderr, "exiting polytointlist0_2\n");fflush(stdout);fflush(stderr);
    return;
  }
  (*r)=e.nops(); // terms in sum
  //cout << "rd= "<< *r << endl;
  (*lst)=(int **)malloc((size_t) (sizeof(int*)*(*r)));
  (*cflst)=(double *)malloc((size_t) (sizeof(double)*(*r)));
  for (size_t k=0; k!=e.nops(); ++k){
    monotointlist0(&((*cflst)[k]),e.op(k),&((*lst)[k]));
  }
  //fprintf(stderr, "exiting polytointlist0_3\n");fflush(stdout);fflush(stderr);
  return;
}




//
// to convert a symbolic matrix to an integer matrix
//

int sumtointlist(ex e, int **lst){
  
  string str;
  (*lst)=NULL; // NULL returned if expression is zero

  if(e==0){
    (*lst)=(int *)malloc((size_t) (sizeof(int)));
    (*lst)[0]=0;
  }
  else if(is_a<add>(e)){ // a sum
    (*lst)=(int *)malloc((size_t) (sizeof(int)*(e.nops()+1)));
    (*lst)[0] = e.nops();
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<symbol>(e.op(k))){
	str=(ex_to<symbol>(e.op(k))).get_name();
	(*lst)[k+1] = atoi(str.c_str()+1);
      }
      else if(is_a<symbol>(-(e.op(k)))){
	str=(ex_to<symbol>(-(e.op(k)))).get_name();
	(*lst)[k+1] = -atoi(str.c_str()+1);
      }
      else{ // something unexpected in the sum
	cout << "failure in sumtointlist: not expecting this:\n";
	cout << e.op(k) << endl;
	exit(0);
      }
    }
    return e.nops()+1;
  }
  else if(is_a<symbol>(e)){ // a positive symbol
    (*lst)=(int *)malloc((size_t) (sizeof(int)*(2)));
    (*lst)[0] = 1;
    str=(ex_to<symbol>(e)).get_name();
    (*lst)[1] = atoi(str.c_str()+1);
    return 2;
  }
  else if(is_a<symbol>(-e)){ // a negative symbol
    (*lst)=(int *)malloc((size_t) (sizeof(int)*(2)));
    (*lst)[0] = 1;
    str=(ex_to<symbol>(-e)).get_name();
    (*lst)[1] = -atoi(str.c_str()+1);
    return 2;
  }
  else{ // something else
    cout << "failure 2 in sumtointlist: not expecting this:\n";
    cout << e << endl;
    exit(0);
  }
 
  return 0;
}

void checkmonotolist(){
  lst l;
  symbol x("x"), y("y");
  ex e=2*x*x*pow(y, 2);
  l=monotolist(e);
  cout << l << endl;
  return;
}

void checkmonotointlist(){
  int *l;
  int r, k;
  symbol x1("x1"), x2("x2");
  ex e=-32*x1*x1*pow(x2, 5);
  monotointlist(&r,e, &l, &k);
  printvec(l, k);
  free ((char *)l);
  return;
}

// A list of monomials into exponential form, with sorting
int **monlisttoexplist(int **mons, int dim, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp(mons[k],dim,expsout[k],numv);
  }

  if(srt)
    qsorti(expsout,numv,0,len-1);
  return expsout;
}

// A list of terms (monomials + coeffs) into exponential form
// The point is that the cfs need to be sorted too
int **moncflisttoexplist(int **mons, int dim, int *cfs, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp(mons[k],dim,expsout[k],numv);
  }

  if(srt)
    qsorti2(expsout,cfs,numv,0,len-1);
  return expsout;
}

//overload to accept double coeffs
int **moncflisttoexplist(int **mons, int dim, double *cfs, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp(mons[k],dim,expsout[k],numv);
  }

  if(srt)
    qsorti2(expsout,cfs,numv,0,len-1);
  return expsout;
}

// don't assume homogeneous polynomials

int **moncflisttoexplist0(int **mons, int *cfs, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp0(mons[k],expsout[k],numv);
  }

  if(srt)
    qsorti2(expsout,cfs,numv,0,len-1);
  return expsout;
}

//overload to accept double cfs

int **moncflisttoexplist0(int **mons, double *cfs, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp0(mons[k],expsout[k],numv);
  }
  if(srt)
    qsorti2(expsout,cfs,numv,0,len-1);
  return expsout;
}

// Extract the variable names from a monomial as strings
void monovars(ex e, char ***vars, int *numvars){
  int j;
  if(is_a<power>(e)){
    for(j=0;j<e.op(1);j++){
      if(is_a<symbol>(e.op(0)))
	(*numvars)=addv1((*numvars), (ex_to<symbol>(e.op(0))).get_name().c_str(), vars);
      else{
	cout << "failure 1 in monovars: expecting a symbol but got:\n" << e.op(0) << endl;
	exit(0);
      }
    }
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    (*numvars)=addv1((*numvars), (ex_to<symbol>(e)).get_name().c_str(), vars);
  }
  else{//must be a product
    for (size_t k=0; k!=e.nops(); ++k){//each part
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  if(is_a<symbol>(e.op(k).op(0)))
	    (*numvars)=addv1((*numvars), (ex_to<symbol>(e.op(k).op(0))).get_name().c_str(), vars);
	  else{
	    cout << "failure 1 in monovars: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	    exit(0);
	  }
	}
      }
      else if(is_a<symbol>(e.op(k))){
	(*numvars)=addv1((*numvars), (ex_to<symbol>(e.op(k))).get_name().c_str(), vars);
      }
      else if(!(is_a<numeric>(e.op(k)))){
	cout << "failure 2 in monovars: expecting a number or symbol but got:\n" << e << endl << e.op(k) << endl;
	exit(0);
      }

    }
  }

  return;
}

// Extract the variable names from a monomial as strings 
// and return the monomial in canonical form
ex monosimp(ex e, char ***vars, int *numvars){
  int p;
  ex ex1=1;
  char str1[10];
  if(is_a<power>(e)){
    if(is_a<symbol>(e.op(0))){
      (*numvars)=addv((*numvars), (ex_to<symbol>(e.op(0))).get_name().c_str(), vars, &p);
      sprintf(str1, "v%d", p+1);
      ex1*=pow(get_possymbol(str1), e.op(1));
      //cout << e.op(0) << " --> " << p+1 << endl;
    }
    else{
      cout << "failure 1 in monosimp: expecting a symbol but got:\n" << e.op(0) << endl;
      exit(0);
    }
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    (*numvars)=addv((*numvars), (ex_to<symbol>(e)).get_name().c_str(), vars, &p);
    sprintf(str1, "v%d", p+1);
    ex1*=get_possymbol(str1);
    //cout << e << " --> " << p+1 << endl;
  }
  else if(is_a<numeric>(e)){
    ex1*=e;
  }
  else{//must be a product
    for (size_t k=0; k!=e.nops(); ++k){//each part
      if(is_a<power>(e.op(k))){
	if(is_a<symbol>(e.op(k).op(0))){
	  (*numvars)=addv((*numvars), (ex_to<symbol>(e.op(k).op(0))).get_name().c_str(), vars, &p);
	  sprintf(str1, "v%d", p+1);
	  ex1*=pow(get_possymbol(str1),e.op(k).op(1));
	  //cout << e.op(k).op(0) << " --> " << p+1 << endl;
	}
	else{
	  cout << "failure 1 in monosimp: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	  exit(0);
	}
      }
      else if(is_a<symbol>(e.op(k))){
	(*numvars)=addv((*numvars), (ex_to<symbol>(e.op(k))).get_name().c_str(), vars,&p);
	sprintf(str1, "v%d", p+1);
	ex1*=get_possymbol(str1);
	//cout << e.op(k) << " --> " << p+1 << endl;
      }
      else if(is_a<numeric>(e.op(k))){
	ex1*=e.op(k);
      }
      else{
	cout << "failure 2 in monosimp: expecting a number or symbol but got:\n" << e << endl << e.op(k) << endl;
	exit(0);
      }

    }
  }
  return ex1;
}


//extract the actual variables (not variable names) from a monomial
//We assume that we already know the variable names (e.g. via monovars)
void monovarss(ex e, ex *vars, char **varnames, int numv){
  int j,tot=0,q;
  if(is_a<power>(e)){
    for(j=0;j<e.op(1);j++){
      if(is_a<symbol>(e.op(0)) && (q=isinarray(varnames, numv, ex_to<symbol>(e.op(0)).get_name().c_str()))!=-1){
	vars[q]=ex_to<symbol>(e.op(0));
	if(tot>numv){
	  cout << "failure in monovarss: too many variables\n";
	  exit(0);
	}
      }
      else{
	cout << "failure 1 in monovarss: expecting a symbol but got:\n" << e.op(0) << endl;
	exit(0);
      }
    }
  }
  else if(is_a<symbol>(e) && (q=isinarray(varnames, numv, ex_to<symbol>(e).get_name().c_str()))!=-1){ // single symbol - not a product
    vars[q]=ex_to<symbol>(e);
    if(tot>numv){
      cout << "failure in monovarss: too many variables\n";
      exit(0);
    }
  }
  else{//must be a product
    for (size_t k=0; k!=e.nops(); ++k){//each part
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  if(is_a<symbol>(e.op(k).op(0)) && (q=isinarray(varnames, numv, ex_to<symbol>(e.op(k).op(0)).get_name().c_str()))!=-1){
	    vars[q]=ex_to<symbol>(e.op(k).op(0));
	    if(tot>numv){
	      cout << "failure in monovarss: too many variables\n";
	      exit(0);
	    }
	  }
	  else{
	    cout << "failure 1 in monovars: expecting a symbol but got:\n" << e.op(k).op(0) << endl;
	    exit(0);
	  }
	}
      }
      else if(is_a<symbol>(e.op(k)) && (q=isinarray(varnames, numv, ex_to<symbol>(e.op(k)).get_name().c_str()))!=-1){
	vars[q]=ex_to<symbol>(e.op(k));
	if(tot>numv){
	  cout << "failure in monovarss: too many variables\n";
	  exit(0);
	}
      }
      else if(!(is_a<numeric>(e.op(k)))){
	cout << "failure 2 in monovars: expecting a number or a symbol but got:\n" << e << endl << e.op(k) << endl;
	exit(0);
      }
    }
  }

  return;
}

// Extract variable names from a polynomial
int polyvars(ex e, char ***vars){
  int numv=0;
  if(e==0){
    (*vars)=NULL;
    return 0;
  }
  if(!(is_a<add>(e)))// a single term
    monovars(e,vars,&numv);
  else{
    for (size_t k=0; k!=e.nops(); ++k)
      monovars(e.op(k),vars, &numv);
  }
  return numv;
}

// Overloading: Extract variable names from a polynomial
// assume an existing list of variables
int polyvars(ex e, char ***vars, int *numv){
  if(e==0)
    return (*numv);
  
  if(!(is_a<add>(e)))// a single term
    monovars(e,vars,numv);
  else{
    for (size_t k=0; k!=e.nops(); ++k)
      monovars(e.op(k),vars, numv);
  }
  return (*numv);
}

// (See monosimp)
ex polysimp(ex e, char ***vars, int *numv){
  ex extot=0;
  if(e==0){
    (*vars)=NULL;
    return 0;
  }
  if(!(is_a<add>(e)))// a single term
    extot=monosimp(e,vars,numv);
  else{
    for (size_t k=0; k!=e.nops(); ++k)
      extot+=monosimp(e.op(k),vars, numv);
  }
  return extot;
}


// Get the variable names where we expect these to be of the form vi, 
// numv is the maximum value of i
int polyvars1(ex e, char ***vars){
  int numv=0, vmax, vmin;
  int i, t,flag=0,vflag=0;
  char str1[5];
  //fprintf(stderr, "entering polyvars1\n");

  if(e==0){
    (*vars)=NULL;
    return 0;
  }
  if(!(is_a<add>(e)))// a single term
    monovars(e,vars,&numv);
  else{
    for (size_t k=0; k!=e.nops(); ++k)
      monovars(e.op(k),vars, &numv);
  }

  for(i=0;i<numv;i++){
    if((*vars)[i][0]=='v' && ispureint((*vars)[i]+1)){
      t=atoi((*vars)[i]+1);
      //fprintf(stderr, "t=%d\n", t);
      if(!vflag){vmin=t;vmax=t;}
      else{
	if(t<vmin)
	  vmin=t;
	if(t>vmax)
	  vmax=t;
      }
      vflag=1;
    }
    else{//variable not of this form
      flag=1;
      break;
    }
  }
  if(vmin==1 && vmax==numv)
    return numv;
  else if(flag){
    fprintf(stderr, "Expecting variables of the form v(integer), but got %s. EXITING.\n", (*vars)[i]);exit(0);
  }
  else if(vmin<=0){
    fprintf(stderr, "Expecting variable indices to be positive, but got v%d. EXITING.\n", vmin);exit(0);
  }

  //Some missing variables
  for(i=0;i<numv;i++)
    free((char*)(*vars)[i]);
  free((char*)(*vars));

  (*vars)=(char**)malloc((size_t) ((vmax)*sizeof(char*)));
  for(i=1;i<=vmax;i++){
    sprintf(str1, "v%d", i);
    (*vars)[i-1]=strdup(str1);
  }
  //fprintf(stderr, "exiting polyvars1\n");
  return vmax;

}

// Get the variables from a polynomial, given variable names
void polyvarss(ex e, ex *vars, char **varnames, int numv){
  if(e==0)
    return;
  if(!(is_a<add>(e)))// a single term
    monovarss(e,vars,varnames,numv);
  else{
    for (size_t k=0; k!=e.nops(); ++k)
      monovarss(e.op(k),vars,varnames,numv);
  }
  return;
}

// From a polynomial in standard form to the same in exponential form
void extolst(ex tmp, int numv, int ***explst, int **cflst, long *r){
  int **lst;
  polytointlist0(tmp, &lst, cflst, r);
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv+1,1);// redefines cflst
  //elongate(explst1,cflst1,r1,numv);
  ifree(lst, *r);
  return;
}

//overload to accept double cfs
void extolst(ex tmp, int numv, int ***explst, double **cflst, long *r){
  int **lst;
  polytointlist0(tmp, &lst, cflst, r);
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv+1,1);// redefines cflst
  //elongate(explst1,cflst1,r1,numv);
  ifree(lst, *r);
  return;
}



void extolst1(ex tmp, int numv, int ***explst, int **cflst, long *r, int *mindegree, int *maxdegree){
  int **lst;
  long s;
  //fprintf(stderr, "gt here1\n");fflush(stdout);fflush(stderr);
  polytointlist0(tmp, &lst, cflst, r);
  //fprintf(stderr, "gt here2\n");fflush(stdout);fflush(stderr);
  for(s=0;s<(*r);s++){
    //printvec(lst[s],numv+1);
    if(s==0){
      (*mindegree)=lst[s][0];
      (*maxdegree)=lst[s][0];
    }
    else{
      (*maxdegree)=lst[s][0]>(*maxdegree)?lst[s][0]:(*maxdegree);
      (*mindegree)=lst[s][0]<(*mindegree)?lst[s][0]:(*mindegree);
    }
  }
  //fprintf(stderr, "gt here3\n");fflush(stdout);fflush(stderr);
  //sorted
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv+1,1);// redefines cflst
  //elongate(explst1,cflst1,r1,numv);
  //fprintf(stderr, "gt here4 r=%ld\n",*r);fflush(stdout);fflush(stderr);
  ifree(lst, *r);
  //fprintf(stderr, "gt here5\n");fflush(stdout);fflush(stderr);
  return;
}


//overload to accept double cfs

void extolst1(ex tmp, int numv, int ***explst, double **cflst, long *r, int *mindegree, int *maxdegree){
  int **lst;
  long s;
  polytointlist0(tmp, &lst, cflst, r);
  for(s=0;s<(*r);s++){
    if(s==0){
      (*mindegree)=lst[s][0];(*maxdegree)=lst[s][0];
    }
    else{
      (*maxdegree)=lst[s][0]>(*maxdegree)?lst[s][0]:(*maxdegree);
      (*mindegree)=lst[s][0]<(*mindegree)?lst[s][0]:(*mindegree);
    }
  }
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv+1,1);// redefines cflst
  ifree(lst, *r);
  return;
}

ex intltomono(int cf, int *l, int n){
  char str1[10];
  ex ex1=1;
  int i;
  ex1*=cf;
  for(i=0;i<n;i++){
    sprintf(str1, "v%d", l[i]);
    ex1*=get_possymbol(str1);
  }
  return ex1;
}

ex expltomono(int cf, int *l, int n){
  char str1[10];
  ex ex1=1;
  int i;
  ex1*=cf;
  for(i=0;i<n;i++){
    if(l[i]!=0){
      sprintf(str1, "v%d", i+1);
      ex1*=pow(get_possymbol(str1),l[i]);
    }
  }
  return ex1;
}

void expltomonoprint(int cf, int *l, int n){
  int i;
  if(cf!=1)
    printf("%d ", cf);
  for(i=0;i<n;i++){
    if(l[i]!=0){
      if(l[i]==1)
	printf("v%d ", i+1);
      else
	printf("v%d^%d ", i+1, l[i]);
    }
  }
  printf("\n");
}

//Simplifies polynomial removing absent variables, and by grouping variables
int polylistsimp(int **lst, int numv, long r, int ***lstnew, int mkhom, int tmppolymindeg, int tmppolymaxdeg, int *polymindeg, int *polymaxdeg){
  int i,j,sameflag,zeroflag,flg=0, mindeg=1000,maxdeg=0;
  long s;
  int deg[r];
  int newlst[numv];
  inittoone(newlst, numv);
  for(i=0;i<numv;i++){
    for(s=0;s<r;s++){
      zeroflag=1;
      if(lst[s][i]){
	zeroflag=0;
	break;
      }
    }
    if(zeroflag){//zero column
      newlst[i]=0;
      flg=1;
    }
  }

  for(i=0;i<numv-1;i++){
    for(j=i+1;j<numv;j++){
      if(newlst[j]){
	sameflag=1;
	for(s=0;s<r;s++){
	  if(lst[s][i]!=lst[s][j]){
	    sameflag=0;
	    break;
	  }
	}
	if(sameflag){//repeated column
	  newlst[j]=0;
	  flg=1;
	}
      }
    }
  }

  (*lstnew)=(int **)malloc((size_t) (r*sizeof(int*)));
  for(s=0;s<r;s++)
    (*lstnew)[s]=(int *)malloc((size_t) ((numv+1)*sizeof(int)));

  j=0;
  for(i=0;i<numv;i++){
    if(newlst[i]){
      for(s=0;s<r;s++)
	(*lstnew)[s][j]=lst[s][i];
      j++;
    }
  }

  //Find mindeg and maxdeg
  if(flg){
    for(s=0;s<r;s++){
      deg[s]=0;
      for(i=0;i<j;i++)
	deg[s]+=(*lstnew)[s][i];
      if(deg[s]>maxdeg)
	maxdeg=deg[s];
      if(deg[s]<mindeg)
	mindeg=deg[s];
    }
  }
  else{//nothing changed, no need to recompute degrees
    for(s=0;s<r;s++){
      deg[s]=0;
      for(i=0;i<j;i++)
	deg[s]+=(*lstnew)[s][i];
    }
    mindeg=tmppolymindeg;
    maxdeg=tmppolymaxdeg;
  }

  if(mkhom){//want homogeneous
    if(mindeg!=maxdeg){//want homogeneous and not already homogeneous
      for(s=0;s<r;s++){
	(*lstnew)[s][j]=maxdeg-deg[s];
	//fprintf(stderr, "here, added %d\n", maxdeg-deg[s]);
      }
      j++;
      (*polymindeg)=maxdeg;
      (*polymaxdeg)=maxdeg;
    }
    else{
      (*polymindeg)=maxdeg;
      (*polymaxdeg)=maxdeg;
    }
  }
  else{
    (*polymindeg)=mindeg;
    (*polymaxdeg)=maxdeg;
  }

  return j;//new number of variables
}

//A list to a polynomial as a symbolic object

ex lrexptopolysimp(int **l, int *cfs, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    ex1=cfs[i];
    for(j=0;j<n;j++){ // each symbol
      if(l[i][j]){
	sprintf(str1, "v%d", j+1);
	ex1*=pow(get_possymbol(str1),l[i][j]);
      }
    }
    extot+=ex1;
  }
  return extot;
}

//overloaded
ex lrexptopolysimp(int **l, double *cfs, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    ex1=cfs[i];
    for(j=0;j<n;j++){ // each symbol
      if(l[i][j]){
	sprintf(str1, "v%d", j+1);
	ex1*=pow(get_possymbol(str1),l[i][j]);
      }
    }
    extot+=ex1;
  }
  return extot;
}

//just the monomials (coeffs all 1)
ex lrexptopolystruct(int **l, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    ex1=1;
    for(j=0;j<n;j++){ // each symbol
      if(l[i][j]){
	sprintf(str1, "v%d", j+1);
	ex1*=pow(get_possymbol(str1),l[i][j]);
      }
    }
    extot+=ex1;
  }
  return extot;
}

//overloaded: adding a tolerance
ex lrexptopolysimp(int **l, double *cfs, double tol, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    if(fabs(cfs[i])>tol){
      ex1=cfs[i];
      for(j=0;j<n;j++){ // each symbol
	if(l[i][j]){
	  sprintf(str1, "v%d", j+1);
	  ex1*=pow(get_possymbol(str1),l[i][j]);
	}
      }
      extot+=ex1;
    }
  }
  return extot;
}


//Simplify any polynomial (put in canonical variables)
//if mkhom, then make homogeneous
//final argument is if the poly has noninteger coeffs

ex polyhomsimp(ex tmp1, int *numvin, int *numvout, int *deg, int mkhom, int dcfs, int debug){
  char **pvars;
  int numv=0,numv1=0;//=polyvars1(tmp, &pvars);
  int **lst=NULL;
  // int or double?
  int *icflst=NULL;
  double *dcflst=NULL;
  int **tmplst;
  long r;
  int tmppolymindeg,tmppolymaxdeg;
  int polymindeg,polymaxdeg;
  ex tmp;
  ex tmpout;
  int debugfull=(debug<=0)?0:debug-1;


  if(tmp1==0){(*numvin)=0;(*numvout)=0;(*deg)=0;return 0;}
  if(is_a<numeric>(tmp1)){(*numvin)=0;(*numvout)=0;(*deg)=0;return tmp1;}

  if(debug){fprintf(stderr, "\n###Entering polyhomsimp.\n");}
  if(debug){cerr << "Poly before simplifying: " << tmp1 << endl;}
  tmp=polysimp(tmp1, &pvars, &numv1);//canonical variables v#
  if(debugfull){
    cerr << "Poly after applying \"polysimp\": " << tmp << endl;
    fprintf(stderr, "variable mapping:\n");
    for(r=0;r<numv1;r++)
      fprintf(stderr, "%s-->v%ld\n",pvars[r], r+1);
  }

  (*numvin)=numv1;

  if(dcfs)
    extolst1(tmp, numv1, &tmplst, &dcflst, &r, &tmppolymindeg, &tmppolymaxdeg);
  else
    extolst1(tmp, numv1, &tmplst, &icflst, &r, &tmppolymindeg, &tmppolymaxdeg);

  if(debugfull){
    fprintf(stderr, "Collapsing variables%s. Before manipulation: number of variables=%d, degree=%d:\n", mkhom?" and homogenising":"", numv1, tmppolymaxdeg);
  }

  //This can reduce the degree of a polynomial by collapsing variables. 
  numv=polylistsimp(tmplst, numv1, r, &lst, mkhom, tmppolymindeg, tmppolymaxdeg, &polymindeg, &polymaxdeg);
  (*deg)=polymaxdeg;

  if(debugfull)
    fprintf(stderr, "After manipulation, number of variables=%d, degree=%d\n",numv,polymaxdeg);

  if(dcfs)
    tmpout=lrexptopolysimp(lst, dcflst, r, numv);
  else
    tmpout=lrexptopolysimp(lst, icflst, r, numv);

  (*numvout)=numv;

  freearraydat(pvars, numv1);
  ifree(lst, r);
  ifree(tmplst,r); 
  if(icflst)
    free((char*)icflst);
  if(dcflst)
    free((char*)dcflst);

  if(debug)
    cerr << "Exiting polyhomsimp. Final polynomial:\n" << tmpout << endl;
  fflush(stderr);

  return tmpout;
}


// takes a polynomial from a string
// expands. Homogenises if desired (mkhom flag)
// puts in standard format

ex polyfromstring(const char *str, int *numvin, int *numv, int *deg, int mkhom, int *dcfs, int debug){
  parser p;
  ex tmp= p(str);
  ex tmp1=expand(tmp);
  *dcfs=polyhasdcfs(tmp1);
  return polyhomsimp(expand(tmp), numvin, numv, deg, mkhom, *dcfs, debug);
}

// polynomial from a file stored in a string
// expands. Homogenises if desired (mkhom flag)
// No error checking

ex polyfromfile(const char filename[], int *numvin, int *numv, int *deg, int mkhom, int *dcfs, int debug){
  ex tmp=0;
  char *str=readfileintostr(filename);
  char *line;
  long pos=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      tmp=polyfromstring(line,numvin,numv,deg,mkhom,dcfs,debug);
      free(line);free(str);
      return tmp;
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);free(str);
  return tmp;
}


// List of polys in a file, one per line, must end with an empty line
// Put all the monomials in a matrix (numv X totmons), one per col
// put all the coeffs in a matrix (totpols X totmons), i.e.,
// col=monomial index, row=poly number
// Convenient for Mixed Volume calculation
void polysfromfile(const char fname[], int ***mons, double ***cfs, int *numv, long *totmons, int *totpols, int debug){
  FILE *fd;
  char oneline[5000];
  ex tmp, tmp1;
  ex structpoly=0;
  parser p;
  char **pvars;
  int i;
  long r, r1, s, t;
  int **lst, **lst1;
  double *cflst, *cflst1;
  (*totpols)=0;
  int totp=0;
  (*numv)=0;

  //Parse 1: get variable names
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "Error in polysfromfile: file %s couldn't be opened. EXITING\n", fname);exit(0);
  }

  while(getline0(fd, oneline, 5000) > 0){
    if(!iscomline(oneline)){
      tmp=expand(p(oneline));
      polyvars(tmp, &pvars, numv);
      (*totpols)++;
      if(debug){
	cerr<< "Poly " << (*totpols) << ": " << tmp << endl;
      }
    }
  }

  if(debug){
    fprintf(stderr, "variables:");
    for(i=0;i<(*numv);i++){fprintf(stderr, " %s",pvars[i]);}
    fprintf(stderr, "\n");
  }
  fclose(fd);

  //second parse
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "Error in polysfromfile: file %s couldn't be opened. EXITING\n", fname);exit(0);
  }
  while(getline0(fd, oneline, 5000) > 0){
    if(!iscomline(oneline)){
      tmp=expand(p(oneline));
      tmp1=polysimp(tmp, &pvars, numv);//standard form
      //cout << tmp << "\n" <<tmp1 << endl;
      extolst(tmp1, (*numv), &lst, &cflst, &r);
      structpoly+=lrexptopolystruct(lst,r,(*numv));
      //cerr<< structpoly << endl;
      ifree(lst,r);free((char*)cflst);
    }
  }
  extolst(structpoly, (*numv), &lst, &cflst, &r);
  fclose(fd);
  *totmons=r;
  (*mons)=imatrix(0,(*numv)-1,0,(*totmons)-1);
  (*cfs)=dmatrix(0,(*totpols)-1,0,(*totmons)-1);

  for(i=0;i<(*numv);i++){//mons is a matrix rather than a list of vectors
    for(s=0;s<(*totmons);s++){
      (*mons)[i][s]=lst[s][i];
    }
  }
  for(i=0;i<(*totpols);i++){//initialise coefficient matrix
    for(s=0;s<(*totmons);s++){
      (*cfs)[i][s]=0.0;
    }
  }
  
  //third parse
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "Error in polysfromfile: file %s couldn't be opened. EXITING\n", fname);exit(0);
  }
  while(getline0(fd, oneline, 5000) > 0){
    if(!iscomline(oneline)){
      tmp=expand(p(oneline));
      tmp1=polysimp(tmp, &pvars, numv);//standard form
      extolst(tmp1, (*numv), &lst1, &cflst1, &r1);
      for(s=0;s<r1;s++){//must be there already
	t=isinarray(lst,(*totmons),lst1[s],(*numv));
	//fprintf(stderr, "(*totmons)=%d,t=%d\n",(int)(*totmons),(int)t);
	(*cfs)[totp][t]=cflst1[s];
      }
      totp++;
      ifree(lst1,r1);free((char*)cflst1);
    }
  }

  fclose(fd);
  ifree(lst,r);free((char*)cflst);
  freearraydat(pvars, (*numv));
  return;
}


// vector of unknowns; for compatibility, we stick with "v"
// and begin numbering at 1
matrix makexvec(int n){
  int i;
  ex x[n];
  char sstr[5];
  matrix sols(n,1);

  for(i=0;i<n;i++){
    sprintf(sstr,"v%d",i+1);
    x[i]=get_possymbol(sstr);
    sols(i,0)=x[i];
  }
  return sols;
}

//overloading: j is initial offset. so makexvec(n)==makexvec(1, n)
matrix makexvec(int j, int n){
  int i;
  ex x[n];
  char sstr[5];
  matrix sols(n,1);

  for(i=0;i<n;i++){
    sprintf(sstr,"v%d",i+j);
    x[i]=get_possymbol(sstr);
    sols(i,0)=x[i];
  }
  return sols;
}

//diagonal matrix instead of vector
matrix makediagxmat(int k, int n, bool minus){
  int i,j;
  ex x[n];
  char sstr[5];
  matrix sols(n,n);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i==j){
	sprintf(sstr,"v%d",i+k);
	x[i]=get_possymbol(sstr);
	if(minus)
	  sols(i,j)=-x[i];
	else
	  sols(i,j)=x[i];
      }
      else
	sols(i,j)=0;
    }
  }
  return sols;
}

matrix makexvect(int j, int n){
  int i;
  ex x[n];
  char sstr[5];
  matrix sols(1,n);

  for(i=0;i<n;i++){
    sprintf(sstr,"v%d",i+j);
    x[i]=get_possymbol(sstr);
    sols(0,i)=x[i];
  }
  return sols;
}



matrix makexvecinv(int n){
  int i;
  ex x[n];
  char sstr[5];
  matrix sols(n,1);

  for(i=0;i<n;i++){
    sprintf(sstr,"v%d",i+1);
    x[i]=get_possymbol(sstr);
    sols(i,0)=1/x[i];
  }
  return sols;
}

//vector of unknowns where each entry is the product of all with one removed
matrix makexveccomp(int n){
  int i;
  ex tmp=1,x[n];
  char sstr[5];
  matrix sols(n,1);

  for(i=0;i<n;i++){
    sprintf(sstr,"v%d",i+1);
    x[i]=get_possymbol(sstr);
    tmp*=x[i];
  }

  for(i=0;i<n;i++)
    sols(i,0)=tmp/x[i];

  return sols;
}
