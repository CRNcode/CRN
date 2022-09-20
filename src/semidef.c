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
#include "symbolic.h"
#include "semidef.h"
#include <limits.h>

// generate all the basic monomials on numv variables
// return the total number of monomials generated (2^numv)

int genmod(int numv, int **mons, int q){
  int i,j;
  int xc[numv];
  int flag=1;
  int tot=0;

  inittozero(mons[tot],numv);
  if(!q){printvec(mons[tot],numv);}
  tot++;
  for(i=1;i<=numv;i++){
    flag=1;
    firstcomb(xc,numv,i);
    while(flag){
      inittozero(mons[tot],numv);
      for(j=0;j<i;j++)
	mons[tot][xc[j]]=1;
      if(!q){printvec(mons[tot],numv);}
      tot++;
      flag=nextcomb(xc,numv,i);
    }
  }
  return tot;
}

/////////////////////////

//SDP return codes: 0: feasible; -1: infeasible; 2: probably feasible (reduced accuracy); -2: something went wrong
int ppp1(int **lst, double *cflst, long r, int numv, int polymindeg, int polymaxdeg, int preprocess, int ldegreemin, int ldegree, int strict, int debug){
  int p,q,j,j1,jc,k,l,blockdeg,blockmindeg,blockmaxdeg,flag;
  // Basic monomials: monomials without repetition on numv variables (including the empty monomial, namely, 1)
  int totblk=(int)pow(2.0, (double)numv);
  int tot=0;//counts the blocks
  long s;
  int **lst2;
  double **lst1;
  int LHtot;
  int **blks;
  double **blksout;
  int **blkinds;
  //stores the basic monomials which figure outside squares
  int **mons=imatrix(0, totblk-1, 0,numv-1);
  int totmons=0;
  int **allmons=NULL;
  int newmon[numv];
  //monomials which can figure in squares
  int t,sqmaxdeg;
  int **sqmons, **lmons;
  int offset;
  int numsqmon;
  int ind;
  int blksz, newblksz;
  int **constraints;
  int numconstraints=0;
  int **newconstraints=NULL;
  int numnewconstraints=0;
  int numtrueblk,numsqblk,numallblk;
  int *trueblksz;
  FILE *fd, *fp;
  int totlmons;
  int **posmons;
  int numposmons;
  int *nondiag;
  int **monrefs;
  int blkind, rowind, monind;
  char *outstr, *datstr;
  double *outdat;
  int datlen;
  int outpolylen=0;
  int **outpoly;
  double *outcfs;
  int outblkpolylen;
  int **outblkpoly;
  double *outblkcfs;
  double tol=0.0;
  double maxval;
  int posleft=0;//1 means LHS multiplier has only positive coefficients
  int retstat=-2;
  char path[2048];
  double maxcf;
  int debugfull=(debug<=0)?0:debug-1;

  if(debugfull){
    fprintf(stderr, "Examining the polynomial: \n");
    cerr << lrexptopolysimp(lst, cflst, r, numv) << endl;
    //exit(0);
  }

  if(!strict)//otherwise possible artifacts
    posleft=1;

  /* fprintf(stderr, "%d\n", NewtonPolyVertSgn(tmp)); */
  /* exit(0); */
  
  //monomials on LHS
  lmons=genmons1(numv,ldegreemin,ldegree,&totlmons);
  /* fprintf(stderr, "totlmons=%d\n", totlmons); */
  /* for(s=0;s<totlmons;s++) */
  /*   printvec(lmons[s],numv); */

  // The LHS poly is the product of original poly P and parameterised poly Q of degrees up to ldegree
  // lst1 stores the coeffs of the LHS poly which are themselves monomials in the additional LHS variables
  // lst2 stores the monomials of the LHS poly
  LHtot=polyprod(lmons, totlmons, cflst, lst, r, numv, &lst1, &lst2);

  //increase the maximum degree and minimum degrees of the polynomial
  polymaxdeg+=ldegree;
  polymindeg+=ldegreemin;
  if(debug){fprintf(stderr, "mindeg=%d, maxdeg=%d\n", polymindeg, polymaxdeg);}

  // posmons stores all monomials involving numv variables
  // in the range between the minimum and maximum degrees of the 
  // LHS polynomial. These are all the possible monomials which could 
  // figure on the RHS. This does not mean that they will.
  posmons=genmons1(numv, polymindeg, polymaxdeg, &numposmons);

  //Vector to hold the number of times that a monomial figures
  //but not as a diagonal entry in a block
  nondiag=(int *)malloc((size_t) (numposmons*sizeof(int)));
  inittozero(nondiag,numposmons);

  // Array of vectors of length 3k+1 which tells us where each monomial
  // occurs in blocks on the RHS (doesn't list occurrence on the LHS)
  monrefs=(int **)malloc((size_t) (numposmons*sizeof(int*)));

  for(s=0;s<numposmons;s++){// add all the monomials in range (this fixes the monomial order)
    totmons=addoneintstra(totmons, posmons[s], numv, &allmons,&ind);
    monrefs[s]=(int *)malloc((size_t) (4*sizeof(int)));
    //total occurrences, block number, then indices.
    monrefs[s][0]=1;monrefs[s][1]=totblk+s;monrefs[s][2]=0;monrefs[s][3]=0;
  }

  if(debugfull){fprintf(stderr, "LHS polynomial:\n");}
  for(s=0;s<LHtot;s++){// add the monomials from the LHS poly, but don't reference (i.e., not in a block)
    totmons=addoneintstra(totmons, lst2[s], numv, &allmons,&ind);
    (nondiag[ind])++;//LHS monomials can't be removed from blocks
    if(debugfull){
      printvec(lst1[s],totlmons);
      printvec(lst2[s],numv);
    }
  }

  //From now on the monomials in allmons are in the same order as those in posmons
  //One block for each basic monomial (i.e., 2^numv monomials where powers are 0 and 1)
  //+ one block for each allowed monomial

  blks=(int **)malloc((size_t) ((totblk+numposmons)*sizeof(int*)));
  blkinds=(int **)malloc((size_t) ((totblk+numposmons)*sizeof(int*)));//stuff that will survive

  sqmaxdeg=polymaxdeg/2;
  t=numv+sqmaxdeg;//largest degree
  sqmons=imatrix(0, choose(t,numv)-1,0,numv-1);
  genmons(numv,sqmaxdeg,sqmons);
  genmod(numv,mons,1);
  //
  // Initial RHS blocks
  // basic monomial of each dimension
  //
  for(blockdeg=0;blockdeg<=numv;blockdeg++){//Each dimension
    jc=choose(numv,blockdeg);
    //mons[blockdeg]=imatrix(0, jc-1, 0,numv-1);//basic monomials of this dimension
    //genmodk(numv, blockdeg, mons[blockdeg],1);
    if(blockdeg>polymaxdeg){//empty blocks
      for(j=0;j<jc;j++){
	blks[tot]=(int *)malloc((size_t) (1*sizeof(int)));
	blkinds[tot]=NULL;
	blks[tot++][0]=0;
	if(debugfull){fprintf(stderr, "dimension: %d, block size: 0\n", blockdeg);}
      }
    }
    else{
      blockmindeg=ceiling(polymindeg-blockdeg, 2);
      blockmaxdeg=(polymaxdeg-blockdeg)/2;
      if(blockmaxdeg<blockmindeg){//empty blocks
	for(j=0;j<jc;j++){
	  blks[tot]=(int *)malloc((size_t) (1*sizeof(int)));
	  blkinds[tot]=NULL;
	  blks[tot++][0]=0;
	  if(debugfull){fprintf(stderr, "dimension: %d, block size: 0\n", blockdeg);}
	}
      }
      else{// nonempty block
	//set offset according to blockmindeg
	//This is to reference the monomials in sqmons
	offset=0;
	for(j=0;j<blockmindeg;j++)
	  offset+=choose(numv+j-1,j);
	numsqmon=0;
	for(j=blockmindeg;j<=blockmaxdeg;j++)//number of monomials
	  numsqmon+=choose(numv+j-1,j);
	if(numsqmon<=1){//ignore these blocks: set their size to zero as they will appear again
	  for(j=0;j<jc;j++){
	    blks[tot]=(int *)malloc((size_t) (1*sizeof(int)));
	    blkinds[tot]=NULL;
	    blks[tot++][0]=0;
	    if(debug){fprintf(stderr, "dimension: %d, block size: 0\n", blockdeg);}
	  }
	}
	else{
	  if(debugfull){fprintf(stderr, "******\nblockdeg=%d\n******\n", blockdeg);}
	  for(j=0;j<jc;j++){// each block
	    blks[tot]=(int *)malloc((size_t) (((numsqmon*(numsqmon+1))/2+1)*sizeof(int)));
	    blks[tot][0]=numsqmon;
	    blkinds[tot]=(int *)malloc((size_t) (blks[tot][0]*sizeof(int)));
	    if(debugfull){fprintf(stderr, "old block %d, size %d\n", tot, numsqmon);}
	    for(k=0;k<numsqmon;k++){
	      for(l=k;l<numsqmon;l++){
		vecsum3(newmon,mons[tot],sqmons[k+offset],sqmons[l+offset],numv);
		//printvec(newmon,numv);
		totmons=addoneintstra(totmons, newmon, numv, &allmons, &ind);
		(monrefs[ind][0])++;
		monrefs[ind]=(int *)realloc(monrefs[ind], (size_t)((3*monrefs[ind][0]+1)*sizeof(int)));
		monrefs[ind][3*monrefs[ind][0]-2]=tot;monrefs[ind][3*monrefs[ind][0]-1]=k;monrefs[ind][3*monrefs[ind][0]]=l;
		if(l!=k)
		  (nondiag[ind])++;
		blks[tot][symmattovec(k,l,numsqmon)+1]=ind;
		//fprintf(stderr, "%d ", ind);
	      }
	    }
	    tot++;
	    if(debugfull){fprintf(stderr, "dimension: %d, block size: %d\n", blockdeg, (numsqmon*(numsqmon+1))/2);}
	  }
	}
      }
    }
  }



  if(totmons!=numposmons){
    fprintf(stderr, "something has gone wrong. totmons=%d, numposmons=%d, EXITING.\n", totmons, numposmons);
    exit(0);
  }

  //Blocks with a single monomial
  for(j=totblk;j<totblk+numposmons;j++){
    blks[j]=(int *)malloc((size_t)(2*sizeof(int)));
    blks[j][0]=1;//block size
    blks[j][1]=j-totblk;//monomial index
    blkinds[j]=(int *)malloc((size_t) (blks[j][0]*sizeof(int)));
    //monomial references already done right at the start
  }

  //
  // If a monomial figures only as a diagonal monomial of RHS blocks, then we can remove it
  // Do this until there is no change
  //

  q=0;
  flag=0;
  if(preprocess)
    flag=1;
  while(flag){
    flag=0;
    // what to do with redundant monomials
    if(debugfull){fprintf(stderr, "q=%d\n", q);}
    for(j=0;j<totmons;j++){//Each monomial
      //!(nondiag[j]) means only figures on diagonals and figures at least once
      if(!(nondiag[j]) && allmons[j][0]){//only diagonal entries in blocks: set whole row-col to zero
	if(debugfull){fprintf(stderr, "Current monomial which doesn't figure off-diagonal:\n");
	  printvec(allmons[j]+2, numv);}

	for(j1=0;j1<monrefs[j][0];j1++){//total references to this monomial
	  //monrefs=(total, block, row, col, block, row, col...)
	  if(monrefs[j][3*j1+2]==monrefs[j][3*j1+3]){//diagonal references to this monomial (off-diagonals must already have been removed)
	    blkind=monrefs[j][3*j1+1];
	    rowind=monrefs[j][3*j1+2];
	    if(debugfull){fprintf(stderr, "*****\nRowcol to remove (%d in block %d; size %d):\n", rowind, blkind, blks[blkind][0]);}
	    //printvec(allmons[j]+2, numv);

	    //remove row and column
	    for(k=0;k<blks[blkind][0];k++){
	      for(l=k;l<blks[blkind][0];l++){
		monind=blks[blkind][symmattovec(k,l,blks[blkind][0])+1];//index of monomial (or -1 if already removed)
		//fprintf(stderr, "monind=%d\n",monind);
		if((k==rowind || l==rowind) && monind!=-1){//not already removed
		  flag=1;//something got removed in this iteration
		  if(debugfull){fprintf(stderr, "removing...(%d, %d) in block %d (monomial %d)\n", k, l, blkind, monind);
		    printvec(allmons[monind]+2, numv);
		    fprintf(stderr, "%d left\n", allmons[monind][0]-1);}
		  //total occurrence of monomial decreases
		  (allmons[monind][0])--;
		  //checks: 
		  if(allmons[monind][0]==0){
		    if(debugfull){fprintf(stderr, "removed all of: ");
		      printvec(allmons[monind]+2, numv);}
		    if(isinarray(lst2, LHtot, allmons[monind]+2, numv)>=0){
		      fprintf(stderr, "Something went wrong: cannot remove a monomial which figures on LHS.\n");exit(0);
		    }
		  }
		  else if(allmons[monind][0]==1 && isinarray(lst2, LHtot, allmons[monind]+2, numv)>=0){
		    fprintf(stderr, "Something went wrong: A monomial which figures on LHS can no longer be cancelled.\n");
		    printvec(allmons[monind]+2, numv);exit(0);
		  }

		  //Indicates that this entry in the block is now empty
		  blks[blkind][symmattovec(k,l,blks[blkind][0])+1]=-1;
		  //Removed a nondiagonal occurrence of the monomial
		  if(l!=k){
		    /* fprintf(stderr, "removing nondiagonal entry...(%d, %d) in block %d (monomial %d)\n", k, l, blkind, monind); */
		    /* printvec(allmons[monind]+2, numv); */
		    (nondiag[monind])--;
		  }
		}
	      }
	    }
	  }
	}
      }

    }

    q++;
  }

  //
  // Enter the main loop
  //

  //
  // A good indicator of whether a solution is artefact is whether it survives with a small positive tol below
  //

  p=0;
  numnewconstraints=0;
  tol=0.005;
  while(p<10){
    //
    // Update the blocks
    //
    for(j=0;j<totblk+numposmons;j++){//each block
      newblksz=0;//what's left in block
      if(blks[j][0]){//potential nonempty block
	//to hold indices of surviving square monomials
	if(debugfull){fprintf(stderr, "\nold block %d, size %d:\n", j, blks[j][0]);}
	for(k=0;k<blks[j][0];k++){
	  if(blks[j][symmattovec(k,k,blks[j][0])+1]!=-1){
	    if(debugfull){
	      fprintf(stderr, "non-empty row-col: block %d, rowcol %d\n", j, k+1);
	      printvec(allmons[blks[j][symmattovec(k,k,blks[j][0])+1]]+2, numv);
	    }
	    blkinds[j][newblksz++]=k;
	  }
	}
	if(newblksz==0){
	  if(debugfull){fprintf(stderr, "Newly empty block: %d\n", j);}
	}
	else if(j<totblk && newblksz==1){//now a redundant (repeated) block
	  if(debugfull){fprintf(stderr, "Block reduced to size 1 (remove) %d\n", j);}
	  (allmons[blks[j][symmattovec(blkinds[j][0],blkinds[j][0],blks[j][0])+1]][0])--;
	  newblksz=0;
	}

	if(debugfull){fprintf(stderr, "new block %d, size %d:\n", j, newblksz);}
	for(k=0;k<newblksz;k++){
	  for(l=k;l<newblksz;l++){
	    //Should be safe: only refers to later monomials
	    blks[j][symmattovec(k,l,newblksz)+1]=blks[j][symmattovec(blkinds[j][k],blkinds[j][l],blks[j][0])+1];
	    if(debugfull){fprintf(stderr, "%d ", blks[j][symmattovec(k,l,newblksz)+1]);
	    }}
	}
	blks[j][0]=newblksz;
      }
    }

    //
    // The surviving monomials
    // One constraint for each surviving monomial
    // Structure of a constraint: total memory allocation; index of monomial; 
    // then blocks of three (block number, row-index, column-index)
    //

    numconstraints=0;
    for(j=0;j<totmons;j++){
      if(allmons[j][0]){//surviving monomials
	newconstraint(&constraints, &numconstraints, j);
	allmons[j][1]=numconstraints-1;//constraint index

	if(debugfull){
	  if(isinarray(lst2, LHtot, allmons[j]+2, numv)>=0)
	    fprintf(stderr, "In LH poly: ");      
	  fprintf(stderr, "monomial to constraint: %d to %d\n", j, numconstraints-1);
	  printvec(allmons[j]+2, numv);}
      }
      else{
	if(debugfull){fprintf(stderr, "monomial deleted: ");
	  printvec(allmons[j]+2,numv);}
      }
    }

    //
    // Output block data and populate the constraints
    //
    numtrueblk=0;
    for(j=0;j<totblk+totmons;j++){
      if(j==totblk)
	numsqblk=numtrueblk;//Number of blocks other than the single entry ones

      if(blks[j][0]){
	if(j<totblk){
	  if(debugfull){fprintf(stderr, "block %d, size %d, corresponding to basic monomial:\n", j, blks[j][0]);
	    printvec(mons[j],numv);}
	}
	else{
	  if(debugfull){fprintf(stderr, "1X1 block %d corresponding to monomial:\n", j);
	    printvec(allmons[j-totblk]+2,numv);}
	}
	/* fprintf(stderr, "square monomials:\n"); */
	/* for(k=0;k<blks[j][0];k++) */
	/* 	printvec(sqmons[blkinds[j][k]],numv); */
	//printvec(blks[j]+1, (blks[j][0]*(blks[j][0]+1))/2);
	for(k=0;k<blks[j][0];k++){
	  for(l=k;l<blks[j][0];l++){
	    ind=allmons[blks[j][symmattovec(k,l,blks[j][0])+1]][1];//index of constraint
	    growconstraint(constraints, ind, numtrueblk, k, l);
	  }
	}
	numtrueblk++;
      }

    }

    p++;
    if(debugfull){fprintf(stderr, "newconstraints=%d\n", numnewconstraints);}
    //
    //Start writing tempfiles/csdp.dat
    //

    if(!(fd=fopen("tempfiles/csdp.dat", "w"))){
      fprintf(stderr, "ERROR: \"tempfiles/csdp.dat\" could not be opened for reading.\n");
      exit(0);
    }

    //Number of constraints
    fprintf(fd, "%d\n", numconstraints+numnewconstraints+1);

    //Number of blocks
    if(ldegree){
      fprintf(fd, "%d\n", numtrueblk+totlmons+1);
      trueblksz=(int*)malloc((size_t) ((numtrueblk+totlmons+1)*sizeof(int)));
    }
    else{
      fprintf(fd, "%d\n", numtrueblk+1);
      trueblksz=(int*)malloc((size_t) ((numtrueblk+1)*sizeof(int)));
    }

    //Block sizes (stored in trueblksz)
    tot=0;
    for(j=0;j<totblk+totmons;j++){
      if(blks[j][0]){
	fprintf(fd, "%d ", blks[j][0]);
	trueblksz[tot++]=blks[j][0];
      }
    }
    if(ldegree){//LHS parameter blocks
      if(posleft){
	for(j=0;j<totlmons;j++){
	  fprintf(fd, "1 ");
	  trueblksz[tot++]=1;
	}
      }
      else{
	for(j=0;j<totlmons;j++){
	  fprintf(fd, "2 ");
	  trueblksz[tot++]=2;
	}
      }
    }
    fprintf(fd, "1\n");//final block
    trueblksz[tot++]=1;

    //Constraint RHSs
    for(j=0;j<numconstraints;j++){//surviving monomials
      if(ldegree<=0 && (s=isinarray(lst, r, allmons[constraints[j][1]]+2, numv))>=0)
	fprintf(fd, "%.10f ", cflst[s]);
      else
	fprintf(fd, "0.0 ");
    }
    for(j=0;j<numnewconstraints;j++)
      fprintf(fd, "%.1f ", (double)(newconstraints[j][1])/1.0);//change if using nearhalfint

    //final constraint RHS
    //This seems to matter: if the poly has large coefficients, then setting this
    //small can make the algorithm fail (roundoff error?)
    //setting this too large can make a degree 0 strict problem fail
    //current heuristic: set it to 1% of largest coeff
    maxcf=0.0;
    for(s=0;s<r;s++){
      if(maxcf<fabs(cflst[s]))
	maxcf=fabs(cflst[s]);
    }
    if(strict && !ldegree)
      fprintf(fd, "%.4f\n",0.01*maxcf);//very heuristic
    else if(strict || ldegree)
      fprintf(fd, "1\n");
    else
      fprintf(fd, "0\n");

    //Objective function

    if(ldegree){
      if(!posleft){
	for(j=0;j<totlmons;j++)
	  fprintf(fd, "0 %d 1 1 -1.0\n0 %d 2 2 -1.0\n", numtrueblk+j+1, numtrueblk+j+1);
      }
      else{
	for(j=0;j<totlmons;j++)
	  fprintf(fd, "0 %d 1 1 -1.0\n", numtrueblk+j+1);
      }
    }
    else{// diagonal
      if(p>0){
	tot=0;
	/* for(j=0;j<totblk;j++){ */
	/*   if(blks[j][0]){ */
	/*     tot++; */
	/*     if(tot==2 || tot==4 || tot>6){ */
	/*       for(k=0;k<blks[j][0];k++) */
	/*       	for(l=k;l<blks[j][0];l++) */
	/*       	  if(l==k) */
	/*       	    fprintf(fd, "0 %d %d %d -1.0\n", tot, k+1, l+1); */
	/*       	  else */
	/*       	    fprintf(fd, "0 %d %d %d -0.5\n", tot, k+1, l+1); */
	/*     } */
	/*     else if(tot==1) */
	/*       for(k=9;k<blks[j][0];k++) */
	/* 	fprintf(fd, "0 %d %d %d 1.0\n", tot, k+1, k+1); */

	/*   } */
	/* } */
	/* for(j=totblk;j<totblk+numposmons;j++){ */
	/* 	if(blks[j][0]){ */
	/* 	  tot++; */
	/* 	  for(k=0;k<blks[j][0];k++) */
	/* 	    fprintf(fd, "0 %d %d %d 1.0\n", tot, k+1, k+1); */
	/* 	} */
	/* } */
	}
    }

    //constraints
    for(j=0;j<numconstraints;j++){
      if(constraints[j][1]>=0){
	if(debugfull){fprintf(stderr, "Constraint on monomial %d: \n", j);
	  printvec(allmons[constraints[j][1]]+2, numv);}
      }
      //    fprintf(stderr, "lines in constraint: %d\n", (constraints[j][0]-2)/3);
      for(k=0;k<(constraints[j][0]-2)/3;k++){
	fprintf(fd, "%d %d %d %d 1.0\n", j+1, constraints[j][3*k+2], constraints[j][3*k+3], constraints[j][3*k+4]);
      }
      if(ldegree>0 && (s=isinarray(lst2, LHtot, allmons[constraints[j][1]]+2, numv))>=0){//LHS parameters
	for(k=0;k<totlmons;k++){
	  if(lst1[s][k]){
	    if(posleft)
	      fprintf(fd, "%d %d 1 1 %.4f\n", j+1, numtrueblk+k+1, -lst1[s][k]);
	    else
	      fprintf(fd, "%d %d 1 2 %.4f\n", j+1, numtrueblk+k+1, -lst1[s][k]/2.0);
	    //fprintf(stderr, "printing %d --> %.1f\n", lst1[s][k], -((double)(lst1[s][k]))/2.0);  
	    //printvec(lst1[s],totlmons);
	  }
	}
      }
    }
    for(j=0;j<numnewconstraints;j++){
      for(k=0;k<(newconstraints[j][0]-2)/3;k++){
	fprintf(fd, "%d %d %d %d 1.0\n", numconstraints+j+1, newconstraints[j][3*k+2], newconstraints[j][3*k+3], newconstraints[j][3*k+4]);
      }
    }

    //final constraint

    if(ldegree && !strict){
      numallblk=numtrueblk+totlmons+1;
      if(posleft){
	for(k=numtrueblk;k<numtrueblk+totlmons;k++)
	  fprintf(fd, "%d %d 1 1 1.0\n", numconstraints+j+1, k+1);
      }
      else{
	for(k=numtrueblk;k<numtrueblk+totlmons;k++)
	  fprintf(fd, "%d %d 1 2 1.0\n", numconstraints+j+1, k+1);
      }
      fprintf(fd, "%d %d 1 1 -1.0\n", numconstraints+j+1, k+1);
    }
    else{
      numallblk=numsqblk;
      for(k=0;k<totmons;k++){
	if(blks[k+totblk][0]){//nonempty
	  numallblk++;
	  fprintf(fd, "%d %d 1 1 1.0\n", numconstraints+j+1, numallblk);
	}
      }
      if(ldegree)//skip the blocks associated with LHS parameters
	numallblk+=totlmons+1;
      else
	numallblk++;
      fprintf(fd, "%d %d 1 1 -1.0\n", numconstraints+j+1, numallblk);
    }

    fclose(fd);

    if(debug){fprintf(stderr, "numallblk=%d\n", numallblk);}

//    if(ldegree==1 && p==2)
//      exit(0);
    //
    //Finished writing tempfiles/csdp.dat, now run
    //

    if(!(fp=popen("csdp tempfiles/csdp.dat tempfiles/tmp.sol", "r"))){
      perror("couldn't run command \"csdp tempfiles/csdp.dat tempfiles/tmp.sol\". EXITING\n");exit(0);
    }

    retstat=-2;
    while(fgets(path, sizeof(path), fp)){
      if(debug){fprintf(stderr, "%s", path);}
      if(strstr(path, "Declaring primal infeasibility"))
	retstat=-1;
      else if(strstr(path, "Partial Success: SDP solved with reduced accuracy"))
	retstat=2;
      else if(strstr(path, "Success: SDP solved"))
	retstat=0;
    }
    pclose(fp);
    if(retstat==-1){
      free((char*)trueblksz);
      for(j=0;j<numconstraints;j++)
	free((char*)constraints[j]);
      free((char*)constraints);
      break;
    }

    //Success or partial success
    outstr=readfileintostr("tempfiles/tmp.sol");
    datstr=strstr(outstr,"\n2");
    datstr++;
    //  fprintf(stderr, "datstr=%s\n", datstr);
    outdat=genfloatvec(datstr, &datlen);
    if(debugfull){fprintf(stderr, "datlen=%d\n", datlen);}
    if(datlen%5!=0){
      fprintf(stderr, "ERROR reading output: data not a multiple of 5\n");
      exit(0);
    }
    // To hold the data
    blksout=(double **)malloc((size_t) (numallblk*sizeof(double*)));
    for(j=0;j<numallblk;j++){
      blksz=trueblksz[j];
      if(debugfull){fprintf(stderr, "j=%d, blksz=%d\n", j, blksz);}
      blksout[j]=(double *)malloc((size_t) (((blksz*(blksz+1))/2+1)*sizeof(double*)));
      blksout[j][0]=blksz;
      for(k=1;k<(blksz*(blksz+1))/2+1;k++)//default is zero
	blksout[j][k]=0.0;
    }

    j=0;
    maxval=0.0;
    if(debugfull){fprintf(stderr, "Output data:\n");}
    while(j<datlen){//reprint output data
      if(debugfull){fprintf(stderr, "%d %d %d %d %f\n", (int)outdat[j], (int)outdat[j+1], (int)outdat[j+2], (int)outdat[j+3], outdat[j+4]);}
      maxval=fabs(outdat[j+4])>maxval?fabs(outdat[j+4]):maxval;
      j+=5;
    }
    if(debugfull){fprintf(stderr, "$$$maxval=%.4f\n$$$\n", maxval);}

    j=0;
    while(j<datlen){
      blksout[(int)(outdat[j+1])-1][symmattovec((int)(outdat[j+2])-1,(int)(outdat[j+3])-1,blksout[(int)(outdat[j+1])-1][0])+1]=outdat[j+4];
      j+=5;
    }


    numallblk=0;outpolylen=0;
    for(j=0;j<totblk+totmons;j++){
      if(blks[j][0]){
	tot=1;
	outblkpolylen=0;
	for(k=0;k<trueblksz[numallblk];k++){
	  for(l=k;l<trueblksz[numallblk];l++){
	    if(l==k){
	      outpolylen=addpolyterm(outpolylen, allmons[blks[j][tot]]+2, blksout[numallblk][tot], numv, &outpoly, &outcfs, &ind);
	      outblkpolylen=addpolyterm(outblkpolylen, allmons[blks[j][tot]]+2, blksout[numallblk][tot], numv, &outblkpoly, &outblkcfs, &ind);
	    }
	    else{
	      outpolylen=addpolyterm(outpolylen, allmons[blks[j][tot]]+2, 2.0*blksout[numallblk][tot], numv, &outpoly, &outcfs, &ind);
	      outblkpolylen=addpolyterm(outblkpolylen, allmons[blks[j][tot]]+2, 2.0*blksout[numallblk][tot], numv, &outblkpoly, &outblkcfs, &ind);
	    }
	    //fprintf(stderr, "%.4f: ", blksout[numallblk][tot]); printvec(allmons[blks[j][tot]]+2,numv);
	    tot++;
	  }
	}//finished dealing with this block

	if(debug){
	  if(j<totblk || fabs(outblkcfs[0])>0.001){
	    fprintf(stderr, "########true block %d\n", numallblk);
	    cerr << lrexptopolysimp(outblkpoly, outblkcfs, 0.001, outblkpolylen, numv) << endl;
	  }
	  for(k=0;k<outblkpolylen;k++){
	    if(fabs(outblkcfs[k])>0.001){
	      fprintf(stderr, "%.3f: ", outblkcfs[k]);
	      printvec(outblkpoly[k],numv);
	    }
	  }
	  if(j<totblk || fabs(outblkcfs[0])>0.001)
	    fprintf(stderr, "########\n");
	}
	ifree(outblkpoly,outblkpolylen);
	free((char*)outblkcfs);

	numallblk++;
      }
    }

    //RHS poly
    if(debug){
      fprintf(stderr, "\n****RHS poly****\n");
      cerr << lrexptopolysimp(outpoly, outcfs, 0.001, outpolylen, numv) << endl;
      for(j=0;j<outpolylen;j++){
	fprintf(stderr, "%.2f: ", outcfs[j]);
	printvec(outpoly[j],numv);
      }
      fprintf(stderr, "\n*********\n");
    }

    if(ldegree){
      outblkpolylen=0;
      for(s=0;s<LHtot;s++){
	if(posleft){
	  for(k=0;k<totlmons;k++)
	    outpolylen=addpolyterm(outpolylen, lst2[s], -1.0*blksout[numallblk+k][1]*lst1[s][k], numv, &outpoly, &outcfs, &ind);
	}
	else{
	  for(k=0;k<totlmons;k++)
	    outpolylen=addpolyterm(outpolylen, lst2[s], -1.0*blksout[numallblk+k][2]*lst1[s][k], numv, &outpoly, &outcfs, &ind);
	}
      }
      if(posleft){
	for(k=0;k<totlmons;k++)
	  outblkpolylen=addpolyterm(outblkpolylen, lmons[k], blksout[numallblk+k][1], numv, &outblkpoly, &outblkcfs, &ind);
      }
      else{
	for(k=0;k<totlmons;k++)
	  outblkpolylen=addpolyterm(outblkpolylen, lmons[k], blksout[numallblk+k][2], numv, &outblkpoly, &outblkcfs, &ind);
      }

      if(debug){
	fprintf(stderr, "########LHS poly\n");
	cerr << lrexptopolysimp(outblkpoly, outblkcfs, 0.001, outblkpolylen, numv) << endl;
	for(k=0;k<outblkpolylen;k++){
	  if(fabs(outblkcfs[k])>0.001){
	    fprintf(stderr, "%.4f: ", outblkcfs[k]);
	    printvec(outblkpoly[k],numv);
	  }
	}
	fprintf(stderr, "########\n");
      }
      ifree(outblkpoly,outblkpolylen);
      free((char*)outblkcfs);
    }
    else{
      for(s=0;s<r;s++)
	outpolylen=addpolyterm(outpolylen, lst[s], (double)(-cflst[s]), numv, &outpoly, &outcfs, &ind);
    }

    if(debug){
      fprintf(stderr, "Final polynomial (to 4 d.p.):\n");
      cerr << lrexptopolysimp(outpoly, outcfs, 0.0001, outpolylen, numv) << endl;
      for(j=0;j<outpolylen;j++){
	if(fabs(outcfs[j])>0.0001){
	  fprintf(stderr, "%.4f: ", outcfs[j]);
	  printvec(outpoly[j],numv);
	}
      }
      fprintf(stderr, "\n***********************************\n\n");
    }

    if(numnewconstraints>0){
      for(j=0;j<numnewconstraints;j++)
	free((char*)newconstraints[j]);
      free((char*)newconstraints);
    }
    numnewconstraints=0;


    // remove rowcols which seem to be zero
    numallblk=0;flag=0;
    for(j=0;j<totblk+numposmons;j++){
      if(blks[j][0]){
	for(k=0;k<trueblksz[numallblk];k++){
	  if(nearint(blksout[numallblk][symmattovec(k,k,blks[j][0])+1], 0, 0, tol)!=-1000){
	    for(l=0;l<trueblksz[numallblk];l++){//rowcol to zero: start at l=0
	      if(blks[j][symmattovec(k,l,blks[j][0])+1]!=-1){//not already set to zero by some other deletion
		flag=1;
		if(debugfull){fprintf(stderr, "setting block %d, entry (%d,%d) = %d, %.5f to zero.\n", j, k,l, symmattovec(k,l,blks[j][0])+1, blksout[numallblk][symmattovec(k,k,blks[j][0])+1]);}
		(allmons[blks[j][symmattovec(k,l,blks[j][0])+1]][0])--;
		blks[j][symmattovec(k,l,blks[j][0])+1]=-1;
	      }
	    }//end of rowcol to zero
	  }
	}
	numallblk++;
      }
    }

    //Any values to be manually fixed? E.g. 
    /* addnewconstraint(&numnewconstraints, &newconstraints, 1, 16, 16, 9); */

    //End of rowcol removal


    ifree(outpoly,outpolylen); 
    free((char*)outcfs);
    free((char*)trueblksz);
    for(j=0;j<numallblk;j++)
      free((char*)blksout[j]);
    free((char*)blksout);
    free((char*)outdat);
    free(outstr);

    for(j=0;j<numconstraints;j++)
      free((char*)constraints[j]);
    free((char*)constraints);

    if(!flag)//nothing changed on this iteration
      break;

  }

  if(numnewconstraints>0){
    for(j=0;j<numnewconstraints;j++)
      free((char*)newconstraints[j]);
    free((char*)newconstraints);
  }

  ifree(allmons,totmons);
  ifree(blks,totblk+numposmons);
  for(j=0;j<totblk+numposmons;j++)
    if(blkinds[j])
      free ((char *)(blkinds[j]));
  free((char *)blkinds);
  free_imatrix(mons,0,totblk-1,0,numv-1);
  free_imatrix(sqmons,0,choose(t,numv)-1,0,numv-1);
  ifree(lst1,LHtot);ifree(lst2,LHtot);
  ifree(monrefs,numposmons);
  free_imatrix(lmons, 0, totlmons-1,0,numv-1);
  free_imatrix(posmons,0,numposmons-1,0,numv-1);
  free((char*)nondiag);

  return retstat;
}

/////////////////

//SDP return codes are the same as for ppp1. 0: feasible; -1: infeasible; 2: probably feasible (reduced accuracy); -2: something went wrong
int pppfull1(int **lst, double *cflst, long r, int numv, int polymindeg, int polymaxdeg, int preprocess, int ldegreemin, int ldegree, int strict, int debug){
  int p,q,j,j1,k,l,blockmindeg,blockmaxdeg,flag;
  // Only the trivial monomial
  int totblk=1;
  int tot=0;//counts the blocks
  long s;
  int **lst2;
  double **lst1;
  int LHtot;
  int **blks;
  double **blksout;
  int **blkinds;
  int totmons=0;
  int **allmons=NULL;
  int newmon[numv];
  //monomials which can figure in squares
  int sqmaxdeg;
  int **sqmons, **lmons;
  int offset;
  int numsqmon;
  int ind;
  int blksz, newblksz;
  int **constraints;
  int numconstraints=0;
  int **newconstraints=NULL;
  int numnewconstraints=0;
  int numtrueblk,numsqblk,numallblk;
  int *trueblksz;
  FILE *fd, *fp;
  int totlmons;
  int **posmons;
  int numposmons;
  int **possqmons;
  int numpossqmons;
  int *nondiag;
  int **monrefs;
  int blkind, rowind, monind;
  char *outstr, *datstr;
  double *outdat;
  int datlen;
  int outpolylen=0;
  int **outpoly;
  double *outcfs;
  int outblkpolylen;
  int **outblkpoly;
  double *outblkcfs;
  double tol=0.0;
  double maxval;
  int posleft=0;//1 means LHS multiplier has only positive coefficients
  int retstat=-2;
  char path[2048];
  double maxcf;
  int debugfull=(debug<=0)?0:debug-1;

  if(polymaxdeg%2!=0 || polymindeg%2!=0){
    fprintf(stderr, "Error in pppfull1: maximum and minimum degrees need to be even.\n");
    exit(0);
  }

  if(polymindeg>polymaxdeg){
    fprintf(stderr, "Error in pppfull1: maximum degree must be >= minimum degree.\n");
    exit(0);
  }


  if(ldegreemin%2!=0 || ldegree%2!=0){
    fprintf(stderr, "Error in pppfull1: maximum and minimum degrees of LHS need to be even.\n");
    exit(0);
  }

  if(ldegreemin>ldegree){
    fprintf(stderr, "Error in pppfull1: maximum ldegree must be >= minimum ldegree.\n");
    exit(0);
  }


  if(!strict)//otherwise possible artifacts
    posleft=1;

  /* fprintf(stderr, "%d\n", NewtonPolyVertSgn(tmp)); */
  /* exit(0); */
  
  //monomials on LHS (all monomials between degree ldegreemin and ldegree)

  lmons=genmons1(numv,ldegreemin,ldegree,&totlmons);

  // The LHS poly is the product of original poly P and parameterised poly Q of degrees up to ldegree
  // lst1 stores the coeffs of the LHS poly which are themselves monomials in the additional LHS variables
  // lst2 stores the monomials of the LHS poly
  LHtot=polyprod(lmons, totlmons, cflst, lst, r, numv, &lst1, &lst2);

  //increase the maximum degree and minimum degrees of the polynomial
  polymaxdeg+=ldegree;
  polymindeg+=ldegreemin;
  if(debugfull){fprintf(stderr, "mindeg=%d, maxdeg=%d\n", polymindeg, polymaxdeg);}

  // posmons stores all monomials involving numv variables
  // in the range between the minimum and maximum degrees of the 
  // LHS polynomial. These are all the possible monomials which could 
  // figure on the RHS. This does not mean that they will.
  posmons=genmons1(numv, polymindeg, polymaxdeg, &numposmons);

  //Just the squares
  //possqmons=gensqmons(numv, polymindeg, polymaxdeg, &numpossqmons);
  if(strict)
    possqmons=genevenpows(numv, polymindeg, polymaxdeg, &numpossqmons);
  else
    possqmons=genemptymonslist(numv, polymindeg, polymaxdeg, &numpossqmons);
  //Vector to hold the number of times that a monomial figures
  //but not as a diagonal entry in a block
  nondiag=(int *)malloc((size_t) (numposmons*sizeof(int)));
  inittozero(nondiag,numposmons);

  // Array of vectors of length 3k+1 which tells us where each monomial
  // occurs in blocks on the RHS (doesn't list occurrence on the LHS)
  monrefs=(int **)malloc((size_t) (numposmons*sizeof(int*)));

  for(s=0;s<numposmons;s++){// add all the monomials in range (this fixes the monomial order)
    totmons=addoneintstra(totmons, posmons[s], numv, &allmons,&ind);
    monrefs[ind]=(int *)malloc((size_t) (sizeof(int)));
    //total occurrences, block number, then indices.
    monrefs[ind][0]=0;
  }

  for(s=0;s<numpossqmons;s++){// add all the square monomials and reference
    totmons=addoneintstra(totmons, possqmons[s], numv, &allmons,&ind);
    (monrefs[ind][0])++;
    monrefs[ind]=(int *)realloc(monrefs[ind], (size_t)((3*monrefs[ind][0]+1)*sizeof(int)));
    //total occurrences, block number, then indices.
    monrefs[ind][3*monrefs[ind][0]-2]=totblk+s;monrefs[ind][3*monrefs[ind][0]-1]=0;monrefs[ind][3*monrefs[ind][0]]=0;
  }

  if(debugfull){fprintf(stderr, "LHS polynomial:\n");}
  for(s=0;s<LHtot;s++){// add the monomials from the LHS poly, but don't reference (i.e., not in a block)
    totmons=addoneintstra(totmons, lst2[s], numv, &allmons,&ind);
    (nondiag[ind])++;//LHS monomials can't be removed from blocks
    if(debugfull){
      printvec(lst1[s],totlmons);
      printvec(lst2[s],numv);
    }
  }

  //From now on the monomials in allmons are in the same order as those in posmons
  //One block for each basic monomial (i.e., 2^numv monomials where powers are 0 and 1)
  //+ one block for each allowed monomial

  blks=(int **)malloc((size_t) ((totblk+numpossqmons)*sizeof(int*)));
  blkinds=(int **)malloc((size_t) ((totblk+numpossqmons)*sizeof(int*)));//stuff that will survive

  // Overkill in generating the monomials which will figure in the squares (unnecessary low degrees)
  sqmaxdeg=polymaxdeg/2;
  sqmons=imatrix(0, choose(numv+sqmaxdeg,numv)-1,0,numv-1);
  genmons(numv,sqmaxdeg,sqmons);

  //
  // Initial RHS block
  //

  blockmindeg=polymindeg/2;
  blockmaxdeg=polymaxdeg/2;
  //set offset according to blockmindeg
  //This is to reference the monomials in sqmons
  offset=0;
  for(j=0;j<blockmindeg;j++)
    offset+=choose(numv+j-1,j);
  numsqmon=0;
  for(j=blockmindeg;j<=blockmaxdeg;j++)//number of monomials
    numsqmon+=choose(numv+j-1,j);

  if(numsqmon<=1){//ignore these blocks: set their size to zero as they will appear again
    blks[tot]=(int *)malloc((size_t) (1*sizeof(int)));
    blkinds[tot]=NULL;
    blks[tot++][0]=0;
    if(debugfull){fprintf(stderr, "dimension: 0, block size: 0\n");}
  }
  else{
    blks[tot]=(int *)malloc((size_t) (((numsqmon*(numsqmon+1))/2+1)*sizeof(int)));
    blks[tot][0]=numsqmon;
    blkinds[tot]=(int *)malloc((size_t) (blks[tot][0]*sizeof(int)));
    if(debug){fprintf(stderr, "old block %d, size %d\n", tot, numsqmon);}
    for(k=0;k<numsqmon;k++){
      for(l=k;l<numsqmon;l++){
	vecsum2(newmon,sqmons[k+offset],sqmons[l+offset],numv);
	//printvec(newmon,numv);
	totmons=addoneintstra(totmons, newmon, numv, &allmons, &ind);
	(monrefs[ind][0])++;
	monrefs[ind]=(int *)realloc(monrefs[ind], (size_t)((3*monrefs[ind][0]+1)*sizeof(int)));
	monrefs[ind][3*monrefs[ind][0]-2]=tot;monrefs[ind][3*monrefs[ind][0]-1]=k;monrefs[ind][3*monrefs[ind][0]]=l;
	if(l!=k)
	  (nondiag[ind])++;
	blks[tot][symmattovec(k,l,numsqmon)+1]=ind;
	//fprintf(stderr, "%d ", ind);
      }
    }
    tot++;
    if(debug){fprintf(stderr, "dimension: 0, block size: %d\n", (numsqmon*(numsqmon+1))/2);}
  }


  if(totmons!=numposmons){
    fprintf(stderr, "something has gone wrong in pppfull1. totmons=%d, numposmons=%d, EXITING.\n", totmons, numposmons);
    exit(0);
  }


  //Blocks with a single monomial
  for(j=totblk;j<totblk+numpossqmons;j++){
    blks[j]=(int *)malloc((size_t)(2*sizeof(int)));
    blks[j][0]=1;//block size

    totmons=addoneintstra(totmons, possqmons[j-totblk], numv, &allmons,&ind);

    blks[j][1]=ind;//monomial index
    blkinds[j]=(int *)malloc((size_t) (blks[j][0]*sizeof(int)));
    //monomial references already done right at the start
  }

  //
  // If a monomial figures only as a diagonal monomial of RHS blocks, then we can remove it
  // Do this until there is no change
  //

  q=0;
  flag=0;
  if(preprocess)
    flag=1;
  while(flag){
    flag=0;
    // what to do with redundant monomials
    if(debug){fprintf(stderr, "q=%d\n", q);}
    for(j=0;j<totmons;j++){//Each monomial
      //!(nondiag[j]) means only figures on diagonals and figures at least once
      if(!(nondiag[j]) && allmons[j][0]){//only diagonal entries in blocks: set whole row-col to zero
	if(debug){fprintf(stderr, "Current monomial which doesn't figure off-diagonal:\n");
	  printvec(allmons[j]+2, numv);}

	if(debug){fprintf(stderr, "total references = %d\n", monrefs[j][0]);
	  printvec(monrefs[j], 3*monrefs[j][0]+1);}
	for(j1=0;j1<monrefs[j][0];j1++){//total references to this monomial
	  //monrefs=(total, block, row, col, block, row, col...)
	  if(monrefs[j][3*j1+2]==monrefs[j][3*j1+3]){//diagonal references to this monomial (off-diagonals must already have been removed)
	    blkind=monrefs[j][3*j1+1];
	    rowind=monrefs[j][3*j1+2];
	    if(debugfull){fprintf(stderr, "*****\nRowcol to remove (%d in block %d; size %d):\n", rowind, blkind, blks[blkind][0]);}
	    //printvec(allmons[j]+2, numv);

	    //remove row and column
	    for(k=0;k<blks[blkind][0];k++){
	      for(l=k;l<blks[blkind][0];l++){
		monind=blks[blkind][symmattovec(k,l,blks[blkind][0])+1];//index of monomial (or -1 if already removed)
		//fprintf(stderr, "monind=%d\n",monind);
		if((k==rowind || l==rowind) && monind!=-1){//not already removed
		  flag=1;//something got removed in this iteration
		  if(debugfull){fprintf(stderr, "removing...(%d, %d) in block %d (monomial %d)\n", k, l, blkind, monind);
		    printvec(allmons[monind]+2, numv);
		    fprintf(stderr, "%d left\n", allmons[monind][0]-1);}
		  //total occurrence of monomial decreases
		  (allmons[monind][0])--;
		  //checks: 
		  if(allmons[monind][0]==0){
		    if(debug){fprintf(stderr, "removed all of: ");
		      printvec(allmons[monind]+2, numv);}
		    if(isinarray(lst2, LHtot, allmons[monind]+2, numv)>=0){
		      fprintf(stderr, "Something went wrong: cannot remove a monomial which figures on LHS.\n");exit(0);
		    }
		  }
		  else if(allmons[monind][0]==1 && isinarray(lst2, LHtot, allmons[monind]+2, numv)>=0){
		    fprintf(stderr, "Something went wrong: A monomial which figures on LHS can no longer be cancelled.\n");
		    printvec(allmons[monind]+2, numv);exit(0);
		  }

		  //Indicates that this entry in the block is now empty
		  blks[blkind][symmattovec(k,l,blks[blkind][0])+1]=-1;
		  //Removed a nondiagonal occurrence of the monomial
		  if(l!=k){
		    /* fprintf(stderr, "removing nondiagonal entry...(%d, %d) in block %d (monomial %d)\n", k, l, blkind, monind); */
		    /* printvec(allmons[monind]+2, numv); */
		    (nondiag[monind])--;
		  }
		}
	      }
	    }
	  }
	}
      }

    }

    q++;
  }

  //
  // Enter the main loop
  //

  //
  // A good indicator of whether a solution is artefact is whether it survives with a small positive tol below
  //

  p=0;
  numnewconstraints=0;
  tol=0.001;
  while(p<1){
    //
    // Update the blocks
    //
    for(j=0;j<totblk+numpossqmons;j++){//each block
      newblksz=0;//what's left in block
      if(blks[j][0]){//potential nonempty block
	//to hold indices of surviving square monomials
	if(debug){fprintf(stderr, "\nold block %d, size %d:\n", j, blks[j][0]);}
	for(k=0;k<blks[j][0];k++){
	  if(blks[j][symmattovec(k,k,blks[j][0])+1]!=-1){
	    if(debug){fprintf(stderr, "non-empty row-col: block %d, rowcol %d\n", j, k+1);
	      printvec(allmons[blks[j][symmattovec(k,k,blks[j][0])+1]]+2, numv);}
	    blkinds[j][newblksz++]=k;
	  }
	}
	if(newblksz==0){
	  if(debug){fprintf(stderr, "Newly empty block: %d\n", j);}
	}
	else if(j<totblk && newblksz==1){//now a redundant (repeated) block
	  if(debug){fprintf(stderr, "Block reduced to size 1 (remove) %d\n", j);}
	  (allmons[blks[j][symmattovec(blkinds[j][0],blkinds[j][0],blks[j][0])+1]][0])--;
	  newblksz=0;
	}

	if(debug){fprintf(stderr, "new block %d, size %d:\n", j, newblksz);}
	for(k=0;k<newblksz;k++){
	  for(l=k;l<newblksz;l++){
	    //Should be safe: only refers to later monomials
	    blks[j][symmattovec(k,l,newblksz)+1]=blks[j][symmattovec(blkinds[j][k],blkinds[j][l],blks[j][0])+1];
	    if(debug){fprintf(stderr, "%d ", blks[j][symmattovec(k,l,newblksz)+1]);
	    }}
	}
	blks[j][0]=newblksz;
      }
    }



    //
    // The surviving monomials
    // One constraint for each surviving monomial
    // Structure of a constraint: total memory allocation; index of monomial; 
    // then blocks of three (block number, row-index, column-index)
    //

    numconstraints=0;
    for(j=0;j<totmons;j++){
      if(allmons[j][0]){//surviving monomials
	newconstraint(&constraints, &numconstraints, j);
	allmons[j][1]=numconstraints-1;//constraint index

	if(debugfull){
	  if(isinarray(lst2, LHtot, allmons[j]+2, numv)>=0)
	    fprintf(stderr, "In LH poly: ");      
	  fprintf(stderr, "monomial to constraint: %d to %d\n", j, numconstraints-1);
	  printvec(allmons[j]+2, numv);}
      }
      else{
	if(debugfull){fprintf(stderr, "monomial deleted: ");
	  printvec(allmons[j]+2,numv);}
      }
    }

    //
    // Output block data and populate the constraints
    //
    numtrueblk=0;
    for(j=0;j<totblk+numpossqmons;j++){
      if(j==totblk)
	numsqblk=numtrueblk;//Number of blocks other than the single entry ones

      if(blks[j][0]){
	if(j<totblk){
	  if(debugfull){fprintf(stderr, "block %d, size %d, corresponding to basic monomial 1:\n", j, blks[j][0]);}
	}
	else{
	  if(debugfull){fprintf(stderr, "1X1 block %d corresponding to monomial:\n", j);
	    printvec(allmons[blks[j][1]]+2,numv);}
	}
	/* fprintf(stderr, "square monomials:\n"); */
	/* for(k=0;k<blks[j][0];k++) */
	/* 	printvec(sqmons[blkinds[j][k]],numv); */
	//printvec(blks[j]+1, (blks[j][0]*(blks[j][0]+1))/2);
	for(k=0;k<blks[j][0];k++){
	  for(l=k;l<blks[j][0];l++){
	    ind=allmons[blks[j][symmattovec(k,l,blks[j][0])+1]][1];//index of constraint
	    growconstraint(constraints, ind, numtrueblk, k, l);
	  }
	}
	numtrueblk++;
      }
    }
    if(numpossqmons==0)
      numsqblk=numtrueblk;//Number of blocks other than the single entry ones

    p++;
    if(debugfull){fprintf(stderr, "newconstraints=%d\n", numnewconstraints);}
    //
    //Start writing tempfiles/csdp.dat
    //

    if(!(fd=fopen("tempfiles/csdp.dat", "w"))){
      fprintf(stderr, "ERROR: \"tempfiles/csdp.dat\" could not be opened for reading.\n");
      exit(0);
    }
    //Number of constraints
    fprintf(fd, "%d\n", numconstraints+numnewconstraints+1);

    //Number of blocks
    if(ldegree){
      fprintf(fd, "%d\n", numtrueblk+totlmons+1);
      trueblksz=(int*)malloc((size_t) ((numtrueblk+totlmons+1)*sizeof(int)));
    }
    else{
      fprintf(fd, "%d\n", numtrueblk+1);
      trueblksz=(int*)malloc((size_t) ((numtrueblk+1)*sizeof(int)));
    }

    //Block sizes (stored in trueblksz)
    tot=0;
    for(j=0;j<totblk+numpossqmons;j++){
      if(blks[j][0]){
	fprintf(fd, "%d ", blks[j][0]);
	trueblksz[tot++]=blks[j][0];
      }
    }
    if(ldegree){//LHS parameter blocks
      if(posleft){
	for(j=0;j<totlmons;j++){
	  fprintf(fd, "1 ");
	  trueblksz[tot++]=1;
	}
      }
      else{
	for(j=0;j<totlmons;j++){
	  fprintf(fd, "2 ");
	  trueblksz[tot++]=2;
	}
      }
    }
    fprintf(fd, "1\n");//final block
    trueblksz[tot++]=1;


    //Constraint RHSs
    for(j=0;j<numconstraints;j++){//surviving monomials
      if(ldegree<=0 && (s=isinarray(lst, r, allmons[constraints[j][1]]+2, numv))>=0)
	fprintf(fd, "%.10f ", cflst[s]);
      else
	fprintf(fd, "0.0 ");
    }
    for(j=0;j<numnewconstraints;j++)
      fprintf(fd, "%.1f ", (double)(newconstraints[j][1])/1.0);//change if using nearhalfint

    //final constraint RHS
    //This seems to matter: if the poly has large coefficients, then setting this
    //small can make the algorithm fail (roundoff error?)
    //setting this too large can make a degree 0 strict problem fail
    //current heuristic: set it to 1% of largest coeff
    maxcf=0.0;
    for(s=0;s<r;s++){
      if(maxcf<fabs(cflst[s]))
	maxcf=fabs(cflst[s]);
    }
    if(strict && !ldegree)
      fprintf(fd, "%.4f\n",0.01*maxcf);//very heuristic
    else if(strict || ldegree)
      fprintf(fd, "1\n");
    else
      fprintf(fd, "0\n");

    //Objective function

    if(ldegree){
      if(!posleft){
	for(j=0;j<totlmons;j++)
	  fprintf(fd, "0 %d 1 1 -1.0\n0 %d 2 2 -1.0\n", numtrueblk+j+1, numtrueblk+j+1);
      }
      else{
	for(j=0;j<totlmons;j++)
	  fprintf(fd, "0 %d 1 1 -1.0\n", numtrueblk+j+1);
      }
    }
    else{// diagonal
      if(p>0){
	tot=0;
	/* for(j=0;j<totblk;j++){ */
	/*   if(blks[j][0]){ */
	/*     tot++; */
	/*     if(tot==2 || tot==4 || tot>6){ */
	/*       for(k=0;k<blks[j][0];k++) */
	/*       	for(l=k;l<blks[j][0];l++) */
	/*       	  if(l==k) */
	/*       	    fprintf(fd, "0 %d %d %d -1.0\n", tot, k+1, l+1); */
	/*       	  else */
	/*       	    fprintf(fd, "0 %d %d %d -0.5\n", tot, k+1, l+1); */
	/*     } */
	/*     else if(tot==1) */
	/*       for(k=9;k<blks[j][0];k++) */
	/* 	fprintf(fd, "0 %d %d %d 1.0\n", tot, k+1, k+1); */

	/*   } */
	/* } */
	/* for(j=totblk;j<totblk+numpossqmons;j++){ */
	/* 	if(blks[j][0]){ */
	/* 	  tot++; */
	/* 	  for(k=0;k<blks[j][0];k++) */
	/* 	    fprintf(fd, "0 %d %d %d 1.0\n", tot, k+1, k+1); */
	/* 	} */
	/* } */
	}
    }


    //constraints
    for(j=0;j<numconstraints;j++){
      if(constraints[j][1]>=0){
	if(debugfull){fprintf(stderr, "Constraint on monomial %d: \n", j);
	  printvec(allmons[constraints[j][1]]+2, numv);}
      }
      //    fprintf(stderr, "lines in constraint: %d\n", (constraints[j][0]-2)/3);
      for(k=0;k<(constraints[j][0]-2)/3;k++){
	fprintf(fd, "%d %d %d %d 1.0\n", j+1, constraints[j][3*k+2], constraints[j][3*k+3], constraints[j][3*k+4]);
      }
      if(ldegree>0 && (s=isinarray(lst2, LHtot, allmons[constraints[j][1]]+2, numv))>=0){//LHS parameters
	for(k=0;k<totlmons;k++){
	  if(lst1[s][k]){
	    if(posleft)
	      fprintf(fd, "%d %d 1 1 %.4f\n", j+1, numtrueblk+k+1, -lst1[s][k]);
	    else
	      fprintf(fd, "%d %d 1 2 %.4f\n", j+1, numtrueblk+k+1, -lst1[s][k]/2.0);
	    //fprintf(stderr, "printing %d --> %.1f\n", lst1[s][k], -((double)(lst1[s][k]))/2.0);  
	    //printvec(lst1[s],totlmons);
	  }
	}
      }
    }
    for(j=0;j<numnewconstraints;j++){
      for(k=0;k<(newconstraints[j][0]-2)/3;k++){
	fprintf(fd, "%d %d %d %d 1.0\n", numconstraints+j+1, newconstraints[j][3*k+2], newconstraints[j][3*k+3], newconstraints[j][3*k+4]);
      }
    }

    //final constraint


    if(ldegree && !strict){
      numallblk=numtrueblk+totlmons+1;
      if(posleft){
	for(k=numtrueblk;k<numtrueblk+totlmons;k++)
	  fprintf(fd, "%d %d 1 1 1.0\n", numconstraints+j+1, k+1);
      }
      else{
	for(k=numtrueblk;k<numtrueblk+totlmons;k++)
	  fprintf(fd, "%d %d 1 2 1.0\n", numconstraints+j+1, k+1);
      }
      fprintf(fd, "%d %d 1 1 -1.0\n", numconstraints+j+1, k+1);
    }
    else{
      numallblk=numsqblk;
      for(k=0;k<numpossqmons;k++){
	if(blks[k+totblk][0]){//nonempty
	  numallblk++;
	  fprintf(fd, "%d %d 1 1 1.0\n", numconstraints+j+1, numallblk);
	}
      }
      if(ldegree)//skip the blocks associated with LHS parameters
	numallblk+=totlmons+1;
      else
	numallblk++;
      fprintf(fd, "%d %d 1 1 -1.0\n", numconstraints+j+1, numallblk);
    }

    fclose(fd);

    if(debug){fprintf(stderr, "numallblk=%d, totlmons=%d\n", numallblk, totlmons);}

    //
    //Finished writing tempfiles/csdp.dat, now run
    //

    if(!(fp=popen("csdp tempfiles/csdp.dat tempfiles/tmp.sol", "r"))){
      perror("couldn't run command \"csdp tempfiles/csdp.dat tempfiles/tmp.sol\". EXITING\n");exit(0);
    }



    retstat=-2;
    while(fgets(path, sizeof(path), fp)){
      if(debug){fprintf(stderr, "%s", path);}
      if(strstr(path, "Declaring primal infeasibility"))
	retstat=-1;
      else if(strstr(path, "Partial Success: SDP solved with reduced accuracy"))
	retstat=2;
      else if(strstr(path, "Success: SDP solved"))
	retstat=0;
    }
    pclose(fp);
    if(retstat==-1){
      free((char*)trueblksz);
      for(j=0;j<numconstraints;j++)
	free((char*)constraints[j]);
      free((char*)constraints);
      break;
    }

    outstr=readfileintostr("tempfiles/tmp.sol");
    datstr=strstr(outstr,"\n2");
    datstr++;

    //  fprintf(stderr, "datstr=%s\n", datstr);
    outdat=genfloatvec(datstr, &datlen);
    if(debug){fprintf(stderr, "datlen=%d\n", datlen);}
    if(datlen%5!=0){
      fprintf(stderr, "ERROR reading output: data not a multiple of 5\n");
      exit(0);
    }
    // To hold the data
    blksout=(double **)malloc((size_t) (numallblk*sizeof(double*)));
    for(j=0;j<numallblk;j++){
      blksz=trueblksz[j];
      if(debugfull){fprintf(stderr, "j=%d, blksz=%d\n", j, blksz);}
      blksout[j]=(double *)malloc((size_t) (((blksz*(blksz+1))/2+1)*sizeof(double)));
      blksout[j][0]=blksz;
      for(k=1;k<(blksz*(blksz+1))/2+1;k++)//default is zero
	blksout[j][k]=0.0;
    }

    j=0;
    maxval=0.0;
    while(j<datlen){//reprint output data
      if(debug){fprintf(stderr, "%d %d %d %d %f\n", (int)outdat[j], (int)outdat[j+1], (int)outdat[j+2], (int)outdat[j+3], outdat[j+4]);}
      maxval=fabs(outdat[j+4])>maxval?fabs(outdat[j+4]):maxval;
      j+=5;
    }
    if(debug){fprintf(stderr, "$$$maxval=%.4f\n$$$\n", maxval);}

    j=0;
    while(j<datlen){
      blksout[(int)(outdat[j+1])-1][symmattovec((int)(outdat[j+2])-1,(int)(outdat[j+3])-1,blksout[(int)(outdat[j+1])-1][0])+1]=outdat[j+4];
      j+=5;
    }

    numallblk=0;outpolylen=0;
    for(j=0;j<totblk+numpossqmons;j++){
      if(blks[j][0]){
	tot=1;
	outblkpolylen=0;
	for(k=0;k<trueblksz[numallblk];k++){
	  for(l=k;l<trueblksz[numallblk];l++){
	    if(l==k){
 	      outpolylen=addpolyterm(outpolylen, allmons[blks[j][tot]]+2, blksout[numallblk][tot], numv, &outpoly, &outcfs, &ind);
	      outblkpolylen=addpolyterm(outblkpolylen, allmons[blks[j][tot]]+2, blksout[numallblk][tot], numv, &outblkpoly, &outblkcfs, &ind);
	    }
	    else{
	      outpolylen=addpolyterm(outpolylen, allmons[blks[j][tot]]+2, 2.0*blksout[numallblk][tot], numv, &outpoly, &outcfs, &ind);
	      outblkpolylen=addpolyterm(outblkpolylen, allmons[blks[j][tot]]+2, 2.0*blksout[numallblk][tot], numv, &outblkpoly, &outblkcfs, &ind);
	    }
	    //fprintf(stderr, "%.4f: ", blksout[numallblk][tot]); printvec(allmons[blks[j][tot]]+2,numv);
	    tot++;
	  }
	}//finished dealing with this block

	if(debug){
	  if(j<totblk || fabs(outblkcfs[0])>0.001){
	    fprintf(stderr, "########true block %d\n", numallblk);
	    cerr << lrexptopolysimp(outblkpoly, outblkcfs, 0.001, outblkpolylen, numv) << endl;
	  }
	  for(k=0;k<outblkpolylen;k++){
	    if(fabs(outblkcfs[k])>0.001){
	      fprintf(stderr, "%.3f: ", outblkcfs[k]);
	      printvec(outblkpoly[k],numv);
	    }
	  }
	  if(j<totblk || fabs(outblkcfs[0])>0.001)
	    fprintf(stderr, "########\n");
	}
	ifree(outblkpoly,outblkpolylen);
	free((char*)outblkcfs);

	numallblk++;
      }
    }

    //RHS poly
    if(debug){
      fprintf(stderr, "\n****RHS poly****\n");
      cerr << lrexptopolysimp(outpoly, outcfs, 0.001, outpolylen, numv) << endl;
      for(j=0;j<outpolylen;j++){
	fprintf(stderr, "%.2f: ", outcfs[j]);
	printvec(outpoly[j],numv);
      }
      fprintf(stderr, "\n*********\n");
    }

    if(ldegree){
      outblkpolylen=0;
      for(s=0;s<LHtot;s++){
	if(posleft){
	  for(k=0;k<totlmons;k++)
	    outpolylen=addpolyterm(outpolylen, lst2[s], -1.0*blksout[numallblk+k][1]*lst1[s][k], numv, &outpoly, &outcfs, &ind);
	}
	else{
	  for(k=0;k<totlmons;k++)
	    outpolylen=addpolyterm(outpolylen, lst2[s], -1.0*blksout[numallblk+k][2]*lst1[s][k], numv, &outpoly, &outcfs, &ind);
	}
      }
      if(posleft){
	for(k=0;k<totlmons;k++)
	  outblkpolylen=addpolyterm(outblkpolylen, lmons[k], blksout[numallblk+k][1], numv, &outblkpoly, &outblkcfs, &ind);
      }
      else{
	for(k=0;k<totlmons;k++)
	  outblkpolylen=addpolyterm(outblkpolylen, lmons[k], blksout[numallblk+k][2], numv, &outblkpoly, &outblkcfs, &ind);
      }

      if(debug){
	fprintf(stderr, "########LHS poly\n");
	cerr << lrexptopolysimp(outblkpoly, outblkcfs, 0.001, outblkpolylen, numv) << endl;
	for(k=0;k<outblkpolylen;k++){
	  if(fabs(outblkcfs[k])>0.001){
	    fprintf(stderr, "%.3f: ", outblkcfs[k]);
	    printvec(outblkpoly[k],numv);
	  }
	}
	fprintf(stderr, "########\n");
      }
      ifree(outblkpoly,outblkpolylen);
      free((char*)outblkcfs);
    }
    else{
      for(s=0;s<r;s++)
	outpolylen=addpolyterm(outpolylen, lst[s], (double)(-cflst[s]), numv, &outpoly, &outcfs, &ind);
    }

    if(debug){
      fprintf(stderr, "Final polynomial (to 4 d.p.):\n");
      cerr << lrexptopolysimp(outpoly, outcfs, 0.0001, outpolylen, numv) << endl;
      for(j=0;j<outpolylen;j++){
	if(fabs(outcfs[j])>0.0001){
	  fprintf(stderr, "%.4f: ", outcfs[j]);
	  printvec(outpoly[j],numv);
	}
      }
      fprintf(stderr, "\n***********************************\n\n");
    }

    if(numnewconstraints>0){
      for(j=0;j<numnewconstraints;j++)
	free((char*)newconstraints[j]);
      free((char*)newconstraints);
    }
    numnewconstraints=0;


    // remove rowcols which seem to be zero
    numallblk=0;flag=0;
    for(j=0;j<totblk+numpossqmons;j++){
      if(blks[j][0]){
	for(k=0;k<trueblksz[numallblk];k++){
	  if(nearint(blksout[numallblk][symmattovec(k,k,blks[j][0])+1], 0, 0, tol)!=-1000){
	    for(l=0;l<trueblksz[numallblk];l++){//rowcol to zero: start at l=0
	      if(blks[j][symmattovec(k,l,blks[j][0])+1]!=-1){//not already set to zero by some other deletion
		flag=1;
		if(debugfull){fprintf(stderr, "setting block %d, entry (%d,%d) = %d to zero.\n", j, k,l, symmattovec(k,l,blks[j][0])+1);}
		(allmons[blks[j][symmattovec(k,l,blks[j][0])+1]][0])--;
		blks[j][symmattovec(k,l,blks[j][0])+1]=-1;
	      }
	    }//end of rowcol to zero
	  }
	}
	numallblk++;
      }
    }

    //Any values to be manually fixed? E.g. 
    /* addnewconstraint(&numnewconstraints, &newconstraints, 1, 16, 16, 9); */

    //End of rowcol removal

    ifree(outpoly,outpolylen); 
    free((char*)outcfs);
    free((char*)trueblksz);
    for(j=0;j<numallblk;j++)
      free((char*)blksout[j]);
    free((char*)blksout);
    free((char*)outdat);
    free(outstr);

    for(j=0;j<numconstraints;j++)
      free((char*)constraints[j]);
    free((char*)constraints);

    if(!flag)//nothing changed on this iteration
      break;

  }

  if(numnewconstraints>0){
    for(j=0;j<numnewconstraints;j++)
      free((char*)newconstraints[j]);
    free((char*)newconstraints);
  }

  ifree(allmons,totmons);
  ifree(blks,totblk+numpossqmons);
  for(j=0;j<totblk+numpossqmons;j++)
    if(blkinds[j])
      free ((char *)(blkinds[j]));
  free((char *)blkinds);

  free_imatrix(sqmons,0,choose(numv+sqmaxdeg,numv)-1,0,numv-1);
  ifree(lst1,LHtot);ifree(lst2,LHtot);
  ifree(monrefs,numposmons);
  free_imatrix(lmons, 0, totlmons-1,0,numv-1);
  free_imatrix(posmons,0,numposmons-1,0,numv-1);
  if(possqmons)
    free_imatrix(possqmons,0,numpossqmons-1,0,numv-1);
  free((char*)nondiag);

  return retstat;
}

//get maximum value, roughly to get a sense of the size of coeffs
double getmaxrough(double *cflst, long r){
  long i;
  double mx=0.0;
  double tol=1.0;//don't need to be very accurate
  for(i=0;i<r;i++){
    if((fabs(cflst[i])-mx)>tol){mx=fabs(cflst[i]);}
  }
  return mx;
}

//Make sure no coefficient much exceeds "maxallowed"
double scalecoeffs(double *cflst, long r, double maxallowed){
  double scl;
  long i;
  double mx=getmaxrough(cflst,r);//roughly max value
  if(mx<maxallowed)
    return 1;
  scl=maxallowed/mx;

  for(i=0;i<r;i++)
    cflst[i]*=scl;

  return scl;
}

//
// Wrapper function which extracts the variables, etc.
//

int ppp(ex tmp1, int preprocess, int maxldegree, int strict, int *ldeg, int debug){
  char **pvars;
  int numv=0,numv1=0;//=polyvars1(tmp, &pvars);
  int **lst=NULL;
  double *cflst;
  int **tmplst;
  long r,s;
  int tmppolymindeg,tmppolymaxdeg;
  int polymindeg,polymaxdeg;
  int ret=-1;
  int ldegree=0;
  int ldegreemin=0;
  ex tmp=polysimp(tmp1, &pvars, &numv1);
  int mkhom=1;//work only with homogeneous polynomials
  double maxallowed=1000.0;//largest allowed coefficient (for polynomial scaling)
  double scl;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering ppp.\n");}
  //cerr << "entering ppp\n";
  //cerr << "numv = " << numv << endl << tmp << endl << tmp1<< endl;

  if(debugfull){
    fprintf(stderr, "new variable order:\n");
    for(r=0;r<numv1;r++)
      fprintf(stderr, "%s\n", pvars[r]);
  }

  if(debug)
    cerr << "checking polynomial:\n" << tmp << endl;

  //maxallowed=-1;//to stop scaling
  extolst1(tmp, numv1, &tmplst, &cflst, &r, &tmppolymindeg, &tmppolymaxdeg);
  if(strict && maxallowed>0){//only scale with "strict"?
    scl=scalecoeffs(cflst,r,maxallowed);
    if(debugfull && scl<0.9999){
      fprintf(stderr, "Coefficients scaled. Scaling factor: %.4f\n", scl);
    }
  }
  //lst=imatrix(0,r-1,0,numv-1);
  numv=polylistsimp(tmplst, numv1, r, &lst, mkhom, tmppolymindeg, tmppolymaxdeg, &polymindeg, &polymaxdeg);
  if(debug && numv!=numv1){
    fprintf(stderr, "After manipulation, number of variables=%d. Poly:\n",numv);
    for(s=0;s<r;s++)
      printvec(lst[s], numv);
  }

  while(ret!=0 && ldegree<=maxldegree){
    if(polymaxdeg==polymindeg)//homogeneous poly: choose LHS multiplier to be homogeneous too: any theoretical justification?
      ldegreemin=ldegree;
    ret=ppp1(lst, cflst, r, numv, polymindeg, polymaxdeg, preprocess, ldegreemin, ldegree, strict, debugfull);
    ldegree++;
  }

  freearraydat(pvars, numv1);
  ifree(lst, r);
  ifree(tmplst,r); 
  free((char*)cflst);
  if(!ret || ret==2)
    *ldeg=ldegree-1;
  else
    *ldeg=-1;

  if(debug){//return codes from ppp1
    fprintf(stderr, "Exiting ppp. Return code = %d. This means that ", ret);
    if(ret==0) 
      fprintf(stderr, "the problem was feasible.\n");
    else if(ret==-1) 
      fprintf(stderr, "the problem was infeasible.\n");
    else if(ret==2) 
      fprintf(stderr, "the problem was probably feasible, but was solved with reduced accuracy.\n");
    else if(ret==-2) 
      fprintf(stderr, "something went wrong.\n");
  }

  return ret;
}


//
// Wrapper function which extracts the variables, etc.
// Assume homogeneous
//

int pppfull(ex tmp, int numv, int preprocess, int maxldegree, int strict, int *ldeg, int debug){
  int **lst=NULL;
  double *cflst;
  long r;
  int polymindeg,polymaxdeg;
  int ret=-1;
  int ldegree=0;
  double maxallowed=1000.0;//largest allowed coefficient (for polynomial scaling: set to -ve for no scaling)
  double scl;
  int debugfull=(debug<=0)?0:debug-1;

  if(debug){fprintf(stderr, "\n###Entering pppfull. ");}
  //cerr << "numv = " << numv << endl << tmp << endl << tmp1<< endl;

  if(debug)
    cerr << "Checking polynomial:\n" << tmp << endl;

  extolst1(tmp, numv, &lst, &cflst, &r, &polymindeg, &polymaxdeg);
  if(strict && maxallowed>0){//only scale with "strict"?
    scl=scalecoeffs(cflst,r,maxallowed);
    if(debugfull && scl<0.9999){
      fprintf(stderr, "Coefficients scaled. Scaling factor: %.4f\n", scl);
    }
  }

  if(polymindeg!=polymaxdeg){
    fprintf(stderr, "ERROR in pppfull: polynomial is not homogeneous.\n");
    exit(0);
  }

  while(ret!=0 && ldegree<=maxldegree){
    ret=pppfull1(lst, cflst, r, numv, polymindeg, polymaxdeg, preprocess, ldegree, ldegree, strict, debugfull);
    ldegree+=2;
  }

  ifree(lst, r);
  free((char*)cflst);
  if(!ret || ret==2)
    *ldeg=ldegree-2;
  else
    *ldeg=-1;

  if(debug){fprintf(stderr, "Exiting pppfull: return code %d.\n", ret);}
  return ret;
}

// Generate monomials upto degree maxdeg in numv variables
// revlex?
// mons must have dimensions numv+maxdeg choose maxdeg X numv
//

int genmons(int numv, int maxdeg, int **mons){
  int i,j,k,t=numv+maxdeg;
  int xc[maxdeg],yc[numv];
  int flag=1;
  int prev;
  int tot=0;

  firstcomb(xc,t,maxdeg);

  while(flag){
    //printvec(xc,maxdeg);
    //complement
    k=0;j=0;i=0;
    while(i<maxdeg){
      while(xc[i]>k){
	yc[j]=k;
	j++;k++;
      }
      i++;k++;
    }
    while(j<numv){
      yc[j]=k;
      j++;k++;
    }
    //printvec(xc,numv);
    //printvec(yc,numv);
    prev=t;
    for(i=0;i<numv;i++){
      mons[tot][numv-i-1]=prev-yc[numv-i-1]-1;
      prev=yc[numv-i-1];
    }
    //printvec(mons[tot],numv);
    tot++;
    flag=nextcomb(xc, t, maxdeg);
  }
  return tot;
}


int **genmons1(int numv, int mindeg, int maxdeg, int *tot){
  int i,j;
  int maxmons=choose(numv+maxdeg,numv);
  int minmons=choose(numv+mindeg-1,numv);
  int **mons1=imatrix(0, maxmons-1,0,numv-1);
  int **mons;
  (*tot)=maxmons-minmons;
  mons=imatrix(0, (*tot)-1,0,numv-1);
  genmons(numv,maxdeg,mons1);
  for(i=minmons;i<maxmons;i++){
    for(j=0;j<numv;j++)
      mons[i-minmons][j]=mons1[i][j];
    //printvec(mons1[i], numv);
  }
  free_imatrix(mons1, 0, maxmons-1,0,numv-1);
  return mons;

}



// The first polynomial has coefficients which are themselves unknowns; these are added as extra variables
int polyprod(int **lst1, int len1, double *cfs2, int **lst2, int len2, int numv, double ***lst1out, int ***lst2out){
  int i,j,k;
  double **tmpout;
  int **tmpout1;
  int tot=0,tot1;
  int pos;

  if(len1<=0){
    fprintf(stderr, "Error in polyprod. First poly must not be empty. Exiting.\n");exit(0);
  }
  //
  //To hold the product as a polynomial (with repetition) with the unknown coefficients as the first len1 unknown variables
  //
  tmpout=dmatrix(0, len1*len2-1, 0,len1-1);
  tmpout1=imatrix(0, len1*len2-1, 0,numv-1);
  for(i=0;i<len1;i++){
    for(j=0;j<len2;j++){
      for(k=0;k<len1;k++){
	if(k==i)
	  tmpout[len2*i+j][k]=cfs2[j];
	else
	  tmpout[len2*i+j][k]=0.0;
      }     
      for(k=0;k<numv;k++){
	tmpout1[len2*i+j][k]=lst1[i][k]+lst2[j][k];
      }
    }
  }

  //Now get rid of repeated coefficients, and split into two polys in one step

  for(i=0;i<len1*len2;i++){
    tot1=addvec1(tot, tmpout1[i], numv, lst2out,1,&pos);
    if(tot1>tot)//new monomial
      addvec1(tot,tmpout[i],len1,lst1out,0,&pos);
    else//existing monomial; update coefficient (itself a monomial)
      vecadd((*lst1out)[pos],tmpout[i],len1);
    tot=tot1;
  }

  free_dmatrix(tmpout, 0, len1*len2-1, 0,len1-1);
  free_imatrix(tmpout1, 0, len1*len2-1, 0,+numv-1);
  return tot;
}


// In this version, the first entry is the number of times the int string occurs
// The second entry is used to reference the constraint associated with the monomial
// The remaining entries are the monomial itself

int addoneintstra(int k, int *s, int len1, int ***t, int *ind)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. */
{
  int i=0,j=0, flag=0;
  while(j<k){
    flag=1;
    for(i=2;i<2+len1;i++){
      if((*t)[j][i]!=s[i-2]){
	flag=0;
	break;
      }
    }
    if(flag){ // found (update where)
      (*ind)=j;
      ((*t)[j][0])++;
      return k;
    }
    j++;
  }
  //not found
  //    fprintf(stderr, "not found: %d\n", k);
  //    printvec(s, len1);printvec(s1,len2);
  if(!k)
    (*t) = (int**) malloc(sizeof(int*)*1);
  else
    (*t)=(int**) realloc((*t), sizeof(int*)*(k+1));
  (*t)[k] = (int*) malloc(sizeof(int)*(len1+2));
  (*t)[k][0]=1;
  (*t)[k][1]=0;
  for(i=2;i<2+len1;i++)
    (*t)[k][i] = s[i-2];
  (*ind)=k;
  return k+1;
}

// ceiling of max(x,0)/y
// y >0

int ceiling(int x, int y){
  if(x<=0)
    return 0;
  return 1+((x-1)/y);
}

//encoding for a symmetric matrix as a vector

int symmattovec(int i, int j, int N){
  if (i<=j)
    return i*N-(i-1)*i/2+j-i;
  else
    return j*N-(j-1)*j/2+i-j;
}


void newconstraint(int ***constraints, int *numconstraints, int mon_ind){
  (*numconstraints)++;
  if((*numconstraints)==1)
    (*constraints)=(int **)malloc((size_t) ((*numconstraints)*sizeof(int*)));
  else
    (*constraints)=(int **)realloc((*constraints), (size_t)((*numconstraints)*sizeof(int*)));

  (*constraints)[(*numconstraints)-1]=(int *)malloc((size_t)(2*sizeof(int)));
  (*constraints)[(*numconstraints)-1][0]=2;//first entry in constraint is storage size
  (*constraints)[(*numconstraints)-1][1]=mon_ind;//second entry in constraint is index of monomial
}


void growconstraint(int **constraints, int ind, int blk, int row, int col){
  constraints[ind]=(int *)realloc(constraints[ind], (size_t)((constraints[ind][0]+3)*sizeof(int)));
  constraints[ind][(constraints[ind][0])++]=blk+1;
  constraints[ind][(constraints[ind][0])++]=row+1;
  constraints[ind][(constraints[ind][0])++]=col+1;
}

double *genfloatvec(char *line, int *num){
  int i=0,j,k;
  char *wd;
  double *dvec=NULL;
  int numt=0;
  int wdmax=0;
  while(line[i]){
    while(line[i] && isspace(line[i]))
      i++;
    j=i;
    while(line[i] && !(isspace(line[i])))
      i++;
    if(i>j){
      numt++;
      wdmax=(i-j)>wdmax?i-j:wdmax;
    }
  }

  dvec=(double *)malloc((size_t) ((numt)*sizeof(double)));
  wd=(char *)malloc((size_t) ((wdmax+1)*sizeof(wd)));

  numt=0;i=0;
  while(line[i]){
    while(line[i] && isspace(line[i]))
      i++;
    j=i;k=0;
    while(line[i] && !(isspace(line[i])))
      wd[k++]=line[i++];
    wd[k++]=0;
    if(i>j)
      dvec[numt++]=atof(wd);
  }

  free(wd);
  *num=numt;
  return dvec;
}

//Add the new term (cf)(newvec) into the poly (cflst)(veclst)
int addpolyterm(int lstlen, int *newvec, double cf, int veclen, int ***veclst, double **cflst, int *pos){
  int i=0;
  (*pos)=-1;

  if(lstlen==0){
    (*veclst) = (int**) malloc(sizeof(int*) * 1);
    (*cflst) = (double*) malloc(sizeof(double) * 1);
    (*pos)=0;
  }
  //not already present
  else if(((*pos)=isinarray((*veclst), lstlen, newvec, veclen))<0){
    (*veclst)=(int**) realloc((*veclst), sizeof(int*) *(lstlen+1));
    (*cflst) = (double*) realloc((*cflst), sizeof(double)*(lstlen+1));
    (*pos)=lstlen;
  }
  else{
    (*cflst)[*pos]+=cf;
    return lstlen;
  }

  (*veclst)[lstlen] = (int*) malloc(sizeof(int) * (veclen));
  for(i=0;i<veclen;i++)
    (*veclst)[lstlen][i] = newvec[i];
  (*cflst)[*pos]=cf;
  return lstlen+1;

}

int nearint(double p, int imin, int imax, double tol){
  int j;
  for(j=imin;j<=imax;j++){
    if(fabs(p-(double)j)<tol){
      //fprintf(stderr, "%f is near %d\n", p, j);
      return j;
    }
  }
  return -1000;
}

//generate all even powers of variables between degree mindeg and maxdeg

int **genevenpows(int numv, int mindeg, int maxdeg, int *tot){
  int i,j,k;
  int **mons;
  (*tot)=0;
  if(maxdeg%2!=0 || mindeg%2!=0){
    fprintf(stderr, "Error in genevenpows: maximum and minimum degrees need to be even.\n");
    exit(0);
  }
  if(maxdeg<mindeg){
    fprintf(stderr, "Error in genevenpows: maximum degree needs to be >= minimum degree.\n");
    exit(0);
  }

  mons=imatrix(0, numv*(1+(maxdeg-mindeg)/2)-1,0,numv-1);
  for(i=mindeg;i<=maxdeg;i+=2){
    for(j=0;j<numv;j++){//even powers of each variable
      for(k=0;k<numv;k++){
	if(k==j)
	  mons[(*tot)][k]=i;
	else
	  mons[(*tot)][k]=0;
      }
      (*tot)++;
    } 
  }
  return mons;

}

int **genemptymonslist(int numv, int mindeg, int maxdeg, int *tot){
  (*tot)=0;
  return NULL;
}


