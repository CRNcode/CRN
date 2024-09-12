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
#include "analysereacs.h"
#include "isomorph.h"
#include "convex.h"

//Generate all at most bimolecular half-reaction vectors on n species

int genallbicomplexes(int n, int N, int **l){
  int i,j,k=1;

  inittozero(l,N,n);//initialise to zero

  //unicomplexes
  for(i=0;i<n;i++){
    l[k][i]=1;k++;
  }

  //bicomplexes
  for(i=0;i<n;i++){
    l[k][i]=2;k++;
  }
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      l[k][i]=1;l[k][j]=1;k++;
    }
  }

  return k;
}

//Generate all at most trimolecular half-reaction vectors on n species
//important that we keep the order consistent with genallbicomplexes

int genalltricomplexes(int n, int N, int **l){
  int i,j,j1,k=1;//First is empty complex

  inittozero(l,N,n);//initialise to zero

  //unicomplexes
  for(i=0;i<n;i++){
    l[k][i]=1;k++;
  }
  //bicomplexes
  for(i=0;i<n;i++){
    l[k][i]=2;k++;
  }
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      l[k][i]=1;l[k][j]=1;k++;
    }
  }

  //tricomplexes
  for(i=0;i<n;i++){
    for(j=i;j<n;j++){
      for(j1=j;j1<n;j1++){
	l[k][i]++;l[k][j]++;l[k][j1]++;k++;
      }
    }
  }

  return k;
}


//Generate all at most trimolecular half-reaction (0,1) vectors on n species
int genallzeroonetricomplexes(int n, int N, int **l){
  int i,j,j1,k=1;//First is empty complex

  inittozero(l,N,n);//initialise to zero

  //unicomplexes
  for(i=0;i<n;i++){
    l[k][i]=1;k++;
  }
  //bicomplexes
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      l[k][i]=1;l[k][j]=1;k++;
    }
  }

  //tricomplexes
  for(i=0;i<n-2;i++){
    for(j=i+1;j<n-1;j++){
      for(j1=j+1;j1<n;j1++){
	l[k][i]=1;l[k][j]=1;l[k][j1]=1;k++;
      }
    }
  }

  return k;
}


//Generate all at most tetramolecular half-reaction vectors on n species
//important that we keep the order consistent with genallbicomplexes and genalltricomplexes

int genalltetracomplexes(int n, int N, int **l){
  int i,j,j1,j2,k=1;//First is empty complex

  inittozero(l,N,n);//initialise to zero

  //unicomplexes
  for(i=0;i<n;i++){
    l[k][i]=1;k++;
  }
  //bicomplexes
  for(i=0;i<n;i++){
    l[k][i]=2;k++;
  }
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      l[k][i]=1;l[k][j]=1;k++;
    }
  }

  //tricomplexes
  for(i=0;i<n;i++){
    for(j=i;j<n;j++){
      for(j1=j;j1<n;j1++){
	l[k][i]++;l[k][j]++;l[k][j1]++;k++;
      }
    }
  }

  //tetracomplexes
  for(i=0;i<n;i++){
    for(j=i;j<n;j++){
      for(j1=j;j1<n;j1++){
	for(j2=j1;j2<n;j2++){
	  l[k][i]++;l[k][j]++;l[k][j1]++;l[k][j2]++;k++;
	}
      }
    }
  }

  return k;
}


//Generate all at most k-molecular complexes on n species
//important that we keep the order consistent with genallbicomplexes
//All complexes up to order molecularity allowed on the RHS
int genallmulticomplexes(int n, int N, int **l, int molecularity){
  int i,j,k=1;//First is empty complex
  int vec[n];
  int vecout[n];
  int mol;
  int flag;
  inittozero(l,N,n);//initialise to zero

  //unicomplexes
  for(i=0;i<n;i++){
    l[k][i]=1;k++;
  }

  //bicomplexes
  for(i=0;i<n;i++){
    l[k][i]=2;k++;
  }
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      l[k][i]=1;l[k][j]=1;k++;
    }
  }

  //higher order complexes (order >=2)
  for(mol=3;mol<=molecularity;mol++){

    flag=1;
    firstcomb(vec,mol+n-1,n-1);
    while(flag==1){
      starbar(vec,mol+n-1,n-1, vecout);
      for(i=0;i<n;i++){
	l[k][i]=vecout[i];
      }
      k++;
      flag=nextcomb(vec,mol+n-1,n-1);
    }
  }
  return k;
}


//are the total degrees of the species in ascending order
bool goodorder(int **L, int **R, int n, int m){
  int i, j, p=0, q=0;
  if(n<=0)
    return 1;

  for(j=0;j<m;j++)//zeroth species
    p+=L[0][j]+R[0][j];

  for(i=1;i<n;i++){
    q=0;
    for(j=0;j<m;j++)
      q+=L[i][j]+R[i][j];
    if(q<p)
      return 0;
    else
      p=q;
  }
  return 1;
}



// a binary number of length 6 to a base 10 integer
int bin6toint(bool *vec){
  return 32*vec[0]+16*vec[1]+8*vec[2]+4*vec[3]+2*vec[4]+vec[5];
}

// To avoid using the NAUTY function amtog and cut down on file IO
// Only for graphs with fewer than 63 vertices (n+m<31), r=2*(n+m)
// V already padded to length 6*l > n*n
// out must have length l+3;
// No error checking

void amtodig(bool *V, int l, char out[]){
  int j;
  //  out[0]=38;
  //  out[1]=63+r;
  for(j=0;j<l;j++){
    out[2+j]=(char)(bin6toint(V)+63);
    V+=6;
  }
  //  out[2+l]=0;
  return;
}

// The number of shortg processes currently running: very crude!
int numshortg(){
  FILE *fpp=popen("pgrep -c shortg", "r");
  char path[10];
  fgets(path, sizeof(path), fpp);
  pclose(fpp);
  return atoi(path);
}

//concatenate all files beginning with fname and followed by a number
//into a single file \"fname\", and then get the short command ready
//and return file pointer to this command (command to be queued)
FILE *mergeandshort(const char fname[], const char ptn[], int numparts){
  int i;
  char systemcom[1024];
  FILE *fd;
  FILE *fp;
  if(!(fd=fopen(fname, "w"))){//make sure it's empty
    fprintf(stderr, "ERROR in mergeandshort: could not open file \"%s\" for writing.\n", fname);exit(0);
  }
  fclose(fd);

  for(i=0;i<numparts;i++){
    sprintf(systemcom, "cat %s_%d >> %s && rm %s_%d", fname, i, fname, fname, i);
    system(systemcom);
  }
  //sprintf(systemcom, "../../nauty/nauty26r5/shortg -q -f%s %s", ptn, fname);
  sprintf(systemcom, "nauty-shortg -q -f%s %s", ptn, fname);
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  return fp;
}

// Is any species trivial? Returns 1 if all species participate in reacs
int has_triv_spec(int **imatL, int **imatR, int n, int m){
  int flag,i=0,j;  

  while(i<n){//each species
    flag=0;//potentially trivial
    for(j=0;j<m;j++){//each reaction
      if(imatL[i][j] || imatR[i][j]){//nontrivial species
	flag=1;
	break;
      }
    }
    if(!flag)//trivial species
      return 1;
    i++;
  }
  return 0;
}


// Generate at most bimolecular CRNs involving n species and m reactions
// the openswtch decides whether we allow 0 --> X, X --> 0, 
// type reactions or not. 
// open means create only fully open CRNs
// rev means reversible CRNs only
// creates lots of temporary files: there must be a directory \"tempfiles\"
unsigned long genbimol(int n, int m, char *outfile, int open, int rev, const char procstr[], int debug){
  int i,j,k,tot=0,tot1=0,flag=1;
  //long milcount=0;
  long alltot=0;
  int m0=rev?(2*m):m;//number of reactions
  int r1=n+m0;
  int r2=4*r1*r1;//minimum dimension of V
  int N=(n+1)*(n+2)/2;//=1+n+n*(n+1)/2
  int Rvec[N*N][2];//all reactions (LHS index, RHS index)
  int **H=imatrix(0, N-1, 0, n-1);//complexes
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  int clen=Vlen/6;
  char out[clen+3];
  int fnum=rev?(2*m+1):(m+1)*(m+1);//equiv classes
  FILE *fd, *fp[fnum][3*m+1];
  FILE *fptemp[4]; // for temporary shortg processes
  FILE *fd1[fnum][3*m+1];
  int oldest_short=0;
  int fptempflag=0;
  char *ptn, *str;
  char **array;
  char systemcom[1024];
  int **L=imatrix(0,n-1,0,m0-1);
  int **R=imatrix(0,n-1,0,m0-1);
  //Transpose of irreversible stoich. mat. (only for rev.)
  int **Gt=imatrix(0,m0-1,0,n-1); 
  int xc[m];
  //Total degree between m and 4*m
  int num2l, num2r, inv, totdeg;
  char current_file[30];
  int fsz[fnum][3*m+1];
  int part[fnum][3*m+1];
  int maxproc=4;//maximum shortg processes to run simultaneously
  FILE *fp0;
  char path[1024];
  unsigned long totnetworks;

  //clear file
  if(!(fd=fopen("tempfiles/_tmp2", "w"))){
    fprintf(stderr, "ERROR in genbimol: could not open file \"tempfiles/_tmp2\" for writing.\n");exit(0);
  }
  fclose(fd);

  if(r1>30){
    fprintf(stderr, "ERROR in genbimol: not sure we can handle the size. Exiting.\n");
    exit(0);
  }

  //initialise
  for(i=0;i<fnum;i++){
    for(j=0;j<=3*m;j++){
      fsz[i][j]=0;part[i][j]=0;
    }
  }
  k=genallbicomplexes(n, N, H);
  if(debug){fprintf(stderr, "total half reacs = %d\n", k);}

  //generate the reactions as index-pairs
  for(i=0;i<k;i++){
    for(j=0;j<k;j++){
      if((!rev && j!=i && (!open || (i!=0 && j!=0) || (i==0 && j>n) || (j==0 && i>n))) || (rev && j>i && (!open || i!=0 || (i==0 && j>n)))){
	Rvec[tot][0]=i;Rvec[tot][1]=j;
	tot++;
      }
    }
  }
  if(debug){fprintf(stderr, "total reacs = %d\n", tot);}
  if(m>tot){
    fprintf(stderr, "Cannot choose %d reactions from %d. Exiting.\n", m, tot);
    exit(0);
  }

  ptn=(char *)malloc((size_t) ((2*r1+1)*sizeof(char)));
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<r1;i++)
    ptn[i]='1';
  for(i=r1;i<r1+n;i++)
    ptn[i]='2';
  for(i=r1+n;i<2*r1;i++)
    ptn[i]='3';
  ptn[i]=0;

  //To store the di6 string
  //Set the initial and final values of the graph string
  //These do not change as long as the dimension is constant

  out[0]=38;
  out[1]=63+2*r1;
  out[2+clen]=0;

  //clear files and then open lots of temporary files
  system("rm -f tempfiles/__tmp*");

  for(inv=0;inv<fnum;inv++){//each value of the invariant
    //    fprintf(stderr, "(reverse): %d/%d\n", fnum-inv, fnum);
    for(j=m;j<=4*m;j++){//open files (total degree at least m)
      sprintf(current_file, "tempfiles/__tmp%d_%d_0", inv, j); 
      if(!(fd1[inv][j-m]=fopen(current_file, "w"))){
	fprintf(stderr, "ERROR in genbimol: could not open file \"%s\" for writing.\n", current_file);exit(0);
      }
    }
  }

  flag=1;
  firstcomb(xc, tot, m);
  while(flag){
    //construct left and right complexes
    for(i=0;i<m;i++){//Each complex pair
      for(j=0;j<n;j++){//Each species
	if(!rev){
	  L[j][i]=H[Rvec[xc[i]][0]][j];
	  R[j][i]=H[Rvec[xc[i]][1]][j];
	  Gt[i][j]=R[j][i]-L[j][i];//transpose of stoichiometric matrix
	}
	else{
	  L[j][2*i]=H[Rvec[xc[i]][0]][j];
	  L[j][2*i+1]=H[Rvec[xc[i]][1]][j];
	  R[j][2*i]=H[Rvec[xc[i]][1]][j];
	  R[j][2*i+1]=H[Rvec[xc[i]][0]][j];
	  Gt[2*i][j]=R[j][2*i]-L[j][2*i];//transpose of stoichiometric matrix
	  Gt[2*i+1][j]=R[j][2*i+1]-L[j][2*i+1];
	}
      }
    }

    if(goodorder(L,R,n,m0)){//degrees in ascending order
      alltot++;
      if(alltot==1000000){//some sense of progress
	fprintf(stderr, "*");
	alltot=0;
      }

      //processing
      if((!strstr(procstr, "connect") || connected(L, R, n, m0)) && (!strstr(procstr, "genuine") || !has_triv_spec(L, R, n, m0)) && (rev || !strstr(procstr, "DN") || !hasposimvec(Gt,m,n)) && (rev || !strstr(procstr, "WR") || CRNweakrev2(Rvec, xc, m))){

	//compute the invariants
	num2l=0;num2r=0;//# of times stoichiometry "2" appears on LHS and RHS
	totdeg=0;//total degree
	for(i=0;i<m;i++){
	  if(Rvec[xc[i]][0]){//Not the zero complex
	    if(Rvec[xc[i]][0]<=n)//LH complex is unimolecular: X_i
	      totdeg++;
	    else if(Rvec[xc[i]][0]<=2*n){//LH complex = 2X_i
	      totdeg+=2;
	      num2l++;
	    }
	    else//LH complex is X_i+X_j
	      totdeg+=2;
	  }
	  if(Rvec[xc[i]][1]){
	    if(Rvec[xc[i]][1]<=n)//RH complex is unimolecular: X_i
	      totdeg++;
	    else if(Rvec[xc[i]][1]<=2*n){//RH complex = 2X_i
	      totdeg+=2;
	      num2r++;
	    }
	    else//RH complex is X_i+X_j
	      totdeg+=2;
	  }
	}
	if(rev)
	  inv=num2r+num2l;
	else
	  inv=(m+1)*num2r+num2l;
	//End of compute invariants


	CRNPN2(L,R,n,m0,V,Vlen);
	amtodig(V,clen,out);
	fprintf(fd1[inv][totdeg-m], "%s\n", out);
	(fsz[inv][totdeg-m])++;

	//If (part)-file got large, then (i) wait from a free process; 
	//(ii) shut (part)-file; (iii) short the (part)-file; open next part;
	if(fsz[inv][totdeg-m]>2000000){
	  fsz[inv][totdeg-m]=0;
	  /* fprintf(stderr, "fptempflag=%d\n", fptempflag); */

	  //wait
	  while(1){if(numshortg()<maxproc){break;}sleep(1);}

	  //finish last temporary shortg
	  if(fptempflag==maxproc){// all running, close one
	    pclose(fptemp[oldest_short]);
	    //fptemp[oldest_short]=NULL;
	    fptempflag--;
	  }

	  sprintf(current_file, "tempfiles/__tmp%d_%d", inv, totdeg);
	  fclose(fd1[inv][totdeg-m]);

	  sprintf(systemcom, "nauty-shortg -q -f%s %s_%d", ptn, current_file, part[inv][totdeg-m]);
	  if(debug){fprintf(stderr, "%s\n", systemcom);}
	  if(!(fptemp[oldest_short]=popen(systemcom, "r"))){
	    perror("couldn't open pipe.\n");exit(0);
	  }
	  (part[inv][totdeg-m])++;

	  if(oldest_short<maxproc-1)//increment available counter for temporary shortg processes
	    oldest_short++;
	  else//used all four
	    oldest_short=0;
	  if(debug){fprintf(stderr, "oldest_short = %d\n", oldest_short);}

	  //fptemp[oldest_short]=shortfilepart(current_file, ptn, &(part[inv][totdeg-m]));
	  fptempflag++;

	  // open next temporary file
	  sprintf(systemcom, "%s_%d", current_file, part[inv][totdeg-m]);
	  if(debug){fprintf(stderr, "opening %s\n", systemcom);}
	  if(!(fd1[inv][totdeg-m]=fopen(systemcom, "w"))){
	    fprintf(stderr, "ERROR in shortfile: could not open file \"%s\" for writing.\n", systemcom);exit(0);
	  }

	}//End of file got too large

	//fprintvec(fd, V, 4*r1*r1);
	//milcount++;
      }//End of check for special conditions

    }//End of "if goodorder..."

    flag=nextcomb(xc, tot, m);
  }

  fprintf(stderr, "\n");
  //close temporary shortgs
  if(debug){fprintf(stderr, "closing temporary shortgs\n");
    fprintf(stderr, "number of open shortgs=%d\n", fptempflag);}
  while(fptempflag){
    pclose(fptemp[oldest_short]);
    if(!oldest_short)
      oldest_short=maxproc-1;
    else
      oldest_short--;
    fptempflag--;
  }

  // Final merge and short all the files
  if(debug){fprintf(stderr, "merge and short all files\n");}
  for(i=0;i<fnum;i++){
    for(j=m;j<=4*m;j++){
      fclose(fd1[i][j-m]);

      //wait
      while(1){if(numshortg()<maxproc){break;}sleep(1);}

      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);
      fp[i][j-m]=mergeandshort(current_file, ptn, part[i][j-m]+1);
    }
  }
      
  //close the final streams and join the files
  if(debug){fprintf(stderr, "Joining and sorting files...\n");}
  for(i=0;i<fnum;i++){//each value of the invariant
    for(j=m;j<=4*m;j++){
      pclose(fp[i][j-m]);
      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);  
      sprintf(systemcom, "cat %s >> tempfiles/_tmp2 && rm %s", current_file, current_file);
      system(systemcom);
    }
  }

  sprintf(systemcom, "LC_COLLATE=C sort tempfiles/_tmp2 > %s && rm tempfiles/_tmp2", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  system(systemcom);

  //risky in terms of memory?
  if(debug==2){
    array=genarraydat(outfile, &tot1);
    for(i=0;i<tot1;i++){
      str=di6toreacstr(array[i], n,m,0);
      fprintf(stderr, "********\n%s\n", str);
      free(str);
    }
    freearraydat(array, tot1);
  }

  free(ptn);
  free_imatrix(H,0,N-1,0,n-1);
  free_imatrix(L,0,n-1,0,m0-1);
  free_imatrix(R,0,n-1,0,m0-1);
  free_imatrix(Gt,0,m0-1,0,n-1);

  sprintf(systemcom, "wc -l %s", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp0=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp0);
  totnetworks=atol(path);
  pclose(fp0);

  printf("Generated %ld networks, and stored them in di6 format in %s.\n", totnetworks, outfile);
  return totnetworks;
}

// Generate CRNs with at most bimolecular reactant complexes and 
// at most trimolecular product complexes, on n species and m reactions
// the openswtch decides whether we allow 0 --> X, X --> 0, 
// type reactions or not. 
// open means create only fully open CRNs
// creates lots of temporary files: there must be a directory \"tempfiles\"
unsigned long genbitrimol(int n, int m, char *outfile, int open, const char procstr[], int debug){
  int i,j,k1,k2,tot=0,tot1=0,flag=1;
  //long milcount=0;
  long alltot=0;
  int m0=m;
  int r1=n+m0;
  int r2=4*r1*r1;//minimum dimension of V
  int N1=1+n+n*(n+1)/2;
  int N2=N1 + n*(n+1)*(n+2)/6; //last=n+2 choose 3
  int Rvec[N1*N2][2];//all reactions (LHS index, RHS index, overcount)
  int **H1=imatrix(0, N1-1, 0, n-1);//at most bi complexes
  int **H2=imatrix(0, N2-1, 0, n-1);//at most tri complexes
  int H1deg[N1];
  int H2deg[N2];
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  int clen=Vlen/6;
  char out[clen+3];
  int fnum=(m+1)*(m+1);//equiv classes
  FILE *fd, *fp[fnum][3*m+1];
  FILE *fptemp[4]; // for temporary shortg processes
  FILE *fd1[fnum][3*m+1];
  int oldest_short=0;
  int fptempflag=0;
  char *ptn, *str;
  char **array;
  char systemcom[1024];
  int **L=imatrix(0,n-1,0,m0-1);
  int **R=imatrix(0,n-1,0,m0-1);
  //Transpose of irreversible stoich. mat.
  int **Gt=imatrix(0,m0-1,0,n-1); 
  int xc[m];
  //Total degree between m and 4*m
  int num2l, num2r, inv, totdeg;
  char current_file[30];
  int fsz[fnum][3*m+1];
  int part[fnum][3*m+1];
  int maxproc=4;//maximum shortg processes to run simultaneously
  FILE *fp0;
  char path[1024];
  unsigned long totnetworks;

  //clear file
  if(!(fd=fopen("tempfiles/_tmp2", "w"))){
    fprintf(stderr, "ERROR in genbitrimol: could not open file \"tempfiles/_tmp2\" for writing.\n");exit(0);
  }
  fclose(fd);

  if(r1>30){
    fprintf(stderr, "ERROR in genbitrimol: not sure we can handle the size. Exiting.\n");
    exit(0);
  }

  //initialise
  for(i=0;i<fnum;i++){
    for(j=0;j<=3*m;j++){
      fsz[i][j]=0;part[i][j]=0;
    }
  }
  k1=genallbicomplexes(n, N1, H1);
  k2=genalltricomplexes(n, N2, H2);
  if(debug){fprintf(stderr, "total bimol/trimol complexes = %d, %d\n", k1,k2);}

  for(i=0;i<k1;i++)
    H1deg[i]=rowsum(H1,k1,n,i);
  for(i=0;i<k2;i++)
    H2deg[i]=rowsum(H2,k2,n,i);

  //generate the reactions as index-pairs
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(j!=i && (!open || (i!=0 && j!=0) || (i==0 && j>n) || (j==0 && i>n))){
      	Rvec[tot][0]=i;Rvec[tot][1]=j;
      	tot++;
      }
    }
  }
  if(debug){fprintf(stderr, "total reacs = %d\n", tot);}
  if(m>tot){
    fprintf(stderr, "Cannot choose %d reactions from %d. Exiting.\n", m, tot);
    exit(0);
  }

  ptn=(char *)malloc((size_t) ((2*r1+1)*sizeof(char)));
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<r1;i++)
    ptn[i]='1';
  for(i=r1;i<r1+n;i++)
    ptn[i]='2';
  for(i=r1+n;i<2*r1;i++)
    ptn[i]='3';
  ptn[i]=0;

  //To store the di6 string
  //Set the initial and final values of the graph string
  //These do not change as long as the dimension is constant

  out[0]=38;
  out[1]=63+2*r1;
  out[2+clen]=0;

  //clear files and then open lots of temporary files
  system("rm -f tempfiles/__tmp*");

  for(inv=0;inv<fnum;inv++){//each value of the invariant
    //    fprintf(stderr, "(reverse): %d/%d\n", fnum-inv, fnum);
    for(j=m;j<=4*m;j++){//open files (total degree at least m)
      sprintf(current_file, "tempfiles/__tmp%d_%d_0", inv, j); 
      if(!(fd1[inv][j-m]=fopen(current_file, "w"))){
	fprintf(stderr, "ERROR in genbitrimol: could not open file \"%s\" for writing.\n", current_file);exit(0);
      }
    }
  }

  flag=1;
  firstcomb(xc, tot, m);
  while(flag){
    //construct left and right complexes
    for(i=0;i<m;i++){//Each complex pair
      for(j=0;j<n;j++){//Each species
	L[j][i]=H1[Rvec[xc[i]][0]][j];
	R[j][i]=H2[Rvec[xc[i]][1]][j];
	Gt[i][j]=R[j][i]-L[j][i];//transpose of stoichiometric matrix
      }
    }

    if(goodorder(L,R,n,m0)){//degrees in ascending order
      alltot++;
      if(alltot==1000000){//some sense of progress
	fprintf(stderr, "*");
	alltot=0;
      }

      //processing
      if((!strstr(procstr, "connect") || connected(L, R, n, m0)) && (!strstr(procstr, "genuine") || !has_triv_spec(L, R, n, m0)) && (!strstr(procstr, "DN") || !hasposimvec(Gt,m,n)) ){

	//compute the invariants
	num2l=0;num2r=0;//# of times stoichiometry "2" appears on LHS and RHS
	totdeg=0;//total degree
	for(i=0;i<m;i++){
	  if(Rvec[xc[i]][0]>n && Rvec[xc[i]][0]<=2*n){//LH complex = 2X_i
	    num2l++;
	  }
	  totdeg+=H1deg[Rvec[xc[i]][0]];
   
	  if(Rvec[xc[i]][1]>n && Rvec[xc[i]][1]<=2*n){//RH complex = 2X_i
	    num2r++;
	  }
	  totdeg+=H2deg[Rvec[xc[i]][1]];
	}
	inv=(m+1)*num2r+num2l;
	//End of compute invariants


	CRNPN2(L,R,n,m0,V,Vlen);
	amtodig(V,clen,out);
	fprintf(fd1[inv][totdeg-m], "%s\n", out);
	(fsz[inv][totdeg-m])++;

	//If (part)-file got large, then (i) wait from a free process; 
	//(ii) shut (part)-file; (iii) short the (part)-file; open next part;
	if(fsz[inv][totdeg-m]>2000000){
	  fsz[inv][totdeg-m]=0;
	  /* fprintf(stderr, "fptempflag=%d\n", fptempflag); */

	  //wait
	  while(1){if(numshortg()<maxproc){break;}sleep(1);}

	  //finish last temporary shortg
	  if(fptempflag==maxproc){// all running, close one
	    pclose(fptemp[oldest_short]);
	    //fptemp[oldest_short]=NULL;
	    fptempflag--;
	  }

	  sprintf(current_file, "tempfiles/__tmp%d_%d", inv, totdeg);
	  fclose(fd1[inv][totdeg-m]);

	  sprintf(systemcom, "nauty-shortg -q -f%s %s_%d", ptn, current_file, part[inv][totdeg-m]);
	  if(debug){fprintf(stderr, "%s\n", systemcom);}
	  if(!(fptemp[oldest_short]=popen(systemcom, "r"))){
	    perror("couldn't open pipe.\n");exit(0);
	  }
	  (part[inv][totdeg-m])++;

	  if(oldest_short<maxproc-1)//increment available counter for temporary shortg processes
	    oldest_short++;
	  else//used all four
	    oldest_short=0;
	  if(debug){fprintf(stderr, "oldest_short = %d\n", oldest_short);}

	  //fptemp[oldest_short]=shortfilepart(current_file, ptn, &(part[inv][totdeg-m]));
	  fptempflag++;

	  // open next temporary file
	  sprintf(systemcom, "%s_%d", current_file, part[inv][totdeg-m]);
	  if(debug){fprintf(stderr, "opening %s\n", systemcom);}
	  if(!(fd1[inv][totdeg-m]=fopen(systemcom, "w"))){
	    fprintf(stderr, "ERROR in shortfile: could not open file \"%s\" for writing.\n", systemcom);exit(0);
	  }

	}//End of file got too large

	//fprintvec(fd, V, 4*r1*r1);
	//milcount++;
      }//End of check for special conditions

    }//End of "if goodorder..."

    flag=nextcomb(xc, tot, m);
  }

  fprintf(stderr, "\n");
  //close temporary shortgs
  if(debug){fprintf(stderr, "closing temporary shortgs\n");
    fprintf(stderr, "number of open shortgs=%d\n", fptempflag);}
  while(fptempflag){
    pclose(fptemp[oldest_short]);
    if(!oldest_short)
      oldest_short=maxproc-1;
    else
      oldest_short--;
    fptempflag--;
  }

  // Final merge and short all the files
  if(debug){fprintf(stderr, "merge and short all files\n");}
  for(i=0;i<fnum;i++){
    for(j=m;j<=4*m;j++){
      fclose(fd1[i][j-m]);

      //wait
      while(1){if(numshortg()<maxproc){break;}sleep(1);}

      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);
      fp[i][j-m]=mergeandshort(current_file, ptn, part[i][j-m]+1);
    }
  }
      
  //close the final streams and join the files
  if(debug){fprintf(stderr, "Joining and sorting files...\n");}
  for(i=0;i<fnum;i++){//each value of the invariant
    for(j=m;j<=4*m;j++){
      pclose(fp[i][j-m]);
      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);  
      sprintf(systemcom, "cat %s >> tempfiles/_tmp2 && rm %s", current_file, current_file);
      system(systemcom);
    }
  }

  sprintf(systemcom, "LC_COLLATE=C sort tempfiles/_tmp2 > %s && rm tempfiles/_tmp2", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  system(systemcom);

  //risky in terms of memory?
  if(debug==2){
    array=genarraydat(outfile, &tot1);
    for(i=0;i<tot1;i++){
      str=di6toreacstr(array[i], n,m,0);
      fprintf(stderr, "********\n%s\n", str);
      free(str);
    }
    freearraydat(array, tot1);
  }

  free(ptn);
  free_imatrix(H1,0,N1-1,0,n-1);
  free_imatrix(H2,0,N2-1,0,n-1);
  free_imatrix(L,0,n-1,0,m0-1);
  free_imatrix(R,0,n-1,0,m0-1);
  free_imatrix(Gt,0,m0-1,0,n-1);

  sprintf(systemcom, "wc -l %s", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp0=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp0);
  totnetworks=atol(path);
  pclose(fp0);

  printf("Generated %ld networks, and stored them in di6 format in %s.\n", totnetworks, outfile);
  return totnetworks;
}


// Generate (0,1) CRNs with at most trimolecular complexes
// on n species and m reactions
// the openswtch decides whether we allow 0 --> X, X --> 0, 
// type reactions or not. 
// open means create only fully open CRNs
// creates lots of temporary files: there must be a directory \"tempfiles\"
unsigned long genzeroonetritrimol(int n, int m, char *outfile, int open, const char procstr[], int debug){
  int i,j,k1,tot=0,tot1=0,flag=1;
  //long milcount=0;
  long alltot=0;
  int m0=m;
  int r1=n+m0;
  int r2=4*r1*r1;//minimum dimension of V
  int N1=1+n+n*(n-1)/2 + n*(n-1)*(n-2)/6;
  int Rvec[N1*N1][2];//all reactions (LHS index, RHS index, overcount)
  int **H1=imatrix(0, N1-1, 0, n-1);//at most tri (0,1) complexes
  int H1deg[N1];
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  int clen=Vlen/6;
  char out[clen+3];
  int totl=3*m+1;//total left molecularity 0--3*m
  int totr=3*m+1;//total right molecularity 0--3*m
  int degl,degr;
  FILE *fd, *fp[totl][totr];
  FILE *fptemp[4]; // for temporary shortg processes
  FILE *fd1[totl][totr];
  int oldest_short=0;
  int fptempflag=0;
  char *ptn, *str;
  char **array;
  char systemcom[1024];
  int **L=imatrix(0,n-1,0,m0-1);
  int **R=imatrix(0,n-1,0,m0-1);
  //Transpose of irreversible stoich. mat.
  int **Gt=imatrix(0,m0-1,0,n-1); 
  int xc[m];
  int inv;
  char current_file[30];
  int fsz[totl][totr];
  int part[totl][totr];
  int maxproc=4;//maximum shortg processes to run simultaneously
  FILE *fp0;
  char path[1024];
  unsigned long totnetworks;

  //clear file
  if(!(fd=fopen("tempfiles/_tmp2", "w"))){
    fprintf(stderr, "ERROR in genzeroonetritrimol: could not open file \"tempfiles/_tmp2\" for writing.\n");exit(0);
  }
  fclose(fd);

  if(r1>30){
    fprintf(stderr, "ERROR in genzeroonetritrimol: not sure we can handle the size. Exiting.\n");
    exit(0);
  }

  //initialise
  for(i=0;i<totl;i++){
    for(j=0;j<totr;j++){
      fsz[i][j]=0;part[i][j]=0;
    }
  }
  k1=genallzeroonetricomplexes(n, N1, H1);
  if(debug){fprintf(stderr, "total trimol complexes = %d\n", k1);}

  for(i=0;i<k1;i++)
    H1deg[i]=rowsum(H1,k1,n,i);

  //generate the reactions as index-pairs
  for(i=0;i<k1;i++){
    for(j=0;j<k1;j++){
      if(j!=i && (!open || (i!=0 && j!=0) || (i==0 && j>n) || (j==0 && i>n))){
      	Rvec[tot][0]=i;Rvec[tot][1]=j;
      	tot++;
      }
    }
  }
  if(debug){fprintf(stderr, "total reacs = %d\n", tot);}
  if(m>tot){
    fprintf(stderr, "Cannot choose %d reactions from %d. Exiting.\n", m, tot);
    exit(0);
  }

  ptn=(char *)malloc((size_t) ((2*r1+1)*sizeof(char)));
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<r1;i++)
    ptn[i]='1';
  for(i=r1;i<r1+n;i++)
    ptn[i]='2';
  for(i=r1+n;i<2*r1;i++)
    ptn[i]='3';
  ptn[i]=0;

  //To store the di6 string
  //Set the initial and final values of the graph string
  //These do not change as long as the dimension is constant

  out[0]=38;
  out[1]=63+2*r1;
  out[2+clen]=0;

  //clear files and then open lots of temporary files
  system("rm -f tempfiles/__tmp*");

  for(inv=0;inv<totl;inv++){//each value of the invariant
    //    fprintf(stderr, "(reverse): %d/%d\n", totl-inv, totl);
    for(j=0;j<totr;j++){//open files (total degree at least m)
      sprintf(current_file, "tempfiles/__tmp%d_%d_0", inv, j); 
      if(!(fd1[inv][j]=fopen(current_file, "w"))){
	fprintf(stderr, "ERROR in genzeroonetritrimol: could not open file \"%s\" for writing.\n", current_file);exit(0);
      }
    }
  }

  flag=1;
  firstcomb(xc, tot, m);
  while(flag){
    //construct left and right complexes
    for(i=0;i<m;i++){//Each complex pair
      for(j=0;j<n;j++){//Each species
	L[j][i]=H1[Rvec[xc[i]][0]][j];
	R[j][i]=H1[Rvec[xc[i]][1]][j];
	Gt[i][j]=R[j][i]-L[j][i];//transpose of stoichiometric matrix
      }
    }

    if(goodorder(L,R,n,m0)){//degrees in ascending order
      alltot++;
      if(alltot==1000000){//some sense of progress
	fprintf(stderr, "*");
	alltot=0;
      }

      //processing
      if((!strstr(procstr, "connect") || connected(L, R, n, m0)) && (!strstr(procstr, "genuine") || !has_triv_spec(L, R, n, m0)) && (!strstr(procstr, "DN") || !hasposimvec(Gt,m,n)) ){

	//compute the invariants
	degl=0;
	degr=0;
	for(i=0;i<m;i++){
	  degl+=H1deg[Rvec[xc[i]][0]];
   	  degr+=H1deg[Rvec[xc[i]][1]];
	}
	//End of compute invariants


	CRNPN2(L,R,n,m0,V,Vlen);
	amtodig(V,clen,out);
	fprintf(fd1[degl][degr], "%s\n", out);
	(fsz[degl][degr])++;

	//If (part)-file got large, then (i) wait from a free process; 
	//(ii) shut (part)-file; (iii) short the (part)-file; open next part;
	if(fsz[degl][degr]>2000000){
	  fsz[degl][degr]=0;
	  /* fprintf(stderr, "fptempflag=%d\n", fptempflag); */

	  //wait
	  while(1){if(numshortg()<maxproc){break;}sleep(1);}

	  //finish last temporary shortg
	  if(fptempflag==maxproc){// all running, close one
	    pclose(fptemp[oldest_short]);
	    //fptemp[oldest_short]=NULL;
	    fptempflag--;
	  }

	  sprintf(current_file, "tempfiles/__tmp%d_%d", degl, degr);
	  fclose(fd1[degl][degr]);

	  sprintf(systemcom, "nauty-shortg -q -f%s %s_%d", ptn, current_file, part[degl][degr]);
	  if(debug){fprintf(stderr, "%s\n", systemcom);}
	  if(!(fptemp[oldest_short]=popen(systemcom, "r"))){
	    perror("couldn't open pipe.\n");exit(0);
	  }
	  (part[degl][degr])++;

	  if(oldest_short<maxproc-1)//increment available counter for temporary shortg processes
	    oldest_short++;
	  else//used all four
	    oldest_short=0;
	  if(debug){fprintf(stderr, "oldest_short = %d\n", oldest_short);}

	  fptempflag++;

	  // open next temporary file
	  sprintf(systemcom, "%s_%d", current_file, part[degl][degr]);
	  if(debug){fprintf(stderr, "opening %s\n", systemcom);}
	  if(!(fd1[degl][degr]=fopen(systemcom, "w"))){
	    fprintf(stderr, "ERROR in shortfile: could not open file \"%s\" for writing.\n", systemcom);exit(0);
	  }

	}//End of file got too large

	//fprintvec(fd, V, 4*r1*r1);
	//milcount++;
      }//End of check for special conditions

    }//End of "if goodorder..."

    flag=nextcomb(xc, tot, m);
  }

  fprintf(stderr, "\n");
  //close temporary shortgs
  if(debug){fprintf(stderr, "closing temporary shortgs\n");
    fprintf(stderr, "number of open shortgs=%d\n", fptempflag);}
  while(fptempflag){
    pclose(fptemp[oldest_short]);
    if(!oldest_short)
      oldest_short=maxproc-1;
    else
      oldest_short--;
    fptempflag--;
  }

  // Final merge and short all the files
  if(debug){fprintf(stderr, "merge and short all files\n");}
  for(i=0;i<totl;i++){
    for(j=0;j<totr;j++){
      fclose(fd1[i][j]);

      //wait
      while(1){if(numshortg()<maxproc){break;}sleep(1);}
      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);
      fp[i][j]=mergeandshort(current_file, ptn, part[i][j]+1);
    }
  }
  //close the final streams and join the files
  if(debug){fprintf(stderr, "Joining and sorting files...\n");}
  for(i=0;i<totl;i++){//each value of the invariant
    for(j=0;j<totr;j++){
      pclose(fp[i][j]);
      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);  
      sprintf(systemcom, "cat %s >> tempfiles/_tmp2 && rm %s", current_file, current_file);
      system(systemcom);
    }
  }

  sprintf(systemcom, "LC_COLLATE=C sort tempfiles/_tmp2 > %s && rm tempfiles/_tmp2", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  system(systemcom);

  //risky in terms of memory?
  if(debug==2){
    array=genarraydat(outfile, &tot1);
    for(i=0;i<tot1;i++){
      str=di6toreacstr(array[i], n,m,0);
      fprintf(stderr, "********\n%s\n", str);
      free(str);
    }
    freearraydat(array, tot1);
  }

  free(ptn);
  free_imatrix(H1,0,N1-1,0,n-1);
  free_imatrix(L,0,n-1,0,m0-1);
  free_imatrix(R,0,n-1,0,m0-1);
  free_imatrix(Gt,0,m0-1,0,n-1);

  sprintf(systemcom, "wc -l %s", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp0=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp0);
  totnetworks=atol(path);
  pclose(fp0);

  printf("Generated %ld networks, and stored them in di6 format in %s.\n", totnetworks, outfile);
  return totnetworks;
}


//n+k choose k
int nkchoosek(int n, int k){
  int j;
  int t=1;
  for(j=0;j<k;j++){
    t*=(n+j)/(j+1);
  }
  return t;
}

// Experimental, and not checked
// Generate CRNs with at most bimolecular reactant complexes and 
// at most k-molecular product complexes ("molecularity" = k<=7 or update layers), 
// on n species and m reactions
// the openswtch decides whether we allow 0 --> X, X --> 0, 
// type reactions or not. "open" means create only fully open CRNs
// creates lots of temporary files: there must be a directory \"tempfiles\"
unsigned long genbimultimol(int n, int m, char *outfile, int open, const char procstr[], int molecularity, int debug){
  int i,j,k1,k2,tot=0,tot1=0,flag=1;
  //long milcount=0;
  long alltot=0;
  int m0=m;
  int r1=n+m0;
  //int layers=3;
  int layers=(int)(log((double)molecularity)/log(2.0)+1.1);//to be safe
  int r2=layers*layers*r1*r1;//minimum dimension of V: l^2*(n+m)*(n+m)
  int N1=1+n+n*(n+1)/2;
  int N2=nkchoosek(n+molecularity, molecularity);
  //fprintf(stderr, "N2=%d\n",N2);exit(0);
  int Rvec[N1*N2][2];//all reactions (LHS index, RHS index, overcount)
  int **H1=imatrix(0, N1-1, 0, n-1);//at most bi complexes
  int **H2=imatrix(0, N2-1, 0, n-1);//at most tetra complexes
  int H1deg[N1];
  int H2deg[N2];
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool *V;
  int clen=Vlen/6;
  char out[clen+3];
  int fnum=m*(1+molecularity/2)+1;//total number of twos
//(m+1)*(m*molecularity/2+1)+1;//number of 2s on each side
  int totmol=m*(2+molecularity)+1;//total molecularity
  FILE *fd, *fp[fnum][totmol];
  FILE *fptemp[4]; // for temporary shortg processes
  FILE *fd1[fnum][totmol];
  int oldest_short=0;
  int fptempflag=0;
  char *ptn, *str;
  char **array;
  char systemcom[1024];
  int **L=imatrix(0,n-1,0,m0-1);
  int **R=imatrix(0,n-1,0,m0-1);
  //Transpose of irreversible stoich. mat.
  int **Gt=imatrix(0,m0-1,0,n-1); 
  int xc[m];
  int num2l, num2r, inv, totdeg;
  char current_file[30];
  int fsz[fnum][totmol];
  int part[fnum][totmol];
  int maxproc=4;//maximum shortg processes to run simultaneously
  FILE *fp0;
  char path[1024];
  unsigned long totnetworks;
  int l;

  //clear file
  if(!(fd=fopen("tempfiles/_tmp2", "w"))){
    fprintf(stderr, "ERROR in genbimultimol: could not open file \"tempfiles/_tmp2\" for writing.\n");exit(0);
  }
  fclose(fd);

  if(r1>30){
    fprintf(stderr, "ERROR in genbimultimol: not sure we can handle the size. Exiting.\n");
    exit(0);
  }

  //initialise
  for(i=0;i<fnum;i++){
    for(j=0;j<totmol;j++){
      fsz[i][j]=0;part[i][j]=0;
    }
  }
  k1=genallbicomplexes(n, N1, H1);
  k2=genallmulticomplexes(n, N2, H2, molecularity);
  if(debug){fprintf(stderr, "total bimol/multimol complexes = %d, %d\n", k1,k2);}

  for(i=0;i<k1;i++)
    H1deg[i]=rowsum(H1,k1,n,i);
  for(i=0;i<k2;i++)
    H2deg[i]=rowsum(H2,k2,n,i);
  //generate the reactions as index-pairs
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(j!=i && (!open || (i!=0 && j!=0) || (i==0 && j>n) || (j==0 && i>n))){
      	Rvec[tot][0]=i;Rvec[tot][1]=j;
      	tot++;
      }
    }
  }
  if(debug){fprintf(stderr, "total reacs = %d\n", tot);}
  if(m>tot){
    fprintf(stderr, "Cannot choose %d reactions from %d. Exiting.\n", m, tot);
    exit(0);
  }

  ptn=(char *)malloc((size_t) ((2*r1+1)*sizeof(char)));
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<r1;i++)
    ptn[i]='1';
  for(i=r1;i<r1+n;i++)
    ptn[i]='2';
  for(i=r1+n;i<2*r1;i++)
    ptn[i]='3';
  ptn[i]=0;

  //To store the di6 string
  //Set the initial and final values of the graph string
  //These do not change as long as the dimension is constant

  out[0]=38;
  out[1]=63+layers*r1;
  out[2+clen]=0;

  //clear files and then open lots of temporary files
  system("rm -f tempfiles/__tmp*");

  for(inv=0;inv<fnum;inv++){//each value of the invariant
    for(j=0;j<totmol;j++){//open files (total degree at least m)
      sprintf(current_file, "tempfiles/__tmp%d_%d_0", inv, j); 
      if(!(fd1[inv][j]=fopen(current_file, "w"))){
	fprintf(stderr, "ERROR in genbimultimol: could not open file \"%s\" for writing.\n", current_file);exit(0);
      }
    }
  }


  flag=1;
  firstcomb(xc, tot, m);
  while(flag){
    //construct left and right complexes
    for(i=0;i<m;i++){//Each complex pair
      for(j=0;j<n;j++){//Each species
	L[j][i]=H1[Rvec[xc[i]][0]][j];
	R[j][i]=H2[Rvec[xc[i]][1]][j];
	Gt[i][j]=R[j][i]-L[j][i];//transpose of stoichiometric matrix
      }
    }

    if(goodorder(L,R,n,m0)){//degrees in ascending order
      alltot++;
      if(alltot==1000000){//some sense of progress
	fprintf(stderr, "*");
	alltot=0;
      }

      //processing
      if((!strstr(procstr, "connect") || connected(L, R, n, m0)) && (!strstr(procstr, "genuine") || !has_triv_spec(L, R, n, m0)) && (!strstr(procstr, "DN") || !hasposimvec(Gt,m,n)) ){

	//compute the invariants
	num2l=0;num2r=0;//# of times stoichiometry "2" appears on both sides
	totdeg=0;//total degree
	for(i=0;i<m;i++){
	  for(j=0;j<n;j++){
	    if(L[j][i]==2)
	      num2l++;
	    if(R[j][i]==2)
	      num2r++;
	  }
	  totdeg+=H1deg[Rvec[xc[i]][0]];
	  totdeg+=H2deg[Rvec[xc[i]][1]];
	}
	inv=num2l+num2r;
	//End of compute invariants
	V=CRNPN3(L,R,n,m0,layers,&l);
	amtodig(V,clen,out);
	free((char *)V);

	fprintf(fd1[inv][totdeg], "%s\n", out);
	(fsz[inv][totdeg])++;

	//If (part)-file got large, then (i) wait from a free process; 
	//(ii) shut (part)-file; (iii) short the (part)-file; open next part;
	if(fsz[inv][totdeg]>2000000){
	  fsz[inv][totdeg]=0;
	  /* fprintf(stderr, "fptempflag=%d\n", fptempflag); */

	  //wait
	  while(1){if(numshortg()<maxproc){break;}sleep(1);}

	  //finish last temporary shortg
	  if(fptempflag==maxproc){// all running, close one
	    pclose(fptemp[oldest_short]);
	    //fptemp[oldest_short]=NULL;
	    fptempflag--;
	  }

	  sprintf(current_file, "tempfiles/__tmp%d_%d", inv, totdeg);
	  fclose(fd1[inv][totdeg]);

	  sprintf(systemcom, "nauty-shortg -q -f%s %s_%d", ptn, current_file, part[inv][totdeg]);
	  if(debug){fprintf(stderr, "%s\n", systemcom);}
	  if(!(fptemp[oldest_short]=popen(systemcom, "r"))){
	    perror("couldn't open pipe.\n");exit(0);
	  }
	  (part[inv][totdeg])++;

	  if(oldest_short<maxproc-1)//increment available counter for temporary shortg processes
	    oldest_short++;
	  else//used all four
	    oldest_short=0;
	  if(debug){fprintf(stderr, "oldest_short = %d\n", oldest_short);}

	  //fptemp[oldest_short]=shortfilepart(current_file, ptn, &(part[inv][totdeg-m]));
	  fptempflag++;

	  // open next temporary file
	  sprintf(systemcom, "%s_%d", current_file, part[inv][totdeg]);
	  if(debug){fprintf(stderr, "opening %s\n", systemcom);}
	  if(!(fd1[inv][totdeg]=fopen(systemcom, "w"))){
	    fprintf(stderr, "ERROR in shortfile: could not open file \"%s\" for writing.\n", systemcom);exit(0);
	  }

	}//End of file got too large

	//fprintvec(fd, V, 4*r1*r1);
	//milcount++;
      }//End of check for special conditions

    }//End of "if goodorder..."

    flag=nextcomb(xc, tot, m);
  }

  fprintf(stderr, "\n");
  //close temporary shortgs
  if(debug){fprintf(stderr, "closing temporary shortgs\n");
    fprintf(stderr, "number of open shortgs=%d\n", fptempflag);}
  while(fptempflag){
    pclose(fptemp[oldest_short]);
    if(!oldest_short)
      oldest_short=maxproc-1;
    else
      oldest_short--;
    fptempflag--;
  }

  // Final merge and short all the files
  if(debug){fprintf(stderr, "merge and short all files\n");}
  for(i=0;i<fnum;i++){
    for(j=0;j<totmol;j++){
      fclose(fd1[i][j]);

      //wait
      while(1){if(numshortg()<maxproc){break;}sleep(1);}

      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);
      fp[i][j]=mergeandshort(current_file, ptn, part[i][j]+1);
    }
  }
      
  //close the final streams and join the files
  if(debug){fprintf(stderr, "Joining and sorting files...\n");}
  for(i=0;i<fnum;i++){//each value of the invariant
    for(j=0;j<totmol;j++){
      pclose(fp[i][j]);
      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);  
      sprintf(systemcom, "cat %s >> tempfiles/_tmp2 && rm %s", current_file, current_file);
      system(systemcom);
    }
  }

  sprintf(systemcom, "LC_COLLATE=C sort tempfiles/_tmp2 > %s && rm tempfiles/_tmp2", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  system(systemcom);

  //risky in terms of memory?
  if(debug==2){
    array=genarraydat(outfile, &tot1);
    for(i=0;i<tot1;i++){
      str=di6toreacstr(array[i], n,m,0);
      fprintf(stderr, "********\n%s\n", str);
      free(str);
    }
    freearraydat(array, tot1);
  }

  free(ptn);
  free_imatrix(H1,0,N1-1,0,n-1);
  free_imatrix(H2,0,N2-1,0,n-1);
  free_imatrix(L,0,n-1,0,m0-1);
  free_imatrix(R,0,n-1,0,m0-1);
  free_imatrix(Gt,0,m0-1,0,n-1);

  sprintf(systemcom, "wc -l %s", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp0=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp0);
  totnetworks=atol(path);
  pclose(fp0);

  printf("Generated %ld networks, and stored them in di6 format in %s.\n", totnetworks, outfile);
  return totnetworks;
}



// Experimental, and not checked
// Generate source-only CRNs with at most k-molecular sources
// ("molecularity" = k<=7 or update layers), 
// on n species and m reactions
// creates lots of temporary files: there must be a directory \"tempfiles\"
unsigned long genmultisrc(int n, int m, const char outfile[], int molecularity, int debug){
  int i,j,k1,tot=0,tot1=0,flag=1;
  //long milcount=0;
  long alltot=0;
  int m0=m;
  int r1=n+m0;
  //int layers=3;
  int layers=(int)(log((double)molecularity)/log(2.0)+1.1);//to be safe
  int r2=layers*layers*r1*r1;//minimum dimension of V: l^2*(n+m)*(n+m)
  int N1=nkchoosek(n+molecularity, molecularity);
  //fprintf(stderr, "N2=%d\n",N2);exit(0);
  int Rvec[N1][2];//all reactions (LHS index, RHS index, overcount)
  int **H1=imatrix(0, N1-1, 0, n-1);
  int H1deg[N1];
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool *V;
  int clen=Vlen/6;
  char out[clen+3];
  int fnum=m*molecularity/2+1;//number of 2s on LHS
  int totmol=m*molecularity+1;//total molecularity
  FILE *fd, *fp[fnum][totmol];
  FILE *fptemp[4]; // for temporary shortg processes
  FILE *fd1[fnum][totmol];
  int oldest_short=0;
  int fptempflag=0;
  char *ptn, *str;
  char **array;
  char systemcom[1024];
  int **L=imatrix(0,n-1,0,m0-1);
  int **R=imatrix(0,n-1,0,m0-1);
  //Transpose of irreversible stoich. mat.
  int **Gt=imatrix(0,m0-1,0,n-1); 
  int xc[m];
  //Total degree between m and 4*m
  int num2l, inv, totdeg;
  char current_file[30];
  int fsz[fnum][totmol];
  int part[fnum][totmol];
  int maxproc=4;//maximum shortg processes to run simultaneously
  FILE *fp0;
  char path[1024];
  unsigned long totnetworks;
  int l;

  //clear file
  if(!(fd=fopen("tempfiles/_tmp2", "w"))){
    fprintf(stderr, "ERROR in genmultisrc: could not open file \"tempfiles/_tmp2\" for writing.\n");exit(0);
  }
  fclose(fd);

  if(r1>30){
    fprintf(stderr, "ERROR in genmultisrc: not sure we can handle the size. Exiting.\n");
    exit(0);
  }

  //initialise
  for(i=0;i<fnum;i++){
    for(j=0;j<totmol;j++){
      fsz[i][j]=0;part[i][j]=0;
    }
  }
  k1=genallmulticomplexes(n, N1, H1, molecularity);
  if(debug){fprintf(stderr, "total multimol complexes = %d\n", k1);}

  for(i=0;i<k1;i++)
    H1deg[i]=rowsum(H1,k1,n,i);

  //generate the reactions as index-pairs
  for(i=0;i<k1;i++){
    Rvec[tot][0]=i;Rvec[tot][1]=0;
    tot++;
  }
  if(debug){fprintf(stderr, "total reacs = %d\n", tot);}

  ptn=(char *)malloc((size_t) ((2*r1+1)*sizeof(char)));
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<r1;i++)
    ptn[i]='1';
  for(i=r1;i<r1+n;i++)
    ptn[i]='2';
  for(i=r1+n;i<2*r1;i++)
    ptn[i]='3';
  ptn[i]=0;

  //To store the di6 string
  //Set the initial and final values of the graph string
  //These do not change as long as the dimension is constant

  out[0]=38;
  out[1]=63+layers*r1;
  out[2+clen]=0;

  //clear files and then open lots of temporary files
  system("rm -f tempfiles/__tmp*");
  //fprintf(stderr, "trying to open %d files.\n", fnum*(3*m+1));
  for(inv=0;inv<fnum;inv++){//each value of the invariant
    //    fprintf(stderr, "(reverse): %d/%d\n", fnum-inv, fnum);
    for(j=0;j<totmol;j++){//open files (total degree at least m)
      sprintf(current_file, "tempfiles/__tmp%d_%d_0", inv, j); 
      if(!(fd1[inv][j]=fopen(current_file, "w"))){
	fprintf(stderr, "ERROR in genmultisrc: could not open file \"%s\" for writing.\n", current_file);exit(0);
      }
    }
  }

  flag=1;
  inittozero(xc,m);
  //firstcomb(xc, tot, m);
  while(flag){
    //printvec(xc,m);
    //construct left and right complexes
    for(i=0;i<m;i++){//Each complex pair
      for(j=0;j<n;j++){//Each species
	L[j][i]=H1[Rvec[xc[i]][0]][j];
	R[j][i]=H1[Rvec[xc[i]][1]][j];
	Gt[i][j]=R[j][i]-L[j][i];//transpose of stoichiometric matrix
      }
    }

    if(goodorder(L,R,n,m0)){//degrees in ascending order
      alltot++;
      if(alltot==1000000){//some sense of progress
	fprintf(stderr, "*");
	alltot=0;
      }

      //compute the invariants
      num2l=0;//# of times stoichiometry "2" appears on LHS
      totdeg=0;//total degree
      for(i=0;i<m;i++){
	for(j=0;j<n;j++){
	  if(L[j][i]==2)
	    num2l++;
	}
	totdeg+=H1deg[Rvec[xc[i]][0]];
   	totdeg+=H1deg[Rvec[xc[i]][1]];
      }
      inv=num2l;
      //End of compute invariants

      V=CRNPN3(L,R,n,m0,layers,&l);
      amtodig(V,clen,out);
      free((char *)V);
      fprintf(fd1[inv][totdeg], "%s\n", out);
      (fsz[inv][totdeg])++;

      //If (part)-file got large, then (i) wait from a free process; 
      //(ii) shut (part)-file; (iii) short the (part)-file; open next part;
      if(fsz[inv][totdeg]>2000000){
	fsz[inv][totdeg]=0;
	/* fprintf(stderr, "fptempflag=%d\n", fptempflag); */

	//wait
	while(1){if(numshortg()<maxproc){break;}sleep(1);}

	//finish last temporary shortg
	if(fptempflag==maxproc){// all running, close one
	  pclose(fptemp[oldest_short]);
	  //fptemp[oldest_short]=NULL;
	  fptempflag--;
	}

	sprintf(current_file, "tempfiles/__tmp%d_%d", inv, totdeg);
	fclose(fd1[inv][totdeg]);

	sprintf(systemcom, "nauty-shortg -q -f%s %s_%d", ptn, current_file, part[inv][totdeg]);
	if(debug){fprintf(stderr, "%s\n", systemcom);}
	if(!(fptemp[oldest_short]=popen(systemcom, "r"))){
	  perror("couldn't open pipe.\n");exit(0);
	}
	(part[inv][totdeg])++;

	if(oldest_short<maxproc-1)//increment available counter for temporary shortg processes
	  oldest_short++;
	else//used all four
	  oldest_short=0;
	if(debug){fprintf(stderr, "oldest_short = %d\n", oldest_short);}

	//fptemp[oldest_short]=shortfilepart(current_file, ptn, &(part[inv][totdeg]));
	fptempflag++;

	// open next temporary file
	sprintf(systemcom, "%s_%d", current_file, part[inv][totdeg]);
	if(debug){fprintf(stderr, "opening %s\n", systemcom);}
	if(!(fd1[inv][totdeg]=fopen(systemcom, "w"))){
	  fprintf(stderr, "ERROR in shortfile: could not open file \"%s\" for writing.\n", systemcom);exit(0);
	}

      }//End of file got too large

      //fprintvec(fd, V, 4*r1*r1);
      //milcount++;

    }//End of "if goodorder..."
    flag=nextnum(xc, m, tot);

  }

  fprintf(stderr, "\n");
  //close temporary shortgs
  if(debug){fprintf(stderr, "closing temporary shortgs\n");
    fprintf(stderr, "number of open shortgs=%d\n", fptempflag);}
  while(fptempflag){
    pclose(fptemp[oldest_short]);
    if(!oldest_short)
      oldest_short=maxproc-1;
    else
      oldest_short--;
    fptempflag--;
  }

  // Final merge and short all the files
  if(debug){fprintf(stderr, "merge and short all files\n");}
  for(i=0;i<fnum;i++){
    for(j=0;j<totmol;j++){
      fclose(fd1[i][j]);

      //wait
      while(1){if(numshortg()<maxproc){break;}sleep(1);}

      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);
      fp[i][j]=mergeandshort(current_file, ptn, part[i][j]+1);
    }
  }
      
  //close the final streams and join the files
  if(debug){fprintf(stderr, "Joining and sorting files...\n");}
  for(i=0;i<fnum;i++){//each value of the invariant
    for(j=0;j<totmol;j++){
      pclose(fp[i][j]);
      sprintf(current_file, "tempfiles/__tmp%d_%d", i, j);  
      sprintf(systemcom, "cat %s >> tempfiles/_tmp2 && rm %s", current_file, current_file);
      system(systemcom);
    }
  }

  sprintf(systemcom, "LC_COLLATE=C sort tempfiles/_tmp2 > %s && rm tempfiles/_tmp2", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  system(systemcom);

  //risky in terms of memory?
  if(debug==2){
    array=genarraydat(outfile, &tot1);
    for(i=0;i<tot1;i++){
      str=di6toreacstr(array[i], n,m,0);
      fprintf(stderr, "********\n%s\n", str);
      free(str);
    }
    freearraydat(array, tot1);
  }

  free(ptn);
  free_imatrix(H1,0,N1-1,0,n-1);
  free_imatrix(L,0,n-1,0,m0-1);
  free_imatrix(R,0,n-1,0,m0-1);
  free_imatrix(Gt,0,m0-1,0,n-1);

  sprintf(systemcom, "wc -l %s", outfile);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp0=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp0);
  totnetworks=atol(path);
  pclose(fp0);

  printf("Generated %ld networks, and stored them in di6 format in %s.\n", totnetworks, outfile);
  return totnetworks;
}

// (New version, l layers)
// Create a canonical "dynamically isomorphic" CRN (i.e. same MA ODEs) 
// (same source complexes, each column of stoichiometric matrix is 
// the smallest integer vector collinear with original)
// Assuming no more than trimolecular
char *dynamic_isomorphMA(char *di6, int n, int m){
  int entries,i,j;
  int **L=imatrix(0,n-1,0,m-1);//left stoich mat
  int **R=imatrix(0,n-1,0,m-1);//right stoich mat
  int **G=imatrix(0,n-1,0,m-1);//stoich mat
  int totV;
  int **AM=di6toCRNam1(di6, n, m, &totV, &entries);
  int layers=numlayers(di6, n, m);
  int layers2;
  int r1=n+m;
  int r2=layers*layers*r1*r1;//minimum dimension of V
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool *V;
  char *out=(char *)malloc((size_t) ((Vlen/6+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+layers*r1;
  out[2+Vlen/6]=0;

  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      R[j][i]=AM[i+n][j];L[j][i]=AM[j][i+n];
      G[j][i]=R[j][i]-L[j][i];
    }
  }

  if(reduce_mat(G,n,m,1)){//reduction occurred
    //redefine R
    for(i=0;i<m;i++){
      for(j=0;j<n;j++)
	R[j][i]=G[j][i]+L[j][i];
    }
    //Make di6
    V=CRNPN3(L,R,n,m,layers,&layers2);
    amtodig(V,Vlen/6,out);
    free((char*)V);
  }
  else//no reduction. size should be right
    strcpy(out,di6);

  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  free_imatrix(L,0,n-1,0,m-1);
  free_imatrix(R,0,n-1,0,m-1);
  free_imatrix(G,0,n-1,0,m-1);
  return out;

}

//old version: 2 layers only
char *dynamic_isomorphMAold(char *di6, int n, int m){
  int entries,i,j;
  int **L=imatrix(0,n-1,0,m-1);//left stoich mat
  int **R=imatrix(0,n-1,0,m-1);//right stoich mat
  int **G=imatrix(0,n-1,0,m-1);//stoich mat
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix
  int r1=n+m;
  int r2=4*r1*r1;//minimum dimension of V (layers here)
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  char *out=(char *)malloc((size_t) ((Vlen/6+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+2*r1;//(layers here)
  out[2+Vlen/6]=0;

  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      R[j][i]=AM[i+n][j];L[j][i]=AM[j][i+n];
      G[j][i]=R[j][i]-L[j][i];
    }
  }

  if(reduce_mat(G,n,m,1)){//reduction occurred
    //redefine R
    for(i=0;i<m;i++){
      for(j=0;j<n;j++)
	R[j][i]=G[j][i]+L[j][i];
    }
    //Make di6
    CRNPN2(L,R,n,m,V,Vlen);//(layers here)
    amtodig(V,Vlen/6,out);

  }
  else//no reduction. size should be right
    strcpy(out,di6);

  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  free_imatrix(L,0,n-1,0,m-1);
  free_imatrix(R,0,n-1,0,m-1);
  free_imatrix(G,0,n-1,0,m-1);

  return out;

}

// Create a canonical "dynamically isomorphic" CRN (i.e. same GK ODEs)
// Stoich mat reduces as in dynamic_isomorph, Left stoich mat reduced to 
// sign pattern
char *dynamic_isomorphGK(char *di6, int n, int m){
  int entries,i,j;
  int **L=imatrix(0,n-1,0,m-1);//left stoich mat
  int **R=imatrix(0,n-1,0,m-1);//right stoich mat
  int **G=imatrix(0,n-1,0,m-1);//stoich mat
  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix

  int r1=n+m;
  int r2=4*r1*r1;//minimum dimension of V
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  char *out=(char *)malloc((size_t) ((Vlen/6+3)*sizeof(char)));
  int flag=0;


  out[0]=38;
  out[1]=63+2*r1;
  out[2+Vlen/6]=0;

  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      R[j][i]=AM[i+n][j];L[j][i]=AM[j][i+n];
      G[j][i]=R[j][i]-L[j][i];
    }
  }

  if(reduce_mat(G,n,m,1))//transpose
    flag=1;
  //  sgn_mat(L, n, m);

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if(L[j][i] && L[j][i]!=max(1,-G[j][i])){
	L[j][i]=max(1,-G[j][i]);//ensures R is nonnegative
	flag=1;
      }
      R[j][i]=G[j][i]+L[j][i];
    }
  }

  //Make di6
  if(flag){
    CRNPN2(L,R,n,m,V,Vlen);
    amtodig(V,Vlen/6,out);
  }
  else//no change
    strcpy(out,di6);

  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  free_imatrix(L,0,n-1,0,m-1);
  free_imatrix(R,0,n-1,0,m-1);
  free_imatrix(G,0,n-1,0,m-1);

  return out;

}

//order reactions by net production of species: highest (inflows) first
void orderreacsbynetprod(int **S, int n, int m, int *neword){
  int i;
  int netprod[m];
  //column sums
  for(i=0;i<m;i++){
    neword[i]=i;
    netprod[i]=colsum(S,n,m,i);
  }
  qsort2(netprod,neword,0,m-1);
  reverse(netprod,m);
  reverse(neword,m);
  return;
}

//order reaction by net species production
char *isomorph_bynetprod(char *di6, int n, int m){
  int entries,i,j;
  int neword[m];
  int **G=imatrix(0,n-1,0,m-1);
  int **G1=imatrix(0,n-1,0,m-1);
  int **L=imatrix(0,n-1,0,m-1);//left stoich mat
  int **R=imatrix(0,n-1,0,m-1);//right stoich mat
  int **L1=imatrix(0,n-1,0,m-1);//left stoich mat
  int **R1=imatrix(0,n-1,0,m-1);//right stoich mat

  int **AM=di6toCRNam(di6,n,m,&entries);//Petri Net adjacency matrix

  int r1=n+m;
  int r2=4*r1*r1;//minimum dimension of V
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  char *out=(char *)malloc((size_t) ((Vlen/6+3)*sizeof(char)));

  out[0]=38;
  out[1]=63+2*r1;
  out[2+Vlen/6]=0;

  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      R[j][i]=AM[i+n][j];L[j][i]=AM[j][i+n];
      G[j][i]=R[j][i]-L[j][i];
    }
  }

  //reorder
  orderreacsbynetprod(G,n,m,neword);
  for(j=0;j<m;j++){// reactions
    for(i=0;i<n;i++){// species
      G1[i][j]=G[i][neword[j]];
      L1[i][j]=L[i][neword[j]];
      R1[i][j]=R[i][neword[j]];
    }
  }


  //Make di6
  CRNPN2(L1,R1,n,m,V,Vlen);
  amtodig(V,Vlen/6,out);


  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  free_imatrix(L,0,n-1,0,m-1);
  free_imatrix(R,0,n-1,0,m-1);
  free_imatrix(G,0,n-1,0,m-1);
  free_imatrix(L1,0,n-1,0,m-1);
  free_imatrix(R1,0,n-1,0,m-1);
  free_imatrix(G1,0,n-1,0,m-1);

  return out;

}

//More complex reduced form
char *reducenn1(char *di6, int n, int m){
  int entries,i,j,rk;
  int **AM=di6toCRNam(di6, n, m, &entries);
  int **S, **Sl, **Sn;
  matrix J;
  int totbasis;
  int **basis;
  int **basisT;
  int **Sr=imatrix(0,n-1,0,m-1);//right stoich mat
  int r1=n+m;
  int r2=4*r1*r1;//minimum dimension of V
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];
  char *out=(char *)malloc((size_t) ((Vlen/6+3)*sizeof(char)));
  int flag=0;

  AMtoSSl(AM, n, m, 0, &S, &Sl);

  basisT=poskerbasis(S,n,m,&totbasis,&rk);
  if(totbasis!=1){
    fprintf(stderr, "ERROR in reducenn1: expecting the stoichiometric matrix to have a 1D kernel. EXITING.\n");
    exit(0);
  }

  basis=transposemat(basisT,totbasis,m);

  Sn=scalecols(S, n, m, basisT[0]);//we need to use basisT here (row vec)
  reduce_mat(Sn,n,m,0);//reduce by rows
  printmat(S,n,m);
  printmat(Sn,n,m);


  for(i=0;i<m;i++){// reactions
    for(j=0;j<n;j++){// species
      Sr[j][i]=Sn[j][i]+Sl[j][i];
      if(Sr[j][i]<0){
	fprintf(stderr, "WARNING in reducenn1: cannot have a negative stoichiometry %d, %d, %d. EXITING.\n", Sn[j][i], Sr[j][i], Sl[j][i]);
	strcpy(out,di6);
	flag=1;//no change
      }
    }
  }

  if(!flag){
    out[0]=38;
    out[1]=63+2*r1;
    out[2+Vlen/6]=0;
    CRNPN2(Sl,Sr,n,m,V,Vlen);
    amtodig(V,Vlen/6,out);
  }

  free_imat(basisT,totbasis);
  free_imatrix(basis, 0, m-1, 0, totbasis-1);
  free_imatrix(AM, 0, n+m-1, 0, n+m-1);
  free_imatrix(Sn,0,n-1,0,m-1);
  free_imatrix(Sr,0,n-1,0,m-1);
  return out;
}

// output only one member in each equivalence class of dynamically 
// isomorphic CRNs (MA). First rewrite each CRN in canonical form
// Then short the list. The third argument is 1 for MA and 0 for GK
int dynamicshortfile(const char *infile, const char *outfile, const char *flag, int n, int m, bool q){
  int linelen,ret=0;
  char *line;
  unsigned long numl=0;
  int maxl=0;
  FILE *fd, *fd1;
  char *str, *str1;
  char systemcom[1024];
  char path[100];
  FILE *fp;
  int i,r1=n+m;
  char *ptn;
  int minlayers=2;
  int totcmplx=0;

  if((numl=numflines(infile, &maxl))<=0){//no reactions
    return 0;
  }
  line = (char*) malloc(sizeof(char) * (maxl));

  if(!(fd=fopen(infile, "r"))){
    fprintf(stderr, "ERROR in dynamicshortfile - File: %s could not be opened...\n", infile);
    exit(0);
  }
  if(!(fd1=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in dynamicshortfile - File: %s could not be opened...\n", outfile);
    exit(0);
  }

  i=1;
  while((linelen = gtline(fd, line, maxl)) > 0){
    str=getnthwd(line, 1);
    if(!(strcmp(flag, "netprod")))
       str1=isomorph_bynetprod(str, n, m);
    else if(!(strcmp(flag, "sourcecomplexes")))//source complex CRN
      str1=cmplxCRN(str, n, m, minlayers, 1, &totcmplx);
    else if(!(strcmp(flag, "sourcecomplexesfull")))//source complex CRN with multiplicity
      str1=cmplxCRN(str, n, m, minlayers, 2, &totcmplx);
    else if(!(strcmp(flag, "allcomplexes")))//all complex CRN
      str1=cmplxCRN(str, n, m, minlayers, 0, &totcmplx);
    else if(!(strcmp(flag, "NN1")))
      str1=reducenn1(str,n,m);
    else if(!(strcmp(flag, "MA")))
      str1=dynamic_isomorphMA(str, n, m);
    else if(!(strcmp(flag, "GK")))
      str1=dynamic_isomorphGK(str, n, m);
    else if(!(strcmp(flag, "nofilter")))//Just remove isomorphs
      str1=strdup(str);
    else{
      fprintf(stderr, "ERROR in dynamicshortfile - filter \"%s\" not recognised.\n", flag);
      exit(0);
    }


    fprintf(fd1, "%s\n", str1);
    //fprintf(stderr, "%s, %s, %d\n", str, str1, i++);fflush(stderr);
    free(str);free(str1);
    ret++;
  }
  fclose(fd);
  fclose(fd1);
  free(line);

  if(!(strcmp(flag, "netprod"))) // no isomorphism checking
    return ret;

  if(!(strcmp(flag, "sourcecomplexes")) || !(strcmp(flag, "allcomplexes")))//complex CRN
    r1=n+totcmplx;
  ptn=(char *)malloc((size_t) ((2*r1+1)*sizeof(char)));
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<r1;i++)
    ptn[i]='1';
  for(i=r1;i<r1+n;i++)
    ptn[i]='2';
  for(i=r1+n;i<2*r1;i++)
    ptn[i]='3';
  ptn[i]=0;

  //The flag "k" means keep original labelling, "v" writes to stderr the correspondence list
  if(!q)
    sprintf(systemcom, "nauty-shortg -k -v -q -f%s %s", ptn, outfile);
  else
    sprintf(systemcom, "nauty-shortg -k -q -f%s %s", ptn, outfile);
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  pclose(fp);
  sprintf(systemcom, "wc -l %s", outfile);
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp);
  ret=atoi(path);
  pclose(fp);

  free(ptn);
  return ret;
}


//create all the bimolecular CRNs with either one more species (r=0) or one more reaction (r=1) than the given one (while remaining bimolecular) and containing the given one as an induced subnetwork. 
//The switch allowopen (only relevant if r=1) means treat all reactions; otherwise exclude 0-->X_i and X_i --> 0. 
// Very important note: all (n,m+1)-CRNs can be generated from (n,m) ones using buildsupCRN; however all (n+1,m)-CRNs cannot be generated from all (n,m) one as some CRNs become invalid after removal of *any* species. E.g., A0+A1-->A0, A1-->A0+A1
//The rankchange flag is negative if we want all CRNs; 0 if we want to preserve rank; and 1 if we want it to increase by 1
unsigned long buildsupCRN(char *infile, char *outfile, int n, int m, unsigned long *totin, bool allowopen, bool r, int rankchange, int debug){
  int i,j,k,k1,nH,flag=0,flagl,flagr;
  int n1,m1;
  int **AM;
  int r1=n+m+1;//either a species or reaction to be added
  int r2=4*r1*r1;//minimum dimension of V
  int rem = r2%6;
  int Vlen=(rem==0)?r2:r2+6-rem;
  bool V[Vlen];//to store output adjmat
  int clen=Vlen/6;
  char out[clen+3];
  int N=(n+1)*(n+2)/2;
  int **H=imatrix(0, N-1, 0, n-1);//complexes
  int Ltot, Rtot, addrowl[m], addrowr[m];
  int xcl[m], xcr[m];
  FILE *fd, *fd1;
  char *ptn;
  char systemcom[1000];
  //left and right matrices
  int **L=imatrix(0,n,0,m);
  int **R=imatrix(0,n,0,m);
  char *str;
  int linelen=clen+4;//+1for the newline
  char line[linelen];
  int entries;
  unsigned long totout=0;
  unsigned long tot=0;
  int rk0;

  (*totin)=0;

  if(r){
    if(debug){fprintf(stderr, "Adding reactions.\n");}
    nH=genallbicomplexes(n, N, H);
  }
  else if(debug)
    fprintf(stderr, "Adding species.\n");

  if(!(fd=fopen("tempfiles/_tmp", "w"))){
    fprintf(stderr, "ERROR in buildsupCRN: could not open file \"tempfiles/_tmp\" for writing.\n");exit(0);
  }

  if(!(fd1=fopen(infile, "r"))){
    fprintf(stderr, "ERROR in buildsupCRN: could not open file \"%s\" for reading.\n", infile);exit(0);
  }

  out[0]=38;
  out[1]=63+2*r1;
  out[2+clen]=0;

  while(gtline(fd1, line, linelen) > 0){
    if(linelen>0 && line[0]!=38)//not di6
	continue;
    //generate L, R matrices
    AM=di6toCRNam(line, n, m, &entries);
    for(j=0;j<n;j++){
      for(k=0;k<m;k++){
	L[j][k]=AM[j][k+n];
	R[j][k]=AM[k+n][j];
      }
    }

    if(rankchange>=0)
      rk0=checkrank(AM,n,m);

    free_imatrix(AM, 0, r1-2, 0, r1-2);

    (*totin)++;

    if(debug){
      str=di6toreacstr(line, n,m,0);
      fprintf(stderr, "******%ld******\n%s", (*totin), str);
      free(str);
    }

    if(r){//add a reaction
      n1=n;m1=m+1;
      for(i=0;i<nH;i++){//left halves
	for(j=0;j<nH;j++){//right halves
	  if(j!=i && (allowopen || (i!=0 && j!=0) || (i==0 && j>n) || (j==0 && i>n))){//valid reactions
	    //check for repetition
	    flag=0;
	    for(k1=0;k1<m;k1++){//each existing react
	      flag=1;
	      for(k=0;k<n;k++){//each species
		if((H[i][k]!=L[k][k1]) || (H[j][k]!=R[k][k1])){//(i,j) not same as k1th react
		  flag=0;
		  break;
		}
	      }
	      if(flag==1)//repeated reaction
		break;
	    }
	    if(!flag){//not repeated, add to L and R
	      for(k=0;k<n;k++){
		L[k][m]=H[i][k];
		R[k][m]=H[j][k];
	      }
	      if(rankchange<0 || (rankchange==0 && checkrankLR(L,R,n,m+1)==rk0) || (rankchange==1 && checkrankLR(L,R,n,m+1)==rk0+1)){
		CRNPN2(L,R,n,m+1,V,Vlen);
		amtodig(V,clen,out);
		fprintf(fd, "%s\n", out);
		tot++;
		if((tot+1)%1000000==0){
		  printf("*");fflush(stdout);
		}
	      }
	    }
	  }
	}
      }
    }//end of add a reaction
    else{//add a species
      n1=n+1;m1=m;
      //Get column totals
      for(k=0;k<m;k++){
	Ltot=colsum(L, n, m, k);
	Rtot=colsum(R, n, m, k);
	addrowl[k]=3-Ltot;
	addrowr[k]=3-Rtot;
      }
      inittozero(xcl, m);
      flagl=1;
      while(flagl){//each left vector
	for(k=0;k<m;k++)
	  L[n][k]=xcl[k];
	flagr=1;inittozero(xcr, m);
	while(flagr){//each right vector
	  for(k=0;k<m;k++)
	    R[n][k]=xcr[k];
	  if(rankchange<0 || (rankchange==0 && checkrankLR(L,R,n+1,m)==rk0) || (rankchange==1 && checkrankLR(L,R,n+1,m)==rk0+1)){
	    CRNPN2(L,R,n+1,m,V,Vlen);
	    amtodig(V,clen,out);
	    fprintf(fd, "%s\n", out);
	    tot++;
	    if((tot+1)%1000000==0){
	      printf("*");fflush(stdout);
	    }
	  }
	  flagr=nextnumb(xcr, m, addrowr);
	}
	flagl=nextnumb(xcl, m, addrowl);
      }

    }//end of add species

  }
  fclose(fd);
  fclose(fd1);

  if(r)
    free_imatrix(H,0,N-1,0,n-1);
  free_imatrix(L,0,n,0,m);
  free_imatrix(R,0,n,0,m);


  ptn=(char *)malloc((size_t) ((2*(n1+m1)+1)*sizeof(char)));
  for(i=0;i<n1;i++)
    ptn[i]='0';
  for(i=n1;i<n1+m1;i++)
    ptn[i]='1';
  for(i=n1+m1;i<2*n1+m1;i++)
    ptn[i]='2';
  for(i=2*n1+m1;i<2*(n1+m1);i++)
    ptn[i]='3';
  ptn[i]=0;

  printf("\n");
  if(tot>5000000){
    printf("removing isomorphic reactions (could take a while).\n");fflush(stdout);
  }
  //sprintf(systemcom, "../../nauty/nauty26r5/shortg -q -f%s tempfiles/_tmp %s", ptn,outfile); /* && rm tempfiles/_tmp */
  sprintf(systemcom, "nauty-shortg -q -f%s tempfiles/_tmp %s", ptn,outfile); /* && rm tempfiles/_tmp */
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  system(systemcom);
  free(ptn);



  if(!(fd1=fopen(outfile, "r"))){
    fprintf(stderr, "ERROR in buildsupCRN: could not open file \"%s\" for reading.\n", outfile);exit(0);
  }
  while(gtline(fd1, line, linelen) > 0 && !iscomline(line)){
    totout++;
    if(debug){
      str=di6toreacstr(line, r?n:n+1,r?m+1:m,0);
      fprintf(stderr, "******%ld******\n%s", totout, str);
      free(str);
    }
  }
  fclose(fd1);

  return totout;
}

//check if a set of CRNs in one file is a subset of those in the other.
//assume systems in di6 format and all of the same dimension.
//do not assume canonical labelling
bool CRNsubset(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5, int n, int m, int debug){
  int i, ret=0;
  char systemcom[1024];
  char ptn[2*(n+m)+1];
  FILE *fp, *fp1;
  unsigned long sz1, sz2, sz3, sz4;
  char path[1035];
  int maxl=0;
  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<n+m;i++)
    ptn[i]='1';
  for(i=n+m;i<2*n+m;i++)
    ptn[i]='2';
  for(i=2*n+m;i<2*(n+m);i++)
    ptn[i]='3';
  ptn[i]=0;

  sz1=numflines(file1, &maxl);
  sz2=numflines(file2, &maxl);
  if(sz1>1000000 || sz2>1000000)
    printf("Large list: this could take a while...\n");

  //sprintf(systemcom, "../../nauty/nauty26r5/shortg -q -f%s tempfiles/_tmp %s", ptn,outfile); /* && rm tempfiles/_tmp */
  sprintf(systemcom, "nauty-labelg -q -f%s %s tempfiles/_tmp1", ptn,file1);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }

  sprintf(systemcom, "nauty-labelg -q -f%s %s tempfiles/_tmp2", ptn,file2);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp1=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  pclose(fp);pclose(fp1);

  sprintf(systemcom, "sort -u tempfiles/_tmp1 > tempfiles/_tmp1a |wc -l");
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  pclose(fp);

  sprintf(systemcom, "sort -u tempfiles/_tmp2 > tempfiles/_tmp2a |wc -l");
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  pclose(fp);


  sprintf(systemcom, "comm -23 tempfiles/_tmp1a tempfiles/_tmp2a |wc -l");
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp);
  sz3=atoi(path);
  pclose(fp);

  if(sz3 && file3){//write difference
    sprintf(systemcom, "comm -23 tempfiles/_tmp1a tempfiles/_tmp2a > %s", file3);
    if(debug){fprintf(stderr, "%s\n", systemcom);}
    if(!(fp=popen(systemcom, "r"))){
      perror("couldn't open pipe.\n");exit(0);
    }
    pclose(fp);
    fprintf(stderr, "networks in \"%s\" but not in \"%s\" written to \"%s\"\n", file1, file2, file3);
  }

  sprintf(systemcom, "comm -23 tempfiles/_tmp2a tempfiles/_tmp1a |wc -l");
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  fgets(path, sizeof(path), fp);
  sz4=atoi(path);
  pclose(fp);

  if(sz4 && file4){//write difference
    sprintf(systemcom, "comm -23 tempfiles/_tmp2a tempfiles/_tmp1a > %s", file4);
    if(debug){fprintf(stderr, "%s\n", systemcom);}
    if(!(fp=popen(systemcom, "r"))){
      perror("couldn't open pipe.\n");exit(0);
    }
    pclose(fp);
    fprintf(stderr, "networks in \"%s\" but not in \"%s\" written to \"%s\"\n", file2, file1, file4);
  }

  if(file5){
    sprintf(systemcom, "comm -12 tempfiles/_tmp1a tempfiles/_tmp2a > %s", file5);
    if(debug){fprintf(stderr, "%s\n", systemcom);}
    if(!(fp=popen(systemcom, "r"))){
      perror("couldn't open pipe.\n");exit(0);
    }
    pclose(fp);
    fprintf(stderr, "networks in the intersection of \"%s\" and \"%s\" written to \"%s\"\n", file1, file2, file5);
  }

  fprintf(stderr, "\t%s: %ld\n\t%s: %ld\n\tnetworks in \"%s\" but not in \"%s\": %ld\n\tnetworks in \"%s\" but not in \"%s\": %ld\n", file1, sz1, file2, sz2, file1, file2, sz3, file2, file1, sz4);

  if(sz3==0 && sz4==0){
    fprintf(stdout, "\"%s\" and \"%s\" contain the same reaction networks.\n\n", file1, file2);
    ret=1;
  }
  else{
    if(sz3==0){
      fprintf(stdout, "all networks in \"%s\" are also in \"%s\".\n\n", file1, file2);
      ret=1;
    }
    if(sz4==0){
      fprintf(stdout, "all networks in \"%s\" are also in \"%s\".\n\n", file2, file1);
      ret=1;
    }
  }

  if(!ret)  
    fprintf(stdout, "neither file is a subset of the other.\n");
  return ret;
}


//check if a set of CRNs in one file has LHS from another set of
//LHS-only CRNs, and write to fileout
//assume systems in di6 format and all of the same dimension.
//do not assume canonical labelling
int CRNwithsources(const char *filein, const char *filesources, const char *fileout, int n, int m, int debug){
  int i;
  char systemcom[1024];
  char ptn[2*(n+m)+1];
  FILE *fp, *fdin, *fdin2, *fdout;
  int sz1, maxl=0;
  char **sourcearray;
  int numsources;
  char *str;
  //int minlayers=2;
  int totcmplx=0;
  char oneline[1024];
  char twoline[1024];
  int tot=0;
  int numl;


  for(i=0;i<n;i++)
    ptn[i]='0';
  for(i=n;i<n+m;i++)
    ptn[i]='1';
  for(i=n+m;i<2*n+m;i++)
    ptn[i]='2';
  for(i=2*n+m;i<2*(n+m);i++)
    ptn[i]='3';
  ptn[i]=0;

  if((sz1=numflines(filein, &maxl))>1000000)
    printf("Large list: this could take a while...\n");

  //Canonically label given list of source CRNs
  sprintf(systemcom, "nauty-labelg -q -f%s %s tempfiles/_sourcelistlabelled", ptn,filesources);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  pclose(fp);


  //generate source arrays
  sourcearray=genarraydat("tempfiles/_sourcelistlabelled", &numsources);
  if(numsources==0){
    fprintf(stderr, "No sources found. EXITING.\n");exit(0);
  }
  //Assumes all the source CRNs have the same number of layers
  numl=numlayers(sourcearray[0], n, m);//get the number of layers

  //generate source CRNs for CRNs in filein
  if(!(fdin = fopen(filein, "r"))){
    fprintf(stderr, "ERROR in CRNwithsources: \"%s\" could not be opened for reading.\n", filein);
    exit(0);
  }
  if(!(fdout = fopen("tempfiles/_tmpsources", "w"))){
    fprintf(stderr, "ERROR in CRNwithsources: \"tempfiles/_tmpsources\" could not be opened for writing.\n");
    exit(0);
  }
  while(getline0(fdin, oneline, 1024) > 0){
    //    numl=numlayers(oneline, n, m);//maintain the layers (bad idea?)
    //use same number of layers as sources, if possible
    str=cmplxCRN(oneline, n, m, numl, 2, &totcmplx);
    fprintf(fdout, "%s\n", str);
    free((char*)str);
  }
  fclose(fdin);
  fclose(fdout);

  //Canonically label source CRNs
  sprintf(systemcom, "nauty-labelg -q -f%s tempfiles/_tmpsources tempfiles/_tmpsourceslabelled", ptn);
  if(debug){fprintf(stderr, "%s\n", systemcom);}
  if(!(fp=popen(systemcom, "r"))){
    perror("couldn't open pipe.\n");exit(0);
  }
  pclose(fp);

  //Write only those with found sources to output file
  if(!(fdin = fopen(filein, "r"))){
    fprintf(stderr, "ERROR in CRNwithsources: \"%s\" could not be opened for reading.\n", filein);
    exit(0);
  }
  if(!(fdin2 = fopen("tempfiles/_tmpsourceslabelled", "r"))){
    fprintf(stderr, "ERROR in CRNwithsources: \"tempfiles/_tmpsourceslabelled\" could not be opened for reading.\n");
    exit(0);
  }
  if(!(fdout = fopen(fileout, "w"))){
    fprintf(stderr, "ERROR in CRNwithsources: \"%s\" could not be opened for writing.\n", fileout);
    exit(0);
  }
  while(getline0(fdin, oneline, 1024) > 0){
    getline0(fdin2, twoline, 1024);
    if(isinarray(sourcearray, numsources, twoline)>=0){
      fprintf(fdout, "%s\n", oneline);
      tot++;
    }
    
  }
  fclose(fdin);
  fclose(fdin2);
  fclose(fdout);

  freearraydat(sourcearray, numsources);
  return tot;
}


//read an isomorphism table in the format output by NAUTY and store
//as an array of unsigned long integers
unsigned long **readisomorphtable(char *fname, unsigned long *tot){
  FILE *fd;
  int maxl=0;
  unsigned long numl=numflines(fname, &maxl);
  char oneline[maxl];
  char *block, *str;
  unsigned long **out = (unsigned long**) malloc(sizeof(unsigned long*)*numl);
  int linelen;
  if(!(fd = fopen(fname, "r"))){
    fprintf(stderr, "ERROR in readisomorphtable: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  (*tot)=0;
  //cout << "here\n";
  linelen=getline0(fd, oneline, maxl);
  while(linelen > 0){
    //fprintf(stderr, "here: %s\n", oneline);
    if(!(iscomline(oneline))){//could be a broken line
      block=strdup(oneline);
      while((linelen=getline0(fd, oneline, maxl))>0 && !strstr(oneline, ":")){
	block=(char*) realloc(block, sizeof(char*) *(strlen(block)+strlen(oneline)+2));
	strcat(block, " ");strcat(block, oneline);
      }
      str=block;//needed as getulvec shifts pointer
      out[(*tot)]=getulvec(str);
      if(out[(*tot)][0]<2 || out[(*tot)][1]!=(*tot)+1){
	fprintf(stderr, "Something went wrong in readisomorphtable \"%s\". EXITING.\n",block);printvec(out[(*tot)]+1,out[(*tot)][0]);exit(0);
      }
      (*tot)++;
      free(block);
    }
    else
      linelen=getline0(fd, oneline, maxl);
  }
  fclose(fd);
  return out;
}

//Remember: NAUTY produces output offset at 1, not 0.
// uses the unsigned long version of isinarray where indices begin at
//1, not 0.
unsigned long getallisomorphs(char *isotablefile, char *infilepre, char *infileDI, char *infile, char *outfile, unsigned long *numin){
  unsigned long **table;
  unsigned long p=0, q, k, tot=0;
  char **inarray, **prearray, **DIarray;
  unsigned long numpre, numDI;
  FILE *fd;
  unsigned long totout=0;

  //generate arrays
  table=readisomorphtable(isotablefile, &tot);
  inarray=genarraydat(infile, numin);
  prearray=genarraydat(infilepre, &numpre);
  DIarray=genarraydat(infileDI, &numDI);
  if(tot!=numDI){
    fprintf(stderr, "Total in table (%ld) should be equal to total in input file (%ld). Something went wrong in getallisomorphs. EXITING.\n", tot, numDI);exit(0);
  }

  if(!(fd=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in getallisomorphs: \"%s\" could not be opened for writing.\n", outfile);
    exit(0);
  }

  for(p=0;p<(*numin);p++){//Each input network
    if((q=isinarray(DIarray,numDI,inarray[p]))){//lots of searching, so can be slow: how to speed up? Assume ordered?
      fprintf(stderr, "q = %ld\n", q);
      for(k=0;k<table[q-1][0]-1;k++){
	if(table[q-1][k+2] > numpre){
	  fprintf(stderr, "Index too large in getallisomorphs: something went wrong. EXITING.\n");exit(0);
	}
	fprintf(fd, "%s\n", prearray[table[q-1][k+2]-1]);
	totout++;
	//fprintf(stderr, "totout=%ld\n",totout);
      }
    }
    else{//don't allow failure: must find input network
      fprintf(stderr, "ERROR in getallisomorphs: couldn't find the input network %s. EXITING.\n", inarray[p]);exit(0);

    }
  }

  fclose(fd);
  freearraydat(inarray, (*numin));
  freearraydat(prearray, numpre);
  freearraydat(DIarray, numDI);
  ulfree_ul(table,tot);
  return totout;
}

