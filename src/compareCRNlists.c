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
#include "isomorph.h"
#include "analysereacs.h"


int main(int argc, char *argv[]){
  int n, m;
  int debug=0;
  if(argc<5){
    fprintf(stderr, "You must supply at least four arguments to \"compareCRNlists\":\n\tfile 1 in di6 format;\n\tfile 2 in di6 format;\n\tnumber of species;\n\tnumber of reactions.\n\tOptionally, in addition you may specify up to three output files to hold the file differences and intersection. EXITING.\n");
    exit(0); 
  }
  if(guessreacformat(argv[1])!=1 || guessreacformat(argv[2])!=1){
    fprintf(stderr, "The files do not appear to be in di6 format. You can first convert them to di6 format using \"convertformat\". EXITING.\n");
    exit(0);
  }

  if(!ispureint(argv[3]) || !ispureint(argv[4])){
    fprintf(stderr, "Third, fourth and fifth arguments must be integers. EXITING.\n");exit(0);
  }
  else{
    n=atoi(argv[3]);m=atoi(argv[4]);
  }

  if(argc>=8)
    CRNsubset(argv[1], argv[2], argv[5], argv[6], argv[7], n, m, debug);
  else if(argc>=7)
    CRNsubset(argv[1], argv[2], argv[5], argv[6], NULL, n, m, debug);
  else if(argc>=6)
    CRNsubset(argv[1], argv[2], argv[5], NULL, NULL, n, m, debug);
  else
    CRNsubset(argv[1], argv[2], NULL, NULL, NULL, n, m, debug);
  return 0;
}

