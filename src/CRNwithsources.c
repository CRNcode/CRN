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
#include "isomorph.h"
#include "analysereacs.h"


int main(int argc, char *argv[]){
  int opt, n, m;
  int debug=0;
  char *file1=NULL, *file2=NULL, *file3=NULL;
  int mainargs=0;
  int tot;

  char HelpStr[2000];
  sprintf(HelpStr, "\nExample command:\n\t%s file1.d6 file2.d6 file3.d6 3 4\n\nMeaning:\n\tFind CRNs in \"file1.d6\" with LHS's from \"file2.d6\" all assumed\n\tto have 3 species and 4 irreversible reactions. \n\tOutput the found CRNs to\"file3.d6\"\n\tYou could generate the sources file with, for example,\n\t\t\"./bin/filterCRNs <filein> <fileout> 3 4 isomorphnewton1\"\n\nArguments: all five are mandatory.\n\tinput file in digraph6 format\n\tLHS CRN file in digraph6 format\n\tutput file in digraph6 format\n\tnumber of species\n\tnumber of reactions\n\nNote: this program relies on calls to nauty.\n\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "h")) != -1) {
    switch (opt) {
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "%s", HelpStr);
      exit(EXIT_FAILURE);
    }
  }

  while(optind < argc){
    switch (mainargs) {
    case 0:
      file1=argv[optind];
      break;
    case 1:
      file2=argv[optind];
      break;
    case 2:
      file3=argv[optind];
      break;
    case 3:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The fourth argument \"%s\" must be the number of species.\n%s", argv[optind], HelpStr);exit(EXIT_FAILURE);
      }
      n=atoi(argv[optind]);
      break;
    case 4:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The fifth argument \"%s\" must be the number of reactions.\n%s", argv[optind], HelpStr);exit(EXIT_FAILURE);
      }
      m=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<5){
    fprintf(stderr, "You need to specify five arguments.\n%s", HelpStr);
    exit(EXIT_FAILURE);
  }

  if(guessreacformat(argv[1])!=1 || guessreacformat(argv[2])!=1){
    fprintf(stderr, "The files do not appear to be in di6 format. You can first convert them to di6 format using \"convertformat\".\n%s", HelpStr);
    exit(0);
  }


  tot=CRNwithsources(file1, file2, file3, n, m, debug);
  fprintf(stderr, "Wrote %d CRNs to file \"%s\"\n", tot, file3);
  return 0;
}

