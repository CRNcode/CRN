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
  int opt, n, m;
  int debug=0;
  char *file1=NULL, *file2=NULL;
  char *a=NULL, *b=NULL, *c=NULL;
  int mainargs=0;

  char HelpStr[2000];
  sprintf(HelpStr, "\nExample command:\n\t%s file1.d6 file2.d6 3 4 a.d6 b.d6 c.d6\n\nMeaning:\n\tCompare CRNs in \"file1.d6\" and \"file2.d6\" all assumed to have\n\t3 species and 4 irreversible reactions. Identify common reactions,\n\tup to isomorphism. Print one representative of each CRN in\n\t\"file1.d6\" but not \"file2.d6\" to \"a.d6\"; one representative\n\tof each CRN in \"file2.d6\" but not \"file1.d6\" to \"b.d6\";\n\tand one representative of each CRN common to both \"file1.d6\"\n\tand \"file2.d6\" to \"c.d6\". If some of these sets are empty or\n\tfilenames are not given, then output files will not be created\n\tbut a report will still be given.\n\nArguments: the first four are mandatory.\n\tinput file in digraph6 format\n\toutput file in digraph6 format\n\tnumber of species\n\tnumber of reactions\n\toptional: a file to hold CRNs in \"file1.d6\" but not \"file2.d6\"\n\toptional: a file to hold CRNs in \"file2.d6\" but not \"file1.d6\"\n\toptional: a file to hold CRNs in both \"file1.d6\" and \"file2.d6\".\n\nNote: this program relies on calls to nauty and makes heavy use of Linux\ntools such as \"sort\" and \"comm\".\n\n", argv[0]);

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
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The third argument \"%s\" must be the number of species.\n%s", argv[optind], HelpStr);exit(EXIT_FAILURE);
      }
      n=atoi(argv[optind]);
      break;
    case 3:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The fourth argument \"%s\" must be the number of reactions.\n%s", argv[optind], HelpStr);exit(EXIT_FAILURE);
      }
      m=atoi(argv[optind]);
      break;
    case 4:
      a=argv[optind];
      break;
    case 5:
      b=argv[optind];
      break;
    case 6:
      c=argv[optind];
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<4){
    fprintf(stderr, "You need to specify at least four arguments.\n%s", HelpStr);
    exit(EXIT_FAILURE);
  }

  if(guessreacformat(argv[1])!=1 || guessreacformat(argv[2])!=1){
    fprintf(stderr, "The files do not appear to be in di6 format. You can first convert them to di6 format using \"convertformat\".\n%s", HelpStr);
    exit(0);
  }


  CRNsubset(file1, file2, a, b, c, n, m, debug);
  return 0;
}

