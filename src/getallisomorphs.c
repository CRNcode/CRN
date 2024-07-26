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


int main(int argc, char *argv[]){
  int opt, mainargs=0;
  char HelpStr[2000];
  char *isotablefile, *base_file, *reduced_file, *infile, *outfile;
  unsigned long totin,totout;

  sprintf(HelpStr, "\nExample command:\n\t%s isotable base.d6 reduced.d6 infile.d6 outfile.d6\n\nMeaning:\n\tThe idea is that we previously examined the CRNs in \"base.d6\",\n\tkept only one representative of each class in \"reduced.d6\", and\n\tstored a list of equivalence classes - in the format output by\n\tnauty - in \"isotable\". We now have a list of CRNs in \"infile.d6\",\n\tand wish to find all the CRNs from the original list in \"base.d6\"\n\twhich are isomorphic to CRNs in \"infile.d6\"; but just using the\n\tlookup table, rather than doing isomorphism checking. The\n\toutput is stored in \"outfile.d6\".\n\nFour input files and one output file are needed, in this order:\n\t1. The isomorphism table in nauty format\n\t2. The file used to generate the isomorphism table in digraph6 format\n\t3. The output from isomorphism generation in digraph6 format\n\t4. The input file whose isomorphs we want to find in digraph6 format\n\t5. The output file to hold isomorphs in digraph6 format.\n\nThe lookup could be generated, for example, using \"filterCRNs\", e.g. via:\n\n./bin/filterCRNs base.d6 reduced.d6 <numspec> <numreacs> dynisomorphMA 2> isotable\n\nAll files must be in digraph6 format.\n\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "i:hve:d:f:")) != -1) {
    switch (opt) {
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "%s", HelpStr);
      exit(EXIT_FAILURE);
    }
  }


  while(mainargs < argc+1){
    switch (mainargs) {
    case 1:
      isotablefile=argv[mainargs];
      break;
    case 2:
      base_file=argv[mainargs];
      break;
    case 3:
      reduced_file=argv[mainargs];
      break;
    case 4:
      infile=argv[mainargs];
      break;
    case 5:
      outfile=argv[mainargs];
      break;
    }
    mainargs++;
  }


  if(mainargs<6){
    fprintf(stderr, "You must supply at least five arguments.\n%s", HelpStr);
    exit(0);
  }

  totout=getallisomorphs(isotablefile, base_file, reduced_file, infile, outfile, &totin);
  printf("A total of %ld reactions went in and %ld came out.\n", totin, totout);
  return 0;
}

