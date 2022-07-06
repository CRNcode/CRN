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


int main(int argc, char *argv[]){
  int mainargs=0;
  char HelpStr[2000];
  char *isotablefile, *base_file, *reduced_file, *infile, *outfile;
  unsigned long totin,totout;

  sprintf(HelpStr, "\nExample command:\n\t%s isotablefile base_file.d6 reduced_file.d6 infile.d6 outfile.d6\n\nThe files are:\n\t1) The isomorphism table\n\t2)The file used to generate the isomorphism table\n\t3)The output from isomorphism generation\n\t4)The input file whose isomorphs we want to find\n\t5)The output file to hold isomorphs.\n\nThe lookup is generated using \"filterCRNs\", e.g. via:\n\n./bin/filterCRNs base_file.d6 reduced_file.d6 <numspec> <numreacs> dynisomorphMA 2> isotablefile\n\n", argv[0]);


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

