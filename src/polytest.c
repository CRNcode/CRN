/* Copyright (C) 2010-2022, Murad Banaji
 *
 * This file is part of QUALMAT
 *
 * QUALMAT is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * QUALMAT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUALMAT: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 


 */

#include "basics.h"
#include "symbolic.h"
#include "pos.h"

int main(int argc, char *argv[]){
  int opt;
  int debug=0;
  int effort=0;//default
  int mainargs=0;
  char *infile;
  const char *filter=NULL;
  char HelpStr[2000];

  sprintf(HelpStr, "\nExample command:\n\t%s -v -fnewton -e0 polyfile 2> tmp.log\n\nOptions:\n\t\"-v\" means full debugging output;\n\t\"-f[filter]\" means a filter (\"factor\", \"newton\", \"heuristic\", \"posorth\", \"posall\", or \"all\" for all tests (default)).\n\t\"-e[effort]\" means effort (integer from 0 to 5, default is 0).\n\nAfter the options, the remaining argument is the input file which holds the polynomial. The polynomial is homogenised (and may be simplified) before any tests.\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "vf:e:hd:")) != -1) {
    switch (opt) {
    case 'v': //verbose
      debug = 1;
      break;
    case 'd': //verbose
      debug = atoi(optarg);
      break;
    case 'f'://filter
      filter = optarg;
      break;
    case 'e'://effort
      if(!ispureint(optarg) || ((effort = atoi(optarg))<0 || effort >5)){
	fprintf(stderr, "The argument of -e (\"%s\") should be the effort: an integer from 0 to 5. EXITING.\n", optarg);exit(EXIT_FAILURE);
      }
      break;
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-v] [-f filter] [-e effort] polyfile\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  //fprintf(stderr, "optind=%d\n", optind);
  while(optind < argc){
    switch (mainargs) {
    case 0:
      infile=argv[optind];
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<1){
    fprintf(stderr, "You need to supply an input file.\n%s", HelpStr);
    exit(0); 
  }

  if(!filter)
    filter="all";

  //fprintf(stderr, "infile=%s, outfilehopf=%s, outfilenohopf=%s, species=%d, reactions=%d, filter=%s, effort=%d, debug=%d\n", infile, outfilehopf, outfilenohopf, numspec, numreac, filter, effort, debug);

  printf("The result of polytest is %d\n", polytest(infile, filter, effort-1, debug));
  return 0;

}
