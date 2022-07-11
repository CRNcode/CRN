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

  sprintf(HelpStr, "\nExample command:\n\t%s -v -fposorth -e2 polyfile 2> tempfiles/tmp.log\n\nMeaning:\n\tExamine the file \"polyfile\" assumed to contain a single multivariate\n\tpolynomial, and attempt to determine if this polynomial is positive\n\tor nonnegative on the positive orthant. Put in a medium level of\n\teffort: \"-e2\". Send debugging output to \"tempfiles/tmp.log\"\n\nOptions:\n\t\"-v\" means verbose output;\n\t\"-f[filter]\" means a filter - see below for a list.\n\t\"-e[effort]\" means effort: an integer from 0 to 5, the default is 0.\n\nFilters:\n\tThe most important filters are \"posorth\": attempt to determine\n\tpositivity on the positive orthant; and \"posall\": attempt to\n\tdetermine positivity on all of R^n, except possibly at 0. Other\n\tfilters are \"newton\": examine the Newton polytope of the\n\tpolynomial; \"heuristic\": do only heuristic tests to determine\n\tpositivity on the positve orthant, i.e., avoid semidefinite\n\tprogramming; and \"factor\": attempt to factor using GiNaC.\n\t\n\nNotes: The polynomial is homogenised, and may be simplified in many ways,\nbefore any tests. If it has many terms or is of high degree, then\nsome tests will likely fail, and the programme may need to be killed.\nThe effort parameter is crucial if semidefinite programming is used\nin any of the tests. If you really want to track what is being done,\nthen increase the debug level, by swapping the option \"-v\" for,\nsay, \"-d5\": but be warned much of this debugging output will be\nhard to understand!\n\n", argv[0]);

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
      fprintf(stderr, "%s", HelpStr);
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

  polytest(infile, filter, effort-1, debug);
  return 0;

}
