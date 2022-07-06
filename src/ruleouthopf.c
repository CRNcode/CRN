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

#include "hopf.h"
#include "basics.h"

// arguments: 
// input file with one bimolecular crn in di6 format per line
// output file for crns with possible Hopf bifurcation (or NULL)
// output file for crns forbidding Hopf bifurcation (or NULL)
// number of species
// number of reactions


// Filters: -f[filtername]
// most likely to want: MAonly, GKonly

// Effort: (0-3) but most likely to want 0 or 1. 0 (the default) means 
// don't use SDP which can really slow things down

// -v means full debugging output

// Returns 1 if we can be confident that Hopf bifurcation
// is forbidden, and returns 0 otherwise

// Example command: ./bin/ruleouthopf -v -fMAonly -e0 crndatDNG/s3r4DNGrank3.d6 crndatDNG/s3r4DNGrank3possiblehopf.d6 crndatDNG/s3r4DNGrank3nohopf.d6 3 4 2> crndatDNG/s3r4hopfcheck.log


int main(int argc, char *argv[]){
  int opt;
  int debug=0;
  int effort=0;//default
  int mainargs=0;
  char *infile, *outfilehopf, *outfilenohopf, *filter=NULL;
  int numspec, numreac;
  char HelpStr[2000];

  sprintf(HelpStr, "\nExample command:\n\t%s -v -fMAonly -e0 infile.d6 outfilehopf.d6 outfilenohopf.d6 3 4 2> tmp.log\n\nOptions:\n\t\"-v\" means full debugging output;\n\t\"-f[filter]\" means a filter (e.g., \"MAonly\" for only mass action, \"GKonly\" for only general kinetics, of \"ALL\" (default) for both).\n\t\"-e[effort]\" means effort (integer from 0 to 3, default is 0).\n\nAfter the options, the remaining arguments are mandatory:\n\tinput file\n\toutput file for networks where H-bif is not ruled out\n\toutput file for networks where H-bif is ruled out\n\tnumber of species\n\tnumber of reactions.\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "vf:e:hd:")) != -1) {
    switch (opt) {
    case 'v': //verbose
      debug = 1;
      break;
    case 'd': //very verbose
      debug = atoi(optarg);
      break;
    case 'f'://filter
      //No need to check: the algorithm fails for unknown filters
      filter = optarg;
      break;
    case 'e'://filter
      if(!ispureint(optarg) || ((effort = atoi(optarg))<0 || effort >3)){
	fprintf(stderr, "The argument of -e (\"%s\") should be the effort: an integer from 0 to 3. EXITING.\n", optarg);exit(EXIT_FAILURE);
      }
      break;
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-v] [-f filter] [-e effort] infile outfilehopf outfilenohopf numspecies numreacs\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  //fprintf(stderr, "optind=%d\n", optind);
  while(optind < argc){
    switch (mainargs) {
    case 0:
      infile=argv[optind];
      break;
    case 1:
      outfilehopf=argv[optind];
      break;
    case 2:
      outfilenohopf=argv[optind];
      break;
    case 3:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The fourth argument \"%s\" must be the number of species. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      numspec=atoi(argv[optind]);
      break;
    case 4:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The fifth argument \"%s\" must be the number of reactions. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      numreac=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<5){
    fprintf(stderr, "You must supply at least five arguments.\n%s", HelpStr);
    exit(0); 
  }

  //fprintf(stderr, "infile=%s, outfilehopf=%s, outfilenohopf=%s, species=%d, reactions=%d, filter=%s, effort=%d, debug=%d\n", infile, outfilehopf, outfilenohopf, numspec, numreac, filter, effort, debug);

  ruleouthopf(infile, outfilehopf, outfilenohopf, numspec, numreac, filter, effort, debug);
  return 0;
}

