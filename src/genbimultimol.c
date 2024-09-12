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

//At most bimolecular reactant complexes and at most tetramolecular
//product complexes

#include "isomorph.h"
#include "basics.h"

int main(int argc, char *argv[]){
  int opt;
  int debug=0;
  int mainargs=0;
  char *outfile;
  const char *filter=NULL;
  int open=0;
  int numspec, numreac, outmol;
  char HelpStr[2000];

  sprintf(HelpStr, "\nExample command:\n\t%s s3r4.d6 3 4 3\n\nMeaning:\n\tGenerate all (nonisomorphic) CRNs on 3 species and\n\t4 reactions with at most bimolecular reactant complexes,\n\tand at most trimolecular product complexes. Write these in\n\tdigraph6 format to the file \"s3r4.d6\"\n\nOptions:\n\t\"-v\" means verbose output;\n\t\"-f[filter]\" means a filter. E.g., \"DN\" means dynamically\n\tnontrivial, \"connect\" means reactions\n\twith connected Petri Net graphs,\n\t\"ALL\", which is the default. Some filters can be strung together\n\t - e.g. \"DNgenuine\" means dynamically nontrivial and genuine.\n\nArguments (after the options, three mandatory, one optional):\n\toutput file\n\tnumber of species\n\tnumber of (irreversible) reactions\n\tmaximum molecularity of product complexes\n\tFinally an optional argument \"open\" means forbid reactions of the\n\tform 0 --> X and X --> 0. You can construct all fully open networks\n\tby adding in all the flows to these networks.\n\nNote that you can always filter reactions later - this may be better than\nfiltering as they are generated. You can also later convert to other formats.\n\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "vf:hd:")) != -1) {
    switch (opt) {
    case 'v': //verbose
      debug = 1;
      break;
    case 'd': //very verbose
      debug = atoi(optarg);
      break;
    case 'f'://filter
      filter = optarg;
      break;
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-v] [-f filter] outfile numspecies numreacs\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  //fprintf(stderr, "optind=%d\n", optind);
  while(optind < argc){
    switch (mainargs) {
    case 0:
      outfile=argv[optind];
      break;
    case 1:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The second argument \"%s\" must be the number of species. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      numspec=atoi(argv[optind]);
      break;
    case 2:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The third argument \"%s\" must be the number of reactions. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      numreac=atoi(argv[optind]);
      break;
    case 3:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The third argument \"%s\" must be the maximum molecularity of the product complexes. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      outmol=atoi(argv[optind]);
      break;
    case 4:
      if(!strcmp(argv[optind], "open"))
	open=1;
      break;
    }
    mainargs++;
    optind++;
  }
  if(mainargs<4){
    fprintf(stderr, "You must supply at least four arguments.\n%s", HelpStr);
    exit(0); 
  }

  if(!filter)
    filter="ALL";

  genbimultimol(numspec,numreac,outfile,open,filter,outmol,debug);

  return 0;
}

