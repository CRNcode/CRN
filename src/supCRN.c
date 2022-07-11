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
  int opt, numspec,numreac;
  int debug=0;
  int mainargs=0;
  char HelpStr[2000];
  char *infile, *outfile;
  bool allowopen=1, r=1;
  char *addspecreac=NULL;
  unsigned long totin,totout;
  int rankchange=-1;

  sprintf(HelpStr, "\nExample command:\n\t%s infile.d6 outfile.d6 3 4 spec\n\nMeaning:\n\tTake all the reactions in infile.d6: these have 3 species and\n\t4 irreversible reactions, and are assumed to be bimolecular.\n\tEnlarge in all possible ways by adding a new species, while\n\tpreserving bimolecularity. Return the results - a list of\n\tnonisomorphic CRNs - in \"outfile.d6\". \n\nOptions:\n\t\"-v\" means full debugging output\n\t\"-o\" means forbid reactions of the form 0-->X and X-->0.\n\t\"-r\" followed by 0 means preserve rank; 1 means increase the rank by 1.\n\tThe default is no checking of rank of enlarged network.\n\nAfter the options, the remaining arguments are mandatory:\n\tinput file\n\toutput file\n\tnumber of species\n\tnumber of reactions\n\t\"spec\" to add species or \"reac\" to add reactions.\n\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "vd:hor:")) != -1) {
    switch (opt) {
    case 'v': //verbose
      debug = 1;
      break;
    case 'd': //very verbose
      debug = atoi(optarg);
      break;
    case 'o': //fully open reacs
      allowopen=0;
      break;
    case 'r': //preserve rank?
      rankchange=atoi(optarg);
      if(rankchange!=0 && rankchange!=1){
	fprintf(stderr, "\"-r\" (currently followed by %s) must be followed by \"0\" (to preserve rank) or \"1\" (to increase the rank by 1. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      break;
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-v] [-o] [-r<rankchange>] infile outfile numspecies numreacs [spec/reac]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  while(optind < argc){
    switch (mainargs) {
    case 0:
      infile=argv[optind];
      break;
    case 1:
      outfile=argv[optind];
      break;
    case 2:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The third argument \"%s\" must be the number of species. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      numspec=atoi(argv[optind]);
      break;
    case 3:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The fourth argument \"%s\" must be the number of reactions. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      numreac=atoi(argv[optind]);
      break;
    case 4:
      addspecreac=argv[optind];
      if(!strcmp(addspecreac, "species") || !strcmp(addspecreac, "spec"))
	r=0;
      else if(!strcmp(addspecreac, "reactions") || !strcmp(addspecreac, "reac"))
	r=1;
      else{
	fprintf(stderr, "The final argument must be \"spec\" (to add species) or \"reac\" (to add reactions).\n");
	exit(EXIT_FAILURE);
      }
      break;
    }
    mainargs++;
    optind++;
  }


  if(mainargs<5){
    fprintf(stderr, "You must supply at least five arguments.\n%s", HelpStr);
    exit(0);
  }


  totout=buildsupCRN(infile, outfile, numspec, numreac, &totin, allowopen, r, rankchange, debug);
  printf("A total of %ld reactions went in %ld and came out.\n", totin, totout);
  return 0;
}

