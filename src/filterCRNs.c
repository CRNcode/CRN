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
#include "analysereacs.h"
int main(int argc, char *argv[]){
  int opt, numspec,numreac;
  int debug=0;
  unsigned long ret;
  int mainargs=0;
  char HelpStr[2000];
  char *infile, *outfile, *filter=NULL;
  int filterarg=0;

  sprintf(HelpStr, "\nExample command:\n\t%s -v infile.d6 outfile.d6 3 4 rank 3 2> logfiles/tmp.log\n\nMeaning:\n\tExamine the file \"infile.d6\" containing 3-species, 4-reaction CRNs,\n\tin digraph6 format, and write to \"outfile.d6\" only those CRNs\n\tof rank 3. Write debugging output to \"logfiles/tmp.log\"\n\nOptions:\n\t\"-v\" means verbose output.\n\nThe following arguments are mandatory:\n\tinput file in d6 format\n\toutput file in d6 format\n\tnumber of species\n\tnumber of reactions\n\tfilter\n\t[for some filters] an additional numerical argument.\n\nThe list of supported filters is long and can be found by examining\n\"allfilters_annotated\".\n\n", argv[0]);


  //options?
  while ((opt = getopt(argc, argv, "vhd:")) != -1) {
    switch (opt) {
    case 'v': //verbose
      debug = 1;
      break;
    case 'd': //very verbose (keep hidden - only for real debugging)
      debug = atoi(optarg);
      break;
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
      infile=argv[optind];
      break;
    case 1:
      if(!strcmp(argv[optind],"NULL"))
	outfile=NULL;//pseudotest
      else
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
      filter=argv[optind];
      break;
    case 5:
      if(!ispureint(argv[optind])){
	fprintf(stderr, "The sixth argument \"%s\" must be a number passed to the filter. EXITING.\n", argv[optind]);exit(EXIT_FAILURE);
      }
      filterarg=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }


  if(mainargs<5){
    fprintf(stderr, "You must supply at least five arguments.\n%s", HelpStr);
    exit(0);
  }


  if(!strcmp(filter,"rank")){
    if(mainargs<6){
      fprintf(stderr, "The sixth argument must be the rank. EXITING.\n");
      exit(0);
    }
  }
  if(outfile && !strcmp(filter,"deficiency")){
    if(mainargs<6){
      fprintf(stderr, "The sixth argument must be the deficiency. EXITING.\n");
      exit(0);
    }
  }

  if(strstr(filter,"concord1") || strstr(filter,"MAeqaccord") || strstr(filter,"MAeqconcord") || strstr(filter,"JadmitsIpair") || strstr(filter,"JMAadmitsIpair") || strstr(filter,"JsquaredisP0") || strstr(filter,"JMAconcord") || strstr(filter,"Jnonsing") || strstr(filter,"JMAnonsing") || strstr(filter,"Jcomp2nonsing") || strstr(filter,"Jcomp2detsigned") || strstr(filter,"JMAcomp2nonsing") || strstr(filter,"JMAcomp2detsigned") || strstr(filter,"JMAcomp2detnonstationary") || strstr(filter,"J2pIdetpos") || strstr(filter,"J2pIMAdetpos") || strstr(filter,"QMAnegdef") || strstr(filter, "QMAnegsemidef") || strstr(filter, "QMApossemidef") || strstr(filter, "JMAisP0") || strstr(filter, "Jcomp2isP0") || strstr(filter, "JMAsquaredisP0") || strstr(filter, "JMAcomp2isP0") || strstr(filter, "MABTbif") || strstr(filter, "MAeqdegen2") || strstr(filter, "JMArealspec")){
    if(mainargs<6){
      fprintf(stderr, "The sixth argument must be effort: 0 means don't use SDP; 1 means use SDP to order 0; etc. If in doubt, try \"1\". EXITING.\n");
      exit(0);
    }
  }

  ret=filterCRNs(infile,outfile,numspec,numreac,filter,filterarg,debug);

  if(outfile)//not a pseudotest
    printf("%ld satisfy the condition.\n", ret);
  return 0;
}

