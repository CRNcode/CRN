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

//Produces a technical report on a CRN: summary of several computations

#include "analysereacs.h"
#include "basics.c"

int main(int argc, char *argv[]){
  int t, opt, numspec=0,numreac=0;
  int mainargs=0;
  char HelpStr[2000];
  char *infile=NULL, *intype=NULL;
  char *filter=NULL;
  long numout;
  int debug=0;
  int effort=0;

  sprintf(HelpStr, "\nExample command:\n\t%s -v -idi6 infile.di6 3 4\n\nOptions:\n\t\"-i[inputformat]\" tells us the input file format (\"di6\" or \"reacstr\" or \"Sauro\" or \"simpstr\"). If not specified, the algorithm tries to determine the format from the file, but may fail.\n\t\"-e[effort]\": an integer from 0 9default) to 5. Only affects some computations.\n\t\"-f[filter]\": to perform only some computations. E.g., \"basic\", \"endo\", \"concord\", \"hopf\".\n\nThe arguments are:\n\tinput file\n\tnumber of species (needed only if input format is \"di6\" or \"simpstr\")\n\tnumber of reactions (needed only if input format is \"di6\" or \"simpstr\").\n\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "i:hve:d:f:")) != -1) {
    switch (opt) {
    case 'i': //input format
      intype = optarg;
      break;
    case 'f': //filter
      filter = optarg;
      break;
    case 'e': //input format
      effort = atoi(optarg);
      break;
    case 'v': //verbose output
      debug=1;
      break;
    case 'd': //very verbose output
      debug=atoi(optarg);
      break;
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-v] -e[effort] -i[inputformat] infile numspecies numreacs\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }


  while(optind < argc){
    switch (mainargs) {
    case 0:
      infile=argv[optind];
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
    }
    mainargs++;
    optind++;
  }

  if(mainargs<1){
    fprintf(stderr, "You need to specify an input file.\n%s", HelpStr);
    exit(EXIT_FAILURE);
  }


  if(!intype){//not given
    t=guessreacformat(infile);
    if(t==1)
      intype=(char*)"di6";
    else if(t==2)
      intype=(char*)"reacstr";
    else if(t==3)
      intype=(char*)"Sauro";
    else if(t==4)
      intype=(char*)"simpstr";
    else if(t==5){
      fprintf(stderr, "Could not be sure of the input format. Could be Sauro format in which case add the option \"-iSauro\", or simpstr format, in which case add the option \"-isimpstr\". EXITING.\n");
      exit(EXIT_FAILURE);
    }
    else{
      fprintf(stderr, "Could not be sure of the input format. EXITING.\n%s", HelpStr);
      exit(EXIT_FAILURE);
    }
  }

  if((!strcmp(intype,"di6") || !strcmp(intype,"d6") || !strcmp(intype,"simpstr") || !strcmp(intype,"simp")) && mainargs<3){
    fprintf(stderr, "With input in digraph 6 or simpstr format, the final two arguments must be the number of species and the number of reactions.\n");
    exit(0);
  }

  if(effort<0 || effort >5){
    fprintf(stderr, "Effort must be a number between 0 and 5.\n");
    exit(0);
  }


  numout=reacreportfile(infile, intype, filter, numspec, numreac, effort-1, debug);

  printf("Analysed %ld %s.\n", numout,numout==1?"network":"networks");
  return 0;
}

