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

//Convert an input file in di6 format (at most trimolecular) to various other formats, including human readable. 

#include "analysereacs.h"
#include "basics.c"

int main(int argc, char *argv[]){
  int opt, numspec=0,numreac=0;
  int mainargs=0;
  char HelpStr[2000];
  char *infile=NULL, *intype=NULL, *outfile=NULL, *outtype=NULL;
  int intyp=0,outtyp=0;//integer versions
  long numout;

  sprintf(HelpStr, "\nExample command:\n\t%s -idi6 infile.di6 outfile.reacs reacstr 3 4\n\nOptions:\n\t\"-i[inputformat]\" tells us the input file format (\"di6\" or \"reacstr\" or \"Sauro\" or \"simpstr\"). If not specified, the algorithm tries to determine the format from the file, but may fail.\n\t\"-o[outformat]\" tells us the output file format (\"di6\" or \"reacstr\" or \"Sauro\" or \"simpstr\" or \"dot\"). If not specified, the algorithm tries to determine the format from the file ending, but may fail.\n\nThe arguments are:\n\tinput file\n\toutput file\n\tnumber of species (needed only if input format is \"di6\" or \"simpstr\")\n\tnumber of reactions (needed only if input format is \"di6\" or \"simpstr\")\n\nOne of input or output files must be in digraph 6 format.\n\n", argv[0]);

  //options?
  while ((opt = getopt(argc, argv, "i:o:h")) != -1) {
    switch (opt) {
    case 'i': //input format
      intype = optarg;
      if((!strcmp(intype, "di6") || !strcmp(intype, "d6")))
	intyp=1;
      else if((!strcmp(intype, "reacstr") || !strcmp(intype, "reacs")))
	intyp=2;
      else if((!strcmp(intype, "Sauro") || !strcmp(intype, "sauro")))
	intyp=3;
      else if((!strcmp(intype, "simpstr") || !strcmp(intype, "simp")))
	intyp=4;
      break;
    case 'o': //output format
      outtype = optarg;
      if((!strcmp(outtype, "di6") || !strcmp(outtype, "d6")))
	outtyp=1;
      else if((!strcmp(outtype, "reacstr") || !strcmp(outtype, "reacs")))
	outtyp=2;
      else if((!strcmp(outtype, "Sauro") || !strcmp(outtype, "sauro")))
	outtyp=3;
      else if((!strcmp(outtype, "simpstr") || !strcmp(outtype, "simp")))
	outtyp=4;
      else if((!strcmp(outtype, "dot")))
	outtyp=5;
      break;
    case 'h': //help
      fprintf(stderr, "%s", HelpStr);
      exit(0);
    default: /* '?' */
      fprintf(stderr, "Usage: %s -i[inputformat] -o[outputformat] infile outfile numspecies numreacs\n", argv[0]);
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
    }
    mainargs++;
    optind++;
  }


  if(!intype){//not given
    intyp=guessreacformat(infile);
    if(intyp==1)
      intype=(char*)"di6";
    else if(intyp==2)
      intype=(char*)"reacstr";
    else if(intyp==3)
      intype=(char*)"Sauro";
    else if(intyp==4)
      intype=(char*)"simpstr";
    else if(intyp==5){
      fprintf(stderr, "Could not be sure of the input format. Could be Sauro format in which case add the option \"-iSauro\", or simpstr format, in which case add the option \"-isimpstr\". EXITING.\n");
      exit(EXIT_FAILURE);
    }
    else{
      fprintf(stderr, "Could not be sure of the input format. EXITING.\n%s", HelpStr);
      exit(EXIT_FAILURE);
    }
  }

  if(!outtype){//not given
    outtyp=guessoutputformat(outfile);
    if(outtyp==1)
      outtype=(char*)"di6";
    else if(outtyp==2)
      outtype=(char*)"reacstr";
    else if(outtyp==3)
      outtype=(char*)"Sauro";
    else if(outtyp==4)
      outtype=(char*)"simpstr";
    else if(outtyp==5)
      outtype=(char*)"dot";
    else{
      fprintf(stderr, "Could not be sure of the output format. EXITING.\n%s", HelpStr);
      exit(EXIT_FAILURE);
    }
  }

  if(intyp==1 && (outtyp<=0 || outtyp>=6)){
    fprintf(stderr, "Given digraph6 input, output type must be one of \"reacstr\", \"simpstr\", \"Sauro\" and \"dot\". EXITING.\n");exit(EXIT_FAILURE);
  }
  if(outtyp==1 && (intyp<=0 || intyp>=5)){
    fprintf(stderr, "For digraph6 output, input type must be one of \"reacstr\", \"simpstr\", or \"Sauro\". EXITING.\n");exit(EXIT_FAILURE);
  }
  if(outtyp>=2 && outtyp<=4 && intyp!=1){
    fprintf(stderr, "For %s output, input type can only be digraph 6. EXITING.\n", outtype);exit(EXIT_FAILURE);
  }
  if(intyp>=2 && intyp<=4 && outtyp!=1){
    fprintf(stderr, "Given %s input, output type can only be digraph 6. EXITING.\n", outtype);exit(EXIT_FAILURE);
  }

  if((intyp==1 || intyp==4) && mainargs<4){
    fprintf(stderr, "With input in digraph 6 or simpstr format, the final two arguments must be the number of species and the number of reactions.\n");
    exit(0);
  }

  numout=convertformat(infile,intyp,outfile,outtyp,&numspec,&numreac);
  printf("Converted %ld reactions.\n", numout);

  return 0;
}

