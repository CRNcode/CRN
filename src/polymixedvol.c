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
#include "symbolic.h"
#include "pos.h"

int main(int argc, char *argv[]){
  int opt;
  int debug=0;
  //int effort=0;//default
  int mainargs=0;
  char *infile;
  const char *filter=NULL;
  char HelpStr[2000];
  char *str;


  sprintf(HelpStr, "\nExample command:\n\t%s -v polyfile 2> tempfiles/tmp.log\n\nMeaning:\n\tExamine the file \"polyfile\" assumed to contain a set of multivariate\n\tpolynomials, and attempt to compute their mixed volume.\n\t(By Bernstein's theorem, this gives a bound on the number\n\tof isolated solutions in (C\\{0})^n, where \"n\" is the number\n\tof variables.)\n\nNot well tested: treat with caution!\n\n", argv[0]);

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

  str=readfileintostr(infile);
  fprintf(stderr, "%s\n", str);
  fprintf(stderr, "Mixed volume of polynomials: %d\n", polymixedvol(infile, debug));
  free(str);
  return 0;

}
