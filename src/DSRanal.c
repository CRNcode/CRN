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
#include "analysereacs.h"

void checkTarjanCircuits1();
void checkDSRCondStar(const char *fname, int debug);

int main(int argc, char *argv[]){

  int i, n, m;
  int debug=0;
  char **array;
  char *str;
  int tot;
  int report;
  int rev=2;

  //checkDSRCondStar("datfiles/matpairDSR", debug);exit(0);
  //checkTarjanCircuits1();exit(0);
  if(argc<4){
    fprintf(stderr, "You need the name of a file containing a reaction systems in digraph6 format, followed by the number of species and the number of reactions. You can also add \"debug\" at the end to see verbose output. EXITING.\n");
    exit(0);
  }
  n=atoi(argv[2]);
  m=atoi(argv[3]);
  if(argc>=5 && (!strcmp(argv[4],"debug") || !strcmp(argv[4],"verbose"))){
    debug=1;
  }

  array=genarraydat(argv[1], &tot);

  for(i=0;i<tot;i++){
    fprintf(stderr,"System %d:\n",i+1);
    str=di6toreacstr(array[i], n, m, 0);
    fprintf(stderr, "%s\n", str);
    if(debug){fprintf(stderr, "%s\n", str);}
    if(DSRCondStar(array[i],n,m,&report,rev,debug)){
      fprintf(stderr, "Condition * holds%s\n", report==1?": no e-cycles found":"");
    }
    else{
      fprintf(stderr, "Condition * fails%s\n", report==-1?": es-cycles found":": e-cycles with odd intersection");
    } 
    free(str);
  }

  freearraydat(array, tot);
  
  return 0;
}
