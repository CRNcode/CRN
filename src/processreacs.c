/* Copyright (C) 2013-2022, Murad Banaji
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

#include "analysereacs.h"

int main(int argc, char *argv[])
{
  int i,q=1; // q=1 is quick, q=0 is slow but more info
  char *fname=NULL;
  bool htswtch=0;
  bool statswtch=0;
  for(i=1;i<argc;i++){
    if(strcmp(argv[i], "--html")==0)
      htswtch=1;
    else if(strcmp(argv[i], "--stats")==0)
      statswtch=1;
    else fname=argv[i];
  }
  if(!fname){
    fprintf(stderr, "ERROR: No data file supplied.\n\n");
    exit(0);
  }
  analysereacs(fname, q, htswtch,statswtch);
  return 0;
}

