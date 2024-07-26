#####################################################
#
#  Makefile for the whole system:
#  Note that certain features are for readability
#  and not for compactness of code
#  
######################################################
CC = gcc -Wall
CXX = g++ -Wall
CFLAGS = -g -O2
CXXFLAGS= -g -O2
VERSION = 1.0
MODELNAME = CRNLIB$(VERSION)

HEADERS =	src/crn.h src/pos.h src/basics.h src/analysereacs.h src/symbolic.h src/hopf.h src/convex.h src/isomorph.h src/semidef.h

MAINLIB = 	src/lib/libcrn.a

MAINLINK =	-lcrn -lcln -lginac -lglpk

OBJDIR = 	src/obj

AR      =       ar crs

OBJECTS =	src/obj/pos.o src/obj/basics.o src/obj/convex.o src/obj/symbolic.o src/obj/semidef.o src/obj/analysereacs.o src/obj/hopf.o src/obj/isomorph.o

project: bin/polytest bin/polymixedvol bin/ruleouthopf bin/processreacs bin/reacreport bin/convertformat bin/filterCRNs bin/genbimol bin/genbitrimol bin/genzeroonetritrimol bin/genbimultimol bin/DSRanal bin/supCRN bin/compareCRNlists bin/CRNwithsources bin/getallisomorphs

clean: 
	rm bin/* src/lib/libcrn.a src/obj/*.o

src/obj/%.o: src/%.c $(HEADERS)
	$(CXX) -c -o $@ $<

src/lib/libcrn.a: $(HEADERS) $(OBJECTS)
	$(AR) src/lib/libcrn.a $(OBJECTS)



##################


bin/polytest: src/polytest.c $(MAINLIB)
	$(CXX) -o bin/polytest src/polytest.c -lm -Lsrc/lib $(MAINLINK)

bin/polymixedvol: src/polymixedvol.c $(MAINLIB)
	$(CXX) -o bin/polymixedvol src/polymixedvol.c -lm -Lsrc/lib $(MAINLINK)

bin/processreacs: src/processreacs.c $(MAINLIB)
	$(CXX) -o bin/processreacs src/processreacs.c -lm -Lsrc/lib $(MAINLINK)

bin/reacreport: src/reacreport.c $(MAINLIB)
	$(CXX) -o bin/reacreport src/reacreport.c -lm -Lsrc/lib $(MAINLINK)

bin/ruleouthopf: src/ruleouthopf.c $(MAINLIB)
	$(CXX) -o bin/ruleouthopf src/ruleouthopf.c -lm -Lsrc/lib $(MAINLINK)

bin/convertformat: src/convertformat.c $(MAINLIB)
	$(CXX) -o bin/convertformat src/convertformat.c -lm -Lsrc/lib $(MAINLINK)

bin/filterCRNs: src/filterCRNs.c $(MAINLIB)
	$(CXX) -o bin/filterCRNs src/filterCRNs.c -lm -Lsrc/lib $(MAINLINK)

bin/genbimol: src/genbimol.c $(MAINLIB)
	$(CXX) -o bin/genbimol src/genbimol.c -lm -Lsrc/lib $(MAINLINK)

bin/genbitrimol: src/genbitrimol.c $(MAINLIB)
	$(CXX) -o bin/genbitrimol src/genbitrimol.c -lm -Lsrc/lib $(MAINLINK)

bin/genzeroonetritrimol: src/genzeroonetritrimol.c $(MAINLIB)
	$(CXX) -o bin/genzeroonetritrimol src/genzeroonetritrimol.c -lm -Lsrc/lib $(MAINLINK)

bin/genbimultimol: src/genbimultimol.c $(MAINLIB)
	$(CXX) -o bin/genbimultimol src/genbimultimol.c -lm -Lsrc/lib $(MAINLINK)

bin/DSRanal: src/DSRanal.c $(MAINLIB)
	$(CXX) -o bin/DSRanal src/DSRanal.c -lm -Lsrc/lib $(MAINLINK)

bin/supCRN: src/supCRN.c $(MAINLIB)
	$(CXX) -o bin/supCRN src/supCRN.c -lm -Lsrc/lib $(MAINLINK)

bin/compareCRNlists: src/compareCRNlists.c $(MAINLIB)
	$(CXX) -o bin/compareCRNlists src/compareCRNlists.c -lm -Lsrc/lib $(MAINLINK)

bin/CRNwithsources: src/CRNwithsources.c $(MAINLIB)
	$(CXX) -o bin/CRNwithsources src/CRNwithsources.c -lm -Lsrc/lib $(MAINLINK)

bin/getallisomorphs: src/getallisomorphs.c $(MAINLIB)
	$(CXX) -o bin/getallisomorphs src/getallisomorphs.c -lm -Lsrc/lib $(MAINLINK)

bin/relconc: src/relconc.c $(MAINLIB)
	$(CXX) -o bin/relconc src/relconc.c -lm -Lsrc/lib $(MAINLINK)

bin/DF: src/DF.c $(MAINLIB)
	$(CXX) -o bin/DF src/DF.c -lm -Lsrc/lib $(MAINLINK)

