CPP := mpicxx
CPPFLAGS := -O3 -g -std=gnu++0x -fpermissive

SRCDIR := ./src
OBJDIR := ./obj
MAIN := pairwise
INCLUDES := -I $(SRCDIR) \
			-I ../../software/Healpix_3.20/src/cxx/generic_gcc/include

OBJECTS := $(OBJDIR)/$(MAIN).o 

#linking to the libraries
# cfitios/lib
xtract : $(OBJECTS) 
	$(CPP) $(CPPFLAGS) $(OBJECTS) \
	-L../../software/Healpix_3.20/src/cxx/generic_gcc/lib \
	-L../../software/cfitsio \
	-lhealpix_cxx \
	-lcxxsupport \
	-lcfitsio \
	-L/usr/local/lib \
    -lgsl -lgslcblas  -lm \
	-o $(MAIN)

#compilation
$(OBJDIR)/$(MAIN).o: $(SRCDIR)/$(MAIN).cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/$(MAIN).cpp -o $(OBJDIR)/$(MAIN).o

clean:
	rm $(OBJDIR)/*.o 