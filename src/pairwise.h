#include <math.h>
#include "healpix_map.h"
#include "fitsio.h"

using namespace std;

vec3 radec_to_vec(float ra, float dec){ // converts RA/DEC to a vec3 element
	float deg2rad = M_PI/180.0;
	float theta_healpix = (90.0 - dec)*deg2rad;
	float phi_healpix = ra * deg2rad;
	pointing point(theta_healpix,phi_healpix);
	return point.to_vec3();
}

int read_column_from_fitsfile(fitsfile *fptr, char* column_name, double* column_data, long int N_clus){
	// some definitions for cfitsio
	int status=0;
	int casesen = 0; // case sensitivity off
	int colnum;
	int typecode;

	// read in the data
	fits_get_colnum(&fptr[0], casesen, column_name, &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, NULL, NULL, &status);
	if (typecode==82){
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, column_data, NULL, &status);
	}
	else if (typecode==42){
		cout<<column_name<<" is a float, not a double - but I will cast it for you.\n";
		float *array = new float[N_clus];
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, array, NULL, &status);
		for (int i=0;i<N_clus;i++){column_data[i]=(double)array[i];}
		delete [] array;
	}
	else cout<<"Error: "<<column_name<<" is neither a double nor float!\n"; return -1;
	return 0;
}

int read_column_from_fitsfile(fitsfile *fptr, char* column_name, float* column_data, long int N_clus){
	// some definitions for cfitsio
	int status=0;
	int casesen = 0; // case sensitivity off
	int colnum;
	int typecode;

	// read in the data
	fits_get_colnum(&fptr[0], casesen, column_name, &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, NULL, NULL, &status);
	if (typecode==42){
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, column_data, NULL, &status);
	}
	else if (typecode==82){
		cout<<column_name<<" is a double, not a float - but I will cast it for you.\n";
		double *array = new double[N_clus];
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, array, NULL, &status);
		for (int i=0;i<N_clus;i++){column_data[i]=(float)array[i];}
		delete [] array;
	}
	else cout<<"Error: "<<column_name<<" is neither a double nor float!\n"; return -1;
	return 0;
}

int read_column_from_fitsfile(fitsfile *fptr, char* column_name, int* column_data, long int N_clus){
	// some definitions for cfitsio
	int status=0;
	int casesen = 0; // case sensitivity off
	int colnum;
	int typecode;

	// read in the data
	fits_get_colnum(&fptr[0], casesen, column_name, &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, NULL, NULL, &status);
	if (typecode==31){
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, column_data, NULL, &status);
	}
	else if (typecode==21){
		cout<<column_name<<" is a short, not an int - but I will cast it for you.\n";
		short *array = new short[N_clus];
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, array, NULL, &status);
		for (int i=0;i<N_clus;i++){column_data[i]=(int)array[i];}
		delete [] array;
	}
	else if (typecode==20){
		cout<<column_name<<" is an unsigned short, not an int - but I will cast it for you.\n";
		unsigned short *array = new unsigned short[N_clus];
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, array, NULL, &status);
		for (int i=0;i<N_clus;i++){column_data[i]=(int)array[i];}
		delete [] array;
	}
	else cout<<"Error: "<<column_name<<" is not the type you think it is -- it is typecode"<<typecode<<"\n"; return -1;
	return 0;
}

int read_column_from_fitsfile(fitsfile *fptr, char* column_name, short* column_data, long int N_clus){
	// some definitions for cfitsio
	int status=0;
	int casesen = 0; // case sensitivity off
	int colnum;
	int typecode;

	// read in the data
	fits_get_colnum(&fptr[0], casesen, column_name, &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, NULL, NULL, &status);
	if (typecode==21){
		fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, column_data, NULL, &status);
	}
	else cout<<"Error: "<<column_name<<" is not the type you think it is!\n"; return -1;
	return 0;
}