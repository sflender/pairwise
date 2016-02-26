#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <cmath>
#include <math.h>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <random>
#include "cosmo.h"
#include "fitshandle.h"
#include "arr.h"
#include "fitsio.h"
#include "pairwise.h"
#include "matrix.h"
#include "config.h"

using namespace std;
int main(int argc, char** argv){

	if (argc != 2) {cout<<"Usage: ./pairwise config_file\n"; return -1;}

	// --- read parameters from the config file

	ConfigFile cfg(argv[1]);

	int N_bootstraps = cfg.getValueOfKey<int>("N_bootstraps");
	float zmin = cfg.getValueOfKey<float>("zmin");
	float zmax = cfg.getValueOfKey<float>("zmax");
	float massproxy_min = cfg.getValueOfKey<float>("massproxy_min");
	float massproxy_max = cfg.getValueOfKey<float>("massproxy_max");
	float sep_min = cfg.getValueOfKey<float>("sep_min");
	float sep_max = cfg.getValueOfKey<float>("sep_max");
	int Nbins = cfg.getValueOfKey<int>("Nbins");
	float h = cfg.getValueOfKey<float>("h");
	float Omega_M = cfg.getValueOfKey<float>("Omega_M");
	float Omega_k = cfg.getValueOfKey<float>("Omega_k");
	float wt = cfg.getValueOfKey<float>("wt");
	bool randomize_redshifts = cfg.getValueOfKey<bool>("randomize_redshifts");
	bool redshift_dependend_sigmaz = cfg.getValueOfKey<bool>("redshift_dependend_sigmaz");
	float sigma_z = cfg.getValueOfKey<float>("sigma_z");
	bool correct_z_evolution = cfg.getValueOfKey<bool>("correct_z_evolution");
	float sigma_K = cfg.getValueOfKey<float>("sigma_K");
	bool convert_units_from_K_to_mikroK = cfg.getValueOfKey<bool>("convert_units_from_K_to_mikroK");
	string input_name = cfg.getValueOfKey<string>("input_name");
	string output_name = cfg.getValueOfKey<string>("output_name");
	string filter_keyword = cfg.getValueOfKey<string>("filter_keyword");
	string massproxy_keyword = cfg.getValueOfKey<string>("massproxy_keyword");
	string redshift_keyword = cfg.getValueOfKey<string>("redshift_keyword");
	bool read_spt_mask = cfg.getValueOfKey<bool>("read_spt_mask");
	bool read_des_mask = 0; read_des_mask = cfg.getValueOfKey<bool>("read_des_mask");
	bool read_scaleval = cfg.getValueOfKey<bool>("read_scaleval");
	bool write_output = cfg.getValueOfKey<bool>("write_output");
	bool do_redshift_nulltest = cfg.getValueOfKey<bool>("do_redshift_nulltest");
	int scatter_in_masscuts = cfg.getValueOfKey<int>("scatter_in_masscuts");
	float scatter_in_masscuts_sigma = cfg.getValueOfKey<float>("scatter_in_masscuts_sigma");
	bool bin_in_redshift = cfg.getValueOfKey<bool>("bin_in_redshift");
	bool compute_distance_from_redshift = cfg.getValueOfKey<bool>("compute_distance_from_redshift");

	// --- parameters derived from the input
	float Delta_sep = (sep_max-sep_min)/Nbins;
	float Delta_zbin = (zmax-zmin)/Nbins;

	float H0 = 100*h;
	cosmo model(H0, Omega_M, Omega_k, wt);
	// --- parameters that will never change
	float deg2rad = M_PI/180.0;

	int rank=0;
	int nranks=1;

	if (rank==0) cout<<"...total number of bootstraps: "<<N_bootstraps*nranks<<endl;
	if (rank==0 and write_output) cout<<"...output file: "<<output_name<<endl;
	if (rank==0) cout<<"...reading input\n";

	fitsfile *fptr;
	int status=0;
	int hdutype;
	long int N_clus;

	fits_open_file(&fptr, (char*)input_name.c_str(), READONLY, &status);
   	if (status) {printf("unable to read files\n"); return -1;}
    fits_movabs_hdu(&fptr[0], 2, &hdutype, &status); //move to HDU 2 (that's where the data sits)
	fits_get_num_rows(&fptr[0], &N_clus, &status); // get the total number of rows (e.g. clusters) in the catalog

	double *ra = new double[N_clus];
   	read_column_from_fitsfile(fptr, (char*)"RA", ra, N_clus); 

	double *dec = new double[N_clus];
   	read_column_from_fitsfile(fptr, (char*)"DEC", dec, N_clus); 

	float *massproxy = new float[N_clus];
	read_column_from_fitsfile(fptr, (char*)massproxy_keyword.c_str(), massproxy, N_clus);

	float *redshift = new float[N_clus];
	read_column_from_fitsfile(fptr, (char*)redshift_keyword.c_str(), redshift, N_clus);

	float *comov_dist_fromfits = new float[N_clus];
	read_column_from_fitsfile(fptr, (char*)"COMOV_DIST_MPC", comov_dist_fromfits, N_clus);

	double *Tfil = new double[N_clus];
	read_column_from_fitsfile(fptr, (char*)filter_keyword.c_str(), Tfil, N_clus);

	float *scaleval = new float[N_clus]; fill(scaleval, scaleval+N_clus,1); // set default value of scaleval to 1
   	if (read_scaleval) read_column_from_fitsfile(fptr, (char*)"SCALEVAL", scaleval, N_clus); 

	short *cmb_valid = new short[N_clus]; fill(cmb_valid,cmb_valid+N_clus,1); // by default everyting is valid
   	if (read_spt_mask) read_column_from_fitsfile(fptr, (char*)"CMB_VALID", cmb_valid, N_clus); 

	short *cmbmap_psflag = new short[N_clus]; fill(cmbmap_psflag, cmbmap_psflag+N_clus,0);
   	if (read_spt_mask) read_column_from_fitsfile(fptr, (char*)"CMBMAP_PSFLAG", cmbmap_psflag, N_clus); 

   	float *desmask = new float[N_clus]; fill(desmask, desmask+N_clus, 1.0);
   	if(read_des_mask) read_column_from_fitsfile(fptr, (char*)"DESMASK", desmask, N_clus);

	if(convert_units_from_K_to_mikroK==1){
		for (int i=0; i<N_clus; i++) Tfil[i] *= 1e6;
	}

	// generate a random number seed (for bootstrapping and for nulltests)
	int seed = time(NULL)*(1+rank); // this is a random unique seed for each rank
	srand(seed); // each time the code is run it uses different random numbers on each rank. call this only once.


	// --- applying selection cuts
	vector<int> sel; // integer list of selected clusters
	int N_selected_clus=0;
	if (scatter_in_masscuts==0){ // select clusters without scatter in the mass-cut
		for (int i=0; i<N_clus; i++){
			if (cmb_valid[i]==1 and cmbmap_psflag[i]==0 and desmask[i]==1 and
				redshift[i]<zmax and redshift[i]>zmin and 
				massproxy[i]<massproxy_max and // upper masscut in lambda
				massproxy[i]/scaleval[i]>massproxy_min){ // more conservative lower cut in lambda/scaleval
				sel.push_back(i);
				N_selected_clus++;
			}
		}
	}

	else if (scatter_in_masscuts==1){ //scatter in lower masscut
		// we have logMtilde = logM + Delta, i.e. Mtilde = M * exp(Delta)
		std::default_random_engine generator(random_device{}());
	  	std::normal_distribution<float> distribution(0.0,scatter_in_masscuts_sigma); 

		for (int i=0; i<N_clus; i++){
			float Delta = distribution(generator); // generates a random scatter with width scatter_..._sigma
			if (cmb_valid[i]==1 and cmbmap_psflag[i]==0 and desmask[i]==1 and
				redshift[i]<zmax and redshift[i]>zmin and 
				massproxy[i]<massproxy_max and 
				massproxy[i] * exp(Delta) > massproxy_min){ 
				sel.push_back(i);
				N_selected_clus++;
			}
		}
	}

	else if (scatter_in_masscuts==2){ //scatter in both upper and lower cut
		std::default_random_engine generator(random_device{}());
	  	std::normal_distribution<float> distribution(0.0,scatter_in_masscuts_sigma); // generates a random scatter with width scatter_..._sigma
		for (int i=0; i<N_clus; i++){
			float Delta = distribution(generator); 
			if (cmb_valid[i]==1 and cmbmap_psflag[i]==0 and desmask[i]==1 and
				redshift[i]<zmax and redshift[i]>zmin and 
				massproxy[i] * exp(Delta) < massproxy_max and 
				massproxy[i] * exp(Delta) > massproxy_min){
				sel.push_back(i);
				N_selected_clus++;
			}
		}
	}
	

	if(rank==0) cout << "...N_clus = "<<N_clus<<endl;
	if(rank==0) cout << "...N_selected_clus = "<<N_selected_clus<<endl;
	
	// --- assign random errors to redshifts
	if (randomize_redshifts){
		float z_dependence = 1; // default: no z-dependence
		if(rank==0) cout<<"...randomizing the redshifts"<<endl;
		for (int i=0; i<N_selected_clus; i++){

			if (redshift_dependend_sigmaz){
				z_dependence = (1+redshift[sel[i]]); // in this case sigma_z will scale like ~(1+z)
			}

			std::default_random_engine generator(random_device{}());
	  		std::normal_distribution<float> distribution(0.0,sigma_z*z_dependence);

			float Delta_z = distribution(generator);
			redshift[sel[i]] += Delta_z;
		}
	}

/*
	//---obsolete content
	redshift_error_hybrid_model = cfg.getValueOfKey<bool>("redshift_error_hybrid_model");
	if (redshift_error_hybrid_model){
		cout<<"hybrid redshift error model!\n";
		int lowz_sample_size=0;
		int highz_sample_size=0;
		for (int i=0; i<N_selected_clus; i++){

			if (redshift[sel[i]]<0.45){
				std::default_random_engine generator(random_device{}());
	  			std::normal_distribution<float> distribution(0.0,0.003*(1+redshift[sel[i]]));
	  			float Delta_z = distribution(generator);
	  			redshift[sel[i]] += Delta_z;
	  			lowz_sample_size++;
				//cout<<i<<" "<<redshift[sel[i]]<<" "<<Delta_z<<endl;
			}

			else if (redshift[sel[i]]>=0.45){
				std::default_random_engine generator(random_device{}());
	  			std::normal_distribution<float> distribution(0.0,0.01*(1+redshift[sel[i]]));
	  			float Delta_z = distribution(generator);
	  			redshift[sel[i]] += Delta_z;
	  			highz_sample_size++;	
			    //cout<<i<<" "<<redshift[sel[i]]<<" "<<Delta_z<<endl;
			}
		}		
		cout<<N_selected_clus<<" "<<lowz_sample_size<<" "<<highz_sample_size<<endl;
		cout<<(float)lowz_sample_size/((float)N_selected_clus)<<" "<<(float)highz_sample_size/((float)N_selected_clus)<<endl;
	}
*/


	// --- destroying the redshift information (redsfhift null test) 
	if (do_redshift_nulltest){
		cout<<"...redshift null test!\n";
		for (int i=0; i<N_selected_clus; i++){
			redshift[sel[i]] = zmin + (float)(rand()%( (int)((zmax-zmin)*10000) ))/10000;
		}
	}


	// --- compute comov. distance to each cluster
	float *comov_dist = new float[N_clus];
	if(compute_distance_from_redshift){ // translate the input redshift into comoving distance
		cout << "...compute comov. distance to each selected cluster"<<endl;
		for (int i=0; i<N_selected_clus; i++){
			comov_dist[sel[i]] = model.comov_dist(redshift[sel[i]]); //comoving distance in units Mpc
		}
	}
	else{ // take the input distance directly
		for (int i=0; i<N_selected_clus; i++){
			comov_dist[sel[i]] = comov_dist_fromfits[sel[i]];
		}
	}


	// --- correct bias from z-dependent noise
	if(correct_z_evolution==1){
		float *curlyT = new float[N_selected_clus];
		if(rank==0) cout << "...correct for bias from redshift evolution"<<endl;
		for (int i=0; i<N_selected_clus;i++){
			float numerator=0;
			float denominator=0;
			for (int j=0; j<N_selected_clus; j++){
				float weight = exp(-0.5*pow( ((redshift[sel[i]]-redshift[sel[j]])/sigma_K) ,2));
				numerator += weight*Tfil[sel[j]];
				denominator += weight;
			}
			curlyT[i] = numerator/denominator;
		}

		for (int i=0; i<N_selected_clus; i++){
			Tfil[sel[i]] -= curlyT[i];
		}
		delete [] curlyT;
	}

	if(rank==0) cout << "...applying pairwise estimator"<<endl;


	cout<<"...computing the bootstrap weights\n";
	int **bootstrap_weight = new int*[N_bootstraps];
	for (int i=0; i<N_bootstraps; i++){
		bootstrap_weight[i] = new int[N_selected_clus]; fill(bootstrap_weight[i], bootstrap_weight[i]+N_selected_clus,0);
	}

	float **numerator = new float*[N_bootstraps]; 
	float **denominator = new float*[N_bootstraps]; 
	int **Npairs = new int*[N_bootstraps];
	float **That = new float*[N_bootstraps];
	for (int i=0; i<N_bootstraps; i++){
		numerator[i] = new float[Nbins]; fill(numerator[i],numerator[i]+Nbins,0.0);
		denominator[i] = new float[Nbins]; fill(denominator[i],denominator[i]+Nbins,0);
		Npairs[i] = new int[Nbins]; fill(Npairs[i],Npairs[i]+Nbins,0);
		That[i] = new float[Nbins]; fill(That[i], That[i]+Nbins,0.0);
	}


	for (int b=0; b<N_bootstraps; b++){
		for (int i=0; i<N_selected_clus; i++){
			int r = rand()%N_selected_clus; // pick a random index 
			bootstrap_weight[b][r]++; // increment bootstrap_weight at that index - cluster r is counted 'bootstrap_weight[b][r]-times' in bootstrap sample  b
		}
	}

	vector<int> pair_bin, partner_number;
	vector<float> pair_geom_weight, pair_Tdiff;

	int this_bin;

	for (int i=1; i<N_selected_clus; i++){
		int N_partners=0; // N_partners counts the number of valid partners (within the distance boundaries) of cluster i
		for (int j=0; j<i; j++){

			float costheta = dotprod(radec_to_vec(ra[sel[i]],dec[sel[i]]),radec_to_vec(ra[sel[j]],dec[sel[j]]));
			float rirjcostheta = dotprod(radec_to_vec(ra[sel[i]],dec[sel[i]])*comov_dist[sel[i]],radec_to_vec(ra[sel[j]],dec[sel[j]])*comov_dist[sel[j]]);
			float comov_sep = sqrt(pow(comov_dist[sel[i]],2) + pow(comov_dist[sel[j]],2) -2*comov_dist[sel[i]]*comov_dist[sel[j]]*costheta); // comoving 

			if (comov_sep>sep_min and comov_sep<sep_max){ // check if the pair is within the boundaries

				if(bin_in_redshift==0) {//binning in separation
					this_bin = floor((comov_sep-sep_min)/Delta_sep);
				}
				else { //binning in redshift
					this_bin = floor((redshift[sel[i]]+redshift[sel[j]]-2.0*zmin)/(2.0*Delta_zbin));
				}

				float pair_Tdiff = Tfil[sel[i]]-Tfil[sel[j]];		
				float geom_weight = (comov_dist[sel[i]]-comov_dist[sel[j]])*(1+costheta)/(2*comov_sep);

				for (int b=0; b<N_bootstraps; b++){
					numerator[b][this_bin] += geom_weight*pair_Tdiff*bootstrap_weight[b][i]*bootstrap_weight[b][j];
					denominator[b][this_bin] += geom_weight*geom_weight*bootstrap_weight[b][i]*bootstrap_weight[b][j];
					Npairs[b][this_bin] += bootstrap_weight[b][i]*bootstrap_weight[b][j];
				}
			}
		}
		cout<<"percent complete: "<<(float)i/((float)N_selected_clus)*100<<"\r";
		
	}
	cout<<endl;

	float *That_bar = new float[Nbins];
	float *mu_npairs = new float[Nbins];

	for (int k=0; k<Nbins; k++){
		That_bar[k]=0.0;
		for (int b=0; b<N_bootstraps; b++){
			That[b][k] = numerator[b][k]/denominator[b][k];
			That_bar[k]+=That[b][k];
			mu_npairs[k] += Npairs[b][k];
		}

		mu_npairs[k] /= (float)N_bootstraps;
		That_bar[k] /= N_bootstraps;
	}


	cout<<"...compute the covariance matrix\n";


	Matrix cov_matrix = Matrix(Nbins, Nbins);

	for (int k=0; k<Nbins; k++){
		for (int l=0; l<=k; l++){
			cov_matrix(k+1,l+1)=0;
			for (int b=0; b<N_bootstraps; b++){
				cov_matrix(k+1,l+1) += (That[b][k] - That_bar[k])*(That[b][l] - That_bar[l]);
			}
		cov_matrix(k+1,l+1) /= N_bootstraps;
		cov_matrix(l+1,k+1) = cov_matrix(k+1,l+1); // the covariance matrix is symmetric
		}
	}

	cout<<endl;
	for (int k=0; k<Nbins; k++){
		for (int l=0; l<Nbins; l++){
			cout<<cov_matrix(k+1,l+1)<<" ";
		}
		cout<<endl;
	}

	float chisq_null=0;
	bool compute_SNR = cfg.getValueOfKey<bool>("compute_SNR");
	if (compute_SNR){
		cout<<"...computing the inverse cov. matrix\n";
		Matrix inv_cov_matrix = Inv(cov_matrix);


		cout<<"...computing chisq_null \n";
		for (int k=0; k<Nbins; k++){
			for (int l=0; l<Nbins; l++){
				chisq_null += inv_cov_matrix(k+1,l+1)*That_bar[k]*That_bar[l];
			}
		}	
		cout<<"-------------------------\n";
		cout<<"   chisq_null = "<<chisq_null<<endl;
		cout<<"   sqrt(chisq_null) = "<<sqrt(chisq_null)<<endl;
		cout<<"-------------------------\n";
	}



	ofstream outfile;
	if (rank==0 and write_output){
		outfile.open(output_name);
	}

	cout<<endl;
	if(rank==0 and write_output) {outfile<<"# --- chisq_null = "<<chisq_null<<endl;}
	// --- output
	if(bin_in_redshift==0){
		for (int k=0; k<Nbins; k++){
			if (rank==0){
				cout<<sep_min+Delta_sep*k<<" "<<sep_min+Delta_sep*(k+1)<<" "<<" "<<That_bar[k]<<" "<<sqrt(cov_matrix(k+1,k+1))<<"   "<<mu_npairs[k]<<endl;
				if (write_output){
				outfile<<sep_min+Delta_sep*k<<" "<<sep_min+Delta_sep*(k+1)<<" "<<" "<<That_bar[k]<<" "<<sqrt(cov_matrix(k+1,k+1))<<"   "<<mu_npairs[k]<<endl;
				}
			}		
		}
	}
	else{
		for (int k=0; k<Nbins; k++){
			if (rank==0){
				cout<<zmin+Delta_zbin*k<<" "<<zmin+Delta_zbin*(k+1)<<" "<<" "<<That_bar[k]<<" "<<sqrt(cov_matrix(k+1,k+1))<<"   "<<mu_npairs[k]<<endl;
				if (write_output){
				outfile<<zmin+Delta_zbin*k<<" "<<zmin+Delta_zbin*(k+1)<<" "<<" "<<That_bar[k]<<" "<<sqrt(cov_matrix(k+1,k+1))<<"   "<<mu_npairs[k]<<endl;
				}
			}		
		}		
	}

	if (rank==0 and write_output){
		outfile.close();
	}

	cout<<endl;

	// --- free up the memory
	
//	delete [] ra, dec, redshift, Tfil;
//	for (int b=0; b<N_bootstraps; b++) delete [] That[b];
//	delete [] That;
	

	//MPI_Finalize();

	return 0;
}
