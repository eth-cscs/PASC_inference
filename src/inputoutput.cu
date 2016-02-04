#include "inputoutput.h"

void InputOutput::saveVTK(std::string name_of_file, DataVector data_vec, GammaVector gamma_vec, int dim, int T, int K)
{
	Timer timer_saveVTK; 
	timer_saveVTK.restart();
	timer_saveVTK.start();
	
	int t,k;
		
	std::ostringstream oss_name_of_file;
    std::ofstream myfile;
	
	/* write to the name of file */
	oss_name_of_file << name_of_file;

	/* open file to write */
	myfile.open(oss_name_of_file.str().c_str());

	/* write header to file */
	myfile << "# vtk DataFile Version 3.1\n";
	myfile << "this is my kmeans data with solution\n";
	myfile << "ASCII\n";
	myfile << "DATASET UNSTRUCTURED_GRID\n";

	/* points - coordinates */
	myfile << "POINTS " << T << " FLOAT\n";
	for(t=0;t < T;t++){
		myfile << data_vec(t) << " "; /* x */
		myfile << data_vec(T+t) << " "; /* y */
		myfile << " 0.0\n"; /* z */
	}
	myfile << "\n";

	/* values in points */
	myfile << "POINT_DATA " <<  T << "\n";
	/* prepare vector with idx of max values */
	HostVector<int> gamma_max_idx(T);
	gamma_max_idx(all) = 0;
	for(k=0;k<K;k++){
		/* write gamma_k */
		myfile << "SCALARS gamma_" << k << " float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for(t=0;t<T;t++){
			myfile << gamma_vec(k*T + t) << "\n";

			/* update maximum */
			if(gamma_vec(k*T+t) > gamma_vec(gamma_max_idx(t)*T+t)){
				gamma_vec(k*T+t) = gamma_vec(gamma_max_idx(t)*T+t);
				gamma_max_idx(t) = k;
			}
		}
	}


	/* store gamma values */
	myfile << "SCALARS gamma_max_id float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for(t=0;t<T;t++){
		myfile << gamma_max_idx(t) << "\n";
	}

	/* close file */
	myfile.close();

	timer_saveVTK.stop();
	if(DEBUG_MODE >= 2) Message_info_time(" - problem saved to VTK in: ",timer_saveVTK.get_value_sum());
	

}
