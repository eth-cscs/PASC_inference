#include "savevtk.h"

void save_VTK(Data data, Gamma gamma) 
{
	int t,k;
		
	int K = gamma.get_K();
	int T = gamma.get_T();
	int dim = data.get_dim();
		
	std::ostringstream oss_name_of_file;
    std::ofstream myfile;
	
	/* write to the name of file */
	oss_name_of_file << EXPORT_SAVEVTK_filename;

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
		myfile << data.data_vecs[0](t) << " "; /* x */
		myfile << data.data_vecs[1](t) << " "; /* y */
		myfile << " 0.0\n"; /* z */
	}
	myfile << "\n";

	/* values in points */
	myfile << "POINT_DATA " <<  T << "\n";
	/* prepare vector with idx of max values */
	GammaVector<int> gamma_max_idx(gamma.get_T());
	gamma_max_idx(all) = 0;
	for(k=0;k<K;k++){
		/* write gamma_k */
		myfile << "SCALARS gamma_" << k << " float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for(t=0;t<T;t++){
			myfile << gamma.gamma_vec(k*T + t) << "\n";

			/* update maximum */
			if(gamma.gamma_vec(k*T+t) > gamma.gamma_vec(gamma_max_idx(t)*T+t)){
				gamma.gamma_vec(k*T+t) = gamma.gamma_vec(gamma_max_idx(t)*T+t);
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

}
