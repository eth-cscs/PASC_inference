
namespace pascinference {

void InputOutput::saveVTK(std::string name_of_file, DataVector data_vec, GammaVector gamma_vec, int dim, int T, int K)
{
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - generating VTK" << std::endl;
 	
	Timer timer_saveVTK; 
	timer_saveVTK.restart();
	timer_saveVTK.start();
	
	int t,k;
		
	std::ostringstream oss_name_of_file;
    std::ofstream myfile;
	
	/* write to the name of file */
	oss_name_of_file << name_of_file;

	/* open file to write */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - open file" << std::endl;
	myfile.open(oss_name_of_file.str().c_str());

	/* write header to file */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - write header" << std::endl;
	myfile << "# vtk DataFile Version 3.1\n";
	myfile << "this is my kmeans data with solution\n";
	myfile << "ASCII\n";
	myfile << "DATASET UNSTRUCTURED_GRID\n";

	/* points - coordinates */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - write points - coordinates" << std::endl;
	myfile << "POINTS " << T << " FLOAT\n";
	for(t=0;t < T;t++){
		myfile << data_vec(t) << " "; /* x */
		myfile << data_vec(T+t) << " "; /* y */
		myfile << " 0.0\n"; /* z */
	}
	myfile << "\n";

	/* values in points */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - write values in points" << std::endl;
	myfile << "POINT_DATA " <<  T << "\n";
	/* prepare vector with idx of max values */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - prepare gamma_max_idx" << std::endl;

	GammaVector gamma_max(T);
	GammaVector  temp;
	gamma_max(petscvector::all) = 0.0; // TODO: deal with all
	GammaVector gamma_max_idx(T); // TODO: use general host vecotr
	
	gamma_max_idx(gall) = 0; // TODO: deal with all
	for(k=0;k<K;k++){
		/* write gamma_k */
		myfile << "SCALARS gamma_" << k << " float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for(t=0;t<T;t++){
			myfile << gamma_vec(k*T + t) << "\n";

			/* update maximum */
			if(gamma_vec(k*T+t) > gamma_max(t)){
				DEBUG_MODE = 101;
				temp = gamma_vec(k*T+t);
				gamma_max(t) = temp;
				
				DEBUG_MODE = 10;
//				gamma_vec(k*T+t) = gamma_vec(gamma_max_idx(t)*T+t);
				gamma_max_idx(t) = k;
			}
		}
	}


	/* store gamma values */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - store gamma values" << std::endl;
	myfile << "SCALARS gamma_max_id float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for(t=0;t<T;t++){
		myfile << gamma_max_idx(t) << "\n";
	}

	/* close file */
	if(DEBUG_MODE >= 11) coutMaster << offset <<" - close file" << std::endl;
	myfile.close();

	timer_saveVTK.stop();
	if(DEBUG_MODE >= 3) Message_info_time(" - problem saved to VTK in: ",timer_saveVTK.get_value_sum());
	

}


} /* end of namespace */
