
void generate_sample(std::string filename, int T, int K){
	std::cout << " - generate sample to: " << filename << std::endl;

	/* open file */
	std::ofstream myfile(filename.c_str(), std::ios::out | std::ios::binary);

	/* header */
	write_int_to_file(myfile, K);
	write_int_to_file(myfile, T);

	int t,k;
	double value;
	
	for(k=0;k<K;k++){
		for(t=0;t<T;t++){
			value = 0.5*cos(t*PI/180.0 + k*PI/5.0)+0.5;
			write_double_to_file(myfile, value);
		}
	}
	
	/* close file */
    myfile.close();	
	
}

