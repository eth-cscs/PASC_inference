
int read_int_from_file(std::ifstream &myfile) {
	int value;
	myfile.read((char *)&value, sizeof(int)); /* read block of memory */
//	value = __builtin_bswap32(value);
	return value;
}

double read_double_from_file(std::ifstream &myfile) {
	double value;
	myfile.read((char *)&value, sizeof(double)); /* read block of memory */
//	value = __builtin_bswap64(value);
	return value;
}

double *read_double_array_from_file(std::ifstream &myfile, const int size) {
	double *return_array = new double[size];
	int i;
	for(i = 0; i < size; i++){
		return_array[i] = read_double_from_file(myfile);
	}
	
	return return_array; /* deal with MATLAB 'ieee-be' */
}


void write_int_to_file(std::ofstream &myfile, int &value) {
//	int new_value = __builtin_bswap32(value);;
	int new_value = value;
	myfile.write((char *)&new_value, sizeof(int)); /* read block of memory */
}

void write_double_to_file(std::ofstream &myfile, double &value) {
//	double new_value = __builtin_bswap64(value);
	double new_value = value;
	myfile.write((char *)&new_value, sizeof(double)); /* read block of memory */
}

double *white_double_array_to_file(std::ofstream &myfile, double *darray, const int size) {
	int i;
	for(i = 0; i < size; i++){
		write_double_to_file(myfile, darray[i]);
	}
}
