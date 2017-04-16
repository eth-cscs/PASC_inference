#include "algebra/arrayoperation.h"

namespace pascinference {
namespace algebra {

template<class MyType>
std::string print_array(MyType *my_array, int my_size){
	std::ostringstream out;
	out << "[";
	for(int i=0;i<my_size;i++){
		out << my_array[i];
		if(i<my_size-1) out << ",";
	}	
	out << "]";
	return out.str();
}
template std::string print_array<int>(int *, int);
template std::string print_array<double>(double *, int);


template<class MyType>
std::string print_array(MyType *my_array, int my_size1, int my_size2){
	std::ostringstream out;
	out << "{";
	for(int i=0;i<my_size1;i++){
		out << "[";
		for(int j=0;j<my_size2;j++){
			out << my_array[i*my_size2+j];
			if(j<my_size2-1) out << ",";
		}
		out << "]";
		if(i<my_size1-1) out << ",";
	}	
	out << "}";
	return out.str();
}
template std::string print_array<int>(int *, int, int);
template std::string print_array<double>(double *, int, int);


bool parse_strings_to_doubles(int K, int Km, std::vector<std::string> Theta_list, double *Theta_solution) {
	std::string token; 
	size_t pos;
	int counter=0;
	
	for(int k=0;k < K;k++){
		pos = 0;
		
		while ((pos = Theta_list[k].find(",")) != std::string::npos) {
			token = Theta_list[k].substr(0, pos);

			if(counter >= K*Km) return false;
			Theta_solution[counter] = atof(token.c_str());
			counter++;

			Theta_list[k].erase(0, pos + 1);
		}

		if(counter >= K*Km) return false;
		Theta_solution[counter] = atof(Theta_list[k].c_str());
		counter++;
	}

	if(counter != K*Km){
		return false;
	} else {
		return true;
	}
}


template<class MyType>
std::string print_vector(std::vector<MyType> &my_vector){
	std::ostringstream out;
	out << "[";
	for(int i=0;i<my_vector.size();i++){
		out << my_vector[i];
		if(i<my_vector.size()-1) out << ",";
	}	
	out << "]";
	return out.str();
}
template std::string print_vector<int>(std::vector<int>&);
template std::string print_vector<double>(std::vector<double>&);


template<class MyType>
void print_array(std::ostream &output, int my_size, MyType *my_array){
	output << "[";
	for(int i=0;i<my_size;i++){
		output << my_array[i];
		if(i<my_size-1) output << ",";
	}	
	output << "]";
}
template void print_array<int>(std::ostream &, int, int *);
template void print_array<double>(std::ostream &, int, double *);


template<class MyType>
MyType sum_array(int my_size, const MyType *my_array){
	MyType sum = 0;
	for(int i=0;i<my_size;i++){
		sum += my_array[i];
	}	
	return sum;
}
template int sum_array<int>(int, const int*);
template double sum_array<double>(int, const double*);


template<class MyType>
MyType max_array(int my_size, const MyType *my_array){
	MyType return_value = my_array[0];
	for(int i=1;i<my_size;i++){
		if(my_array[i] > return_value){
			return_value = my_array[i];
		}
	}	
	return return_value;
}
template int max_array<int>(int, const int*);
template double max_array<double>(int, const double*);


template<class MyType>
int max_arg_array(int my_size, const MyType *my_array){
	MyType max_value = my_array[0];
	int return_value;
	for(int i=1;i<my_size;i++){
		if(my_array[i] > return_value){
			max_value = my_array[i];
			return_value = i;
		}
	}	
	return return_value;
}
template int max_arg_array<int>(int, const int*);
template int max_arg_array<double>(int, const double*);


template<class MyType>
void set_value_array(int my_size, MyType *my_array, MyType my_value){
	for(int i=0;i<my_size;i++){
		my_array[i] = my_value;
	}	
}
template void set_value_array<int>(int, int *, int );
template void set_value_array<double>(int, double *, double );


template<class MyType>
MyType max_diff_array(int my_size, const MyType *my_array){
	MyType return_value = my_array[1] - my_array[0];
	for(int i=2;i<my_size;i++){
		if((my_array[i]-my_array[i-1]) > return_value){
			return_value = my_array[i]-my_array[i-1];
		}
	}	
	return return_value;
}
template int max_diff_array<int>(int, const int *);
template double max_diff_array<double>(int, const double *);


template<class MyType>
MyType sum_subarray(int start, int end, const MyType *my_array){
	MyType sum = 0;
	for(int i=start;i<=end;i++){
		sum += my_array[i];
	}	
	return sum;
}
template int sum_subarray<int>(int, int, const int *);
template double sum_subarray<double>(int, int, const double *);


template<class MyType>
MyType dot_arrays(int size, const MyType *my_array1, const MyType *my_array2){
	MyType sum = 0;
	for(int i=0;i<=size;i++){
		sum += my_array1[i]*my_array2[i];
	}	
	return sum;
}
template int dot_arrays<int>(int, const int *, const int *);
template double dot_arrays<double>(int, const double *, const double *);


template<class MyType>
MyType max(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 > a1){
		return_value = a2;
	}
	return return_value;
}
template int max<int>(const int, const int);
template double max<double>(const double, const double);


template<class MyType>
MyType min(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 < a1){
		return_value = a2;
	}
	return return_value;
}
template int min<int>(const int, const int);
template double min<double>(const double, const double);


template<class MyType>
void mult_pw_array(int my_size, MyType *my_array, const MyType *my_array1, const MyType *my_array2){
	for(int i=0;i<my_size;i++){
		my_array[i] = my_array1[i]*my_array2[i];
	}	
}
template void mult_pw_array<int>(int, int *, const int *, const int *);
template void mult_pw_array<double>(int, double *, const double *, const double *);


}
} /* end of namespace */
