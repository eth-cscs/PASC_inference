#ifndef PASC_COMMON_ARRAYOPERATION_H
#define	PASC_COMMON_ARRAYOPERATION_H

namespace pascinference {
	
template<class MyType>
void print_array(std::ostream &output, int my_size, MyType *my_array){
	output << "[";
	for(int i=0;i<my_size;i++){
		output << my_array[i];
		if(i<my_size-1) output << ",";
	}	
	output << "]";
}

template<class MyType>
MyType sum_array(int my_size, const MyType *my_array){
	MyType sum = 0;
	for(int i=0;i<my_size;i++){
		sum += my_array[i];
	}	
	return sum;
}

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

template<class MyType>
MyType sum_subarray(int start, int end, const MyType *my_array){
	MyType sum = 0;
	for(int i=start;i<=end;i++){
		sum += my_array[i];
	}	
	return sum;
}

template<class MyType>
MyType dot_arrays(int size, const MyType *my_array1, const MyType *my_array2){
	MyType sum = 0;
	for(int i=0;i<=size;i++){
		sum += my_array1[i]*my_array2[i];
	}	
	return sum;
}

template<class MyType>
MyType max(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 > a1){
		return_value = a2;
	}
	return return_value;
}

template<class MyType>
MyType min(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 < a1){
		return_value = a2;
	}
	return return_value;
}


/* my_array[i] = my_array1[i]*my_array2[i] */
template<class MyType>
void mult_pw_array(int my_size, MyType *my_array, const MyType *my_array1, const MyType *my_array2){
	for(int i=0;i<my_size;i++){
		my_array[i] = my_array1[i]*my_array2[i];
	}	
}


} /* end of namespace */

#endif
