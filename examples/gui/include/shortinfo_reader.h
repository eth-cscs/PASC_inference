#ifndef SHORTINFO_READER_H
#define	SHORTINFO_READER_H

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

void remove_extra_whitespaces(const std::string &input, std::string &output){
	output.clear();
	unique_copy(input.begin(), input.end(), std::back_insert_iterator<std::string>(output),
		[](char a,char b){ return std::isspace(a) && std::isspace(b);});  
}

class Shortinfo_reader {
	private:
		class col_node {
			private:
				std::string name; /* header */
				std::vector<double> values_double;
				std::vector<std::string> values_string;
			
			public:
				col_node(std::string new_name){
					remove_extra_whitespaces(new_name, name);
				}

				void print(){
					std::cout << "header: " << name << std::endl;
					std::cout << "values:" << std::endl;

					for(int i=0;i<values_double.size();i++){
						std::cout << "   " << i << ". " << values_double[i] << " (double)" << std::endl;
					}
					
					for(int i=0;i<values_string.size();i++){
						std::cout << "   " << i << ". " << values_string[i] << " (string)" << std::endl;
					}
				}
				
				void add_value(std::string value_string){
					double value_double;
					try {
						value_double = std::stod(value_string);
						values_double.push_back(value_double);
					} catch (const std::invalid_argument&) {
						values_string.push_back(value_string);
					}
				}

				bool is_me(std::string find_name){
					bool return_value;

					std::string find_name_pure;
					remove_extra_whitespaces(find_name, find_name_pure);

					if(find_name_pure == name){
						return_value = true;
					} else {
						return_value = false;
					}
					
					return return_value;
				}
				
				std::vector<double> get_values_double(){
					return values_double;
				}
				
				std::vector<std::string> get_values_string(){
					return values_string;
				}

				std::string get_header(){
					return name;
				}

		};

		bool is_opened;
		std::ifstream *myfile;						/**< the file reader */
		int line_counter; 							/**< number of lines in the file */
		std::vector<col_node> cols;					/**< vector of founded columns with values */
	
		void process_header(std::string myline){
			std::vector<std::string> words;
			boost::split(words, myline, boost::is_any_of(","), boost::token_compress_on);

			for(int i=0;i<words.size();i++){
				cols.push_back(col_node(words[i]));
			}
		}

		void process_line(std::string myline){
			std::vector<std::string> words;
			boost::split(words, myline, boost::is_any_of(","), boost::token_compress_on);

			for(int i=0;i<words.size();i++){
				cols[i].add_value(words[i]);
			}
		}
		
	public:
		Shortinfo_reader(std::string filename){
			myfile = new std::ifstream(filename);
			line_counter = 0;
			
			if (myfile->is_open()){ /* if the file can be opened */
				is_opened = true;
			} else {
				is_opened = false;
			}
		}
		
		void process(){
			if(is_opened){
				int line_counter = 0;
				
				std::string myline; /* one loaded line from the file */
				while ( getline (*myfile,myline) ){ /* we are reading one line */
					/* process the line */
					if(line_counter == 0){
						process_header(myline);
					} else {
						process_line(myline);
					}

					/* increase counter of lines */
					line_counter++;
				}
			}

		}
		
		void print(){
			for(int i=0;i<cols.size();i++){
				cols[i].print();
				std::cout << "-------------------------------------------------------" << std::endl;
			}
		}

		void print_headers(){
			for(int i=0;i<cols.size();i++){
				std::cout << cols[i].get_header() << std::endl;
			}
		}
		
		std::vector<double> get_values_double(std::string find_name){
			for(int i=0;i<cols.size();i++){
				if(cols[i].is_me(find_name)){
					return cols[i].get_values_double();
				}
			}
		}

		std::vector<std::string> get_values_string(std::string find_name){
			for(int i=0;i<cols.size();i++){
				if(cols[i].is_me(find_name)){
					return cols[i].get_values_string();
				}
			}
		}

		std::vector<std::string> get_headers(){
			std::vector<std::string>  return_vector;
			for(int i=0;i<cols.size();i++){
				return_vector.push_back(cols[i].get_header());
			}
			return return_vector;
		}


};

#endif
