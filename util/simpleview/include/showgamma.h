


void show_gamma(std::string filename){
	std::cout << " - show gamma from file " << filename << std::endl;

	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	int K = read_int_from_file(myfile);
	int T = read_int_from_file(myfile);
	int t,k;

	Gnuplot gp;

	/* header */
	gp << "set terminal qt size 1000,400 \n"; 
	gp << "set terminal x11 size 1000,400 \n"; 
	gp << "set title \"Gamma results\"\n";
	gp << "set xlabel \"t\"\n";
	gp << "set ylabel \"gamma_k(t)\"\n";
	gp << "set xrange [0:"<< T << "]\nset yrange [0:1]\n";

	/* set colors */
	gp << "set linetype 1 lc rgb 'red' 			lw 4\n";
	gp << "set linetype 2 lc rgb 'blue'			lw 4\n";
	gp << "set linetype 3 lc rgb 'green'		lw 4\n";
	gp << "set linetype 4 lc rgb 'orange'		lw 4\n";
	gp << "set linetype 5 lc rgb 'dark-violet'	lw 4\n";
	gp << "set linetype 6 lc rgb 'sea-green' 	lw 4\n";
	gp << "set linetype 7 lc rgb 'cyan' 		lw 4\n";
	gp << "set linetype 8 lc rgb 'dark-red'		lw 4\n";
	gp << "set linetype 9 lc rgb 'blue'			lw 4\n";
	gp << "set linetype 10 lc rgb 'dark-orange'	lw 4\n";
	gp << "set linetype 11 lc rgb 'goldenrod'	lw 4\n";
	gp << "set linetype cycle 11\n";
	

	/* plot header */
	gp << "plot ";
	for(k=0;k<K;k++){
		gp << "'-' using 1:2 title 'gamma_" << k << "' with lines";
		if(k < K-1) gp << ", ";
	}
	gp << "\n";	

	for(k=0;k<K;k++){
		for(t=0;t<T;t++){
			gp << "\t" << t << " " << read_double_from_file(myfile) << "\n";
		}
	
		gp << "e\n";
	}

    myfile.close();
	
}


