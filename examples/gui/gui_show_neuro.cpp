#include "pascinference.h"

#ifndef USE_PETSC
 #error 'this util is for PETSC'
#endif

#ifndef USE_DLIB
 #error 'this util is for DLIB'
#endif

#define ENABLE_ASSERTS

#include <dlib/gui_widgets.h>
#include <sstream>
#include <string>

#define NUMBER_OF_LABELS 14
#define DEFAULT_RADIUS_RATIO 30
#define DRAW_SPACE 0.1
#define PIXEL_SIZE 10.0

#define DEFAULT_TYPE 0

#define COLOR_MIN_R 0
#define COLOR_MIN_G 0
#define COLOR_MIN_B 255

#define COLOR_MAX_R 0
#define COLOR_MAX_G 255
#define COLOR_MAX_B 0

#define VALUE_MIN -200
#define VALUE_MAX 200

using namespace dlib;
using namespace pascinference;


class eegimageplotter : public drawable {
	private:
		std::string data_filename;
		Vec *data_Vec;
		bool data_loaded;

		std::string graph_filename;
		Vec *graph_Vec;
		bool graph_loaded;

		void draw (const canvas& c) const;
		void plot_image(const canvas& c) const;
		void plot_graph(const canvas& c) const;
		void plot_data(const canvas& c) const;
		void plot_graph_data(const canvas& c) const;

		int type;

		int T;
		int t;
		int npoints;

		double radius_ratio;

		void compute_mapping();
		double ax, bx, ay, by; /* cooeficients of linear mapping */
	public: 
		eegimageplotter(drawable_window& w, int type, double radius_ratio);
		~eegimageplotter();	
		
		void set_plotting_area(rectangle area);
		
		void load_data( const std::string& filename );
		bool get_data_loaded();
		std::string get_data_filename() const;
		Vec *get_data_Vec();

		void load_graph( const std::string& filename );
		bool get_graph_loaded();
		std::string get_graph_filename() const;
		Vec *get_graph_Vec();
	
		int get_T() const;
		int get_npoints() const;
		
		int set_t(int new_t);
};

class show_image_window : public drawable_window {
private:
    menu_bar mbar; /* gui: menubar */
	eegimageplotter myeegimageplotter; /* gui: canvas for drawing */
 
    label **labels_properties;
	scroll_bar *timescroll;
	label t_label;

    void load_data( const std::string& filename );
    void load_graph( const std::string& filename );

	/* menu events */
    void on_menu_file_open_data ();
    void on_menu_file_open_graph ();
    void on_menu_file_quit ();
    void on_menu_help_about();
    void on_timescroll();
	
	/* general window events */
	void on_window_resized();
	
	template<class ValueType> void set_label_properties(int label_idx, const std::string &text, ValueType value);
	void set_label_properties(int label_idx, const std::string &text);
	std::string cut_filename(const std::string &input);
	void fill_labels();

	int t;

public:
	show_image_window(int type, double radius_ratio);
	show_image_window(std::string data_filename, std::string graph_filename, std::string title, int type, double radius_ratio);

	~show_image_window();

};

//  ----------------------------------------------------------------------------

int main( int argc, char *argv[] ) {
	/* add local program options */
	boost::program_options::options_description opt_problem("GUI_SHOW_NEURO", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("data_filename", boost::program_options::value<std::string>(), "name of input file [string]")
		("graph_filename", boost::program_options::value<std::string>(), "name of input file with coordinates [string]")
		("type", boost::program_options::value<int>(), "type of plotting [0=color vertices,1=interpolated]")
		("radius_ratio", boost::program_options::value<double>(), "how large is radius in comparison with max(width,height) [double]")
		("title", boost::program_options::value<std::string>(), "title of window [string]")	;
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* get values provided as console parameters */
	std::string data_filename, graph_filename, title;
	int type;
	double radius_ratio;

	consoleArg.set_option_value("data_filename", &data_filename, "");
	consoleArg.set_option_value("graph_filename", &graph_filename, "");
	consoleArg.set_option_value("title", &title, "");
	consoleArg.set_option_value("radius_ratio", &radius_ratio, DEFAULT_RADIUS_RATIO);
	consoleArg.set_option_value("type", &type, DEFAULT_TYPE);


    // create our window
    show_image_window my_window(data_filename, graph_filename, title, type, radius_ratio);

    // wait until the user closes this window before we let the program 
    // terminate.
    my_window.wait_until_closed();

	/* finalize the library */
	Finalize<PetscVector>();
	
    return 0;
}


/* ---------------- implementation -------------- */

show_image_window::show_image_window(int type, double radius_ratio) : /* All widgets take their parent window as an argument to their constructor. */
        mbar(*this),
        t_label(*this),
        myeegimageplotter(*this,type,radius_ratio)
{

	/* set the size of window */
    set_size(600,350);
		
	/* allocate labels */
	labels_properties = new label*[NUMBER_OF_LABELS];
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		labels_properties[i] = new label(*this);
	}

	/* main timer */
	t = 0;
	myeegimageplotter.set_t(t);
    t_label.set_pos(10,50);

	/* prepare time scroll_bar */
	timescroll = new scroll_bar(*this, scroll_bar::HORIZONTAL);
	timescroll->set_scroll_handler( *this, &show_image_window::on_timescroll );
	timescroll->set_length(180);
	timescroll->set_pos(10,80);
	timescroll->disable();

    /* prepare position of labels of vector properties */
    labels_properties[0]->set_pos(10,110);
	for(int i=1; i < NUMBER_OF_LABELS; i++){
		labels_properties[i]->set_pos(labels_properties[i-1]->left(),labels_properties[i-1]->bottom()+20);
	}

	/* window title */
	set_title("Show PETSc Neuro utility");
        
	/* create menu bar */
    mbar.set_number_of_menus(2);
    mbar.set_menu_name(0,"File",'F');
    mbar.set_menu_name(1,"Help",'H');

    /* add the entries to the File menu. */
    mbar.menu(0).add_menu_item(menu_item_text("Open Graph",   *this, &show_image_window::on_menu_file_open_graph,    'G'));
    mbar.menu(0).add_menu_item(menu_item_text("Open Data",   *this, &show_image_window::on_menu_file_open_data,    'O'));
    mbar.menu(0).add_menu_item(menu_item_separator());
    mbar.menu(0).add_menu_item(menu_item_text("Quit",   *this, &show_image_window::on_menu_file_quit,    'Q'));

    /* Add the entries to the Help menu. */
    mbar.menu(1).add_menu_item(menu_item_text("About",  *this, &show_image_window::on_menu_help_about,   'A'));

	/* arrange the window */
	on_window_resized();

	fill_labels();
	show();
} 

show_image_window::show_image_window(std::string data_filename, std::string graph_filename, std::string title, int type, double radius_ratio) : show_image_window(type,radius_ratio) {
	if(title != ""){
		set_title(title);		
	}

	if(data_filename != ""){
		load_data( data_filename );
	}

	if(graph_filename != ""){
		load_graph( graph_filename );
	}
}

show_image_window::~show_image_window(){

	delete timescroll;
	
	/* destroy labels */
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		free(labels_properties[i]);
	}
	free(labels_properties);

	/* close window */
	close_window();
}

void show_image_window::on_menu_help_about(){
     message_box("About","This application is for PETSc EEG visualisation\n");
}

void show_image_window::on_menu_file_open_data(){
    /* display a file chooser window and when the user choses a file */
    open_existing_file_box(*this, &show_image_window::load_data);
}

void show_image_window::on_menu_file_open_graph(){
    /* display a file chooser window and when the user choses a file */
    open_existing_file_box(*this, &show_image_window::load_graph);
}

void show_image_window::on_menu_file_quit(){
	close_window();
}

void show_image_window::on_window_resized() {

	/* set new plotting area */
	unsigned long width,height;
    get_size(width,height);
	myeegimageplotter.set_plotting_area(rectangle(200,mbar.bottom(),width,height));
	drawable_window::on_window_resized();

	
}

void show_image_window::on_timescroll() {
	t = timescroll->slider_pos();
	myeegimageplotter.set_t(t);

	std::ostringstream sout;
	sout << "t = " << t << " (" << myeegimageplotter.get_T() << ")";
	t_label.set_text(sout.str());
}

void show_image_window::load_data( const std::string& filename ) {
	myeegimageplotter.load_data(filename);
	fill_labels();
}

void show_image_window::load_graph( const std::string& filename ) {
	myeegimageplotter.load_graph(filename);
	fill_labels();
}

void show_image_window::fill_labels() {
	/* compute basic properties of loaded vector */
	int graph_size;
	int data_size;
	double data_norm2;
	double data_max;
	double data_min;

	/* print properties of vectors to labels */
	set_label_properties(0, "GRAPH");
	set_label_properties(1, " loaded:    ", print_bool(myeegimageplotter.get_graph_loaded()));
	if(myeegimageplotter.get_graph_loaded()){
		Vec *graph_Vec = myeegimageplotter.get_graph_Vec();
		set_label_properties(2, " name:      ", cut_filename(myeegimageplotter.get_graph_filename()));
		set_label_properties(3, " npoints:   ", myeegimageplotter.get_npoints());

		TRYCXX( VecGetSize(*graph_Vec, &graph_size) );
		set_label_properties(4, " size:   ", graph_size);
	}
	set_label_properties(5, "");


	set_label_properties(6, "DATA");
	set_label_properties(7, " loaded:    ", print_bool(myeegimageplotter.get_data_loaded()));
	if(myeegimageplotter.get_data_loaded()){
		Vec *data_Vec = myeegimageplotter.get_data_Vec();
		set_label_properties(8, " name:      ", cut_filename(myeegimageplotter.get_data_filename()));
		set_label_properties(9, " T:         ", myeegimageplotter.get_T());

		TRYCXX( VecGetSize(*data_Vec, &data_size) );
		set_label_properties(10, " size:      ", data_size);

		TRYCXX( VecNorm(*data_Vec, NORM_2, &data_norm2) );
		set_label_properties(11, " norm2:     ", data_norm2);

		TRYCXX( VecMax(*data_Vec, NULL, &data_max) );
		set_label_properties(12, " max:       ", data_max);

		TRYCXX( VecMin(*data_Vec, NULL, &data_min) );
		set_label_properties(13, " min:       ", data_min);
	}
	
	if(myeegimageplotter.get_graph_loaded() && myeegimageplotter.get_data_loaded()){
		if(t > myeegimageplotter.get_T()-1){
			t = myeegimageplotter.get_T()-1;
		}

		timescroll->set_max_slider_pos(myeegimageplotter.get_T()-1);
		timescroll->set_slider_pos(t);
		timescroll->set_jump_size(timescroll->max_slider_pos()/10.0);
		timescroll->enable();

		std::ostringstream sout;
		sout << "t = " << t << " (" << myeegimageplotter.get_T() << ")";
		t_label.set_text(sout.str());
//		myeegimageplotter.set_t(t);

	} else {
		t = 0;
		t_label.set_text("");
		timescroll->disable();
//		myeegimageplotter.set_t(0);
	}

}

template<class ValueType>
void show_image_window::set_label_properties(int label_idx, const std::string &text, ValueType value){
	std::ostringstream sout;
	sout << std::setprecision(5);
    sout << text << std::setw(20) << value;
    labels_properties[label_idx]->set_text(sout.str());	
}

void show_image_window::set_label_properties(int label_idx, const std::string &text){
    labels_properties[label_idx]->set_text(text);	
}

std::string show_image_window::cut_filename(const std::string &input){
	std::ostringstream sout;
	boost::filesystem::path p(input);
	sout << p.filename();

	return sout.str();
}




/* ------------------------- eeg image plotter -------------- */

void eegimageplotter::draw(const canvas& c) const {
	/* draw background */
	fill_rect(c,rect,rgb_pixel(255,255,255));
	
	/* plot image */
	plot_image(c);

}

eegimageplotter::eegimageplotter(drawable_window& w, int type, double radius_ratio): 
			drawable(w)
{
	/* create empty vectors */
	data_Vec = new Vec();
	TRYCXX( VecCreate(PETSC_COMM_WORLD,data_Vec) );
	data_loaded = false;

	graph_Vec = new Vec();
	TRYCXX( VecCreate(PETSC_COMM_WORLD,graph_Vec) );
	graph_loaded = false;
	
	T = 0;
	npoints = 0;
	
	this->type = type;
	this->radius_ratio = radius_ratio;

	enable_events();
}

eegimageplotter::~eegimageplotter(){
	free(data_Vec);
	free(graph_Vec);

	disable_events();
	parent.invalidate_rectangle(rect);
}

void eegimageplotter::set_plotting_area(rectangle area){
	rect = area;
	compute_mapping();
}

void eegimageplotter::load_data( const std::string& filename ) {
	this->data_filename = filename;

	if(data_loaded){
		/* destroy existing vector */
		TRYCXX( VecDestroy(data_Vec) );
		TRYCXX( VecCreate(PETSC_COMM_WORLD,data_Vec) );
	}

	/* load the data from file */
	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_SELF, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_SELF ,this->data_filename.c_str(), FILE_MODE_READ, &mviewer) );

	/* load vector from viewer */
	TRYCXX( VecLoad(*data_Vec, mviewer) );

	/* destroy the viewer */
	TRYCXX( PetscViewerDestroy(&mviewer) );

	data_loaded = true;

	/* recompute T */
	int data_size;
	TRYCXX( VecGetSize(*data_Vec, &data_size) );
	if(graph_loaded){
		T = (int)(data_size/(double)get_npoints());
	} else {
		T = data_size;
	}

	compute_mapping();

	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);
}

void eegimageplotter::load_graph( const std::string& filename ) {
	this->graph_filename = filename;

	if(graph_loaded){
		/* destroy existing vector */
		TRYCXX( VecDestroy(graph_Vec) );
		TRYCXX( VecCreate(PETSC_COMM_WORLD, graph_Vec) );
	}

	/* load the data from file */
	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_SELF, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_SELF ,this->graph_filename.c_str(), FILE_MODE_READ, &mviewer) );

	/* load vector from viewer */
	TRYCXX( VecLoad(*graph_Vec, mviewer) );

	/* destroy the viewer */
	TRYCXX( PetscViewerDestroy(&mviewer) );

	graph_loaded = true;

	/* recompute npoints */
	int graph_size;
	TRYCXX( VecGetSize(*graph_Vec, &graph_size) );
	npoints = (int)(graph_size/2.0);
	if(data_loaded){
		int data_size;
		TRYCXX( VecGetSize(*data_Vec, &data_size) );
		T = (int)(data_size/(double)get_npoints());
	}

	compute_mapping();
	
	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);
}

Vec * eegimageplotter::get_data_Vec(){
	return data_Vec;
}

Vec * eegimageplotter::get_graph_Vec(){
	return graph_Vec;
}

bool eegimageplotter::get_data_loaded(){
	return data_loaded;
}

bool eegimageplotter::get_graph_loaded(){
	return graph_loaded;
}

std::string eegimageplotter::get_data_filename() const{
	return this->data_filename;
}

std::string eegimageplotter::get_graph_filename() const{
	return this->graph_filename;
}

int eegimageplotter::get_T() const{
	return this->T;
}

int eegimageplotter::get_npoints() const{
	return this->npoints;
}

void eegimageplotter::plot_image(const canvas& c) const{
	/* color graph */
	if(type==0){
		if(graph_loaded && data_loaded){
			plot_graph_data(c);
		}
	}

	/* interpolation */
	if(type==1){
		/* plot data */
		if(graph_loaded && data_loaded){
			plot_data(c);
		}

		/* plot graph */
		if(graph_loaded){
			plot_graph(c);
		}
	}
}

void eegimageplotter::plot_graph(const canvas& c) const{
	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();

	unsigned long x_size = this->width();
	unsigned long y_size = this->height();

	/* define radius */
	double radius;
	if(x_size < y_size){
		radius = x_size/(double)(this->radius_ratio);
	} else {
		radius = y_size/(double)(this->radius_ratio);
	}
	
	point center_point;
	double *coordinates_arr; /* [x_0, x_1, ... x_{n-1}, y_0, y_1, ... y_{n-1} ] */
	TRYCXX( VecGetArray(*graph_Vec, &coordinates_arr) );
	for(int i=0;i<npoints;i++){
		center_point.x() = x_begin + ax*(coordinates_arr[i]) + bx;
		center_point.y() = y_begin + (y_size - (ay*(coordinates_arr[npoints + i]) + by));
		
		draw_circle (c, center_point, radius, rgb_pixel(0,0,0));
	}
	TRYCXX( VecRestoreArray(*graph_Vec, &coordinates_arr) );

}

void eegimageplotter::plot_data(const canvas& c) const{
	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();

	unsigned long x_size = this->width();
	unsigned long y_size = this->height();

	/* get the values on vertices */
	int npoints = get_npoints();
	double *data_arr;
	TRYCXX( VecGetArray(*data_Vec, &data_arr) );

	/* maybe it is not necessary to plot all pixels */
	double pixel_size = PIXEL_SIZE;

	/* go through image and plot pixels */
	double *coordinates_arr; /* [x_0, x_1, ... x_{n-1}, y_0, y_1, ... y_{n-1} ] */
	TRYCXX( VecGetArray(*graph_Vec, &coordinates_arr) );
	
	rectangle rect_pixel;
	rgb_pixel color_pixel;
	double coeff;
	double x1, x2, y1, y2;
	double center_x, center_y;
	double myvalue;
	
	double value_min, value_max; //VALUE_MIN, VALUE_MAX
	value_min = VALUE_MIN;
	value_max = VALUE_MAX;
//	TRYCXX( VecMin(*data_Vec, NULL, &value_min) );
//	TRYCXX( VecMax(*data_Vec, NULL, &value_max) );

	point center_point;
	double maxdistance = (x_size*x_size + y_size*y_size);
		
	for(int y_coor=0; y_coor < y_size/pixel_size; y_coor++){
		for(int x_coor=0; x_coor < x_size/pixel_size; x_coor++){
			x1 = x_begin + x_coor*pixel_size;
			x2 = x_begin + (x_coor+1)*pixel_size;
			y1 = y_begin + y_coor*pixel_size;
			y2 = y_begin + (y_coor+1)*pixel_size;
			
			rect_pixel.set_left(x1);
			rect_pixel.set_top(y1);
			rect_pixel.set_right(x2);
			rect_pixel.set_bottom(y2);

			center_x = (x1+x2)/2.0;
			center_y = (y1+y2)/2.0;

			/* compute coeffs as distance from nodes */
			myvalue = 0;
			for(int r=0;r<npoints;r++){
				center_point.x() = x_begin + ax*(coordinates_arr[r]) + bx;
				center_point.y() = y_begin + (y_size - (ay*(coordinates_arr[npoints + r]) + by));

				coeff = maxdistance - ((center_x - center_point.x())*(center_x - center_point.x()) + (center_y - center_point.y())*(center_y - center_point.y()));
				coeff = sqrt(coeff)/sqrt(maxdistance);

				myvalue += coeff*data_arr[t*npoints + r];

			}

			/* my value to [0,1] */
			myvalue = (myvalue - value_min)/(value_max - value_min);
			if(myvalue > 1.0) myvalue = 1.0;
			if(myvalue < 0.0) myvalue = 0.0;

			color_pixel.red =   (int)(COLOR_MIN_R + myvalue*(COLOR_MAX_R - COLOR_MIN_R));
			color_pixel.green = (int)(COLOR_MIN_G + myvalue*(COLOR_MAX_G - COLOR_MIN_G));
			color_pixel.blue =  (int)(COLOR_MIN_B + myvalue*(COLOR_MAX_B - COLOR_MIN_B));
			
			fill_rect(c,rect_pixel,color_pixel);	
		}
	}

	//TODO: temp
	coutMaster << "t = " << t << " - computed" << std::endl;
	coutMaster << "data: ";
	for(int r=0;r<npoints;r++){
		 coutMaster << data_arr[t*npoints + r];
		 if(r < npoints-1) coutMaster << ", ";
	}
	coutMaster << std::endl;

	TRYCXX( VecRestoreArray(*data_Vec, &data_arr) );
	TRYCXX( VecRestoreArray(*graph_Vec, &coordinates_arr) );

}

void eegimageplotter::plot_graph_data(const canvas& c) const{

	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();

	unsigned long x_size = this->width();
	unsigned long y_size = this->height();

	/* define radius */
	double radius;
	if(x_size < y_size){
		radius = x_size/(double)(this->radius_ratio);
	} else {
		radius = y_size/(double)(this->radius_ratio);
	}
	
	point center_point;
	rgb_pixel color_pixel;
	double myvalue;

	double value_min, value_max; //VALUE_MIN, VALUE_MAX
	value_min = VALUE_MIN;
	value_max = VALUE_MAX;

	double *data_arr;
	TRYCXX( VecGetArray(*data_Vec, &data_arr) );

	double *coordinates_arr; /* [x_0, x_1, ... x_{n-1}, y_0, y_1, ... y_{n-1} ] */
	TRYCXX( VecGetArray(*graph_Vec, &coordinates_arr) );
	for(int r=0;r<npoints;r++){
		center_point.x() = x_begin + ax*(coordinates_arr[r]) + bx;
		center_point.y() = y_begin + (y_size - (ay*(coordinates_arr[npoints + r]) + by));
		
		myvalue = data_arr[t*npoints + r];
		myvalue = (myvalue - value_min)/(value_max - value_min);
		if(myvalue > 1.0) myvalue = 1.0;
		if(myvalue < 0.0) myvalue = 0.0;
		
		color_pixel.red =   (int)(COLOR_MIN_R + myvalue*(COLOR_MAX_R - COLOR_MIN_R));
		color_pixel.green = (int)(COLOR_MIN_G + myvalue*(COLOR_MAX_G - COLOR_MIN_G));
		color_pixel.blue =  (int)(COLOR_MIN_B + myvalue*(COLOR_MAX_B - COLOR_MIN_B));
		
		draw_solid_circle(c, center_point, radius, color_pixel);
	}

	TRYCXX( VecRestoreArray(*graph_Vec, &coordinates_arr) );
	TRYCXX( VecRestoreArray(*data_Vec, &data_arr) );

}

int eegimageplotter::set_t(int new_t){
	this->t = new_t;

	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);
}

void eegimageplotter::compute_mapping(){

	if(graph_loaded){
		double *coordinates_arr;

		TRYCXX( VecGetArray(*graph_Vec, &coordinates_arr) );
		
		unsigned long x_size = this->width();
		unsigned long y_size = this->height();

		/* find minimal coordinates in graph nodes */
		double coorx_min = std::numeric_limits<double>::max();
		double coorx_max = - std::numeric_limits<double>::max();
		double coory_min = std::numeric_limits<double>::max();
		double coory_max = - std::numeric_limits<double>::max();
		for(int i=0;i<npoints;i++){
			if(coordinates_arr[i] < coorx_min) coorx_min = coordinates_arr[i];
			if(coordinates_arr[i] > coorx_max) coorx_max = coordinates_arr[i];
			if(coordinates_arr[npoints + i] < coory_min) coory_min = coordinates_arr[npoints + i];
			if(coordinates_arr[npoints + i] > coory_max) coory_max = coordinates_arr[npoints + i];
		}

		/* define linear mapping between coordinates and plotting area */
		/* free space parameters */
		double px_min = DRAW_SPACE*x_size;
		double px_max = (1 - DRAW_SPACE)*x_size;
		double py_min = DRAW_SPACE*y_size;
		double py_max = (1 - DRAW_SPACE)*y_size;

		/* coefficients of mapping t to x */
		ax = (px_max - px_min)/(double)(coorx_max - coorx_min);
		bx = px_max - ax*coorx_max;

		/* coefficients of mapping value to y */
		ay = (py_max - py_min)/(double)(coory_max - coory_min);
		by = py_max - ay*coory_max;

		TRYCXX( VecRestoreArray(*graph_Vec, &coordinates_arr) );

	} else {
		ax = 0;
		bx = 0;
		ay = 0;
		by = 0;
	}
}

