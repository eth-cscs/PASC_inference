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

#define NUMBER_OF_LABELS 13

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

		int T;
		int npoints;
	public: 
		eegimageplotter(drawable_window& w);
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
    show_image_window();
    ~show_image_window();

};

//  ----------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

    // create our window
    show_image_window my_window;

    // wait until the user closes this window before we let the program 
    // terminate.
    my_window.wait_until_closed();

	/* finalize the library */
	Finalize<PetscVector>();
	
    return 0;
}


/* ---------------- implementation -------------- */

show_image_window::show_image_window() : /* All widgets take their parent window as an argument to their constructor. */
        mbar(*this),
        t_label(*this),
        myeegimageplotter(*this)
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
    t_label.set_pos(10,90);

    /* prepare position of labels of vector properties */
    labels_properties[0]->set_pos(10,160);
	for(int i=1; i < NUMBER_OF_LABELS; i++){
		labels_properties[i]->set_pos(labels_properties[i-1]->left(),labels_properties[i-1]->bottom()+20);
	}

	/* window title */
	set_title("Show PETSc EEG utility");
        
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

	/* prepare time scroll_bar */
	timescroll = new scroll_bar(*this, scroll_bar::HORIZONTAL);
	timescroll->set_scroll_handler( *this, &show_image_window::on_timescroll );
	timescroll->set_length(180);
	timescroll->set_pos(10,120);
	timescroll->disable();

	/* arrange the window */
	on_window_resized();

	fill_labels();
	show();
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

	set_label_properties(5, "DATA");
	set_label_properties(6, " loaded:    ", print_bool(myeegimageplotter.get_data_loaded()));
	if(myeegimageplotter.get_data_loaded()){
		Vec *data_Vec = myeegimageplotter.get_data_Vec();
		set_label_properties(7, " name:      ", cut_filename(myeegimageplotter.get_data_filename()));
		set_label_properties(8, " T:         ", myeegimageplotter.get_T());

		TRYCXX( VecGetSize(*data_Vec, &data_size) );
		set_label_properties(9, " size:      ", data_size);

		TRYCXX( VecNorm(*data_Vec, NORM_2, &data_norm2) );
		set_label_properties(10, " norm2:     ", data_norm2);

		TRYCXX( VecMax(*data_Vec, NULL, &data_max) );
		set_label_properties(11, " max:       ", data_max);

		TRYCXX( VecMin(*data_Vec, NULL, &data_min) );
		set_label_properties(12, " min:       ", data_min);
	}
	
	if(myeegimageplotter.get_graph_loaded() && myeegimageplotter.get_data_loaded()){
		std::ostringstream sout;
		sout << "t = " << t << " (" << myeegimageplotter.get_T() << ")";
		t_label.set_text(sout.str());


		if(t > myeegimageplotter.get_T()-1){
			t = myeegimageplotter.get_T()-1;
		}
		timescroll->set_max_slider_pos(myeegimageplotter.get_T()-1);
		timescroll->set_slider_pos(t);
		timescroll->set_jump_size(timescroll->max_slider_pos()/10.0);
		timescroll->enable();
	} else {
		t_label.set_text("");
		timescroll->disable();
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
	
	/* plot vector */
	if(data_loaded && graph_loaded){
		plot_image(c);
	}

}

eegimageplotter::eegimageplotter(drawable_window& w): 
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
//	unsigned long x_begin = this->left();
//	unsigned long y_begin = this->top();
	
//	unsigned long x_size = this->width();
//	unsigned long y_size = this->height();

	/* get the length of the signal */
//	int myvector_size;
//	TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );

	/* coefficients of mapping t to x */
//	double ax = (px_max - px_min)/(double)(myvector_size - 1);
//	double bx = px_max - ax*myvector_size;

	/* coefficients of mapping value to y */
//	double ay = (py_max - py_min)/(double)(myvector_max - myvector_min);
//	double by = py_max - ay*myvector_max;

	/* get some plotting variables */
//	rectangle rect_pixel;
//	int x_coor, y_coor;
//	double pixel_size_x = x_size/(double)image_width;
//	double pixel_size_y = y_size/(double)image_height;
//	double pixel_size;
//	if(pixel_size_x < pixel_size_y){
//		pixel_size = pixel_size_x;
//	} else {
//		pixel_size = pixel_size_y;
//	}

	
//	double *values;
//	TRYCXX( VecGetArray(*myvector_Vec,&values) );
//	for(int t=1;t<myvector_size;t++){
//		y_coor = t/(double)image_width;
//		x_coor = t - y_coor*image_width;

//		rect_pixel.set_left(x_begin + x_coor*pixel_size);
//		rect_pixel.set_top(y_begin + y_coor*pixel_size);
//		rect_pixel.set_right(x_begin + (x_coor+1)*pixel_size);
//		rect_pixel.set_bottom(y_begin + (y_coor+1)*pixel_size);

//		fill_rect(c,rect_pixel,rgb_pixel(values[t-1]*255,values[t-1]*255,values[t-1]*255));		
//	}
//	TRYCXX( VecRestoreArray(*myvector_Vec,&values) );
	
}

