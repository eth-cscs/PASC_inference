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

#define NUMBER_OF_LABELS 9

using namespace dlib;
using namespace pascinference;

class imageplotter : public drawable {
	private:
		Vec *myvector_Vec;
		bool myvector_loaded;

		void draw (const canvas& c) const;
		void plot_vector(const canvas& c) const;

	public: 
		imageplotter(drawable_window& w);
		~imageplotter();	
		
		void set_plotting_area(rectangle area);
		void load_image( const std::string& file_name );
		bool get_myvector_loaded();
		Vec *get_myvector_Vec();
};

class show_image_window : public drawable_window {
private:
    menu_bar mbar; /* gui: menubar */
	imageplotter myimageplotter; /* gui: canvas for drawing */
 
    label **labels_myvector_properties;

    void load_image( const std::string& file_name );

	/* menu events */
    void on_menu_file_open ();
    void on_menu_file_quit ();
    void on_menu_help_about();
	
	/* general window events */
	void on_window_resized();
	
	template<class ValueType>
	void set_label_myvector_properties(int label_idx, const std::string &text, ValueType value);
	std::string cut_filename(const std::string &input);

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
        myimageplotter(*this)
{

	/* set the size of window */
    set_size(600,350);
		
	/* allocate labels */
	labels_myvector_properties = new label*[NUMBER_OF_LABELS];
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		labels_myvector_properties[i] = new label(*this);
	}

    /* prepare position of labels of vector properties */
    labels_myvector_properties[0]->set_pos(10,30);
	for(int i=1; i < 9; i++){
		labels_myvector_properties[i]->set_pos(labels_myvector_properties[i-1]->left(),labels_myvector_properties[i-1]->bottom()+20);
	}

	/* window title */
	set_title("Show PETSc Vector utility");
        
	/* create menu bar */
    mbar.set_number_of_menus(2);
    mbar.set_menu_name(0,"File",'F');
    mbar.set_menu_name(1,"Help",'H');

    /* add the entries to the File menu. */
    mbar.menu(0).add_menu_item(menu_item_text("Open",   *this, &show_image_window::on_menu_file_open,    'O'));
    mbar.menu(0).add_menu_item(menu_item_separator());
    mbar.menu(0).add_menu_item(menu_item_text("Quit",   *this, &show_image_window::on_menu_file_quit,    'Q'));

    /* Add the entries to the Help menu. */
    mbar.menu(1).add_menu_item(menu_item_text("About",  *this, &show_image_window::on_menu_help_about,   'A'));

	/* arrange the window */
	on_window_resized();

	show();
} 

show_image_window::~show_image_window(){

	/* destroy labels */
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		free(labels_myvector_properties[i]);
	}
	free(labels_myvector_properties);

	/* close window */
	close_window();
}

void show_image_window::on_menu_help_about(){
     message_box("About","This application is for PETSc vector visualisation\n");
}

void show_image_window::on_menu_file_open(){
    /* display a file chooser window and when the user choses a file */
    open_existing_file_box(*this, &show_image_window::load_image);
}

void show_image_window::on_menu_file_quit(){
	close_window();
}

void show_image_window::on_window_resized() {

	/* set new plotting area */
	unsigned long width,height;
    get_size(width,height);
	myimageplotter.set_plotting_area(rectangle(200,mbar.bottom(),width,height));

	drawable_window::on_window_resized();

}

void show_image_window::load_image( const std::string& file_name ) {
	myimageplotter.load_image(file_name);

	Vec *myvector_Vec = myimageplotter.get_myvector_Vec();

	/* compute basic properties of loaded vector */
	int myvector_size;
	double myvector_norm2;
	double myvector_norm1;
	double myvector_normInf;
	double myvector_sum;
	double myvector_max;
	double myvector_min;
	TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );
	TRYCXX( VecNorm(*myvector_Vec, NORM_2, &myvector_norm2) );
	TRYCXX( VecNorm(*myvector_Vec, NORM_1, &myvector_norm1) );
	TRYCXX( VecNorm(*myvector_Vec, NORM_INFINITY, &myvector_normInf) );
	TRYCXX( VecSum(*myvector_Vec, &myvector_sum) );
	TRYCXX( VecMax(*myvector_Vec, NULL, &myvector_max) );
	TRYCXX( VecMin(*myvector_Vec, NULL, &myvector_min) );

	/* print properties of vectors to labels */
	set_label_myvector_properties(0, "name:    ", cut_filename(file_name));
	set_label_myvector_properties(1, "size:    ", myvector_size);
	set_label_myvector_properties(2, "norm2:   ", myvector_norm2);
	set_label_myvector_properties(3, "norm1:   ", myvector_norm1);
	set_label_myvector_properties(4, "normInf: ", myvector_normInf);
	set_label_myvector_properties(5, "sum:     ", myvector_sum);
	set_label_myvector_properties(6, "max:     ", myvector_max);
	set_label_myvector_properties(7, "min:     ", myvector_min);
	set_label_myvector_properties(8, "mean:    ", myvector_sum/(double)myvector_size);

}

template<class ValueType>
void show_image_window::set_label_myvector_properties(int label_idx, const std::string &text, ValueType value){
	std::ostringstream sout;
	sout << std::setprecision(5);
    sout << text << std::setw(20) << value;
    labels_myvector_properties[label_idx]->set_text(sout.str());	
}

std::string show_image_window::cut_filename(const std::string &input){
	std::ostringstream sout;
	boost::filesystem::path p(input);
	sout << p.filename();

	return sout.str();
}

/* ------------------------- graph plotter -------------- */

void imageplotter::draw(const canvas& c) const {
	/* draw background */
	fill_rect(c,rect,rgb_pixel(255,255,255));
	
	/* plot vector */
	if(myvector_loaded){
		plot_vector(c);
	}

}

imageplotter::imageplotter(drawable_window& w): 
			drawable(w)
{
	/* create empty vector */
	myvector_Vec = new Vec();
	TRYCXX( VecCreate(PETSC_COMM_WORLD,myvector_Vec) );
	myvector_loaded = false;
	
	enable_events();
}

imageplotter::~imageplotter(){
	free(myvector_Vec);

	disable_events();
	parent.invalidate_rectangle(rect);
}

void imageplotter::set_plotting_area(rectangle area){
	rect = area;
}

void imageplotter::load_image( const std::string& file_name ) {
	if(myvector_loaded){
		/* destroy existing vector */
		TRYCXX( VecDestroy(myvector_Vec) );
		TRYCXX( VecCreate(PETSC_COMM_WORLD,myvector_Vec) );
	}

	/* load the data from file */
	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_SELF, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_SELF ,file_name.c_str(), FILE_MODE_READ, &mviewer) );

	/* load vector from viewer */
	TRYCXX( VecLoad(*myvector_Vec, mviewer) );

	/* destroy the viewer */
	TRYCXX( PetscViewerDestroy(&mviewer) );

	myvector_loaded = true;

	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);

}


Vec * imageplotter::get_myvector_Vec(){
	return myvector_Vec;
}

bool imageplotter::get_myvector_loaded(){
	return myvector_loaded;
}

void imageplotter::plot_vector(const canvas& c) const{
	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();
	
	unsigned long x_size = this->width();
	unsigned long y_size = this->height();

	/* get the length of the signal */
	int    myvector_size;
	double myvector_max;
	double myvector_min;
	TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );
	TRYCXX( VecMax(*myvector_Vec, NULL, &myvector_max) );
	TRYCXX( VecMin(*myvector_Vec, NULL, &myvector_min) );

	/* free space parameters */
	double px_min = 0.01*x_size;
	double px_max = 0.99*x_size;
	double py_min = 0.1*y_size;
	double py_max = 0.9*y_size;

	/* coefficients of mapping t to x */
	double ax = (px_max - px_min)/(double)(myvector_size - 1);
	double bx = px_max - ax*myvector_size;

	/* coefficients of mapping value to y */
	double ay = (py_max - py_min)/(double)(myvector_max - myvector_min);
	double by = py_max - ay*myvector_max;

	/* I will use these points */
	point mypoint1;
	point mypoint2;

	double *values;
	TRYCXX( VecGetArray(*myvector_Vec,&values) );
	for(int t=1;t<myvector_size;t++){
		mypoint1(0) = x_begin + ax*(t-1) + bx;
		mypoint1(1) = y_begin + (y_size - (ay*values[t-1] + by));

		mypoint2(0) = x_begin + ax*t + bx;
		mypoint2(1) = y_begin + (y_size - (ay*values[t] + by));
		
		draw_line(c,mypoint1,mypoint2, rgb_pixel(0,0,255));

	}
	TRYCXX( VecRestoreArray(*myvector_Vec,&values) );
	
	/* plot a line */
	
}
