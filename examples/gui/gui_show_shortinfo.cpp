#include "pascinference.h"

#ifndef USE_PETSC
 #error 'this util is for PETSC'
#endif

#ifndef USE_DLIB
 #error 'this util is for DLIB'
#endif

#ifndef USE_BOOST
 #error 'this util is for BOOST'
#endif

#define ENABLE_ASSERTS

#include <iostream>
#include <fstream>
#include <string>

#include <dlib/gui_widgets.h>
#include <sstream>
#include <string>

#include "shortinfo_reader.h"

using namespace dlib;
using namespace pascinference;

class graphplotter : public drawable {
	private:
		std::vector<double> x_values;
		std::vector<double> y_values;

		bool values_loaded;

		void draw (const canvas& c) const;
		void plot_graph(const canvas& c) const;

	public: 
		graphplotter(drawable_window& w);
		~graphplotter();	
		
		void set_plotting_area(rectangle area);

		void set_values(std::vector<double> x_values, std::vector<double> y_values);
		bool get_values_loaded();
		std::vector<double> get_x_values();
		std::vector<double> get_y_values();
};

class show_shortinfo_window : public drawable_window {
private:
    menu_bar mbar; /* gui: menubar */
	popup_menu submenu_epssqr_vs; /* gui: submenu */

	graphplotter mygraphplotter; /* gui: canvas for drawing */
	std::string xlabel;
	std::string ylabel;	
 
	Shortinfo_reader *myreader;
    void load_shortinfo( const std::string& filename );

	/* menu events */
    void on_menu_file_open ();
    void on_menu_file_quit ();
    void on_menu_help_about();
	
	/* general window events */
	void on_window_resized();
	
	class menu_item_whattoplot : public menu_item_text	{
		private:
			std::string xlabel;
			std::string ylabel;
			show_shortinfo_window *mywindow;
			
		public:
//			template <typename T>
//			menu_item_whattoplot(const std::string &name, const std::string &xlabel, const std::string &ylabel, T& object ) : menu_item_text(name, object ) {
			menu_item_whattoplot(const std::string &name, const std::string &xlabel, const std::string &ylabel, show_shortinfo_window& mywindow, void (show_shortinfo_window::*on_click_handler)() ) : menu_item_text(name, mywindow, on_click_handler) {
				this->xlabel = xlabel;
				this->ylabel = ylabel;
				this->mywindow = &mywindow;
			}

			virtual void on_click () const {
				mywindow->set_value_names(this->xlabel, this->ylabel);
			}	
	};

	void set_value_names(std::string xlabel, std::string ylabel);
	void empty_void(){}

public:
    show_shortinfo_window();
    show_shortinfo_window(std::string filename, std::string title, std::string xlabel, std::string ylabel);

    ~show_shortinfo_window();

};

//  ----------------------------------------------------------------------------



int main( int argc, char *argv[] ) {
	/* add local program options */
	boost::program_options::options_description opt_problem("GUI_SHOW_SHORTINFO", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename", boost::program_options::value<std::string>(), "shortinfo to load [string]")
		("title", boost::program_options::value<std::string>(), "title of window [string]")
		("xlabel", boost::program_options::value<std::string>(), "label in shortinfo to use as x [string]")
		("ylabel", boost::program_options::value<std::string>(), "label in shortinfo to use as y [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* get values provided as console parameters */
	std::string filename, title, plot_epssqrvs, xlabel, ylabel;

	consoleArg.set_option_value("filename", &filename, "");
	consoleArg.set_option_value("title", &title, "");
	consoleArg.set_option_value("xlabel", &xlabel, "");
	consoleArg.set_option_value("ylabel", &ylabel, "");

	/* create our window */
	show_shortinfo_window my_window(filename, title, xlabel, ylabel);

	my_window.wait_until_closed();

	/* finalize the library */
	Finalize<PetscVector>();
	
    return 0;
}


/* ---------------- implementation -------------- */

show_shortinfo_window::show_shortinfo_window() : /* All widgets take their parent window as an argument to their constructor. */
        mbar(*this),
        mygraphplotter(*this)
{

	/* set the size of window */
    set_size(600,350);
		
	/* window title */
	set_title("Show PETSc Shortinfo utility");
        
	/* create menu bar */
    mbar.set_number_of_menus(3);
    mbar.set_menu_name(0,"File",'F');
	mbar.set_menu_name(1,"Plot",'W');
    mbar.set_menu_name(2,"Help",'H');

    /* add the entries to the File menu. */
    mbar.menu(0).add_menu_item(menu_item_text("Open",   *this, &show_shortinfo_window::on_menu_file_open,    'O'));
    mbar.menu(0).add_menu_item(menu_item_separator());
    mbar.menu(0).add_menu_item(menu_item_text("Quit",   *this, &show_shortinfo_window::on_menu_file_quit,    'Q'));

	/* add plotting menus */
	mbar.menu(1).add_submenu(menu_item_submenu("plot epssqr vs",'e'), submenu_epssqr_vs);

    /* Add the entries to the Help menu. */
    mbar.menu(2).add_menu_item(menu_item_text("About",  *this, &show_shortinfo_window::on_menu_help_about,   'A'));

	/* arrange the window */
	on_window_resized();

	show();
} 

show_shortinfo_window::show_shortinfo_window(std::string filename, std::string title, std::string xlabel, std::string ylabel) : show_shortinfo_window(){
	if(title != ""){
		set_title(title);		
	}

	if(xlabel != "" && ylabel != ""){
		set_value_names(xlabel, ylabel);
	}

	if(filename != ""){
		load_shortinfo( filename );
	}
}

show_shortinfo_window::~show_shortinfo_window(){

	/* close window */
	close_window();
}

void show_shortinfo_window::on_menu_help_about(){
     message_box("About","This application is for PETSc shortinfo visualisation\n");
}

void show_shortinfo_window::on_menu_file_open(){
    /* display a file chooser window and when the user choses a file */
    open_existing_file_box(*this, &show_shortinfo_window::load_shortinfo);
}

void show_shortinfo_window::on_menu_file_quit(){
	close_window();
}

void show_shortinfo_window::on_window_resized() {

	/* set new plotting area */
	unsigned long width,height;
    get_size(width,height);
	mygraphplotter.set_plotting_area(rectangle(0,mbar.bottom(),width,height));

	drawable_window::on_window_resized();

}

void show_shortinfo_window::load_shortinfo( const std::string& filename ) {

	/* prepare shortinfo reader */
	myreader = new Shortinfo_reader(filename);
	myreader->process();

	/* clear menu */
	submenu_epssqr_vs.clear();

	/* get headers */
	std::vector<std::string> headers = myreader->get_headers();
	for(int i=0; i < headers.size();i++){
		submenu_epssqr_vs.add_menu_item(menu_item_whattoplot(headers[i], "epssqr", headers[i], *this, &show_shortinfo_window::empty_void) );
	}

}


void show_shortinfo_window::set_value_names(std::string xlabel, std::string ylabel){
	this->xlabel = xlabel;
	this->ylabel = ylabel;
	std::cout << "xvalue: " << xlabel << std::endl;
	std::cout << "yvalue: " << ylabel << std::endl;
	
//	std::vector<double> myvalues = myreader.get_values_double("epssqr");
}


/* ------------------------- graph plotter -------------- */

void graphplotter::draw(const canvas& c) const {
	/* draw background */
	fill_rect(c,rect,rgb_pixel(255,255,255));
	
	/* plot vector */
	if(values_loaded){
		plot_graph(c);
	}

}

graphplotter::graphplotter(drawable_window& w): 
			drawable(w)
{
	values_loaded = false;
	enable_events();
}

graphplotter::~graphplotter(){
	disable_events();
	parent.invalidate_rectangle(rect);
}

void graphplotter::set_plotting_area(rectangle area){
	rect = area;
}

void graphplotter::set_values(std::vector<double> x_values, std::vector<double> y_values){
	this->x_values = x_values;
	this->y_values = y_values;
	values_loaded = true;

	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);
}

bool graphplotter::get_values_loaded(){
	return values_loaded;
}

std::vector<double> graphplotter::get_x_values(){
	return x_values;
}

std::vector<double> graphplotter::get_y_values(){
	return y_values;
}

void graphplotter::plot_graph(const canvas& c) const{

	
}
