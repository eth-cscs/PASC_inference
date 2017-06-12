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

#include "graphplotter.h"
#include "shortinfo_reader.h"

using namespace dlib;
using namespace pascinference;

class show_shortinfo_window : public drawable_window {
private:
    menu_bar mbar; /* gui: menubar */
	popup_menu submenu_xlabel; /* gui: submenu */
	popup_menu submenu_ylabel; /* gui: submenu */

	graphplotter mygraphplotter; /* gui: canvas for drawing */
	std::string xlabel;
	std::string ylabel;	
 
	int debug;
 
	Shortinfo_reader *myreader;
    void load_shortinfo( const std::string& filename );

	/* menu events */
    void on_menu_file_open ();
    void on_menu_file_quit ();
    void on_menu_help_about();
	
	/* general window events */
	void on_window_resized();
	
	class menu_item_xlabel : public menu_item_text	{
		private:
			std::string xlabel;
			show_shortinfo_window *mywindow;
		public:
			menu_item_xlabel(const std::string &name, const std::string &xlabel, show_shortinfo_window& mywindow, void (show_shortinfo_window::*on_click_handler)() ) : menu_item_text(name, mywindow, on_click_handler) {
				this->xlabel = xlabel;
				this->mywindow = &mywindow;
			}

			virtual void on_click () const {
				mywindow->set_xlabel(this->xlabel);
			}
	};

	class menu_item_ylabel : public menu_item_text	{
		private:
			std::string ylabel;
			show_shortinfo_window *mywindow;
		public:
			menu_item_ylabel(const std::string &name, const std::string &ylabel, show_shortinfo_window& mywindow, void (show_shortinfo_window::*on_click_handler)() ) : menu_item_text(name, mywindow, on_click_handler) {
				this->ylabel = ylabel;
				this->mywindow = &mywindow;
			}

			virtual void on_click () const {
				mywindow->set_ylabel(this->ylabel);
			}	
	};


	void set_xlabel(std::string xlabel);
	void set_ylabel(std::string ylabel);
	void empty_void(){}
	void set_title_labels();

public:
    show_shortinfo_window(int debug = 0);
    show_shortinfo_window(std::string filename, std::string title, std::string xlabel, std::string ylabel, int debug = 0);

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
		("ylabel", boost::program_options::value<std::string>(), "label in shortinfo to use as y [string]")
		("debug", boost::program_options::value<int>(), "debug mode [int]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* get values provided as console parameters */
	std::string filename, title, plot_epssqrvs, xlabel, ylabel;
	int debug;

	consoleArg.set_option_value("filename", &filename, "");
	consoleArg.set_option_value("title", &title, "");
	consoleArg.set_option_value("xlabel", &xlabel, "");
	consoleArg.set_option_value("ylabel", &ylabel, "");
	consoleArg.set_option_value("debug", &debug, 0);

	/* create our window */
	show_shortinfo_window my_window(filename, title, xlabel, ylabel, debug);

	my_window.wait_until_closed();

	/* finalize the library */
	Finalize<PetscVector>();
	
    return 0;
}


/* ---------------- implementation -------------- */

show_shortinfo_window::show_shortinfo_window(int debug) : /* All widgets take their parent window as an argument to their constructor. */
        mbar(*this),
        mygraphplotter(*this, debug)
{
	this->debug = debug;

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
	mbar.menu(1).add_submenu(menu_item_submenu("X axis",'x'), submenu_xlabel);
	mbar.menu(1).add_submenu(menu_item_submenu("Y axis",'y'), submenu_ylabel);

    /* Add the entries to the Help menu. */
    mbar.menu(2).add_menu_item(menu_item_text("About",  *this, &show_shortinfo_window::on_menu_help_about,   'A'));

	/* arrange the window */
	on_window_resized();

	show();
} 

show_shortinfo_window::show_shortinfo_window(std::string filename, std::string title, std::string xlabel, std::string ylabel, int debug) : show_shortinfo_window(debug){
	if(title != ""){
		set_title(title);		
	}

	if(xlabel != ""){
		set_xlabel(xlabel);
	}

	if(ylabel != ""){
		set_ylabel(ylabel);
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

	/* clear menus */
	submenu_xlabel.clear();
	submenu_ylabel.clear();

	/* get headers */
	std::vector<std::string> headers = myreader->get_headers();
	for(int i=0; i < headers.size();i++){
		submenu_xlabel.add_menu_item(menu_item_xlabel(headers[i], headers[i], *this, &show_shortinfo_window::empty_void) );
		submenu_ylabel.add_menu_item(menu_item_ylabel(headers[i], headers[i], *this, &show_shortinfo_window::empty_void) );
	}

}

void show_shortinfo_window::set_title_labels() {

	if(mygraphplotter.get_xvalues_loaded() && mygraphplotter.get_yvalues_loaded()){
		std::ostringstream sout;
		sout << mygraphplotter.get_xlabel() << " x " << mygraphplotter.get_ylabel();

		set_title(sout.str());
	}
	
}

void show_shortinfo_window::set_xlabel(std::string xlabel){
	this->xlabel = xlabel;

	//TODO: categorical (string) values? 
	std::vector<double> xvalues = myreader->get_values_double(this->xlabel);

	mygraphplotter.set_xlabel(xlabel);
	mygraphplotter.set_xvalues(xvalues);

	set_title_labels();
}

void show_shortinfo_window::set_ylabel(std::string ylabel){
	this->ylabel = ylabel;

	//TODO: categorical (string) values? 
	std::vector<double> yvalues = myreader->get_values_double(this->ylabel);
	
	mygraphplotter.set_ylabel(ylabel);
	mygraphplotter.set_yvalues(yvalues);

	set_title_labels();
}

