#include "pascinference.h"

#ifndef USE_PETSC
 #error 'this util is for PETSC'
#endif

#ifndef USE_DLIB
 #error 'this util is for DLIB'
#endif

#include <dlib/gui_widgets.h>
#include <sstream>
#include <string>

//using namespace std;
using namespace dlib;
using namespace pascinference;


class show_vec_window : public drawable_window 
{
public:
    show_vec_window();
    ~show_vec_window();

private:
	Vec *myvector_Vec;
	GeneralVector<PetscVector> *myvector;


    menu_bar mbar; /* gui: menubar */

    void load_vector( const std::string& file_name );

    /* Event handlers */
    void on_menu_file_open ();
    void on_menu_file_quit ();
//   void on_menu_file_save ();
//   void on_menu_file_save_as ();
    void on_menu_help_about ();
//    void on_menu_help_help ();
    void on_window_resized ();

    // Member data
/*    const rgb_pixel color_non_evidence;
    const rgb_pixel color_default_bg;
    const rgb_pixel color_evidence;
    const rgb_pixel color_error;
    const rgb_pixel color_gray;
    bool graph_modified_since_last_recalc;

    button btn_calculate;
    check_box sel_node_is_evidence;
    directed_graph_drawer<directed_graph_type> graph_drawer;
    label sel_node_index;
    label sel_node_num_values_label; 
    label sel_node_text_label;
    label sel_node_evidence_label;

    named_rectangle selected_node_rect;
    tabbed_display tables;
    text_field sel_node_num_values;
    text_field sel_node_text;
    text_field sel_node_evidence;
    text_grid cpt_grid;
    text_grid ppt_grid;
    unsigned long selected_node_index;
    bool node_is_selected;
    widget_group cpt_group;
    widget_group ppt_group;

    scoped_ptr<bayesian_network_join_tree> solution;
    join_tree_type join_tree;
    // The std_vector_c is an object identical to the std::vector except that it checks
    // all its preconditions and throws a dlib::fatal_error if they are violated.
    std_vector_c<assignment> cpt_grid_assignments;
    std::string graph_file_name;
*/ 
};




int main( int argc, char *argv[] )
{
	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 


    /* create our window */
    show_vec_window my_window;

    /* tell our window to put itself on the screen */
    my_window.show();

    /* wait until the user closes this window before we let the program */
    my_window.wait_until_closed();		


	/* finalize the library */
	Finalize<PetscVector>();

    return 0;
}



show_vec_window::show_vec_window() : mbar(*this) {
	
	/* window title */
	set_title("Show PETSc Vector utility");

	/* create menu bar */
    mbar.set_number_of_menus(2);
    mbar.set_menu_name(0,"File",'F');
    mbar.set_menu_name(1,"Help",'H');

    /* add the entries to the File menu. */
    mbar.menu(0).add_menu_item(menu_item_text("Open",   *this, &show_vec_window::on_menu_file_open,    'O'));
    mbar.menu(0).add_menu_item(menu_item_separator());
    mbar.menu(0).add_menu_item(menu_item_text("Quit",   *this, &show_vec_window::on_menu_file_quit,    'Q'));

    /* Add the entries to the Help menu. */
    mbar.menu(1).add_menu_item(menu_item_text("About",  *this, &show_vec_window::on_menu_help_about,   'A'));

	/* set the size of window */
    set_size(600,350);
	

	/* prepare vectors */
	myvector_Vec = new Vec();
	TRYCXX( VecCreate(PETSC_COMM_WORLD,myvector_Vec) );
	myvector = new GeneralVector<PetscVector>(*myvector_Vec);

	/* arrange the window */
	on_window_resized();
}

show_vec_window::~show_vec_window() {
    close_window();
}

void show_vec_window::load_vector( const std::string& file_name ) {
	myvector->load_local(file_name);

	/* print properties of vectors */
	coutMaster << std::setprecision(17);	
	coutMaster << std::endl;
	coutMaster << " size      = " << std::setw(30) << myvector->size() << std::endl;
	coutMaster << " norm      = " << std::setw(30) << norm(*myvector) << std::endl;
	coutMaster << " sum       = " << std::setw(30) << sum(*myvector) << std::endl;
	coutMaster << " max       = " << std::setw(30) << max(*myvector) << std::endl;
	coutMaster << " min       = " << std::setw(30) << min(*myvector) << std::endl;

}

void show_vec_window::on_menu_file_open () {
    /* display a file chooser window and when the user choses a file call the on_open_file_selected() function */
    open_existing_file_box(*this, &show_vec_window::load_vector);	
}


void show_vec_window::on_menu_file_quit () {
	close_window();	
}

void show_vec_window::on_menu_help_about () {
	message_box("About","This application is for PETSc vector visualisation\n");
}

void show_vec_window::on_window_resized () {
	
//	drawable_window::on_window_resized();
	
	show();
}
