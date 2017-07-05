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
#define DLIB_JPEG_SUPPORT

#include <dlib/gui_widgets.h>
#include <sstream>
#include <string>
#include <dlib/image_io.h>
//#include <dlib/image_transforms.h>

#define NUMBER_OF_LABELS 9

using namespace dlib;
using namespace pascinference;

std::string cut_filename(const std::string &input);
std::string cut_extension(const std::string &input);

class movieplotter : public drawable {
	private:
		int type;
		Vec *myvector_Vec;
		bool myvector_loaded;

		void draw (const canvas& c) const;
		void plot_image(const canvas& c) const;

		int movie_width;
		int movie_height;
		int movie_T;
		double movie_scale;
		int movie_xdim;
		std::string filename;

	public: 
		movieplotter(drawable_window& w, int type);
		~movieplotter();	
		
		void set_plotting_area(rectangle area);
		void load_movie( const std::string& file_name, int new_xdim );
		void save_image(const std::string& file_name);
		
		bool get_myvector_loaded();
		Vec *get_myvector_Vec();
		
		void set_size(int new_width, int new_xdim, int new_T);
		
		int get_width() const;
		int get_height() const;
		int get_T() const;
		double get_scale() const;
		int get_xdim() const;

		std::string get_filename() const;
		
		void recompute_scale();
};

class show_image_window : public drawable_window {
private:
    menu_bar mbar; /* gui: menubar */
	movieplotter mymovieplotter; /* gui: canvas for drawing */
 
    label **labels_myvector_properties;
	text_field select_width;
	label label_width;
	text_field select_height;
	label label_height;
	text_field select_T;
	label label_T;
	text_field select_scale;
	label label_scale;
	text_field select_xdim;
	label label_xdim;

	button button_apply_size;

    void load_movie( const std::string& file_name );
    void save_image( const std::string& file_name );
    
	/* menu events */
    void on_menu_file_open ();
	void on_menu_file_save();
    void on_menu_file_quit ();
    void on_menu_help_about();
	
	/* general window events */
	void on_window_resized();
	
	template<class ValueType>
	void set_label_myvector_properties(int label_idx, const std::string &text, ValueType value);
	void on_button_apply_size();

	void fill_labels();

public:
    show_image_window(int type);
    show_image_window(std::string filename, int type, std::string title, int width = 0, int T = 1, int xdim = 0);
    ~show_image_window();

};

//  ----------------------------------------------------------------------------

int main( int argc, char *argv[] ) {
	/* add local program options */
	boost::program_options::options_description opt_problem("GUI_SHOW_IMAGE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename", boost::program_options::value<std::string>(), "image to load [string]")
		("width", boost::program_options::value<int>(), "width of movie [int]")
		("T", boost::program_options::value<int>(), "number of frames [int]")
		("xdim", boost::program_options::value<int>(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("type", boost::program_options::value< int >(), "type of output vector [0=Rn, 1=nR]")
		("title", boost::program_options::value<std::string>(), "title of window [string]")	;
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* get values provided as console parameters */
	std::string filename, title;
	int width, xdim, type, T;

	consoleArg.set_option_value("filename", &filename, "");
	consoleArg.set_option_value("title", &title, "");
	consoleArg.set_option_value("width", &width, 0);
	consoleArg.set_option_value("T", &T, 1);
	consoleArg.set_option_value("xdim", &xdim, 1);
	consoleArg.set_option_value("type", &type, 1);

    /* create our window */
    show_image_window my_window(filename,type,title,width,T,xdim);

    /* wait until the user closes this window before we let the program terminate */
    my_window.wait_until_closed();

	/* finalize the library */
	Finalize<PetscVector>();
	
    return 0;
}


/* ---------------- implementation -------------- */

show_image_window::show_image_window(int type) : /* All widgets take their parent window as an argument to their constructor. */
        mbar(*this),
        select_width(*this),
		label_width(*this),
        select_height(*this),
		label_height(*this),
        select_T(*this),
		label_T(*this),
        select_scale(*this),
		label_scale(*this),
        select_xdim(*this),
		label_xdim(*this),
        button_apply_size(*this),
        mymovieplotter(*this, type)
{

	/* set the size of window */
    set_size(600,350);
		
	/* allocate labels */
	labels_myvector_properties = new label*[NUMBER_OF_LABELS];
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		labels_myvector_properties[i] = new label(*this);
	}

    /* prepare position of labels of vector properties */
    labels_myvector_properties[0]->set_pos(10,220);
	for(int i=1; i < 9; i++){
		labels_myvector_properties[i]->set_pos(labels_myvector_properties[i-1]->left(),labels_myvector_properties[i-1]->bottom()+20);
	}

	/* window title */
	set_title("Show PETSc movie utility");
      
	/* create menu bar */
    mbar.set_number_of_menus(2);
    mbar.set_menu_name(0,"File",'F');
    mbar.set_menu_name(1,"Help",'H');

    /* add the entries to the File menu. */
    mbar.menu(0).add_menu_item(menu_item_text("Open",   *this, &show_image_window::on_menu_file_open,    'O'));
    mbar.menu(0).add_menu_item(menu_item_text("Save",   *this, &show_image_window::on_menu_file_save,    'S'));
    mbar.menu(0).add_menu_item(menu_item_separator());
    mbar.menu(0).add_menu_item(menu_item_text("Quit",   *this, &show_image_window::on_menu_file_quit,    'Q'));

    /* Add the entries to the Help menu. */
    mbar.menu(1).add_menu_item(menu_item_text("About",  *this, &show_image_window::on_menu_help_about,   'A'));

	/* register input for size */
	select_width.set_pos(120,30);
	select_width.set_width(70);
	select_width.set_text("0");
	label_width.set_pos(10,35);
	label_width.set_text("width  :");
	
	select_height.set_pos(120,60);
	select_height.set_width(70);
	select_height.set_text("0");
	label_height.set_pos(10,65);
	label_height.set_text("height :");
	select_height.disable();

	select_T.set_pos(120,90);
	select_T.set_width(70);
	select_T.set_text("1");
	label_T.set_pos(10,95);
	label_T.set_text("T      :");

	select_scale.set_pos(120,120);
	select_scale.set_width(70);
	select_scale.set_text("1.0");
	label_scale.set_pos(10,125);
	label_scale.set_text("scale  :");
	select_scale.disable();

	select_xdim.set_pos(120,150);
	select_xdim.set_width(70);
	select_xdim.set_text("1");
	label_xdim.set_pos(10,155);
	label_xdim.set_text("xdim   :");

	/* register button for changing the size */
	button_apply_size.set_click_handler(*this, &show_image_window::on_button_apply_size);
	button_apply_size.set_name("apply size");
	button_apply_size.set_pos(120,180);

	/* arrange the window */
	on_window_resized();

	show();
} 

show_image_window::show_image_window(std::string filename, int type, std::string title, int width, int T, int xdim) : show_image_window(type) {

	if(title != ""){
		set_title(title);		
	}

	if(filename != ""){
		load_movie( filename );
	}

	select_xdim.set_text(std::to_string(xdim));
	select_width.set_text(std::to_string(width));
	select_T.set_text(std::to_string(T));
	
	on_button_apply_size();
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
     message_box("About","This application is for PETSc movie visualisation\n");
}

void show_image_window::on_menu_file_open(){
    /* display a file chooser window and when the user choses a file */
    open_existing_file_box(*this, &show_image_window::load_movie);
}

void show_image_window::on_menu_file_save(){
    /* display a file chooser window and when the user choses a file */
    save_file_box(*this, &show_image_window::save_image);
}

void show_image_window::on_menu_file_quit(){
	close_window();
}

void show_image_window::on_window_resized() {

	/* set new plotting area */
	unsigned long width,height;
    get_size(width,height);
	mymovieplotter.set_plotting_area(rectangle(200,mbar.bottom(),width,height));
	mymovieplotter.recompute_scale();
	select_scale.set_text(std::to_string(mymovieplotter.get_scale()));

	drawable_window::on_window_resized();

	
}

void show_image_window::load_movie( const std::string& file_name ) {
	mymovieplotter.load_movie(file_name, std::stoi(trim(select_xdim.text())) );

	if(mymovieplotter.get_myvector_loaded()){
		/* get recomputed height and set value */
		select_width.set_text(std::to_string(mymovieplotter.get_width()));
		select_height.set_text(std::to_string(mymovieplotter.get_height()));

		fill_labels();
	}

	mymovieplotter.recompute_scale();
	select_scale.set_text(std::to_string(mymovieplotter.get_scale()));

}

void show_image_window::save_image( const std::string& file_name ) {
	mymovieplotter.save_image(file_name);
}

void show_image_window::fill_labels() {

	Vec *myvector_Vec = mymovieplotter.get_myvector_Vec();

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
	set_label_myvector_properties(0, "name:    ", cut_filename(mymovieplotter.get_filename()));
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

void show_image_window::on_button_apply_size(){
	int new_width = std::stoi(trim(select_width.text()));
	int new_xdim = std::stoi(trim(select_xdim.text()));
	int new_T = std::stoi(trim(select_T.text()));
	
	mymovieplotter.set_size(new_width, new_xdim, new_T);
	mymovieplotter.recompute_scale();

	/* get recomputed height and set value */
	select_width.set_text(std::to_string(mymovieplotter.get_width()));
	select_height.set_text(std::to_string(mymovieplotter.get_height()));
	select_T.set_text(std::to_string(mymovieplotter.get_T()));
	select_scale.set_text(std::to_string(mymovieplotter.get_scale()));

	if(mymovieplotter.get_myvector_loaded()){
		fill_labels();
	}

}

/* ------------------------- graph plotter -------------- */

void movieplotter::draw(const canvas& c) const {
	/* draw background */
	fill_rect(c,rect,rgb_pixel(255,255,255));
	
	/* plot vector */
	if(myvector_loaded){
		plot_image(c);
	}

}

movieplotter::movieplotter(drawable_window& w, int type): 
			drawable(w)
{
	/* create empty vector */
	myvector_Vec = new Vec();
	TRYCXX( VecCreate(PETSC_COMM_WORLD,myvector_Vec) );
	myvector_loaded = false;
	
	this->type = type;
	movie_width = 0;
	movie_height = 0;
	movie_scale = 1.0;
	movie_xdim = 1;
	
	enable_events();
}

movieplotter::~movieplotter(){
	free(myvector_Vec);

	disable_events();
	parent.invalidate_rectangle(rect);
}

void movieplotter::set_plotting_area(rectangle area){
	rect = area;
}

void movieplotter::load_movie( const std::string& file_name, int new_xdim ) {
	this->filename = file_name;
	this->movie_xdim = new_xdim;

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

	/* compute some implicit width */
	int myvector_size;
	TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );
	int x_size = rect.width();
	if(myvector_size < x_size){ 
		movie_width = myvector_size;
	} else {
		movie_width = x_size;
	}
	movie_height = (int)(myvector_size/(double)movie_width);
	movie_scale = 1.0;

	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);

}

void movieplotter::save_image( const std::string& file_name ) {
	coutMaster << "saving image : " << file_name << std::endl;
	
	/* get filename extension */
	std::string extension = cut_extension(file_name);
	
	array2d<rgb_pixel> image_dlib(movie_height,movie_width);

	double *x_arr;
	TRYCXX( VecGetArray(*myvector_Vec, &x_arr) );
	for(int i=0;i<movie_height;i++){
		for(int j=0;j<movie_width;j++){
			
			/* greyscale */
			if(movie_xdim==1){
				image_dlib[i][j].red = x_arr[i*movie_width + j]*255;
				image_dlib[i][j].green = x_arr[i*movie_width + j]*255;
				image_dlib[i][j].blue = x_arr[i*movie_width + j]*255;
			}
			
			/* rgb - Rn */
			if(type==0){
				if(movie_xdim==3){
					image_dlib[i][j].red   = x_arr[i*movie_width + j*movie_xdim + 0]*255;
					image_dlib[i][j].green = x_arr[i*movie_width + j*movie_xdim + 1]*255;
					image_dlib[i][j].blue  = x_arr[i*movie_width + j*movie_xdim + 2]*255;
				}
			}

			/* rgb - nR */
			if(type==1){
				if(movie_xdim==3){
					image_dlib[i][j].red   = x_arr[0*movie_width*movie_height + i*movie_width + j]*255;
					image_dlib[i][j].green = x_arr[1*movie_width*movie_height + i*movie_width + j]*255;
					image_dlib[i][j].blue  = x_arr[2*movie_width*movie_height + i*movie_width + j]*255;
				}
			}
			
		}
	}
	TRYCXX( VecRestoreArray(*myvector_Vec, &x_arr) );

	bool saved;
	if(extension == ".jpg"){
		save_jpeg(image_dlib , file_name, 75);
		coutMaster << "jpeg saved: " << file_name << std::endl;
		saved = true;
	}

	if(!saved){
		coutMaster << "ERROR: extension of file '" << extension << "' was not recognized." << std::endl;
	}
}

Vec * movieplotter::get_myvector_Vec(){
	return myvector_Vec;
}

bool movieplotter::get_myvector_loaded(){
	return myvector_loaded;
}

void movieplotter::plot_image(const canvas& c) const{
	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();
	
	unsigned long x_size = this->width();
	unsigned long y_size = this->height();

	/* get the length of the signal */
	int myvector_size;
	TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );

	/* get some plotting variables */
	rectangle rect_pixel;
	int x_coor, y_coor;
	double pixel_size_x = x_size/(double)movie_width;
	double pixel_size_y = y_size/(double)movie_height;
	double pixel_size;
	if(pixel_size_x < pixel_size_y){
		pixel_size = pixel_size_x;
	} else {
		pixel_size = pixel_size_y;
	}

	
	double *values;
	TRYCXX( VecGetArray(*myvector_Vec,&values) );

	/* greyscale */
	if(movie_xdim==1){
		for(int t=0;t<myvector_size;t++){
			y_coor = t/(double)movie_width;
			x_coor = t - y_coor*movie_width;

			rect_pixel.set_left(x_begin + x_coor*pixel_size);
			rect_pixel.set_top(y_begin + y_coor*pixel_size);
			rect_pixel.set_right(x_begin + (x_coor+1)*pixel_size);
			rect_pixel.set_bottom(y_begin + (y_coor+1)*pixel_size);
		
			fill_rect(c,rect_pixel,rgb_pixel(values[t]*255,values[t]*255,values[t]*255));		
		}
	}

	/* rgb */
	if(movie_xdim==3){
		int R = movie_width*movie_height;

		for(int t=0;t<R;t++){
			y_coor = t/(double)movie_width;
			x_coor = t - y_coor*movie_width;

			rect_pixel.set_left(x_begin + x_coor*pixel_size);
			rect_pixel.set_top(y_begin + y_coor*pixel_size);
			rect_pixel.set_right(x_begin + (x_coor+1)*pixel_size);
			rect_pixel.set_bottom(y_begin + (y_coor+1)*pixel_size);

			/* Rn */
			if(type == 0){
				fill_rect(c,rect_pixel,rgb_pixel(values[t*movie_xdim+0]*255,values[t*movie_xdim+1]*255,values[t*movie_xdim+2]*255));
			}

			/* nR */
			if(type == 1){
				fill_rect(c,rect_pixel,rgb_pixel(values[0*R + t]*255,values[1*R + t]*255,values[2*R + t]*255));
			}
		}
	}

	TRYCXX( VecRestoreArray(*myvector_Vec,&values) );
	
}

void movieplotter::set_size(int new_width, int new_xdim, int new_T){

	if(new_width > 0){
		this->movie_width = new_width;
	}

	if(new_xdim==1 | new_xdim == 3){
		this->movie_xdim = new_xdim;
	}

	if(new_T > 0){
		this->movie_T = new_T;
	}

	/* compute height */
	if(myvector_loaded){
		int    myvector_size;
		TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );
		
		this->movie_height = (int)((myvector_size/(double)movie_xdim)/(double)(new_width*new_T));
	
	} else {
		this->movie_height = 0;
	}

	/* repaint ! */
	parent.invalidate_rectangle(rect);
}

int movieplotter::get_width() const {
	return movie_width;
}

int movieplotter::get_height() const {
	return movie_height;
}

int movieplotter::get_T() const {
	return movie_T;
}

double movieplotter::get_scale() const {
	return movie_scale;
}

int movieplotter::get_xdim() const {
	return movie_xdim;
}

std::string movieplotter::get_filename() const{
	return this->filename;
}

void movieplotter::recompute_scale() {
	if(!myvector_loaded){
		movie_scale = 1.0;
	} else {
		unsigned long x_size = this->width();
		unsigned long y_size = this->height();

		double pixel_size_x = x_size/(double)movie_width;
		double pixel_size_y = y_size/(double)movie_height;
		
		if(pixel_size_x < pixel_size_y){
			movie_scale = pixel_size_x;
		} else {
			movie_scale = pixel_size_y;
		}
		
	}
}

std::string cut_filename(const std::string &input){
	std::ostringstream sout;
	boost::filesystem::path p(input);
	sout << p.filename();

	return sout.str();
}

std::string cut_extension(const std::string &input){
	std::ostringstream sout;
	boost::filesystem::path p(input);
	sout << p.extension();

	return sout.str();
}
