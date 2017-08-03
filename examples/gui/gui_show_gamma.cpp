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

#define NUMBER_OF_LABELS 11
#define MAX_K 20
#define DEFAULT_TYPE 1
#define DEFAULT_K 1

using namespace dlib;

class gammaplotter : public drawable {
	private:
		int type;
		Vec *myvector_Vec;
		bool myvector_loaded;
		int K;
		int T;
		std::string filename;
		int myvector_size;
		long cut;

		void draw (const canvas& c) const;
		void plot_vector(const canvas& c) const;

		template <typename image_type, typename pixel_type>
		void draw_dotted_line(image_type& c, const point& p1, const point& p2, const pixel_type& val, double size) const;

		template <typename image_type>
		void draw_marker(image_type& c, long x, long y, double width) const;
		
		void fill_label_value(long x, long y, double value, int k) const;
		void fill_label_t(long x, long y, double value) const;
		
		
		label **value_labels;
		label *t_label;
		
	public: 
		gammaplotter(drawable_window& w, int type);
		~gammaplotter();	
		
		void set_plotting_area(rectangle area);
		void load_vector( const std::string& file_name );
		bool get_myvector_loaded();
		Vec *get_myvector_Vec();
		
		void set_K(int new_K);
		int get_K() const;
		int get_T() const;
		void set_newcut(long cut);
		
		std::string get_filename() const;

};

class mymouse_tracker : public draggable {
    private:
		gammaplotter *mygammaplotter;
    
		const long offset;
		std::ostringstream sout;

	public:
		mymouse_tracker(drawable_window& w, gammaplotter& mygammaplotter);
		~mymouse_tracker();

		void show();
		void hide();

		void enable();
		void disable();

		void set_pos(long x, long y);
		void set_main_font(const std::shared_ptr<font>& f);

    protected:
		void on_mouse_move (unsigned long state, long x, long y);

		void draw (const canvas& c) const;
};

class show_gamma_window : public drawable_window {
private:
	mymouse_tracker mtracker; /* gui: for tracking mouse position */

    menu_bar mbar; /* gui: menubar */
	gammaplotter mygammaplotter; /* gui: canvas for drawing */
 
    label **labels_myvector_properties;
	text_field select_K;
	button button_K;

    void load_vector( const std::string& file_name );

	/* menu events */
    void on_menu_file_open ();
    void on_menu_file_quit ();
    void on_menu_help_about();
	
	/* general window events */
	void on_window_resized();
	
	template<class ValueType>
	void set_label_myvector_properties(int label_idx, const std::string &text, ValueType value);
	std::string cut_filename(const std::string &input);
	void on_button_K();
	
	void fill_labels();

public:
    show_gamma_window(int type);
    show_gamma_window(std::string filename, int type, std::string title, int K = 1);
    ~show_gamma_window();

};

//  ----------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

	/* add local program options */
	boost::program_options::options_description opt_problem("GUI_SHOW_GAMMA", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename", boost::program_options::value<std::string>(), "gamma file to load [string]")
		("K", boost::program_options::value<int>(), "number of clusters [int]")
		("type", boost::program_options::value< int >(), "type of input vector [0=TRn, 1=TnR, 2=nTR]")
		("title", boost::program_options::value<std::string>(), "title of window [string]")	;
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* get values provided as console parameters */
	std::string filename, title;
	int K, type;

	consoleArg.set_option_value("filename", &filename, "");
	consoleArg.set_option_value("title", &title, "");
	consoleArg.set_option_value("type", &type, DEFAULT_TYPE);
	consoleArg.set_option_value("K", &K, DEFAULT_K);

    // create our window
    show_gamma_window my_window(filename,type,title,K);

    // wait until the user closes this window before we let the program 
    // terminate.
    my_window.wait_until_closed();

	/* finalize the library */
	Finalize<PetscVector>();
	
    return 0;
}


/* ---------------- implementation -------------- */

show_gamma_window::show_gamma_window(std::string filename, int type, std::string title, int K) : show_gamma_window(type) {
	if(title != ""){
		set_title(title);		
	}

	if(filename != ""){
		load_vector( filename );
	}

	select_K.set_text(std::to_string(K));
	
	on_button_K();

	
}

show_gamma_window::show_gamma_window(int type) : /* All widgets take their parent window as an argument to their constructor. */
        mbar(*this),
        mtracker(*this, mygammaplotter),
        select_K(*this),
        button_K(*this),
        mygammaplotter(*this,type)
{

	/* set the size of window */
    set_size(600,350);
		
	/* allocate labels */
	labels_myvector_properties = new label*[NUMBER_OF_LABELS];
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		labels_myvector_properties[i] = new label(*this);
	}

    /* prepare position of labels of vector properties */
    labels_myvector_properties[0]->set_pos(10,100);
	for(int i=1; i < 11; i++){
		labels_myvector_properties[i]->set_pos(labels_myvector_properties[i-1]->left(),labels_myvector_properties[i-1]->bottom()+20);
	}

	/* window title */
	set_title("Show PETSc Vector utility");
        
	/* create menu bar */
    mbar.set_number_of_menus(2);
    mbar.set_menu_name(0,"File",'F');
    mbar.set_menu_name(1,"Help",'H');

    /* add the entries to the File menu. */
    mbar.menu(0).add_menu_item(menu_item_text("Open",   *this, &show_gamma_window::on_menu_file_open,    'O'));
    mbar.menu(0).add_menu_item(menu_item_separator());
    mbar.menu(0).add_menu_item(menu_item_text("Quit",   *this, &show_gamma_window::on_menu_file_quit,    'Q'));

    /* Add the entries to the Help menu. */
    mbar.menu(1).add_menu_item(menu_item_text("About",  *this, &show_gamma_window::on_menu_help_about,   'A'));

	/* register input for K */
	select_K.set_pos(10,30);
	select_K.set_width(50);
	select_K.set_text("1");
//	select_K.disable();

	/* register button for changing the K */
	button_K.set_click_handler(*this, &show_gamma_window::on_button_K);
	button_K.set_name("apply K");
	button_K.set_pos(10,60);

	mtracker.set_pos(10,320);

	/* arrange the window */
	on_window_resized();

	show();
} 

show_gamma_window::~show_gamma_window(){

	/* destroy labels */
	for(int i=0; i < NUMBER_OF_LABELS; i++){
		free(labels_myvector_properties[i]);
	}
	free(labels_myvector_properties);

	/* close window */
	close_window();
}

void show_gamma_window::on_menu_help_about(){
     message_box("About","This application is for PETSc vector visualisation\n");
}

void show_gamma_window::on_menu_file_open(){
    /* display a file chooser window and when the user choses a file */
    open_existing_file_box(*this, &show_gamma_window::load_vector);
}

void show_gamma_window::on_menu_file_quit(){
	close_window();
}

void show_gamma_window::on_window_resized() {

	/* set new plotting area */
	unsigned long width,height;
    get_size(width,height);
	mygammaplotter.set_plotting_area(rectangle(200,mbar.bottom(),width,height));

	drawable_window::on_window_resized();

}

void show_gamma_window::load_vector( const std::string& file_name ) {
	mygammaplotter.load_vector(file_name);

	if(mygammaplotter.get_myvector_loaded()){
		fill_labels();
	}

	mygammaplotter.enable();
	mtracker.enable();
}

void show_gamma_window::fill_labels() {
	Vec *myvector_Vec = mygammaplotter.get_myvector_Vec();

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
	set_label_myvector_properties(0, "name:    ", cut_filename(mygammaplotter.get_filename()));
	set_label_myvector_properties(1, "K:       ", mygammaplotter.get_K());
	set_label_myvector_properties(2, "T:       ", mygammaplotter.get_T());
	set_label_myvector_properties(3, "size:    ", myvector_size);
	set_label_myvector_properties(4, "norm2:   ", myvector_norm2);
	set_label_myvector_properties(5, "norm1:   ", myvector_norm1);
	set_label_myvector_properties(6, "normInf: ", myvector_normInf);
	set_label_myvector_properties(7, "sum:     ", myvector_sum);
	set_label_myvector_properties(8, "max:     ", myvector_max);
	set_label_myvector_properties(9, "min:     ", myvector_min);
	set_label_myvector_properties(10, "mean:    ", myvector_sum/(double)myvector_size);	
}	

template<class ValueType>
void show_gamma_window::set_label_myvector_properties(int label_idx, const std::string &text, ValueType value){
	std::ostringstream sout;
	sout << std::setprecision(5);
    sout << text << std::setw(20) << value;
    labels_myvector_properties[label_idx]->set_text(sout.str());	
}

std::string show_gamma_window::cut_filename(const std::string &input){
	std::ostringstream sout;
	boost::filesystem::path p(input);
	sout << p.filename();

	return sout.str();
}

void show_gamma_window::on_button_K(){
	int newK = std::stoi(trim(select_K.text()));
	mygammaplotter.set_K(newK);
	
	if(mygammaplotter.get_myvector_loaded()){
		fill_labels();
	}

}


/* ------------------------- graph plotter -------------- */

void gammaplotter::draw(const canvas& c) const {
	/* draw background */
	fill_rect(c,rect,rgb_pixel(255,255,255));
	
	/* plot vector */
	if(myvector_loaded){
		plot_vector(c);
	}

}

gammaplotter::gammaplotter(drawable_window& w, int type): drawable(w) {
	/* create empty vector */
	myvector_Vec = new Vec();
	TRYCXX( VecCreate(PETSC_COMM_WORLD,myvector_Vec) );
	myvector_loaded = false;
	
	/* default number of clusters */
	K = 1;
	T = 0;
	cut = -1;
	this->type = type;
	
	t_label = new label(w);
	t_label->hide();
	t_label->set_z_order(1000000);
	t_label->set_text_color(rgb_pixel(100,100,255));
	
	value_labels = new label*[MAX_K];
	for(int k=0; k < MAX_K; k++){
		value_labels[k] = new label(w);
		value_labels[k]->set_text("");
		value_labels[k]->set_pos(10,50);
		value_labels[k]->set_z_order(1000000);
		value_labels[k]->set_text_color(rgb_pixel(255,0,255));

		value_labels[k]->hide();
	}
	
	enable_events();
}

gammaplotter::~gammaplotter(){
	free(myvector_Vec);

	for(int k=0; k < MAX_K; k++){
		free(value_labels[k]);
	}
	free(value_labels);

	disable_events();
	parent.invalidate_rectangle(rect);
}

void gammaplotter::set_plotting_area(rectangle area){
	rect = area;
}

void gammaplotter::load_vector( const std::string& file_name ) {
	this->filename = file_name;
	
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

	TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );
	this->T = myvector_size/(double)K;


	/* the whole rectangle with plotted graph will be repainted */
	parent.invalidate_rectangle(rect);

}

std::string gammaplotter::get_filename() const{
	return this->filename;
}

Vec * gammaplotter::get_myvector_Vec(){
	return myvector_Vec;
}

bool gammaplotter::get_myvector_loaded(){
	return myvector_loaded;
}

void gammaplotter::plot_vector(const canvas& c) const{
	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();
	
	unsigned long x_size = this->width();
	unsigned long y_size = this->height();

	/* get the length of the signal */
	double myvector_max = 1;
	double myvector_min = 0;

	/* free space parameters */
	double px_min = 0.01*x_size;
	double px_max = 0.99*x_size;
	double py_min = 0.1*y_size;
	double py_max = 0.9*y_size;
	double py_space = 0.1*y_size;

	double py_step = (py_max-py_min - (K-1)*py_space)/(double)K;

	double *values;
	TRYCXX( VecGetArray(*myvector_Vec,&values) );
	for(int k=0;k<K;k++){
		/* coefficients of mapping t to x */
		double ax = (px_max - px_min)/(double)(T - 1);
		double bx = px_max - ax*T;

		/* coefficients of mapping value to y */
		double ay = (py_step)/(double)(myvector_max - myvector_min);
		double by = py_step - ay*myvector_max;

		/* I will use these points */
		point mypoint1;
		point mypoint2;

		/* plot cut */
		double cut_t = (cut-x_begin - bx)/ax;
		if(cut_t >= 0 && cut_t < T){
			mypoint1(0) = cut;
			mypoint1(1) = y_begin + py_min + k*py_step + (py_step -(ay*1 + by)) - 5;

			mypoint2(0) = cut;
			mypoint2(1) = y_begin + py_min + k*py_step + (py_step - (ay*0 + by)) + 5;
		
			if(k>=1){
				mypoint1(1) += k*py_space;
				mypoint2(1) += k*py_space;
			}

			draw_line(c,mypoint1,mypoint2, rgb_pixel(100,100,255));
			
			/* interpolate values */
			int left_t = (int)cut_t;
			int right_t = left_t + 1;
			double av, bv;
			if(type == 0 | type == 1){
				av = (values[right_t*K + k] - values[left_t*K + k])/(double)(right_t - left_t);
				bv = values[right_t*K + k] - av*right_t;
			}
			if(type == 2){
				av = (values[k*T+right_t] - values[k*T+left_t])/(double)(right_t - left_t);
				bv = values[k*T+right_t] - av*right_t;
			}
			double v = av*cut_t + bv;

			double yv = y_begin + py_min + k*py_step + (py_step -(ay*v + by));
			if(k>=1){
				yv += k*py_space;
			}			
			
			this->draw_marker(c, cut, yv, 5);
			
			this->fill_label_value(cut, yv, v, k);
			if(k == K-1){
				this->fill_label_t(cut, y_begin + y_size - 20, cut_t);
			}
		} else {
			/* hide values */
			for(int k=0; k < MAX_K; k++){
				value_labels[k]->hide();
			}
			t_label->hide();
		}

		/* plot x-axis */
		mypoint1(0) = x_begin + ax*0 + bx - 5;
		mypoint1(1) = y_begin + py_min + k*py_step + (py_step -(ay*0 + by));
		mypoint2(0) = x_begin + ax*T + bx + 5;
		mypoint2(1) = y_begin + py_min + k*py_step + (py_step - (ay*0 + by));
		if(k>=1){
			mypoint1(1) += k*py_space;
			mypoint2(1) += k*py_space;
		}
		draw_dotted_line(c,mypoint1,mypoint2, rgb_pixel(0,0,0), 10);

		/* plot y-axis */
		mypoint1(0) = x_begin + ax*0 + bx;
		mypoint1(1) = y_begin + py_min + k*py_step + (py_step -(ay*0 + by)) + 5;
		mypoint2(0) = x_begin + ax*0 + bx;
		mypoint2(1) = y_begin + py_min + k*py_step + (py_step - (ay*1 + by)) - 5;
		if(k>=1){
			mypoint1(1) += k*py_space;
			mypoint2(1) += k*py_space;
		}
		draw_dotted_line(c,mypoint1,mypoint2, rgb_pixel(0,0,0), 10);


		/* plot gamma */
		for(int t=1;t<T;t++){
			mypoint1(0) = x_begin + ax*(t-1) + bx;
			mypoint2(0) = x_begin + ax*t + bx;
			if(type == 0 | type == 1){
				mypoint1(1) = y_begin + py_min + k*py_step + (py_step -(ay*values[(t-1)*K + k] + by));
				mypoint2(1) = y_begin + py_min + k*py_step + (py_step - (ay*values[t*K + k] + by));
			}
			if(type == 2){
				mypoint1(1) = y_begin + py_min + k*py_step + (py_step -(ay*values[k*T+t-1] + by));
				mypoint2(1) = y_begin + py_min + k*py_step + (py_step - (ay*values[k*T+t] + by));
			}
		
			if(k>=1){
				mypoint1(1) += k*py_space;
				mypoint2(1) += k*py_space;
			}
		
			draw_line(c,mypoint1,mypoint2, rgb_pixel(255,0,0));

		}

	}
	TRYCXX( VecRestoreArray(*myvector_Vec,&values) );
	
	
}

template <typename image_type, typename pixel_type>
void gammaplotter::draw_dotted_line(image_type& c, const point& p1, const point& p2, const pixel_type& val, double size) const{
	point diff_vector = p2 - p1;
	double diff_size = sqrt( diff_vector(0)*diff_vector(0) + diff_vector(1)*diff_vector(1) );
	point step_vector = size*(diff_vector/diff_size);
	
	double pass_size = 0;
	point pass_point = p1;
	int i = 0;
	while(pass_size < diff_size){
		if(i % 2 == 0){
			if(pass_size + size > diff_size){
				step_vector = p2 - pass_point;
			}
			draw_line(c,pass_point,pass_point+step_vector,val);
			pass_point += step_vector;
			pass_size += size;
		} else {
			pass_point += 0.5*step_vector;
			pass_size += 0.5*size;
		}
		i++;
	}
}

template <typename image_type>
void gammaplotter::draw_marker(image_type& c, long x, long y, double width) const {
	point mypoint1;
	point mypoint2;	
	
	mypoint1(0) = x-width;
	mypoint1(1) = y;
	mypoint2(0) = x;
	mypoint2(1) = y+width;
	draw_line(c,mypoint1,mypoint2, rgb_pixel(255,0,255));

	mypoint1(0) = x+width;
	mypoint1(1) = y;
	draw_line(c,mypoint1,mypoint2, rgb_pixel(255,0,255));

	mypoint2(0) = x;
	mypoint2(1) = y-width;
	draw_line(c,mypoint1,mypoint2, rgb_pixel(255,0,255));

	mypoint1(0) = x-width;
	mypoint1(1) = y;
	draw_line(c,mypoint1,mypoint2, rgb_pixel(255,0,255));

}

void gammaplotter::fill_label_value(long x, long y, double value, int k) const {
	std::ostringstream sout;
	sout << std::setprecision(5);
    sout << value;
	value_labels[k]->set_text(sout.str());
	value_labels[k]->set_pos(x+8,y-16);
	value_labels[k]->show();
}

void gammaplotter::fill_label_t(long x, long y, double value) const {
	std::ostringstream sout;
	sout << std::setprecision(5);
    sout << value;
	t_label->set_text(sout.str());
	t_label->set_pos(x,y);
	t_label->show();
}

void gammaplotter::set_newcut(long cut){
	unsigned long x_begin = this->left();
	unsigned long y_begin = this->top();
	unsigned long x_size = this->width();
	unsigned long y_size = this->height();


	/* repaint old cut */
	rectangle cut_rect(this->cut-10,y_begin,this->cut+50, y_begin+y_size);
	parent.invalidate_rectangle(cut_rect);

	this->cut = cut;
	
	/* repaint new cut */
	cut_rect.left() = cut-10;
	cut_rect.right() = cut+50;
	parent.invalidate_rectangle(cut_rect);
}

void gammaplotter::set_K(int new_K){
	this->K = new_K;

	if(myvector_loaded){
		TRYCXX( VecGetSize(*myvector_Vec, &myvector_size) );
		this->T = myvector_size/(double)K;
	}
		
	/* repaint ! */
	parent.invalidate_rectangle(rect);
}

int gammaplotter::get_K() const {
	return K;
}

int gammaplotter::get_T() const {
	return T;
}

/* ------------------------- mymouse_tracker -------------- */

mymouse_tracker::mymouse_tracker(drawable_window& w, gammaplotter& mygammaplotter) :
	draggable(w),
	offset(18)
{
	this->mygammaplotter = &mygammaplotter;
	
	set_draggable_area(rectangle(0,0,1,1));
//	auto_mutex M(m); 
//	area = rectangle(0,0,500,500);

	enable_events();
}

mymouse_tracker::~mymouse_tracker(){ 
	disable_events(); 
	parent.invalidate_rectangle(rect); 
}

void mymouse_tracker::set_main_font (const std::shared_ptr<font>& f){
	auto_mutex M(m);
	mfont = f;
}

void mymouse_tracker::set_pos(long x, long y){
	draggable::set_pos(x,y);
}

void mymouse_tracker::show(){
	draggable::show();
}

void mymouse_tracker::hide(){
	draggable::hide();
}

void mymouse_tracker::enable(){
	draggable::enable();
}

void mymouse_tracker::disable(){
	draggable::disable();
}

void mymouse_tracker::on_mouse_move(unsigned long state, long x, long y){
	if (!hidden && enabled){
		parent.invalidate_rectangle(rect);
		draggable::on_mouse_move(state,x,y);

		mygammaplotter->set_newcut(x);
	}
}

void mymouse_tracker::draw(const canvas& c) const{ 
//	fill_rect(c, rect,rgb_pixel(232,228,230));
//	draw_pixel(c, point(click_x,click_y),rgb_pixel(255,0,0));
}

