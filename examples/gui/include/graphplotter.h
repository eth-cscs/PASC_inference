#ifndef GRAPHPLOTTER_H
#define	GRAPHPLOTTER_H

#ifndef USE_DLIB
 #error 'this util is for DLIB'
#endif

#include <dlib/gui_widgets.h>
#include <sstream>
#include <string>

using namespace dlib;

class graphplotter : public drawable {
	private:
		std::string xlabel;
		std::string ylabel;
	
		std::vector<double> xvalues;
		std::vector<double> yvalues;

		bool xvalues_loaded;
		bool yvalues_loaded;

		void draw (const canvas& c) const;
		void plot_graph(const canvas& c) const;
		void sort_values();

		int xsize;
		double xmax;
		double xmin;

		int ysize;
		double ymax;
		double ymin;

	public: 
		graphplotter(drawable_window& w);
		~graphplotter();	
		
		void set_plotting_area(rectangle area);

		void set_xvalues(std::vector<double> xvalues);
		void set_yvalues(std::vector<double> yvalues);

		std::vector<double> get_xvalues() const;
		std::vector<double> get_yvalues() const;

		void set_xlabel(std::string xlabel);
		void set_ylabel(std::string ylabel);

		std::string get_xlabel() const;
		std::string get_ylabel() const;

		bool get_xvalues_loaded() const;
		bool get_yvalues_loaded() const;

};

/* -- implementation -- */

void graphplotter::draw(const canvas& c) const {
	/* draw background */
	fill_rect(c,rect,rgb_pixel(255,255,255));
	
	/* plot vector */
	if(xvalues_loaded && yvalues_loaded){
		plot_graph(c);
	}

}

graphplotter::graphplotter(drawable_window& w): 
			drawable(w)
{
	xvalues_loaded = false;
	yvalues_loaded = false;

	enable_events();
}

graphplotter::~graphplotter(){
	disable_events();
	parent.invalidate_rectangle(rect);
}

void graphplotter::set_plotting_area(rectangle area){
	rect = area;
}

void graphplotter::set_xvalues(std::vector<double> xvalues){
	this->xvalues = xvalues;
	xvalues_loaded = true;

	/* get the properties of input values */
	xsize = xvalues.size();
	xmax = *std::max_element(xvalues.begin(), xvalues.end());
	xmin = *std::min_element(xvalues.begin(), xvalues.end());

	if(xvalues_loaded && yvalues_loaded){
		sort_values();	
		parent.invalidate_rectangle(rect);
	}
}

void graphplotter::set_yvalues(std::vector<double> yvalues){
	this->yvalues = xvalues;
	yvalues_loaded = true;

	/* get the properties of input values */
	ysize = yvalues.size();
	ymax = *std::max_element(yvalues.begin(), yvalues.end());
	ymin = *std::min_element(yvalues.begin(), yvalues.end());

	if(xvalues_loaded && yvalues_loaded){
		sort_values();	
		parent.invalidate_rectangle(rect);
	}
}

void graphplotter::sort_values(){
	//TODO: sort vectors with respect to x!

}

std::vector<double> graphplotter::get_xvalues() const {
	return xvalues;
}

std::vector<double> graphplotter::get_yvalues() const {
	return yvalues;
}

void graphplotter::set_xlabel(std::string xlabel) {
	this->xlabel = xlabel;
}

void graphplotter::set_ylabel(std::string ylabel) {
	this->ylabel = ylabel;
}

std::string graphplotter::get_xlabel() const {
	return this->xlabel;
}

std::string graphplotter::get_ylabel() const {
	return this->ylabel;
}

bool graphplotter::get_xvalues_loaded() const {
	return xvalues_loaded;
}

bool graphplotter::get_yvalues_loaded() const {
	return yvalues_loaded;
}


void graphplotter::plot_graph(const canvas& c) const{
	unsigned long wx_begin = this->left();
	unsigned long wy_begin = this->top();
	
	unsigned long wx_size = this->width();
	unsigned long wy_size = this->height();

	/* free space parameters */
	double px_min = 0.1*wx_size;
	double px_max = 0.9*wx_size;
	double py_min = 0.1*wy_size;
	double py_max = 0.9*wy_size;

	/* coefficients of mapping t to x */
	double ax = (px_max - px_min)/(double)(xmax - xmin);
	double bx = px_max - ax*xmax;

	/* coefficients of mapping value to y */
	double ay = (py_max - py_min)/(double)(ymax - ymin);
	double by = py_max - ay*ymax;

	/* I will use these points */
	point mypoint1;
	point mypoint2;

	for(int t=0;t<xsize-1;t++){
		mypoint1(0) = wx_begin + ax*xvalues[t] + bx;
		mypoint1(1) = wy_begin + (wy_size - (ay*yvalues[t] + by));

		mypoint2(0) = wx_begin + ax*xvalues[t+1] + bx;
		mypoint2(1) = wy_begin + (wy_size - (ay*yvalues[t+1] + by));
		
		draw_line(c,mypoint1,mypoint2, rgb_pixel(0,0,255));

	}

	
}


#endif
