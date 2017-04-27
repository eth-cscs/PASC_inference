#include "general/solver/generalsolver.h"

namespace pascinference {
namespace solver {

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralSolver &solver){
	output << solver.get_name();
	return output;
}

void GeneralSolver::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN
	
	output_global << this->get_name() << std::endl;
	output_global.synchronize();

	LOG_FUNC_END
}

void GeneralSolver::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": status" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": status" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": content" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": timer" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printshort(std::ostringstream &header, std::ostringstream &values) const{
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

void GeneralSolver::printshort_sum(std::ostringstream &header, std::ostringstream &values) const{
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

std::string GeneralSolver::get_name() const {
	return "GeneralSolver";
}

GeneralData *GeneralSolver::get_data() const{
	return this->data;
}

int GeneralSolver::get_it() const{
	return 0;
}

int GeneralSolver::get_debugmode() const {
	return this->debugmode;
}

void GeneralSolver::set_debugmode(int debugmode){
	this->debugmode = debugmode;
}

int GeneralSolver::get_maxit() const {
	return this->maxit;
}

void GeneralSolver::set_maxit(int maxit) {
	this->maxit = maxit;
}

double GeneralSolver::get_eps() const {
	return this->eps;
}

void GeneralSolver::set_eps(double eps) {
	this->eps = eps;
}



}
} /* end of namespace */

