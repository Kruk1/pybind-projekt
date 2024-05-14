#define _USE_MATH_DEFINES
#include <cmath>
#include <list>
#include <complex>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <matplot/matplot.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace matplot;
using namespace std;

int add(int i, int j) {
    return i + j;
}

std::list<double>* generate_signal(double frequency, double freq_samp, std::list<int> samp, int choose, float b = 0)
{
    int s = size(samp);
    list<double>* res = new list<double>[s];
    for (int i = 0; i < s; i++)
    {
        int curr_sample = samp.front();
        double phase = fmod(curr_sample * 2 * M_PI * frequency / freq_samp, 2 * M_PI);
        if (choose == 1)
        {
            res->push_back(sin(phase));
        }
        else if (choose == 2)
        {
            res->push_back(cos(phase));
        }
        else if (choose == 3)
        {
            res->push_back(fmod(curr_sample * 2 * frequency / freq_samp + 1, 2) - 1);
        }
        else if (choose == 4)
        {
            double con = fmod(curr_sample * 2 * frequency / freq_samp + 1, 2) + (b - 1);
            if (con >= 0.0)
                res->push_back(1);
            else
                res->push_back(0);
        }
        samp.pop_front();
    }
    return res;
}

void show_plot(std::vector<double> x, std::vector<double> y)
{
    tiledlayout(1, 1);
    auto ax1 = nexttile();
    plot(ax1, x, y);
    title(ax1, "Wartosc sygnalu");

    show();
}

std::list<complex<double>> *dft(std::list<complex<double>> input) {
	double s = size(input);
	std::list<complex<double>>* result = new std::list<complex<double>>[s];
	complex<double> ret,roundret;
	complex<double> pow = -2i * M_PI;
	
	int nr = 0;
	double j = 0;
	for (double i = 0; i < s;i++) {
		j = 0;
		ret = 0;
		roundret = 0;
		for (complex<double> z : input) {
			ret += z * exp(pow * i * j / s);
			j++;
		}
		roundret = (round(real(ret))+round(imag(ret))*1i);
		result->push_back(roundret);
		nr++;

	}
	return result;
}

std::list<complex<double>> *inverse_transform(std::list<complex<double>> input) {
	double s = size(input);
	list<complex<double>>* result = new list<complex<double>>[s];
	complex<double> ret, roundret;
	complex<double> pow = 2i * M_PI;

	int nr = 0;
	double j = 0;
	for (double i = 0; i < s;i++) {
		j = 0;
		ret = 0;
		roundret = 0;
		for (complex<double> z : input) {
			ret += 1 / s * z * exp(pow * i * j / s);
			j++;
		}
		roundret = (round(real(ret)) + round(imag(ret)) * 1i);
		result->push_back(roundret);
		nr++;
	}
	return result;
}

std::list<double>* generate_signal_noise(int num, double seed)
{ 
    list<double>* res = new list<double>[num];
    double m = pow(2, 31); 
    double a = 1103515245; 
    double c = 12243; 
    
    double x = fmod(seed, m); 

    for (int i = 0; i < num; i++)
    {
        x = fmod(a * x + c, m);
        res->push_back(x / m);
    }
   
    return res; 
}

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("generate_signal", &generate_signal, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("show_plot", &show_plot, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("dft", &dft);

    m.def("inverse_transform", &inverse_transform);

    m.def("generate_signal_noise", &generate_signal_noise);
   
    m.attr("__version__") = "dev";
}
