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

int add(int i, int j) {
    return i + j;
}

double generate_signal(double frequency, double freq_samp, int samp, int choose, float b = 0)
{
    double res;
    double phase = fmod(samp * 2 * M_PI * frequency / freq_samp, 2 * M_PI);
    if (choose == 1)
    {
        res = sin(phase);
        return res;
    }
    else if(choose == 2)
    {
        res = cos(phase);
        return res;
    }
    else if (choose == 3)
    {
        res = fmod(samp * 2 * frequency / freq_samp + 1, 2) - 1;
        return res;
    }
    else if (choose == 4)
    {
        res = fmod(samp * 2 * frequency / freq_samp + 1, 2) + (b - 1);
        if (res >= 0.0)
            res = 1;
        else
            res = 0;
        return res;
    }
}

void show_plot(std::vector<double> x, std::vector<double> y)
{
    tiledlayout(1, 1);
    auto ax1 = nexttile();
    plot(ax1, x, y);
    title(ax1, "Wartosc sygnalu");

    show();
}

const double PI = atan(1.0) * 4;

std::list<complex<double>> *dft(std::list<complex<double>> input) {
	double s = size(input);
	std::list<complex<double>>* result = new std::list<complex<double>>[s];
	complex<double> ret,roundret;
	complex<double> pow = -2i * PI;
	
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
	complex<double> pow = 2i * PI;

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

   
    m.attr("__version__") = "dev";
}
