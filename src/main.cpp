#define _USE_MATH_DEFINES
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <matplot/matplot.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace matplot;

int add(int i, int j) {
    return i + j;
}

double generate_signal(double frequency, double freq_samp, int samp, int choose)
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
        res = cos(phase);
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
   
    m.attr("__version__") = "dev";
}
