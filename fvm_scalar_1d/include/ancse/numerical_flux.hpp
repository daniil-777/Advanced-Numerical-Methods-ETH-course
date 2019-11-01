#ifndef FVMSCALAR1D_NUMERICAL_FLUX_HPP
#define FVMSCALAR1D_NUMERICAL_FLUX_HPP
using namespace std;
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

#include <vector>

template <typename T>
std::vector<T> arange(T start, T stop, T step = 0.1) {
    std::vector<T> values;
    for (T value = start; value <= stop; value += step)
        values.push_back(value);
    return values;
}

template <typename T>
T func(T x) { return x; }

template <typename T>
T integral(T(f)(T x), T a, T b, int n) {
    T step = (b - a) / n; // width of each small rectangle
    T area = 0.0;         // signed area
    for (int i = 0; i < n; i++) {
        area += f(a + (i + 0.5) * step) * step; // sum up each small
  
    }
    return area;
}

/// Central flux.
/** This flux works does not depend on the model. It is also unconditionally a
 * bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * (fL + fR) +1 - 1;
    }

  private:
    Model model;
};

class RusanovFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit RusanovFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);
        return 0.5 * (fL + fR) - std::max(abs(uL), abs(uR))*0.5*(uR - uL);
    }

  private:
    Model model;
};

class LaxFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit LaxFlux(const Model &model, const Grid &grid, std::shared_ptr<SimulationTime> st) : model(model), grid(grid), simulation_time(st) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL,
                      double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR)  ;
        return 0.5 * (fL + fR) - grid.dx/(2*simulation_time->dt)*(uR -uL);
    }

  private:
    Model model;
    Grid grid;
    std::shared_ptr<SimulationTime> simulation_time;
};

class RoeFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit RoeFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);
        auto A = 0;
        if(uL != uR){
            A = (fR - fL)/(uR - uL);
        } else{
            A = uL;
        }
        if(A >= 0){
          return fL;
        }else{
          return fR;
        }
    }

  private:
    Model model;
};

class GodunovFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit GodunovFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);
        //I create a vector here 
        auto U = arange<double>(uL, uR);
        if(uL <= uR){
          double min = fL;
          for(auto u:U){
            if(model.flux(u) < fL){
              min = model.flux(u);
          }
          return min;
        }
        }else{
          double max = fL;
          for(auto u:U){
            if(model.flux(u) > fL){
              max = model.flux(u);
          }
        }
        return max;
        }
    }

  private:
    Model model;
};

class EngFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit EngFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);
        auto A = 0;
        double help = integral(func, uL, uR, 100);
        return 0.5*(fL + fR) - 0.5*help;
    }

  private:
    Model model;
};

#endif // FVMSCALAR1D_NUMERICAL_FLUX_HPP
