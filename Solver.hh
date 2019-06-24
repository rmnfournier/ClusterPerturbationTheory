/**
 * \file Solver.hh
 * \author Romain Fournier
 * \date 25.10.18
 **/

#ifndef INC_2DTOOLS_SOLVER_HH
#define INC_2DTOOLS_SOLVER_HH

#include <string>
#include <ostream>
#include "Material.hh"

using namespace std;

/** \class Solver
   *\brief Abstract Class for computing the Green's function of a Material.
**/
class Solver {
public:
    Solver(const string &simulation_name_, const string &filename_, const Material &mat_, int k_resolution,
           double omega_max, double omega_min, double omega_resolution_, const vector<string> &directions_name_,
           const vector<Vector2d> &starting_k_points_, const vector<Vector2d> &finishing_k_points_) : simulation_name_(
            simulation_name_), filename_(filename_), mat_(mat_), k_resolution(k_resolution), omega_max(omega_max),
                                                                                                      omega_min(
                                                                                                              omega_min),
                                                                                                      omega_resolution_(
                                                                                                              omega_resolution_),
                                                                                                      directions_name_(
                                                                                                              directions_name_),
                                                                                                      starting_k_points_(
                                                                                                              starting_k_points_),
                                                                                                      finishing_k_points_(
                                                                                                              finishing_k_points_) {}


        /**
         * \method get_Green_function
         * \brief This method saves the Green's function of mat_ in the file filename_.
         * */
        virtual void get_Green_function()=0;
        /**
         * Method used to compatibility problems
         * @param x
         * @return exp(i*x)
         */
        complex<double> compl_exp(double x){return cos(x)+complex<double>(0,1)*sin(x);}

    protected:
        string simulation_name_; /**< Name of the current simulation */
        string filename_; /**< Name of file in which the results will be saved */
        Material mat_; /**< Material for which the Green's function must be computed */
        int k_resolution; /**< Number of k points */
        double omega_max; /**< Maximum energy range*/
        double omega_min; /**< Minimum of the energy range */
        double omega_resolution_;
protected:
    /**< Number of energy points */

        vector<string> directions_name_;/**< Name of the direction in which the ith Green's function is compute. */
        vector<Vector2d> starting_k_points_; /**< First k_point of the computation. */
        vector<Vector2d> finishing_k_points_; /**< Last k_point of the computation. */
};


#endif //INC_2DTOOLS_SOLVER_HH
