/**
 * \file MeanFieldSolver.hh
 * \author Romain Fournier
 * \date 30.10.18
 **/

#ifndef INC_2DTOOLS_MEANFIELDSOLVER_HH
#define INC_2DTOOLS_MEANFIELDSOLVER_HH

#include "Solver.hh"
#include<eigen3/Eigen/Dense>

/** \class MeanFieldSolver
 *  \brief Implements a Mean Field approach for computing the Green's function of a Material.
**/

class MeanFieldSolver: public Solver {
public:
    MeanFieldSolver(const string &simulation_name_, const string &filename_, const Material &mat_, int k_resolution,
                    double omega_max, double omega_min, double omega_resolution_,
                    const vector<string> &directions_name_, const vector<Vector2d> &starting_k_points_,
                    const vector<Vector2d> &finishing_k_points_);

    void get_Green_function() override;

    /**
     * \brief Compute the Hamiltonian for a given wave vector k
     * @param k
     */
    void compute_H(const VectorXd& k );
private:
    MatrixXcd H_;
public:
    const MatrixXcd &getH_() const;
    /**< Hamiltonian of the system **/
};


#endif //INC_2DTOOLS_MEANFIELDSOLVER_HH
