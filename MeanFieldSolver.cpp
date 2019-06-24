/**
 * \file MeanFieldSolver.cpp
 * \author Romain Fournier
 * \date 30.10.18
 **/
#include "MeanFieldSolver.hh"
#include <fstream>
#include <iostream>
using namespace std;
MeanFieldSolver::MeanFieldSolver(const string &simulation_name_, const string &filename_, const Material &mat_,
                                 int k_resolution, double omega_max, double omega_min, double omega_resolution_,
                                 const vector<string> &directions_name_, const vector<Vector2d> &starting_k_points_,
                                 const vector<Vector2d> &finishing_k_points_) : Solver(simulation_name_, filename_,
                                                                                       mat_, k_resolution, omega_max,
                                                                                       omega_min, omega_resolution_,
                                                                                       directions_name_,
                                                                                       starting_k_points_,
                                                                                       finishing_k_points_)
{
    H_.setZero(mat_.get_number_elements_in_basis()*mat_.nb_orbitals(),mat_.get_number_elements_in_basis()*mat_.nb_orbitals());
}

void MeanFieldSolver::get_Green_function() {
    // For each chosen direction
    for(int ii(0);ii<directions_name_.size();++ii){
        // Compute the dk
        Vector2d dk;
        for(int jj(0);jj<2;jj++){
            dk(jj)=(finishing_k_points_[ii](jj)-starting_k_points_[ii](jj))/(k_resolution+0.0);
        }

        // For each dk in the direction
        Vector2d k= starting_k_points_[ii];
        ofstream myfile;
        myfile.open(directions_name_[ii]+".csv");
        for(int kk(0);kk<=k_resolution;kk++){
            compute_H(k);
            SelfAdjointEigenSolver <MatrixXcd> mysolver;
            mysolver.compute(H_);
            myfile<<mysolver.eigenvalues().transpose()<<endl;
            k+=dk;
        }
        myfile.close();
    }
}
template <typename T> int my_sign(T val) {
    return (T(0) < val) - (val < T(0));
}
void MeanFieldSolver::compute_H(const VectorXd &k) {
    // For each column of the Hamiltonian
    H_.setZero(mat_.get_number_elements_in_basis()*mat_.nb_orbitals(),mat_.get_number_elements_in_basis()*mat_.nb_orbitals());

    for (unsigned int ii(0);ii<mat_.get_number_elements_in_basis();++ii){
       // cout<<" ii = "<<ii<<endl;
        unsigned int counter(0);
        for(const auto& neighbor : mat_.get_neighborhood(ii)){
            // Hopping from the basis element ii to the one of the neighbor, with amplitude depending on the distance
            // Intra orbitals terms
            for(unsigned int jj(0);jj<mat_.nb_orbitals();++jj){
              //   std::cout<<" New neighbor "<<endl<<neighbor<<endl;
              //   std::cout<<"added phase : "<<endl<<neighbor[0]*mat_.get_lattice().getA1_()+neighbor[1]*mat_.get_lattice().getA2_()<<endl;

                H_(neighbor[2]+jj*mat_.get_number_elements_in_basis(),ii+jj*mat_.get_number_elements_in_basis())+=-mat_.get_t(neighbor[3]-1)*compl_exp((neighbor[0]*mat_.get_lattice().getA1_()+neighbor[1]*mat_.get_lattice().getA2_()).dot(k));
            }
            // inter orbitals terms
            /** todo : make it work for arbitrary number of orbitals, now only for 2 **/
            int hop_from(neighbor[2]+mat_.get_number_elements_in_basis());
            int hop_to(ii);
            Vector2d translation(neighbor[0]*mat_.get_lattice().getA1_()+neighbor[1]*mat_.get_lattice().getA2_());
            double amplitude(mat_.getCross_orbital_(neighbor[3]-1));
            /** not very general */
            if ((neighbor[0]+neighbor[1])==2 || (neighbor[0]+neighbor[1])==-1 ){
                amplitude*=-1;
                cout<<"changing"<<endl;
            }
            if(abs(amplitude)>0.0001) {
                counter+=1;
                cout<<"adding new phase ("<<counter<<"th time) with translation "<<translation.transpose()<<endl;
            }
            H_(hop_to,hop_from)+=-complex<double>(0,2)*amplitude*sin(k.dot(translation));
            H_(hop_from,hop_to)+=complex<double>(0,2)*amplitude*sin(k.dot(translation));
         }
    }
}

const MatrixXcd &MeanFieldSolver::getH_() const {
    return H_;
}
