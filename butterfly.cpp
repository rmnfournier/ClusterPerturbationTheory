//
// Created by romain on 20/11/18.
//

/**
 * \file test_solver.cpp
 * \author Romain Fournier
 * \date 26.10.18
 **/


#include "CPTSolver.hh"
#include <iostream>
using namespace std;
#include <omp.h>
#include <mpi/mpi.h>

int main(int argc, char * argv[]) {
    MPI_Init(&argc,&argv);

    int prank,psize;

    MPI_Comm_rank(MPI_COMM_WORLD,&prank);
    MPI_Comm_size(MPI_COMM_WORLD,&psize);
    //Define the Material
    Graphene graphene;

    // Construct the CPTSolver
    vector<string>names;
    vector<Vector2d>starting_k;
    vector<Vector2d>finishing_k;
    Vector2d k_f;
    Vector2d k_s;
    //Define the positions of interest
    Vector2d G2 = graphene.get_lattice().getG2_();
    //cout<<G2<<endl;
    //cout<<graphene.get_lattice().getG1_()<<endl;
    double norme_G2=G2.norm();
    double K_Y=norme_G2/2.;
    double K_X =tan(30./360.*2*M_PI)*K_Y;

    //1st direction
    names.emplace_back("Gamma_K")       ;
    k_s.setZero();
    k_s<<0,0;
    k_f.setZero();
    k_f<<2*M_PI/3./0.142,2*M_PI/3./0.142/sqrt(3);
    starting_k.push_back(k_s);
    finishing_k.push_back(k_f);

    //2nd direction
    /*names.emplace_back("K_M");
    k_s.setZero();
    k_s=k_f;
    k_f.setZero();
    k_f<<2*M_PI/3/0.142,0;
    starting_k.push_back(k_s);
    finishing_k.push_back(k_f);

    //3rd direction
    names.emplace_back("M_Gamma");
    k_s.setZero();
    k_s=k_f;
    k_f.setZero();
    k_f<<0.,0.;
    starting_k.push_back(k_s);
    finishing_k.push_back(k_f);*/
  //  cout<<"Starting Butterflies"<<endl;

    int q = 104;
    //cout<< " process "<<prank<<" goes from "<<prank*q/psize+1<<" to "<<((prank+1)*q/psize)<<endl;

    for(int p(prank*q/psize+1);p<=((prank+1)*q/psize);p++){
        cout<<"starting p = "<<p<<" from thread "<<prank<<endl;
        CPTSolver solver("butterfly_2c_q_"+to_string(q), "/home/romain/Documents/EPFL/C3MP/2DTools/Results/2cluster_u1_t1/butterfly_2c_q_"+to_string(q)+"_p_"+to_string(p), graphene, 1,4, -4, 200, names,starting_k,finishing_k,2,1,2, p,q);
        solver.get_Green_function();
    }
    MPI_Finalize();
    return 0;
}