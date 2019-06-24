/**
 * \file Lattice.cpp
 * \author Romain Fournier
 * \date 30.10.18
 **/
#include "Lattice.hh"
#include <iostream>
#include <algorithm>

using namespace std;
void Lattice::set_next_nearest_neigbors(unsigned int number_neighbors) {
    // We start by creating "pseudo_neighbors", which are neighbors with an extra dimension representing the euclidian distance
    vector<Matrix<double,5,1>> pseudo_neighbors;
    vector<vector<Matrix<double,5,1>>> pseudo_neighborhood;

    for(size_t i(0);i<get_number_elements_in_basis();i++){
        pseudo_neighbors.clear();
       // cout<<"starting "<<i<<endl;
        for(const auto& neighbor : neighbors_ [i]){
            //cout<<"new Vector"<<endl<<neighbor<<endl;
            Matrix<double,5,1> handler;
            //copy the first 4 dimensions
            for(size_t j(0);j<4;j++){
                handler(j)=neighbor(j);
            }
            // the fifth dimension is the distance
            handler(4)=compute_distance(neighbor,i);
            //cout<<"new distance "<<compute_distance(neighbor,i)<<endl;
            pseudo_neighbors.emplace_back(handler);
        }
        pseudo_neighborhood.push_back(pseudo_neighbors);
    }



    for(int basis(0);basis<get_number_elements_in_basis();++basis){
        // For each ramification
        for(int ii(2);ii<=number_neighbors;ii++){
            // We look for the neighbors located at a distance i-1
            for (const auto neighbor : pseudo_neighborhood[basis]){
                if(abs(neighbor[3]-(ii-1))<1e-5){

                    // Now we have to follow the 1st neighbors of this atom
                    for(const auto& next_neighbor : pseudo_neighborhood[neighbor[2]]){
                        //cout<<"My next neighbor : "<<next_neighbor.transpose()<<endl;

                        if(abs(next_neighbor[3]-1)<1e-5){
                            // Find the absolute path to reach this atom
                            Matrix<double,5,1> potential_neighbor;
                            potential_neighbor.setZero();
                            potential_neighbor=neighbor+next_neighbor;
                           // cout<<" sum : "<<potential_neighbor.transpose()<<endl;
                            //the position in basis is not additive
                            potential_neighbor[2]=next_neighbor[2];
                            // the range is the current loop
                            potential_neighbor[3]=ii;
                            //cout<<"Neighbor under investigation : "<<potential_neighbor.transpose()<<endl;
                            potential_neighbor[4]=compute_distance(potential_neighbor,basis);
                            // Now we have to make sure that this atom is not in the current list (no close loops)
                            bool found(false);
                            //Check if this is not the origine
                           if (!( abs(potential_neighbor[0])<1e-5 && abs(potential_neighbor[1])<1e-5 && abs(potential_neighbor[2]-basis)<1e-5) ) {
                                for ( auto& present : pseudo_neighborhood[basis]){
                                    if(abs(potential_neighbor[0]-present[0])<1e-5 &&abs(potential_neighbor[1]-present[1])<1e-5 && abs(potential_neighbor[2]-present[2])<1e-5){
                                      //should never be used, but nevermind
                                      if(present[3]>ii){
                                          present[3]=ii;
                                      }
                                      found=true;
                                      break;
                                    }
                                }
                                if(!found){
                                    pseudo_neighborhood[basis].push_back(potential_neighbor);

                                }
                            }

                        }
                    }
                }
            }
        }
        // Now we must add the neighbors to the lattice
        // First we make a list of all distances
        vector<double> distances;
        distances.empty();
        for(const auto& n : pseudo_neighborhood[basis]){
            bool already_here(false);
            for(const auto& d : distances){
                if(abs(d-n[4])<1e-9){
                    already_here=true;
                    break;
                }
            }
            if(!already_here){
                distances.push_back(n[4]);
            }
        }



        // Then We sort that list
        sort(distances.begin(),distances.end());

        // for each pseudo_neighbor of the current basis
        cout<<"basis "<<basis<<endl;
        for(const auto& n : pseudo_neighborhood[basis] ){
            //go through the list and find the position
            unsigned int counter(0);
            for(const auto d : distances){
                //if the position is less than the number of neighbor, add it
                counter++;
                if(abs(n(4)-d)<1e-5 && counter<=number_neighbors){
                    Neighbor new_neighbor;
                    for(size_t jn(0);jn<3;jn++){
                        new_neighbor(jn)=int(n(jn));
                    }
                    new_neighbor(3)=counter;
                    neighbors_[basis].push_back(new_neighbor);
                }
            }
        }
    }


}
double Lattice::compute_distance(const Neighbor& n,int b){
    return (n(0)*a1_+n(1)*a2_+basis_[n(2)]-basis_[b]).norm();
}
double Lattice::compute_distance(const Matrix<double,5,1>& n,int b){
    return (n(0)*a1_+n(1)*a2_+basis_[n(2)]-basis_[b]).norm();
}
