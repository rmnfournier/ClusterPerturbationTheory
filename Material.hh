/**
 * \file Material.hh
 * \author Romain Fournier
 * \date 25.10.18
 * \include eigen3
 **/

#ifndef INC_2DTOOLS_MATERIAL_HH
#define INC_2DTOOLS_MATERIAL_HH

#include <eigen3/Eigen/Dense>
#include <vector>
#include "Lattice.hh"

using namespace std;

/** \class Material
   *\brief Class representing a 2D Material
   *  This is nothing but a lattice with physical properties
**/
class Material {
    public:
        /**
         * \brief returns the number of elements in the basis of the lattice
         * @return number of atoms
         */
        unsigned int get_number_elements_in_basis() const;
        /**
         * \brief Get the neighborhood of the ith site
         * @param i
         * @return neighbors of the ith site
         */
        vector<Neighbor> get_neighborhood(unsigned long i) const {return lat_.get_neighborhood(i);};
        /**
         * \Brief first n.n. hopping amplitude
         * @return first n.n. hopping amplitude
         */
        double get_t(unsigned int i=0) const {return t_[i];}
        double getCross_orbital_(unsigned int i=0) const { return cross_orbital_[i]; }
        unsigned int nb_orbitals() const {return orbitals_;}
        double get_u() { return U_;}
        const Lattice& get_lattice() const {return lat_;}
    protected:
        vector<double> t_; /**< list of the hopping amplitudes (1st element corresponds to first nearest neighbors, 2nd to next nearest neig. etc. */
        vector<double> cross_orbital_;/**< hopping term between orbitals **/
        double U_; /**< force of the electron-electron repulsion for doubly occupied sites */
        Lattice lat_;/**< underlying lattice */

        unsigned int orbitals_;/** number of orbitals per site**/
};

/** \class Graphene
   *\brief Class representing Graphene
   *  First nearest neighbors only, no interaction U
**/
class Graphene:public Material{
public:
    Graphene(){
        /*t_.push_back(2.97);
        t_.push_back(0.073);
        t_.push_back(0.33);*/
        t_.push_back(1.0);
        U_=1.0;
        orbitals_=1;
        cross_orbital_.push_back(0);
        cross_orbital_.push_back(0);
        cross_orbital_.push_back(0);

        lat_=HoneyCombLattice(0.142);

    }
};

class TwistedBilayerGraphene:public Material{
public:
    explicit TwistedBilayerGraphene(double twist_angle){
        /** Maximally Localized Wannier Orbitals and the Extended Hubbard Model for Twisted Bilayer Graphene (Koshino 2018) **/
//        t_.push_back(0.331);t_.push_back(0.016); t_.push_back(0.036); t_.push_back(0.119); t_.push_back(-0.01); // meV;
        t_.push_back(0.331);t_.push_back(0.0); t_.push_back(0.0); t_.push_back(0.119); t_.push_back(-0.01); // meV;

        //t_.push_back(2);t_.push_back(0.02);t_.push_back(0.0);t_.push_back(0.0); t_.push_back(0.00); // meV;

        /** We have two different orbitals(px,py) **/
        orbitals_=2;
        /** Hopping occure only for second nn **/
        cross_orbital_.push_back(0);cross_orbital_.push_back(0.0);cross_orbital_.push_back(0);cross_orbital_.push_back(0);cross_orbital_.push_back(0.097);
       // cross_orbital_.push_back(0);cross_orbital_.push_back(0.0);cross_orbital_.push_back(0);cross_orbital_.push_back(0);cross_orbital_.push_back(2);

        /** No interaction at the moment **/
        U_=0;
        /** Set up a lattice with 5 nearest neighbors **/
        lat_=HoneyCombLattice(0.246/2./sin(twist_angle/2.)); //nm
        lat_.set_next_nearest_neigbors(5);
    }
};

class ToyMaterial:public Material{
public:
    ToyMaterial(){
        t_.push_back(1.);
        U_=0;
        orbitals_=1;
        lat_=TriangularLattice(1.);
    }
};

#endif //INC_2DTOOLS_MATERIAL_HH
