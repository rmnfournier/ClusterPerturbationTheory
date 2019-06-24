/**
 * \file Lattice.hh
 * \author Romain Fournier
 * \date 25.10.18
 * \include eigen3
 **/

#ifndef INC_2DTOOLS_LATTICE_HH
#define INC_2DTOOLS_LATTICE_HH

#include <eigen3/Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;

/**
 * \typedef Neighbor
 * \brief represent a relation between two atoms
 *  This typedef is used with the following convention in mind:
 *  <#times a1,#times a2, #position in the basis,#distance from first atom (in terms of jumps)>
 *  for example, <1,-1,3,3> means that to reach our neighbor, we have to go into a site with a a1-a2 translation, and choose the 3rd element of the basis, and that this atom is 3 jumps away from the origine
 */
typedef Vector4i Neighbor;
/**
 * \typedef Neighborhood
 * \brief represent the relation between all atoms
 *  This typedef is used with the following convention in mind:
 *  Neighborhood (i,j) represents the jth neighbor of the ith element in the basis.
 */
typedef vector<vector<Neighbor>> Neighborhood ;

/** \class Lattice
   *\brief Class representing a 2D lattice
   *  This class handles interactions with a general 2D lattice
   **/
class Lattice {
    public:
        Lattice (Vector2d a1, Vector2d a2,vector<Vector2d> basis,Neighborhood neighbors) : a1_(std::move(a1)),a2_(std::move(a2)),basis_(std::move(basis)),neighbors_(std::move(neighbors)){}
        Lattice(){}
        void initialize(Vector2d a1, Vector2d a2,vector<Vector2d> basis,Neighborhood neighbors){
            a1_=(std::move(a1));a2_=(std::move(a2));basis_=(std::move(basis));neighbors_=(std::move(neighbors));
        }
        /**
         * returns the number of atoms in the basis
         * @return number of atoms in the basis
         */
        unsigned int get_number_elements_in_basis() const {return basis_.size();}

        /**
         * \brief Get the neighborhood of the ith site
         * @param i
         * @return neighbors of the ith site
         */
        vector<Neighbor> get_neighborhood(unsigned long i) const {return neighbors_[i];};
    /**
     * \brief starting with already set-up first neighbors, this function add the references to the number_neighbors next nearest neighbors
     * @param number_neighbors
     */
        void set_next_nearest_neigbors(unsigned int number_neighbors);

        /**
         * Compute the Euclidian distance
         * @param n
         * @param b
         * @return
         */
        double compute_distance(const Neighbor& n,int b);
        double compute_distance(const Matrix<double,5,1>& n,int b);

        const Vector2d &getA1_() const {
        return a1_;
    }

    const Vector2d &getA2_() const {
        return a2_;
    }

    Vector2d getG1_() const {
        Matrix2d R ;
        R<<0.,-1.,
            1,0;
        return 2*M_PI*R*a2_/a1_.dot(R*a2_);
    }

    Vector2d getG2_() const {
        Matrix2d R ;
        R<<0.,-1.,
                1,0;

        return 2*M_PI*R*a1_/a2_.dot(R*a1_);
    }

    Vector2d getBasis(const unsigned int i) const {return basis_[i];}

    double get_surface() const { return  abs(a1_(0)*a2_(1)-a1_(1)*a2_(0));}
private:

        Vector2d a1_; /**< first lattice vector */
        Vector2d a2_; /**< second lattice vector */
        vector<Vector2d> basis_;  /**< list of all vectors describing the position of the atoms in one lattice site */
        Neighborhood neighbors_; /**< list of all vectors describing the neighbors of the atoms in the basis */


};
/** \class HoneyCombLattice
   *\brief HoneyCombLattice, with a the interatomic distance

   **/
class HoneyCombLattice :public Lattice {
public:
    explicit HoneyCombLattice(double a){
        // Lattice vectors
        Vector2d a1,a2;
        a1<<3,sqrt(3);
        a2<<0,-2*sqrt(3);
        a1*=a/2.;
        a2*=a/2.;

        vector<Vector2d> basis;
        //Sublattice A
        Vector2d sub_a;
        sub_a<<0.,0.;
        //Sublattice B
        Vector2d sub_b;
        sub_b<<a,0;
        basis.push_back(sub_a);
        basis.push_back(sub_b);

        // Neighbors of sublattice A
        vector<Neighbor> n_a;
        Neighbor neighbor;
        neighbor<<0,0,1,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<-1,-1,1,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<-1,0,1,1;
        n_a.push_back(neighbor);

        // Neighbors of sublattice B
        vector<Neighbor> n_b;
        neighbor.setZero();
        neighbor<<0,0,0,1;
        n_b.push_back(neighbor);

        neighbor.setZero();
        neighbor<<1,1,0,1;
        n_b.push_back(neighbor);

        neighbor.setZero();
        neighbor<<1,0,0,1;
        n_b.push_back(neighbor);

        Neighborhood nh;
        nh.push_back(n_a);
        nh.push_back(n_b);

        initialize(a1,a2,basis,nh);
    }
};

class TriangularLattice : public Lattice{
public:
    explicit TriangularLattice(double a){
        // Lattice vectors
        Vector2d a1,a2;
        a1<<3,sqrt(3);
        a2<<3,-sqrt(3);

        vector<Vector2d> basis;

        Vector2d sub_a;
        sub_a<<0.,0.;
        basis.push_back(sub_a);

        // Neighbors of sublattice A
        vector<Neighbor> n_a;
        Neighbor neighbor;
        neighbor<<1,0,0,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<-1,0,0,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<0,1,0,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<0,-1,0,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<1,-1,0,1;
        n_a.push_back(neighbor);

        neighbor.setZero();
        neighbor<<-1,1,0,1;
        n_a.push_back(neighbor);

        Neighborhood nh;
        nh.push_back(n_a);

        initialize(a1,a2,basis,nh);

    }
};

#endif //INC_2DTOOLS_LATTICE_HH
