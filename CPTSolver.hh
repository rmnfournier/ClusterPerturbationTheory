/**
 * \file CPTSolver.hh
 * \author Romain Fournier
 * \date 25.10.18
 **/

#ifndef INC_2DTOOLS_CPTSOLVER_HH
#define INC_2DTOOLS_CPTSOLVER_HH

#include "Solver.hh"
#include "BitOperationsHelper.hh"
#include <algorithm>
#define A 0
#define B 1
/**
 * \typedef TridiagonalMatrix
 * \brief first element is the diagonal terms, second the offdiagonals
 */
typedef vector<vector<double>> TridiagonalMatrix;
/**
 * \typedef ClusterPosition
 * \brief convention : position (m,n,b) corresponds to site b+m*b +b*n*m
 */
typedef Matrix<int,3,1> ClusterPosition;

/** \class CPTSolver
   *\brief Implements a Cluster Perturbation theory approach for computing the Green's function of a Material.
**/
class CPTSolver :public Solver{
public:
    CPTSolver(const string &simulation_name_, const string &filename_, const Material &mat_, int k_resolution,
              double omega_max, double omega_min, double omega_resolution_, const vector<string> &directions_name_,
              const vector<Vector2d> &starting_k_points_, const vector<Vector2d> &finishing_k_points_, unsigned int M_, unsigned int N_, unsigned long filling_,int p=0,int q=1,
              unsigned int magnetic_cluster=0);
        void get_Green_function() override;
        /**
         * \brief Update the coefficients representing the intracluster Green's function
         */
        void update_Green_intra_coefficients();
        /**
         * \brief Update V_ to take into acccount the impact of the outsideworld on the cluster for wavevector k
         * @param k vector for which we compute the perturbation
         */
        void update_perturbation(const Vector2d& k);
        /**
         * \brief Potential part of update_perturbation method
         * @param k
         */
        void update_outside_potential(const Vector2d& k);


        /**
         * \brief Return the dimension of our vector space
         * @return the dimension of the Fock Space
         */
        unsigned long get_Fock_dim() const;
        /**
         * \brief Return a vector of integer going from 0 to number of sites-1
         * @return a list of all sites
         */
        vector<int> get_sites_list() const;
        /**
         * \brief Return the number of different sites in one cluster
         * @return number of sites
         */
        unsigned long get_number_of_sites() const;

        /**
         * \brief returns the dimensionality of the Geen's function, which is N_*M_*elements in basis * 2 spins
         * @return
         */
        unsigned int get_occupation_states_dim() const{ return mat_.get_number_elements_in_basis()*M_*N_*2;}
        /**
         * \brief Compute the internal Green's function from its coefficients stored in memory
         * @param omega frequency
         * @return Matr
         */
        void compute_G_intern(double omega);
        /**
         * \brief Helps the function compute_G_intern to compute the elements in the diagonal or Gsym_ab  otherwise
         * @param i ith line element
         * @param j jth column element
         * @param omega frequency
         * @return G_intern_ii or Gsym_ij
         */
        complex<double> compute_G_intern_diag_element(int i,int j, complex<double> omega);
        complex<double> compute_G_intern_diag_element_helper(int i,int j,complex<double> omega,bool electron);
        /**
         * \brief Helps the function compute_G_intern to compute the offdiagonal elements. Need compute_G_intern_diag_element to be run before.
         * @param i line
         * @param j column
         * @param omega frequency
         * @return G_intern_ij
         */
        complex<double> compute_G_intern_offdiag_element(int i,int j, complex<double> omega);

        /**
         * \brief Compute the tridiagonal form of intracluster Hamiltonian.
         * @param phi_0 starting point of lanczos procedure
         * @param a storage for the diagonal terms
         * @param b storage for the offdiagonals terms
         * @param tolerance tolerance for the maximal error of the eigenvalues
         * @param full_orthogonalization if true, ensures the full orthogonalization. Can be memory/time consuming
         * @param verbose output some informations for debugging purposes
         */
        void get_lanczos_coefficients_from_states(const VectorXcd& phi_0 ,vector<double>& a,vector<double>& b,double tolerance=1e-15,bool full_orthogonalization=true,bool reconstruct=false,bool early_stop=true);
        /**
         * \brief returns \f$ H_{intra}|q>\f$
         * @param q vector on which we whant to apply the intracluster Hamiltonian
         * @return \f$ H_{intra}|q>\f$
         */
        VectorXcd apply_intra_H(const VectorXcd& q);
        VectorXcd apply_intra_H(const VectorXcd& q,const VectorXd& k);

    /**
         * \brief give all states that can be reached by one hopping, starting from state i
         * @param i starting site
         * @return all states "one hopping" away form i
         */
        vector<unsigned long> get_hopping_from_i_states(unsigned long i);
        /**
         * \brief return true if the sites i and j are conected within the cluster
         * @param i first site
         * @param j second site
         * @return true if the sites i and j are conected within the cluster
         */
        bool are_connected(unsigned long i, unsigned long j);

        /**
         * \brief convert site i into the representation position
         * @param i site
         * @return ClusterPosition
         */
        ClusterPosition get_cluster_position(unsigned long i);
        /**
         * \brief get the neighbors of site i outside the cluster as well as the vector relating the two clusters
         * @param i site of interest
         * @return vector of pairs, first element being the site number of the neighbor in the other cluster, the second one the translation vector to reach the other cluster
         */
        vector<std::pair<unsigned long,Vector2d>> get_translation_for_neighbors_from_other_cluster(unsigned long i,bool magnetic_cluster=false);
        /**
         * \brief convert the representation position into site representation
         * @param pos
         * @return
         */
        int get_site_from_position (ClusterPosition pos);
        /**
         * \brief add the hopping term hopping to all <phi|V|psi> such that phi and psi differ by the sites i and j, with psi having a 1 at site i and at site j
         * @param i position of the electron before hopping
         * @param j position of the electron after hopping
         * @param hopping phase acquire by the electron in the fourier space
         */
        void add_perturbative_hopping_between(int i,int j,complex<double> hopping);


        /**
         * \brief Give back the Green's function with wave vector corresponding to the real lattice, not the cluster
         * Require that G_ has been already computed for the cluster. Note that the k for which we compute G_ is the k corresponding to the real lattice
         */
        MatrixXcd compute_G_periodic(const Vector2d& k);

        /**
         * \brief initialize the list of available states
         */
        void init_states();

        /**
         * \brief compute the Ground state of the Hamiltonian
         */
        void compute_ground_state();
        Vector2d get_vector_from_cluster_position(const ClusterPosition& cp) const ;
            const vector<unsigned long> &getStates_() const {
                    return states_;
            }
        /**
         * \brief Return the dimension of the hamiltonian (number of states with current filling)
         * @return
         */
        unsigned long working_dim(){ return states_.size();}
        /**
         * Correspondance between the integer representation of a state and its numerotation
         * @param s
         * @return on which dimension is state s
         */
       unsigned long get_dim_from_state(unsigned long s){
            for(unsigned long ii(0);ii<working_dim();++ii){
                if(s==states_[ii]){
                    return ii;
                }
            }
        }
        /**
         * \brief Compute adagger_i |groundstate>
         * @param i
         * @return
         */
        VectorXcd creator(unsigned long i);
        /**
         * \brief Compute a_i |groundstate>
         * @param i
         * @return
         */
        VectorXcd destructor(unsigned long i);
        /**
         * \brief return the real space coordinate of the site
         * @param site site of interest
         * @return
         */
        Vector2d get_absolute_position(unsigned int site);
        /**
         * \brief compute exp(-ip/q integral (A.dr))
         * @param i first site
         * @param j second site
         * @return
         */
        complex<double> get_peierl_substitution(unsigned int i, unsigned int j);
        complex<double> get_peierl_substitution(const Vector2d& start, const Vector2d& finish) ;
        /**
         * \brief returns the number of clusters in one magnetic cell
         * @return
         */
        unsigned int number_clusters_in_magnetic_cell() const {return std::max(int(q_/M_),1);}

        void update_v_magnetic(const VectorXd & k);
        void update_current_cluster(unsigned int current);


        private:

        unsigned int M_; /**< Width of the cluster (number of times the basis is repeated in the a1_ direction) */
        unsigned int N_; /**< Height of the cluster (number of times the basis is repeated in the a2_ direction) */

        unsigned long filling_; /**< Number of electrons of spin up (=spin down) in one cluster */
        vector<unsigned long> states_; /**< list of the states having filling_ electrons */
        unsigned long nb_states_; /**< number of states */
        VectorXcd ground_state_; /**< ground state of the system */
        double E0_;/**< energy of the ground state */
        MatrixXcd V_;/**< Perturbation of the system (intercluster terms) */
        MatrixXcd V_magnetic_;/**< Perturbation of the system (inter magnetic cells terms) */

        vector<vector<vector<TridiagonalMatrix>>> intra_G_parameters_e_; /**< parameters allowing to compute \f$ G^' (\omega \f$ for each cluster*/
        vector<vector<vector<TridiagonalMatrix>>> intra_G_parameters_h_; /**< parameters allowing to compute \f$ G^' (\omega \f$ for each cluster*/

        MatrixXcd G_;/**< Full Green's function */
        MatrixXcd G_magnetic_; /**< Green's function for the magnetic cell (without intercluster coupling in the direction of the magnetic potential) */
        MatrixXcd G_intern_;/**< Intracluster Green's function */

        Vector2d corner_position_;/**< absolute coordinate of the first point of the cluster (needed for magnetic field computation) */

        unsigned int magnetic_cluster_; /**< number of the magnetic clusters, in range[0,q_/M_) */
        int p_,q_; /** magnetic field */
};


#endif //INC_2DTOOLS_CPTSOLVER_HH
