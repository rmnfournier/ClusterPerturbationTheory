/**
 * \file CPTSolver.hh
 * \author Romain Fournier
 * \date 25.10.18
 **/

#include "CPTSolver.hh"
#include <fstream>

/**
 * @tparam T
 * @param a starting point
 * @param b ending point
 * @param N number of steps
 * @return equivalent to python's linspace
 */
template <typename T>
vector<T> linspace(T a, T b, int N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

CPTSolver::CPTSolver(const string &simulation_name_, const string &filename_, const Material &mat_, int k_resolution,
                     double omega_max, double omega_min, double omega_resolution_,
                     const vector<string> &directions_name_, const vector<Vector2d> &starting_k_points_,
                     const vector<Vector2d> &finishing_k_points_, unsigned int M_, unsigned int N_, unsigned long filling_,int p,int q,
                     unsigned int  magnetic_cluster) : Solver(simulation_name_, filename_,
                                                                                                                                     mat_, k_resolution,
                                                                                                                                     omega_max, omega_min,
                                                                                                                                     omega_resolution_,
                                                                                                                                     directions_name_,
                                                                                                                                     starting_k_points_,
                                                                                                                                     finishing_k_points_), M_(M_),
                                                                                                                              N_(N_),filling_(filling_),magnetic_cluster_(magnetic_cluster),p_(p),q_(q) {
    V_.setZero(get_occupation_states_dim(),get_occupation_states_dim());
    G_intern_=V_;
    G_magnetic_.setZero(get_occupation_states_dim()*number_clusters_in_magnetic_cell(),get_occupation_states_dim()*number_clusters_in_magnetic_cell());
    V_magnetic_=G_magnetic_;
    G_=G_magnetic_;

    vector<vector<vector<TridiagonalMatrix>>> temp (number_clusters_in_magnetic_cell(),vector<vector<TridiagonalMatrix>>(get_occupation_states_dim(),vector<TridiagonalMatrix>(get_occupation_states_dim())));
    intra_G_parameters_e_=std::move(temp);
    intra_G_parameters_h_=intra_G_parameters_e_;
    init_states();
    nb_states_=states_.size();
    ground_state_=VectorXcd::Random(nb_states_);
    corner_position_=mat_.get_lattice().getA1_()*magnetic_cluster_*M_;
}
void CPTSolver::get_Green_function()  {
    // For each cluster in the magnetic cell
    for (unsigned int c(0); c<number_clusters_in_magnetic_cell();c++) {
        //Initialize the process by choosing the appropriate cluster
        update_current_cluster(c);
        // Compute the Ground state for the current filling
        compute_ground_state();
        //cout<<"Ground state : "<<endl<<ground_state_<<endl;
        // Compute the intracluster Green's function parameter, that allows to compute G-1
        update_Green_intra_coefficients();
        //cout<<"Intra Green's coefficients ok"<<endl;
    }

    // For each chosen direction
    for(int ii(0);ii<directions_name_.size();++ii) {
        // Compute the dk
        Vector2d dk;
        for (int jj(0); jj < 2; jj++) {
            dk(jj) = (finishing_k_points_[ii](jj) - starting_k_points_[ii](jj)) / (k_resolution + 0.0);
        }

        // For each dk in the direction
        Vector2d k = starting_k_points_[ii];
        ofstream my_file;
        my_file.open(filename_+"_green_"+directions_name_[ii]+".csv");
        //cout<<filename_+"_green_"+directions_name_[ii]+".csv"<<endl;
        //cout<<"starting direction "<<directions_name_[ii]<<endl;
        for (int kk(0); kk <= k_resolution; kk++) {
            //cout<<kk/(k_resolution+0.0)*100<<"%"<<endl;
            vector<double> energy_range = linspace(omega_min, omega_max, omega_resolution_);
            // Now we must add the perturbation between cluster in the same magnetic cell
            if(p_!=0)
                update_v_magnetic(k);
            for (const auto &omega : energy_range) {
                //Here, we compute the block diagonal matrix corresponding to the intracluster green's function+ hopping between equivalent cluster
                for (unsigned int c(0); c < number_clusters_in_magnetic_cell(); c++) {
                    update_current_cluster(c);
                    compute_G_intern(omega);
                    update_perturbation(k);
                    cout<<"G intern -1 : "<<endl<<G_intern_.inverse()<<endl<<" V : "<<endl<<V_<<endl<<"vmagn = "<<endl<<V_magnetic_<<endl;
                    G_magnetic_.block(get_occupation_states_dim() * magnetic_cluster_,
                                      get_occupation_states_dim() * magnetic_cluster_, get_occupation_states_dim(),
                                      get_occupation_states_dim()) = ((G_intern_).inverse() - V_);
                }
                    // Finally we can get the green's function
              /*  cout<<"************ G_Magnetic "<<endl<<G_magnetic_<<endl;
                cout<<"************ V_Magnetic "<<endl<<V_magnetic_<<endl;*/

                G_ = (G_magnetic_- V_magnetic_).inverse();
                    if(p_==0){
                       // cout<<-1/M_PI*imag(compute_G_periodic(k).trace())<<",";
                        my_file<<-1/M_PI*imag(compute_G_periodic(k).trace())<<",";
                    }
                    else{
                       // cout<<-1/M_PI*imag(compute_G_periodic(k).trace())<<",";
                        my_file<<-1/M_PI*imag(G_.trace())<<",";
                    }

            }
            my_file<<endl;
            //Increment k
            k += dk;
        }
        my_file.close();
    }
}


void CPTSolver::update_v_magnetic(const VectorXd& k){
    V_magnetic_.setZero();
    int tmp(magnetic_cluster_);
    for(unsigned int cluster (0);cluster<number_clusters_in_magnetic_cell();cluster++){
        update_current_cluster(cluster);
        for(unsigned int ii(0);ii<get_occupation_states_dim();ii++){
            int i_no_spin_cluster_offset = ((cluster*get_occupation_states_dim()+ii)/get_number_of_sites())/2;
            i_no_spin_cluster_offset*=get_number_of_sites();
            int i_in_cluster = (cluster*get_occupation_states_dim()+ii)%get_number_of_sites();
            int spin_offset = (((cluster*get_occupation_states_dim()+ii)/get_number_of_sites())%2)*get_number_of_sites();
            for(const auto& j : get_translation_for_neighbors_from_other_cluster(i_no_spin_cluster_offset+i_in_cluster,true)){
                // the tricks is to represent the state in the correct way for the get_absolute_position(ii), which does not take the spin into account
                int j_no_spin_cluster_offset = ((j.first+spin_offset)/get_number_of_sites())/2;
                j_no_spin_cluster_offset*=get_number_of_sites();
                int j_in_cluster = (j.first+spin_offset)%get_number_of_sites();

              /*  cout<<"**************"<<endl;
                cout<<"Considering i="<<cluster*get_occupation_states_dim()+ii<<", located at "<<get_absolute_position(i_no_spin_cluster_offset+i_in_cluster).transpose()<<endl;
                cout<<"Neighbor at j="<<j.first+spin_offset<<", cooresponding to clusteroffset "<<j_no_spin_cluster_offset<<" and position "<<j_in_cluster<<endl;
                cout<<" reaching j with translation "<<j.second.transpose()<<endl;
                 cout<<"Vector from i to j "<<(j.second+get_absolute_position(j_no_spin_cluster_offset+j_in_cluster)-get_absolute_position(i_no_spin_cluster_offset+i_in_cluster)).transpose()<<endl;
                 cout<<"location of the neighbor is "<<(j.second+get_absolute_position(j_no_spin_cluster_offset+j_in_cluster)).transpose()<<endl;
                cout<<" translation to reach the magnetic cell "<<j.second.transpose()<<endl;
                cout<<" peierl term = "<<get_peierl_substitution(get_absolute_position(i_no_spin_cluster_offset+i_in_cluster),j.second+get_absolute_position(j_no_spin_cluster_offset+j_in_cluster))<<endl;*/
                V_magnetic_(cluster*get_occupation_states_dim()+ii,j.first+spin_offset)+=-mat_.get_t(0)*compl_exp(k.dot(j.second))*get_peierl_substitution(get_absolute_position(i_no_spin_cluster_offset+i_in_cluster),j.second+get_absolute_position(j_no_spin_cluster_offset+j_in_cluster));
            }
        }
    }
    update_current_cluster(tmp);
}
void CPTSolver::update_current_cluster(unsigned int current){
    magnetic_cluster_=current;
    corner_position_=mat_.get_lattice().getA1_()*current;
}

void CPTSolver::update_Green_intra_coefficients() {
    //Computation of G'_ii,jj
    for(int ii(0);ii<get_occupation_states_dim();++ii){
            for(int jj(0);jj<get_occupation_states_dim();++jj){
                    // Electron part
                    vector<double>a,b;
                    a.clear();
                    b.clear();
                    VectorXcd phi_0(creator(ii)+creator(jj));
                    if(ii==jj) phi_0/=2.;
                if(phi_0.norm()>1e-8){
                        get_lanczos_coefficients_from_states(phi_0,a,b,1e-8,true,false,true);
                    }
                    TridiagonalMatrix current;
                    current.push_back(a);
                    current.push_back(b);
                    intra_G_parameters_e_[magnetic_cluster_][ii][jj]=current;

                    // Hole part
                    a.clear();
                    b.clear();
                    phi_0.setZero();
                    phi_0=(destructor(ii)+destructor(jj));
                    if(ii==jj) phi_0/=2.;

                if(phi_0.norm()>1e-8) {
                    get_lanczos_coefficients_from_states(phi_0,a,b,1e-8,true,false,true);
                }
                    current.clear();
                    current.push_back(a);
                    current.push_back(b);
                    intra_G_parameters_h_[magnetic_cluster_][ii][jj]=current;
            }
        }
}

void CPTSolver::update_perturbation(const Vector2d &k) {
    V_.setZero(get_occupation_states_dim(),get_occupation_states_dim());
    // Kinetic part
    for ( long ii : get_sites_list() ){
        // get all "outside" neighbors of ii in the form (site,vector translation)
        for(const auto& jj : get_translation_for_neighbors_from_other_cluster(ii)){
            complex<double> hopping(-mat_.get_t()*compl_exp(-k.dot(jj.second)));
            //Spin up
            V_(ii,jj.first)+=hopping;
            //Spin down
            V_(ii+get_number_of_sites() ,jj.first+get_number_of_sites())+=hopping;
         }
    }
    // Potential part
    update_outside_potential(k);
}

unsigned long CPTSolver::get_Fock_dim() const{
    // 2 ^ the number of sites * 2, because of the spin
    return 2<<(get_occupation_states_dim()-1);
}

vector<int> CPTSolver::get_sites_list() const{
    vector<int>list(get_number_of_sites());
    std::iota(std::begin(list),std::end(list),0);
    return list;
}

unsigned long CPTSolver::get_number_of_sites() const{
    return M_*N_*mat_.get_number_elements_in_basis();
}

void CPTSolver::compute_G_intern(double omega){
    // We start with the diagonal elements, because we use them to compute the off diagonal elements
    for(int ii(0);ii<get_occupation_states_dim();++ii) {
        complex<double> tmp = compute_G_intern_diag_element(ii,ii,omega);
        G_intern_(ii,ii)=compute_G_intern_diag_element(ii,ii,omega);
    }
    // Then we switch to the non-diagonal elements
    for(int ii(0);ii<get_occupation_states_dim();++ii){
        for(int jj(0);jj<get_occupation_states_dim();++jj){
            if(ii!=jj){
                complex<double> tmp = compute_G_intern_offdiag_element(ii,jj,omega);
                G_intern_(ii,jj)=compute_G_intern_offdiag_element(ii,jj,omega);
            }
        }
    }
}

complex<double> CPTSolver::compute_G_intern_diag_element_helper(int i,int j,complex<double> omega,bool electron){
    //spread parameter
    complex<double>etha(0,0.01);

    complex<double> diag_e,diag_h; //holes and electrons diagonal contributions
    diag_e=0;
    diag_h=0;
    //number of terms in the continuous fraction

     long size_e=intra_G_parameters_e_[magnetic_cluster_][i][j][0].size();
     long size_h=intra_G_parameters_h_[magnetic_cluster_][i][j][0].size();

    //Compute the fraction recusively (electrons)

    for( long jj(size_e-1);jj>0;jj--){
        diag_e= intra_G_parameters_e_[magnetic_cluster_][i][j][B][jj-1]*intra_G_parameters_e_[magnetic_cluster_][i][j][B][jj-1]/(omega-etha+E0_-intra_G_parameters_e_[magnetic_cluster_][i][j][A][jj]-diag_e);
    }
    for( long jj(size_h-1);jj>0;jj--){
        diag_h= intra_G_parameters_h_[magnetic_cluster_][i][j][B][jj-1]*intra_G_parameters_h_[magnetic_cluster_][i][j][B][jj-1]/(omega+etha-E0_+intra_G_parameters_h_[magnetic_cluster_][i][j][A][jj]-diag_h);
    }
    //coefficient
    double b0_2(0);
    if(i!=j){
        b0_2=pow((creator(i)+creator(j)).norm(),2);
    }else{
        b0_2=pow((creator(i)).norm(),2);
    }

    if(size_e>=1) diag_e=(b0_2/(omega-etha+E0_-intra_G_parameters_e_[magnetic_cluster_][i][j][A][0]-diag_e));

    //coefficient
    if(i!=j){
        b0_2=pow((destructor(i)+destructor(j)).norm(),2);
    }else{
        b0_2=pow((destructor(i)).norm(),2);
    }
    if(size_h>=1) diag_h=(b0_2/(omega+etha-E0_+intra_G_parameters_h_[magnetic_cluster_][i][j][A][0]-diag_h));

    if(electron)
        return diag_e;
    else
        return diag_h;
}

complex<double> CPTSolver::compute_G_intern_diag_element(int i,int j,complex<double> omega){
    return compute_G_intern_diag_element_helper(i,j,omega,true)+compute_G_intern_diag_element_helper(i,j,omega,false);
}

complex<double> CPTSolver::compute_G_intern_offdiag_element(int i,int j, complex<double> omega){
    complex<double> sym_e(compute_G_intern_diag_element_helper(i,j,omega,true)),sym_h(compute_G_intern_diag_element_helper(i,j,omega,false)); // contains the term for the electron and the one for the hole
    return 0.5*(sym_e-compute_G_intern_diag_element_helper(i,i,omega,true)-compute_G_intern_diag_element_helper(j,j,omega,true) +sym_h-compute_G_intern_diag_element_helper(i,i,omega,false)-compute_G_intern_diag_element_helper(j,j,omega,false));
}

void CPTSolver::update_outside_potential(const Vector2d& k){

}

VectorXcd CPTSolver::creator(unsigned long i){
    VectorXcd phi=VectorXcd::Zero(get_Fock_dim());
    for( long ii(0);ii<get_Fock_dim();++ii){
        if(!CHECK_BIT(ii,i)){
            unsigned long target = toggled(ii,i);
            double invert=(1);
            for(unsigned int jj(i+1);jj<get_occupation_states_dim();jj++){
                if(CHECK_BIT(ii,jj))
                    invert*=-1;
            }
            phi(target)+=ground_state_(ii)*invert;
        }
    }
    //cout<<phi.transpose()<<endl;
    return phi;
}

VectorXcd CPTSolver::destructor(unsigned long i){
    VectorXcd phi=VectorXcd::Zero(get_Fock_dim());
    //cout<<"desetructin at "<<i<<endl;
    // go through the whole Fock Space
    for( long ii(0);ii<get_Fock_dim();++ii){
        //If the ith site is occupied, we can destroy it
        if(CHECK_BIT(ii,i)){
            //find at which dimension we end up after destroying the ith bit of the current state
            unsigned long target = toggled(ii,i);
            // Since c and cdag anticommute, we need to take care of the sign
            double invert=(1);
            for(unsigned int jj(i+1);jj<get_occupation_states_dim();jj++){
                if(CHECK_BIT(ii,jj))
                    invert*=-1;
            }
           // cout<<"for ii="<<ii<<" going to "<<target<<" we have sign "<<invert<<endl;
            phi(target)+=ground_state_(ii)*invert;
        }
    }
    //cout<<phi.transpose()<<endl;
    return phi;
}

void CPTSolver::get_lanczos_coefficients_from_states(const VectorXcd& phi_0,vector<double>& a,vector<double>& b,double tolerance,bool full_orthogonalization,bool reconstruct,bool early_stop){
    //Make sure the storage vectors are empty
    a.clear();
    b.clear();

    //Lanczos vectors for spanning the Krilov space
    VectorXcd q=phi_0;
    VectorXcd q_last=VectorXcd::Zero(get_Fock_dim(),1);
    // List of the q vectors for full_orthogonalization
    vector<VectorXcd> listq;


    q/=q.norm();
    VectorXcd gs;
    if(reconstruct){
        gs.setZero(get_Fock_dim());
    }
    double beta(0),alpha(0);

    int jj(0);
    for ( jj=0; jj<get_Fock_dim() ;++jj){
        //cout<<jj<<" beta = "<<beta<<endl;
        // save the vector if we need it later for the reorthogonalization
        if(full_orthogonalization || !early_stop){
            listq.emplace_back(q);
        }
        //compute ground_state
        if(reconstruct){
            for(size_t kk(0);kk<working_dim();kk++){
                gs(states_[kk])+=q(states_[kk])*ground_state_(jj);
            }
        }
        //Compute H|q>
        VectorXcd z = apply_intra_H(q);
        //Compute the next diagonal element
        alpha=real(q.dot(z));

        // Now we can look for the next vector
        if(full_orthogonalization){
            //for some computational reason, Gram-Schmidt must be applied twice
            for (int mystery(0);mystery<2;mystery++){
                for(int kk=0;kk<=jj;kk++){
                    z=z-listq[kk].dot(z)*listq[kk];

                }
            }
        }else {
            z=z-alpha*q-beta*q_last;
        }
        //Update beta
        beta=z.norm();
        //saves the elements
        a.push_back(alpha);
        b.push_back(beta);

        if(beta>tolerance){
            q_last=q;
            q=z/beta;
        }
        else if((jj!= get_Fock_dim()-1) && !early_stop){
            //random restart (we were in an invariant subspace, we need to go out of it)
            do{
                q=VectorXcd::Random(get_Fock_dim(),1);
                for (int mystery(0);mystery<2;mystery++){
                    for(int kk=0;kk<=jj;kk++){
                        q=q-listq[kk].dot(q)*listq[kk];
                    }
                }
            }while(q.norm()<tolerance);
            q=q/q.norm();
        }
        else {
            break;
        }
    }
    //drop the last element, since we do not need it
    b.pop_back();
    if(reconstruct)
        ground_state_=gs;
}

VectorXcd CPTSolver::apply_intra_H(const VectorXcd& q){
  //  cout<<"Applying H on "<<endl<<q.transpose()<<endl;
    VectorXcd results;
    results.setZero(get_Fock_dim());
    // Kinetic part
    // todo : add other neighbors, for now there are only first nearest neighbors

    for(unsigned long ii(0); ii<get_Fock_dim();++ii){
        for(unsigned long j : get_hopping_from_i_states(ii)){
            //We need to be carefull about an eventual change of sign
            pair<unsigned long,unsigned long>xchange=get_diff_index(ii,j);
            //we need to count the number of occupied sites between the hopping terms
            double sign(1);
            for(unsigned long i(min(xchange.first,xchange.second)+1);i<max(xchange.second,xchange.first);i++){
                if(CHECK_BIT(ii,i)){
                    sign*=-1;
                }
            }
            results(ii)+=-mat_.get_t()*q(j)*sign*get_peierl_substitution(xchange.first,xchange.second);
        }
    }
   // cout<<"gives "<<endl<<results.transpose()<<endl;
    // Potential part
    int states_spin=sqrt(get_Fock_dim());
    for(unsigned long ii(0); ii<get_Fock_dim();++ii) {
        // First we solve exactly for the cluster
        unsigned long left = ii % states_spin;
        unsigned long right = ii / states_spin;
        unsigned long doubly = bit_sum(left & right);
        results(ii) += mat_.get_u() * doubly * q(ii);

        //Then, add the mean field term
/*        double mean_field(0.);
        for(unsigned int i(0);i< get_number_of_sites();++i){
            mean_field+= CHECK_BIT(left,i)*mean_occupation_down_[i]+CHECK_BIT(right,i)*mean_occupation_up_[i]-mean_occupation_up_[i]*mean_occupation_down_[i];
        }
        results(ii)+=mat_.get_u()*mean_field;*/
    }
    return results;
}

vector<unsigned long> CPTSolver::get_hopping_from_i_states(unsigned long i){
    vector<unsigned long> list;
    int separator(sqrt(get_Fock_dim()));
    for(unsigned long current(0);current<get_Fock_dim();current++){
        //we must first of all make sure that we have the same occupation for the two subspace (up and down)
        if( (bit_sum(current%separator)== bit_sum(i%separator) ) && (bit_sum(current/separator)== bit_sum(i/separator) )){
            //if spin up differ at 2 sites but down are the same or otherwise
            if( (number_of_different_bits_for_same_occupation(current%separator,i%separator)==2 &&number_of_different_bits_for_same_occupation(current/separator,i/separator)==0) || (number_of_different_bits_for_same_occupation(current%separator,i%separator)==0 &&number_of_different_bits_for_same_occupation(current/separator,i/separator)==2) ){
                std::pair<unsigned long,unsigned long> coef = get_diff_index(i,current);
                if(are_connected(coef.first%get_number_of_sites() ,coef.second % get_number_of_sites())){
                    list.push_back(current);
                }
            }
        }

    }

    return list;
}

bool CPTSolver::are_connected(unsigned long i, unsigned long j){
    // Convert the number of the site into the m,n,b representation
    ClusterPosition pos_i(get_cluster_position(i)),pos_j(get_cluster_position(j));
  //  cout<<pos_i<<endl<<" and "<<endl<<pos_j<<endl<<"gives "<<endl;
    // for each neighbor of i
    for(const auto& neighbor_i : mat_.get_neighborhood(pos_i(2))){
   //     cout<<"for "<<neighbor_i.transpose()<<endl ;
        //Return true if j is a neighbor of i
        if( (pos_i(0)+neighbor_i(0) == pos_j(0)) && (pos_i(1)+neighbor_i(1) == pos_j(1)) && (neighbor_i(2)==pos_j(2)) ) {
            return true;
        }
    }
    // if we didn't find j in the neighborhood of i, return false
    return false;
}

ClusterPosition CPTSolver::get_cluster_position(unsigned long i){
    //Initialize the position
    ClusterPosition pos;
    pos.setZero();

    //start with the basis element
    pos(2)=i % mat_.get_number_elements_in_basis();
    // Then the a1 translation m
    pos(0)= (i / mat_.get_number_elements_in_basis())%M_;
    //Finally the a2 translation n
    pos(1)= (i/mat_.get_number_elements_in_basis()/M_);
    return pos;
}

vector<std::pair<unsigned long,Vector2d>> CPTSolver::get_translation_for_neighbors_from_other_cluster(unsigned long i,bool magnetic_cluster){
    // Prepare the list of neighbors-translation pair that will be returned
    vector<std::pair<unsigned long,Vector2d>> list;

        // Get the representation in Cluster coord. (a1,a2,basis)
    ClusterPosition pos_i(get_cluster_position(i));
    if(magnetic_cluster){
        pos_i(1)-=magnetic_cluster_;
    }
    // Now, loop over each neighbor of an element located at basis pos_i(2)
    for(const auto& neighbor_i : mat_.get_neighborhood(pos_i(2))){
            // prepare markers for later
            bool outside(false); //True if we leave the cluster
            bool flag(false); // True if we leave via A1

            Vector2d translation;
            translation.setZero();

            //find the position of the neighbor position_j
            ClusterPosition pos_j=pos_i;
            pos_j(0)+=neighbor_i(0);
            pos_j(1)+=neighbor_i(1);
            pos_j(2)=neighbor_i(2); //the basis is not additive, we have to take the one of the neighbor

            // Check if we left via A1
            if(pos_j(0)>=int(M_)){
                flag=true;
                if(p_==0){
                    outside=true;
                    pos_j(0)%=M_;
                    translation+=mat_.get_lattice().getA1_()*M_;
                }
                //If we are in a simulation with a magnetic field, we only have a fourier term if we jump from the last cluster
                else if(magnetic_cluster){
                    outside=true;
                    //if it is the last cell
                    if( magnetic_cluster_==(number_clusters_in_magnetic_cell()-1)){
                        translation+=mat_.get_lattice().getA1_()*M_*number_clusters_in_magnetic_cell();
                        pos_j(0)%=M_;
                    }
                    //othervise we put it in the next cluster, without translation term since we are in the same magnetic cell
                    else{
                        pos_j(0)=pos_j(0)%M_+2*(magnetic_cluster_+1)*M_*N_;
                    }
                }

            }else if(pos_j(0)<0){
                flag=true;
                if(p_==0){
                    outside=true;
                    pos_j(0)+=M_;
                    translation-=mat_.get_lattice().getA1_()*M_;
                }
                else if(magnetic_cluster){
                    outside=true;
                    //if it is the first cluster, it hopes to the last one in the fourier space
                    if( magnetic_cluster_==0){
                        pos_j(0)+=M_+2*M_*N_*(number_clusters_in_magnetic_cell()-1);
                        translation-=mat_.get_lattice().getA1_()*M_*number_clusters_in_magnetic_cell();
                    }
                    //othervise just join the previous one
                    else{
                        pos_j(0)+=M_+2*M_*N_*(magnetic_cluster_-1);
                    }
                }
            }

            // We consider now the case when we leave via A2
            if(pos_j(1)>=int(N_)){
                //in case of the magnetic_cluster, we only consider the case where we jump to the other cluster previously
                if (p_==0){
                    outside=true;
                    pos_j(1)%=N_;
                    translation+=mat_.get_lattice().getA2_()*N_;
                }
                else if( (magnetic_cluster) || ( (!magnetic_cluster) &&(!flag)) ){
                    outside=true;
                    pos_j(1)%=N_;
                    translation+=mat_.get_lattice().getA2_()*N_;
                }
            }else if(pos_j(1)<0){
                if(p_==0){
                    outside=true;
                    pos_j(1)+=N_;
                    translation-=mat_.get_lattice().getA2_()*N_;
                }
                else if((magnetic_cluster) || ( (!magnetic_cluster) &&(!flag))){
                    outside=true;
                    pos_j(1)+=N_;
                    translation-=mat_.get_lattice().getA2_()*N_;
                }
            }
            // If we left the cluster, we add the info to the list
            if(outside){
                int j(get_site_from_position(pos_j));
                //If we need to go back to the other spin
                list.emplace_back(std::make_pair(j,translation));
            }
    }
    return list;
}

int CPTSolver::get_site_from_position (ClusterPosition pos){
    return pos(2)+mat_.get_number_elements_in_basis()*pos(0)+M_*mat_.get_number_elements_in_basis()*pos(1);
}

void CPTSolver::add_perturbative_hopping_between(int i,int j,complex<double> hopping){
    //We iterate over all states
    for (unsigned long ii(0);ii<get_Fock_dim();++ii){
        //if the integer does not contain a 1 on the ith and the jth place
        if( (!CHECK_BIT(ii,i)) &&(!CHECK_BIT(ii,j)) ){
            //we use this number by flipping the jth bit to one for the row and the ith for the column
            unsigned long line(toggled(ii,j)) ,column(toggled(ii,i));
            V_(line,column)=hopping;
        }
    }
}

MatrixXcd CPTSolver::compute_G_periodic(const Vector2d& k) {
    /* We need to perform a residual fourier transform */
    /* Now, G_ is a sxMxN square matrix. we have to reduce it into a sxs matrix */
    // Initialize the matrix
    MatrixXcd rft;
    rft.setZero(mat_.get_number_elements_in_basis(),mat_.get_number_elements_in_basis() );
    //loop over all position in the cluster
    for(unsigned int a : get_sites_list()){
        for(unsigned int b : get_sites_list()){
            ClusterPosition ca(get_cluster_position(a)), cb(get_cluster_position(b));
         //   cout<<"ca "<<endl<<ca<<endl;
         //   cout<<"cb "<<endl<<cb<<endl;

            //We only need to perform a F.T. between the different repetitions of the basis
            Vector2d ra=get_vector_from_cluster_position(ca);
            Vector2d rb=get_vector_from_cluster_position(cb);
           // cout<<"ra "<<endl<<ra<<endl;
           // cout<<"rb "<<endl<<rb<<endl;



            rft(ca(2),cb(2))+= compl_exp(-k.dot(ra-rb))*G_(ca(0)*mat_.get_number_elements_in_basis()+ca(1)*M_*mat_.get_number_elements_in_basis()+ca(2) , cb(0)*mat_.get_number_elements_in_basis()+cb(1)*M_*mat_.get_number_elements_in_basis()+cb(2));
        }
    }
//    cout<<G_<<endl;
  //  cout<<" vs "<<endl<<rft<<endl;
    return rft/N_/M_;
    //return G_;
}

Vector2d CPTSolver::get_vector_from_cluster_position(const ClusterPosition& cp) const {
    return cp(0)*mat_.get_lattice().getA1_()+cp(1)*mat_.get_lattice().getA2_();
}

void CPTSolver::init_states() {
    // states with spin up and states with spin down must have a number of one equal to filling
    int separator(    2<<(get_occupation_states_dim()/2-1));
    for (unsigned long ii(0);ii<get_Fock_dim();++ii){
        //Check both side
        if(bit_sum(ii %separator)==filling_ && bit_sum(ii/separator)==filling_){
            states_.push_back(ii);
        }
    }
}

void CPTSolver::compute_ground_state() {
    vector<double> a,b;
    a.clear();b.clear();

    VectorXcd phi_0;
    phi_0.setZero(get_Fock_dim());
    for(unsigned int ii=0;ii<nb_states_;ii++){
        phi_0(states_[0])=1;
    }
    phi_0/=phi_0.norm();
    //cout<<"Starting lanczos for ground state "<<endl;
    get_lanczos_coefficients_from_states(phi_0,a,b,1e-8,false,false);
   // cout<<"coefficients found"<<endl;

    SelfAdjointEigenSolver<MatrixXcd> solver;
    VectorXd diag;
    VectorXd subdiag;

    diag.setZero(a.size());
    subdiag.setZero(b.size());

    for(int ii=0;ii<b.size();ii++){
        diag(ii)=a[ii];
        subdiag(ii)=b[ii];
    }
    diag(a.size()-1)=a[a.size()-1];
    solver.computeFromTridiagonal(diag,subdiag);
    MatrixXcd ev = solver.eigenvectors();
    ground_state_=ev.col(0);
    get_lanczos_coefficients_from_states(phi_0,a,b,1e-8,false,true);
    E0_=solver.eigenvalues()(0);
    //cout<<"energy and eigenstates found"<<endl;
}

Vector2d CPTSolver::get_absolute_position(unsigned int site){
    int offset=site/get_number_of_sites();
    site%=get_number_of_sites();

    ClusterPosition cp(get_cluster_position(site));
    return (offset+cp(0))*mat_.get_lattice().getA1_()+cp(1)*mat_.get_lattice().getA2_()+mat_.get_lattice().getBasis(cp(2));
}

complex<double> CPTSolver::get_peierl_substitution(unsigned int i, unsigned int j){
    //Compute the direction towards which we integrate
    Vector2d starting_point (get_absolute_position(i));
    Vector2d finishing_point (get_absolute_position(j));
    return get_peierl_substitution(starting_point,finishing_point);
}

complex<double> CPTSolver::get_peierl_substitution(const Vector2d & start,const Vector2d & finish) {
    // dx and dy
    double dy=(finish(1)-start(1));
    double dx=(finish(0)-start(0));

    return compl_exp(2*M_PI* p_/(q_+0.0)*(dy)*dx/mat_.get_lattice().get_surface()*(start(0)/dx+0.5) );

}