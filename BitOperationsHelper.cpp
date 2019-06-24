/**
 * \file BitOperationsHelper.cpp
 * \author Romain Fournier
 * \date 19.10.18
 **/

#include "BitOperationsHelper.hh"
unsigned long number_of_different_bits_for_same_occupation(unsigned long i,unsigned long j){
    unsigned int v (i^j);

    if (bit_sum(i)==bit_sum(j))
        return bit_sum(v);
    else
        return 0;

}
unsigned long bit_sum(unsigned long v){
    unsigned int c ( ((v & 0xfff) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f);
    c += (((v & 0xfff000) >> 12) * 0x1001001001001ULL & 0x84210842108421ULL) %
         0x1f;
    c += ((v >> 24) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f;
    //sum for
    return c;
}

std::pair<unsigned long,unsigned long> get_diff_index(unsigned long i, unsigned long j){
    int counter(0);
    std::pair<unsigned int,unsigned int> coeffs=std::make_pair(-1,-1);
    for (unsigned int k(0);k< sizeof(i)*8;k++){
        int i_kth = i & (1<<k);
        int j_kth = j & (1<<k);
        if (i_kth!=j_kth){
            if(i_kth!=0){
                coeffs.first=k;
                counter++;
            }else{
                coeffs.second=k;
                counter++;
            }
        }
        if(counter==2) break;
    }
    return coeffs;
}
unsigned long toggled(unsigned long i,int pos){
    return (i ^= 1UL << pos);
}
Eigen::VectorXd int_to_state(unsigned long l, unsigned int size){
    Eigen::VectorXd state;
    state.setZero(size);
    for(unsigned int ii(0);ii<size;ii++){
        if(l&((1 << (ii )))){
            state(ii)=1;
        }
    }
    return state;
}
