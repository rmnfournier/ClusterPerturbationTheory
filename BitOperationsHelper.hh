/**
 * \file BitOperationsHelper.hh
 * \author Romain Fournier
 * \date 19.10.18
 **/

#ifndef INC_2DCPT_BITOPERATIONSHELPER_HH
#define INC_2DCPT_BITOPERATIONSHELPER_HH

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

#include <iostream>
#include <eigen3/Eigen/Dense>
/**
 * \brief Compute the sum of all bits of unsigned integer
 * @param i
 * @return sum of all bits of unsigned integer i
 */
unsigned long bit_sum(unsigned long i);
/**
 * \brief if i and j have the same number of 1 in their binary representation, returns the sum of the positions at which their bit differ
 * @param i first state
 * @param j second state
 * @return 0 if i and j do not have the same number of atoms, number of differences otherwise
 */
unsigned long number_of_different_bits_for_same_occupation(unsigned long i, unsigned long j);

/**
 * \brief Compute the first two index at which the bits from i differ from the one of j
 * @param i first state
 * @param j second state
 * @return the first two index at which the bits from i differ from the one of j
 */
std::pair<unsigned long,unsigned long> get_diff_index(unsigned long i, unsigned long j);
/**
 * Flips the pos th bit ( starting counting at 0) of i
 * @param i
 * @param pos
 * @return i with a bit flipped at position pos
 */
unsigned long toggled(unsigned long i,int pos);

Eigen::VectorXd int_to_state(unsigned long l, unsigned int size);
/**
 * Count the number of occupied sites between positions pos1 and pos2
 * @param state
 * @param pos1
 * @param pos2
 * @return
 */
unsigned long occupancy_between_position(unsigned long state, unsigned long pos1, unsigned long pos2);
#endif //INC_2DCPT_BITOPERATIONSHELPER_HH
