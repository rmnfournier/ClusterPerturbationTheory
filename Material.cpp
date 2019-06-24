//
// Created by romain on 25/10/18.
//

#include "Material.hh"
unsigned int Material::get_number_elements_in_basis() const {
    return lat_.get_number_elements_in_basis();
}
