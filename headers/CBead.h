/*
 * CPolymers.h
 *
 *  Created on: Aug 9, 2013
 *      Author: jh
 */

#ifndef CBEAD_H_
#define CBEAD_H_

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "boost/random.hpp"
#include <Eigen/Dense>

using namespace std;


class CBead {//TODO include CPolymers into this here, by setting random sign!
private:
    //...

public:
    Eigen::Vector3d pos; 
    Eigen::Vector3d f_mob;   //store mobility and stochastic force
    Eigen::Vector3d f_sto;
    int sign = 1;
    
    CBead();
    CBead(int signx, Eigen::Vector3d posx){
        sign = signx;
        pos = posx;
    }
};



#endif /* CBEAD_H_ */
