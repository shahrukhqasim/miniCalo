//
// Created by Shah Rukh Qasim on 10/8/21.
//
#include "json/json.h"
#include <fstream>
#include <iostream>
#include <CLHEP/Vector/ThreeVector.h>


int main(int argc, char*argv[]) {

    CLHEP::Hep3Vector v1(1500,0,3207);

    std::cout<<"Eta: "<<v1.eta()<<std::endl;
    std::cout<<"Phi: "<<v1.phi()<<std::endl;
    std::cout<<"Theta: "<<v1.theta()<<std::endl;
    std::cout<<"Z: "<<v1.z()<<std::endl;


    return EXIT_SUCCESS;
}