//
// Created by tomlucas on 05.11.20.
//

#ifndef ZAVI_BOULDERN_TRANSFORMIT_H
#define ZAVI_BOULDERN_TRANSFORMIT_H

#include <Eigen/Core>
template<typename Iterable,typename Unary>
std::vector<Eigen::Vector3d> transformToPoints(const Iterable & states, Unary transform,int stride=1 ){
    std::vector<Eigen::Vector3d> points;
    int i=rand()%stride;
    for(auto state: states){
        if (i++%stride==0)
        points.push_back(transform(state));
    }
    return points;
}

#endif //ZAVI_BOULDERN_TRANSFORMIT_H
