#pragma once

/* Class definition for Intersection */
#ifndef INTERSECTION_H
#define INTERSECTION_H
class Shape;
//#include "Header.h"
#include <vector>
class Intersection {
public:
	Intersection(){ primative = NULL; };
	Intersection(std::vector<Shape*>&, Ray&);
	Shape* primative;
	vec3 point;	//交点
	vec3 sourceDirection;	//从交点出发的方向，单位向量//sourceDirection = -ray.direction;

	double travel;
};

#endif