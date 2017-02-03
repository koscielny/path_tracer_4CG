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
	vec3 point;	//����
	vec3 sourceDirection;	//�ӽ�������ķ��򣬵�λ����//sourceDirection = -ray.direction;

	double travel;
};

#endif