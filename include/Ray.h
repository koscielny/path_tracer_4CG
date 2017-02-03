#pragma once


/* Class description for Ray */
#ifndef RAY_H
#define RAY_H
#ifndef GLMH_H
#endif
#include "glmHeader.h"

class Ray {
public:
	Ray(){};
	Ray(const vec3&, const vec3&);

	vec3 getPoint(double);
	vec3 origin;
	vec3 direction;
};


#endif