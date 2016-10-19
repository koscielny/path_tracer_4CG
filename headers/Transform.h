// Transform header file to define the interface. 
// The class is all static for simplicity
// You need to implement left, up and lookAt
// Rotate is a helper function

// Include the helper glm library, including matrix transform extensions
#pragma once
#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "glmHeader.h"


const double pi = 3.14159265358979 ; // For portability across platforms

class Transform	 
{
public:
	Transform();
	virtual ~Transform();
	static mat4 rotate(double degrees, const vec3& axis) ;
	static mat4 scale(double sx, double sy, double sz) ; 
	static mat4 translate(double tx, double ty, double tz);
};

#endif

