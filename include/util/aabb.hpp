#pragma once

// from https://github.com/brandonpelfrey/Fast-BVH

#ifndef BBox_h
#define BBox_h

//#include "Ray.h"

//#ifndef TRANSFORM_H
//#endif
class Ray;
class Vector;
#include "../Ray.h"
#include "vec3.hpp"
#include "glmHeader.h"
#include <stdint.h>
struct AABB2 {
	Vector min, max, extent;
	AABB2() { }
	AABB2(const Vector& min, const Vector& max);
	AABB2( vec3& min,  vec3& max);
	AABB2(const Vector& p);

	bool intersect(const Ray& ray, float *tnear, float *tfar) const;
 void expandToInclude(const Vector& p);
 void expandToInclude(const AABB2& b);
 uint32_t maxDimension() const;
 float surfaceArea() const;
};

#endif
