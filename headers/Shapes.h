#pragma once

/* Class definition for a shape */
#ifndef SHAPE_H
#define SHAPE_H

//#include "Transform.h"
//#include <util/vec3.hpp>
//#include <util/aabb.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

#include <algorithm>
#include "glmHeader.h"
class Ray;
class Intersection;
class BVH;
class AABB2;
class Material;
class Vector;

#include "util\aabb.hpp"
#include "util\vec3.hpp"
class AABB {
public:
	AABB(){};
	AABB(vec3&,vec3&);
	vec3 aabbmax;
	vec3 aabbmin;
	double intersect(Ray&);
	vec3 center;
};


class Shape {
	public:
		virtual double intersect(Ray&)=0;
		virtual vec3 getNormal(const vec3&)=0;
		virtual vec3 getNormal(const Intersection&) = 0;//比getNormal(const vec3&)能防止法向弄反
		virtual double getSubtendedAngle(const vec3&)=0;
		virtual vec3 shade(const Intersection&, BVH*, vec3 s_norm, bool) = 0;
	//	virtual vec3 getTexture(vec3&){};
		//virtual vec3 gettrilog();


		/*! This method returns the axis-aligned bounding box of the primitive.
		\return The primitive's bounding box.
		\remark The bounding box need not be ideal, but the bounding volume hierarchy is more efficient if the
		returned bounding box tightly fits the primitive. The same bounding box must always be returned. */
		virtual AABB2 BoundingBox() = 0;

		/*! This method returns the centroid of the primitive.
		\return The primitive's centroid.
		\remark If the centroid is not well-defined, pass the best one and the bounding volume hierarchy will do its
		best to handle it. The same centroid must always be returned. */
		virtual Vector Centroid() = 0;

//		AABB aabb;
		//bool hasTexture;
		
		/* material properties */
		vec3 diffuse;
		vec3 specular;
		double shininess;
		vec3 emission;
		double indexofrefraction;
		double refractivity;
		//假如全归进material类
		Material *material;

		bool assertnorm;//如为NormTriangle则true；如仅为Triangle，norm不一定在正确的那一方向


};

class Sphere : public Shape {
	public:
		Sphere(mat4);
		double intersect(Ray&);
		vec3 getNormal(const vec3&);
		vec3 getNormal(const Intersection&);
		double getSubtendedAngle(const vec3&);
		vec3 shade(const Intersection&, BVH*, vec3 s_norm, bool);
		/* This function returns the bounding box of the sphere. */
		virtual AABB2 BoundingBox(){ return this->boundingBox; }

		/* This function returns the centroid of the sphere. */
		virtual Vector Centroid(){ return this->center; }

		AABB2 boundingBox;

		/* The sphere's center. */
		Vector center;

		mat4 mv;
		mat4 inv;	
};

class Triangle : public Shape {
	public:
		Triangle(){};
		Triangle(vec3,vec3,vec3);
		double intersect(Ray&);
		virtual vec3 getNormal(const vec3&);
		vec3 getNormal(const Intersection&);////比getNormal(const vec3&)能防止法向弄反

		double getSubtendedAngle(const vec3&);
		vec3 shade(const Intersection&, BVH*, vec3 s_norm, bool);
	//	vec3 getTexture(vec3&);
//		vec3 gettrilog(){ return ; };//辅助输出面片信息到log文件

		/* This function returns the bounding box of the triangle. */
		virtual AABB2 BoundingBox(){ return this->boundingBox; }

		/* This function returns the centroid of the triangle. */
		virtual Vector Centroid(){ return this->centroid; }


		AABB2 boundingBox;
		Vector centroid;
		vec3 p0;
		vec3 p1;
		vec3 p2;
		vec3 n0;
};

class NormTriangle : public Triangle {
public:
	NormTriangle(vec3,vec3,vec3,vec3,vec3,vec3);
	vec3 getNormal(const vec3&);
	vec3 getNormal(const Intersection&);

	vec3 n1;
	vec3 n2;

	/* This function returns the bounding box of the sphere. */
//	virtual AABB2 BoundingBox();

	/* This function returns the centroid of the sphere. */
//	virtual Vector Centroid();

};

#endif 