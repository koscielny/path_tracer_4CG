#pragma once

/* Class definition for lights */
#ifndef LIGHT_H
#define LIGHT_H
#include "glmHeader.h"
//#include "KDTree.h"
//#include <sceneaccel/bvh.hpp>
class Intersection;
class BVH;
class Shape;
//class Material;
class Light {
public:
	//Light();
	virtual ~Light(){};
	virtual vec3 shade(const Intersection&, BVH* tree, vec3 s_norm, bool single_ray) = 0;
	vec3 color;
};

class AreaLight : public Light 
{
public:
	AreaLight(Shape* s) : shapelight(s){};//干脆不用color属性 用Shape里的emission代替以实现shade()
	~AreaLight(){};//shapelight指针交给Vector<Shapes*> shapes 释放

	vec3 shade(const Intersection&, BVH* tree, vec3 s_norm, bool single_ray);
	Shape* shapelight;
};

class DirectionalLight : public Light {
public:
	~DirectionalLight(){};
	DirectionalLight(const vec3& color,const vec3& dir);
	vec3 shade(const Intersection& hit, BVH* tree, vec3 s_norm, bool single_ray);
	bool isVisible(const vec3& point, BVH* tree);
	vec3 direction;
};

class PointLight : public Light {
public:
	~PointLight(){};
	PointLight(const vec3& color,const vec3& p, double,double,double);
	vec3 shade(const Intersection& hit, BVH* tree, vec3 s_norm, bool single_ray);
	bool isVisible(const vec3& point, BVH* tree);
	vec3 point;
	double constant;
	double linear;
	double quadratic;
	double lightradius=0;
	int shadowrays=1;//
};

#endif