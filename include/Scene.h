#pragma once


#ifndef SCENE_H  
#define SCENE_H
#include "glmHeader.h"
//#include "KDTree.h"
class Ray;
class BVH;
class AABB2;
class Material;
class Vector;
class Shape;
class Light;
#include "util\aabb.hpp"
class Scene {
public:
	vec3 w;
	vec3 u;
	vec3 v;
	int width;
	int height;
	int maxdepth;
	std::string filename;
	double fovy;
	double fovx;
	vec3 eye;
	glm::vec3 camera_eye;
	glm::vec3 camera_lookat;
	glm::vec3 camera_up;
	double camera_fovy;

	std::vector<Shape*> objects;
	std::vector<Material*> materials;
	//		std::vector<diffuseMaterial*> diffmaterials;
	//		AABB sceneAABB;
	AABB2 sceneAABB;
	//TreeNode* KDTree;
	BVH* bvh;

	vec3 diffuse;
	vec3 specular;
	double shininess;
	vec3 emission;
	double indexofrefraction;
	double refractivity;

	double antialias;
	int shadowrays;
	double lightradius;

	std::vector<Light*> lights;
	//std::vector<Shape*> lights;
	bool isLight;

	Scene(char*);
	~Scene();
	Ray castEyeRay(double, double);
	vec3 castpoint(double i, double j);
	void setCoordinateFrame(vec3&, vec3&);
	void parseLine(std::string, std::stack<mat4>&,
		std::vector<vec3>&, std::vector<vec3>&, std::vector<vec3>&);
	void parse(char*);
	void updateAABB(vec3&);
};


#endif