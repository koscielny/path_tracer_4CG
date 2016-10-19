#pragma once


//#include "Transform.h"
#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>
#include "glmHeader.h"
#define EPSILON1 0.0000001

using namespace std;
class Material{

public:
	//����Ҫ�ԣ��������䷽�� vec3 Exitance//hitpoint Ҫ�ܱ��޸ģ��ڷ�����΢Сλ�ƣ��Է�ֹ��һ���Խ�������artifact//indirection�����hitpointָ�����䷽��ķ���
	Material(){};
	Material(stringstream& line)
	{
		line >> extinct_inside;
		line >> extinct_outside;
	}//һ��Ҫ��ϸ��������Ĺ��캯���ĸ��ȱ����á�stringstream������ָ��


	//indirection��ʵ�����䷽�����Խ���Ϊ���
	virtual vec3 genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside)=0;
	virtual vec3 getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled) = 0;

	double extinct_inside;
	double extinct_outside;
};


class SmoothGlassMaterial : public Material
{
public:
	SmoothGlassMaterial(stringstream& line);
	//����Ҫ�ԣ��������䷽�� vec3 Exitance
	 vec3 genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside);

	 vec3 getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled);

	double refractiveIndex;
};

class diffuseMaterial : public Material
{
public:
	diffuseMaterial(stringstream& line);
		vec3 genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside);
	
		vec3 getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled);

public:
	vec3 kd;//diffuse Spectrum index������ϵ��

};

class specularMaterial : public Material
{
public:
	specularMaterial(stringstream& line);
		vec3 genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside);

		vec3 getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled);

	vec3 ks;
};


class CookTorranceMaterial : public Material
{
public:
	CookTorranceMaterial(stringstream& line);
	vec3 genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside);

	vec3 getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled);


	/* The index of the reflectance, and refractive index, distributions. */
	vec3 reflectance;
	double	refractiveIndex;
	/* The roughness of the surface. */
	double roughness;

};


class FrostedGlassMaterial : public Material
{
public:
	FrostedGlassMaterial(stringstream& line);
	virtual	vec3 genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside);

	virtual vec3 getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled);


	double	refractiveIndex;
	double roughness;

};

inline bool vec3delta(double x){ return (std::abs(x) <= 0.001); }
inline vec3 spherical(double phi, double theta)
{
	return vec3{ cosf(phi) * sinf(theta), cosf(theta), sinf(phi) * sinf(theta) };
}

/* Rotates a unit vector around a normal vector - more accurately, transforms the vector from world coordinates to
* normal coordinates defined by the orthonormal basis described by the normal vector as the upwards axis. */
vec3 rotate(vec3 &a, vec3 &n);
