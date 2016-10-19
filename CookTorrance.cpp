
#include "Material.h"
#include "includeALL.h"





CookTorranceMaterial::CookTorranceMaterial(stringstream &line) :Material(line)
{
	line >> reflectance[0];
	line >> reflectance[1];
	line >> reflectance[2];

	line >> refractiveIndex;
	line >> roughness;
}

//indirection其实是入射方向反向，以交点为起点
vec3 CookTorranceMaterial::genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside)
{
	/* Move the origin outside the surface slightly. */
	(*hitpoint) = (*hitpoint) + normal * EPSILON1;

	/* Generate a random microfacet normal based on the Beckmann distribution with the given roughness. */
	double r1 = ((double)rand() / ((double)RAND_MAX + 1));
	double r2 = ((double)rand() / ((double)RAND_MAX + 1));
	double theta = atan(-pow(this->roughness, 2.0f) * log(1.0f - r1));
	double phi = 2.0f * PI * r2;
	vec3 m = spherical(phi, theta);
	m = rotate(m, normal);

	/* Reflect the incident vector accordingly. */
	return glm::reflect(-indirection, m);
}

vec3 CookTorranceMaterial::getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled)
{
	/* Compute the half-angle vector. */
	vec3 H = glm::normalize(outdirection + indirection);

	/* If the ray was not importance-sampled, we need to take into account the distribution. */
	double D = 1.0;
	if (!sampled)
 	{
		/* Get the half angle vector's angle with the normal. */
		double alpha = acos(glm::dot(H , normal));

		/* Compute the Beckmann distribution. */
		D = exp(-pow(tanf(alpha) / this->roughness, 2.0f));
	}

	/* Compute the refractive indices. */
	double n2 = this->refractiveIndex;//	double n2 = this->refractiveIndex->Lookup(wavelength);
	double n1 = 1.0;

	/* Compute the theoretical reflected and refracted angles. */
	double cosI = std::abs(glm::dot(-indirection , normal));
	double cosT = sqrtf(1.0 - pow(n1 / n2, 2.0) * (1.0 - pow(cosI, 2.0)));

	/* Compute the Fresnel term for the incident vector. */
	double F = (pow((n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT), 2.0f) + pow((n2 * cosI - n1 * cosT) / (n1 * cosT + n2 * cosI), 2.0f)) * 0.5f;

	/* Compute the geometric attenuation term. */
	double NdL = std::abs(glm::dot(normal , outdirection));
	double VdH = std::abs(glm::dot(-indirection , H));
	double NdH = std::abs(glm::dot(normal , H));
	double NdV = cosI;
	double G = std::min(1.0, std::min(2.0 * NdH * NdV / VdH, 2.0 * NdH * NdL / VdH));

	/* Compute the microfacet normalization term. */
	double norm = 1.0f / (PI * pow(this->roughness, 2.0) * pow(NdH, 4.0));

	/* Compute the reflectance (note the lambertian term cancels a dot product out). */
	return norm * this->reflectance * (F * D * G) / (NdV);

}



