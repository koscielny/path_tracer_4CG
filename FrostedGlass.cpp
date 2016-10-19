#include "Material.h"
#include "includeALL.h"
/* Creates the material from a scene file. */


FrostedGlassMaterial::FrostedGlassMaterial(stringstream& line) : Material(line)
{
	line >> refractiveIndex;
	line >> roughness;
}

vec3 FrostedGlassMaterial::genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside)
{
	/* Generate a random microfacet normal based on the Beckmann distribution with the given roughness. */
	double r1 = ((double)rand() / ((double)RAND_MAX + 1));
	double r2 = ((double)rand() / ((double)RAND_MAX + 1));
	double theta = atan(-pow(this->roughness, 2.0f) * log(1.0f - r1));
	double phi = 2.0f * PI * r2;
	vec3 m = spherical(phi, theta);

	/* Rotate the microfacet normal according to the actual surface normal. */
	m = rotate(m, normal);
	double n1, n2;
	/* Work out the correct n1 and n2 depending on the incident vector's direction relative to the normal. */
	if (isInside)
	{
		/* 光线在物体内部相交 */
		n1 = refractiveIndex;
		n2 = 1.0f;
		/* Flip the microfacet normal around. */
		//m = -m;
	}
	else
	{
		/* Incident and normal have opposite directions, so the ray is outside the material. */
		n2 = refractiveIndex;
		n1 = 1.0f;

	}
	double cosI = glm::dot(indirection, normal);



	/* Calculate the refracted angle's cosine. */
	float cosT = 1.0f - pow(n1 / n2, 2.0f) * (1.0f - pow(cosI, 2.0f));

	/* Check for total internal reflection. */
	if (cosT < 0.0f)
	{
		/* Total internal reflection occurred. */
		(*hitpoint) = (*hitpoint) + m * EPSILON1;
		return glm::reflect(-indirection, m);
	}

	/* Otherwise, finish computing the angle. */
	cosT = sqrt(cosT);

	/* Now, compute the Fresnel coefficients for reflection and refraction for a randomly polarized ray. */
	float R = (pow((n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT), 2.0f) + pow((n2 * cosI - n1 * cosT) / (n1 * cosT + n2 * cosI), 2.0f)) * 0.5f;

	/* Perform a random trial to decide whether to reflect or refract the ray. */
	if (((double)rand() / ((double)RAND_MAX + 1)) < R)
	{
		/* Reflection. */
		(*hitpoint) = (*hitpoint) + m * EPSILON1;
		return glm::reflect(-indirection, m);
	}
	else
	{
		/* Refraction. */
		(*hitpoint) = (*hitpoint) - m * EPSILON1;
		return -indirection * (n1 / n2) + m * ((n1 / n2) * cosI - cosT);
	}
}

/* This returns the reflectance for an incident and exitant vector. */
vec3 FrostedGlassMaterial::getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled)
{
	/* Work out the refractive indices. */
	double n1, n2;
	if (isInside)
	{
		/* Incident and normal have the same direction, ray is inside the material. */
		n1 = this->refractiveIndex;
		n2 = 1.0f;

	}
	else
	{
		/* Incident and normal have opposite directions, so the ray is outside the material. */
		n2 = this->refractiveIndex;
		n1 = 1.0f;
	}

	/* Check whether reflection or refraction occurred. */
	vec3 H;
	double D = 1.0f;
	if (glm::dot(-indirection , outdirection) < 0.0f)
	{
		/* Reflection occurred, find the half-angle vector. */
		H = glm::normalize(outdirection + indirection);

		/* If the ray was not importance-sampled, we need to take into account the distribution. */
		if (!sampled)
		{
			/* Get the half angle vector's angle with the normal. */
			double alpha = acos(glm::dot(H , normal));

			/* Compute the Beckmann distribution. */
			D = exp(-pow(tanf(alpha) / this->roughness, 2.0f));
		}
	}
	else
	{
		/* Refraction occurred, we have to find the microfacet normal. */
		double cI = std::abs(glm::dot(-indirection , normal));
		double cT = 1.0f - pow(n1 / n2, 2.0f) * (1.0f - pow(cI, 2.0f));
		H = (-indirection * (n1 / n2) - outdirection) / ((n1 / n2) * cI - cT);

		/* If the ray was not importance-sampled, we need to take into account the distribution. */
		if (!sampled)
		{
			/* Get the half angle vector's angle with the normal. */
			float alpha = acos(glm::dot(H , normal));

			/* Compute the Beckmann distribution. */
			D = exp(-pow(tanf(alpha) / this->roughness, 2.0f));
		}
	}

	/* Compute the geometric attenuation term. */
	double NdV = std::abs(glm::dot(-indirection , normal));
	double NdL = std::abs(glm::dot(normal , outdirection));
	double VdH = std::abs(glm::dot(-indirection , H));
	double NdH = std::abs(glm::dot(normal , H));
	double G = std::min(1.0, std::min(2.0 * NdH * NdV / VdH, 2.0 * NdH * NdL / VdH));

	/* Compute the microfacet normalization term. */
	double norm = 1.0f / (PI * pow(this->roughness, 2.0f) * pow(NdH, 4.0f));

	/* Compute the reflectance (note the lambertian term cancels a dot product out).
	* Also note we do NOT use the fresnel term if the ray was importance-sampled,
	* since we were already weighting the probability of reflection and refraction
	* with it when sampling the BTDF. */
	return norm * (D * G) / (NdV)*vec3(1.0,1.0,1.0);
}


