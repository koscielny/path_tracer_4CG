#include "Material.h"
#include "includeALL.h"

SmoothGlassMaterial::SmoothGlassMaterial(stringstream& line) :Material(line)
{
	line >> refractiveIndex;
}
vec3 SmoothGlassMaterial::genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside)
{
	/* Work out the correct n1 and n2 depending on the incident vector's direction relative to the normal. */
	double n1, n2;
	if (isInside)
	{
		/* 光线在物体内部相交 */
		n1 = refractiveIndex;
		n2 = 1.0f;

	}
	else
	{
		/* Incident and normal have opposite directions, so the ray is outside the material. */
		n2 = refractiveIndex;
		n1 = 1.0f;

	}
	double cosI = glm::dot(indirection, normal);

	/* Calculate the refracted angle's cosine. */
	double cosT = 1.0f - pow(n1 / n2, 2.0f) * (1.0f - pow(cosI, 2.0f));

	/* Check for total internal reflection. */
	if (cosT < 0.0f)
	{
		/* 只发生镜面反射Total internal reflection occurred. */
		(*hitpoint) += normal * EPSILON1;
		return glm::reflect(-indirection, normal);
	}

	/* Otherwise, finish computing the angle. */
	cosT = sqrt(cosT);

	/* Now, compute the Fresnel coefficients for reflection and refraction for a randomly polarized ray. */
	float R = (pow((n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT), 2.0f) + pow((n2 * cosI - n1 * cosT) / (n1 * cosT + n2 * cosI), 2.0f)) * 0.5f;

	/* Perform a random trial to decide whether to reflect or refract the ray. */
	if (((double)rand() / ((double)RAND_MAX + 1)) < R)
	{
		/* Reflection. */
		(*hitpoint) += normal * EPSILON1;
		return glm::reflect(-indirection, normal);
	}
	else
	{
		/* Refraction. */
		(*hitpoint) -= normal * EPSILON1;
		return -indirection * (n1 / n2) + normal * ((n1 / n2) * cosI - cosT);
	}
}

vec3 SmoothGlassMaterial::getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled)
{
	if (!sampled)
	{
		/* Otherwise, find the expected reflected vector. */
		vec3 expected = glm::reflect(-indirection, normal);

		/* Return full reflectance if the exitant vector is the expected vector, otherwise zero. */
		if (!delta(1.0f - glm::dot(expected, outdirection))) return vec3(0, 0, 0);
		else return vec3(1.0, 1.0, 1.0);
	}
	else {
		/* The reflectance will here be constant, because we've already weighted the
		* reflection/refraction probability according to the Fresnel equations. */

		return vec3(1.0, 1.0, 1.0);
	}
}