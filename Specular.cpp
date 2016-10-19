#include "Material.h"
#include "includeALL.h"

specularMaterial::specularMaterial(stringstream& line) :Material(line)
{
	line >> ks[0];
	line >> ks[1];
	line >> ks[2];
}
vec3 specularMaterial::genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside)
{


	/* Move the origin outside the surface slightly. */
	(*hitpoint) = (*hitpoint) + normal * EPSILON1;

	/* Just return the reflected angle. */
	return glm::reflect(-indirection, normal);
}

vec3 specularMaterial::getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled)
{//if importance sampled
	if (!sampled)    
	{
		/* Otherwise, find the expected reflected vector. */
		vec3 expected = glm::reflect(-indirection, normal);

		/* Return full reflectance if the exitant vector is the expected vector, otherwise zero. */
		if (!delta(1.0f - glm::dot(expected , outdirection))) return vec3(0,0,0);
		else return this->ks;
	}
	else return ks;
}

//#define delta(x) (float)(std::abs(x) <= 1e-3f)