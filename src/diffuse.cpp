#include "Material.h"
#include "includeALL.h"

diffuseMaterial::diffuseMaterial(stringstream& line) :Material(line)
{
	line >> kd[0];
	line >> kd[1];
	line >> kd[2];
}

vec3 diffuseMaterial::genExitSample(vec3 *hitpoint, vec3 indirection, vec3 normal, bool isInside)
{//Œ¥—È÷§
	*hitpoint += normal*EPSILON1;

	double u1 = ((double)rand() / (double)RAND_MAX);
	double u2 = ((double)rand() / (double)RAND_MAX);

	vec3 y = vec3(normal);
	vec3 h = vec3(normal);
	double theta = acos(sqrt(1.0 - u1));
	double phi = 2.0 * M_PI * u2;
	double xs = sin(theta) * cos(phi);
	double ys = cos(theta);
	double zs = sin(theta) * sin(phi);
	if ((abs(h[0]) <= abs(h[1])) && (abs(h[0]) <= abs(h[2])))
		h[0] = 1.0;
	else if ((abs(h[1]) <= abs(h[0])) && (abs(h[1]) <= abs(h[2])))
		h[1] = 1.0;
	else
		h[2] = 1.0;
	vec3 x = glm::cross(h, y);
	vec3 z = glm::cross(x, y);

	vec3 direction = xs * x + ys * y + zs * z;
	return direction;
}

vec3 diffuseMaterial::getSpectrumIndex(vec3 hitpoint, vec3 indirection, vec3 outdirection, vec3 normal, bool isInside, bool sampled)
{
	if (!sampled)
	{
		/* Otherwise, use the uniform sampling formulation. */
		return 2.0 * this->kd * std::abs(glm::dot(outdirection, normal));
	}
	else
	{
		return kd;
	}
}