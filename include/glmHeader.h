#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/epsilon.hpp>

// glm provides vector, matrix classes like glsl
// Typedefs to make code more readable 

typedef glm::dmat2 mat2;
typedef glm::dmat3 mat3;
typedef glm::dmat4 mat4;
typedef glm::dvec2 vec2;
typedef glm::dvec3 vec3;
typedef glm::dvec4 vec4;
