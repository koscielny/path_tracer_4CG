
#include "Material.h"
#include "includeALL.h"

vec3 rotate(vec3 &a, vec3 &n)
{
	Vector av(a[0], a[1], a[2]);
	Vector nv(n[0], n[1], n[2]);

	Vector rv = rotate(av, nv);

	return vec3(rv.x, rv.y, rv.z);
}

