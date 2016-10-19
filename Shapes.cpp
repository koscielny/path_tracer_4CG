/* Class definitions for Triangles and Spheres */
#include <iostream>
#include <vector>

#include "Shapes.h"
//#include "KDTree.h"
#include "includeALL.h"

//浮点数误差1e-6f影响面片求交、和面片大小也有关
#define EPSILON3 0.0000001

using namespace std;

/***  RAY  ***/
Ray::Ray(const vec3& o,const vec3& d){
	origin = o;
	direction = glm::normalize(d);
}

vec3 Ray::getPoint(double t) {
	return origin + t*direction;
}

/***  INTERSECTION  ***/
Intersection::Intersection(vector<Shape*>& objects, Ray& ray) {
	double min_t = DBL_MAX;
	primative = NULL;
	for(vector<Shape*>::iterator prim=objects.begin(); prim!=objects.end(); prim++) {
		double t = (*prim)->intersect(ray);
		if (t>0 && t<min_t){
			primative = *prim;
			min_t = t;
		}
	}
	point = ray.getPoint(min_t);
	sourceDirection = -ray.direction;

	travel = min_t;
}


/***  TRIANGLE  ***/
Triangle::Triangle(vec3 point0, vec3 point1, vec3 point2) {
	p0 = point0;
	p1 = point1;
	p2 = point2;
	n0 = glm::normalize(glm::cross(p1-p0,p2-p0));
	double minx = min(min(p0[0],p1[0]),p2[0]);
	double maxx = max(max(p0[0],p1[0]),p2[0]);
	double miny = min(min(p0[1],p1[1]),p2[1]);
	double maxy = max(max(p0[1],p1[1]),p2[1]);
	double minz = min(min(p0[2],p1[2]),p2[2]);
	double maxz = max(max(p0[2],p1[2]),p2[2]);
	vec3 minvec = vec3(minx,miny,minz);
	vec3 maxvec = vec3(maxx,maxy,maxz);
	//aabb = AABB(minvec, maxvec);
	this->boundingBox = AABB2(minvec, maxvec);

	this->centroid = Vector((this->p0 + this->p1 + this->p2)[0] / 3.0f, (this->p0 + this->p1 + this->p2)[1] / 3.0f, (this->p0 + this->p1 + this->p2)[2] / 3.0f);

	assertnorm = false;
}


double Triangle::intersect(Ray& ray){
	vec3 col1 = p1-p0;
	vec3 col2 = p2-p0;
	mat3 M = mat3(col1, col2, -ray.direction);
	double det = glm::determinant(M);
	if (det<EPSILON3 && det>-EPSILON3) return -1.0;
	
	M[0] = ray.origin-p0;
	double alpha = glm::determinant(M)/det;
	if (alpha<0.0 || alpha>1.0) return -1.0;
	M[0] = col1;
	M[1] = ray.origin-p0;

	double beta = glm::determinant(M)/det;
	if (beta<0.0 || beta+alpha>1.0) return -1.0;
	M[1] = col2;
	M[2] = ray.origin-p0;
	return glm::determinant(M)/det;
}

vec3 Triangle::getNormal(const vec3& hitpoint){
	return n0;//不安全 ，可能反向
}

vec3 Triangle::getNormal(const Intersection& hit){
	if (glm::dot(n0, hit.sourceDirection) < 0)
		return -n0;
	else return n0;
}

double Triangle::getSubtendedAngle(const vec3& origin) {
	vec3 a = p0 - origin;
	vec3 b = p1 - origin;
	vec3 c = p2 - origin;
	
	double det = glm::determinant( mat3(a,b,c) );
	
	double a1 = glm::length(a);
	double b1 = glm::length(b);
	double c1 = glm::length(c);
	
	double result = a1*b1*c1 + glm::dot(a,b)*c1 + 
				glm::dot(a,c)*b1 + glm::dot(b,c)*a1;
					
	result = atan2(det, result);
	//if (result < 0.0) result += M_PI;
	result *= 2.0;
	return abs(result);
}

inline vec3 genSample(const vec3& p0, const vec3& p1, const vec3& p2) {
	double u1 = double(rand()) / double(RAND_MAX);
	double u2 = double(rand()) / double(RAND_MAX);
	u1 = sqrt(u1);
	return (1.0-u1)*p0 + u1*(1.0-u2)*p1 + u1*u2*p2;
}

inline vec3 genUniSphereSample(const vec3& p0, const vec3& p1, const vec3& p2) {
	double u1 = double(rand()) / double(RAND_MAX);
	double u2 = double(rand()) / double(RAND_MAX);
//	u1 = sqrt(u1);
//	return (1.0 - u1)*p0 + u1*(1.0 - u2)*p1 + u1*u2*p2;
}



vec3 getShadeIndex(Intersection it, vec3 exitent, vec3 hitnormal, bool isInside)
{

	return  it.primative->material->getSpectrumIndex(it.point, it.sourceDirection, exitent, hitnormal, isInside, false);
}


/*
returns shade
*/
vec3 Triangle::shade(const Intersection& hit, BVH* tree, vec3 s_norm, bool single_ray) {//single_ray:对光源单采样，计算间接光照用，计算直接光照不用
	//vec3 s_norm = hit.primative->getNormal(hit);
//	if (!hit.primative->assertnorm)
//	{
//		if (glm::dot(s_norm, hit.sourceDirection) < 0){ s_norm = -s_norm; }
//	}
	//vec3 s_diffuse = hit.primative->diffuse;
	//vec3 s_specular = hit.primative->specular;
	//double s_shininess = hit.primative->shininess;



	if (single_ray){
		vec3 light_samp = genSample(p0, p1, p2);
		
		vec3 dir = glm::normalize(light_samp - hit.point);
		Ray ray = Ray(hit.point + EPSILON3*s_norm, dir);//+ EPSILON * s_norm

		Intersection light_hit;
		tree->getIntersection(ray, &light_hit, false);//shadow_ray
		vec3 l_norm = this->getNormal(light_hit);

		if (light_hit.primative != this) {
			return vec3(0,0,0); // is a shadow
		}
		else return getShadeIndex(hit, dir,s_norm,false)*(this->emission);
		/* Calculate shading */
//		vec3 shade = max(0.0,glm::dot(s_norm, dir)) * s_diffuse;
//		vec3 half = glm::normalize(hit.sourceDirection+dir);
//		double phong = pow( max(0.0,glm::dot(half,s_norm)) , s_shininess);
//		shade += phong * s_specular;
//		shade *= this->emission;

		/* weigh by dA */
		double dist = glm::distance(hit.point, light_hit.point);
		dist *= dist; // square distance
		double cos_weight = glm::dot(s_norm, dir);
		cos_weight *= glm::dot(l_norm, -dir);
		/* Multiply by dA */
		cos_weight *= 0.5 * glm::dot(p1-p0, p2-p0);
		cos_weight /= dist;
		return getShadeIndex(hit, dir, s_norm, false)*(this->emission)*cos_weight; //cos_weight * shade;
	}	
	
	
	
	/* Setup array to rotate through parts of triangle */
	/* Slightly hardcoded for improved performace */
	vec3 centroid = p0 + p1 + p2;
	centroid /= 3.0;
	vec3 corners[] = {p0, (p0+p1)/2.0, p1, (p1+p2)/2.0, p2, (p2+p0)/2.0, p0};
	
	const int num_samples = 6;
	vec3 color = vec3(0,0,0);
	for (int i=0; i<num_samples; i+=1){
		
		vec3 light_samp = genSample(corners[i], corners[i+1], centroid);
		
		vec3 dir = glm::normalize(light_samp - hit.point);
		Ray ray = Ray(hit.point + EPSILON3*s_norm, dir);//+ EPSILON * s_norm//s_norm可能朝背面
		Intersection light_hit;
		tree->getIntersection(ray, &light_hit, false);//shadow_ray
		vec3 l_norm = this->getNormal(light_hit);

		if (light_hit.primative != this) {
			continue; // is a shadow
		}
		else return getShadeIndex(hit, dir, s_norm, false)*(this->emission);


		/* Calculate shading */
//		vec3 shade = max(0.0,glm::dot(s_norm, dir)) * s_diffuse;
//		vec3 half = glm::normalize(hit.sourceDirection+dir);
//		double phong = pow( max(0.0,glm::dot(half,s_norm)) , s_shininess);
//		shade += phong * s_specular;
//		shade *= this->emission;

		/* weigh by dA */
		double dist = glm::distance(hit.point, light_hit.point);
		dist *= dist; // square distance
		double cos_weight = glm::dot(s_norm, dir);
		cos_weight *= glm::dot(l_norm, -dir);
		cos_weight /= dist;
		
		color += getShadeIndex(hit, dir, s_norm, false)*(this->emission)*cos_weight;
	}
	color /= double(num_samples);
	return color * 0.5 * glm::dot(p1-p0, p2-p0);
}


/*
vec3 Triangle::getTexture(vec3& hit){
	mat2 M = mat2(p1[0]-p0[0], p1[1]-p0[1], p2[0]-p0[0], p2[1]-p0[1]);
	double det = glm::determinant(M);
	M[0][0] = hit[0]-p0[0];
	M[0][1] = hit[1]-p0[1];
	double alpha = glm::determinant(M)/det;//p0 to p1
	M[0][0] = p1[0]-p0[0];
	M[0][1] = p1[1]-p0[1];
	M[1][0] = hit[0]-p0[0];
	M[1][1] = hit[1]-p0[1];
	double beta = glm::determinant(M)/det;//p0 to p2
}

*/
/***  NORMTRIANGLE  ***/
NormTriangle::NormTriangle(vec3 point0, vec3 point1, vec3 point2,
							vec3 norm0, vec3 norm1, vec3 norm2) {
	p0 = point0;
	p1 = point1;
	p2 = point2;
	n0 = glm::normalize(norm0);
	n1 = glm::normalize(norm1);
	n2 = glm::normalize(norm2);

	assertnorm = true;
}

vec3 NormTriangle::getNormal(const vec3& hitpoint){
	mat2 M = mat2(p1[0]-p0[0], p1[1]-p0[1], p2[0]-p0[0], p2[1]-p0[1]);
	double det = glm::determinant(M);
	M[0][0] = hitpoint[0] - p0[0];
	M[0][1] = hitpoint[1] - p0[1];
	double beta = glm::determinant(M)/det;
	M[0][0] = p1[0]-p0[0];
	M[0][1] = p1[1]-p0[1];
	M[1][0] = hitpoint[0] - p0[0];
	M[1][1] = hitpoint[1] - p0[1];
	double gamma = glm::determinant(M)/det;
	return glm::normalize((1-beta-gamma)*n0 + beta*n1 + gamma*n2);
}
vec3 NormTriangle::getNormal(const Intersection& hit){
	mat2 M = mat2(p1[0] - p0[0], p1[1] - p0[1], p2[0] - p0[0], p2[1] - p0[1]);
	double det = glm::determinant(M);
	M[0][0] = hit.point[0] - p0[0];
	M[0][1] = hit.point[1] - p0[1];
	double beta = glm::determinant(M) / det;
	M[0][0] = p1[0] - p0[0];
	M[0][1] = p1[1] - p0[1];
	M[1][0] = hit.point[0] - p0[0];
	M[1][1] = hit.point[1] - p0[1];
	double gamma = glm::determinant(M) / det;
	return glm::normalize((1 - beta - gamma)*n0 + beta*n1 + gamma*n2);	
	//	if (glm::dot(n0, hit.sourceDirection) < 0)
	//		return -n0;
	//	else return n0;
}

/***  SPHERE  ***/
Sphere::Sphere(mat4 trans){
	mv = trans;
	inv = glm::inverse(mv);
	
	/*Bounding Box*/
	mat4 S = mat4(1.0);
	S[3][3] = -1.0;
	mat4 R = mv*S*glm::transpose(mv);
	double minx = (R[0][3]+sqrt(R[0][3]*R[0][3] -R[3][3]*R[0][0]))/R[3][3];
	double maxx = (R[0][3]-sqrt(R[0][3]*R[0][3] -R[3][3]*R[0][0]))/R[3][3];
	double miny = (R[1][3]+sqrt(R[1][3]*R[1][3] -R[3][3]*R[1][1]))/R[3][3];
	double maxy = (R[1][3]-sqrt(R[1][3]*R[1][3] -R[3][3]*R[1][1]))/R[3][3];
	double minz = (R[2][3]+sqrt(R[2][3]*R[2][3] -R[3][3]*R[2][2]))/R[3][3];
	double maxz = (R[2][3]-sqrt(R[2][3]*R[2][3] -R[3][3]*R[2][2]))/R[3][3];
	
	vec3 minvec = vec3(minx,miny,minz);
	vec3 maxvec = vec3(maxx,maxy,maxz);
	//aabb = AABB(minvec, maxvec);
	this->boundingBox = AABB2(minvec, maxvec);
	this->center = Vector(R[0][3], R[1][3], R[2][3]);
	assertnorm = true;
}


/*
double Sphere::intersect(Ray& ray) {//妈呀，搞个坐标的debug吧@！！！
	vec3 direction = glm::normalize(vec3(inv * vec4(ray.direction, 0.0)));
	vec3 origin =vec3(inv * vec4(ray.origin,1.0));
#ifdef DEBUG
	cout << "原方向" << " " << ray.direction.x << " " << ray.direction.y << " " << ray.direction.z << endl;
	cout << "射线起点" << " " << origin.x << " " << origin.y << " " << origin.z << endl;
	cout << "方向" << " " << direction.x << " " << direction.y << " " << direction.z << endl;
	if (glm::distance(glm::vec3(0, 0, 0), glm::vec3(origin)) <= 1.0)
		cout << "球内" << glm::distance(glm::vec3(0, 0, 0), glm::vec3(origin)) << endl;

#endif // DEBUG

	double b = 2.0f * glm::dot(direction, origin);
	double c = glm::dot(origin,origin) - 1.0;
	double det = b*b - 4.0f*c;
	if (det < 0.0) {
//		cout << "举手";
		return -1.0;
	}
	det = sqrt(det);
#ifdef DEBUG
	cout << "det" << det << endl;
	cout << "b" << b << endl;
	cout << "c" << c << endl;
#endif // DEBUG


	double t1 = (-b+det)*0.5;
	double t2 = (-b-det)*0.5;
	
	if (t1<=0.0 && t2<=0.0) return -1.0;
	
	if (t2>0.0f) {
		vec4 hit = mv * vec4(origin + t2*direction, 1.0);
		cout << "举手t2";
		return glm::distance(ray.origin,vec3(hit));

	} else { //t1 is closer
		vec4 hit = mv * vec4(origin+t1*direction,1.0);
		cout << "举手t1" << " " << t1 <<"!!!"<< (origin + t1*direction).x << " " << (origin + t1*direction).y << " " << (origin + t1*direction).z<<endl;
		return glm::distance(ray.origin, vec3(hit));
	}
}
*/

double Sphere::intersect(Ray& ray) {//妈呀，搞个坐标的debug吧@！！！
	vec3 direction = glm::normalize(vec3(inv * vec4(ray.direction, 0.0)));
	vec3 origin = vec3(inv * vec4(ray.origin, 1.0));
#ifdef DEBUG
	cout << "原方向" << " " << ray.direction.x << " " << ray.direction.y << " " << ray.direction.z << endl;
	cout << "射线起点" << " " << origin.x << " " << origin.y << " " << origin.z << endl;
	cout << "方向" << " " << direction.x << " " << direction.y << " " << direction.z << endl;
	if (glm::distance(glm::vec3(0, 0, 0), glm::vec3(origin)) <= 1.0)
		cout << "球内" << glm::distance(glm::vec3(0, 0, 0), glm::vec3(origin)) << endl;

#endif // DEBUG
	vec3 vpc = -origin;
	if (glm::dot(direction, vpc) < 0)
	{//圆心在射线后
		if (glm::length(vpc) > 1)
		{
			return -1.0;
			//无交点
		}
		else if (glm::length(vpc) == 1){
			double t = 0;
			vec4 hit = mv * vec4(origin + t*direction, 1.0);
			return glm::distance(ray.origin, vec3(hit));
		}
			
		else
		{//球心在射线后，射线起点在球内
			//求球心向射线投影点 pc
			double t_pc = glm::dot(direction, vpc) / glm::length(direction);
			vec3 pc = origin + t_pc* direction;//t<0//direction是单位向量

			double dist = sqrt(1 - glm::dot(pc, pc));
			double t = (dist + t_pc);//减t的绝对值即加t 
			vec4 hit = mv * vec4(origin + t*direction, 1.0);
			return glm::distance(ray.origin, vec3(hit));
		}
	}
	else
	{//球心能投影在射线正方向			
		//求球心向射线投影点 pc
		double t_pc = glm::dot(direction, vpc) / glm::length(direction);
		vec3 pc = origin + t_pc* direction;//t<0//direction是单位向量

		if (glm::length(pc) > 1)
		{
			//无交点
			return -1.0;
		}
		else 
		{
			double dist = sqrt(1 - glm::dot(pc,pc));
			if (glm::length(vpc) > (1))
			{
//				cout << "举手t1";
				double t = (t_pc - dist);
				vec4 hit = mv * vec4(origin + t*direction, 1.0);
				return glm::distance(ray.origin, vec3(hit));
			}
			else //if (glm::length(vpc) < (1-EPSILON))
			{
//				cout << "举手t2";
				double t = (t_pc + dist);
				vec4 hit = mv * vec4(origin + t*direction, 1.0);
				return glm::distance(ray.origin, vec3(hit));
			}
		}


	}

#ifdef DEBUG
	cout << "det" << det << endl;
	cout << "b" << b << endl;
	cout << "c" << c << endl;
#endif // DEBUG

}
vec3 Sphere::getNormal(const vec3& hitpoint){
	return glm::normalize(vec3(glm::transpose(inv)*inv*vec4(hitpoint, 1.0)));
}
vec3 Sphere::getNormal(const Intersection& hit){
	return glm::normalize(vec3(glm::transpose(inv)*inv*vec4(hit.point, 1.0)));
	//	if (glm::dot(n0, hit.sourceDirection) < 0)
//		return -n0;
//	else return n0;
}


double Sphere::getSubtendedAngle(const vec3& origin) {
	cout << "THIS IS NOT IMPLEMENTED" << endl;
	exit(1);
}

vec3 Sphere::shade(const Intersection& hit, BVH* tree,vec3 s_norm,bool derp) {
	cout << "THIS IS NOT IMPLEMENTED" << endl;
	exit(1);
}


/***  AABB  ***/
AABB::AABB(vec3& minarg, vec3& maxarg) {
	this->aabbmin = minarg;
	this->aabbmax = maxarg;
	this->center = (minarg+maxarg)/2.0;
	
}

inline
bool intersect1D(double start, double dir, double axisMin, double axisMax, double& near, double& far){
	// Parallel
	if (dir<EPSILON3 && dir>-EPSILON3){
		return (start>axisMin) && (start<axisMax);
	}
	
	//intersection parameters
	double t0 = (axisMin-start)/dir;
	double t1 = (axisMax-start)/dir;
	
	if(t0>t1){
		double temp = t1;
		t1 = t0;
		t0 = temp;
	}
	
	near = max(t0,near);
	far = min(t1,far);
	
	if(near>far || far<0.0) return false;
	
	return true;
}

bool
operator<(const vec3 &vecA, const vec3 &vecB){ 
	if (vecA[0] > vecB[0]) return false;
	if (vecA[1] > vecB[1]) return false;
	if (vecA[2] > vecB[2]) return false;
	return true;
}

bool
operator>(const vec3 &vecA, const vec3 &vecB){ 
	if (vecA[0] < vecB[0]) return false;
	if (vecA[1] < vecB[1]) return false;
	if (vecA[2] < vecB[2]) return false;
	return true;
}

double AABB::intersect(Ray& ray){
	
//	if (ray.origin < aabbmax && ray.origin > aabbmin){
//		return EPSILON; // always first if inside	
//	}
	
	
	double far = DBL_MAX;
	double near = DBL_MIN;
	
	if (!intersect1D(ray.origin[0],ray.direction[0],aabbmin[0],aabbmax[0],near,far)) return false;
	if (!intersect1D(ray.origin[1],ray.direction[1],aabbmin[1],aabbmax[1],near,far)) return false;
	if (!intersect1D(ray.origin[2],ray.direction[2],aabbmin[2],aabbmax[2],near,far)) return false;
	
	return near;
}






