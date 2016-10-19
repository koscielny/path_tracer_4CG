/* Deals with the rendered scene */
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <string>
#include <vector>

#include "Scene.h"
#include "includeALL.h"
using namespace std;


/***  SCENE  ***/

// sets default values
Scene::Scene(char* file) {
	filename = "OUTPUT";
	diffuse = vec3(0,0,0);
	specular = vec3(0,0,0);
	shininess = 0;
	emission = vec3(0,0,0);
	indexofrefraction = 1;
	refractivity = 0;
	antialias = 2;
	shadowrays = 1;
	lightradius = 0;
	isLight = false;
	
	vec3 maxvec = vec3(DBL_MAX,DBL_MAX,DBL_MAX);
	vec3 minvec = vec3(DBL_MIN,DBL_MIN,DBL_MIN);
	sceneAABB = AABB2(maxvec, minvec);
	parse(file);
	clog << "Constructing KDTree... ";
//	KDTree = new TreeNode(objects,sceneAABB);
	bvh = new BVH(&objects);
	clog << "done"<<endl;
}

Scene::~Scene() {
	vector<Light*>::iterator lit;
	for (lit = lights.begin(); lit != lights.end(); lit++){
		delete *lit;
	}
	vector<Shape*>::iterator it;
	for(it=objects.begin(); it!=objects.end(); it++){
		delete *it;
	}
	vector<Material*>::iterator mit;
	for (mit = materials.begin(); mit != materials.end(); mit++){
		delete *mit;
	}
}
 
Ray Scene::castEyeRay(double i, double j){
	double alpha = (2.0*i-width)/width;

	alpha *=tan(fovx/2.0);
	double beta = (2.0*j-height)/height;
	beta *= tan(fovy/2.0);

	Ray ray(eye, alpha*u + beta*v - w);
	return ray;
}
vec3 Scene::castpoint(double i, double j){
	double d = glm::distance(glm::vec3(eye) , (camera_lookat));
	double alpha = (2.0*i-width)/width;
	alpha *=tan(fovx/2.0)*d;
	double beta = (2.0*j-height)/height;
	beta *= tan(fovy / 2.0)*d;
	
	vec3 planepoint(eye + alpha*u + beta*v - d*w);
	cout << "distance eye to plane" <<d<<endl ;
	return planepoint;
}

/***********************************************
** Parsing input and setting global variables **
***********************************************/

void Scene::setCoordinateFrame(vec3& lookat, vec3& up){
	w = glm::normalize(eye - lookat);
	u = glm::normalize(glm::cross(up, w));
	v = glm::cross(w,u);
}

void Scene::updateAABB(vec3& point){
	Vector p(point[0], point[1], point[2]);
	sceneAABB.min[0] = min(sceneAABB.min[0],p[0]);
	sceneAABB.max[0] = max(sceneAABB.max[0],p[0]);
	sceneAABB.min[1] = min(sceneAABB.min[1], p[1]);
	sceneAABB.max[1] = max(sceneAABB.max[1],p[1]);
	sceneAABB.min[2] = min(sceneAABB.min[2], p[2]);
	sceneAABB.max[2] = max(sceneAABB.max[2],p[2]);
}

void Scene::parseLine(string l, stack<mat4>& mv, vector<vec3>& verts, 
				vector<vec3>& normverts, vector<vec3>& norms){
					
	stringstream line(l);
	string cmd;
	line >> cmd;
	if (cmd[0]=='#' || cmd=="") { //comment or blank line
		return;
	} else if (cmd == "size") {
		line >> width >> height;
	} else if (cmd == "maxdepth") {
		line >> maxdepth;
	} else if (cmd == "output") {
		line >> filename;
		filename.erase(filename.end()-4,filename.end());
	} else if (cmd == "camera") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		eye = vec3(arg1,arg2,arg3);
		camera_eye = eye;
		

		line >> arg1 >> arg2 >> arg3;
		vec3 lookat = vec3(arg1,arg2,arg3);
		camera_lookat = lookat;

		line >> arg1 >> arg2 >> arg3;
		vec3 up = vec3(arg1,arg2,arg3);
		camera_up = up;

		setCoordinateFrame(lookat,up);
		line >> fovy;
		camera_fovy = fovy;
		fovy *= pi / 180.0;
		double d = height / (2.0 * tan( fovy * 0.5 ) );
		fovx = 2.0 * atan( width / (2.0 * d) );
	} else if (cmd == "SmoothGlassMaterial"){
		SmoothGlassMaterial *m = new SmoothGlassMaterial(line);
		materials.push_back(m);
	}
	else if (cmd == "diffuseMaterial"){
		cout << "yeahyeahyeah";
		
		diffuseMaterial *m = new diffuseMaterial(line);
		materials.push_back(m);
		cout << "m" << m->extinct_inside;
	}
	else if (cmd == "CookTorranceMaterial"){
		Material *m = new CookTorranceMaterial(line);
		materials.push_back(m);
	}
	else if (cmd == "FrostedGlassMaterial"){
		FrostedGlassMaterial *m = new FrostedGlassMaterial(line);
		materials.push_back(m);
	}
	else if (cmd == "specularMaterial"){
		specularMaterial *m = new specularMaterial(line);
		materials.push_back(m);
	}
	else if (cmd == "sphere") {
		double arg1, arg2, arg3, arg4;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		line >> arg4;
		mat4 trans = mv.top();
		trans *= Transform::translate(arg1,arg2,arg3);
		trans *= Transform::scale(arg4,arg4,arg4);
		Sphere* s = new Sphere(trans);
		sceneAABB.min[0] = min(s->boundingBox.min[0], sceneAABB.min[0]);
		sceneAABB.max[0] = max(s->boundingBox.max[0], sceneAABB.max[0]);
		sceneAABB.min[1] = min(s->boundingBox.min[1], sceneAABB.min[1]);
		sceneAABB.max[1] = max(s->boundingBox.max[1], sceneAABB.max[1]);
		sceneAABB.min[2] = min(s->boundingBox.min[2], sceneAABB.min[2]);
		sceneAABB.max[2] = max(s->boundingBox.max[2], sceneAABB.max[2]);
		//s->ambient = ambient;
		s->diffuse = diffuse;
		s->specular = specular;
		s->shininess = shininess;
		s->emission = emission;	
		s->indexofrefraction = indexofrefraction;
		s->refractivity = refractivity;
//		s->material = *diffmaterials.end();
		s->material = materials.back();
		objects.push_back(s);
		if (isLight){
			AreaLight* a = new AreaLight(s);
			lights.push_back(a);
		}
	} else if (cmd == "maxverts") {
		int maxverts;
		line >> maxverts;
		verts.reserve(maxverts);
	} else if (cmd == "maxvertnorms") {
		int maxvertnorms;
		line >> maxvertnorms;
		verts.reserve(maxvertnorms);
		verts.reserve(maxvertnorms);
	} else if (cmd == "vertex") {
		double arg1, arg2, arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		vec3 v(arg1,arg2,arg3);
		verts.push_back(v);
	} else if (cmd == "vertexnormal") {
		double arg1, arg2, arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		vec3 v(arg1,arg2,arg3);
		normverts.push_back(v);
		line >> arg1;
		line >> arg2;
		line >> arg3;
		vec3 n(arg1,arg2,arg3);
		norms.push_back(n);
	} else if (cmd == "tri") {
		int a1, a2, a3;
		line >> a1 >> a2 >> a3;
		mat4 top =  mv.top();
		vec3 v1 = vec3(top * vec4(verts[a1],1));
		vec3 v2 = vec3(top * vec4(verts[a2],1));
		vec3 v3 = vec3(top * vec4(verts[a3],1));
		Triangle* t = new Triangle(v1,v2,v3);
		updateAABB(v1);
		updateAABB(v2);
		updateAABB(v3);
		//t->ambient = ambient;
		t->diffuse = diffuse;
		t->specular = specular;
		t->shininess = shininess;
		t->emission = emission;
		t->indexofrefraction = indexofrefraction;
		t->refractivity = refractivity;
		t->material = materials.back();
//		cout << "t" << t->material->extinct_outside << t->material->extinct_inside << t->material->kd[0];

		objects.push_back(t);
		if (isLight){
			AreaLight* a = new AreaLight(t);
			lights.push_back(a);
		}
	} else if(cmd == "trinormal") {
		int a1,a2,a3;
		line >> a1 >> a2 >> a3;
		mat4 top =  mv.top();
		vec3 v1 = vec3(top * vec4(normverts[a1],1));
		vec3 v2 = vec3(top * vec4(normverts[a2],1));
		vec3 v3 = vec3(top * vec4(normverts[a3],1));
		top = glm::transpose(glm::inverse(top));
		vec3 n1 = vec3(top * vec4(norms[a1],0));
		vec3 n2 = vec3(top * vec4(norms[a2],0));
		vec3 n3 = vec3(top * vec4(norms[a3],0));
		NormTriangle* t = new NormTriangle(v1,v2,v3,n1,n2,n3);
		//t->ambient = ambient;
		t->diffuse = diffuse;
		t->specular = specular;
		t->shininess = shininess;
		t->emission = emission;
		t->indexofrefraction = indexofrefraction;
		t->refractivity = refractivity;
		t->material = materials.back();
		objects.push_back(t);
		if (isLight){
			AreaLight* a = new AreaLight(t);
			lights.push_back(a);
		}
	} else if(cmd == "translate") {
		double arg1,arg2,arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		mv.top() *= Transform::translate(arg1, arg2, arg3);
	} else if(cmd == "rotate") {
		double arg1,arg2,arg3,arg4;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		line >> arg4;
		mv.top() *= Transform::rotate(arg4,vec3(arg1,arg2,arg3));
	} else if (cmd=="scale") {
		double arg1,arg2,arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		mv.top() *= Transform::scale(arg1, arg2, arg3);
	} else if (cmd == "pushTransform") {
		mv.push(mv.top());
	} else if (cmd == "popTransform"){
		mv.pop();
	} else if (cmd == "diffuse") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		diffuse = vec3(arg1,arg2,arg3);
	} else if (cmd == "specular") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		specular = vec3(arg1,arg2,arg3);
	} else if (cmd == "shininess") {
		line >> shininess;
	} else if (cmd == "emission") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		emission = vec3(arg1,arg2,arg3);
		if (arg1 || arg2 || arg3){
			isLight = true;
		} else {
			isLight = false;
		}
	} else if (cmd == "indexofrefraction") {
		line >> indexofrefraction;
	} else if (cmd == "refractivity") {
		line >> refractivity;
	} else if (cmd == "antialias") {
		line >> antialias;
	}
	else if (cmd == "directionallight") {

		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		vec3 c(arg1, arg2, arg3);
		line >> arg1 >> arg2 >> arg3;
		vec3 dir(arg1, arg2, arg3);

		DirectionalLight *dlight = new DirectionalLight(c, dir);
		lights.push_back(dlight);
	}
	else if (cmd == "pointlight") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		vec3 colour(arg1, arg2, arg3);
		line >> arg1 >> arg2 >> arg3;
		vec3 poi(arg1, arg2, arg3);
		double con;
		double lin;
		double quad;
		line >> con >> lin >> quad;

		PointLight *plight = new PointLight(colour, poi, con, lin, quad);
		lights.push_back(plight);

	}
	//cout << cmd << endl;
}

void Scene::parse(char* filename) {
	ifstream file(filename, ifstream::in);
	string line;
	if (file.is_open()) {
		stack<mat4> modelview;
		modelview.push(mat4(1.0));
		vector<vec3> verts;
		vector<vec3> normverts;
		vector<vec3> norms;
		while (file.good()) {
			getline(file, line);
			parseLine(line, modelview, verts, normverts, norms);
		}

		//log
		FILE *fp;
		fp = fopen("logvec.txt", "w+");
		if (fp == NULL) { cout << "Ê²Ã´£¿"; return; }
		
		ofstream log;
		log.open("logvec.txt");
		if (log.is_open())
		{
			int i=0;
			for (vector<vec3>::iterator it = verts.begin(); it != verts.end(); it++)
			{
	//			log << it->x << " " << it->y << " " << it->z << " " << ++i << endl;
			}
			log.close();
		}
		fclose(fp);
		//
	} else {
		cerr << "Unable to open file " << filename << endl;
		exit(1);
	}
	file.close();
}
