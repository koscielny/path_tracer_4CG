/* This is the main file for the raytracer */

#include "includeALL.h"
#define draw_debug
#define BPP 24
#define EPSILON4 0.00001
#define BUFFER_OFFSET(i) (reinterpret_cast<void*>(i))
#define air_RIndex 1.0
using namespace std;

/* Paramaters */
double rays_cast;
Scene* scene;
int height;
int width;
int numFrames = 0;

/* Shaders */
GLuint vertexshader;
GLuint fragmentshader;
GLuint shaderprogram;
GLuint texture[2];
FIBITMAP* bitmap;
FIBITMAP* colortexel;
vec3 * pixels;
vec3 * direct_pixels;
int update;

const vec3 X = vec3(1, 0, 0);
const vec3 Y = vec3(0, 1, 0);
const vec3 Z = vec3(0, 0, 1);

void drawray();
void draw_ray_debug(int, int);
//？？？
vec3 genpositiveNormal(Intersection it, bool *insideflag);

static int mousex, mousey;
static int mcx, mcy;
/*GLUT提供鼠标动作检测能力。有两种GLUT处理的motion：active motion和passive motion。Active motion是指鼠标移动并且有一个鼠标键被按下。Passive motion是指当鼠标移动时，并有没鼠标键按下。如果一个程序正在追踪鼠标，那么鼠标移动期间，每一帧将产生一个结果。*/
void mouseclick(int x, int y)
{
	mcx = x;
	mcy = glutGet(GLUT_WINDOW_HEIGHT) - y;
	glutPostRedisplay();

}
//鼠标移动回调函数
void mousemove(int x, int y)
{
	mousex = x;
	mousey = glutGet(GLUT_WINDOW_HEIGHT) - y;
	glutPostRedisplay();
}

/*-------------------------------------------------------------------*//*!
																	   * \摘要 	在窗口中指定位置指定大小的矩形内显示鼠标的当前位置
																	   * \参数  	x:           指定矩形的左小角的x坐标
																	   y:           指定矩形的左小角的y坐标
																	   width:       指定矩形的宽度
																	   height:      指定矩形的高度
																	   win_width:   窗口的宽度
																	   win_height:  窗口的高度
																	   * \返回值
																	   * \标注    	当前光标的坐标为：mousex mousey
																	   *//*-------------------------------------------------------------------*/
void put_coordinate(float x, float y, float width, float height, float win_width, float win_height)
{
	GLfloat w_ratio = 0.1*win_width / 330;
	GLfloat h_ratio = 0.1*win_height / 220;

	//正投影
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-win_width / 2, win_width / 2, -win_height / 2, win_height / 2, 1.0, -1.0);
	//glViewport(0, 0, width, height);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glLoadIdentity();

	glTranslatef(x - win_width / 2, y - win_height / 2 + 110 * h_ratio, 0);           //将要显示的文字沿y轴上移110个单位(像素)
	glScalef(w_ratio, h_ratio, 1);                                        //将要显示的文字的x,y坐标大小缩放至原来的0.1倍

	//生成文字的固有大小为(330,220)
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 88);                   //输入要显示文字“X”的ASC码值
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 58);                   //输入要显示文字“:”的ASC码值
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 48 + (mousex / 100));      //计算当前位置的x坐标的百位数字
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 48 + (mousex % 100) / 10);   //计算当前位置的x坐标的十位数字
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 48 + (mousex % 10));       //计算当前位置的x坐标的个位数字

	glLoadIdentity();
	glTranslatef(x - win_width / 2, y - win_height / 2, 0);
	glScalef(w_ratio, h_ratio, 1);                                        //将要显示的文字的x,y坐标大小缩放至原来的0.1倍
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 89);                   //输入要显示文字“Y”的ASC码值
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 58);                   //输入要显示文字“:”的ASC码值
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 48 + (mousey / 100));      //计算当前位置的y坐标的百位数字
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 48 + (mousey % 100) / 10);   //计算当前位置的y坐标的十位数字
	glutStrokeCharacter(GLUT_STROKE_ROMAN, 48 + (mousey % 10));       //计算当前位置的y坐标的个位数字

	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glm::mat4 mv0 = glm::lookAt(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	glLoadMatrixf(&mv0[0][0]);
}

//将用户坐标映射到窗口坐标；当窗口的大小改变时，该函数会被调用
void ChangeSize(GLsizei w, GLsizei h)
{
	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//正投影
	GLfloat aspectRatio = (GLfloat)w / (GLfloat)(h + 1);           //h+1防止除数为零
	if (w <= h)
	{
		int windowWidth = 100;
		int windowHeight = 100 / aspectRatio;
		glOrtho(-100.0, 100.0, -windowHeight, windowHeight, 1.0, -1.0);
	}
	else
	{
		int windowWidth = 100 * aspectRatio;
		int windowHeight = 100;
		glOrtho(-windowWidth, windowWidth, -100.0, 100.0, 1.0, -1.0);
	}
}
//？？？

inline double average(vec3& v){
	return (v[0] + v[1] + v[2]) / 3.0;
}

void
print_vector(const vec3& v) {
	cout << v[0] << " " << v[1] << " " << v[2] << endl;
}

void
print_matrix(const mat3& m) {
	cout << m[0][0] << "|" << m[0][1] << "| " << m[0][2] << endl;
	cout << m[1][0] << "|" << m[1][1] << "| " << m[1][2] << endl;
	cout << m[2][0] << "|" << m[2][1] << "| " << m[2][2] << endl;
}

/** rotate the Z vector in the direction of norm */
mat3 rotate_axis(const vec3& sample, const vec3& reflected_dir) {
	double u1 = ((double)rand() / (double)RAND_MAX);
	double u2 = ((double)rand() / (double)RAND_MAX);
	double u3 = ((double)rand() / (double)RAND_MAX);
	vec3 randompoint = vec3(u1, u2, u3);
	vec3 u = glm::normalize(glm::cross(randompoint, reflected_dir));
	vec3 v = glm::normalize(glm::cross(reflected_dir, u));
	mat3 rot = mat3(u, v, reflected_dir);
	rot = glm::transpose(rot);
	return rot;
}

mat3 rotate_around_axis(vec3& rotation_axis, double theta, vec3& vec_in) {
	rotation_axis = glm::normalize(rotation_axis);
	double common_factor = sin(theta*0.5);
	double a = cos(theta*0.5);
	double b = rotation_axis[0] * common_factor;
	double c = rotation_axis[1] * common_factor;
	double d = rotation_axis[2] * common_factor;

	mat3 m = mat3(a*a + b*b - c*c - d*d, 2 * (b*c - a*d), 2 * (b*d + a*c),
		2 * (b*c + a*d), a*a - b*b + c*c - d*d, 2 * (c*d - a*b),
		2 * (b*d - a*c), 2 * (c*d + a*b), a*a - b*b - c*d + d*d);
	return m;

}

/* Sample a hemisphere f(r diffuse ray */
vec3 cos_weighted_hem(vec3& norm){
	double u1 = ((double)rand() / (double)RAND_MAX);
	double u2 = ((double)rand() / (double)RAND_MAX);

	vec3 y = vec3(norm);
	vec3 h = vec3(norm);
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

/* Sample a hemispehre for specular ray */
vec3 specular_weighted_hem(vec3& reflection, const vec3& norm, double n){
	double u1 = ((double)rand() / (double)RAND_MAX);
	double u2 = ((double)rand() / (double)RAND_MAX);
	vec3 y = vec3(reflection);
	vec3 h = vec3(reflection);
	double alpha = acos(pow(u1, 1.0 / (n + 1.0)));
	double phi = 2.0 * M_PI * u2;
	double xs = sin(alpha) * cos(phi);
	double ys = cos(alpha);
	double zs = sin(alpha) * sin(phi);
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

vec3 uniform_sample_hem(vec3& norm) {
	double u1 = ((double)rand() / (double)RAND_MAX);
	double u2 = ((double)rand() / (double)RAND_MAX);
	double phi = 2.0* M_PI * u2;
	double r = sqrt(1.0 - (u1 * u1));
	vec3 dir = vec3(cos(phi) * r, sin(phi) * r, u1);
	if (glm::dot(dir, norm) < 0)
		return -dir;
	return dir;
	//dir = rotate_axis(dir, norm);
	//return dir;
}




/* Shoot ray at scene */
vec3 findColor2(Scene* scene, Ray& ray, double weight, vector<Intersection> * pointlist) {


	/* Intersect scene */
	Intersection hit;
	scene->bvh->getIntersection(ray, &hit, false);
	if (!hit.primative) {
		return vec3(0, 0, 0); //background color
	}

	if (average(hit.primative->emission) > EPSILON4){
		return hit.primative->emission;
	}

#ifdef draw_debug

	if (pointlist != NULL)
	{
		pointlist->push_back(hit);
	}
#endif // draw_debug

	/* Russian Roulette */
	double russian = 1.0;
	const double cutoff = 0.25;
	if (weight < 0.001) {
		return vec3(0, 0, 0);
		double u1 = ((double)rand() / (double)RAND_MAX);
		if (u1 > cutoff) {
			return vec3(0, 0, 0);
		}
		else {
			russian = 1.0 / cutoff;
		}
	}

	vec3 color = vec3(0, 0, 0);
	bool insideflag = false;
	vec3 hitnormal = genpositiveNormal(hit, &insideflag);
	/*********************************************
	Add direct lighting contribution
	*********************************************/
	if (weight < 1.0 && weight > 0.01){
		int numLights = scene->lights.size();
		color += scene->lights[rand() % numLights]->shade(hit, scene->bvh, hitnormal, true);
	}
	weight *= 0.95; // make sure it doesnt go forever


	/*********************************************
	Importance sample global illumination
	*********************************************/
	double diffWeight = average(hit.primative->diffuse);
	double specWeight = average(hit.primative->specular);
	double refractWeight = hit.primative->refractivity;
	double threshold1 = diffWeight / (diffWeight + specWeight + refractWeight);
	double threshold2 = specWeight / (specWeight + refractWeight);

	vec3 normal = hit.primative->getNormal(hit);
	double rindex = hit.primative->indexofrefraction;

	if (refractWeight == 0)
	{
		threshold2 = 1.0;
	}
	else
	{
		double normal_test = glm::dot(normal, -ray.direction);
		if (normal_test < 0)
		{
			normal = -normal;
			rindex = 0.66;
			//cout << refractWeight<<endl;
		}
	}
	double u1 = ((double)rand() / ((double)RAND_MAX + 1));
	double u2 = ((double)rand() / ((double)RAND_MAX + 1));

	/* Importance sample on macro level to choose diffuse or specular */


	Ray newRay;
	if (u1 < threshold1) {//计算漫反射//反射率为零则不搞反射，折射率为零则不搞折射 || glm::length(hit.primative->specular)<EPSILON
		vec3 newDirection = cos_weighted_hem(normal);
		newRay = Ray(hit.point + 3 * EPSILON4*normal, newDirection);
		double prob = 1.0 / threshold1;
		color += prob * hit.primative->diffuse * findColor2(scene, newRay, weight*threshold1, pointlist);
	}
	else if (u2 <= threshold2){//计算镜面反射
		vec3 reflect = glm::reflect(ray.direction, normal);
		//newRay = Ray(hit.point+EPSILON*normal, reflect);
		//color += findColor(scene, newRay, specWeight*weight);

		vec3 newDirection = specular_weighted_hem(reflect, normal, hit.primative->shininess);
		//vec3 newDirection = uniform_sample_hem(normal);
		newRay = Ray(hit.point + 3 * EPSILON4*normal, newDirection);

		double dot = glm::dot(normal, newDirection);
		if (dot < 0.0) {
			return vec3(0, 0, 0);
		}

		double n = hit.primative->shininess;

		double multiplier = (n + 2.0) / (n + 1.0);
		multiplier *= 1.0 / ((1.0 - threshold1)*threshold2);
		//multiplier *= dot;	//封印
		color += multiplier * hit.primative->specular * findColor2(scene, newRay, (1.0 - threshold1)*threshold2*weight, pointlist);

	}
	else{//计算折射//失败原因：在球内部求不到交点，交点交不到本球
		//cout << "pig" << endl;

		//if (rindex!=1.5)
		//cout << rindex<<endl;
		assert(hit.primative->refractivity > 0);
		vec3 refract = glm::refract(ray.direction, normal, rindex);
		//cout << refract.x;
		{
			newRay = Ray(hit.point - 3 * EPSILON4*normal, refract);
			double multiplier = 1.0 / ((1.0 - threshold1)*(1 - threshold2));
			//vec3 absorbance = hit.primative->diffuse * hit.primative->absorb_constant * dist;
			//vec3 transparency = 
			double transparency = 0.9;
			color += transparency * multiplier * findColor2(scene, newRay, weight*(1.0 - threshold1)*(1 - threshold2), pointlist);
		}
		//color += findColor(scene, newRay, specWeight*weight);


	}

	return russian*color;
}

vec3 genExitSample(Intersection* it, Material* m, vec3 hitnormal, bool isInside)
{
	return m->genExitSample(&it->point, it->sourceDirection, hitnormal, isInside);//sourceDirection起点在交点， 是入射反向

}


vec3 genSpectrumIndex(Intersection it, vec3 exitent, vec3 hitnormal, bool isInside)
{

	return it.primative->material->getSpectrumIndex(it.point, it.sourceDirection, exitent, hitnormal, isInside, true);
}

vec3 genpositiveNormal(Intersection it, bool *insideflag)
{
	vec3 N = it.primative->getNormal(it);
	if (glm::dot(N, it.sourceDirection) < 0)
	{
		*insideflag = true;
		return -N;
	}
	else return N;
}

double Transmittance(Intersection hit, bool isInside)
{
	if (isInside) return exp(-hit.travel * hit.primative->material->extinct_inside);
	else return exp(-hit.travel * hit.primative->material->extinct_outside);
}
int times, nulltimes, errortimes = 0;
/* Shoot ray at scene */
vec3 findColor(Scene* scene, Ray& ray, double weight, vector<vec3> * pointlist, int recurdepth) {
	times++;
	//cout << "次数" << times<<endl;
	//已知进入场景的ray， 求交，返回颜色
	/* Intersect scene */
	Intersection hit;
	scene->bvh->getIntersection(ray, &hit, false);
	if (!hit.primative) {
		//cout << "无交正常返回次数" << ++nulltimes<<endl;
		return vec3(0, 0, 0); //background color
	}
	else{
		if (average(hit.primative->emission) > EPSILON4){
			return hit.primative->emission;
		}

		//----镜面反射：未采样到交点，也加光源

#ifdef draw_debug

		if (pointlist != NULL)
		{
			pointlist->push_back(hit.point);
		}
#endif // draw_debug

		/* Russian Roulette */
		double russian = 1.0;
		const double cutoff = 0.25;
		if (weight < 0.01) {
			return vec3(0, 0, 0);
			double u1 = ((double)rand() / (double)RAND_MAX);
			if (u1 > cutoff) {
				return vec3(0, 0, 0);
			}
			else {
				russian = 1.0 / cutoff;
			}
		}

		vec3 color = vec3(0, 0, 0);
		bool insideflag = false;
		vec3 hitnormal = genpositiveNormal(hit, &insideflag);
		/*********************************************
		Add direct lighting contribution
		*********************************************/
		//		if (weight == 1.0)
		//		{
		//			for (unsigned int i = 0; i < scene->lights.size(); i++) {
		//				color += scene->lights[i]->shade(hit, scene->bvh, hitnormal, false);
		//			}
		//		}
		if (weight < 1.0 && weight > 0.01){
			int numLights = scene->lights.size();
			color += scene->lights[rand() % numLights]->shade(hit, scene->bvh, hitnormal, true);
		}

		weight *= 0.95; // make sure it doesnt go forever
		//重要性采样交给material做
		//        Vector exitant = intersection.primitive->material->Sample(&point, incident, normal, wavelength, prng);

		if (hit.primative->material == NULL)
			cout << "无交错误返回次数" << ++errortimes << endl;
		vec3 exitent = genExitSample(&hit, hit.primative->material, hitnormal, insideflag);
		//反射的颜色值也交给material做
		//		float radiance = intersection.primitive->material->Reflectance(incident, exitant, normal, wavelength, true);
		vec3 SpectrumIndex = genSpectrumIndex(hit, exitent, hitnormal, insideflag);

		/*********************************************
		Add indirect lighting contribution
		*********************************************/
		if (pointlist != NULL)
		{
			pointlist->push_back(hit.point);
		}
		//发出下一条ray给递归//滚雪球过程//而pbrt（非递归）是个多光路求和过程，每次pathThroughput最终表示长度不等的一条连向光源的光路
		if (recurdepth <= 0){
			int numLights = scene->lights.size();
			return scene->lights[rand() % numLights]->shade(hit, scene->bvh, hitnormal, true);
		}
		else
		{
			color += Transmittance(hit, insideflag)* SpectrumIndex * findColor(scene, Ray(hit.point, exitent), weight, pointlist, recurdepth - 1);

			//absorbance 或称 体传播 或称 transmission
			//color *= Transmittance(hit, insideflag);		//render->Transmittance();//Intersection加flag表示在什么介质中？
			//返回颜色值

			return russian*color;
		}

	}
}

/* Main raytracing function. Shoots ray for each pixel with anialiasing */
/* ouputs bitmap to global variable*/
void raytrace(double rayscast) {
	double subdivisions = scene->antialias;
	double subdivide = 1.0 / subdivisions;

	double old_weight = rayscast / (rayscast + 1.0);
	double new_weight = 1.0 - old_weight;

#pragma omp parallel for
	for (int j = 0; j < scene->height; j++){
		int tid = omp_get_thread_num();
		if (tid == 0) {
			clog << "Progress: " << (j * 100 * omp_get_num_threads()) / scene->height << "%" << "\r";
		}
		RGBQUAD rgb;
		for (int i = 0; i < scene->width; i++) {
			vec3 color;
			for (double a = 0; a < subdivisions; a += 1) {
				for (double b = 0; b < subdivisions; b += 1) {
					double randomNum1 = ((double)rand() / (double)RAND_MAX) * subdivide;
					double randomNum2 = ((double)rand() / (double)RAND_MAX) * subdivide;
					Ray ray = scene->castEyeRay(i + (a*subdivide) + randomNum1, j + (b*subdivide) + randomNum2);
					color += findColor(scene, ray, 1.0, NULL, 5);
				}
			}
			color *= (subdivide * subdivide);

			color[0] = min(color[0], 1.0);
			color[1] = min(color[1], 1.0);
			color[2] = min(color[2], 1.0);

			pixels[i + scene->width*j] *= old_weight;
			pixels[i + scene->width*j] += new_weight*color;
			color = pixels[i + scene->width*j] + direct_pixels[i + scene->width*j];
			rgb.rgbRed = min(color[0], 1.0)*255.0;
			rgb.rgbGreen = min(color[1], 1.0)*255.0;
			rgb.rgbBlue = min(color[2], 1.0)*255.0;
			FreeImage_SetPixelColor(bitmap, i, j, &rgb);
		}
	}
}

/* Only the direct lighting */
vec3 direct_lighting(Scene* scene, Ray& ray){

	/* Intersect scene */
	Intersection hit;
	scene->bvh->getIntersection(ray, &hit, false);

	if (!hit.primative) {
		return vec3(0, 0, 0); //background color
	}
	if (glm::length(hit.primative->emission) > EPSILON4){
		return hit.primative->emission;
	}

	vec3 color = vec3(0, 0, 0);
	bool insideflag = false;
	vec3 hitnormal = genpositiveNormal(hit, &insideflag);
	for (unsigned int i = 0; i < scene->lights.size(); i++) {
		color += scene->lights[i]->shade(hit, scene->bvh, hitnormal, false);
	}
	return color;
}


/* Calculate the direct lighting seperatly */
void direct_raytrace() {
	double subdivisions = scene->antialias;
	double subdivide = 1.0 / subdivisions;

#pragma omp parallel for
	for (int j = 0; j < scene->height; j++){
		int tid = omp_get_thread_num();
		if (tid == 0) {
			clog << "Direct Lighting: " << (j * 100 * omp_get_num_threads()) / scene->height << "%" << "\r";
		}
		RGBQUAD rgb;
		for (int i = 0; i < scene->width; i++) {
			vec3 color;
			for (double a = 0; a < subdivisions; a += 1) {
				for (double b = 0; b < subdivisions; b += 1) {
					double randomNum1 = ((double)rand() / (double)RAND_MAX) * subdivide;
					double randomNum2 = ((double)rand() / (double)RAND_MAX) * subdivide;
					Ray ray = scene->castEyeRay(i + (a*subdivide) + randomNum1, j + (b*subdivide) + randomNum2);
					color += direct_lighting(scene, ray);
				}
			}
			color *= (subdivide * subdivide);
			color[0] = min(color[0], 1.0);
			color[1] = min(color[1], 1.0);
			color[2] = min(color[2], 1.0);


			direct_pixels[i + scene->width*j] = color;
			rgb.rgbRed = min(color[0], 1.0)*255.0;
			rgb.rgbGreen = min(color[1], 1.0)*255.0;
			rgb.rgbBlue = min(color[2], 1.0)*255.0;
			FreeImage_SetPixelColor(bitmap, i, j, &rgb);
		}
	}
	clog << "\n	done\n";
}

/* Everything below here is openGL boilerplate */

void reshape(int w, int h){
	width = w;
	height = h;
	glViewport(0, 0, w, h);
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	stringstream ss;
	switch (key){
	case 'l':
		update = 1;
		break;
	case 's':
		ss << scene->filename << "_" << int(rays_cast) << ".png";
		FreeImage_Save(FIF_PNG, bitmap, ss.str().c_str(), 0);
		cout << "Image saved!" << endl;
		break;
	case 'f':
		update = 10000;
		break;
	case 'h':
		update = 0;
		break;
	case 'r':
		glutReshapeWindow(scene->width, scene->height);
		break;
	case 27:  // Escape to quit
		FreeImage_DeInitialise();
		delete scene;
		delete pixels;
		delete direct_pixels;
		exit(0);
		break;

	}
	glutPostRedisplay();
}

void init(char* filename) {
	numFrames = 0;
	scene = new Scene(filename);
	rays_cast = 0.0;
	update = numFrames ? numFrames : false;
	width = scene->width;
	height = scene->height;
	pixels = new vec3[width*height];
	direct_pixels = new vec3[width*height];
	memset(pixels, 0, sizeof(vec3)*width*height);
	memset(direct_pixels, 0, sizeof(vec3)*width*height);

	FreeImage_Initialise();
	bitmap = FreeImage_Allocate(width, height, BPP);

	if (!numFrames) {
		vertexshader = initshaders(GL_VERTEX_SHADER, "shaders/vert.glsl");
		fragmentshader = initshaders(GL_FRAGMENT_SHADER, "shaders/frag.glsl");
		shaderprogram = initprogram(vertexshader, fragmentshader);

		glGenTextures(2, texture);

		//glEnable(GL_TEXTURE_2D) ;
		glBindTexture(GL_TEXTURE_2D, texture[0]);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		BYTE* bits = FreeImage_GetBits(bitmap);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene->width, scene->height,
			0, GL_BGR, GL_UNSIGNED_BYTE, (GLvoid*)bits);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-1, 1, -1, 1, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glm::mat4 mv = glm::lookAt(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
		glLoadMatrixf(&mv[0][0]);
		//glDisable(GL_TEXTURE_2D);
	}

}

void display(){
	if (!numFrames)
		glClear(GL_COLOR_BUFFER_BIT);

	int repetitions = numFrames;
	if (update){
		do {
			cout << "Iterations left: " << update << endl;
			time_t seconds = time(NULL);
			if (!rays_cast) {
				//direct_raytrace();
			}
			//direct_raytrace();
			raytrace(rays_cast);
			BYTE* bits = FreeImage_GetBits(bitmap);
			if (!numFrames){
				glBindTexture(GL_TEXTURE_2D, texture[0]);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene->width, scene->height,
					0, GL_BGR, GL_UNSIGNED_BYTE, (GLvoid*)bits);

			}
			rays_cast += 1.0;
			cout << "Number of Samples: " << rays_cast <<
				"\tTime: " << time(NULL) - seconds << " seconds" << endl;
			update -= 1;
			if (!numFrames) glutPostRedisplay();
			else repetitions -= 1;
		} while (repetitions);
	}
	if (numFrames && !update) {
		stringstream saved;
		saved << scene->filename << "_" << int(rays_cast) << ".png";
		FreeImage_Save(FIF_PNG, bitmap, saved.str().c_str(), 0);
		cout << "Image saved to " << saved.str() << "!" << endl;
		exit(0);
	}


	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glBegin(GL_QUADS);
	glTexCoord2d(0, 0); glVertex3d(-1, -1, 0);
	glTexCoord2d(0, 1); glVertex3d(-1, 1, 0);
	glTexCoord2d(1, 1); glVertex3d(1, 1, 0);
	glTexCoord2d(1, 0); glVertex3d(1, -1, 0);
	glEnd();

	//debug
	//drawray();
	put_coordinate(0, 0, 60, 40, scene->width, scene->height);
	draw_ray_debug(mousex, mousey);


	glutSwapBuffers();
}

void parse_command_line(int argc, char* argv[]) {
	for (int i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-f") == 0) {
			numFrames = atoi(argv[i + 1]);
			cout << "Running " << numFrames << " times" << endl;
			break;
		}
	}
}

int main(int argc, char* argv[]){
	if (argc < 2) {
		cerr << "You need at least 1 scene file as the argument" << endl;
		exit(1);
	}
	srand(time(0));
	parse_command_line(argc, argv);
	if (!numFrames) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutCreateWindow("Path Tracer");
		GLenum err = glewInit();
		if (GLEW_OK != err)
		{  /* Problem: glewInit failed, something is seriously wrong. */
			fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		}
	}
	init(argv[argc - 1]);
	if (numFrames)
		display();
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutReshapeWindow(width, height);
	//设置鼠标移动回调函数
	glutPassiveMotionFunc(mousemove);
	glutMotionFunc(mouseclick);
	glutMainLoop();
	system("pause");
	return 0;
}

void drawray()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(scene->camera_fovy, 1, 0, 100);//(GLdouble fovy, GLdouble aspect, GLdouble zNear//距离, GLdouble zFar//距离);
	glMatrixMode(GL_MODELVIEW);
	//	glm::mat4 mv1 = glm::lookAt(glm::vec3(0, 0.25, 0.2), glm::vec3(0, 0.25, 0), glm::vec3(0, 1, 0));
	glm::mat4 mv1 = glm::lookAt(scene->camera_eye, scene->camera_lookat, scene->camera_up);

	colortexel = FreeImage_Allocate(1, 1, BPP);
	RGBQUAD color;
	color.rgbBlue = 1.0 * 255;
	color.rgbRed = 0.0 * 255;
	color.rgbGreen = 0.0 * 255;
	FreeImage_SetPixelColor(colortexel, 0, 0, &color);
	BYTE* c = FreeImage_GetBits(colortexel);
	glBindTexture(GL_TEXTURE_2D, texture[1]);//texture0 是 颜色纹元，只用于设置颜色
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1,
		0, GL_BGR, GL_UNSIGNED_BYTE, c);

	glLoadMatrixf(&mv1[0][0]);
	glBegin(GL_LINE_LOOP);
	//	glVertex3d(-0.5, -0.1, -0.5);
	//	glVertex3d(-0.5, -0.1, 0.2);
	//	glVertex3d(0.5, -0.1, 0.2);
	//	glVertex3d(0.5, -0.1, -0.5);
	glVertex3d(552.8, 0.0, 0.0);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, 559.2);
	glVertex3d(549.6, 0.0, 559.2);

	glEnd();

	//没有找到 push pop方法
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glm::mat4 mv0 = glm::lookAt(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	glLoadMatrixf(&mv0[0][0]);

}

void draw_ray_debug(int i, int j)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(scene->camera_fovy, 1, 0, 100);//(GLdouble fovy, GLdouble aspect, GLdouble zNear//距离, GLdouble zFar//距离);
	glMatrixMode(GL_MODELVIEW);
	//	glm::mat4 mv1 = glm::lookAt(glm::vec3(0, 0.25, 0.2), glm::vec3(0, 0.25, 0), glm::vec3(0, 1, 0));
	glm::mat4 mv1 = glm::lookAt(scene->camera_eye, scene->camera_lookat, scene->camera_up);
	glLoadMatrixf(&mv1[0][0]);
	colortexel = FreeImage_Allocate(1, 1, BPP);

	vector<vec3> pointlist;
	cout << "坐标" << " " << i << " " << j << endl;
	//for (){
	RGBQUAD rgb;
	vec3 color;
	double randomNum1 = ((double)rand() / (double)RAND_MAX);
	double randomNum2 = ((double)rand() / (double)RAND_MAX);
	//Ray ray = scene->castEyeRay(i + randomNum1, j + randomNum2);
	Ray ray = scene->castEyeRay(i, j);
	color += findColor(scene, ray, 1.0, &pointlist, 20);

	//color[0] = min(color[0], 1.0);
	//color[1] = min(color[1], 1.0);
	//color[2] = min(color[2], 1.0);

	rgb.rgbRed = min(color[0], 1.0)*255.0;
	rgb.rgbGreen = min(color[1], 1.0)*255.0;
	rgb.rgbBlue = min(color[2], 1.0)*255.0;

	rgb.rgbBlue = 1.0 * 255;
	rgb.rgbRed = 0.0 * 255;
	rgb.rgbGreen = 0.0 * 255;
	FreeImage_SetPixelColor(colortexel, 0, 0, &rgb);

	BYTE* c = FreeImage_GetBits(colortexel);
	glBindTexture(GL_TEXTURE_2D, texture[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1,
		0, GL_BGR, GL_UNSIGNED_BYTE, c);

	glBegin(GL_LINE_LOOP);
	//	glVertex3d(-0.5, -0.1, -0.5);
	//	glVertex3d(-0.5, -0.1, 0.2);
	//	glVertex3d(0.5, -0.1, 0.2);
	//	glVertex3d(0.5, -0.1, -0.5);
	glVertex3d(552.8, 0.0, 0.0);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, 559.2);
	glVertex3d(549.6, 0.0, 559.2);

	glEnd();
	glBegin(GL_LINE_STRIP);
	//	glBegin(GL_LINE_LOOP);
	vec3 p = scene->castpoint(i, j);
	glVertex3d(p.x, p.y, p.z);
	cout << "开始";
	cout << p.x << " " << p.y << " " << p.z << endl;
	cout << pointlist.size() << endl;
	for (int ii = 0; ii < pointlist.size(); ii++, ii++)
	{
		glVertex3d(pointlist[ii].x, pointlist[ii].y, pointlist[ii].z);
		cout << "hit:" << pointlist[ii].x << " " << pointlist[ii].y << " " << pointlist[ii].z
			<< "origin:" << pointlist[ii + 1].x << " " << pointlist[ii + 1].y << " " << pointlist[ii + 1].z
			<< endl;
	}
	glEnd();


	//没有找到 push pop方法
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glm::mat4 mv0 = glm::lookAt(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	glLoadMatrixf(&mv0[0][0]);

}