/******************************************************************
*
* PathTracing.cpp
*
* Description: This program demonstrates global illumination rendering
* based on the path tracing method. The intergral in the rendering
* equation is approximated via Monte-Carlo integration; explicit 
* direct lighting is included to improve quality; the rendered image 
* is saved in PPM format.
*
* The code is largely based on the software smallpt by Kevin Beason,
* released under the MIT License.
*
* Advanced Computer Graphics Proseminar WS 2015
* 
* Interactive Graphics and Simulation Group
* Institute of Computer Science
* University of Innsbruck
*
*******************************************************************/

//object loader
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* source
#include "tiny_obj_loader.h"

/* Standard includes */
#include <cmath>   
#include <cstdlib> 
#include <iostream>
#include <fstream>
#ifdef _WIN32
#include "drand48.h"
#define M_PI    3.14159265358979323846
#define M_1_PI  0.31830988618379067154
#define M_4_PI  1.27323954473516268615107010698011488
#endif

using namespace std;

/*------------------------------------------------------------------
| Struct for standard Vector operations in 3D 
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
struct Vector 
{        
    double x, y, z;           /* Position XYZ or color RGB */

    Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
    Vector(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    
    Vector operator+(const Vector &b) const 
    {
        return Vector(x + b.x, y + b.y, z + b.z);
    }

    Vector operator-(const Vector &b) const
    {
        return Vector(x - b.x, y - b.y, z - b.z);
    }

    Vector operator/(double c) const 
    {
        return Vector(x / c, y / c, z / c);
    }

    Vector operator*(double c) const 
    {
        return Vector(x * c, y * c, z * c);
    }

    friend Vector operator*(double c, const Vector &b) 
    { 
        return b * c; 
    }

    Vector MultComponents(const Vector &b) const
    {
        return Vector(x * b.x, y * b.y, z * b.z);
    }

    const Vector Normalized() const
    {
        return Vector(x, y, z) / sqrt(x*x + y*y + z*z);
    }

    const double Dot(const Vector &b) const 
    {
        return x * b.x + y * b.y + z * b.z;
    }

    const Vector Cross(const Vector &b) const
    {
        return Vector((y * b.z) - (z * b.y), 
                      (z * b.x) - (x * b.z), 
                      (x * b.y) - (y * b.x));
    }
     
    const double Max() 
    {
        return fmax(x, fmax(x, y));
    }

    Vector& clamp() 
    {
        x = x<0 ? 0.0 : x>1.0 ? 1.0 : x;
        y = y<0 ? 0.0 : y>1.0 ? 1.0 : y;
        z = z<0 ? 0.0 : z>1.0 ? 1.0 : z;
        return *this;   
    }
};

typedef Vector Color;
const Color BackgroundColor(0.0, 0.0, 0.0);

/*------------------------------------------------------------------
| Structure for rays (e.g. viewing ray, path tracing)
------------------------------------------------------------------*/
struct Ray 
{
    Vector org, dir;    /* Origin and direction */
    Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};


/*------------------------------------------------------------------
| Struct holds pixels/colors of rendered image
------------------------------------------------------------------*/
struct Image 
{
    int width, height;
    Color *pixels;

    Image(int _w, int _h) : width(_w), height(_h) 
    {
        pixels = new Color[width * height];        
    }

    Color getColor(int x, int y) 
    {
        int image_index = (height-y-1) * width + x; 
        return pixels[image_index];
    }

    void setColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x; 
        pixels[image_index] = c;
    }

    void addColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x; 
        pixels[image_index] = pixels[image_index] + c;
    }

    int toInteger(double x)
    { 
        /* Clamp to [0,1] */
        if (x<0.0)
            x = 0.0;        
        
        if (x>1.0)
            x = 1.0;             

        /* Apply gamma correction and convert to integer */
        return int(pow(x,1/2.2)*255+.5); 
    }

    void Save(const string &filename) 
    {
        /* Save image in PPM format */
        FILE *f = fopen(filename.c_str(), "wb");
        fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
        for (int i = 0; i < width * height; i++)
            fprintf(f,"%d %d %d ", toInteger(pixels[i].x), 
                                   toInteger(pixels[i].y), 
                                   toInteger(pixels[i].z));
        fclose(f);
    }
};



/*------------------------------------------------------------------
| Scene objects are spheres or triangle spheres; material either 
| perfectly diffuse, specular (mirror reflection), transparent 
| (refraction/reflection), glossy or translucent
| (DIFFuse, SPECular, REFRactive, GLOSSY, TRANSlucent)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR, GLOSSY, TRANS }; 

struct Geom 
{
    Geom(Vector position_, Vector emission_, 
           Vector color_, Refl_t refl_, double glossy_factor, double trans_factor): position(position_), emission(emission_), 
           color(color_), refl(refl_), glossiness_factor(glossy_factor), translucency_factor(trans_factor) {}
    Vector position;
    Color emission, color;      
    Refl_t refl;
    double glossiness_factor; // between 0.0 (very glossy) and 1.0 (mirror)
	double translucency_factor; // between 0.0 (translucent) and 1.0 (not translucent)

    virtual Vector getNormalAtPoint(Vector point) const = 0;
    virtual double Intersect(const Ray &ray) const = 0;
};

struct Sphere : Geom
{
    double radius;
    
    Sphere(double radius_, Vector position_, Vector emission_, 
           Vector color_, Refl_t refl_, double glossy_factor = 0.5, double trans_factor = 0.5): 
		   Geom(position_, emission_, color_, refl_, glossy_factor, trans_factor), radius(radius_) {}

    Vector getNormalAtPoint(Vector p) const {
        return (p - position).Normalized();
    }

    double Intersect(const Ray &ray) const 
    { 
        /* Check for ray-sphere intersection by solving for t:
            t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
        Vector op = position - ray.org; 
        double eps = 1e-4;
        double b = op.Dot(ray.dir);
        double radicant = b*b - op.Dot(op) + radius*radius;
        if (radicant < 0.0) 
            return 0.0;      /* No intersection */
        else   
            radicant = sqrt(radicant);
        
        double t;
        t = b - radicant;    /* Check smaller root first */
        if(t > eps)
            return t;
        
        t = b + radicant;
        if(t > eps)          /* Check second root */
            return t;
        
        return 0.0;          /* No intersection in ray direction */  
    }
};

struct Triangle : Geom
{
	Vector p0; //origin               
	Vector edge_a, edge_b, edge_c; //edges making up the triangle
    Vector normal;

	//create a triangle, diagonal edge is calculated from 2 given edges
	Triangle(Vector p0_, Vector a_, Vector b_, Vector normal_, Color emission_, Color color_, Refl_t refl_, double glossy_factor = 0.5, double trans_factor = 0.5):
    Geom((p0 + (p0+edge_a) + (p0+edge_b)) / 3, emission_, color_, refl_, glossy_factor, trans_factor), p0(p0_), edge_a(a_), edge_b(b_), normal(normal_)
    {
		edge_c = edge_a - edge_b; //diagonal edge from point b to point a
		/*normal = edge_a.Cross(edge_b);
		normal = normal.Normalized();*/
	}

    Vector getNormalAtPoint(Vector p) const {
        return normal;
    }

	/* Triangle-ray intersection 
	 * Moeller-Trumbore approach
	*/
	double Intersect(const Ray &ray) const 
    {
		double epsilon = 0.0000001;

		Vector d = ray.dir;
		Vector p = ray.org;
		Vector s = p - p0;

		//ray and triangle are parallel if determinant is close to 0
		double detT = d.Cross(edge_b).Dot(edge_a);

		if (fabs(detT) < epsilon) {
			return 0.0;
		}

		//check ray-plane instersection
		double inv_detT = 1 / detT;
		double t = (s.Cross(edge_a).Dot(edge_b)) * inv_detT;

		if (t <= epsilon)
			return 0.0;

		//check if intersection is within triangle
		double bary_x = -1*(d.Cross(s).Dot(edge_b)) * inv_detT;

		if (bary_x < 0 || bary_x > 1) {
			return 0.0;
		}

		double bary_y = -1*(d.Cross(edge_a).Dot(s)) * inv_detT;

		if (bary_y < 0 || bary_x + bary_y > 1) {
			return 0.0;
		}

		return t;
	}
};

/******************************************************************
* Hard-coded scene definition: the geometry is composed of spheres
* (i.e. Cornell box walls are part of very large spheres). 
* These are defined by:
* - radius, center 
* - emitted light (light sources), surface reflectivity (~color), 
*   material
*******************************************************************/
Sphere spheres[] = 
{
    Sphere( 1e5, Vector( 1e5  +1,      40.8,      81.6),  Vector(), Vector(.75,.25,.25), DIFF), /* Left wall */
    Sphere( 1e5, Vector(-1e5 +99,      40.8,      81.6),  Vector(), Vector(.25,.25,.75), DIFF), /* Right wall */
    Sphere( 1e5, Vector(      50,      40.8,       1e5),  Vector(), Vector(.75,.75,.75), DIFF), /* Back wall */
    Sphere( 1e5, Vector(      50,      40.8, -1e5 +170),  Vector(), Vector(),            DIFF), /* Front wall */
    Sphere( 1e5, Vector(      50,       1e5,      81.6),  Vector(), Vector(.75,.75,.75), DIFF), /* Floor */
    Sphere( 1e5, Vector(      50,-1e5 +81.6,      81.6),  Vector(), Vector(.75,.75,.75), DIFF), /* Ceiling */

    Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1,1,1)*.999,  SPEC), /* Mirror sphere */
    Sphere(10.5, Vector(50, 16.5, 105), Vector(), Vector(1,1,1)*.999,  GLOSSY, 0.8), /* Glossy sphere */
	Sphere(10.5, Vector(10.5, 50.0, 105), Vector(), Vector(1,1,1)*.999,  TRANS, 1, 0.8), /* Transluscent sphere */
    Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1,1,1)*.999,  REFR), /* Glass sphere */

    Sphere( 1.5, Vector(50, 81.6-16.5, 81.6), Vector(4,4,4)*100, Vector(), DIFF), /* Light */
};

vector<Triangle> triangles;

/* Load the triangle geometry */
void loadObjects() {
    string inputfile = "pyramid_rotated.obj";
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;
    string err;
    bool ret = tinyobj::LoadObj(shapes, materials, err, inputfile.c_str());
    if (!err.empty()) { // `err` may contain warning message.
      std::cerr << err << std::endl;
    }
    if (!ret) {
      exit(1);
    }
    
    vector<unsigned int> indices = shapes[0].mesh.indices;

    //obj loader has trouble to parse both vertex and normal indices at the same time, hence we provided 2 separate obj files
    inputfile = "pyramid_rotated_with_normals.obj";
    ret = tinyobj::LoadObj(shapes, materials, err, inputfile.c_str());
    if (!err.empty()) { // `err` may contain warning message.
      std::cerr << err << std::endl;
    }
    if (!ret) {
      exit(1);
    }

    Vector obj_center = Vector(48, 36.5, 37);

    for (int i = 0; i < int(shapes.size()); i++) {
        for (int j = 0; j < int(shapes[i].mesh.indices.size()); j+=3) {
            int i1 = indices[j];
            int i2 = indices[j+1];
            int i3 = indices[j+2];
            Vector v1 = Vector(shapes[i].mesh.positions[i1], shapes[i].mesh.positions[i1+1], shapes[i].mesh.positions[i1+2])*16.5 + obj_center;
            Vector v2 = Vector(shapes[i].mesh.positions[i2], shapes[i].mesh.positions[i2+1], shapes[i].mesh.positions[i2+2])*16.5 + obj_center; 
            Vector v3 = Vector(shapes[i].mesh.positions[i3], shapes[i].mesh.positions[i3+1], shapes[i].mesh.positions[i3+2])*16.5 + obj_center;
            Vector normal = Vector(shapes[i].mesh.normals[j*3], shapes[i].mesh.normals[j*3+1], shapes[i].mesh.normals[j*3+2]);
            triangles.push_back(Triangle(v1, v2-v1, v3-v1, normal, Vector(), Vector(.25,.75,.25), DIFF)); 
        }
    }
}

vector<Geom *> geoms;

void init_geom_array() {
    int ns = sizeof(spheres) / sizeof(Sphere);
    int ts = triangles.size();
    for(int i = 0; i < ns; i++) {
        geoms.push_back(&spheres[i]);
    }
    for(int i = 0; i < ts; i++) {
        geoms.push_back(&triangles[i]); 
    }
}

/******************************************************************
* Check for closest intersection of a ray with the scene;
* returns true if intersection is found, as well as ray parameter
* of intersection and id of intersected object
*******************************************************************/
bool Intersect(const Ray &ray, double &t, int &id)
{
    const int n = geoms.size();
    t = 1e20;

    for (int i = 0; i < n; i++) 
    {
        double d = geoms[i]->Intersect(ray);
        if (d > 0.0  && d < t) 
        {
            t = d;
            id = i;
        }
    }
    
    return t < 1e20;
}

/*
 * return a random vector around axis
 *
 * n (0.0 to 1.0) determines how close the vector is to the axis vector (1.0 is the same)
 */
Vector getSample(Vector axis, double n) {
	double e1 = drand48();
	double e2 = drand48();
	
	// Set up local orthogonal coordinate system u,v,w on surface
	Vector w = axis; 
	Vector u;
	Vector r (drand48(), drand48(), drand48());
	while (fabs(r.Dot(w)) < 0.000000001) {
		r = Vector(drand48(), drand48(), drand48());
	}
	u = (r.Cross(w)).Normalized();

	Vector v = w.Cross(u);
		
	// calc sample vector
	double phi = 2.0 * M_PI * e1;
	double cosT = pow(e2, 1/(n+1));

	double z = cosT + (1 - cosT) * n; // determine angle spread depended on n
	double sinT = sqrt(1 - z * z);
	double x = cos(phi) * sinT;
	double y = sin(phi) * sinT;

	Vector sample (x,y,z);

	//translate sample into correct space
	Vector utrans(u.x, v.x, w.x);
	Vector vtrans(u.y, v.y, w.y);
	Vector wtrans(u.z, v.z, w.z); 
	Vector sampleDirection(utrans.x * sample.x + utrans.y * sample.y + utrans.z * sample.z,
	                       vtrans.x * sample.x + vtrans.y * sample.y + vtrans.z * sample.z,
	                       wtrans.x * sample.x + wtrans.y * sample.y + wtrans.z * sample.z);
	return sampleDirection;
}

/******************************************************************
* Recursive path tracing for computing radiance via Monte-Carlo
* integration; only considers perfectly diffuse, specular or 
* transparent materials; 
* after 5 bounces Russian Roulette is used to possibly terminate rays;  
* emitted light from light source only included on first direct hit 
* (possibly via specular reflection, refraction), controlled by 
* parameter E = 0/1;  
* on diffuse surfaces light sources are explicitely sampled;
* for transparent objects, Schlick's approximation is employed;
* for first 3 bounces obtain reflected and refracted component,
* afterwards one of the two is chosen randomly   
*
* samples should be kept low to avoid segfault (stack full, too much recursions)
*******************************************************************/
Color Radiance(const Ray &ray, int depth, int E)
{
    depth++;

    int numSpheres = int(sizeof(spheres) / sizeof(Sphere));

    double t;                               
    int id = 0; 
                             
    if (!Intersect(ray, t, id))   /* No intersection with scene */
        return BackgroundColor; 

    const Geom* obj = geoms[id];

    Vector hitpoint = ray.org + ray.dir * t;    /* Intersection point */
    Vector normal = obj->getNormalAtPoint(hitpoint);  /* Normal at intersection */ 
    Vector nl = normal;

    /* Obtain flipped normal, if object hit from inside */
    if (normal.Dot(ray.dir) >= 0) 
        nl = nl*-1.0;

    Color col = obj->color; 

    /* Maximum RGB reflectivity for Russian Roulette */
    double p = col.Max();

    if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
    {
        if (drand48() < p)            /* Russian Roulette */
            col = col * (1/p);        /* Scale estimator to remain unbiased */
        else 
            return obj->emission * E;  /* No further bounces, only return potential emission */
    }

    if (obj->refl == DIFF)
    {                  
        /* Compute random reflection vector on hemisphere */
        double r1 = 2.0 * M_PI * drand48(); 
        double r2 = drand48(); 
        double r2s = sqrt(r2); 
        
        /* Set up local orthogonal coordinate system u,v,w on surface */
        Vector w = nl; 
        Vector u;
        
        if(fabs(w.x) > .1)
            u = Vector(0.0, 1.0, 0.0);
        else
            u = (Vector(1.0, 0.0, 0.0).Cross(w)).Normalized(); 

        Vector v = w.Cross(u);  

        /* Random reflection vector d */
        Vector d = (u * cos(r1) * r2s + 
                    v * sin(r1) * r2s + 
                    w * sqrt(1 - r2)).Normalized();  

        /* Explicit computation of direct lighting */
        Vector e;
        for (int i = 0; i < numSpheres; i ++)
        {
            const Sphere &sphere = spheres[i];
            if (sphere.emission.x <= 0 && sphere.emission.y <= 0 && sphere.emission.z <= 0) 
                continue; /* Skip objects that are not light sources */
      
            /* Randomly sample spherical light source from surface intersection */

            /* Set up local orthogonal coordinate system su,sv,sw towards light source */
            Vector sw = sphere.position - hitpoint;
            Vector su;
            
            if(fabs(sw.x) > 0.1)
                su = Vector(0.0, 1.0, 0.0);
            else
                su = Vector(1.0, 0.0, 0.0);

            su = (su.Cross(w)).Normalized();
            Vector sv = sw.Cross(su);

            /* Create random sample direction l towards spherical light source */
            double cos_a_max = sqrt(1.0 - sphere.radius * sphere.radius / 
                                    (hitpoint - sphere.position).Dot(hitpoint-sphere.position));
            double eps1 = drand48();
            double eps2 = drand48();
            double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
            double sin_a = sqrt(1.0 - cos_a * cos_a);
            double phi = 2.0*M_PI * eps2;
            Vector l = su * cos(phi) * sin_a + 
                       sv * sin(phi) * sin_a + 
                       sw * cos_a;
            l = l.Normalized();

            /* Shoot shadow ray, check if intersection is with light source */
            if (Intersect(Ray(hitpoint,l), t, id) && id==i)
            {  
                double omega = 2*M_PI * (1 - cos_a_max);

                /* Add diffusely reflected light from light source; note constant BRDF 1/PI */
                e = e + col.MultComponents(sphere.emission * l.Dot(nl) * omega) * M_1_PI; 
            }
        }
   
        /* Return potential light emission, direct lighting, and indirect lighting (via
           recursive call for Monte-Carlo integration */      
        return obj->emission * E + e + col.MultComponents(Radiance(Ray(hitpoint,d), depth, 0));
    } 
    else if (obj->refl == SPEC) 
    {  
        /* Return light emission mirror reflection (via recursive call using perfect
           reflection vector) */
        return obj->emission + 
            col.MultComponents(Radiance(Ray(hitpoint, ray.dir - normal * 2 * normal.Dot(ray.dir)),
                               depth, 1));
    }
	else if (obj->refl == GLOSSY) 
    { 
		Vector perfectReflectionDirectionN = (ray.dir - normal * 2 * normal.Dot(ray.dir)).Normalized();

        // shoot secondary rays, perturbed around perfect reflection ray
		// make sure that stack has enough space for recursions by limiting depth
		if (depth < 3) {
			int num_samples = 4;
			Color avrg(0.0,0.0,0.0);
			for (int i = 0; i < num_samples; i++) {
				Color rad (Radiance(Ray(hitpoint, getSample(perfectReflectionDirectionN, obj->glossiness_factor)), depth, 1));
				avrg.x += rad.x;
				avrg.y += rad.y;
				avrg.z += rad.z;
			}
			avrg = avrg/num_samples;
			return obj->emission + col.MultComponents(avrg);
		}
		//could shoot more rays, but the result is not worth the additional rendering time
		return obj->emission;
    }

    /* Otherwise object transparent or translucent, i.e. assumed dielectric glass material */
    Ray reflRay (hitpoint, ray.dir - normal * 2 * normal.Dot(ray.dir));  /* Prefect reflection */  
    bool into = normal.Dot(nl) > 0;       /* Bool for checking if ray from outside going in */
    double nc = 1;                        /* Index of refraction of air (approximately) */  
    double nt = 1.5;                      /* Index of refraction of glass (approximately) */
    double nnt;

    if(into)      /* Set ratio depending on hit from inside or outside */
        nnt = nc/nt;
    else
        nnt = nt/nc;

    double ddn = ray.dir.Dot(nl);
    double cos2t = 1 - nnt * nnt * (1 - ddn*ddn);

    /* Check for total internal reflection, if so only reflect */
    if (cos2t < 0)
        return obj->emission + col.MultComponents( Radiance(reflRay, depth, 1));

    /* Otherwise reflection and/or refraction occurs */
    Vector tdir;

    /* Determine transmitted ray direction for refraction */
    if(into)
        tdir = (ray.dir * nnt - normal * (ddn * nnt + sqrt(cos2t))).Normalized();
    else
        tdir = (ray.dir * nnt + normal * (ddn * nnt + sqrt(cos2t))).Normalized();

    /* Determine R0 for Schlick's approximation */
    double a = nt - nc;
    double b = nt + nc;
    double R0 = a*a / (b*b);
  
    /* Cosine of correct angle depending on outside/inside */
    double c;
    if(into)
        c = 1 + ddn;
    else
        c = 1 - tdir.Dot(normal);

    /* Compute Schlick's approximation of Fresnel equation */ 
    double Re = R0 + (1 - R0) *c*c*c*c*c;   /* Reflectance */
    double Tr = 1 - Re;                     /* Transmittance */

    /* Probability for selecting reflectance or transmittance */
    double P = .25 + .5 * Re;
    double RP = Re / P;         /* Scaling factors for unbiased estimator */
    double TP = Tr / (1 - P);

	double translucency = obj->translucency_factor;
    if (depth < 3) {   /* Initially both reflection and transmission */
		//shoot multiple rays perturbed around refraction ray if object is translucent
		if (obj->refl == TRANS) {
			int num_samples = 2;
			Color avrg(0.0,0.0,0.0);
			for (int i = 0; i < num_samples; i++) {
				Color rad (Radiance(Ray(hitpoint, getSample(tdir, translucency)), depth, 1)*Tr);
				avrg.x += rad.x;
				avrg.y += rad.y;
				avrg.z += rad.z;
			}
			avrg = avrg/num_samples;
			return obj->emission + col.MultComponents(Radiance(reflRay, depth, 1) * Re + avrg);
		}
        return obj->emission + col.MultComponents(Radiance(reflRay, depth, 1) * Re + 
		                                          Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
	}
	//either shoot reflective or transmissive ray, decided by russian roulette
    else {
        if (drand48() < P) {
            return obj->emission + col.MultComponents(Radiance(reflRay, depth, 1) * RP);
		}
        else {
			//shoot multiple rays perturbed around refraction ray if object is translucent
			if (obj->refl == TRANS) {
				int num_samples = 2;
				Color avrg(0.0,0.0,0.0);
				for (int i = 0; i < num_samples; i++) {
					Color rad (Radiance(Ray(hitpoint, getSample(tdir, translucency)), depth, 1)*TP);
					avrg.x += rad.x;
					avrg.y += rad.y;
					avrg.z += rad.z;
				}
				avrg = avrg/num_samples;
				return obj->emission + col.MultComponents(avrg);
			}
			return obj->emission + col.MultComponents(Radiance(Ray(hitpoint,tdir), depth, 1) * TP);
		}
	}
            
}


/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height 
* - Number of samples per subpixel (non-uniform filtering): samples 
* Rendered result saved as PPM image file
*******************************************************************/
 
int main(int argc, char *argv[]) 
{
    // PathTracing [samples per subpixel] [(sqrt of) samples on lens] [focal distance|cm] [aperture|f-stops]
    
    //assumption: values are given in cm (somewhat realistic size of box)   
    /* Default values */
    int width = 1024;
    int height = 768;
    int samples = 4;
    int lens_samples = 1;
    double focal_distance = 248.6;  //focused on metal sphere (217.6 would focus on glass sphere)   
    double image_distance = 1;  //for our virtual camera the image sensor is 1 unit away from the lens   
    double fstops = .2; // smaller f-stop -> shallower depth of field (here: ridiculously low value to demonstrate depth of field)

    /* User input */
    if(argc >= 2)
        samples = atoi(argv[1]);
    if(argc >= 3)
        lens_samples = atoi(argv[2]);
    if(argc == 5) {
        focal_distance = atof(argv[3]);
        fstops = atof(argv[4]);
    }
    
    /* Computed values */
    double focal_length = 1 / (1/image_distance + 1/focal_distance);    //focus length for given distances of focal & image plane to the lens
     /* Adapt behaviour of fstops to those of a real camera with focal length of 100mm (where distance to image sensor is not exactly 10mm but dependent on the focal distance) */
    fstops *= focal_length / 10;
    double aperture = 0/*focal_length/fstops*/ /*temporarily disabled depth of field for development convenience*/;
     
    /* Set camera origin and viewing direction (negative z direction) */
    Ray camera(Vector(50.0, 52.0, 170 + 125.6), Vector(0.0, -0.042612, -1.0).Normalized());

    /* Image edge vectors for pixel sampling */
    Vector cx = Vector(width * 0.5135 / height);    // cx = aspect ratio * .5, e.g. 1.2 * .5 (-> so that cx/cy reach end of column/row after the same number of steps)
    Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135; // cy = 1 * .5

    /* Final rendering */
    Image img(width, height);

    loadObjects();
    init_geom_array();

    /* Loop over image rows */
    for (int y = 0; y < height; y ++) 
    { 
        cout << "\rRendering (" << samples * 4 << " spp) " << (100.0 * y / (height - 1)) << "%     ";
        srand(y * y * y);
 
        /* Loop over row pixels */
        for (int x = 0; x < width; x ++)  
        {
            img.setColor(x, y, Color());
 
            /* 2x2 subsampling per pixel */
            for (int sy = 0; sy < 2; sy ++) 
            {
                for (int sx = 0; sx < 2; sx ++) 
                {
                    Color accumulated_radiance = Color();

                    /* Compute radiance at subpixel using multiple samples */
                    for (int s = 0; s < samples; s ++) 
                    {
                        const double r1 = 2.0 * drand48();
                        const double r2 = 2.0 * drand48();

                        /* Transform uniform into non-uniform filter samples */
                        double dx;
                        if (r1 < 1.0)
                            dx = sqrt(r1) - 1.0;
                        else
                            dx = 1.0 - sqrt(2.0 - r1);

                        double dy;
                        if (r2 < 1.0)
                            dy = sqrt(r2) - 1.0;
                        else
                            dy = 1.0 - sqrt(2.0 - r2);
        
                        /* Intersection with focal plane from camera center */
                        Vector focal_point = camera.org +
                                   (cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                                   cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) + 
                                   camera.dir * image_distance) * focal_distance;   //focal "plane" is actually curved                       
                        
                        /* Stratified camera lens sampling, using a concentric map */
                        for(int i = 0; i < lens_samples; i++)
                        {
                            for(int j = 0; j < lens_samples; j++)
                            { 
                                const double lsx = (drand48()+i)/lens_samples - 0.5;
                                const double lsy = (drand48()+j)/lens_samples - 0.5;
                                double r, phi;
                                if(lsx > lsy && lsx > -lsy) {
                                    r = lsx;
                                    phi = lsy/lsx * M_PI_4;
                                }
                                else if(lsx > lsy && lsx < -lsy) {
                                    r = -lsy;
                                    if (lsy == 0)
                                        phi = 0;
                                    else
                                        phi = (6 - lsx/lsy) * M_PI_4;
                                }
                                else if(lsx < lsy && lsx > -lsy) {
                                    r = lsy;
                                    phi = (2 - lsx/lsy) * M_PI_4;
                                }
                                else {
                                    r = -lsx;
                                    phi = (4 + lsy/lsx) * M_PI_4;
                                }
                                //lens is slightly elliptical to match image aspect ratio (use [cx/cy].Normalized() for a circular lens)
                                Vector offset = r * cos(phi) * cx + r * sin(phi) * cy;
                                offset = offset * aperture;                        
            
                                Vector start = camera.org + offset;
                                Vector dir = focal_point - start;
                                dir = dir.Normalized();
                                start = start + dir * ((camera.org.z-169)/abs(dir.z)); // Extend camera ray to start inside box (camera positioned at z~300, front wall at z~170)

                                /* Accumulate radiance */
                                accumulated_radiance = accumulated_radiance + 
                                    Radiance( Ray(start, dir), 0, 1) / samples / (lens_samples*lens_samples);
                            }
                        }
                    } 
                    
                    accumulated_radiance = accumulated_radiance.clamp() * 0.25;

                    img.addColor(x, y, accumulated_radiance);
                }
            }
        }
    }
    cout << endl;

    img.Save(string("image.ppm"));

}
