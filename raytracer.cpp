#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <thread>
// #include <string>

//exclude before submitting 
// #include <random>
// #include <ctime>
// #include <iomanip>
//exclude before submitting 

#define MAX_T 99999.0f
#define MAX_I 999999

// #include <fstream>


using namespace parser;

enum Obj {SPHERE, TRIANGLE, MESH_TRIANGLE};
struct IntersectionResult
{
    Material m;
    Vec3f n;
    int minI;
    int minJ;
    double minT;
    char obj_type;
    bool intersected;
    

};
struct Ray
{
    Vec3f o, d;
    Ray(Vec3f o, Vec3f d)
    {
        this-> o = o;
        this-> d = d;
    }
    Vec3f computeColor(const Scene &scene, int depth, bool is_recursive)
    {
        Vec3f color = scene.background_color;
        
        IntersectionResult normal_intersection = (this->intersectObjects(scene,false,-scene.shadow_ray_epsilon/100));
        
        if(normal_intersection.intersected)
        {
            
            color = (normal_intersection.m.ambient).multVecs(scene.ambient_light);
            Vec3f w_o = (d * (-1)).normalize();
            Vec3f x = o + d*normal_intersection.minT;
            Vec3f n = normal_intersection.n;
            Material m = normal_intersection.m;
            Vec3f L_d;
            Vec3f L_s;
            // calculate color
            for(int i = 0; i < scene.point_lights.size(); i++)
            {
                PointLight light = scene.point_lights[i];
                Vec3f w_i = (light.position - x).normalize();
                float length_r = sqrt((light.position - x) * (light.position - x));
               
                Ray shadow_ray(x + n * scene.shadow_ray_epsilon , w_i);
                IntersectionResult shadow_intersection = (shadow_ray.intersectObjects(scene,true,0));
                
                if(shadow_intersection.intersected)
                {
                    if((shadow_intersection.minT > scene.shadow_ray_epsilon) && (shadow_intersection.minT < length_r))
                        continue;
                }
                
                //diffuse
                float cos_alfa;
                float inv_r2 = 1.0/(length_r*length_r); 
                cos_alfa = std::max(w_i*n, 0.0);
                L_d = ((m.diffuse).multVecs(light.intensity) * (cos_alfa * inv_r2)); 
                
                //specular
                Vec3f h = (w_i + w_o).normalize();
                if((h * w_i) > 0)
                {
                    cos_alfa = std::max(n*h,0.0);
                    L_s = (m.specular).multVecs(light.intensity) * ((pow(cos_alfa, m.phong_exponent) * inv_r2));
                }
                else 
                    L_s.x = L_s.y = L_s.z = 0;

                color = color + L_d + L_s;
                
                
            }
            if(m.is_mirror && depth > 0)
                {
                    Vec3f L_m;
                    Vec3f w_r = (w_o * (-1)) + n*(n*w_o) * 2;
                    Ray reflectance(x + n* scene.shadow_ray_epsilon, w_r);
                    L_m = reflectance.computeColor(scene, depth - 1,true).multVecs(m.mirror);
                    color = color + L_m;;
                }
        }
        else if(is_recursive)
            color = Vec3f();
        
        return color;
    }

    IntersectionResult intersectObjects(const Scene& scene, bool shadow, double epsilon)
    {
        Vec3f baryo_param; // triangle
        IntersectionResult res;
        
        int minJ = MAX_I;

        double minT1 = MAX_T;
        double minT2 = MAX_T;
        double minT3 = MAX_T;
       
        for(int i = 0; i < scene.spheres.size(); i++)
        {
            double t = this->intersectSphere(
            scene.vertex_data[scene.spheres[i].center_vertex_id - 1], 
            scene.spheres[i]
            );
            if(t < minT1 && t >= 0)
            {
                res.minI = i;
                minT1 = t;
                res.obj_type = SPHERE;
            }
            if(shadow &&( minT1 != MAX_T)) 
            {
                res.minT = minT1;
                res.intersected = true;
                return res;
            }
        }

        //render meshes
        for(int i = 0; i < scene.meshes.size(); i++)
        {
            Mesh mesh = scene.meshes[i];
            for(int j = 0; j < mesh.faces.size(); j++)
            {
                if(!shadow && (mesh.faces[j].normal * d.normalize()) >0)
                    continue;
                baryo_param = this->intersectTriangle(mesh.faces[j], scene.vertex_data);
                if(baryo_param.y + baryo_param.z - 1 <=  epsilon
                && baryo_param.y >= epsilon
                && baryo_param.z >= epsilon)
                {
                    if(baryo_param.x < minT2 && baryo_param.x >= epsilon)
                    {
                        minT2 = baryo_param.x;
                        if(minT2 < minT1) 
                        {
                            res.minI= i;
                            res.obj_type = MESH_TRIANGLE;
                            minJ = j;

                        }
                        
                    }
                   if(shadow && (minT2 != MAX_T)) 
                    {
                        res.minT = minT2;
                        res.intersected = true;
                        return res;
                    }
                    
                }
            }
            
        }
        //render triangles
        for(int i = 0; i < scene.triangles.size(); i++)
        {
            Face tri = scene.triangles[i].indices;
            if(!shadow && (tri.normal * d.normalize()) >0)
                    continue;
            baryo_param = this->intersectTriangle(tri, scene.vertex_data);
            
            if(baryo_param.y + baryo_param.z  - 1<= epsilon  
                && baryo_param.y >= epsilon 
                && baryo_param.z >= epsilon)
            {
                
                if(baryo_param.x < minT3 && baryo_param.x >= epsilon)
                {
                    minT3 = baryo_param.x;
                    if(minT3 < minT1 || minT3 < minT2) 
                    {
                        res.minI= i;
                        res.obj_type = TRIANGLE;
                    }
                }
                if(shadow && (minT3 != MAX_T)) 
                {
                    res.minT = minT3;
                    res.intersected = true;
                    return res;
                }
            }
                
            
        }
        if((minT3 == MAX_T) && (minT2 == MAX_T) && (minT1 == MAX_T))
            {
                res.intersected = false;
                return res;
            }
            else
            {
                res.intersected = true;

                res.minT = std::min(minT1, minT2);
       		    res.minT = std::min(res.minT, minT3);


                res.minJ = minJ;
                if(res.obj_type == SPHERE)
                {
                    
                    Sphere current_sphere = scene.spheres[res.minI]; 
                    res.m = scene.materials[current_sphere.material_id - 1];
                    Vec3f center = scene.vertex_data[current_sphere.center_vertex_id - 1];
                    res.n = ((o + d*res.minT) - center)* (1.0/current_sphere.radius);
                }
                else if (res.obj_type == TRIANGLE)
                {
                    res.m = scene.materials[scene.triangles[res.minI].material_id - 1];
                    res.n = scene.triangles[res.minI].indices.normal;
                }
                else
                {
                   res.m = scene.materials[scene.meshes[res.minI].material_id - 1];
                   res.n = scene.meshes[res.minI].faces[res.minJ].normal;
                } 

            }
        return res;
    }
    
    float intersectSphere(const Vec3f &center, const Sphere &sphere)
    {
        // parser::Vec3f center ;
        float A, B, C;
        float delta;
        float t;
        
        
        A = (d * d);
        B = (d * (o - center)) * 2;
        C = (o - center) * (o - center) - std::pow(sphere.radius,2);
        delta = B*B - 4*A*C;
        if(delta < 0) return MAX_T;
        // else if (delta < 0.0001)
        // {
            // t = (-B)/(2*A);
        // }
        else 
        {
            //might be irrelevant, since t2 < t1 always ?
            delta = std::sqrt(delta);
            // float t1 = (-B + delta)/(2*A);
            // float t2 = (-B - delta)/(2*A);
            t = (-B - delta)/(2*A);
        }
        return t;
    }
    
    Vec3f intersectTriangle(Face tri, const std::vector<Vec3f> &vertexes)
    {
        Vec3f result;
        
        double t, beta, gamma;
        double 
            ax = vertexes[tri.v0_id - 1].x, 
            bx = vertexes[tri.v1_id - 1].x, 
            cx = vertexes[tri.v2_id - 1].x;
        double 
            ay = vertexes[tri.v0_id - 1].y, 
            by = vertexes[tri.v1_id - 1].y, 
            cy = vertexes[tri.v2_id - 1].y;
        double 
            az = vertexes[tri.v0_id - 1].z, 
            bz = vertexes[tri.v1_id - 1].z, 
            cz = vertexes[tri.v2_id - 1].z;
        double matrix[3][4] = {  
                                {ax-bx, ax-cx, d.x, ax - o.x},
                                {ay-by, ay-cy, d.y, ay - o.y},
                                {az-bz, az-cz, d.z, az - o.z}, 
                            };
        cramer(matrix, t, beta, gamma);

        result.x = t;
        result.y = beta;
        result.z = gamma;

         
        return result;
    }
};

Ray generateRay(int i, int j, const Camera & cam)
{
    //x    y     z       w
    //left right bottom  top
   
    float su, sv;
    Vec3f m, q, s, u, v, w, e;
    Camera my_cam = cam;

    su = (i + 0.5)* (my_cam.near_plane.y - my_cam.near_plane.x)/my_cam.image_width;
    sv = (j + 0.5)* (my_cam.near_plane.w - my_cam.near_plane.z)/my_cam.image_height;

    v = my_cam.up;
    w = my_cam.gaze * (-1);
    u = v & w;
    e = my_cam.position;
    
    
    m = e + (my_cam.gaze * my_cam.near_distance);
    q = m + ((u * my_cam.near_plane.x) + (v * my_cam.near_plane.w));
    s = q + ((u * su)+(v * -sv));
    
    Ray result(e,(s + e * (-1)));
    return result;
}

void precomputeNormals(Scene &scene)
{
    for(int i = 0; i < scene.meshes.size(); i++)
    {
        Mesh &mesh = scene.meshes[i];
        for(int j = 0; j < mesh.faces.size(); j++)
        {
            Face &tri = mesh.faces[j];

            tri.normal = 
            ((scene.vertex_data[tri.v1_id - 1] - scene.vertex_data[tri.v0_id - 1]) & 
            (scene.vertex_data[tri.v2_id - 1] - scene.vertex_data[tri.v0_id - 1])).normalize();
         
        }

    }
    for(int i = 0; i < scene.triangles.size(); i++)
    {
        Face &tri = scene.triangles[i].indices;
        tri.normal = 
            ((scene.vertex_data[tri.v1_id - 1] - scene.vertex_data[tri.v0_id - 1]) & 
            (scene.vertex_data[tri.v2_id - 1] - scene.vertex_data[tri.v0_id - 1])).normalize();
        

    }
}


void render(Scene scene, unsigned char* image, int height, int width, int y_start, int index)
{
    int i = (y_start) * width * 3;
    for (int y = y_start; y < y_start + height/8; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                Ray ray = generateRay(x,y, scene.cameras[index]); 
                // Ray ray = generateRay(509,397, scene.cameras[cam_index]); 

                Vec3f pixel = ray.d + ray.o;
                Vec3f rayColor = ray.computeColor(scene,scene.max_recursion_depth, false);
                rayColor = rayColor.clamp(0.0f, 255.0f);
                image[i++] = (int)(rayColor.x);
                image[i++] = (int)(rayColor.y);
                image[i++] = (int)(rayColor.z);

            }
            
        }
}

int main(int argc, char* argv[])
{
    
    Scene scene;


    scene.loadFromXml(argv[1]);
   
    std::thread threads[8];

    unsigned char * image;
    
    precomputeNormals(scene);

    
    
   
    for(int cam_index = 0; cam_index < scene.cameras.size(); cam_index++)
    {
        
        int width = scene.cameras[cam_index].image_width; 
        int height = scene.cameras[cam_index].image_height;
        image = new unsigned char [width * height * 3];
        std::string file_name = scene.cameras[cam_index].image_name;
        for(int j = 0; j < 8; j++)
        {
            threads[j] = std::thread(render, scene, image, height, width, j*height/8, cam_index);
        }
        for(int j = 0; j < 8; j++)
        {
            threads[j].join();
        }
        
        write_ppm((file_name).c_str(), image, width, height);
        // {std::cout<<std::setw(3)<<"%"<< (int)((float)/(width * height * 3) * 100) << '\r'<<std::flush;}

        delete [] image;
        
    }


}
