// Ray_tracer.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include "rtweekend.h"

#include"camera.h"
//#include "color.h"
//#include "vec3.h"
//#include"ray.h"
//#include "hittable.h"
#include"color.h"
#include"material.h"
#include "hittable_list.h"
#include "sphere.h"

//#include <iostream>

/*
//time when ray hits the closest object
double hit_sphere(const point3& center, double radius, const ray& r) {
    //center is the center of the sphere;
    //auto a = dot(r.direction(), r.direction());
    vec3 oc = r.origin() - center;
    //auto b = 2 * dot(r.direction(), oc);
    //auto c = dot(oc, oc) - radius * radius;
    //auto discriminant = b * b - 4 * a * c;
    //return(discriminant >= 0);
    auto a = r.direction().length_squared();
    auto half_b = dot(r.direction(), oc);
    auto c = oc.length_squared() - radius * radius;
    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-half_b - sqrt(discriminant)) / a;
    }
}
*/
/*
color ray_color(const ray& r,const hittable& world) {
    
    //t is the minimum time hitting certain object
    auto t = hit_sphere(point3(0, 0, -1), 0.5, r);
    if (t > 0.0) {
        //the result of at function is the point the ray hit the sphere at time t
        vec3 N = unit_vector(r.at(t) - vec3(0, 0, -1));
        return 0.5 * color(N.x() + 1, N.y() + 1, N.z() + 1);
    }
    
    if (hit_sphere(point3(0, 0, -1), 0.5, r))
        return color(1, 0, 0);
    
    hit_record rec;
    if (world.hit(r, interval(0,infinity), rec)) {
        return 0.5 * (rec.normal + color(1, 1, 1));
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto a = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
}
*/

int main()
{
    //When start a renderer, we need a way to see the image. 
    //PPM file is a nice format to store data of the image.
    //Image:determine the width and height of the image
    //int image_width = 256;
    //int image_height = 256;

    //Image

    //auto aspect_ratio = 16.0 / 9.0;
    //int image_width = 400;

    //calculate the image height, plus ensure that it's at least 1
    //int image_height = static_cast<int>(image_width / aspect_ratio);
    //image_height = (image_height < 1) ? 1 : image_height;

    //world

    hittable_list world;

    //auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    //auto material_center = make_shared<lambertian>(color(0.7, 0.3, 0.3));
    //auto material_left = make_shared<metal>(color(0.8, 0.8, 0.8),0.3);
    //auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    //auto material_left = make_shared<dielectric>(1.5);
    //auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2),0.0);

    //world.add(make_shared<sphere>(point3(0.0, -100.5, -1.0), 100.0, material_ground));
    //world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));
    //world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    //world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));


    //world.add(make_shared<sphere>(point3(0, 0, -1), 0.5));
    //world.add(make_shared<sphere>(point3(0, -100.5, -1), 100));

    //auto R = cos(pi / 4);
    
    //auto material_left = make_shared<lambertian>(color(0, 0, 1));
    //auto material_right = make_shared<lambertian>(color(1, 0, 0));

    //world.add(make_shared<sphere>(point3(-R, 0, -1), R, material_left));
    //world.add(make_shared<sphere>(point3(R, 0, -1), R, material_right));

    /*
    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    auto material_left = make_shared<dielectric>(1.5);
    auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);

    world.add(make_shared<sphere>(point3(0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));
    world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), -0.4, material_left));
    world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));
    */
    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 1200;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;

    cam.vfov = 20;
    cam.lookfrom = point3(13, 2, 3);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0.6;
    cam.focus_dist = 10.0;


    cam.render(world);

    /*
    //Camera

    auto focal_length = 1.0;
    auto camera_center = point3(0, 0, 0);
    auto viewport_height = 2.0;
    auto viewport_width = viewport_height * (static_cast<double>(image_width) / image_height);
    //calculate the vectors across the horizontal and down the vertical viewport edges
    //Attention, this is to calculate the vector of both horizontal and vertical
    auto viewport_u = vec3(viewport_width, 0, 0);
    auto viewport_v = vec3(0, -viewport_height, 0);

    //calculate the horizontal and vertical delta vectors from pixel to pixel
    auto pixel_delta_u = viewport_u / image_width;
    auto pixel_delta_v = viewport_v / image_height;

    //calculate the location of the upper left pixel
    auto viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u / 2 - viewport_v / 2;
    auto pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

    //Render
    //The pixels are written out in rows.
    //Every row of pixels is written out left to right.
    //each red/green/blue components are represented from 0.0 to 1.0

    //format of PPM file
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = 0; j < image_height; ++j) {
        // \r is enter
        std::clog << "\rScanlines remaining:" << (image_height - j) << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            
            auto r = double(i) / (image_width - 1);
            auto g = double(j) / (image_height - 1);
            auto b = 0;

            // Then trans to [0,255)
            int ir = static_cast<int>(255.999 * r);
            int ig = static_cast<int>(255.999 * g);
            int ib = static_cast<int>(255.999 * b);

            std::cout << ir << ' ' << ig << ' ' << ib << '\n';
           
            //auto pixel_color = color(double(i) / (image_width - 1), double(j) / (image_height - 1), 0);
            auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
            auto ray_direction = pixel_center - camera_center;
            ray r(camera_center, ray_direction);

            color pixel_color = ray_color(r,world);
            write_color(std::cout, pixel_color);
        }
    }
    std::clog << "\rDone.                 \n";
    */


}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
