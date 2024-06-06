//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include <thread>
#include <mutex>

std::mutex mtx;


inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 1e-8; //0.00001

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    //Vector3f eye_pos(0, 0, 11);
    //Vector3f eye_pos(0);
    int m = 0;

    // change the spp value to change sample ammount
    int spp = 256; // 16
    std::cout << "SPP: " << spp << "\n";
    // thread @wmz
    int num_threads = 8;
    std::thread th[num_threads];
    int thread_height = scene.height / num_threads;

    //for (uint32_t j = 0; j < scene.height; ++j) {
    //    for (uint32_t i = 0; i < scene.width; ++i) {
    //        // generate primary ray direction
    //        float x = (2 * (i + 0.5) / (float)scene.width - 1) *
    //                  imageAspectRatio * scale;
    //        float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

    //        Vector3f dir = normalize(Vector3f(-x, y, 1));
    //        for (int k = 0; k < spp; k++){
    //            framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;  
    //        }
    //        m++;
    //    }
    //    UpdateProgress(j / (float)scene.height);
    //}
    float aperture = 0.01f;
    float focus = 8.f;

    int process = 0;
    auto castRayMultithread = [&](uint32_t ymin, uint32_t ymax) {
        for (uint32_t j = ymin; j < ymax; ++j) {
            int m = j * scene.width;
            for (uint32_t i = 0; i < scene.width; ++i) {
                // generate primary ray direction
                /*float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                    imageAspectRatio * scale;
                float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;
                Vector3f dir = normalize(Vector3f(-x, y, 1));*/

                Vector3f offset = aperture * 0.5f * random_in_unit_disk();
                eye_pos += offset;

                for (int k = 0; k < spp; k++){
                    // MSAA
                    float x = (2 * (i + get_random_float()) / (float)scene.width - 1) *
                        imageAspectRatio * scale;
                    float y = (1 - 2 * (j + get_random_float()) / (float)scene.height) * scale;
                    Vector3f dir = normalize(Vector3f(-x, y, 1));
                    framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;  
                }
                m++;
            }
            mtx.lock();
            process++;
            UpdateProgress(1.0 * process / (float)scene.height);
            mtx.unlock();
        }
    };

    // light field camera
    //auto castRayMultithread = [&](uint32_t ymin, uint32_t ymax) {
    //    for (uint32_t j = ymin; j < ymax; ++j) {
    //        int m = j * scene.width;
    //        for (uint32_t i = 0; i < scene.width; ++i) {
    //            Vector3f offset = aperture * 0.5f * random_in_unit_disk();
    //            eye_pos += offset;

    //            for (int k = 0; k < spp; k++) {
    //                // MSAA
    //                float x = (2 * (i + get_random_float()) / (float)scene.width - 1) *
    //                    imageAspectRatio * scale;
    //                float y = (1 - 2 * (j + get_random_float()) / (float)scene.height) * scale;
    //                //Vector3f dir = normalize(Vector3f(-x, y, 1));
    //                Vector3f dir(-x, y, 1);
    //                dir = dir * focus - offset;
    //                dir = normalize(dir);
    //                framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
    //            }
    //            eye_pos = Vector3f(278, 273, -800);
    //            m++;
    //        }
    //        mtx.lock();
    //        process++;
    //        UpdateProgress(1.0 * process / (float)scene.height);
    //        mtx.unlock();
    //    }
    //};

    //auto castRayMultithread = [&](uint32_t ymin, uint32_t ymax) {
    //    for (uint32_t j = ymin; j < ymax; ++j) {  
    //        int m = j * scene.width;
    //        for (uint32_t i = 0; i < scene.width; ++i) {
    //            for (int k = 0; k < spp; k++) {
    //                float x = (2 * (i + get_random_float()) / (float)scene.width - 1) * imageAspectRatio * scale;
    //                float y = (1 - 2 * (j + get_random_float()) / (float)scene.height) * scale;
    //                Vector3f dir = normalize(Vector3f(x, y, -1));
    //                dir.rotateAxis(Vector3f(1, 0, 0), 10);
    //                framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
    //            }
    //            m++;
    //        }
    //        mtx.lock();
    //        process++;
    //        UpdateProgress(1.0 * process / (float)scene.height);
    //        mtx.unlock();
    //    }
    //};

    for (int i = 0; i < num_threads; ++i) 
        th[i] = std::thread(castRayMultithread, i * thread_height, (i + 1) * thread_height);

    for (int i = 0; i < num_threads; ++i)
        th[i].join();

    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("lightfieldcamera_ori.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        //framebuffer[i] = ACESToneMapping(framebuffer[i], 2.f);
        /*color[0] = (unsigned char)(255 * std::pow(ACESToneMapping(clamp(0, 1, framebuffer[i].x), 1.f), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(ACESToneMapping(clamp(0, 1, framebuffer[i].y), 1.f), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(ACESToneMapping(clamp(0, 1, framebuffer[i].z), 1.f), 0.6f));*/
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
