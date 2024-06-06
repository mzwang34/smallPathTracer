#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    //Scene scene(784, 784);
    Scene scene(512, 512);
    //Scene scene(100, 100);
    
    // ----------------- diffuse ----------------- 
    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f + 0.058f, 0.747f + 0.258f, 0.747f) + 15.6f * Vector3f(0.740f + 0.287f, 0.740f + 0.160f, 0.740f) + 18.4f * Vector3f(0.737f + 0.642f, 0.737f + 0.159f, 0.737f)));
    light->Kd = Vector3f(0.65f);

    // ----------------- microfacet ----------------- 
    //Material* red = new Material(MICROFACET, Vector3f(0.0f));
    //red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    //red->roughness = .8f;
    //red->metallic = 0.3f;
    //Material* green = new Material(MICROFACET, Vector3f(0.0f));
    //green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    //green->roughness = .5f;
    //green->metallic = 0.5f;
    //Material* white = new Material(MICROFACET, Vector3f(0.0f));
    //white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    //white->roughness = .2f;
    //white->metallic = 0.2f;
    //Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    //light->Kd = Vector3f(0.65f);
    //Material* microfacet = new Material(MICROFACET, Vector3f(0.0f));
    //microfacet->Kd = Vector3f(0.913f, 0.921f, 0.925f);
    //microfacet->roughness = 0.05f;
    //microfacet->metallic = 1.f;
    //Material* microfacet2 = new Material(MICROFACET, Vector3f(0.0f));
    //microfacet2->Kd = Vector3f(1.000f, 0.766f, 0.336f);
    //microfacet2->roughness = 0.2f; // 0.2f
    //microfacet2->metallic = 0.8f;
    //Material* microfacet3 = new Material(MICROFACET, Vector3f(0.0f));
    //microfacet3->Kd = Vector3f(1.f, 0.86f, 0.57f);
    //microfacet3->roughness = 0.15f;
    //microfacet3->metallic = 0.8f;

    MeshTriangle floor("../models/cornellbox/floor.obj", white);
    //MeshTriangle shortbox("../models/cornellbox/shortbox.obj", white);
    //MeshTriangle tallbox("../models/cornellbox/tallbox.obj", white);
    MeshTriangle left("../models/cornellbox/left.obj", red);
    MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/light.obj", light);
    /*MeshTriangle bunny("../models/bunny/bunny.obj", microfacet2, Vector3f(200, -60, 150), Vector3f(1500, 1500, 1500),
                        Vector3f(-1, 0, 0), Vector3f(0, 1, 0), Vector3f(0, 0, -1));*/
    //Sphere sph1(Vector3f(120, 100, 300), 100, white);
    Sphere sph1(Vector3f(278, 273, 150), 80, white);
    Sphere sph2(Vector3f(400, 100, 400), 50, red);
    Sphere sph3(Vector3f(100, 150, 370), 30, green);

    // ----------------- mirror ----------------- 
    /*Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f + 0.058f, 0.747f + 0.258f, 0.747f) + 15.6f * Vector3f(0.740f + 0.287f, 0.740f + 0.160f, 0.740f) + 18.4f * Vector3f(0.737f + 0.642f, 0.737f + 0.159f, 0.737f)));
    light->Kd = Vector3f(0.65f);
    Material* mirror = new Material(MIRROR, Vector3f(0.0f));
    mirror->ior = 10.f;*/

    //MeshTriangle floor("../models/cornellbox/floor.obj", white);
    ////MeshTriangle tallbox("../models/cornellbox/tallbox.obj", mirror);
    //MeshTriangle left("../models/cornellbox/left.obj", red);
    //MeshTriangle right("../models/cornellbox/right.obj", green);
    //MeshTriangle light_("../models/cornellbox/light.obj", light);
    //MeshTriangle bunny("../models/bunny/bunny.obj", mirror, Vector3f(200, -60, 150), Vector3f(1500, 1500, 1500),
    //                    Vector3f(-1, 0, 0), Vector3f(0, 1, 0), Vector3f(0, 0, -1));
    //Sphere sph1(Vector3f(120, 100, 300), 100, mirror);

    scene.Add(&floor);
    /*scene.Add(&shortbox);
    scene.Add(&tallbox);*/
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);
    //scene.Add(&bunny);
    scene.Add(&sph1);
    scene.Add(&sph2);
    scene.Add(&sph3);

    // ----------------- mis ----------------- 
    /*Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* plate1 = new Material(MICROFACET, Vector3f(0.0f));
    plate1->Kd = Vector3f((0.56, 0.57, 0.58));
    plate1->roughness = 0.001f;
    plate1->F0 = Vector3f((0.56, 0.57, 0.58));
    plate1->metallic = 1.0f;
    Material* plate2 = new Material(MICROFACET, Vector3f(0.0f));
    plate2->Kd = Vector3f((0.56, 0.57, 0.58));
    plate2->roughness = 0.05f;
    plate2->F0 = Vector3f((0.56, 0.57, 0.58));
    plate2->metallic = 1.0f;
    Material* plate3 = new Material(MICROFACET, Vector3f(0.0f));
    plate3->Kd = Vector3f((0.56, 0.57, 0.58));
    plate3->roughness = 0.1f;
    plate3->F0 = Vector3f((0.56, 0.57, 0.58));
    plate3->metallic = 1.0f;
    Material* plate4 = new Material(MICROFACET, Vector3f(0.0f));
    plate4->Kd = Vector3f((0.56, 0.57, 0.58));
    plate4->roughness = 0.2f;
    plate4->F0 = Vector3f((0.56, 0.57, 0.58));
    plate4->metallic = 1.0f;

    Material* light_mis1 = new Material(DIFFUSE, Vector3f(2.0));
    light_mis1->Kd = Vector3f(1.0f);
    Material* light_mis2 = new Material(DIFFUSE, Vector3f(0.0,0.0,10.0));
    light_mis2->Kd = Vector3f(1.0f);
    Material* light_mis3 = new Material(DIFFUSE, Vector3f(0.0,50.0,0.0));
    light_mis3->Kd = Vector3f(1.0f);
    Material* light_mis4 = new Material(DIFFUSE, Vector3f(400.0,0.0,0.0));
    light_mis4->Kd = Vector3f(1.0f);

    MeshTriangle floor("../models/meshes/floor.obj", white);
    MeshTriangle plate1_obj("../models/meshes/plate1.obj", plate1);
    MeshTriangle plate2_obj("../models/meshes/plate2.obj", plate2);
    MeshTriangle plate3_obj("../models/meshes/plate3.obj", plate3);
    MeshTriangle plate4_obj("../models/meshes/plate4.obj", plate4);

    Sphere sphere1(Vector3f(3.5, 0, 0),     1.0, light_mis1);
    Sphere sphere2(Vector3f(1, 0, 0),       0.6, light_mis2);
    Sphere sphere3(Vector3f(-1.7, 0, 0),    0.2, light_mis3);
    Sphere sphere4(Vector3f(-4, 0, 0),      0.05, light_mis4);

    scene.Add(&floor);
    scene.Add(&plate1_obj);
    scene.Add(&plate2_obj);
    scene.Add(&plate3_obj);
    scene.Add(&plate4_obj);
    scene.Add(&sphere1);
    scene.Add(&sphere2);
    scene.Add(&sphere3);
    scene.Add(&sphere4);*/

    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}