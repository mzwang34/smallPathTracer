//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

void Scene::sampleLightSphere(Intersection& pos, float& pdf) const
{
    float emit_areaMultiEmit_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k)
    {
        if (objects[k]->hasEmit())
        {
            //面积和光的强度做比重
            emit_areaMultiEmit_sum += objects[k]->getEmitNorm() * objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_areaMultiEmit_sum;
    emit_areaMultiEmit_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k)
    {
        if (objects[k]->hasEmit())
        {
            emit_areaMultiEmit_sum += objects[k]->getEmitNorm() * objects[k]->getArea();
            //随机找到一个光源面，再在这个光源面中找到一个点
            if (p <= emit_areaMultiEmit_sum)
            {
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Intersection inter = intersect(ray);
    // if no intersection
    if (!inter.happened)
        return Vector3f();
    
    // if ray hits the light
    if (inter.m->hasEmission() && depth == 0) 
        return inter.m->getEmission();

    float pdf_light = 1.0; 
    Intersection light_pos;
    sampleLight(light_pos, pdf_light);

    // object    
    auto p = inter.coords;
    auto N = inter.normal;
    auto wo = -ray.direction; // negative
    // light
    auto x = light_pos.coords;
    auto NN = light_pos.normal;
    auto ws = (x - p).normalized();
    auto emit = light_pos.emit;

    float EPSILON = 0.005f;//1e-5
    float ws_dis = (x - p).norm();

    if (inter.m->getType() == MIRROR) {
        Vector3f L_indir = 0.0;
        if (get_random_float() < RussianRoulette) {
            auto wi = inter.m->sample(wo, N).normalized();
            Ray r(p, wi);
            Intersection wi_inter = intersect(r);
            if (wi_inter.happened) {
                float pdf = inter.m->pdf(wo, wi, N);
                if (pdf > EPSILON)
                     L_indir = castRay(r, depth + 1) * inter.m->eval(wo, wi, N) * dotProduct(wi, N) / inter.m->pdf(wo, wi, N) / RussianRoulette;
            }
        }
        return L_indir;
    }

    // shoot a ray from p to x
    Ray p2x_ray(p, ws);
    Intersection ws_ray_inter = intersect(p2x_ray);
    // if the ray is not blocked in the middle
    Vector3f L_dir;
    if (ws_ray_inter.happened && (ws_ray_inter.distance - ws_dis) > -EPSILON)
        L_dir = emit * inter.m->eval(wo, ws, N) * dotProduct(ws, N) * dotProduct(-ws, NN) / dotProduct(x - p, x - p) / pdf_light; //std::pow(ws_dis, 2)

    Vector3f L_indir = 0.0;
    //if (get_random_float() > RussianRoulette)
    //    return L_dir;
    //    //return Vector3f::Min(Vector3f::Max(L_dir, Vector3f(0.0f)), Vector3f(1.0f));;

    if (get_random_float() < RussianRoulette) {
        //Vector3f wi = (inter.m->sample(wo, N)).normalized();
        Vector3f wi = inter.m->sample(wo, N);
        float pdf = inter.m->pdf(wo, wi, N);
        auto fr = inter.m->eval(wo, wi, N);
        wi = wi.normalized();
        Ray r(p, wi);
        Intersection wi_inter = intersect(r);
        if (wi_inter.happened && !wi_inter.m->hasEmission()) {
            // white denoise
            if (pdf > EPSILON)
                //L_indir = castRay(r, depth + 1) * inter.m->eval(wo, wi, N) * dotProduct(wi, N) / pdf / RussianRoulette;
                L_indir = castRay(r, depth + 1) * fr * dotProduct(wi, N) / pdf / RussianRoulette;
        }
    }

    //return L_dir + L_indir;
    return Vector3f::Min(Vector3f::Max(L_dir + L_indir, Vector3f(0.0f)), Vector3f(1.0f));   // reduce noise but lose energy
}


// sample BRDF
//Vector3f Scene::castRay(const Ray& ray, int depth) const
//{
//    Intersection inter = intersect(ray);
//    // if no intersection
//    if (!inter.happened)
//        return Vector3f();
//
//    Vector3f L;
//    if (inter.m->hasEmission())
//        return inter.m->getEmission();
//    
//    auto p = inter.coords;
//    auto N = inter.normal;
//    auto wo = -ray.direction;
//    auto wi = (inter.m->sample(wo, N)).normalized();
//    auto fr = inter.m->eval(wo, wi, N);
//    Ray r(p, wi);
//    Intersection wi_inter = intersect(r);
//    float pdf = dotProduct(wi, N) / M_PI;   // cos weight pdf
//    if (wi_inter.happened && wi_inter.m->hasEmission()) {
//        L += wi_inter.m->getEmission() * fr * std::max(dotProduct(wi, N), 0.f) / pdf;
//    }
//
//    if (get_random_float() > RussianRoulette) return L;
//
//    if (wi_inter.happened && !wi_inter.m->hasEmission()) {
//        L += castRay(r, depth + 1) * fr * std::max(dotProduct(wi, N), 0.f) / pdf / RussianRoulette;
//    }
//    return L;
//}


// light only sampling for mis test
//Vector3f Scene::castRay(const Ray& ray, int depth) const 
//{
//    Intersection inter = intersect(ray);
//    if (!inter.happened)
//        return Vector3f();
//    
//    if (inter.m->hasEmission() && depth == 0)
//        return inter.m->getEmission();
//    
//    // direct light
//    Vector3f L_dir;
//    auto p = inter.coords;
//    auto N = inter.normal;
//    auto wo = -ray.direction;
//
//    Intersection L_dir_inter;
//    float pdf_light;
//    L_dir_inter.tcoords = inter.coords;
//    sampleLightSphere(L_dir_inter, pdf_light);    // sample sphere
//
//    // light
//    auto x = L_dir_inter.coords;
//    auto ws = (x - p).normalized();
//    Ray p2x_ray(p, ws);
//    Intersection p2x_inter = intersect(p2x_ray);
//    if (p2x_inter.happened && p2x_inter.m->hasEmission()) {
//        if (pdf_light > EPSILON)
//            L_dir = L_dir_inter.emit * inter.m->eval(wo, ws, N) * std::max(dotProduct(ws, N) ,0.f) / pdf_light;
//    }
//
//    if (get_random_float() > RussianRoulette)
//        return L_dir;
//
//    // indirect light
//    Vector3f L_indir;
//    auto wi = inter.m->sample(wo, N);
//    float pdf = inter.m->pdf(wo, wi, N);
//    wi = wi.normalized();
//    Ray L_indir_ray(p, wi);
//    Intersection L_indir_inter = intersect(L_indir_ray);
//    
//    if (L_indir_inter.happened && !L_indir_inter.m->hasEmission()) {
//        if (pdf > EPSILON)
//            L_indir = castRay(L_indir_ray, depth + 1) * inter.m->eval(wo, wi, N) * std::max(dotProduct(wi, N), 0.f) / pdf / RussianRoulette;
//    }
//
//    return L_dir + L_indir;
//}

// BRDF only sampling for mis test
//Vector3f Scene::castRay(const Ray& ray, int depth) const
//{
//    Intersection inter = intersect(ray);
//    // if no intersection
//    if (!inter.happened)
//        return Vector3f();
//
//    
//    if (inter.m->hasEmission() && depth == 0)
//        return inter.m->getEmission();
//
//    Vector3f L_dir(0.f);
//    auto p = inter.coords;
//    auto N = inter.normal;
//    auto wo = -ray.direction;
//
//    auto wi = inter.m->sample(wo, N);
//    float pdf = inter.m->pdf(wo, wi, N);
//    Vector3f fr = inter.m->eval(wo, wi, N);
//    wi = normalize(wi);
//    Ray r(p, wi);
//    Intersection wi_inter = intersect(r);
//    //float pdf = dotProduct(wi, N) / M_PI;   // cos weight pdf
//    if (wi_inter.happened && wi_inter.m->hasEmission()) {
//        L_dir += wi_inter.m->getEmission() * fr * std::max(dotProduct(wi, N), 0.f) / pdf;
//    }
//    if (get_random_float() > RussianRoulette) return L_dir;
//
//    Vector3f L_indir;
//    if (wi_inter.happened && !wi_inter.m->hasEmission()) {
//        if (pdf > EPSILON)
//            L_indir = castRay(r, depth + 1) * inter.m->eval(wo, wi, N) * dotProduct(wi, N) / pdf / RussianRoulette;
//    }
//
//    return L_dir + L_indir;
//}

// mis
//Vector3f Scene::castRay(const Ray& ray, int depth) const
//{
//    Intersection inter = intersect(ray);
//    // if no intersection
//    if (!inter.happened)
//        return Vector3f();
//
//    if (inter.m->hasEmission() && depth == 0)
//        return inter.m->getEmission();
//
//    Vector3f L_dir(0.f);
//    auto p = inter.coords;
//    auto N = inter.normal;
//    auto wo = -ray.direction;
//
//    // BRDF
//    auto wi = inter.m->sample(wo, N);
//    float pdf = inter.m->pdf(wo, wi, N);
//    Vector3f fr = inter.m->eval(wo, wi, N);
//    wi = normalize(wi);
//
//    Ray r(p, wi);
//    Intersection wi_inter = intersect(r);
//    //float pdf = dotProduct(wi, N) / M_PI;   // cos weight pdf
//
//    if (wi_inter.happened && wi_inter.m->hasEmission()) {
//        float pdf_brdf = pdf;
//        float pdf_light;
//        Intersection inter_tmp;
//        inter_tmp.tcoords = inter.coords;
//        wi_inter.obj->Sample(inter_tmp, pdf_light);
//        float weight_brdf = PowerHeuristic(1.f, pdf_brdf, 1.f, pdf_light);
//        if (pdf > EPSILON)
//            L_dir += wi_inter.m->getEmission() * fr * std::max(dotProduct(wi, N), 0.f) / pdf * weight_brdf;
//    }
//    
//    // light
//    Intersection L_dir_inter;
//    float pdf_light;
//    L_dir_inter.tcoords = inter.coords;
//    sampleLightSphere(L_dir_inter, pdf_light);    // sample sphere
//    auto x = L_dir_inter.coords;
//    auto ws = (x - p).normalized();
//    //pdf_brdf = inter.m->pdf(wo, ws, N);
//    Ray p2x_ray(p, ws);
//    Intersection p2x_inter = intersect(p2x_ray);
//    /*float weight_light = PowerHeuristic(1.f, pdf_light, 1.f, pdf_brdf);*/
//    if (p2x_inter.happened && p2x_inter.m->hasEmission()) {
//
//        float pdf_brdf = inter.m->pdf(wo, ws, N);
//        float weight_light = PowerHeuristic(1.f, pdf_light, 1.f, pdf_brdf);
//        if (pdf_light > EPSILON)
//            L_dir += L_dir_inter.emit * inter.m->eval(wo, ws, N) * std::max(dotProduct(ws, N) ,0.f) / pdf_light * weight_light;
//    }
//
//    if (get_random_float() > RussianRoulette) return L_dir;
//
//    Vector3f L_indir;
//    if (wi_inter.happened && !wi_inter.m->hasEmission()) {
//        if (pdf > EPSILON)
//            L_indir = castRay(r, depth + 1) * inter.m->eval(wo, wi, N) * dotProduct(wi, N) / pdf / RussianRoulette;
//    }
//
//    return L_dir + L_indir;
//}