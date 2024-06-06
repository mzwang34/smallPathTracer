//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType { DIFFUSE, MICROFACET, MIRROR};

class Pdf {
public:
    Pdf() :pdfValue(0), halfVec(Vector3f(0)), isReflect(true) {}
    float pdfValue;
    Vector3f halfVec;
    bool isReflect;
};

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f fresnelSchlick(const Vector3f &F0, const Vector3f &V, const Vector3f &H)
    {
        // TODO: To calculate Schlick F here
        float VdotH = std::max(dotProduct(V, H), 0.0f);

        return F0 + (Vector3f(1.0) - F0) * std::pow(1.0 - VdotH, 5.0);
    }

    float DistributionGGX(const Vector3f &N, const Vector3f &H, const float &roughness)
    {
        // TODO: To calculate GGX NDF here
        float a = roughness; // mis test  
        //float a = roughness * roughness; // black around?
        float a2 = a * a;
        float NdotH = std::max(dotProduct(N, H), 0.0f);
        float NdotH2 = NdotH * NdotH;

        float nom = a2;
        float denom = (NdotH2 * (a2 - 1.0) + 1.0);
        denom = M_PI * denom * denom;

        //return nom / std::max(denom, 0.01f);
        return nom / denom;
    }

    float GeometrySchlickGGX(const float &NdotV, const float &roughness)
    {
        // TODO: To calculate Smith G1 here
        float a = roughness + 1.0;
        //float k = (a * a) / 2.f;
        float k = (a * a) / 8.0;

        float nom = NdotV;
        float denom = NdotV * (1.f - k) + k;

        return nom / denom;
    }

    float GeometrySmith(const Vector3f &N, const Vector3f &V, const Vector3f &L, const float &roughness)
    {
        // TODO: To calculate Smith G here
        float NoV = std::max(dotProduct(N, V), 0.0f);
        float NoL = std::max(dotProduct(N, L), 0.0f);
        float ggx2 = GeometrySchlickGGX(NoV, roughness);
        float ggx1 = GeometrySchlickGGX(NoL, roughness);

        return ggx1 * ggx2;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

    Vector3f mix(const Vector3f& a, const Vector3f& b, const float &x) {
        return (1 - x) * a + x * b;
    }

    /*float NDF(const float& Roughness, const Vector3f& h, const Vector3f& N) {
        float cosTheta = dotProduct(h, N);
        float divisor = 1.0 + cosTheta * cosTheta * (Roughness * Roughness - 1);
        float divisor2 = M_PI * divisor * divisor;
        if (divisor2 < EPSILON)
            return (Roughness * Roughness) / EPSILON;
        else
            return (Roughness * Roughness) / divisor2;
    }*/

    Vector3f Fresnel_Schlick(const Vector3f& h, const Vector3f& wo)
    {
        if (ior != 0) {
            F0 = (1 - ior) / (1 + ior);
            F0 = F0 * F0;
        }
        else {
            F0 = mix(F0, Kd, metallic);
        }
        float subItem = 1.0 - dotProduct(h, wo);
        subItem = subItem * subItem * subItem * subItem * subItem;
        return F0 + (Vector3f(1.0f) - F0) * subItem;
    }

    float NDF(const float& Roughness, const Vector3f& h, const Vector3f& N) {
        double cosTheta = dotProduct(h, N);
        double divisor = 1.0 + cosTheta * cosTheta * (Roughness * Roughness - 1);
        double divisor2 = M_PI * divisor * divisor;
        if (divisor2 < EPSILON)
            return (Roughness * Roughness) / EPSILON;
        else
            return (Roughness * Roughness) / divisor2;
        //float a2 = Roughness * Roughness;
        //a2 = a2 * a2;
        //float costheta = dotProduct(h, N);  // compare with 0?
        //float exp = (a2 - 1) * costheta * costheta + 1;
        //float D = a2 / (M_PI * exp * exp);
        //return D;
    }

    float G1G2(const float& Roughness, const Vector3f& wi, const Vector3f& wo, const Vector3f& N) {
        float A_wi, A_wo;
        A_wi = (-1 + sqrt(1 + Roughness * Roughness * pow(tan(acos(dotProduct(wi, N))), 2))) / 2;
        A_wo = (-1 + sqrt(1 + Roughness * Roughness * pow(tan(acos(dotProduct(wo, N))), 2))) / 2;
        float divisor = (1 + A_wi + A_wo);
        if (divisor < EPSILON)
            return 1.0f / EPSILON;
        else
            return 1.0f / divisor;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;
    float roughness;
    float metallic;
    Vector3f F0;

    //Vector3f halfVec;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    //inline Vector3f sample(const Vector3f &wi, const Vector3f &N, Pdf& pdf);
    inline Vector3f sample(const Vector3f& wi, const Vector3f& N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    //inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N, const Pdf& pdf);
    inline Vector3f eval(const Vector3f& wi, const Vector3f& wo, const Vector3f& N);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


//Vector3f Material::sample(const Vector3f &wi, const Vector3f &N, Pdf& pdf){
Vector3f Material::sample(const Vector3f& wi, const Vector3f& N) {
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            
            break;
        }
        case MICROFACET:
        {
        //    // uniform sample on the hemisphere
        //    /*float x_1 = get_random_float(), x_2 = get_random_float();
        //    float z = std::fabs(1.0f - 2.0f * x_1);
        //    float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
        //    Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
        //    return toWorld(localRay, N);*/

            float x_1 = get_random_float(), x_2 = get_random_float();
            float a2 = roughness * roughness * roughness * roughness;   // !!!
            //float theta = std::atan(std::sqrt((a2 * x_1) / (1 - x_1))); // polar angle
            float theta = std::acos(std::sqrt((1 - x_1) / (x_1 * (a2 - 1) + 1)));
            float phi = 2 * M_PI * x_2; // azimuth angle
            Vector3f localRay(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));

            return reflect(-wi, toWorld(localRay, N));

            break;
        }
        case MIRROR :
        {
            Vector3f localRay = reflect(-wi, N);
            return localRay;

            break;
        }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case MICROFACET:
        {
            /*if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;*/

            float a2 = roughness * roughness;
            Vector3f H = normalize(wi + wo);
            float costheta = dotProduct(H, N);  // compare with 0?
            float exp = (a2 - 1) * costheta * costheta + 1;
            float D = a2 / (M_PI * exp * exp);
            return D * costheta / (4.f * dotProduct(H, wo));

            break;
        }
        case MIRROR :
        {
            if (dotProduct(wo, N) > 0.0f)
                return 1.f;
            else
                return 0.0f;
            break;
        }
    }
}

//Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N, const Pdf& pdf){
Vector3f Material::eval(const Vector3f & wi, const Vector3f & wo, const Vector3f & N) {
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MICROFACET:    // Cook-Torrance  BRDF (reflect only)
        {
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                //Vector3f V = -wi;
                //Vector3f L = wo;
                ////Vector3f H = (V + L).normalized();
                //Vector3f H = pdf.halfVec;

                //float NdotV = std::max(dotProduct(N, V), 0.f);
                //float NdotL = std::max(dotProduct(N, L), 0.f);

                ////Vector3f F0(0.04f);
                //F0 = mix(F0, Kd, metallic);

                ////float roughness = 0.f;
                //float NDF = DistributionGGX(N, H, roughness);
                //float G = GeometrySmith(N, V, L, roughness);
                //Vector3f F = fresnelSchlick(F0, V, H);
                ///*float F;
                //fresnel(V, N, ior, F);*/

                //Vector3f numerator = NDF * G * F;
                ////float numerator = NDF * G * F;
                //float denominator = std::max((4.f * NdotL * NdotV), 0.001f);
                ////float specular = numerator / denominator;
                //Vector3f specular = numerator / denominator;

                ////float ks = F;
                ////float kd = 1.f - ks;
                //Vector3f kd = Vector3f(1.f) - F;

                //Vector3f diffuse = Kd * 1.f / M_PI;
                //return kd * diffuse * (1 - metallic) + specular;
                ////return specular;

                Vector3f V = wi;
                Vector3f L = wo;
                Vector3f H = normalize(V + L);

                float NdotV = std::max(dotProduct(N, V), 0.f);
                float NdotL = std::max(dotProduct(N, L), 0.f);

                F0 = mix(F0, Kd, metallic);

                float D = DistributionGGX(N, H, roughness);
                //float D = NDF(roughness, H, N);
                float G = GeometrySmith(N, V, L, roughness);
                Vector3f F = fresnelSchlick(F0, V, H);

                Vector3f numerator = D * G * F;
                float denominator = std::max((4.f * NdotL * NdotV), 0.001f);
                Vector3f specular = numerator / denominator;
                Vector3f kd = Vector3f(1.f) - F;

                Vector3f diffuse = Kd * 1.f / M_PI;
                return kd * diffuse * (1 - metallic) + specular;
                //return kd * diffuse + specular;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MIRROR:
        {
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                float kr;
                fresnel(wi, N, ior, kr);
                Vector3f mirror = 1.f / cosalpha;
                return kr * mirror;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        //case TRANSPARENT:  // BSDF (reflect and refract)
        //{
        //    Vector3f V = -wi;
        //    Vector3f L = wo;
        //    Vector3f H = (V + L).normalized();

        //    float NdotV = std::max(dotProduct(N, V), 0.f);
        //    float NdotL = std::max(dotProduct(N, L), 0.f);

        //    Vector3f F0(0.04f);
        //    F0 = mix(F0, Kd, metallic);

        //    float NDF = DistributionGGX(N, H, roughness);
        //    float G = GeometrySmith(N, V, L, roughness);
        //    Vector3f F = fresnelSchlick(F0, V, H);
        //    /*float F;
        //    fresnel(V, N, ior, F);*/

        //    Vector3f numerator = NDF * G * F;
        //    //float numerator = NDF * G * F;
        //    float denominator = std::max((4.f * std::abs(NdotL) * std::abs(NdotV)), 0.001f);
        //    //float specular = numerator / denominator;
        //    Vector3f specular = numerator / denominator;    // brdf

        //    // bsdf
        //    Vector3f ht = 
        //}
    }
}

#endif //RAYTRACING_MATERIAL_H