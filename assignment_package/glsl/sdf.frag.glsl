
#define FOVY 45 * PI / 180.f
Ray rayCast() {
    vec2 ndc = fs_UV;
    ndc = ndc * 2.f - vec2(1.f);

    float aspect = u_ScreenDims.x / u_ScreenDims.y;
    vec3 ref = u_CamPos + u_Forward;
    vec3 V = u_Up * tan(FOVY * 0.5);
    vec3 H = u_Right * tan(FOVY * 0.5) * aspect;
    vec3 p = ref + H * ndc.x + V * ndc.y;

    return Ray(u_CamPos, normalize(p - u_CamPos));
}


#define MAX_ITERATIONS 300
#define THICKNESS 0.01
#define MAX_DISTANCE 300
MarchResult raymarch(Ray ray) {
    vec3 origin = ray.origin;
    vec3 direction = ray.direction;
    float t = 0.0;
    for(int i = 0; i < MAX_ITERATIONS; ++i) {
        vec3 samplePoint = origin + direction * t;
        vec3 repeatPeriod =vec3(10.0+clamp(cos(samplePoint.z * sin(u_Time*0.01)),-1.0,1.0),10 + clamp(sin(samplePoint.z*sin(u_Time*0.01)),-1.0,1.0),10.0);
//        vec3 repeatPeriod =vec3(10);
        float distance = sceneSDF(samplePoint);
        if (distance < THICKNESS) {
            BSDF bsdf = sceneBSDF(samplePoint, repeatPeriod);
            return MarchResult(t, 1, bsdf);
        }
        t += distance;
        if (t > MAX_DISTANCE) {
            break;
        }
    }
    //if hit nothing
    return MarchResult(-1, 0, BSDF(vec3(0.), vec3(0.), vec3(0.), 0., 0., 0., 0., 0., 0., 0.));
}

#define K_COEFF 2
#define AO_DIST 0.085
float subsurfaceThickness(vec3 pos, vec3 norm, float k) {
    float ao = 0.;
    for(float i = 0.0; i < 5.0; ++i) {
        float coeff = 1.0 / pow(2.0, i);
        ao += max(coeff * (i * AO_DIST - (sceneSDF(pos + norm * i * AO_DIST))), 0.);
    }
    return 1.0 - k * ao;
}

#define AMBIENT 0.f
vec3 subsurfaceColor(vec3 lightDir, vec3 normal, vec3 viewVec, BSDF bsdf) {
        vec3 scatterDir = lightDir + normal * bsdf.subsurfaceDistortion;
        float lightReachingEye = pow(clamp(dot(viewVec, -scatterDir), 0.0, 1.0), bsdf.subsurfaceGlow)
                                 * bsdf.subsurfaceScale;
        float attenuation = max(0.0, dot(normal, lightDir)
                             + dot(viewVec, -lightDir));
        float totalLight = attenuation * (lightReachingEye + AMBIENT) * bsdf.thinness;
        return texture(u_DiffuseIrradianceMap, lightDir).rgb * totalLight;
}



#define DISTORTION 0.2
#define GLOW 6.0
#define SCALE 3.f
//I also add those into BSDF struct
void main()
{
    Ray ray = rayCast();
    MarchResult result = raymarch(ray);
    BSDF bsdf = result.bsdf;
    vec3 pos = ray.origin + result.t * ray.direction;
    bsdf.pos = pos;
    float thinness = subsurfaceThickness(bsdf.pos,bsdf.nor,K_COEFF);
    bsdf.thinness = thinness;
    bsdf.subsurfaceGlow = GLOW;
    bsdf.subsurfaceScale = SCALE;
    bsdf.subsurfaceDistortion = DISTORTION;
    vec3 viewVec = normalize(u_CamPos - pos);
    vec3 reflectedLightDir = reflect(viewVec, normalize(bsdf.nor));
    vec3 lightDir = -normalize(reflectedLightDir);
    vec3 subsurface = subsurfaceColor(lightDir,bsdf.nor,viewVec, bsdf);

    vec3 color = metallic_plastic_LTE(bsdf, -ray.direction) +  subsurface * (1.0 - bsdf.metallic);

    // Reinhard operator to reduce HDR values from magnitude of 100s back to [0, 1]
    color = color / (color + vec3(1.0));
    // Gamma correction
    color = pow(color, vec3(1.0/2.2));

    out_Col = vec4(color, result.hitSomething > 0 ? 1. : 0.);
}

