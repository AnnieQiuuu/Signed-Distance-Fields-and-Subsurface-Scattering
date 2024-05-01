#version 330 core

uniform float u_Time;

uniform vec3 u_CamPos;
uniform vec3 u_Forward, u_Right, u_Up;
uniform vec2 u_ScreenDims;

// PBR material attributes
uniform vec3 u_Albedo;
uniform float u_Metallic;
uniform float u_Roughness;
uniform float u_AmbientOcclusion;
// Texture maps for controlling some of the attribs above, plus normal mapping
uniform sampler2D u_AlbedoMap;
uniform sampler2D u_MetallicMap;
uniform sampler2D u_RoughnessMap;
uniform sampler2D u_AOMap;
uniform sampler2D u_NormalMap;
// If true, use the textures listed above instead of the GUI slider values
uniform bool u_UseAlbedoMap;
uniform bool u_UseMetallicMap;
uniform bool u_UseRoughnessMap;
uniform bool u_UseAOMap;
uniform bool u_UseNormalMap;

// Image-based lighting
uniform samplerCube u_DiffuseIrradianceMap;
uniform samplerCube u_GlossyIrradianceMap;
uniform sampler2D u_BRDFLookupTexture;

// Varyings
in vec2 fs_UV;
out vec4 out_Col;

const float PI = 3.14159f;

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct BSDF {
    vec3 pos;
    vec3 nor;
    vec3 albedo;
    float metallic;
    float roughness;
    float ao;
    float thinness;
    float subsurfaceGlow;
    float subsurfaceScale;
    float subsurfaceDistortion;

};

struct MarchResult {
    float t;
    int hitSomething;
    BSDF bsdf;
};

struct SmoothMinResult {
    float dist;
    float material_t;
};

float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

float sceneSDF(vec3 query);

vec3 SDF_Normal(vec3 query) {
    vec2 epsilon = vec2(0.0, 0.001);
    return normalize( vec3( sceneSDF(query + epsilon.yxx) - sceneSDF(query - epsilon.yxx),
                            sceneSDF(query + epsilon.xyx) - sceneSDF(query - epsilon.xyx),
                            sceneSDF(query + epsilon.xxy) - sceneSDF(query - epsilon.xxy)));
}

float SDF_Sphere(vec3 query, vec3 center, float radius) {
    return length(query - center) - radius;
}

float SDF_Box(vec3 query, vec3 bounds ) {
  vec3 q = abs(query) - bounds;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float SDF_RoundCone( vec3 query, vec3 a, vec3 b, float r1, float r2) {
  // sampling independent computations (only depend on shape)
  vec3  ba = b - a;
  float l2 = dot(ba,ba);
  float rr = r1 - r2;
  float a2 = l2 - rr*rr;
  float il2 = 1.0/l2;

  // sampling dependant computations
  vec3 pa = query - a;
  float y = dot(pa,ba);
  float z = y - l2;
  float x2 = dot2( pa*l2 - ba*y );
  float y2 = y*y*l2;
  float z2 = z*z*l2;

  // single square root!
  float k = sign(rr)*rr*rr*x2;
  if( sign(z)*a2*z2>k ) return  sqrt(x2 + z2)        *il2 - r2;
  if( sign(y)*a2*y2<k ) return  sqrt(x2 + y2)        *il2 - r1;
                        return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float smooth_min( float a, float b, float k ) {
    float h = max(k - abs(a - b), 0.0) / k;
    return min(a, b) - h * h * k * 0.25;
}

SmoothMinResult smooth_min_mix( float a, float b, float k ) {
    float h = max( k-abs(a-b), 0.0 )/k;
    float m = h*h*0.5;
    float s = m*k*0.5;
    if(a < b) {
        return SmoothMinResult(a-s,m);
    }
    return SmoothMinResult(b-s,1.0-m);
}
vec3 repeat(vec3 query, vec3 cell) {
    return mod(query + 0.5 * cell, cell) - 0.5 * cell;
}

float subtract(float d1, float d2) {
    return max(d1, -d2);
}

float opIntersection( float d1, float d2 ) {
    return max(d1,d2);
}

float opOnion(float sdf, float thickness ) {
    return abs(sdf)-thickness;
}

vec3 rotateX(vec3 p, float angle) {
    angle = angle * 3.14159 / 180.f;
    float c = cos(angle);
    float s = sin(angle);
    return vec3(p.x, c * p.y - s * p.z, s * p.y + c * p.z);
}

vec3 rotateZ(vec3 p, float angle) {
    angle = angle * 3.14159 / 180.f;
    float c = cos(angle);
    float s = sin(angle);
    return vec3(c * p.x - s * p.y, s * p.x + c * p.y, p.z);
}

vec3 rotateY(vec3 p, float angle) {
    angle = angle * 3.14159 / 180.f;
    float c = cos(angle);
    float s = sin(angle);
    return vec3(c * p.x + s * p.z, p.y, -s * p.x + c * p.z);
}

float SDF_Stache(vec3 query) {
    float left = SDF_Sphere(query / vec3(1,1,0.3), vec3(0.2, -0.435, 3.5), 0.1) * 0.1;
    left = min(left, SDF_Sphere(query / vec3(1,1,0.3), vec3(0.45, -0.355, 3.5), 0.1) * 0.1);
    left = min(left, SDF_Sphere(query / vec3(1,1,0.3), vec3(0.7, -0.235, 3.5), 0.09) * 0.1);
    left = subtract(left, SDF_Sphere(rotateZ(query, -15) / vec3(1.3,1,1), vec3(0.3, -0.1, 1.), 0.35));

    float right = SDF_Sphere(query / vec3(1,1,0.3), vec3(-0.2, -0.435, 3.5), 0.1) * 0.1;
    right = min(right, SDF_Sphere(query / vec3(1,1,0.3), vec3(-0.45, -0.355, 3.5), 0.1) * 0.1);
    right = min(right, SDF_Sphere(query / vec3(1,1,0.3), vec3(-0.7, -0.235, 3.5), 0.09) * 0.1);
    right = subtract(right, SDF_Sphere(rotateZ(query, 15) / vec3(1.3,1,1), vec3(-0.3, -0.1, 1.), 0.35));

    return min(left, right);
}

float SDF_Wahoo_Skin(vec3 query) {
    // head base
    float result = SDF_Sphere(query / vec3(1,1.2,1), vec3(0,0,0), 1.) * 1.1;
    // cheek L
    result = smooth_min(result, SDF_Sphere(query, vec3(0.5, -0.4, 0.5), 0.5), 0.3);
    // cheek R
    result = smooth_min(result, SDF_Sphere(query, vec3(-0.5, -0.4, 0.5), 0.5), 0.3);
    // chin
    result = smooth_min(result, SDF_Sphere(query, vec3(0.0, -0.85, 0.5), 0.35), 0.3);
    // nose
    result = smooth_min(result, SDF_Sphere(query / vec3(1.15,1,1), vec3(0, -0.2, 1.15), 0.35), 0.05);
    return result;
}


float SDF_Wahoo_Hat(vec3 query) {
    float result = SDF_Sphere(rotateX(query, 20) / vec3(1.1,0.5,1), vec3(0,1.65,0.4), 1.);
    result = smooth_min(result, SDF_Sphere((query - vec3(0,0.7,-0.95)) / vec3(2.5, 1.2, 1), vec3(0,0,0), 0.2), 0.3);
    result = smooth_min(result, SDF_Sphere(query / vec3(1.5,1,1), vec3(0, 1.3, 0.65), 0.5), 0.3);

    float brim = opOnion(SDF_Sphere(query / vec3(1.02, 1, 1), vec3(0, -0.15, 1.), 1.1), 0.02);

    brim = subtract(brim, SDF_Box(rotateX(query - vec3(0, -0.55, 0), 10), vec3(10, 1, 10)));

    result = min(result, brim);

    return result;
}


float SDF_Wahoo(vec3 query) {
    // Flesh-colored parts
    float result = SDF_Wahoo_Skin(query);
    // 'stache parts
    result = min(result, SDF_Stache(query));
    // hat
    result = min(result, SDF_Wahoo_Hat(query));

    return result;
}

BSDF BSDF_Wahoo(vec3 query) {
    // Head base
    //deafult to skin color
    BSDF result = BSDF(query, normalize(query), pow(vec3(239, 181, 148) / 255., vec3(2.2)),
                       0., 0.7, 1., 0.0, 0.0, 0.0, 0.0);

    result.nor = SDF_Normal(query);

    float skin = SDF_Wahoo_Skin(query);
    float stache = SDF_Stache(query);
    float hat = SDF_Wahoo_Hat(query);

    if(stache < skin && stache < hat) {
        //if stache is nearest then brown
        result.albedo = pow(vec3(68,30,16) / 255., vec3(2.2));
    }
    if(hat < skin && hat < stache) {
        //if hat is nearest then red
        result.albedo = pow(vec3(186,45,41) / 255., vec3(2.2));
    }

    return result;
}

//--------------------------Replicate the Shader Toy Scene------------------------------
//reference to https://www.shadertoy.com/view/WdBcRh

#define MAX_STEPS 300
#define SURF_DIST 1e-3
#define MAX_DIST 100.

float hash21(vec2 p) {
    p = fract(p * vec2(233.34, 851.74));
    p += dot(p, p + 23.45);
    return fract(p.x * p.y);
}
vec2 hash22(vec2 p) {
    float k = hash21(p);
    return vec2(k, hash21(p + k));
}
float SD_Sphere(vec3 p, float s)
{
    return length(p) - s;
}

mat3x3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);

    return mat3x3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

float opSmoothUnion(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

float opSmoothSubtraction(float d1, float d2, float k) {
    float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
    return mix(d2, -d1, h) + k * h * (1.0 - h);
}

float SDF_RoundBox(vec3 p, vec3 b, float r)
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - r;
}

float SDF_Pig(vec3 p,float jump) {
    p*= 1.0 + vec3(-0.2,0.2,-0.2)*(0.5+0.5*sin(u_Time*8.0+3.5));
    vec3 j = vec3(.0, -jump, .0);
    p.x = abs(p.x);
    float g = opSmoothUnion(SDF_RoundBox(p+j, vec3(1.), 0.1), SD_Sphere(p + j, 1.2), 1.); //Main Body
    g = min(g,
            opSmoothUnion(
                SDF_RoundBox(p - vec3(0, -0.25, 0.9) + j, vec3(0.4, 0.3, 0.5), 0.1),
                SD_Sphere(p - vec3(0, -0.25, 0.9) + j, .5), .5) //nose
           );
    float s = SDF_RoundBox(p - vec3(.2, -0.25, 1.5) + j, vec3(0.03, 0.13, 0.2), 0.05); //nostrile
    s = min(s, SDF_RoundBox(p - vec3(.4, 0.5, 1.3) + j, vec3(0.05, 0.2, 0.05), 0.05)); //eye
    return opSmoothSubtraction(s, g, 0.02);
}

float SDF_Bridge(vec3 p, float t) {
    float gap = 2.4;
    float tread = min(mod(t, 3.141529 * 2.) / 3.141529, 1.) * gap;
    float backScale = smoothstep(3.141529 * 2., 3.141529, mod(t, 3.141529 * 2.));
    float frontScale = smoothstep(0., 3.141529, mod(t, 3.141529 * 2.));
    float g = min(
        SDF_RoundBox(p - vec3(0., -2.3 - ((1. - backScale) * 3.), gap * -1. - tread), vec3(backScale), 0.1),
        SDF_RoundBox(p - vec3(0., -2.3, 0. - tread), vec3(1.), 0.1)
    );
    g = min(g, SDF_RoundBox(p - vec3(0., -2.3, gap - tread), vec3(1.), 0.1));
    float alternate = mod(floor(t / (3.141529 * 2.)), 2.);
    p = (rotateY(alternate > 0.5 ? (frontScale - 1.) : (1. - frontScale))* p);
    return min(g, SDF_RoundBox(p - vec3(0., -2.3, gap * 2. - tread), vec3(frontScale), 0.1));
}

float SDF_All(vec3 query) {

    //how fast it moving
    float t = u_Time * 8.;
    float pig = SDF_Pig(query, max(sin(u_Time * 8.), .0));
    float bridge = SDF_Bridge(query, t);
    float result = min(pig,bridge);

    return result;
}

float CalculateAO(vec3 p, vec3 n) {
    float d = 0.6;
    return smoothstep(0.,d,SDF_All(p + n*d));
}

//for beautiful color
vec3 increaseSaturation(vec3 color, float saturationFactor) {
    float grey = dot(color, vec3(0.333));
    return mix(vec3(grey), color, saturationFactor);
}


//-----------------------SDF Pingu----------------------------

float SDF_Ellipsoid(vec3 query, vec3 center, vec3 r )
{
    float k0 = length((query - center)/r);
    float k1 = length((query - center)/(r*r));
    return k0*(k0-1.0)/k1;
}

//reference from https://www.shadertoy.com/view/Xds3zN
float SD_RoundCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;

    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot2( pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float SDF_Pingu_Head(vec3 query){
    //head
    float head = SDF_Sphere(query / vec3(1.15,1,1), vec3(0,1.3,0), 1.) * 1.1;
    return head;
}

float SDF_Pingu_Mouth(vec3 query){
    return SDF_Sphere(query / vec3(1.2,1,1), vec3(0, 1.1, 1), 0.35);
}

float SDF_Pingu_Body(vec3 query){
    //main body
    float body = SDF_Ellipsoid(query - vec3(0.0, 0.0, 0.0), vec3(0, -1, 0),vec3(1.3,1.6,1.));
    float button = SDF_Sphere(query / vec3(1.4,1.3,1.), vec3(0, -1., 0), 1.);
    body = smooth_min(body,button,0.08);

    // wing
    float leftWing = SDF_Ellipsoid(rotateZ(query,32) - vec3(-1.5, -0.5, 0.0), vec3(.6, -1., 0), vec3(0.2, 1.2, 0.4));
    float rightWing = SDF_Ellipsoid(rotateZ(query,-32) - vec3(1.5, -0.5, 0.0), vec3(-.6, -1.,0), vec3(0.2, 1.2, 0.4));

    // combine all
    float completeBody = smooth_min(smooth_min(body, leftWing, 0.05), rightWing, 0.1);
    return completeBody;
}

float SDF_Pingu_Belly(vec3 query){
    float belly = SDF_Ellipsoid(query - vec3(0.0, 0.0, 0.0), vec3(0, -1.1, 0.4),vec3(1.05,1.25,0.9));
    return belly;
}

float SDF_Pingu_Belly2(vec3 query){
    float belly = SDF_Ellipsoid(query - vec3(0.0, 0.0, 0.0), vec3(0, -1.1, 0.4),vec3(1.1,1.3,0.8));
    return belly;
}

float SDF_Pingu_Foot(vec3 query){
    //foot
    float leftFoot = SDF_Ellipsoid(query - vec3(-0.2, -0.9, 0.0), vec3(-0.6, -1.5, 0.8),vec3(0.7,0.2,0.7));
    float rightFoot = SDF_Ellipsoid(query - vec3(0.2, -0.9, 0.0), vec3(0.6, -1.5, 0.8),vec3(0.7,0.2,0.7));
    return min(leftFoot,rightFoot);
}

float SDF_Pingu_Eyes(vec3 query) {
    float leftEye = SDF_Sphere(rotateY(query,30)/ vec3(2.9,2.9,0.98), vec3(-0.04, 0.48, 1.01), 0.08);
    float rightEye = SDF_Sphere(rotateY(query,-30)/ vec3(2.9,2.9,0.98), vec3(0.04, 0.48, 1.01), 0.08);

    return min(leftEye, rightEye);
}

float SDF_Pingu_Eyeballs(vec3 query) {
    float expressionTime = sin(u_Time *0.5);
    float leftEye = SDF_Sphere(rotateY(query,30)/ vec3(1.8,1.8,0.98), vec3(-0.04, 0.75, 1.05), 0.08);
    float rightEye = SDF_Sphere(rotateY(query,-30)/ vec3(1.8,1.8,0.98), vec3(0.04, 0.75, 1.05), 0.08);

    if(expressionTime > 0.3 && expressionTime < 0.8){
        leftEye = SDF_Sphere(rotateY(query,30)/ vec3(1.,1.,0.98), vec3(-0.04, 1.5, 1.05), 0.08);
        rightEye = SDF_Sphere(rotateY(query,-30)/ vec3(1.,1.,0.98), vec3(0.04, 1.5, 1.05), 0.08);
    }

    return min(leftEye, rightEye);
}


float SDF_Pingu(vec3 query){
    float drop = 1e10;
    float expressionTime = sin(u_Time * 0.5);
    if(expressionTime > 0.3 && expressionTime < 0.8){
        drop = SD_RoundCone(rotateZ(query,26)-vec3( 0.0,2., 0.8), vec3(0.1,0.0,0.0), vec3(-0.1,0.35,0.1), 0.15, 0.05);
    }
    float result = SDF_Pingu_Head(query);
    result = smooth_min(result,drop,0.01);
    result = smooth_min(result, SDF_Pingu_Eyes(query),0.01);
    result = smooth_min(result, SDF_Pingu_Eyeballs(query),0.001);
    result = smooth_min(result, SDF_Pingu_Mouth(query),0.05);
    result = smooth_min(result, SDF_Pingu_Body(query),0.1);
    result = min(result,SDF_Pingu_Belly(query));
    result = smooth_min(result,SDF_Pingu_Belly2(query),0.05);
    result = smooth_min(result,SDF_Pingu_Foot(query),0.05);
    return result;
}

BSDF BSDF_Pingu(vec3 query) {
    // Head base
    //deafult to black color
    BSDF result = BSDF(query, normalize(query), pow(vec3(0, 0, 0) / 255., vec3(2.2)),
                       0., 0.7, 1., 0.0, 0.0, 0.0, 0.0);

    result.nor = SDF_Normal(query);

    float head = SDF_Pingu_Head(query);
    float eyes = SDF_Pingu_Eyes(query);
    float eyeballs = SDF_Pingu_Eyeballs(query);
    float mouth = SDF_Pingu_Mouth(query);
    float belly = SDF_Pingu_Belly(query);
    float belly2 = SDF_Pingu_Belly2(query);
    float body = SDF_Pingu_Body(query);
    float foot = SDF_Pingu_Foot(query);
    float drop = SD_RoundCone(rotateZ(query,26)-vec3( 0.0,2., 0.8), vec3(0.1,0.0,0.0), vec3(-0.1,0.35,0.1), 0.15, 0.05);

    float expressionTime = sin(u_Time * 0.5);
    if(expressionTime > 0.3 && expressionTime < 0.8){
        if(drop < head && drop < body && drop < belly && drop < foot && drop < eyes && drop < eyeballs) {
            result.albedo = pow(vec3(63, 183, 252) / 255., vec3(2.2));
        }
    }

    if(mouth < head && mouth < body && mouth < belly && mouth < foot && mouth < eyes && mouth < eyeballs && mouth < drop) {

        result.albedo = pow(vec3(184,38,13) / 255., vec3(2.2));
    }

    //blink
    float blinkTime = sin(u_Time * 0.8);
    if(eyes < head && eyes < body && eyes < belly && eyes < foot && eyes < mouth && eyes < eyeballs && eyes < drop) {

        if(blinkTime > 0.1 && blinkTime < 0.2){
             result.albedo = pow(vec3(0,0,0), vec3(2.2));
        }else{
            result.albedo = pow(vec3(1,1,1), vec3(2.2));
        }
    }

    if(eyeballs < head && eyeballs < body && eyeballs < belly && eyeballs < foot && eyeballs < mouth && eyeballs < eyes && eyeballs < drop) {
        result.albedo = pow(vec3(0,0,0), vec3(2.2));
    }

    if(belly < head && belly < body && belly < mouth && belly < foot && belly < eyes && belly < eyeballs && belly < belly2) {

        result.albedo = pow(vec3(1,1,1), vec3(2.2));
    }

    if(belly2 < head && belly2 < body && belly2 < mouth && belly2 < foot && belly2 < eyes && belly2 < eyeballs && belly2 < belly) {

        result.albedo = pow(vec3(252, 230, 63) / 255., vec3(2.2));
    }

    if(foot < head && foot < body && foot < mouth && foot < belly && foot < belly2) {

        result.albedo = pow(vec3(235, 110, 21) / 255., vec3(2.2));
    }

    return result;
}


//--------------------------Duplicate the Scene------------------------------
float opRepetitionPig(vec3 p, vec3 s) {
    vec3 q = p - s * round(p / s);
    return SDF_All(q);
}


float opRepetition(vec3 p, vec3 s) {
    vec3 q = p - s * round(p / s);
//    return SDF_All(q);
    return SDF_Pingu(q);
}


vec3 opRepetition2(vec3 p, vec3 s) {
    return  p - s * round(p / s);
}

vec3 getCellIndex(vec3 p, vec3 s) {
    return floor(p / s);
}

float randomFromVec3(vec3 v) {
    return fract(sin(dot(v, vec3(12.9898, 78.233, 54.145))) * 43758.5453);
}


//------------------------SDF BSDF----------------------------------
// Random material
void assignMaterialProperties(vec3 index, out float metallic, out float roughness) {
    float randValue = randomFromVec3(index);
    metallic = randValue;
    roughness = (1.0 - randValue);
}


//helper function for color changing
vec3 changeColor() {
    float value = fract(u_Time);
    return vec3(value, 1.0 - value, sin(u_Time * 2.0 * 3.14159));
}


BSDF BSDF_Pig(vec3 query, vec3 c) {
    BSDF result = BSDF(query, normalize(query), pow(vec3(0.,0.,0.) / 255., vec3(2.2)),
                      0., 0, 1., 0.0, 0.0, 0.0, 0.0);
    vec3 norm = SDF_Normal(query);
    float ao = CalculateAO(query, norm);
    result.nor = norm;
    result.ao = ao;

    vec3 repetQuery = opRepetition2(query, c);
//    vec3 repetQuery = query;
    vec4 col = vec4(0.);
    vec3 light1Dir = normalize(vec3(.8, 1, .2));
    vec3 light1Color = vec3(1, 0.9, 0.9);
    float ground = smoothstep(-1.18, -1.19, repetQuery.y);
    col.rgb = mix(vec3(1, .7, .8), vec3(0.5, 0.6, 0.9), ground);
    col.rgb += pow(clamp(dot(reflect(normalize(repetQuery - u_CamPos), norm), light1Dir),0.0,0.1), .6) * light1Color * 0.3;
    col.rgb += norm * 0.15;
    col.rgb *= ao * 0.4 + 1;
    col.rgb = increaseSaturation(col.rgb, 1.5);
    result.albedo = pow(vec3(col.rgb),vec3(2.2));

    vec3 spacing = vec3(10.0+clamp(cos(query.z * sin(u_Time*0.01)),-1.0,1.0),10 + clamp(sin(query.z*sin(u_Time*0.01)),-1.0,1.0),10.0);
    vec3 cellIndex = getCellIndex(query+vec3(5), spacing);
    float metallic, roughness;
    assignMaterialProperties(cellIndex, result.metallic, result.roughness);

    return result;
}



float sceneSDF(vec3 query) {
# if 0
    //For Pingu SDF
    return SDF_Pingu(query);
# endif

# if 1
    //For Pig SDF
    vec3 spacing = vec3(10.0+clamp(cos(query.z * sin(u_Time*0.01)),-1.0,1.0),10 + clamp(sin(query.z*sin(u_Time*0.01)),-1.0,1.0),10.0);
    return opRepetitionPig(query, spacing);
# endif
}



BSDF sceneBSDF(vec3 query,vec3 repeatPeriod) {
# if 0
    //For Pingu BSDF
    return BSDF_Pingu (query);
# endif

# if 1
    //For Pig BSDF
    return BSDF_Pig (query, repeatPeriod);
# endif

 # if 0
    //For subsurface Sphere
    return BSDF(query, normalize(query), vec3(1, 1, 1),
                0., 0.2, 1., 0.0, 0.0, 0.0, 0.0);
# endif
}
