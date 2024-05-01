
// TODO add any helper functions you need here
vec3 FresnelSchlickRoughness(float cosTheta, vec3 R, float roughness){
    return R + (max(vec3(1.0 - roughness),R)-R)* pow(clamp(1.0-cosTheta,0.0,1.0),5.0);
}

vec3 metallic_plastic_LTE(BSDF bsdf, vec3 wo) {
    vec3 N = bsdf.nor;
    vec3 albedo = bsdf.albedo;
    float metallic = bsdf.metallic;
    float roughness = bsdf.roughness;
    float ambientOcclusion = bsdf.ao;

    // TODO
    //diffuse illumination combine with albedo
    vec3 diffuseIrradiance = texture(u_DiffuseIrradianceMap, N).rgb;
    vec3 diffuse = diffuseIrradiance * albedo;

    vec3 V = wo;
    vec3 R = reflect(-V, N);

    vec3 Lo = vec3(0.f);

    vec3 F0 = mix(vec3(0.04), albedo, metallic);

    vec3 F = FresnelSchlickRoughness(max(dot(N, V), 0.0), F0, roughness);

    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= (1.0 - metallic);


    const float MAX_REFLECTION_LOD = 4.0;
    vec3 glossyIrradiance = textureLod(u_GlossyIrradianceMap, R,
                                       roughness * MAX_REFLECTION_LOD).rgb;
    vec2 brdf = texture(u_BRDFLookupTexture, vec2(max(dot(N, V), 0.0), roughness)).rg;
    vec3 specular = glossyIrradiance * (F * brdf.x + brdf.y);


    vec3 ambient = (kD * diffuse + specular) * ambientOcclusion;
    vec3 color = ambient + Lo;

    return color;
}
