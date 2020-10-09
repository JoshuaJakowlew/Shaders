const int MAX_MARCHING_STEPS = 60;
const float MIN_DIST = 0.1;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;

// Axis rotation

mat3 rotateX(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(1, 0, 0),
        vec3(0, c, -s),
        vec3(0, s, c)
    );
}

mat3 rotateY(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

mat3 rotateZ(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, -s, 0),
        vec3(s, c, 0),
        vec3(0, 0, 1)
    );
}

// Per-primitive operations

float opUnion(float d1, float d2)
{
    return min(d1, d2);
}

float opSubtraction(float d1, float d2)
{
    return max(d1, -d2);
}

float opIntersection(float d1, float d2)
{
    return max(d1, d2);
}

float opSmoothUnion(float d1, float d2, float k)
{
    float h = clamp( 0.5 + 0.5 * (d2-d1) / k, 0.0, 1.0 );
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

float opSmoothSubtraction(float d1, float d2, float k)
{
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h);
}

float opSmoothIntersection(float d1, float d2, float k)
{
    float h = clamp( 0.5 - 0.5 * (d2 - d1) / k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k * h * (1.0 - h);
}

// Primitives

float sdSphere(vec3 p, float r)
{
    return length(p) - r;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz) - t.x, p.y);
  return length(q) - t.y;
}

// Scene

float sdScene(vec3 p)
{
    p = rotateY(iTime / 2.0) * p;
    
    float res = MAX_DIST;
    
    { // Box
    float box = sdBox(p, vec3(.9) + 0.1 * sin(3. * iTime)) - .1;
    
    float clipSphere = sdSphere(p, 1.5 + 0.1 * sin(3. * iTime));
    res = opUnion(res, opIntersection(box, clipSphere));
    
    float subSphere = sdSphere(p, 1.3);
    res = opSubtraction(res, subSphere);
    }
    
    { // Center ball
    float sphere = sdSphere(p, .4);
    res = opUnion(res, sphere);   
    }
    
    { // Side balls   
    float angle = 6.28 / 4.;
    float sector = round(atan(p.z, p.x) / angle);
    
    vec3 q = p;
    float an = sector * angle;
    q.xz = mat2(cos(an), -sin(an), sin(an), cos(an)) * q.xz;
    
    
    float sphere = sdSphere(q + vec3(-1.5 + 0.5 * sin(3. * iTime), 0, 0), .4 - 0.1 * sin(3. * iTime));
    res = opUnion(res, sphere);
    }
    
    // Wobble
    res *= inversesqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    // Displacement
    res += .3 * sin(sin(1. * iTime) * 25. * p.x) * .3 * sin(sin(1. * iTime) * 25. * p.y) * .3 * sin(sin(1. * iTime) * 25. * p.z);
    return res;
} 

// Raymarching

float rayMarch(vec3 eye, vec3 marchingDirection, float start, float end, int maxSteps)
{
    float depth = start;
    for (int i = 0; i < maxSteps; ++i) {
        float dist = sdScene(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}

vec3 rayDirection(float fieldOfView, vec2 size, vec2 fragCoord) {
    vec2 xy = fragCoord - size / 2.0;
    float z = size.y / tan(radians(fieldOfView) / 2.0);
    return normalize(vec3(xy, -z));
}

mat3 viewMatrix_(vec3 eye, vec3 center, vec3 up) {
    // Based on gluLookAt man page
    vec3 f = normalize(center - eye);
    vec3 s = normalize(cross(f, up));
    vec3 u = cross(s, f);
    return mat3(s, u, -f);
}

// Lightning

vec3 estimateNormal(vec3 p)
{
    float pDist = sdScene(p);
    return normalize(vec3(
        sdScene(vec3(p.x + EPSILON, p.y, p.z)) - pDist,
        sdScene(vec3(p.x, p.y + EPSILON, p.z)) - pDist,
        sdScene(vec3(p.x, p.y, p.z  + EPSILON)) - pDist
    ));
}

vec3 phongIllumination(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity)
{
    vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = clamp(dot(L, N),0.,1.); 
    //float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        // Light not visible from this point on the surface
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        // Light reflection in opposite direction as viewer, apply only diffuse
        // component
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

vec3 lightScene(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(4.0 * sin(iTime),
                          2.0,
                          4.0 * cos(iTime));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongIllumination(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity);
    
    vec3 light2Pos = vec3(2.0 * sin(0.37 * iTime),
                          2.0 * cos(0.37 * iTime),
                          2.0);
    vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongIllumination(k_d, k_s, alpha, p, eye,
                                  light2Pos,
                                  light2Intensity);    
    return color;
}

// Main

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Setup camera
    vec3 viewDir = rayDirection(45.0, iResolution.xy, fragCoord);
	vec3 eye = vec3(8.0, 5.0 * sin(0.2 * iTime), 7.0);    
    
    mat3 viewToWorld = viewMatrix_(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
    
    vec3 worldDir = viewToWorld * viewDir;
    
    // Calc SDF
    float dist = rayMarch(eye, worldDir, MIN_DIST, MAX_DIST, MAX_MARCHING_STEPS);
    
    if (dist > MAX_DIST - EPSILON)
    {
		fragColor = vec4(0);
        return;
    }
    
    // Lightning

    vec3 p = eye + dist * worldDir;
    vec3 K_a = 0.5 + 0.5 * cos(iTime + p.xyx + vec3(0, 2, 4));
    vec3 K_d = K_a;
    vec3 K_s = vec3(1.0, 1.0, 1.0);
    float shininess = 10.0;
    
    vec3 color = lightScene(K_a, K_d, K_s, shininess, p, eye);
    
    // Gamma correction
    color = pow(color, vec3(0.4545));
    fragColor = vec4(color, 1.0);
}