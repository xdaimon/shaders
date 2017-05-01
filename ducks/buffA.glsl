#define PI 3.14159
#define R iResolution

#define IT 36
#define s 2.816

vec2 cmul(vec2 a, vec2 b) {
	return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

vec2 triangle(vec2 x, float d) {
	return abs(fract(x) * d - d / 2.);
}

vec2 rotate(vec2 a, float th) {
	float c = cos(th);
	float ss = sin(th);
	return mat2(c, ss, -ss, c) * a;
}

#define ZOOM 1.1

#define M (iMouse.xy / iResolution.xy)
vec2 formula(in vec2 z, float i) {

	/*  
	float ln = length(z);
    float th = iGlobalTime/10.;
    z = triangle(z, 5.*113./800.);
    return rotate(z-.5, th)/ln;
    //return rotate(z-iMouse.xy/iResolution.xy*2., th)/ln;
//*/

	//float m = dot(z,z)*i;
	//float m = dot(z,z)/i;
	float m = dot(z, z);
	/*
    return abs(z) / m *s-iMouse.xy/iResolution.xy*2.;
    //wow
    //return abs(z) / m * iMouse.xy/iResolution.xy - iMouse.xy/iResolution.xy/s;
/*/

	//*
	//z.y=abs(z.y);
	z = abs(z);
	return vec2(log(length(z)), atan(z.y, z.x)) + M * 2. - 1.;
	/*/ 
    
/*
    // Map points to the first quadrant using imaginary numbers
    vec2 c = iMouse.xy/iResolution.xy*2.;
    float z_phase = atan(z.y/z.x);
    z_phase = acos(cos(z_phase)) - atan(z.y,z.x);
    vec2 u = vec2(cos(z_phase), sin(z_phase));
    return cmul(z, u)/m*s - c;
//*/

	/*
    vec2 c = iMouse.xy/iResolution.xy*2.-1. + vec2(12./800., 11./450.);
    float z_phase = atan(z.y/z.x);
    z_phase = acos(cos(z_phase));
    //return (cmul(cmul(z,z), z) - c*s)/length(z*2.);
    //return 1./((cmul(cmul(z,z), z) - c)/length(z*2.));
    //return (cmul(cmul(z,z), vec2(0.,1.)) - c)/length(z*2.);
    //return (cmul(cmul(z,z),iMouse.xy/iResolution.xy) - c)/length(z*2.);
    //return ((cmul(cmul(z,z)/(cos(iGlobalTime/2.)*sin(iGlobalTime/2.)*.2-1.), vec2(1.)) + vec2(cos(iGlobalTime/2.), sin(iGlobalTime/2.))*s)/length(z));
    //return ((cmul(cmul(z,z)/(cos(iGlobalTime/2.)*sin(iGlobalTime/2.)*.2-1.), vec2(1.)) + c*s)/length(z));
//*/
}

float map(in vec2 p) {
	float f = 0.;
	float ln = 0.;
	float prevln = 0.;

	p /= ZOOM;

	for (int i = 0; i < IT; i++) {
		p = formula(p, float(i + 1));

		//f = length(p)*(34./450.)+f*(1.-34./450.);

		prevln = ln;
		ln = length(p);
		f += exp(-1. / abs(prevln - ln));
	}
	return f;
}

vec4 shade(float t) {
	const vec3 a = vec3(.5);
	const vec3 b = vec3(.5);
	const vec3 cc = vec3(1., 1., 1.25);
	vec3 d = vec3(1., 0., -0.1);

	vec3 color = a + b * cos(cc * t + d);
	return vec4(color, 1.);
}
vec2 hash(vec2 i) {
	return fract(vec2(sin(i.y * 845.589 + iDate.w) * 5865.9, sin(iDate.w + i.x * 8454.589) * 865.9));
}
#define saturate(x) clamp(x, 0., 1.)
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

	vec4 oldColor = texture(iChannel0, fragCoord.xy / iResolution.xy);
	//if (iFrame > 2000) {fragColor = oldColor; return;}
	vec2 uv = fragCoord.xy + 2. * hash(fragCoord.xy);

	//uv = 2.*uv - 1.;
	//uv *= 1.9;
	//uv.x*=iResolution.x/iResolution.y;

	// jitter
	//uv+= 0.001*hash(uv*1000.);

	uv = 2. * (2. * uv - iResolution.xy) / iResolution.y;

	float t = map(uv);
	fragColor -= fragColor;
	//fragColor = shade(t);
	//fragColor = vec4(pow(cos(.05960*t), 1.));
	fragColor = vec4(cos(.05960 * t));
	fragColor = mix(saturate(oldColor), fragColor, iMouse.z > 0. ? 1. : .01);
}
