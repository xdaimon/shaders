#define time iGlobalTime / 5.
#define DX 4.
#define PLOT_RADIUS 3.
// ray marching
const int max_iterations = 1500;
const float stop_threshold = 0.001;
const float clip_far = 30.0;

// math
const float PI = 3.14159265359;
const float DEG_TO_RAD = PI / 180.0;
vec3 rot3(vec3 pos, float th) {
	return mat3(1., 0., 0., 0., cos(th), sin(th), 0., -sin(th), cos(th)) * pos;
}
vec2 rot2(vec2 pos, float th) {
	return mat2(cos(th), sin(th), -sin(th), cos(th)) * pos;
}
float dist_plane(float p) {
	return p + 5.;
}
float dist_sphere(vec4 p) {
	return length(p) - PLOT_RADIUS;
}
float smin(float a, float b) {
	return (a * exp(-a * 32.) + b * exp(-b * 32.)) / (exp(-a * 32.) + exp(-b * 32.));
}
float dist_box(vec3 p) {
	return length(max(abs(p) - vec3(PLOT_RADIUS), 0.));
}
float val(vec4 p) {
	//	p.yw=rot2(p.yw,rotate.x);
	//	p.xw=rot2(p.xw,rotate.y);
	//	p.zw=rot2(p.zw,rotate.z);
	//	p.zx=rot2(p.zx,time);
	//	p.w-=fd;
	//	p.z-=zzz;
	float x = p.x;
	float y = p.y;
	float z = p.z;
	float w = p.w;

	vec3 a = vec3(.3, .43, .1);
	vec3 b = vec3(-.9, .0, .1);
	vec3 c = vec3(.2, .04, .4);
	vec3 d = vec3(-.853, 1.14, .6);

	float score = dot(p.xyz, a);
	//*
	score = max(score, dot(p.xyz, b));
	//score = max(score, dot(p.xyz, c));
	//score = max(score, dot(p.xyz, d));
	/*/
    score  *= dot(p.xyz,b);
    score *= dot(p.xyz,c);
    score *= dot(p.xyz,d);
	*/
	return score;

	// 2D surfaces
	//return exp(8.*z - x*x - 2./3.*z*z) - y;
	//return (x*x+y*y-4.);
	//return (x*x+z-y);

	//return (sin(z)+sin(y)+sin(x+time+w));
	//return (sin(z*x- time)+sin(y+time)+sin(z-time));
	//return (sin(x+y+z));
	//return sin(x*y)/(sin(x*z)+1.2) + y; // tubular
	//return (z*y*sin(x)-2.*x+cos(y));
	//return 2.*z*(pow(w,time) - 3.*pow(x,time))*(1.-pow(y,time))+pow(pow(x,time)+pow(z,time),time) - (9.*pow(y,time)-1.)*(1.-pow(y,time));
	#define sS(X,Z) 1./sqrt((x+X)*(x+X) + (z+Z)*(z+Z) + y*y)
	return sS(-1.,-1.) + sS(1.,1.) + sS(-1.,1.) - sS(1.,-1.) - time;
	//#define SS(X,Z) 1./sqrt((x+X)*(x+X) + (z+Z)*(z+Z))
	//return SS(-1.,-1.) + SS(1.,1.) + SS(-1.,1.) - SS(1.,-1.) - y;
	//return x*x*z*z + x*x*y*y + z*z*y*y + x*y*z;
	//return (z*y*sin(x)-2.*x+cos(y*z));
	//return (4.*x+z*z*z + x*y - 4.)/10.;
	//return sin(z-time)-y;
	//float R=.5;
	//return -R*R*z + x*x*z + z*z*z - 2.*R*x*y - 2.*x*x*y - 2.*z*z*y+z*y*y;
	//float t = dot(p, p) - 2.*z - 1.;
	//return (t+4.*z)*(t*t - 8.*y*y) + 16.*x*y*t;

	//p*=fd;
	//p.zx -= sin(p.zx*1.)*.45;
	//p.xy = rot2(p.xy,p.z*sin(time)*0.3);
	//p.x += sin(p.z*3.)*.35;
	//p.y += sin(p.y*4.)*.3*sin(time);
	//p.z += p.y*sin(time*0.9)*2.+0.5;
	//return length(p)-6.;

	//return z*y*w*sin(y)-2.*x*cos(z*w);
	//return (sin(w*y)/(sin(w+z)+1.2)+y);
	//return sin(z)+sin(y)+sin(x)+sin(w);
}
vec4 v_grad(vec4 p) {
	vec2 e = vec2(0., .01);
	float w = val(p);
	return vec4(val(p + e.yxxx) - w,
	            val(p + e.xyxx) - w,
	            val(p + e.xxyx) - w, 0.) /
	       e.y;
}
float dist_field(vec4 p) {
	float plane = dist_plane(p.y);
	//float box = dist_box(p.xyz);
	float sphere = dist_sphere(p);
	float d0 = val(p);
	p += v_grad(p) / 50.;
	float d1 = val(p);
	// union     : min( d0,  d1 )
	// intersect : max( d0,  d1 )
	// subtract  : max( d1, -d0 )
	return min(plane, max(max(d0, -d1), sphere));
}
float dist_bounding(vec4 p) {
	return min(dist_plane(p.y), dist_sphere(p));
}
//#define BUG
float ray_marching(vec4 origin, vec4 dir, float start, float end) {
	float depth = start;
	float dist;
	int i = 0;
	// bounding box/sphere/plane
	for (int j = 0; j < max_iterations; j++) {
		dist = dist_bounding(origin + dir * depth) / 10.;
		if (dist < stop_threshold) {
			break;
		}
		depth += dist;
		if (depth >= end) {
#ifdef BUG
			return float(i);
#endif
			return end;
		}
		i++;
	}
	i = 0;
	// object to plot
	for (int j = 0; j < max_iterations; j++) {
		dist = dist_field(origin + dir * depth);
		if (dist < stop_threshold) {
#ifdef BUG
			return float(i);
#endif
			return depth;
		}
		//		if (dist_sphere(origin+dir*depth)>rad || depth >= end)
		if (depth >= end) {
#ifdef BUG
			return float(i);
#endif
			return end;
		}
		depth += dist / (DX + 1.);
		i++;
	}
#ifdef BUG
	return float(i);
#endif
	return end;
}
vec4 ray_dir(vec2 pos) {
	float cot_half_fov = tan(67.5 * DEG_TO_RAD);
	float z = .5 * cot_half_fov;
	return normalize(vec4(pos, -z, 0.));
}
vec4 calcNormal(vec4 p) // for shading
{
	vec2 e = vec2(.01, .0);
	float d = dist_field(p);
	vec4 nor = vec4(
	    dist_field(p + e.xyyy) - d,
	    dist_field(p + e.yxyy) - d,
	    dist_field(p + e.yyxy) - d, 0.);
	return normalize(nor);
}
float softshadow(in vec4 ro, in vec4 rd, in float mint, in float tmax) {
	float res = 1.0;
	float t = mint;
	for (int i = 0; i < 32; i++) {
		float h = dist_field(ro + rd * t);
		res = min(res, 8.0 * h / t);
		t += clamp(h, 0.02, 0.10);
		if (h < 0.001 || t > tmax)
			break;
	}
	return clamp(res, 0.0, 1.0);
}
// shading function due to IQ at shadertoy.com
vec3 shading(vec4 eye, float depth, vec4 dir) {
	// shading
	vec4 pos = eye + depth * dir;
	vec4 nor = calcNormal(pos);
	vec4 ref = reflect(dir, nor);
	vec3 col = 0.45 + 0.3 * sin(vec3(.5, .5, .5));
	vec4 lig = normalize(vec4(1., 1., -3., 0.));
	float amb = clamp(0.5 + 0.5 * nor.y, 0.0, 1.0);
	float dif = clamp(dot(nor, lig), 0.0, 1.0);
	float dom = smoothstep(-0.1, 0.1, ref.y);
	float spe = pow(clamp(dot(ref, lig), 0.0, 1.0), 16.0);
	float fre = pow(clamp(1.0 + dot(nor, dir), 0.0, 1.0), 2.0);

	//dif *= softshadow( pos, lig, 0.02, 2.5 );
	//dom *= softshadow( pos, ref, 0.02, 2.5 );

	vec3 lin = vec3(0.);
	lin += 1.20 * dif * vec3(1., .85, .55);
	lin += 1.20 * spe * vec3(1.00, 0.85, 0.55) * dif;
	lin += 0.20 * amb * vec3(0.50, 0.70, 1.00);
	lin += 0.30 * dom * vec3(0.50, 0.70, 1.00);
	lin += 0.40 * fre * vec3(1.00, 1.00, 1.00);

	// xor to add modulus two... cool
	if (pos.y < -4.99)
		col = vec3(.5) * mod(floor(mod(pos.z, 2.)) + floor(mod(pos.x, 2.)), 2.);
	return col * lin;
}
vec3 perf_color(float depth) {
	vec3 col = vec3(1.) * sin(3.141592 / 2. * depth / 1000.);
	return col;
}
void mainImage(out vec4 fragColor, in vec2 c) {
	c /= iResolution.xy;
	c -= .5;
	c.x *= iResolution.x / iResolution.y;
	vec2 m = iMouse.xy / iResolution.xy * 2. - 1.;

	float look_xz = 6.27 * m.x;
	float look_zy = 0.22;
	vec4 dir = ray_dir(c);
	dir.zy = rot2(dir.zy, look_zy);
	dir.xz = rot2(dir.xz, look_xz);
	vec4 eye = normalize(vec4(0., 4., 25., 0.)) * 15.;
	eye.xz = rot2(eye.xz, look_xz);

	float depth = ray_marching(eye, dir, 0., clip_far);
#ifndef BUG
	if (depth >= clip_far) {
		dir = eye + dir * depth + 16.;
		fragColor = vec4(0.17, 0.21, .46, 1.) * dir.y / 10.;
		return;
	}
#endif

#ifdef BUG
	fragColor.rgb = perf_color(depth) * 4.;
	return;
#endif
	fragColor = vec4(shading(eye, depth, dir), 1.);
	return;

	//if (pos.y < -4.99) color = 2.*vec4(1.)/abs(6.+pos.z*pos.z+pos.x*pos.x);
	//fragColor = color*vec4(lin,1.);
}
