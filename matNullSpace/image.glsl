const int MaxNumClasses;
const int MaxNumExamples;

uniform vec4 Examples[MaxNumExamples];
uniform int Labels[MaxNumExamples];
uniform vec4 Weights[MaxNumClasses];
uniform vec3 Colors[MaxNumClasses];
uniform float weightToClassMap[MaxNumClasses];

uniform vec2 translation;
uniform int numberOfClasses;
uniform int numberOfExamples;

const float Gamma = 1.81;
const float Bound = 2.5;
const float Lens = 2.2;

const float Shininess = 23.;
const float PlaneOpacity = .80;
const float ExampleOpacity = .86;
const float ExampleSphereRadius = 0.1;

#define COLORFUL_ISOLINES

vec3 Light0;
vec3 Light1;
const vec3 SKY = vec3(.85);

// not feature toggles
const int EXAMPLE = 0;
const int PLANE = 1;

out vec4 fragColor;

struct Intersection {
	// What type of object did we intersect
	int objectID;

	// Which two weights define the plane we've intersected
	ivec2 ids;

	float t;

	float opacity;
};

void makeRay(out vec3 ro, out vec3 rd) {
	vec2 pixelCoord = gl_FragCoord.xy / resolution.xy;
	pixelCoord = pixelCoord * 2.0 - 1.0;
	pixelCoord.x *= resolution.x / resolution.y;

	rd = vec3(pixelCoord, Lens);
	rd = normalize(rd);

	ro = vec3(0.0, 0.0, -7.0);

	float cx = cos(10. * translation.x);
	float cy = cos(10. * translation.y);
	float sx = sin(10. * translation.x);
	float sy = sin(-10. * translation.y);
	mat3 rotate = mat3(cx, -sx * sy, sx * cy,
	                   0., cy, sy,
	                   -sx, -cx * sy, cx * cy);
	rd *= rotate;
	ro *= rotate;
}

// Plane-Ray intersection
float plane(in vec3 ro, in vec3 rd, in vec3 norm, in vec3 pointInPlane) {
	return dot(pointInPlane - ro, norm) / dot(rd, norm);
}

// Sphere-Ray intersection
vec3 sphere(in vec3 ro, in vec3 rd, in vec3 center, in float radius) {
	vec3 oc = ro - center;
	float b = dot(oc, rd);
	float c = dot(oc, oc) - radius * radius;
	float h = b * b - c;
	return vec3(-b - sqrt(max(h, 0.)), -b + sqrt(max(h, 0.)), h);
}

// TODO inline?
int maxWeightAtPoint(in vec3 p) {
	// Find the class of the point
	int max_weight = -1;
	float score = 0.;
	float max_score = -1e7;
	for (int i = 0; i < numberOfClasses; ++i) {
		score = dot(vec4(p, 1.0), Weights[i]);
		if (score > max_score) {
			max_score = score;
			max_weight = i;
		}
	}
	return max_weight;
}

#define withinBound(t) (t > BoundInterval.x && t < BoundInterval.y)

// More than enough intersections for all the planes and examples
Intersection intersections[MaxNumClasses * (MaxNumClasses - 1) / 2];
int numIntersections = 0;

// Fill the intersections array
bool intersect(in vec3 ro, in vec3 rd) {
	// Intersect bounding sphere
	vec2 BoundInterval;
	vec3 sphereIntersect = sphere(ro, rd, vec3(0.), Bound);
	// float nearestSphereRayDist = sqrt(max(Bound*Bound-sphereIntersect.z, 0.0))-Bound;
	// if (nearestSphereRayDist < 0.01)
	if (sphereIntersect.z > 0.) {
		// intersections[numIntersections].objectID = BOUND;
		// // intersections[numIntersections].opacity = 1.-clamp(nearestSphereRayDist/0.02, 0., 1.);
		// intersections[numIntersections].opacity = abs(sphereIntersect.y - sphereIntersect.x);
		// intersections[numIntersections].t = sphereIntersect.x;
		// numIntersections++;
		BoundInterval = sphereIntersect.xy;
	} else {
		return false;
	}

	// Find intersections with examples
	for (int i = 0; i < 60; ++i) {
		sphereIntersect = sphere(ro, rd, Examples[i].xyz, ExampleSphereRadius);

		// float nearestSphereRayDist = sqrt(max(ExampleSphereRadius*ExampleSphereRadius-sphereIntersect.z, 0.0))-ExampleSphereRadius;
		// const float closeEnough = ExampleSphereRadius/5.;
		//if (nearestSphereRayDist < closeEnough)

		if (sphereIntersect.z > 0.) {
			// checking for positive is needed when tracing the shadows??
			if (withinBound(sphereIntersect.x) && sphereIntersect.x > 0.) {
				intersections[numIntersections].objectID = EXAMPLE;

				intersections[numIntersections].ids.x = i;
				intersections[numIntersections].t = sphereIntersect.x;

				// intersections[numIntersections].opacity = 1.-clamp(nearestSphereRayDist/closeEnough, 0., 1.);
				// intersections[numIntersections].opacity *= ExampleOpacity;
				intersections[numIntersections].opacity = clamp(ExampleOpacity * (sphereIntersect.y - sphereIntersect.x) / ExampleSphereRadius, 0., 1.);

				numIntersections++;
			}
		}
	}

	// Find intersections with planes
	//
	// f(x,y,z) is a coefficient matrix for a 3d function formed by the difference between two weight functions
	// The solutions to f(x,y,z) == 0 form a plane. The plane can be parameterized by solving for z and forming
	// the funtion h(x,y). Similarly for two planes the level set of their difference h0(x,y) - h1(x,y) == 0
	// forms a line which can be parameterized by solving for y and forming the function g(x). Given this
	// function g(x) we can perform antialiasing on the seams where two planes intersect
	vec4 f;
	vec3 P;
	float t;
	int max_weight;
	// Find points on ray where any two scores are equivalent
	for (int i = 0; i < numberOfClasses; ++i) {
		for (int j = numberOfClasses - 1; j > i; --j) {
			// Construct the coefficients for the linear function f(x,y,z)
			f = Weights[i] - Weights[j];
			// Let P be the z intercept of f(x,y,z) == 0
			P = vec3(0., 0., -f.w / f.z);
			// Solve for t such that ro+rd*t is on the plane. normailzation of f.xyz (plane normal) cancels out
			t = dot(P - ro, f.xyz) / dot(rd, f.xyz);
			// We're only interested in intersections where the two score functions are maximum.
			// checking for positive is needed when tracing the shadows??
			if (withinBound(t) && t > 0.) {
				// Try to ensure that intersections[k].ids.x holds the id of the max weight in the space before
				// the ray hits the plane -> so find the maxWeight slightly before the intersection
				max_weight = maxWeightAtPoint(ro + rd * (t - 0.0001));
				if (max_weight == j || max_weight == i) {
					intersections[numIntersections].objectID = PLANE;
					if (max_weight == j) {
						intersections[numIntersections].ids.x = j;
						intersections[numIntersections].ids.y = i;
					} else {
						intersections[numIntersections].ids.x = i;
						intersections[numIntersections].ids.y = j;
					}
					intersections[numIntersections].t = t;
					intersections[numIntersections].opacity = PlaneOpacity;
					numIntersections++;
				}
			}
		}
	}
	// Sort intersections. Intersection closest to ro is in intersections[0]
	for (int i = 0; i < numIntersections - 1; i++)
		for (int j = i + 1; j < numIntersections; j++) {
			if (intersections[j].t < intersections[i].t) {
				Intersection temp = intersections[i];
				intersections[i] = intersections[j];
				intersections[j] = temp;
			}
		}

	return numIntersections > 0;
}

// This shader really needs to take into consideration shadowing. Especially since the specular component
// shows up in places where there should be no light from the sun.
vec3 lighting(in vec3 hit, in vec3 rd, in vec3 norm, in vec3 diffuseColor) {
	vec3 ret;

	vec3 lr = normalize(hit - Light0);
	float diff = max(dot(norm, -lr), 0.0);
	vec3 refL = reflect(lr, norm);
	float spec = pow(max(dot(refL, -rd), 0.0), Shininess);
	vec3 diffComp = diffuseColor * diff;
	vec3 specComp = vec3(spec);
	vec3 ambientComp = diffuseColor * .3;

	// Add a flashlight with the camera
	diffComp = mix(diffComp, abs(dot(rd, norm)) * diffuseColor, .5);
	ret = diffComp + .2 * specComp + ambientComp;

	lr = normalize(hit - Light1);
	diff = max(dot(norm, -lr), 0.0);
	refL = reflect(lr, norm);
	spec = pow(max(dot(refL, -rd), 0.0), Shininess);
	diffComp = diffuseColor * diff;
	specComp = vec3(spec);

	diffComp = mix(diffComp, abs(dot(rd, norm)) * diffuseColor, .5);
	ret = mix(diffComp + .2 * specComp + ambientComp, ret, .5);

	return ret;
}

vec3 getPlaneColor(in vec3 ro, in vec3 rd, in float t, in int i, in int j) {
	vec3 f = (Weights[i] - Weights[j]).xyz;
	f = normalize(f);
	return lighting(ro + rd * t, rd, f.xyz, Colors[int(weightToClassMap[i])]);
}

// k -> kth intersection
vec3 isolines(in vec3 ro, in vec3 rd, in int k, in vec3 baseColor) {
	// AntiAlias the seams where two planes intersect
	// The goal is to find two points on this seam. A line is constructed from these points.
	// The minimum distance between the primary rays intersecting point Q and
	// the line is found and used for color mixing.
	float dist_to_seam = 1e10;
	vec4 f0;
	vec4 f1;
	vec3 g;
	vec3 z0;
	vec3 z1;
	int idOfPlane;
	float s;
	float t;
	vec3 Q;
	vec2 line;
	vec3 lineP0;
	vec3 lineP1;
	vec3 lineDir;
	Q = ro + rd * intersections[k].t;
	f0 = Weights[intersections[k].ids.x] - Weights[intersections[k].ids.y];
	// A point on the z0 plane is then vec3(s, t, dot(z0, vec3(s,t,1.)))
	z0 = -vec3(f0.x, f0.y, f0.w) / f0.z;
	for (int i = 0; i < numberOfClasses; ++i) {
		//TODO project into screen plane and then take distances in order to smooth out the pixel mixing gradients

		if (i == intersections[k].ids.y)
			continue;

		f1 = Weights[intersections[k].ids.x] - Weights[i];

		// Coefficients for equation of plane
		z1 = -vec3(f1.x, f1.y, f1.w) / f1.z;

		g = z0 - z1;

		// A point on the line is vec2(s, dot(line, vec2(s, 1.)))
		line = -vec2(g.x, g.z) / g.y;

		s = 0.;
		t = dot(line, vec2(s, 1.));
		lineP0 = vec3(s, t, dot(z0, vec3(s, t, 1.)));

		s = 1.;
		t = dot(line, vec2(s, 1.));
		lineP1 = vec3(s, t, dot(z0, vec3(s, t, 1.)));

		lineDir = normalize(lineP1 - lineP0);

		// Point on line such that point-Q is orthogonal to line
		lineP1 = lineP0 + lineDir * dot(lineDir, Q - lineP0);
		float D = length(Q - lineP1);
		if (D < dist_to_seam) {
			dist_to_seam = D;
			idOfPlane = i;
		}
	}

	vec3 colToLerp = getPlaneColor(ro, rd, intersections[k].t, intersections[k].ids.x, idOfPlane);
#ifdef COLORFUL_ISOLINES
	// return mix(baseColor, colToLerp, .5-.5*smoothstep(0., 20./resolution.x, dist_to_seam));
	return mix(baseColor, colToLerp, .5 - .5 * smoothstep(0., .022, dist_to_seam));
#else
	// return mix(baseColor, vec3(0.), .5-.5*smoothstep(0., 20./resolution.x, dist_to_seam));
	return mix(baseColor, .2 * colToLerp, .5 - .5 * smoothstep(0., .022, dist_to_seam));
#endif
}

vec3 getColor(in vec3 ro, in vec3 rd) {
	vec3 color = SKY;
	vec3 objColor = vec3(0.);
	for (int i = numIntersections - 1; i >= 0; --i) {
		switch (intersections[i].objectID) {
		case PLANE:
			objColor = getPlaneColor(ro, rd, intersections[i].t, intersections[i].ids.x, intersections[i].ids.y);
			objColor = isolines(ro, rd, i, objColor);
			break;
		case EXAMPLE:
			vec3 hit = ro + rd * intersections[i].t;
			objColor = lighting(hit, rd, normalize(hit - Examples[intersections[i].ids.x].xyz), Colors[Labels[intersections[i].ids.x] - 1]);
			break;
		}
		color = mix(color, objColor, intersections[i].opacity);
	}
	return color;
}
void main() {
	vec3 ro;
	vec3 rd;
	makeRay(ro, rd);

	Light0 = 1.5 * vec3(Bound, .75, Bound);
	Light1 = 1.5 * vec3(-Bound, .75, -Bound);

	fragColor.rgb = SKY;
	if (intersect(ro, rd)) {
		vec3 color = getColor(ro, rd);
		// if (intersect(toward light))
		// {
		// 	color = shadow;
		// }
		fragColor.rgb = color;
	}
	// TODO intersect another ray from hit in the direction toward light
	// if there is no occluders then intersect() returns false

	// I dont really understand gamma i guess.
	fragColor.rgb = pow(fragColor.rgb, vec3(1. / Gamma));
	fragColor.a = 1.;
}
