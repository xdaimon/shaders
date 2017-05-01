#version 100

#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif

uniform sampler2D wave_texture;
varying vec2 tex_pos;

uniform float time;
uniform vec2 resolution;
uniform vec3 mouse;
uniform vec3 color1;
uniform vec3 color2;

uniform vec2 info; // x = line_scalar
                   // y = texture size

float line(vec2 p, vec2 a, vec2 b) {
	vec2 pa = p - a;
	vec2 ba = b - a;

	vec2 w_or_max = ba * clamp(dot(pa, ba) / dot(ba, ba), 0., 1.);
	return length(pa - w_or_max);
}

float f(float x) {
	vec4 t = texture2D(wave_texture, vec2(x, 0.));
	return t.r / 256. + t.a;
}

void main() {
	vec2 cp = tex_pos * 2. - 1.;
	cp.x *= resolution.x / resolution.y;

	// How many pixels does a circle of radius r intersect?
	float pixels = floor(length(resolution * (tex_pos - .5)) + .5);
	// Approximate the theta step to the next radial line
	float indx_step_to_next_radial = 1. / ((2. * 3.141592653589793) * pixels);

	float theta = atan(cp.y, cp.x);
	if (theta < 0.)
		theta += (2. * 3.141592653589793);

	// Radius of polar function f
	float r0 = f(theta / (2. * 3.141592653589793)) + .13;
	float r1 = f(theta / (2. * 3.141592653589793) + indx_step_to_next_radial) + .13;

	// Cartesian coords
	vec2 c0 = r0 * vec2(cos(theta), sin(theta));
	vec2 c1 = r1 * vec2(cos(theta + indx_step_to_next_radial), sin(theta + indx_step_to_next_radial));

	float dist = 1. / line(cp, c0, c1);
	gl_FragColor = vec4(mix(color2, color1, clamp(dist * dist / info.x, 0., 1.)), 1.);
}
