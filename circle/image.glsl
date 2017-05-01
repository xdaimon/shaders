#version 100

#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif

uniform sampler2D fft_texture;
varying vec2 tex_pos;

uniform float time;
uniform vec2 resolution;
uniform vec3 mouse;
uniform vec3 color1;
uniform vec3 color2;

uniform vec2 info; // x : bass amplitude
                   // y : scale

float f(float x) {
	vec4 t = texture2D(fft_texture, vec2(x, 0.));
	return t.r / 256. + t.a;
}

void main() {
	vec2 cp = tex_pos * 2. - 1.;
	cp.x *= resolution.x / resolution.y;

	vec2 m = mouse.xy;
	m.x = smoothstep(0.10, 1., m.x);

	float theta = atan(cp.x, cp.y) + 0.02;
	if (theta < 0.)
		theta += (2. * 3.141592653589793);
	theta /= (2. * 3.141592653589793);

	vec3 fg = color1;
	vec3 bg = color2;

	const float base_win_size = 2. / 24.;
	const float final_win_size = 9. / 24.;
	const float win_shift = 3. / 24.;

	float look = (theta * (base_win_size + final_win_size * m.x) + win_shift * smoothstep(0.1, 1., m.x));
	float v = f(look);

	// Mix the frequencies where theta < 0.13 with the frequencies where theta > 1. so that there is no visible split.
	look = (theta + 1.) * (base_win_size + final_win_size * m.x) + win_shift * smoothstep(0.1, 1., m.x);
	float cross_over_wave = f(look);
	const float cross_over_range = 0.13;
	v = mix(cross_over_wave, v, smoothstep(0., 1., clamp(theta / cross_over_range, 0., 1.)));
	v *= info.y;

	float mouse_scale = 210. * smoothstep(0., 1., m.y - .45 * sqrt(m.x)) + 25.;
	float len = length(cp);
	float ring_density = 2.75 - sqrt(len);
	float ring_speed = -2.2 * time;
	float ring_shift = 3.;
	float zero = len - v;
	float ring = fract(ring_shift + ring_speed + ring_density * zero);
	float center_glow = v / len;
	float lerp = center_glow / abs(smoothstep(0., 1., ring) - .5) / info.x;
	gl_FragColor = vec4(mix(bg, fg, lerp / mouse_scale), 1.);
}
