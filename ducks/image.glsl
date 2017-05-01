vec4 color(vec2 uv) {

	vec4 col = .7 * texture(iChannel0, uv);
	float height = col.r;
	float sound = texture(iChannel1, vec2(height, 0.)).r;
	//fragColor = vec4(smoothstep(0.3,1.3,sound));
	col = vec4(smoothstep(0.1, 1., sound));
	//fragColor = vec4(sound);
	return pow(col, vec4(2.));
}

#define r iResolution
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//fragColor = texture(iChannel0, fragCoord.xy/iResolution.xy);
	//return;

	vec2 uv = fragCoord.xy / iResolution.xy;
	fragColor = color(uv);
	fragColor += color(uv + vec2(.5 / r.x, 0.));
	fragColor += color(uv + vec2(-.5 / r.x, 0.));
	fragColor += color(uv + vec2(.5 / r.x, .5 / r.y));
	fragColor += color(uv + vec2(-.5 / r.x, .5 / r.y));
	fragColor += color(uv + vec2(.5 / r.x, -.5 / r.y));
	fragColor += color(uv + vec2(-.5 / r.x, -.5 / r.y));
	fragColor += color(uv + vec2(0., .5 / r.y));
	fragColor += color(uv + vec2(0., -.5 / r.y));
	fragColor *= 1. / 9.;
}
