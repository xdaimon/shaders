void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	fragColor -= fragColor;
	vec2 uv = fragCoord.yx / iResolution.yx;

	uv = uv * 2. - 1.;

	vec2 m = iMouse.yx / iResolution.yx * 2. - 1.;

	//uv/=atan(uv.y,uv.x)+3.14;
	//uv/=length(uv) + 0.23;

	//float ttt = atan(uv.y, uv.x)+length(5.*uv)+3.141592;
	//uv *= vec2(cos(ttt), sin(ttt));

	vec2 s;

	vec2 sum = vec2(0);

	//*
	uv.y += iGlobalTime / 5.;
	uv.y += 100.;
	if (iMouse.w > 0.)
		uv.y += floor(iGlobalTime * 5.);
	uv.x -= m.x;

	vec2 PlotScale = 3. * vec2(1., 1.);
	uv *= PlotScale;

	s = uv;

	vec2 abln;
	for (float n = 1.; n < 100.; ++n) {
		abln = s * log(n);
		abln.x = exp(abln.x);
		sum.x += cos(abln.y) / abln.x;
		sum.y += -sin(abln.y) / abln.x;
	}

	//float th = atan(sum.y, sum.x);
	//fragColor = 1.6 * vec4(.2*sin(th), .25*cos(1.5*th), .8*cos(th), 1.);

	float th = atan(sum.y, sum.x);
	fragColor += 1.6 * vec4(.2 * sin(th), .25 * cos(1.5 * th), .8 * cos(th), 1.) / length(sum);

	//float l = length(sum);
	//fragColor += vec4(smoothstep(.5-fwidth(l), 1. + fwidth(l), l));

	//float s = abs(sum.x+sum.y);
	//fragColor += vec4(smoothstep(0.-3.*fwidth(s), 1. + 3.*fwidth(s), s));

	//float es = exp(sum.x);
	//fragColor += vec4(smoothstep(0., .1,es));
	//fragColor = clamp(fragColor, 0., 1.);

	float w = abs(.5 - s.x);
	fragColor = mix(fragColor, vec4(1., 0., 0., 0.), 1. - smoothstep(0., .01, w));

	/*/

    vec2 PlotScale = 8.* vec2(1., 1.);
    uv *= PlotScale;
    s=uv.yx;

    for ( float n = 1.; n <= 1000.; ++n)
        sum.x += 1./pow(n, s.x);
    fragColor += vec4(1.-smoothstep(0., 0.1, abs(s.y-sum.x)));

//*/
}
