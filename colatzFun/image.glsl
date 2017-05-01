//#define FAST

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	vec2 p = fragCoord / 40.;

#ifdef FAST
	p.y *= 22.;
#else
	p.y *= 5.;
#endif

	p = floor(p);
	float odd = p.y * 2. + 1.;

	// Pair each odd number with a power of two
	float even = odd * pow(2., p.x);

	// Inverse of 3x + 1
	float link = (even - 1.) / 3.;

	// Color block white if link is integer which is true when even == 3*o + 1
	// for some odd integer o
	fragColor = .75 * vec4(fract(link) < 0.1);

	float t = iGlobalTime;
#ifdef FAST
	t = floor(t * 25.);
#else
	t = floor(t);
#endif
	const float evenSkip = 1.;
	float oddint = t * (1. + evenSkip) + 1.;
	// Lets look at which even number (3*oddint+1) will map to
	if (abs(oddint - link) < 0.01 && fract(link) < 0.01)
		fragColor = vec4(1., 0., 0., 1.);

#ifdef FAST
	return;
#endif

	// I wonder what expression tells when a block in the ith column is lit
	//0==mod(t+2-1,2);			lit at t == 1
	//0==mod(t+4-0,4);			lit at t == 0
	//0==mod(t+8-6,8);			lit at t == 6
	//0==mod(t+16-2,16);		lit at t == 2
	//0==mod(t+32-26,32);		lit at t == 26
	//0==mod(t+64-10,64);		lit at t == 10
	//0==mod(t+128-106,128);	lit at t == 106
	//0==mod(t+256-42,256);		lit at t == 42
	//0==mod(t+512-426,512);	lit at t == 426
	//0==mod(t+1024-170,1024);	lit at t == 170
	//0==mod(t+2048-1706,2048); lit at t == 1706
	//0==mod(t+4096-682,4096);	lit at t == 682
	// 0 6 2 6 0 6 2 6 0 6 2 6 0 6 2

	//0==mod(t+2-1,2);			lit at t == 1
	//0==mod(t+8-6,8);			lit at t == 6
	//0==mod(t+32-26,32);		lit at t == 26
	//0==mod(t+128-106,128);	lit at t == 106
	//0==mod(t+512-426,512);	lit at t == 426
	//0==mod(t+2048-1706,2048); lit at t == 1706
	// some more values

	//				 1706
	//				 6826
	//				27306
	//			   109226
	//			   436906
	//			  1747626
	//			  6990506
	//			 27962026
	//			111848106
	//			447392426
	//		   1789569706
	//		   7158278826
	//		  28633115306
	//		 114532461226
	//		 458129844906
	//		1832519379626
	//		7330077518506
	//	   29320310074026
	//    117281240296106
	//	  469124961184426
	//   1876499844737706
	//   7505999378950826
	//  30023997515803306
	// 120095990063213226
	// 480383960252852906
	//1921535841011411626

	//0==mod(t+4-0,4);			lit at t == 0
	//0==mod(t+16-2,16);		lit at t == 2
	//0==mod(t+64-10,64);		lit at t == 10
	//0==mod(t+256-42,256);		lit at t == 42
	//0==mod(t+1024-170,1024);	lit at t == 170
	//0==mod(t+4096-682,4096);	lit at t == 682
	// this sequence's last digit likes to toggle between 0 and 2

	// hmm
	// This sequence is found here http://oeis.org/A176965
	// 2^(-1 + n) - 1/3 * (2 + (-2)^n)

	float n = p.x;
	float powwow = pow(2., n);
	float sgn = 1. - 2. * step(.000001, fract(n / 2.));
	float lit = pow(2., n - 1.) - 1. / 3. * (2. + sgn * powwow);
	if (mod(t + powwow - lit, powwow) < 0.00001)
		fragColor.rgb += .2;
}
