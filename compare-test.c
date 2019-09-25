/**********************************************************************
  gcc -O2 -I secp256k1/src/ -I secp256k1/ compare-test.c -lgmp
 **********************************************************************/

#include "libsecp256k1-config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _SYS_TIME_H_
#include <sys/time.h>
#endif
#ifndef _SYS_TIME_H_
#include <time.h>
#endif

//#include <unistd.h>
//#include <math.h>

#include "include/secp256k1.h"
#include "secp256k1.c"

/////////////////////////////////////////////////
// user settings

// 0: A + A -> A,  xA ready ; 0.291Mk/s; the fastest
// 1: J + A -> J,  xJ->xA   ; 0.274Mk/s; secp256k1_gej_add_ge_var()
// 2: J + A -> J,  xJ->xA   ; 0.267Mk/s; secp256k1_gej_add_ge()
// 3: 0*P0+k*G->J, xJ->xA   ; 0.091Mk/s; secp256k1_ecmult()
// 4: k*G->J, J+J->J, xJ->xA; 0.031Mk/s; secp256k1_ecmult_gen() + secp256k1_gej_add_var()
#define ALGO_CALC 0


/////////////////////////////////////////////////


void secp256k1_scalar_set_uint64(secp256k1_scalar *r, uint64_t v){
    #if defined(USE_SCALAR_4X64)
	r->d[0] = v;
	r->d[1] = r->d[2] = r->d[3] = 0;
    #elif defined(USE_SCALAR_8X32)
	r->d[0] = (uint32_t)(v);
	r->d[1] = (uint32_t)(v>>32);
	r->d[2]=r->d[3]=r->d[4]=r->d[5]=r->d[6]=r->d[7]=0;
    #endif
}


void passtime(char *s, size_t max, const struct tm *tm){
	snprintf(s, max, "%01iy %01im %01id %02i:%02i:%02is", tm->tm_year-70, tm->tm_mon, tm->tm_mday-1, tm->tm_hour, tm->tm_min, tm->tm_sec);
}


/////////////////////////////////////////////////


int main(int argc, char **argv) {

	setbuf(stdout, NULL); // cancel channel buffering (for printf)

	/////////////////////////////////////////////////

	//secp256k1_context *ctx = secp256k1_context_create(SECP256K1_CONTEXT_NONE);
	secp256k1_context *ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY);

	//ecmult vars
	secp256k1_scalar scalar_K; secp256k1_scalar_set_int(&scalar_K, 123456);
	secp256k1_scalar scalar_0; secp256k1_scalar_set_int(&scalar_0, 0);
	secp256k1_gej point_0gej; secp256k1_gej_set_infinity(&point_0gej);
	//Double multiply: R = na*A + ng*G
	//secp256k1_ecmult(&ctx->ecmult_ctx, &gej, &secp256k1_gej_const_g,	&scalar_K, &scalar_0); //1x
	//secp256k1_ecmult(&ctx->ecmult_ctx, &gej, NULL,			&scalar_0, &scalar_K); //2x
	//secp256k1_ecmult(&ctx->ecmult_ctx, &gej, &point_0gej,			&scalar_0, &scalar_K); //2x

	//for secp256k1_gej_add_var
	secp256k1_gej secp256k1_gej_const_g;
	secp256k1_gej_set_ge(&secp256k1_gej_const_g, &secp256k1_ge_const_g);

	// 1
	secp256k1_scalar scalar_1;
	secp256k1_scalar_set_int(&scalar_1, (unsigned int) 1);

	unsigned char buff_b32[32];

	/////////////////////////////////////////////////

	time_t timetotal, timelast, timenow, timepass;
	timetotal = timelast = timenow = time(NULL);
	struct tm *timetm;
	char timebuff[255];

	struct timeval utimetotal, utimelast, utimenow;
	gettimeofday(&utimetotal, NULL);
	gettimeofday(&utimelast, NULL);
	gettimeofday(&utimenow, NULL);


	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	timenow = time(NULL);
	strftime(timebuff, sizeof(timebuff), "%d %b %Y %H:%M:%S", gmtime(&timenow));
	//strftime(timebuff, sizeof(timebuff), "%X", gmtime(&timenow));
	printf("\n[DATE(utc)] %s", timebuff);
	printf("\n[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]");


	#if	ALGO_CALC == 0 
	printf("\n[algo#0]	A + A -> A,  xA ready ; the fastest! ");
	#elif	ALGO_CALC == 1 
	printf("\n[algo#1]	J + A -> J,  xJ->xA   ; secp256k1_gej_add_ge_var() ");
	#elif	ALGO_CALC == 2 
	printf("\n[algo#2]	J + A -> J,  xJ->xA   ; secp256k1_gej_add_ge() ");
	#elif	ALGO_CALC == 3 
	printf("\n[algo#3]	0*P0+k*G->J, xJ->xA   ; secp256k1_ecmult() ");
	#elif	ALGO_CALC == 4 
	printf("\n[algo#4]	k*G->J, J+J->J, xJ->xA; secp256k1_ecmult_gen() ");
	#endif

	printf("\n[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]");

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

// 2A -> A (1I, 2M, 2S)
void secp256k1_ge_double_var(secp256k1_ge *r, const secp256k1_ge *a){

	secp256k1_fe rX, rY;
	secp256k1_fe tmp1, tmp2;

	tmp1 = a->y;
	if ( (a->infinity == 1) || secp256k1_fe_normalizes_to_zero_var(&tmp1) ) { 
		secp256k1_ge_set_infinity(r);
		return;
	}

	// s = (3x^2 + A)/(2y)				# 1I 2M 1S
	// A=0 B=7 (secp256k1)
	secp256k1_fe_sqr(&tmp2, &a->x);//1S
	tmp1 = tmp2;
	secp256k1_fe_add(&tmp1, &tmp2);
	secp256k1_fe_add(&tmp1, &tmp2);
	secp256k1_fe_normalize_weak(&tmp1);

	tmp2 = a->y;
	secp256k1_fe_add(&tmp2, &a->y);
	secp256k1_fe_normalize_weak(&tmp2);
	secp256k1_fe_inv_var(&tmp2, &tmp2);//1I

	secp256k1_fe_mul(&tmp1, &tmp1, &tmp2);//1M

	// x' = s^2-2x					# 1S
	secp256k1_fe_sqr(&rX, &tmp1);//1S
	tmp2 = a->x;
	secp256k1_fe_add(&tmp2, &a->x);
	secp256k1_fe_negate(&tmp2, &tmp2, 1);
	secp256k1_fe_add(&rX, &tmp2);
	secp256k1_fe_normalize_weak(&rX);

	// y' = s(x-x')-y				# 1M
	secp256k1_fe_negate(&tmp2, &rX, 1);
	rY = a->x;
	secp256k1_fe_add(&rY, &tmp2);
	secp256k1_fe_normalize_weak(&rY);
	secp256k1_fe_mul(&rY, &rY, &tmp1);//1M
	secp256k1_fe_negate(&tmp2, &a->y, 1);
	secp256k1_fe_add(&rY, &tmp2);
	secp256k1_fe_normalize_weak(&rY);

	r->x = rX;
	r->y = rY;
}


// A + A -> A (1I, 2M, 1S)
void secp256k1_ge_add_ge_var(secp256k1_ge *r, const secp256k1_ge *a, const secp256k1_ge *b){
	

	if ( (a->infinity == 1) || (b->infinity == 1) ) { 
		if ( (a->infinity == 1) ^  (b->infinity == 1) ) {
			if (a->infinity == 1) { r->x = b->x; r->y = b->y; return; }
			if (b->infinity == 1) { r->x = a->x; r->y = a->y; return; }
		}else{
			secp256k1_ge_set_infinity(r);
			return;
		}
	}

	secp256k1_fe rX, rY;
	secp256k1_fe tmp1, tmp2;

	secp256k1_fe aXneg;
	secp256k1_fe_negate(&aXneg, &a->x, 1);
	secp256k1_fe aYneg;
	secp256k1_fe_negate(&aYneg, &a->y, 1);

	//dx = B.x - A.x
	tmp1 = b->x;
	secp256k1_fe_add(&tmp1, &aXneg);

	//dy = B.y - A.y
	tmp2 = b->y;
	secp256k1_fe_add(&tmp2, &aYneg);


	if (secp256k1_fe_normalizes_to_zero_var(&tmp1)) {
		if (secp256k1_fe_normalizes_to_zero_var(&tmp2)) {
			secp256k1_ge_double_var(r, a);
			return;
		}else{
			secp256k1_ge_set_infinity(r);
			return;
		}
	}

	secp256k1_fe_normalize_weak(&tmp1);
	secp256k1_fe_normalize_weak(&tmp2);

	secp256k1_fe_inv_var(&tmp1, &tmp1);//1I

	//c = dy * invert(dx, p) % p			# 1I,1M
	secp256k1_fe c;
	secp256k1_fe_mul(&tmp2, &tmp2, &tmp1);//1M

	//R.x = (c**2 - A.x - B.x) % p			# 1S
	secp256k1_fe_sqr(&rX, &tmp2);//1S
	secp256k1_fe_add(&rX, &aXneg);
	secp256k1_fe_negate(&tmp1, &b->x, 1);
	secp256k1_fe_add(&rX, &tmp1);
	secp256k1_fe_normalize_weak(&rX);

	//R.y = (c*(A.x - R.x) - A.y) % p		# 1M
	rY = a->x;
	secp256k1_fe_negate(&tmp1, &rX, 1);
	secp256k1_fe_add(&rY, &tmp1);
	secp256k1_fe_normalize_weak(&rY);
	secp256k1_fe_mul(&rY, &rY, &tmp2);//1M
	secp256k1_fe_add(&rY, &aYneg);
	secp256k1_fe_normalize_weak(&rY);

	r->x = rX;
	r->y = rY;
}

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

// BENCHMARKS
if(1){

	secp256k1_ge TestGe;
	secp256k1_gej TestGej;
	TestGe = secp256k1_ge_const_g;
	secp256k1_gej_set_ge(&TestGej, &TestGe);

	secp256k1_ge TestGe2;
	secp256k1_gej TestGej2;
	secp256k1_gej_add_ge_var(&TestGej2, &TestGej, &secp256k1_ge_const_g, NULL);
	secp256k1_ge_set_gej(&TestGe2, &TestGej2);


	printf("\n");
	#if defined USE_ENDOMORPHISM
	printf("\n[USE_ENDOMORPHISM]     %d", USE_ENDOMORPHISM);
	#else
	printf("\n[USE_ENDOMORPHISM]     0");
	#endif

	printf("\n[ECMULT_WINDOW_SIZE]   %d", WINDOW_G);

	#if defined ECMULT_GEN_PREC_BITS
	printf("\n[ECMULT_GEN_PREC_BITS] %d", ECMULT_GEN_PREC_BITS);
	#else
	printf("\n[ECMULT_GEN_PREC_BITS] 0");
	#endif
	printf("\n");

	uint64_t Ntimeit = 1000000; // 1M
	printf("\n[timeit] %lu keys calculating...", Ntimeit);
	printf("\n");

	gettimeofday(&utimelast, NULL);

	for(uint64_t i = 1; i<Ntimeit ; ++i){

		//// MULTIPLICATIONS
		if(1){ 

			//secp256k1_scalar_set_uint64(&scalar_K, i);		// SEQUENCE
			secp256k1_scalar_set_uint64(&scalar_K, TestGej.x.n[0]);	// RANDOM


			//// JACOBIAN


			//// Multiply: R = q*A (in constant-time)
			//// 1M at 58,0 sec of 1core i7-6820
			//secp256k1_ecmult_const(&TestGej, &secp256k1_ge_const_g, &scalar_K, 256);

#if	ALGO_CALC == 4 
			//// Multiply: R = a*G (with the generator)
			//// To harden against timing attacks .. Break up the multiplicand into groups of..
			//// 1M at 26,5 sec of 1core i7-6820
			secp256k1_ecmult_gen(&ctx->ecmult_gen_ctx, &TestGej, &scalar_K);
#endif

			//// R = na*A + ng*G
			////      k*G + 0*G
			//// 1M at 16,8 sec of 1core i7-6820
			//// secp256k1_ecmult(*ctx, gej *r, const gej *a, const scalar *na, const scalar *ng);
			//secp256k1_ecmult(&ctx->ecmult_ctx, &TestGej, &secp256k1_gej_const_g,	&scalar_K,	NULL); // NULL not recommended!
			//secp256k1_ecmult(&ctx->ecmult_ctx, &TestGej, &secp256k1_gej_const_g,	&scalar_K,	&scalar_0); //k*G + 0*G

#if	ALGO_CALC == 3 
			//// R = na*A + ng*G
			////     0*P0 + k*G
			//// 1M at 10,2 sec of 1core i7-6820
			//secp256k1_ecmult(&ctx->ecmult_ctx ,&TestGej, NULL,			&scalar_0,	&scalar_K); // NULL not recommended!
			secp256k1_ecmult(&ctx->ecmult_ctx, &TestGej, &point_0gej,		&scalar_0,	&scalar_K); // 0*P0 + k*G
#endif
		}

		//// ADDITIONS
		if(1){ 

			//// JACOBIAN

#if	ALGO_CALC == 2 
			//// J + A -> J (7M, 5S)
			//// 1M at 0,33 sec of 1core i7-6820
			secp256k1_gej_add_ge(&TestGej, &TestGej, &secp256k1_ge_const_g);
#endif

#if	ALGO_CALC == 1 
			//// J + A -> J (8M, 3S)
			//// 1M at 0,27 sec of 1core i7-6820
			secp256k1_gej_add_ge_var(&TestGej, &TestGej, &secp256k1_ge_const_g, NULL);
#endif
			//// 2J -> J (4M, 4S)
			//// 1M at 0,17 sec of 1core i7-6820
			//secp256k1_gej_double_var(&TestGej, &TestGej, NULL);

			//// J + J -> J (12M, 4S)
			//// 1M at 0,41 sec of 1core i7-6820
			//secp256k1_gej_add_var(&TestGej, &TestGej, &secp256k1_gej_const_g, NULL);



			////AFFINES
#if	ALGO_CALC == 0 
			//// A + A -> A (1I, 2M, 1S)
			//// 1M at 3,36 sec of 1core i7-6820
			secp256k1_ge_add_ge_var(&TestGe, &TestGe, &secp256k1_ge_const_g);
#endif
			//// 2A -> A (1I, 2M, 2S)
			//// 1M at 3,32 sec of 1core i7-6820
			//secp256k1_ge_double_var(&TestGe, &TestGe);

		}

		//// Jacobian -> Affine
		//secp256k1_ge_set_gej_var(&TestGe, &TestGej);
#if	ALGO_CALC == 0 
		if(0){
#else
		if(1){
#endif
			secp256k1_fe FE;
			secp256k1_fe_storage ST;
			secp256k1_fe Zinv, Zinv2;

			secp256k1_fe_inv_var(&Zinv, &TestGej.z);
			secp256k1_fe_sqr(&Zinv2, &Zinv);

			//// X
			secp256k1_fe_mul(&FE, &TestGej.x, &Zinv2);

			if(1){
			secp256k1_fe_normalize_var(&FE);
			secp256k1_fe_to_storage(&ST, &FE);
			}

			//// Y
			if(0){
				secp256k1_fe Zinv3;
				secp256k1_fe_mul(&Zinv3, &Zinv2, &Zinv);
				secp256k1_fe_mul(&FE, &TestGej.y, &Zinv3);
				if(1){
				secp256k1_fe_normalize_var(&FE);
				secp256k1_fe_to_storage(&ST, &FE);
				}
			}
		}


	}


	/////////////////////////////////////////////////

	gettimeofday(&utimenow, NULL);
	if(utimenow.tv_usec >= utimelast.tv_usec){ 
		printf("\n[passtime] %3lisec %3limsec", utimenow.tv_sec-utimelast.tv_sec, (utimenow.tv_usec-utimelast.tv_usec)/1000 );
	}else{
		printf("\n[passtime] %3lisec %3limsec", utimenow.tv_sec-utimelast.tv_sec-1, (1000000 + utimenow.tv_usec-utimelast.tv_usec)/1000 );
	}
	gettimeofday(&utimelast, NULL);

	exit(1);
}

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////



if(1){
  timenow = time(NULL);
  timepass = difftime(timenow,timetotal);
  passtime(timebuff, sizeof(timebuff), gmtime(&timepass));
  printf("\n[#############################; passtime: %s]", timebuff);
}
if(1){
  gettimeofday(&utimenow, NULL);
  if(utimenow.tv_usec >= utimelast.tv_usec){ 
	printf("\n[####################################; precision: %3lis %3lims]", utimenow.tv_sec-utimelast.tv_sec, (utimenow.tv_usec-utimelast.tv_usec)/1000 );
  }else{
	printf("\n[####################################; precision: %3lis %3lims]", utimenow.tv_sec-utimelast.tv_sec-1, (1000000 + utimenow.tv_usec-utimelast.tv_usec)/1000 );
  }
}

	printf("\n[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]");

	timenow = time(NULL);
	strftime(timebuff, sizeof(timebuff), "%d %b %Y %H:%M:%S", gmtime(&timenow));
	//strftime(timebuff, sizeof(timebuff), "%X", gmtime(&timenow));
	printf("\n[DATE(utc)] %s", timebuff);

	printf("\n[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]");

  printf("\n[x] EXIT\n");exit(EXIT_SUCCESS);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////


