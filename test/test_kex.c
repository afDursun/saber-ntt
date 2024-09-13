#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "../api.h"
#include "../poly.h"
#include "../rng.h"
#include "../SABER_indcpa.h"
#include "../verify.h"
#include <sys/time.h>
// void fprintBstr(char *S, unsigned char *A, unsigned long long L)
// {
// 	unsigned long long  i;

// 	printf("%s", S);

// 	for ( i=0; i<L; i++ )
// 		printf("%02X", A[i]);

// 	if ( L == 0 )
// 		printf("00");

// 	printf("\n");
// }

uint64_t clock1,clock2;
uint64_t clock_kp_mv,clock_cl_mv, clock_kp_sm, clock_cl_sm;

static int test_kem_cca()
{

struct timeval timeval_start, timeval_end;
  	uint8_t pk[SABER_PUBLICKEYBYTES];
  	uint8_t sk[SABER_SECRETKEYBYTES];
 	 uint8_t c[SABER_BYTES_CCA_DEC];	
  	uint8_t k_a[SABER_KEYBYTES], k_b[SABER_KEYBYTES];
	
  	unsigned char entropy_input[48];
	
  	uint64_t i, j, NTESTS;
  	NTESTS=10000;	

	gettimeofday(&timeval_start, NULL);
    for ( i = 0; i < NTESTS; i++)
    {
        crypto_kem_keypair(pk, sk);
    }
    gettimeofday(&timeval_end, NULL);
    printf("The average time of crypto_kem_keypair:\t %.3lf us \n", ((timeval_end.tv_usec + timeval_end.tv_sec * 1000000) - (timeval_start.tv_sec * 1000000 + timeval_start.tv_usec)) / (NTESTS * 1.0));


	gettimeofday(&timeval_start, NULL);
    for ( i = 0; i < NTESTS; i++)
    {
        crypto_kem_enc(c, k_a, pk);
    }
    gettimeofday(&timeval_end, NULL);
    printf("The average time of crypto_kem_enc:\t %.3lf us \n", ((timeval_end.tv_usec + timeval_end.tv_sec * 1000000) - (timeval_start.tv_sec * 1000000 + timeval_start.tv_usec)) / (NTESTS * 1.0));



	gettimeofday(&timeval_start, NULL);
    for ( i = 0; i < NTESTS; i++)
    {
     crypto_kem_dec(k_b, c, sk);
    }
    gettimeofday(&timeval_end, NULL);
    printf("The average time of crypto_crypto_kem_deckem_keypair:\t %.3lf us \n", ((timeval_end.tv_usec + timeval_end.tv_sec * 1000000) - (timeval_start.tv_sec * 1000000 + timeval_start.tv_usec)) / (NTESTS * 1.0));



	return 0;
}

/*
void test_kem_cpa(){

	uint8_t pk[SABER_PUBLICKEYBYTES];
	uint8_t sk[SABER_SECRETKEYBYTES];

	indcpa_kem_keypair(unsigned char *pk, unsigned char *sk);
	indcpa_kem_enc(unsigned char *message_received, unsigned char *noiseseed, const unsigned char *pk, unsigned char *ciphertext)
	indcpa_kem_dec(const unsigned char *sk, const unsigned char *ciphertext, unsigned char message_dec[])
}
*/
int main()
{
	test_kem_cca();
	return 0;
}
