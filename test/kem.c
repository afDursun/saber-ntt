#include "../api.h"
#include "../poly.h"
#include "../rng.h"
#include "../SABER_indcpa.h"
#include "../verify.h"
#include "cpucycles.c"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define NTESTS 10000
static int test_kem_cca()
{


  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t ct[CRYPTO_CIPHERTEXTBYTES];	
  uint8_t ss_a[CRYPTO_BYTES], ss_b[CRYPTO_BYTES];
	
  unsigned char entropy_input[48];
	
  uint64_t i;


   
    	for (i=0; i<48; i++)
        	entropy_input[i] = i;
    	randombytes_init(entropy_input, NULL, 256);

struct timeval timeval_start, timeval_end;

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
     crypto_kem_enc(ct, ss_a, pk);
    }
    gettimeofday(&timeval_end, NULL);
    printf("The average time of crypto_kem_enc:\t %.3lf us \n", ((timeval_end.tv_usec + timeval_end.tv_sec * 1000000) - (timeval_start.tv_sec * 1000000 + timeval_start.tv_usec)) / (NTESTS * 1.0));


    gettimeofday(&timeval_start, NULL);
    for ( i = 0; i < NTESTS; i++)
    {
         crypto_kem_dec(ss_b, ct, sk);
	  
    }
    gettimeofday(&timeval_end, NULL);
    printf("The average time of crypto_kem_dec:\t %.3lf us \n", ((timeval_end.tv_usec + timeval_end.tv_sec * 1000000) - (timeval_start.tv_sec * 1000000 + timeval_start.tv_usec)) / (NTESTS * 1.0));


	    

	   

	  

	    // Functional verification: check if ss_a == ss_b?
	    for(i=0; i<SABER_KEYBYTES; i++)
	    {
		//printf("%u \t %u\n", ss_a[i], ss_b[i]);
		if(ss_a[i] != ss_b[i])
		{
			printf(" ----- ERR CCA KEM ------\n");		
			break;
		}
		else{
			printf(" ----- BAŞARILı ------\n");	
			break;
		}
	    }


  	return 0;
}



int main()
{

	test_kem_cca();
	return 0;
}
