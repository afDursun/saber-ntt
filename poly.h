#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "SABER_params.h"


#define MOD 20972417
#define Mprime 1544799103
#define RmodM -4378189
#define RinvN 4191613
#define R2invN 8870406


typedef struct
{
  uint16_t coeffs[SABER_N];
} poly;

typedef struct{
  poly vec[SABER_L];
} polyvec;

void InnerProd1(uint16_t pkcl[SABER_L][SABER_N],uint16_t skpv[SABER_L][SABER_N],uint16_t mod,uint16_t res[SABER_N]);
void MatrixVectorMul1(polyvec *a, uint16_t skpv[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], uint16_t mod, int16_t transpose);
void MatrixVectorMul(const uint16_t a[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);
void GenMatrix(uint16_t a[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES]);
void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES]);

#endif
