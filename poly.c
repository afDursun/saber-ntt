#include <stdio.h>
#include "api.h"
#include "poly.h"
#include "poly_mul.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"

const int mul_table[64] = {9259893, -9259893, 7725652, -7725652, -741305, 741305, 5430861, -5430861, 8941483, -8941483, -1904165, 1904165, -1510590, 1510590, 4179776, -4179776, 1087096, -1087096, 5477510, -5477510, 4288247, -4288247, -4589585, 4589585, 4406378, -4406378, -962293, 962293, -4903106, 4903106, -3217734, 3217734, -4768008, 4768008, -10114282, 10114282, -1582389, 1582389, -9352320, 9352320, 2044676, -2044676, -1889742, 1889742, -3496807, 3496807, -5815793, 5815793, -1035446, 1035446, 8431160, -8431160, 563090, -563090, -4899008, 4899008, 9979851, -9979851, -8135419, 8135419, -3075030, 3075030, 7867132, -7867132};
const int root_table[64] = {1126591, -5474819, -547396, -5772959, 4437632, 2009459, -8766387, 5132106, -9595122, -3921010, 9259893, 7725652, -741305, 5430861, -3392653, 6496586, -9784786, 8941483, -1904165, -1510590, 4179776, 3590343, -8974179, 1102780, 1087096, 5477510, 4288247, -4589585, -245695, 1446745, -440588, 4406378, -962293, -4903106, -3217734, 9479185, 710543, 7712526, -4768008, -10114282, -1582389, -9352320, 1213433, -9241782, -4253525, 2044676, -1889742, -3496807, -5815793, -2968645, -8075200, 9346672, -1035446, 8431160, 563090, -4899008, -5531124, 8159833, 1229943, 9979851, -8135419, -3075030, 7867132};
const int inv_root_table[65] = {5531124, -1229943, -8159833, -7867132, 3075030, 8135419, -9979851, 2968645, -9346672, 8075200, 4899008, -563090, -8431160, 1035446, -1213433, 4253525, 9241782, 5815793, 3496807, 1889742, -2044676, -9479185, -7712526, -710543, 9352320, 1582389, 10114282, 4768008, 245695, 440588, -1446745, 3217734, 4903106, 962293, -4406378, -3590343, -1102780, 8974179, 4589585, -4288247, -5477510, -1087096, 3392653, 9784786, -6496586, -4179776, 1510590, 1904165, -8941483, -5132106, 3921010, 9595122, -5430861, 741305, -7725652, -9259893, -179895, 547396, 5474819, 8766387, -2009459, -4437632, 5772959, 8870406, RmodM};

void MatrixVectorMul(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose)
{
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc(A[j][i], s[j], res[i]);
			}
			else
			{
				poly_mul_acc(A[i][j], s[j], res[i]);
			}
		}
	}
}

uint16_t modinv(uint16_t x, uint16_t mod)
{
	int16_t t, newt, r, newr, quotient, temp;
	t = 0;
	newt = 1;
	r = mod;
	newr = x;

	while (newr != 0)
	{
		quotient = r / newr;
		temp = newt;
		newt = t - quotient * newt;
		t = temp;
		temp = newr;
		newr = r - quotient * newr;
		r = temp;
	}

	if (t < 0)
	{
		t = t + mod;
	}

	return t;
}

int16_t montgomery_reduce(int32_t a, int32_t M, int32_t M_inv)
{
	int16_t t;

	// Montgomery indirgeme işlemi
	t = (int16_t)((a * M_inv) & 0xFFFF); // lower 16-bit alınarak Montgomery indirgeme yapılır
	t = (a - (int32_t)t * M) >> 16;		 // Yüksek bitler alınarak tekrar mod M ile indirgeme yapılır

	return t; // 16-bit sonuç döndürülür
}

int16_t barrett_reduce(int16_t a, int16_t M)
{
	int16_t t;
	const int16_t v = ((1 << 26) + M / 2) / M;

	t = ((int32_t)v * a + (1 << 25)) >> 26;
	t *= M;
	return a - t;
}

static int16_t fqmul(int16_t a, int16_t b)
{
	return montgomery_reduce((int32_t)a * b, MOD, Mprime);
}

void my_mul1(int16_t *mul_out, int16_t *vector, const int16_t *root_table, int16_t M, int16_t M_inv,
			 int16_t *Matrix_0, int16_t *Matrix_1, int16_t *Matrix_2, int16_t counter)
{

	int32_t lower, upper, tmp1;
	int32_t K0, K1, K2, K3;
	int32_t B0, B1, B2, B3;
	int32_t root;
	int32_t c0_lower, c0_upper, c1_lower, c1_upper, c2_lower, c2_upper, c3_lower, c3_upper;

	while ((int32_t)mul_out < counter)
	{
		// small[0] * Matrix_0
		K0 = vector[0];
		K1 = vector[1];
		K2 = vector[2];
		K3 = vector[3];

		B0 = Matrix_0[0];
		B1 = Matrix_0[1];
		B2 = Matrix_0[2];
		B3 = Matrix_0[3];

		root = *root_table++;

		// c0 = K0B0 + root*(K1B3 + K2B2 + K3B1)
		lower = (int32_t)K1 * B3;
		upper = lower >> 16;
		lower = (int32_t)K2 * B2;
		upper += lower >> 16;
		lower = (int32_t)K3 * B1;
		upper += lower >> 16;
		tmp1 = (upper * M_inv) & 0xFFFF;
		upper = (upper + tmp1 * M) >> 16;
		c0_upper = (int32_t)K0 * B0 + root * upper;

		// c1 = K0B1 + K1B0 + root*(K2B3 + K3B2)
		lower = (int32_t)K2 * B3;
		upper = lower >> 16;
		lower = (int32_t)K3 * B2;
		upper += lower >> 16;
		tmp1 = (upper * M_inv) & 0xFFFF;
		upper = (upper + tmp1 * M) >> 16;
		c1_upper = (int32_t)K0 * B1 + (int32_t)K1 * B0 + root * upper;

		// c2 = K0B2 + K2B0 + K1B1 + root*(K3B3)
		lower = (int32_t)K3 * B3;
		upper = lower >> 16;
		tmp1 = (upper * M_inv) & 0xFFFF;
		upper = (upper + tmp1 * M) >> 16;
		c2_upper = (int32_t)K0 * B2 + (int32_t)K2 * B0 + (int32_t)K1 * B1 + root * upper;

		// c3 = K0B3 + K3B0 + K1B2 + K2B1
		c3_upper = (int32_t)K0 * B3 + (int32_t)K3 * B0 + (int32_t)K1 * B2 + (int32_t)K2 * B1;

		// Sonuçları mul_out'a yaz
		mul_out[0] = c0_upper;
		mul_out[1] = c1_upper;
		mul_out[2] = c2_upper;
		mul_out[3] = c3_upper;

		// Diğer matrisler ile işlemleri devam ettir
		vector += 16;
		Matrix_0 += 16;
		Matrix_1 += 16;
		Matrix_2 += 16;
		mul_out += 16;
	}
}

void my_mul(uint16_t *a, const uint16_t *mul_table, uint16_t mod, uint16_t modprime, uint16_t *matrix0, uint16_t *matrix1, uint16_t *res)
{
	uint64_t temp0, temp1, result;
	uint16_t tmp;

	for (int i = 0; i < 256; i++)
	{
		// Vektör ile matris çarpımı
		temp0 = (uint64_t)a[i] * mul_table[i];	   // Vektör ile çarpım
		temp1 = (uint64_t)matrix0[i] * matrix1[i]; // İki matris arasında çarpım

		// Modüler Montgomery işlemleri
		tmp = (temp0 * modprime) & ((1 << 16) - 1);	 // Montgomery indirgeme
		temp0 = (temp0 + (uint64_t)tmp * mod) >> 16; // Modprime ile çarpım ve indirgeme

		tmp = (temp1 * modprime) & ((1 << 16) - 1);
		temp1 = (temp1 + (uint64_t)tmp * mod) >> 16;

		// Sonuçların mod ile sınırlandırılması
		result = (temp0 + temp1) % mod;
		res[i] = (uint16_t)result; // Sonucu res dizisine kaydet
	}
}

void NTT_forward(int16_t r[256])
{
	unsigned int len, start, j, k;
	int16_t t, zeta;
	k = 1;
	for (len = 128; len >= 2; len >>= 1)
	{
		for (start = 0; start < 256; start = j + len)
		{
			zeta = root_table[k++];
			for (j = start; j < start + len; j++)
			{
				t = fqmul(zeta, r[j + len]);
				r[j + len] = r[j] - t;
				r[j] = r[j] + t;
			}
		}
	}
}
void NTT_inv(int16_t r[256], int16_t M)
{
	unsigned int start, len, j, k;
	int16_t t, zeta;
	const int16_t f = 1441; // mont^2/128

	k = 127;
	for (len = 2; len <= 128; len <<= 1)
	{
		for (start = 0; start < 256; start = j + len)
		{
			zeta = inv_root_table[k--];
			for (j = start; j < start + len; j++)
			{
				t = r[j];
				r[j] = barrett_reduce(t + r[j + len], M);
				r[j + len] = r[j + len] - t;
				r[j + len] = fqmul(zeta, r[j + len]);
			}
		}
	}

	for (j = 0; j < 256; j++)
		r[j] = fqmul(r[j], f);
}

void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta)
{
	r[0] = fqmul(a[1], b[1]);
	r[0] = fqmul(r[0], zeta);
	r[0] += fqmul(a[0], b[0]);
	r[1] = fqmul(a[0], b[1]);
	r[1] += fqmul(a[1], b[0]);
}

void poly_basemul_montgomery(uint16_t *r, const uint16_t *a, const uint16_t *b, const int16_t *root)
{
	unsigned int i;
	for (i = 0; i < SABER_N / 4; i++)
	{
		basemul(&r[4 * i], &a[4 * i], &b[4 * i], root[64 + i]);
		basemul(&r[4 * i + 2], &a[4 * i + 2], &b[4 * i + 2], -root[64 + i]);
	}
}

void invntt(int16_t r[256], int16_t *inv_root_table, u_int16_t M)
{
	unsigned int start, len, j, k;
	int16_t t, zeta;
	const int16_t f = 1441; // mont^2/128

	k = 127;
	for (len = 2; len <= 128; len <<= 1)
	{
		for (start = 0; start < 256; start = j + len)
		{
			zeta = inv_root_table[k++];
			for (j = start; j < start + len; j++)
			{
				t = r[j];
				r[j] = barrett_reduce(t + r[j + len], M);
				r[j + len] = r[j + len] - t;
				r[j + len] = fqmul(zeta, r[j + len]);
			}
		}
	}

	for (j = 0; j < 256; j++)
		r[j] = fqmul(r[j], f);
}

void MatrixVectorMul1(polyvec *a, uint16_t skpv[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], uint16_t mod, int16_t transpose)
{
	int matrix_temp[SABER_L][SABER_L][SABER_N];
	int mul_out[SABER_L][SABER_N];

	for (int i = 0; i < SABER_L; i++)
	{
		for (int j = 0; j < SABER_L; j++)
		{
			NTT_forward(a[i].vec[j].coeffs);
		}
	}

	for (int i = 0; i < SABER_L; i++)
	{
		NTT_forward(skpv[i]);
	}

	if (transpose == 1)
	{
		for (int j = 0; j < SABER_L; j++)
		{
			poly_basemul_montgomery(res[j], a[j].vec, skpv, root_table);
		}
	}
	else
	{
		for (int j = 0; j < SABER_L; j++)
		{
			poly_basemul_montgomery(res[j], a[j].vec, skpv, root_table);
		}
	}
	for (int j = 0; j < SABER_L; j++)
	{
		invntt(res[j], inv_root_table, MOD);
	}
}

void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N])
{
	int j;
	for (j = 0; j < SABER_L; j++)
	{
		poly_mul_acc(b[j], s[j], res);
	}
}

int32_t central_reduce(int32_t target, int32_t Mhalf, int32_t M)
{
	if (target >= Mhalf)
	{
		target -= M;
	}
	if (target < -Mhalf)
	{
		target += M;
	}
	return target;
}

// Add/Sub Fonksiyonu
static inline void add_sub(int16_t *a0, int16_t *b0,
						   int16_t *a1, int16_t *b1,
						   int16_t *a2, int16_t *b2,
						   int16_t *a3, int16_t *b3)
{
	// Toplama
	int16_t temp_a0 = *a0 + *b0;
	int16_t temp_a1 = *a1 + *b1;
	int16_t temp_a2 = *a2 + *b2;
	int16_t temp_a3 = *a3 + *b3;

	// Çıkarma (b = a - 2 * b)
	*b0 = temp_a0 - (int16_t)(*b0 << 1);
	*b1 = temp_a1 - (int16_t)(*b1 << 1);
	*b2 = temp_a2 - (int16_t)(*b2 << 1);
	*b3 = temp_a3 - (int16_t)(*b3 << 1);
}
// Montgomery Multiply Fonksiyonu
static inline int16_t montgomery_mul(int16_t a, int16_t b, int16_t M, int16_t M_inv)
{
	int64_t product = (int64_t)a * b; // 32-bit çarpım sonucu 64-bit
	int32_t lower = (int32_t)(product & 0xFFFFFFFF);
	int32_t upper = (int32_t)(product >> 32);

	int32_t tmp = (int32_t)(lower * M_inv);
	int64_t result = (int64_t)upper + (int64_t)tmp * M;

	return (int16_t)(result & 0xFFFFFFFF);
}

void NTT_inv_inner(int16_t *src, const int16_t *root, int16_t mod, int16_t mod_inv, int16_t *dest)
{
	int16_t temp[8];
	int32_t lower, upper, tmp;
	int32_t r4, r5, r6, r7, r8, r9, r10, r11;
	int32_t mod_half = mod / 2; // Mhalf
	int32_t i, j;

	// İlk dış döngü
	for (i = 0; i < 1024; i += 16)
	{
		// Vektör elemanlarını yükleme
		r4 = src[i];
		r5 = src[i + 16];
		r6 = src[i + 32];
		r7 = src[i + 48];
		r8 = src[i + 64];
		r9 = src[i + 80];
		r10 = src[i + 96];
		r11 = src[i + 112];

		// Add_sub fonksiyonu ile toplama çıkarma
		add_sub(&r4, &r5, &r6, &r7, &r8, &r9, &r10, &r11);

		// Montgomery çarpımları
		for (j = 0; j < 4; j++)
		{
			r5 = montgomery_mul(r5, root[j], mod, mod_inv);
			r7 = montgomery_mul(r7, root[j + 1], mod, mod_inv);
			r9 = montgomery_mul(r9, root[j + 2], mod, mod_inv);
			r11 = montgomery_mul(r11, root[j + 3], mod, mod_inv);
		}

		// Central reduce işlemleri
		r4 = central_reduce(r4, mod_half, mod);
		r5 = central_reduce(r5, mod_half, mod);
		r6 = central_reduce(r6, mod_half, mod);
		r7 = central_reduce(r7, mod_half, mod);
		r8 = central_reduce(r8, mod_half, mod);
		r9 = central_reduce(r9, mod_half, mod);
		r10 = central_reduce(r10, mod_half, mod);
		r11 = central_reduce(r11, mod_half, mod);

		// Sonuçları kaydetme
		dest[i] = r4 & 0x1FFF; // Sonuçları 13 bit ile sınırlandırma
		dest[i + 16] = r5 & 0x1FFF;
		dest[i + 32] = r6 & 0x1FFF;
		dest[i + 48] = r7 & 0x1FFF;
		dest[i + 64] = r8 & 0x1FFF;
		dest[i + 80] = r9 & 0x1FFF;
		dest[i + 96] = r10 & 0x1FFF;
		dest[i + 112] = r11 & 0x1FFF;
	}
}

void NTT_forward_pk(int16_t *src, const int16_t *root, int16_t M, int16_t M_inv, int16_t *dest)
{
	int32_t lower, upper, tmp;
	int32_t K0, K1, K2, K3;
	int32_t B0, B1, B2, B3;
	int32_t root_value;
	int16_t *vector = src;
	int16_t *matrix_0 = src + 1024;
	int16_t *matrix_1 = src + 2048;
	int16_t *matrix_2 = src + 3072;

	while (vector < src + 1536)
	{ // Dış döngü

		// İlk set K0-K3 ve B0-B3 değerlerini yükle
		K0 = vector[0];
		K1 = vector[1];
		K2 = vector[2];
		K3 = vector[3];

		B0 = matrix_0[0];
		B1 = matrix_0[1];
		B2 = matrix_0[2];
		B3 = matrix_0[3];

		root_value = *root++;

		// c0 = K0B0 + root * (K1B3 + K2B2 + K3B1)
		lower = (int64_t)K1 * B3;
		upper = lower >> 16;
		lower = (int64_t)K2 * B2;
		upper += lower >> 16;
		lower = (int64_t)K3 * B1;
		upper += lower >> 16;
		tmp = (upper * M_inv) & 0xFFFF;
		upper = (upper + tmp * M) >> 16;
		upper = (int64_t)K0 * B0 + root_value * upper;
		dest[0] = upper & 0x1FFF; // Sonucu mod ile sınırlıyoruz (13 bit)

		// c1 = K0B1 + K1B0 + root * (K2B3 + K3B2)
		lower = (int64_t)K2 * B3;
		upper = lower >> 16;
		lower = (int64_t)K3 * B2;
		upper += lower >> 16;
		tmp = (upper * M_inv) & 0xFFFF;
		upper = (upper + tmp * M) >> 16;
		upper = (int64_t)K0 * B1 + (int64_t)K1 * B0 + root_value * upper;
		dest[1] = upper & 0x1FFF;

		// c2 = K0B2 + K2B0 + K1B1 + root * (K3B3)
		lower = (int64_t)K3 * B3;
		upper = lower >> 16;
		tmp = (upper * M_inv) & 0xFFFF;
		upper = (upper + tmp * M) >> 16;
		upper = (int64_t)K0 * B2 + (int64_t)K2 * B0 + (int64_t)K1 * B1 + root_value * upper;
		dest[2] = upper & 0x1FFF;

		// c3 = K0B3 + K3B0 + K1B2 + K2B1
		upper = (int64_t)K0 * B3 + (int64_t)K3 * B0 + (int64_t)K1 * B2 + (int64_t)K2 * B1;
		dest[3] = upper & 0x1FFF;

		// Sonraki elemanlara geç
		vector += 4;
		matrix_0 += 4;
		dest += 4;
	}
}



void InnerProd1(uint16_t pkcl[SABER_L][SABER_N], uint16_t skpv[SABER_L][SABER_N], uint16_t mod, uint16_t res[SABER_N])
{

	int pkcl_temp[SABER_L][SABER_N];
	int skpv_temp[SABER_L][SABER_N];
	int res_temp[SABER_N];

	for (int i = 0; i < SABER_L; i++)
	{
		NTT_forward(&pkcl[i][0]);
	}
	for (int i = 0; i < SABER_L; i++)
	{
		NTT_forward(&skpv[i][0]);
	}

	poly_basemul_montgomery(res_temp, pkcl[0], skpv[0], root_table);
	invntt(res_temp, inv_root_table, MOD);
	// NTT_inv_inner(&res_temp[0], inv_root_table, MOD, Mprime, &res[0]);
}

void GenMatrix(uint16_t A[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYVECBYTES];
	int i;

	shake128(buf, sizeof(buf), seed, SABER_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		BS2POLVECq(buf + i * SABER_POLYVECBYTES, A[i]);
	}
}

void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYCOINBYTES];
	size_t i;

	shake128(buf, sizeof(buf), seed, SABER_NOISE_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		cbd(s[i], buf + i * SABER_POLYCOINBYTES);
	}
}
