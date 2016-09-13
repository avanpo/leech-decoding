/* 
 * Author: Alex van Poppelen
 * Email: avanpoppelen@gmail.com
 *
 * A collection of utility functions for debugging
 * Leech lattice vectors, and partial information
 * collected during the different decoding functions.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "leech.h"

/* Leech lattice generator matrix in standard MOG
 * coordinates. Taken from Conway, JH and Sloane, NJA.
 * "Sphere Packings, Lattices and Groups". Springer
 * (1993).
 *
 * Note that the decoder by Vardy and Be'ery uses a
 * different set of coordinates.
 */
uint8_t leech[24][24] = {
	{ 8,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,4,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,0,4,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,0,0,4, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},

	{ 4,0,0,0, 4,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,0,0,0, 0,4,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,0,0,0, 0,0,4,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 2,2,2,2, 2,2,2,2, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},

	{ 4,0,0,0, 0,0,0,0, 4,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,0,0,0, 0,0,0,0, 0,4,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 4,0,0,0, 0,0,0,0, 0,0,4,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 2,2,2,2, 0,0,0,0, 2,2,2,2, 0,0,0,0, 0,0,0,0, 0,0,0,0},

	{ 4,0,0,0, 0,0,0,0, 0,0,0,0, 4,0,0,0, 0,0,0,0, 0,0,0,0},
	{ 2,2,0,0, 2,2,0,0, 2,2,0,0, 2,2,0,0, 0,0,0,0, 0,0,0,0},
	{ 2,0,2,0, 2,0,2,0, 2,0,2,0, 2,0,2,0, 0,0,0,0, 0,0,0,0},
	{ 2,0,0,2, 2,0,0,2, 2,0,0,2, 2,0,0,2, 0,0,0,0, 0,0,0,0},

	{ 4,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 4,0,0,0, 0,0,0,0},
	{ 2,0,2,0, 2,0,0,2, 2,2,0,0, 0,0,0,0, 2,2,0,0, 0,0,0,0},
	{ 2,0,0,2, 2,0,2,0, 2,2,0,0, 2,0,2,0, 2,0,2,0, 0,0,0,0},
	{ 2,2,0,0, 2,0,2,0, 2,0,0,2, 0,0,0,0, 2,0,0,2, 0,0,0,0},

	{ 0,2,2,2, 2,0,0,0, 2,0,0,0, 2,0,0,0, 2,0,0,0, 2,0,0,0},
	{ 0,0,0,0, 0,0,0,0, 2,2,0,0, 2,2,0,0, 2,2,0,0, 2,2,0,0},
	{ 0,0,0,0, 0,0,0,0, 2,0,2,0, 2,0,2,0, 2,0,2,0, 2,0,2,0},
	{-3,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}};

/*******************************************************
 * Computation
 *******************************************************/

/* Copies a point in 24-dimensional space from src
 * to destination.
 */
void copyv(uint32_t *dest, const uint32_t *src)
{
	memcpy(dest, src, sizeof(uint32_t) * 24);
}

/* Compare two points in 24-dimensional space for
 * equality. Returns 0 if equal, 1 otherwise.
 */
int cmpv(const uint32_t *a, const uint32_t *b)
{
	int i;
	for (i = 0; i < 24; ++i) {
		if (a[i] != b[i])
			return 1;
	}
	return 0;
}

/* Return the minimum of two unsigned integers.
 */
uint64_t minu(uint64_t u1, uint64_t u2)
{
	return u1 < u2 ? u1 : u2;
}

/* Return the absolute value of a signed integer.
 */
uint64_t absi(int64_t i)
{
	return i >= 0 ? i : -i;
}

/* Calculates squared Euclidean distance between
 * two points in 24-dimensional space.
 */
uint64_t sed(const uint32_t *a, const uint32_t *b)
{
	int i;
	uint64_t d = 0;
	for (i = 0; i < 24; ++i) {
		uint32_t di = absi((int32_t)a[i] - (int32_t)b[i]);
		di = minu(di, Q - di);
		d += di * di;
	}
	return d;
}

/******************************************************
 * Hexacode
 ******************************************************/

uint8_t projection(uint8_t i1, uint8_t j1, uint8_t i2, uint8_t j2)
{
	if ((i1 && j1 && i2 && j2) || (!i1 && !j1 && !i2 && !j2) ||
		(!i1 && j1 && i2 && j2) || (i1 && !j1 && !i2 && !j2)) {
		return 0;
	} else if ((i1 && j1 && !i2 && !j2) || (!i1 && !j1 && i2 && j2) ||
		(!i1 && j1 && !i2 && !j2) || (i1 && !j1 && i2 && j2)) {
		return 1;
	} else if ((i1 && !j1 && i2 && !j2) || (!i1 && j1 && !i2 && j2) ||
		(!i1 && !j1 && i2 && !j2) || (i1 && j1 && !i2 && j2)) {
		return 2;
	} else {
		return 3;
	}
}

/******************************************************
 * Encodings
 ******************************************************/

/* Decodes the output from a decoder to a point in space.
 */
void decode_pt(uint32_t *out, uint64_t cv, int verbose)
{
	// determine if A-type or if B-type by counting
	// k-parity
	int l, coset_h = 0;
	for (l = 0; l < 12; ++l) {
		coset_h ^= (cv >> 3 * l) & 1;
	}

	int n1 = 0, n2 = 0;
	char line1[256], line2[256];

	for (l = 0; l < 6; ++l) {
		uint8_t col = cv >> (6 * l);
		uint8_t o1 = cv >> (36 + 2 * l);
		uint8_t o2 = cv >> (36 + 2 * l + 1);
		out[4 * l] = p[coset_h][col >> 5 & 1][col >> 4 & 1][col >> 3 & 1][0] + 4 * o1;
		out[4 * l + 1] = p[coset_h][col >> 5 & 1][col >> 4 & 1][col >> 3 & 1][1] + 4 * o1;
		out[4 * l + 2] = p[coset_h][col >> 2 & 1][col >> 1 & 1][col & 1][0] + 4 * o2;
		out[4 * l + 3] = p[coset_h][col >> 2 & 1][col >> 1 & 1][col & 1][1] + 4 * o2;
		if (verbose) {
			n1 += snprintf(line1 + n1, 256 - n1, "%c_%1d%1d%1d   ", coset_h ? 'B' : 'A', col >> 5 & 1, col >> 4 & 1, col >> 3 & 1);
			n2 += snprintf(line2 + n2, 256 - n2, "%c_%1d%1d%1d   ", coset_h ? 'B' : 'A', col >> 2 & 1, col >> 1 & 1, col & 1);
		}
	}

	if (verbose) {
		printf("%s\n", line1);
		printf("%s\n", line2);
	}

	for (l = 0; l < 24; ++l) {
		out[l] = SCALE * out[l] % Q;
	}
}


/******************************************************
 * Printing
 ******************************************************/

/* Print a 24-dimensional vector.
 */
void printvu(const uint32_t *v)
{
	int i;
	for (i = 0; i < 24; ++i) {
		printf("%4u ", v[i]);
		if (i == 11 || i == 23) printf("\n");
	}
}

/* Print a binary representation of a Leech lattice
 * vector. Includes 12 offsets bits + 36 point bits
 * (12 redundant).
 */
void printvb(const uint64_t v, uint8_t spacing)
{
	int i;
	for (i = 47; i >= 0; --i) {
		printf("%1lu", v >> i & 1);
		if (i % spacing == 0) printf(" ");
	}
	printf("\n");
}
