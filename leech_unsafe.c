/*
 * Author: Alex van Poppelen
 * Email: avanpoppelen@gmail.com
 *
 * This is the maximum likelihood decoder of the Leech
 * lattice of Vardy and Be'ery, "Maximum Likelihood
 * Decoding of the Leech Lattice", Trans. on Information
 * Theory, Vol.39, No. 4, July 1993.
 *
 * This is an integer implementation. It assumes the
 * space ZZ mod Q, where Q is divisible by 8. It assumes
 * a scaling of the Leech lattice such that its minimal
 * vectors have norm Q/sqrt(2). In other words, it is a
 * decoder of L/2L, where L is the Leech lattice.
 *
 * This implementation is not intended for cryptographic
 * purposes, as it is not resistant to timing attacks.
 */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "leech.h"

typedef struct Pen {
	uint64_t val;
	uint64_t ijk;
	uint8_t l;
	uint8_t c;
} Pen;

typedef struct Q24 {
	uint64_t (*d_ij)[2][2];
	int64_t (*delta_ij)[2][2];
	uint8_t (*offsets);

	uint8_t i1[6][4];
	uint8_t j1[6][4];
	uint8_t k1[6][4];
	uint8_t i2[6][4];
	uint8_t j2[6][4];
	uint8_t k2[6][4];
	uint64_t mu_x[6][4];
	Pen *pens;
	Pen *pen[3][24];

	uint64_t S_j[3][4][4];

	uint64_t M[64];
	uint64_t pts[64];

	uint64_t cv;
	uint64_t d;

	int coset_h;
	int coset_q;
} Q24;

/* Returns the minimum of two unsigned integers.
 */
static uint64_t minu_unsafe(uint64_t u1, uint64_t u2)
{
	return u1 < u2 ? u1 : u2;
}

/* Same as minu_unsafe() but sets a flag.
 */
static uint64_t minu_b_unsafe(uint64_t u1, uint64_t u2, uint8_t *flag)
{
	if (u1 < u2) {
		*flag = 0;
		return u1;
	} else {
		*flag = 1;
		return u2;
	}
}

/* Returns the absolute value of a signed integer.
 */
static uint64_t abs_unsafe(int64_t i)
{
	return i >= 0 ? i : -i;
}

/* Computes the squared Euclidean distance (SED)
 * between two points in ZZ^2_Q.
 */
static uint64_t dist_unsafe(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2)
{
	uint32_t dx = abs_unsafe((int64_t)x1 - (int64_t)x2);
	uint32_t dy = abs_unsafe((int64_t)y1 - (int64_t)y2);
	dx = minu_unsafe(dx, Q - dx);
	dy = minu_unsafe(dy, Q - dy);
	return dx * dx + dy * dy;
}

/* Sort an array of pointers to penalties, by
 * comparing the val fields. Due to the size of
 * the arrays being small (24), insertion sort is
 * used to minimize overhead.
 */
static void sort_penalties(Pen **pens, int n)
{
	int i;
	for (i = 1; i < n; ++i) {
		Pen *p = pens[i];
		int j = i;
		while (j > 0 && p->val < pens[j - 1]->val) {
			pens[j] = pens[j - 1];
			--j;
		}
		pens[j] = p;
	}
}

/* Decodes a X_ijk subset, for the 32-QAM constellation
 * (see comments above). This function will need to be
 * modified if another constellation is used.
 */
static void decode_subset(const uint32_t *t, int coset_h, int i, int j, int k,
		uint64_t *d, uint8_t *o)
{
	uint64_t d1 = dist_unsafe(t[0], t[1], SCALE * p[coset_h][i][j][k][0],
			SCALE * p[coset_h][i][j][k][1]);
	uint64_t d2 = dist_unsafe(t[0], t[1], SCALE * (p[coset_h][i][j][k][0] + 4) % Q,
			SCALE * (p[coset_h][i][j][k][1] + 4) % Q);

	*d = d1 < d2 ? d1 : d2;
	*o ^= (d1 < d2 ? 0 : 1) << (4 * i + 2 * j + k);
}

/* The precomputation step to calculate the d_ij's and
 * delta_ij's. This step is done for each Leech half-
 * lattice.
 */
static void precomputation_H24(const uint32_t *t, uint64_t d_ij[12][2][2],
		int64_t delta_ij[12][2][2], uint8_t offsets[12], int coset_h)
{
	// d_ij and delta_ij are set explicitly, but
	// the offsets need to be zeroed.
	memset(offsets, 0, sizeof(uint8_t) * 12);

	int n, i, j;
	uint64_t d_ij0, d_ij1;
	for (n = 0; n < 12; ++n) {
		for (i = 0; i < 2; ++i) {
			for (j = 0; j < 2; ++j) {
				decode_subset(t + 2 * n, coset_h, i, j, 0, &d_ij0,
						&offsets[n]);
				decode_subset(t + 2 * n, coset_h, i, j, 1, &d_ij1,
						&offsets[n]);
				d_ij[n][i][j] = minu_unsafe(d_ij0, d_ij1);
				delta_ij[n][i][j] = d_ij1 - d_ij0;
			}
		}
	}
}

/* Initialize Q24 struct.
 */
static Q24 *init_Q24(uint64_t (*d_ij)[2][2], int64_t (*delta_ij)[2][2], uint8_t *offsets, int coset_h, int coset_q)
{
	Q24 *q = calloc(1, sizeof(struct Q24));

	q->d_ij = d_ij;
	q->delta_ij = delta_ij;
	q->offsets = offsets;

	q->coset_h = coset_h;
	q->coset_q = coset_q;

	Pen *pens = calloc(72, sizeof(struct Pen));
	q->pens = pens;

	int i, l, c;
	for (i = 0; i < 3; ++i) {
		for (l = 0; l < 6; ++l) {
			for (c = 0; c < 4; ++c) {
				q->pen[i][4 * l + c] = &pens[24 * i + 4 * l + c];
				q->pen[i][4 * l + c]->l = l;
				q->pen[i][4 * l + c]->c = c;
			}
		}
	}
	return q;
}

/* Computes confidence values, preferable representations,
 * and penalties. c refers to the character of GF(4).
 *
 * The preferable representations are calculated from the
 * (i1,j1,i2,j2) interpretation of c, and its complement.
 */
static void compute_vals(Q24 *q, int l, int c, int i1, int j1, int i2, int j2)
{
	int64_t d1 = q->d_ij[2 * l][i1][j1] + q->d_ij[2 * l + 1][i2][j2];
	int64_t d2 = q->d_ij[2 * l][!i1][!j1] + q->d_ij[2 * l + 1][!i2][!j2];

	// preferable representations
	if (d2 < d1) {
		i1 = !i1;
		j1 = !j1;
		i2 = !i2;
		j2 = !j2;
	}
	q->i1[l][c] = i1;
	q->j1[l][c] = j1;
	q->k1[l][c] = q->delta_ij[2 * l][i1][j1] < 0;
	q->i2[l][c] = i2;
	q->j2[l][c] = j2;
	q->k2[l][c] = q->delta_ij[2 * l + 1][i2][j2] < 0;

	// confidence value
	q->mu_x[l][c] = minu_unsafe(d1, d2);

	// penalties
	uint8_t bit, ck1, ck2;
	ck1 = (q->delta_ij[2 * l][!i1][!j1] < 0) ^ q->k1[l][c];
	ck2 = (q->delta_ij[2 * l + 1][!i2][!j2] < 0) ^ q->k2[l][c];

	q->pen[0][4 * l + c]->val = minu_b_unsafe(abs_unsafe(q->delta_ij[2 * l][i1][j1]),
			abs_unsafe(q->delta_ij[2 * l + 1][i2][j2]), &bit);
	q->pen[0][4 * l + c]->ijk = bit ? 0x01 : 0x08;

	q->pen[1][4 * l + c]->val = q->d_ij[2 * l][!i1][!j1] + q->d_ij[2 * l + 1][!i2][!j2] -
			q->d_ij[2 * l][i1][j1] - q->d_ij[2 * l + 1][i2][j2];
	q->pen[1][4 * l + c]->ijk = 0x36 ^ (ck1 ? 0x08 : 0) ^ (ck2 ? 0x01 : 0);

	q->pen[2][4 * l + c]->val = q->pen[1][4 * l + c]->val +
			minu_b_unsafe(abs_unsafe(q->delta_ij[2 * l][!i1][!j1]),
			abs_unsafe(q->delta_ij[2 * l + 1][!i2][!j2]), &bit);
	q->pen[2][4 * l + c]->ijk = 0x36 ^ (ck1 ? 0x08 : 0) ^ (ck2 ? 0x01 : 0) ^
			(bit ? 0x01 : 0x08);

	// switch last two penalties based on intrinsic
	// change to k-parity
	if (ck1 != ck2) {
		uint64_t tmp_val = q->pen[1][4 * l + c]->val;
		uint64_t tmp_ijk = q->pen[1][4 * l + c]->ijk;
		q->pen[1][4 * l + c]->val = q->pen[2][4 * l + c]->val;
		q->pen[1][4 * l + c]->ijk = q->pen[2][4 * l + c]->ijk;
		q->pen[2][4 * l + c]->val = tmp_val;
		q->pen[2][4 * l + c]->ijk = tmp_ijk;
	}
}

/* Uses offset and Leech lattice point information
 * to extract 48 bits from the decoder, of which
 * 12 are redundant. For constellations other than
 * 32-QAM, this function will need to be modified.
 */
uint64_t output_bits(Q24 *q)
{
	uint64_t o;
	int l;

	// 36 bits of point information (12 redundant)
	o = q->cv;

	// 12 bits of offset information (0 redundant)
	for (l = 0; l < 6; ++l) {
		uint8_t ijk = q->cv >> (6 * l);
		o ^= (uint64_t)(q->offsets[2 * l] >> (ijk >> 3 & 0x7) & 1) << (36 + 2 * l);
		o ^= (uint64_t)(q->offsets[2 * l + 1] >> (ijk & 0x7) & 1) << (36 + 2 * l + 1);
	}

	return o;
}

/* Seek to the first penalty matching a hexacode-
 * word digit in a sorted array of penalties.
 */
static Pen **penalty_seek(Pen **pens, uint8_t *word)
{
	while ((*pens)->c != word[(*pens)->l]) {
		++pens;
	}
	return pens;
}

/* Decodes the Leech quarter-lattice.
 */
static void decoder_Q24(Q24 *q)
{
	// Step 1: Computing the confidence values,
	// preferable representations, and penalties
	int l;
	for (l = 0; l < 6; ++l) {
		if (q->coset_q == 0) {
			compute_vals(q, l, 0, 0, 0, 0, 0);
			compute_vals(q, l, 1, 0, 0, 1, 1);
			compute_vals(q, l, 2, 0, 1, 0, 1);
			compute_vals(q, l, 3, 0, 1, 1, 0);
		} else {
			compute_vals(q, l, 0, 1, 0, 0, 0);
			compute_vals(q, l, 1, 0, 1, 0, 0);
			compute_vals(q, l, 2, 0, 0, 1, 0);
			compute_vals(q, l, 3, 0, 0, 0, 1);
		}
	}

	// Step 2: Sorting the penalties
	for (l = 0; l < 3; ++l) {
		sort_penalties(q->pen[l], 24);
	}

	// Step 3: Computing the confidence values
	// of the blocks
	int x1, x2;
	for (l = 0; l < 3; ++l) {
		for (x1 = 0; x1 < 4; ++x1) {
			for (x2 = 0; x2 < 4; ++x2) {
				q->S_j[l][x1][x2] = q->mu_x[2 * l][x1] +
						q->mu_x[2 * l + 1][x2];
			}
		}
	}

	// Step 4: Finding the images of the hexacodewords
	// and computing their metrics
	int w;
	for (w = 0; w < 64; ++w) {
		uint8_t *word = hexacode[w];

		// calculate parities
		int h_parity = 0;
		int k_parity = 0;
		for (l = 0; l < 6; ++l) {
			h_parity ^= q->i1[l][word[l]];
			k_parity ^= q->k1[l][word[l]] ^ q->k2[l][word[l]];
		}

		// calculate base metric
		q->M[w] = q->S_j[0][word[0]][word[1]] + q->S_j[1][word[2]][word[3]] +
				q->S_j[2][word[4]][word[5]];

		// construct base point
		q->pts[w] = 0;
		for (l = 0; l < 6; ++l) {
			q->pts[w] ^= (uint64_t)q->i1[l][word[l]] << (l * 6 + 5);
			q->pts[w] ^= (uint64_t)q->j1[l][word[l]] << (l * 6 + 4);
			q->pts[w] ^= (uint64_t)q->k1[l][word[l]] << (l * 6 + 3);
			q->pts[w] ^= (uint64_t)q->i2[l][word[l]] << (l * 6 + 2);
			q->pts[w] ^= (uint64_t)q->j2[l][word[l]] << (l * 6 + 1);
			q->pts[w] ^= (uint64_t)q->k2[l][word[l]] << (l * 6);
		}

		// check parity relative to Q24 coset
		h_parity ^= q->coset_q;
		k_parity ^= q->coset_h;

		// find and apply appropriate penalties
		if (h_parity || k_parity) {
			// In the array containing pointers to the penalties,
			// 0 stands for altering k_parity, 1 for h_parity, and
			// 2 for both.
			uint8_t a, b1, b2;
			if (h_parity && !k_parity) {
				a = 1;
				b1 = 0;
				b2 = 2;
			} else if (!h_parity && k_parity) {
				a = 0;
				b1 = 1;
				b2 = 2;
			} else {
				a = 2;
				b1 = 0;
				b2 = 1;
			}

			// check two penalty case
			Pen **pos1 = penalty_seek(q->pen[b1], word);
			Pen **pos2 = penalty_seek(q->pen[b2], word);
			Pen *p1 = *pos1;
			Pen *p2 = *pos2;

			// -> account for coordinate collision
			if (p1->l == p2->l) {
				Pen *p1_next = *penalty_seek(pos1 + 1, word);
				Pen *p2_next = *penalty_seek(pos2 + 1, word);
				if (p1->val + p2_next->val < p2->val + p1_next->val) {
					p2 = p2_next;
				} else {
					p1 = p1_next;
				}
			}

			// compare to single penalty case
			Pen *pa = *penalty_seek(q->pen[a], word);
			if (pa->val < p1->val + p2->val) {
				p1 = pa;
				p2 = NULL;
			}

			// apply penalties
			q->M[w] += p1->val + (p2 == NULL ? 0 : p2->val);
			q->pts[w] ^= ((uint64_t)p1->ijk << (6 * p1->l)) ^
					(p2 == NULL ? 0 : ((uint64_t)p2->ijk << (6 * p2->l)));
		}
	}

	// Step 5: Final minimization
	q->d = UINT64_MAX;
	for (w = 0; w < 64; ++w) {
		if (q->M[w] < q->d) {
			q->d = q->M[w];
			q->cv = q->pts[w];
			q->cv = output_bits(q);
		}
	}
}

/* Decodes the Leech lattice by means of four
 * Leech quarter-lattice decoders.
 *  *t: pointer to target vector
 * *cv: pointer to outputted bits
 *  *d: pointer to outputted SEQ
 */
void decoder_L24_unsafe(const uint32_t *t, uint64_t *cv, uint64_t *d)
{
	uint64_t d_ij[12][2][2];
	int64_t delta_ij[12][2][2];
	uint8_t offsets[12];

	Q24 *q[4];
	Q24 *q_best = NULL;
	uint64_t d_best = UINT64_MAX;

	int coset_h, coset_q;
	for (coset_h = 0; coset_h < 2; ++coset_h) {
		precomputation_H24(t, d_ij, delta_ij, offsets, coset_h);

		for (coset_q = 0; coset_q < 2; ++coset_q) {
			int i = coset_h * 2 + coset_q;
			q[i] = init_Q24(d_ij, delta_ij, offsets,
					coset_h, coset_q);
			decoder_Q24(q[i]);
			if (q[i]->d < d_best) {
				d_best = q[i]->d;
				q_best = q[i];
			}
		}
	}

	*d = d_best;
	*cv = q_best->cv;

	int i;
	for (i = 0; i < 4; ++i) {
		free(q[i]->pens);
		free(q[i]);
	}
}
