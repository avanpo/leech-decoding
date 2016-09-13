/*
 * Author: Alex van Poppelen
 * Email: avanpoppelen@gmail.com
 *
 * This is the maximum likelihood decoder of the Leech
 * lattice of Vardy and Be'ery, "Maximum Likelihood
 * Decoding of the Leech Lattice", Trans. on Information
 * Theory, Vol. 39, No. 4, July 1993.
 *
 ************************************************************
 * I have attempted to make this code timing and cache
 * attack resistant. Therefore there is some deviation from
 * the paper described above. Instead of sorting the
 * penalties once and keeping an ordered list, the
 * relevant penalties will need to be (partially) sorted for
 * each hexacodeword. This is possibly slower. The paper
 * assumes a priori information for each hexacodeword after
 * sorting, which is not realistic in an implementation, as
 * the first applicable penalty may be near the end of the
 * sorted list, and will require significant traversing of
 * said list.
 *
 * Notes (and assumptions):
 *
 * - The constant factor of of 1/sqrt(8) typically associated
 *   with the Leech lattice has been ignored. Additionally,
 *   the order of the coordinates of the construction used in
 *   the paper is not entirely consistent with the generator
 *   matrix typically associated with the Leech lattice. The
 *   resulting lattice is therefore some rotation and/or
 *   reflection of the familiar Leech lattice basis.
 *
 * - The decoder only recognizes 2^24 different "types" of
 *   Leech lattice points. Points offset by a minimal vector
 *   of shape (4,4,0^22) mod Q where the non-zero coordinates
 *   are consecutive and begin in an even coordinate, are
 *   considered to be of the same type. As we are considering
 *   2^36 points, this means there are 2^12 points for every
 *   "type" (in the context of the paper, this is the
 *   difference between using a 16-QAM or 32-QAM
 *   constellation). This gives 12 "offset" bits.
 *
 * - The points are represented in memory using 36 point bits
 *   (of which 12 are redundant) and 12 offset bits. The 12
 *   offset bits are stored first, bitwise little endian. The
 *   point bits are stored next, in groups of six
 *   (i1,j1,k1,i2,j2,k2), where the groups are ordered little
 *   endian.
 *
 *   Example:
 *     000000 000001 110111 000000 111110 000000 000000 000000
 *     +-----------+ +---------------------------------------+
 *        offsets                   point bits
 *
 *   Refers to:
 *     ^A_000   A_000   A_000   A_111   A_000   A_110
 *      A_000   A_000   A_000   A_110   A_000   A_111
 *   where the hat ^ refers to an offset.
 *
 *   In vector notation this is:
 *     (4,4,0,0,0,0,0,0,0,0,0,0,2,6,2,2,0,0,0,0,2,2,2,6).
 *
 *   The points are therefore stored in a uint64_t.
 *
 * - I have assumed arithmetic instructions on 64 bit
 *   integers are constant time. This includes multiplication,
 *   addition, and bit shifts of arbitrary size.
 *
 * - I have assumed a right bit shift of size 63 on a 64 bit
 *   integer casted to an unsigned 64 bit integer will reveal
 *   the sign bit, and that this is portable.
 *
 * - This code is set up for a 32-QAM constellation. This means
 *   it decodes exactly 2^36 Leech lattice points. It can be
 *   reworked to decode 2^(12*k) points for k >= 2, but this
 *   will require different machinery in the precomputation
 *   phase (specifically the decode_subset(...) function).
 *   It will also require the output to include the extra bits
 *   that are extracted. Note that the offset bits do not need
 *   to be carried through the Q_24 decoder. Instead, these
 *   bits are calculated during the precomputation phase and
 *   can be applied to the ouputted point using (assumed
 *   constant) bit shifts.
 *
 * - There is some potential for optimization. Right now, a
 *   penalty stored a 64 bit SED value, and 64 bits of point
 *   information (how it will alter the point). This last part
 *   could easily be packed down to 10 bits. If we are given
 *   a bound on Q, then we can halve the size of a penalty in
 *   memory. This would give significant savings as the
 *   minimization of penalties is a large part of the algorithm,
 *   and requires time complexity linear in this size.
 *
 * - This optimization is ripe for parallelization or
 *   vectorization. The timing and cache safe nature means this
 *   can be done without altering the algorithm at all.
 */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "leech.h"

/* This structure holds a penalty.
 *
 *   val: Change in SEQ when applied
 *   ijk: Change in point representation when applied
 *        (applied by XOR). Also contains the hexacode
 *        coordinate l in the 4 most significant bits.
 */
typedef struct Pen {
	uint64_t val;
	uint64_t ijk;
} Pen;

/* This structure holds the data for a single
 * Leech quarter-lattice coset.
 *
 * The first three fields are pointers to the data
 * calculated per Leech half-lattice coset, and is
 * constant. The fields all follow the naming
 * conventions used in the paper.
 *
 * This structure can be improved in terms of
 * memory footprint, focus was on speed.
 */
typedef struct Q24 {
	uint64_t (*d_ij)[2][2];
	int64_t (*delta_ij)[2][2];
	uint8_t *offsets;

	uint8_t i1[6][4];
	uint8_t j1[6][4];
	uint8_t k1[6][4];
	uint8_t i2[6][4];
	uint8_t j2[6][4];
	uint8_t k2[6][4];
	uint8_t o12[6][4];
	uint64_t mu_x[6][4];
	Pen pen[6][4][3];

	uint64_t S_j[3][4][4];

	uint64_t M[64];
	uint64_t pts[64];

	int coset_h;
	int coset_q;
} Q24;

/* The representations of the 16 subsets of points
 * in D_2, in a 16-QAM constellation, are given.
 *
 * A-type, even    | A-type, odd
 * ----------------+---------------
 * A_000 = (0,0)   | A_010 = (2,0)
 * A_001 = (4,0)   | A_011 = (6,0)
 * A_110 = (2,2)   | A_100 = (4,-2)
 * A_111 = (2,-2)  | A_101 = (4,2)
 *                 |
 * B-type, even    | B-type, odd
 * ----------------+---------------
 * B_000 = (1,1)   | B_010 = (3,1)
 * B_001 = (5,1)   | B_011 = (3,-3)
 * B_110 = (3,3)   | B_100 = (5,-1)
 * B_111 = (3,-1)  | B_101 = (1,-1)
 *
 * This constellation can be scaled by powers of
 * two, and it can be rotated 45 degrees (imposing
 * an additional scaling factor of sqrt(2)). This
 * preserves the integer"ness" of the Leech lattice.
 *
 * Assuming a modulus environment, only positive
 * coordinates are used. Since the decoder is
 * integer, the given representations must be
 * multiplied by SCALE.
 *
 * The representation below is equivalent to a
 * 32-QAM constellation for SCALE = Q / 8. Only
 * Only 16 D_2 points are given, the others can be
 * generated by adding SCALE*(4,4) modulo Q. It
 * represents exactly 2^36 Leech lattice points.
 *
 * To use a (4^k * 32)-QAM constellation for k >= 0,
 * or alternatively to increase the number of Leech
 * lattice points by k * 2^24, simply let
 * SCALE = Q / (8 * 2^k).
 */
uint8_t p[2][2][2][2][2] = {{{
		{{0, 0}, {4, 0}},
		{{2, 0}, {6, 0}}
		}, {
		{{4, 6}, {4, 2}},
		{{2, 2}, {2, 6}}
		}}, {{
		{{1, 1}, {5, 1}},
		{{3, 1}, {3, 5}}
		}, {
		{{5, 7}, {1, 7}},
		{{3, 3}, {3, 7}}}}};

/* The representation below is equivalent to a
 * 16-QAM constellation in an integer modulo r=4
 * environment. It represents exactly 2^24 Leech
 * lattice points.
 *
 * Note that using this representation will require
 * some refactoring in the decode_subset(...)
 * function to account for the additional points in
 * the subset, along with the appropriate value in
 * R.
uint8_t p[2][2][2][2][2] = {{{
		{{0, 0}, {2, 2}},
		{{1, 1}, {3, 3}}
		}, {
		{{3, 1}, {1, 3}},
		{{0, 2}, {2, 0}}
		}}, {{
		{{0, 1}, {2, 3}},
		{{1, 2}, {3, 0}}
		}, {
		{{3, 2}, {1, 0}},
		{{0, 3}, {2, 1}}}}};*/

/* Enumeration of the hexacode over GF(4). The four
 * characters of GF(4) are represented as {0,1,2,3}.
 */
uint8_t hexacode[64][6] = {
	{0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 1, 1},
	{0, 0, 2, 2, 2, 2}, {0, 0, 3, 3, 3, 3},
	{1, 1, 0, 0, 1, 1}, {2, 2, 0, 0, 2, 2},
	{3, 3, 0, 0, 3, 3}, {1, 1, 1, 1, 0, 0},
	{2, 2, 2, 2, 0, 0}, {3, 3, 3, 3, 0, 0},
	{1, 1, 2, 2, 3, 3}, {1, 1, 3, 3, 2, 2},
	{2, 2, 1, 1, 3, 3}, {2, 2, 3, 3, 1, 1},
	{3, 3, 1, 1, 2, 2}, {3, 3, 2, 2, 1, 1},
	{2, 3, 2, 3, 2, 3}, {3, 1, 3, 1, 3, 1},
	{1, 2, 1, 2, 1, 2}, {2, 3, 3, 2, 3, 2},
	{3, 2, 2, 3, 3, 2}, {3, 2, 3, 2, 2, 3},
	{3, 1, 1, 3, 1, 3}, {1, 3, 3, 1, 1, 3},
	{1, 3, 1, 3, 3, 1}, {1, 2, 2, 1, 2, 1},
	{2, 1, 1, 2, 2, 1}, {2, 1, 2, 1, 1, 2}, 
	{0, 1, 0, 1, 2, 3}, {0, 1, 1, 0, 3, 2},
	{1, 0, 0, 1, 3, 2}, {1, 0, 1, 0, 2, 3},
	{0, 1, 2, 3, 0, 1}, {0, 1, 3, 2, 1, 0},
	{1, 0, 2, 3, 1, 0}, {1, 0, 3, 2, 0, 1},
	{2, 3, 0, 1, 0, 1}, {2, 3, 1, 0, 1, 0},
	{3, 2, 0, 1, 1, 0}, {3, 2, 1, 0, 0, 1},
	{0, 2, 0, 2, 3, 1}, {0, 2, 2, 0, 1, 3},
	{2, 0, 0, 2, 1, 3}, {2, 0, 2, 0, 3, 1},
	{0, 2, 3, 1, 0, 2}, {0, 2, 1, 3, 2, 0},
	{2, 0, 3, 1, 2, 0}, {2, 0, 1, 3, 0, 2},
	{3, 1, 0, 2, 0, 2}, {3, 1, 2, 0, 2, 0},
	{1, 3, 0, 2, 2, 0}, {1, 3, 2, 0, 0, 2},
	{0, 3, 0, 3, 1, 2}, {0, 3, 3, 0, 2, 1},
	{3, 0, 0, 3, 2, 1}, {3, 0, 3, 0, 1, 2},
	{0, 3, 1, 2, 0, 3}, {0, 3, 2, 1, 3, 0},
	{3, 0, 1, 2, 3, 0}, {3, 0, 2, 1, 0, 3},
	{1, 2, 0, 3, 0, 3}, {1, 2, 3, 0, 3, 0},
	{2, 1, 0, 3, 3, 0}, {2, 1, 3, 0, 0, 3}};

/* Returns the minimum value among two unsigned
 * integers without branching.
 */
static uint64_t minu_safe(uint64_t d1, uint64_t d2)
{
	int sign = (uint64_t)(d1 - d2) >> 63;
	return d1 * sign + d2 * (sign ^ 1);
}

/* Returns the absolute value of a SED without
 * branching.
 */
static uint64_t abs_safe(int64_t d)
{
	int sign = (uint64_t)d >> 63;
	return d * (sign ^ 1) - d * sign;
}

/* Computes squared euclidean distance (SED) between
 * points in ZZ^2_q.
 */
static uint64_t dist_safe(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2)
{
	uint32_t dx = abs_safe((int32_t)x1 - (int32_t)x2);
	uint32_t dy = abs_safe((int32_t)y1 - (int32_t)y2);
	dx = minu_safe(dx, Q - dx);
	dy = minu_safe(dy, Q - dy);
	return dx * dx + dy * dy;
}

/* Compares two penalties and places the best in
 * first, the other in second.
 *
 * Timing and cache access are independent of input.
 */
static void cmp_place_penalties(Pen *first, Pen *second, Pen *ap, Pen *bp)
{
	int64_t a = ap->val;
	int64_t b = bp->val;
	uint8_t sign = (uint64_t)(a - b) >> 63;

	first->val = a * sign + b * (sign ^ 1);
	second->val = b * sign + a * (sign ^ 1);
	first->ijk = ap->ijk * sign + bp->ijk * (sign ^ 1);
	second->ijk = bp->ijk * sign + ap->ijk * (sign ^ 1);
}

/* Compares the best, second best, and another third
 * penalty. Assumes best and second best are already
 * in order. The third penalty is compared and possibly
 * replaces the second penalty. Then the first and now
 * second penalties are compared and possibly swapped.
 *
 * Timing and cache access are independent of input.
 */
static void cmp_swp_penalties(Pen *first, Pen *second, Pen *ap)
{
	int64_t snd = second->val;
	int64_t a = ap->val;
	uint8_t sign = (uint64_t)(snd - a) >> 63;

	// move third penalty into second
	second->val = snd * sign + a * (sign ^ 1);
	second->ijk = second->ijk * sign + ap->ijk * (sign ^ 1);

	int64_t fst = first->val;
	snd = second->val;
	sign = (uint64_t)(fst - snd) >> 63;

	int64_t fst_ijk = first->ijk;

	// swap first and second
	first->val = fst * sign + snd * (sign ^ 1);
	second->val = snd * sign + fst * (sign ^ 1);
	first->ijk = first->ijk * sign + second->ijk * (sign ^ 1);
	second->ijk = second->ijk * sign + fst_ijk * (sign ^ 1);
}

/* Compares three penalties and places the best
 * penalty at index 0. May destroy the others.
 *
 * Timing and cache access are independent of input.
 */
static void cmp_bestof3_penalties(Pen *ap, Pen *bp, Pen *cp)
{
	int64_t b = bp->val;
	int64_t c = cp->val;
	uint8_t sign = (uint64_t)(b - c) >> 63;

	bp->val = b * sign + c * (sign ^ 1);
	bp->ijk = bp->ijk * sign + cp->ijk * (sign ^ 1);

	int64_t a = ap->val;
	b = bp->val;
	sign = (uint64_t)(a - b) >> 63;

	ap->val = a * sign + b * (sign ^ 1);
	ap->ijk = ap->ijk * sign + bp->ijk * (sign ^ 1);
}

/* Compares two doubles and moves the minimum double
 * forward, along with an associated vectors in a
 * timing and cache access pattern consistent manner.
 * May destroy the larger of the pair (incl. vector).
 */
static void cmp_mov_with_vector(uint64_t *ap, uint64_t *bp, uint32_t *a_arr,
		uint32_t *b_arr, int n_arr)
{
	int64_t a = *ap;
	int64_t b = *bp;
	uint8_t sign = (uint64_t)(a - b) >> 63;

	*ap = a * sign + b * (sign ^ 1);

	int i;
	for (i = 0; i < n_arr; ++i) {
		a_arr[i] = a_arr[i] * sign + b_arr[i] * (sign ^ 1);
	}
}

/* Similar to cmp_mov_with_vector but moves associated
 * point instead of vector. May destroy the larger of
 * the pair (incl. point).
 */
static void cmp_mov_with_pt(uint64_t *ap, uint64_t *bp, uint64_t *a_pt,
		uint64_t *b_pt, int n)
{
	int64_t a = *ap;
	int64_t b = *bp;
	uint8_t sign = (uint64_t)(a - b) >> 63;

	*ap = a * sign + b * (sign ^ 1);
	*a_pt = *a_pt * sign + *b_pt * (sign ^ 1);
}

/* Moves the minimum value to the front of the array,
 * along with its associated array. Destroys both
 * inputs, excepting the minimum value and associated
 * array.
 */
static void min_in_array_with_vector(uint64_t *vals, uint32_t *arrs, int n, int n_arr)
{
	int i;
	for (i = n - 2; i >= 0; --i) {
		cmp_mov_with_vector(&vals[i], &vals[i + 1], &arrs[n_arr * i],
				&arrs[n_arr * (i + 1)], n_arr);
	}
}

/* Moves minimum metric to front of the array, along
 * with the associated point. Destroys both inputs,
 * excepting the minimum metric and associated point.
 */
static void min_metric(uint64_t *metrics, uint64_t *pts, int n)
{
	int i;
	for (i = n - 2; i >= 0; --i) {
		cmp_mov_with_pt(&metrics[i], &metrics[i + 1], &pts[i], &pts[i + 1], n);
	}
}

/* Decodes a X_ijk subset, for the 32-QAM constellation
 * (see comments above). This function will need to be
 * modified if another constellations is used.
 */
static void decode_subset(const uint32_t *t, int coset_h, int i, int j, int k, uint64_t *d, uint8_t *o)
{
	uint32_t offset[2][2] = {{0,0},{4,4}};
	uint64_t distances[2];
	
	int index;
	for (index = 0; index < 2; ++index) {
		distances[index] = dist_safe(t[0], t[1], SCALE * (p[coset_h][i][j][k][0] +
				offset[index][0]) % Q,
				SCALE * (p[coset_h][i][j][k][1] + offset[index][1]) % Q);
	}

	min_in_array_with_vector(distances, offset[0], 2, 2);
	*d = distances[0];
	*o ^= (offset[0][0] / 4) << (4 * i + 2 * j + k);
}

/* The precomputation step to calculate the d_ij's and
 * delta_ij's. This step is done for each Leech half-
 * lattice.
 */
static void precomputation_H24(const uint32_t *t, uint64_t d_ij[12][2][2],
		int64_t delta_ij[12][2][2], uint8_t offsets[12], int coset_h)
{
	// d_ij and delta_ij are set explicitly, but
	// the offsets need to be wiped.
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
				d_ij[n][i][j] = minu_safe(d_ij0, d_ij1);
				delta_ij[n][i][j] = d_ij1 - d_ij0;
			}
		}
	}
}

/* Initializes Q24 struct.
 */
static Q24 *init_Q24(uint64_t (*d_ij)[2][2], int64_t (*delta_ij)[2][2], uint8_t *offsets, int coset_h, int coset_q)
{
	Q24 *q = calloc(1, sizeof(struct Q24));

	q->d_ij = d_ij;
	q->delta_ij = delta_ij;
	q->offsets = offsets;

	q->coset_h = coset_h;
	q->coset_q = coset_q;

	int i, j, k;
	for (i = 0; i < 6; ++i) {
		for (j = 0; j < 4; ++j) {
			for (k = 0; k < 3; ++k) {
				// store penalty coordinate in four
				// most significant bits
				q->pen[i][j][k].ijk = (uint64_t)i << 60;
			}
		}
	}
	return q;
}

/* Computes confidence values, preferable representations,
 * and penalties. x refers to the character of GF(4).
 *
 * The preferable representations are calculated from the
 * (i1,j1,i2,j2) interpretation of c, and its complement.
 */
static void compute_vals(Q24 *q, int l, int x, int i1, int j1, int i2, int j2)
{
	int64_t a = q->d_ij[2 * l][i1][j1] + q->d_ij[2 * l + 1][i2][j2];
	int64_t b = q->d_ij[2 * l][!i1][!j1] + q->d_ij[2 * l + 1][!i2][!j2];
	int sign = (uint64_t)(a - b) >> 63;

	int k1, k2;
	// preferable representations
	i1 = i1 * sign + !i1 * (sign ^ 1);
	j1 = j1 * sign + !j1 * (sign ^ 1);
	k1 = (uint64_t)q->delta_ij[2 * l][i1][j1] >> 63;
	i2 = i2 * sign + !i2 * (sign ^ 1);
	j2 = j2 * sign + !j2 * (sign ^ 1);
	k2 = (uint64_t)q->delta_ij[2 * l + 1][i2][j2] >> 63;

	q->i1[l][x] = i1;
	q->j1[l][x] = j1;
	q->k1[l][x] = k1;
	q->i2[l][x] = i2;
	q->j2[l][x] = j2;
	q->k2[l][x] = k2;

	// confidence value
	q->mu_x[l][x] = a * sign + b * (sign ^ 1);

	// precalculate some penalty stuff
	int64_t delta_i1j1 = q->delta_ij[2 * l][0][0] * !i1 * !j1 +
			q->delta_ij[2 * l][0][1] * !i1 * j1 +
			q->delta_ij[2 * l][1][0] * i1 * !j1 +
			q->delta_ij[2 * l][1][1] * i1 * j1;
	int64_t delta_ni1j1 = q->delta_ij[2 * l][0][0] * i1 * j1 +
			q->delta_ij[2 * l][0][1] * i1 * !j1 +
			q->delta_ij[2 * l][1][0] * !i1 * j1 +
			q->delta_ij[2 * l][1][1] * !i1 * !j1;
	int64_t delta_i2j2 = q->delta_ij[2 * l + 1][0][0] * !i2 * !j2 +
			q->delta_ij[2 * l + 1][0][1] * !i2 * j2 +
			q->delta_ij[2 * l + 1][1][0] * i2 * !j2 +
			q->delta_ij[2 * l + 1][1][1] * i2 * j2;
	int64_t delta_ni2j2 = q->delta_ij[2 * l + 1][0][0] * i2 * j2 +
			q->delta_ij[2 * l + 1][0][1] * i2 * !j2 +
			q->delta_ij[2 * l + 1][1][0] * !i2 * j2 +
			q->delta_ij[2 * l + 1][1][1] * !i2 * !j2;

	uint64_t e = minu_safe(abs_safe(delta_i1j1), abs_safe(delta_i2j2));
	uint64_t f = abs_safe(a - b);
	uint64_t g = minu_safe(abs_safe(delta_ni1j1), abs_safe(delta_ni2j2));
	int parity = ((uint64_t)delta_i1j1 >> 63) ^
			((uint64_t)delta_i2j2 >> 63) ^
			((uint64_t)delta_ni1j1 >> 63) ^
			((uint64_t)delta_ni2j2 >> 63);
	int min_e = (uint64_t)((int64_t)abs_safe(delta_i1j1) -
			(int64_t)abs_safe(delta_i2j2)) >> 63;
	int min_g = (uint64_t)((int64_t)abs_safe(delta_ni1j1) -
			(int64_t)abs_safe(delta_ni2j2)) >> 63;
	int not_i1j1_k = (uint64_t)delta_ni1j1 >> 63;
	int not_i2j2_k = (uint64_t)delta_ni2j2 >> 63;
	
	// calculate penalties. ijk contains information on
	// penalty hexacode coordinate l in the most
	// significant bits, hence the XOR.
	q->pen[l][x][0].val = e;
	q->pen[l][x][0].ijk ^= (uint64_t)(0x08 * min_e ^ 0x01 * (min_e ^ 1)) << l * 6;
	q->pen[l][x][1].val = f + g * parity;
	q->pen[l][x][1].ijk ^= (uint64_t)(0x36 ^ (0x08 * (k1 ^ not_i1j1_k)) ^
				(0x01 * (k2 ^ not_i2j2_k)) ^ (0x08 * min_g ^ 0x01 * (min_g ^ 1)) *
				parity) << l * 6;
	q->pen[l][x][2].val = f + g * (parity ^ 1);
	q->pen[l][x][2].ijk ^= (uint64_t)(0x36 ^ (0x08 * (k1 ^ not_i1j1_k)) ^
				(0x01 * (k2 ^ not_i2j2_k)) ^ (0x08 * min_g ^ 0x01 * (min_g ^ 1)) *
				(parity ^ 1)) << l * 6;
}

/* Decodes the Leech quarter-lattice.
 */
static void decoder_Q24(Q24 *q, uint64_t *cv, uint64_t *d)
{
	// Step 1: Computing confidence values, penalties
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

	// Step 2: Sort the penalties
	// This step is done later due to timing concerns.

	// Step 3: Computing confidence values of block
	int x1, x2;
	for (l = 0; l < 3; ++l) {
		for (x1 = 0; x1 < 4; ++x1) {
			for (x2 = 0; x2 < 4; ++x2) {
				q->S_j[l][x1][x2] = q->mu_x[2 * l][x1] +
						q->mu_x[2 * l + 1][x2];
			}
		}
	}

	// Step 4: Finding images of hexacodewords, calculating
	// metrics
	int w;
	for (w = 0; w < 64; ++w) {
		uint8_t *word = hexacode[w];

		// calculate parities of preferable
		// representations
		int h_parity = 0;
		int k_parity = 0;
		for (l = 0; l < 6; ++l) {
			h_parity ^= q->i1[l][word[l]];
			k_parity ^= q->k1[l][word[l]] ^ q->k2[l][word[l]];
		}

		// parity disparities
		int h_d = h_parity ^ q->coset_q;
		int k_d = k_parity ^ q->coset_h;

		// get first and second best penalties
		// for each parity change
		Pen pens[3][2];
		int uv;
		for (uv = 0; uv < 3; ++uv) {
			// place penalties belonging to first two
			// hexacode digits in first and second best
			// slots
			cmp_place_penalties(&pens[uv][0], &pens[uv][1], &q->pen[0][word[0]][uv],
					&q->pen[1][word[1]][uv]);
			for (l = 2; l < 6; ++l) {
				// swap penalty in hexacode digit l into
				// first and/or second best slots if better
				cmp_swp_penalties(&pens[uv][0], &pens[uv][1],
						&q->pen[l][word[l]][uv]);
			}
		}

		// penalty collisions
		int h_c = (pens[0][0].ijk >> 60) == (pens[2][0].ijk >> 60);
		int k_c = (pens[1][0].ijk >> 60) == (pens[2][0].ijk >> 60);
		int hk_c = (pens[0][0].ijk >> 60) == (pens[1][0].ijk >> 60);

		// cases to fix penalties: 1 coordinate fix,
		// 2 coordinate fix with priority on h or k,
		// 2 coordinate fix with priority on hk or k
		Pen cases[3];
		memset(cases, 0, 3 * sizeof(Pen));
		cases[0].val = h_d * (k_d ^ 1) * pens[1][0].val +
				(h_d ^ 1) * k_d * pens[0][0].val +
				h_d * k_d * pens[2][0].val;
		cases[0].ijk = h_d * (k_d ^ 1) * pens[1][0].ijk ^
				(h_d ^ 1) * k_d * pens[0][0].ijk ^
				h_d * k_d * pens[2][0].ijk;
		cases[1].val = h_d * (k_d ^ 1) * (pens[0][0].val +
				(h_c ^ 1) * pens[2][0].val +
				h_c * pens[2][1].val) + 
				(h_d ^ 1) * k_d * (pens[1][0].val +
				(k_c ^ 1) * pens[2][0].val +
				k_c * pens[2][1].val) +
				h_d * k_d * (pens[0][0].val +
				(hk_c ^ 1) * pens[1][0].val +
				hk_c * pens[1][1].val);
		cases[1].ijk = h_d * (k_d ^ 1) * (pens[0][0].ijk ^
				(h_c ^ 1) * pens[2][0].ijk ^
				h_c * pens[2][1].ijk) ^
				(h_d ^ 1) * k_d * (pens[1][0].ijk ^
				(k_c ^ 1) * pens[2][0].ijk ^
				k_c * pens[2][1].ijk) ^
				h_d * k_d * (pens[0][0].ijk ^
				(hk_c ^ 1) * pens[1][0].ijk ^
				hk_c * pens[1][1].ijk);
		cases[2].val = h_d * (k_d ^ 1) * (pens[2][0].val +
				(h_c ^ 1) * pens[0][0].val +
				h_c * pens[0][1].val) + 
				(h_d ^ 1) * k_d * (pens[2][0].val +
				(k_c ^ 1) * pens[1][0].val +
				k_c * pens[1][1].val) +
				h_d * k_d * (pens[1][0].val +
				(hk_c ^ 1) * pens[0][0].val +
				hk_c * pens[0][1].val);
		cases[2].ijk = h_d * (k_d ^ 1) * (pens[2][0].ijk ^
				(h_c ^ 1) * pens[0][0].ijk ^
				h_c * pens[0][1].ijk) ^
				(h_d ^ 1) * k_d * (pens[2][0].ijk ^
				(k_c ^ 1) * pens[1][0].ijk ^
				k_c * pens[1][1].ijk) ^
				h_d * k_d * (pens[1][0].ijk ^
				(hk_c ^ 1) * pens[0][0].ijk ^
				hk_c * pens[0][1].ijk);

		// get best case at index 0
		cmp_bestof3_penalties(&cases[0], &cases[1], &cases[2]);

		// calculate metric
		q->M[w] = q->S_j[0][word[0]][word[1]] + q->S_j[1][word[2]][word[3]] +
				q->S_j[2][word[4]][word[5]] + cases[0].val;

		// construct point
		q->pts[w] = 0;
		for (l = 0; l < 6; ++l) {
			q->pts[w] ^= (uint64_t)q->i1[l][word[l]] << (l * 6 + 5);
			q->pts[w] ^= (uint64_t)q->j1[l][word[l]] << (l * 6 + 4);
			q->pts[w] ^= (uint64_t)q->k1[l][word[l]] << (l * 6 + 3);
			q->pts[w] ^= (uint64_t)q->i2[l][word[l]] << (l * 6 + 2);
			q->pts[w] ^= (uint64_t)q->j2[l][word[l]] << (l * 6 + 1);
			q->pts[w] ^= (uint64_t)q->k2[l][word[l]] << (l * 6);
		}
		q->pts[w] ^= cases[0].ijk;
	}
	
	// Step 5: Final minimization
	min_metric(q->M, q->pts, 64);

	// Finalize output
	*d = q->M[0];

	// 36 bits of point information (12 redundant)
	*cv = q->pts[0];

	// 12 bits of offset information (0 redundant)
	for (l = 0; l < 6; ++l) {
		// point is represented in groups of six,
		// little endian
		uint8_t ijk = *cv >> (6 * l);
		// for each group of six, output 2 offset bits
		*cv ^= (uint64_t)(q->offsets[2 * l] >> (ijk >> 3 & 0x7) & 1) << (36 + 2 * l);
		*cv ^= (uint64_t)(q->offsets[2 * l + 1] >> (ijk & 0x7) & 1) << (36 + 2 * l + 1);
	}

	/* Note that the actual output will depend on how the
	 * points are encoded. Some of the bits in *cv are
	 * redundant, due to the parities, and hexacode
	 * representation. In addition, some bits are not yet
	 * outputted, such as the Leech half-lattice coset,
	 * and the offsets recorded in the precomputation.
	 *
	 * How these bits are outputted and in what order can
	 * be chosen arbitrarily, and is dependent on the
	 * specifics of the complete implementation. For the
	 * use of large constellations, for example, the number
	 * of outputted bits can vastly exceed 64. For this
	 * proof of concept, we output 12 offset bits, plus 36
	 * point bits (of which 12 are redundant).
	 */
}

/* Decodes the Leech lattice by means of four
 * Leech quarter-lattice decoders.
 *  *t: pointer to target vector
 * *cv: pointer to outputted bits
 *  *d: pointer to outputted SED
 */
void decoder_L24(const uint32_t *t, uint64_t *cv, uint64_t *d)
{
	uint64_t cv_q[4];
	uint64_t d_q[4];

	uint64_t d_ij[12][2][2];
	int64_t delta_ij[12][2][2];
	uint8_t offsets[12];

	int coset_h, coset_q;
	for (coset_h = 0; coset_h < 2; ++coset_h) {
		precomputation_H24(t, d_ij, delta_ij, offsets, coset_h);

		for (coset_q = 0; coset_q < 2; ++coset_q) {
			int coset = 2 * coset_h + coset_q;
			Q24 *q = init_Q24(d_ij, delta_ij, offsets, coset_h, coset_q);
			decoder_Q24(q, &cv_q[coset], &d_q[coset]);
			free(q);
		}
	}

	min_metric(d_q, cv_q, 4);

	*d = d_q[0];
	*cv = cv_q[0] & 0x00ffffffffffffff;
}
