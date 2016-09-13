#ifndef LEECH_H
#define LEECH_H

/* This defines the modulus Q. Since vectors are
 * stored in uint32_t types, and SEQs in uint64_t
 * types, the modulus should be:
 *
 *   floor(log_2(sqrt(2^64 / 24))) = 29
 *
 * bits or less.
 *
 * To ensure that the Leech lattice remains integer,
 * Q should be divisible by 8.
 */
#define Q     1024

/* SCALE is the size of the smallest non-zero
 * coordinate of a vector in the Leech lattice,
 * which is dependent on Q, and the number of
 * Leech lattice layers we wish to use. The number
 * of layers determines how many points of the
 * Leech lattice we used.
 *
 * This code is built around the following:
 *
 *   Q*ZZ^24_Q \subset (Q/sqrt(8))*L \subset ZZ^24_Q
 *
 * where L is the Leech lattice with minimal vector
 * of norm 2. This gives a total of 2^36 Leech
 * lattice points, and means SCALE is defined as
 * follows.
 *
 * For the CVP decoder by Vardy and Be'ery, this
 * means a square 32*QAM constellation is used.
 */
#define SCALE (Q / 8)

extern uint8_t p[2][2][2][2][2];
extern uint8_t hexacode[64][6];

/* Maximum likelihood decoder of the Leech lattice.
 * This function is timing safe.
 *    t: input vector
 *   cv: output bits (corresponding to a point in the
 *       leech lattice)
 *    d: output distance (squared euclidean distance
 *       between input/output vectors)
 */
void decoder_L24(const uint32_t *t, uint64_t *cv, uint64_t *d);

/* Maximum likelihood decoder of the Leech lattice.
 * This function is not cryptographically secure as it
 * leaks timing information.
 *    t: input vector
 *   cv: output bits (corresponding to a point in the
 *       leech lattice)
 *    d: output distance (squared euclidean distance
 *       between input/output vectors)
 */
void decoder_L24_unsafe(const uint32_t *t, uint64_t *cv, uint64_t *d);

#endif
