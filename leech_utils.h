#ifndef LEECH_UTILS_H
#define LEECH_UTILS_H

extern uint8_t leech[24][24];

/* computation */
void copyv(uint32_t *dest, const uint32_t *src);
int cmpv(const uint32_t *a, const uint32_t *b);
uint64_t minu(uint64_t u1, uint64_t u2);
uint64_t absi(int64_t i);
uint64_t sed(const uint32_t *a, const uint32_t *b);

/* hexacode */
uint8_t projection(uint8_t i1, uint8_t j1, uint8_t i2, uint8_t j2);

/* encoding */
void decode_pt(uint32_t *out, uint64_t cv, int verbose);

/* printing */
void printvu(const uint32_t *v);
void printvb(const uint64_t v, uint8_t spacing);

#endif
