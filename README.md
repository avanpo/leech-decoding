# Leech lattice decoding

The Leech lattice is a very remarkable lattice in 24 dimensions. Its density and high degree of structure make it a very suitable candidate for error-correction.

This library aims to efficiently implement several different Leech lattice decoders.

For a discussion of the original, cryptographically safe integer implementation, please see my thesis at the Cryptology ePrint [archive](https://eprint.iacr.org/2016/1050).

## Contents

- *leech.h:* Header file containing constants and declarations of all integer decoders.
- *leech.c:* Non-vectorized, cryptographically safe, integer implementation of [1]. If a tighter bound on the modulus `q` can be assumed, this implementation can be optimized further.
- *leech_unsafe.c:* Non-vectorized integer implementation of [1]. Some minor optimizations are probably possible.

### To do

- Vectorized versions of both the above.
- Floating point versions of the unsafe version.
- Bounded distance decoder, based on [2].

## Maximum likelihood decoding

Maximum likelihood decoding is equivalent to the Closest Vector Problem (CVP) on lattices. Given a point in space, it finds the closest lattice vector to this point. Alternatively, it maps every point in the Voronoi cell of a lattice point to that point, except those points on the boundary.

The most efficient decoder known is described in [1]. This library contains a cryptographically secure implementation, and an unsafe version.

## Bounded distance decoding

Bounded distance decoding is equivalent to maximum likelihood decoding, with the additional restriction that the input point lies within a certain distance to a lattice point. For the Leech lattice, this distance is equal to half the minimal distance. Essentially, it maps each point within a sphere packing of the Leech lattice to the relevant point.

The most efficient bounded distance decoder known is based on the same construction as [1]. It is described in [2].

[1] Vardy and Be'ery, "Maximum Likelihood Decoding of the Leech lattice", Trans. on Information Theory, Vol. 39, No. 4, July 1993.

[2] Vardy, "Even More Efficient Bounded-Distance Decoding of the Hexacode, the Golay Code, and the Leech Lattice", Trans. on Information Theory, Vol. 41, No. 5, September 1995.
