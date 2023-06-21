# TerraciniLoci package for Macaulay2

This package implements the algorithms from Section 8 of the paper
*Terracini loci of interesting surfaces* by F. Galuppi,
P. Santarsiero, D. Torrance, and E. Turatti.

## Examples

### Rational normal curves

```m2
i1 : needsPackage "TerraciniLoci"

o1 = TerraciniLoci

o1 : Package

i2 : rationalNormalCurve = d -> (
    kk := ZZ/32003;
    (s, t, x) := (symbol s, symbol t, symbol x);
    R := kk[s, t];
    S := kk[x_0..x_d];
    map(R, S, apply(d + 1, i -> s^(d - i) * t^i)))

o2 = rationalNormalCurve

o2 : FunctionClosure

i3 : assertEmptyTerracini = (r, f) -> elapsedTime assert(terraciniLocus(r, f) == 1)

o3 = assertEmptyTerracini

o3 : FunctionClosure

i4 : -- twisted cubic
     assertEmptyTerracini(2, rationalNormalCurve 3)
 -- 0.129286 seconds elapsed

i5 : -- rational normal quartic
     assertEmptyTerracini(2, rationalNormalCurve 4)
 -- 0.167151 seconds elapsed

i6 : -- rational normal quintic
     assertEmptyTerracini(2, rationalNormalCurve 5)
 -- 0.194115 seconds elapsed

i7 : assertEmptyTerracini(3, rationalNormalCurve 5)
 -- 0.425239 seconds elapsed

i8 : -- rational normal sextic
     assertEmptyTerracini(2, rationalNormalCurve 6)
 -- 0.223492 seconds elapsed

i9 : assertEmptyTerracini(3, rationalNormalCurve 6)
 -- 1.04454 seconds elapsed

i10 : -- rational normal septic
      assertEmptyTerracini(2, rationalNormalCurve 7)
 -- 0.19556 seconds elapsed

i11 : assertEmptyTerracini(3, rationalNormalCurve 7)
 -- 2.20825 seconds elapsed

i12 : assertEmptyTerracini(4, rationalNormalCurve 7)
 -- 1028.94 seconds elapsed
 ```
