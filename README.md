# TerraciniLoci package for Macaulay2

This package implements the algorithms from Section 8 of the paper
*Geometry of first nonempty Terracini loci* by F. Galuppi,
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

### Rational quintic in P^4 (Example 4.5)
```m2
i3 : kk = ZZ/32003;

i4 : R = kk[x, y];

i5 : S = kk[z_0..z_4];

i6 : f = map(R, S, {x^5, x^4*y, x^3*y^2, x*y^4, y^5});

o6 : RingMap R <--- S

i7 : assertEmptyTerracini(2, f)
 -- 0.167679 seconds elapsed
 ```

### Elliptic normal quintic

```m2
i3 : kk = ZZ/32003;

i4 : S = kk[y_0..y_9];

i5 : A = genericSkewMatrix(S, 5);

             5       5
o5 : Matrix S  <--- S

i6 : R = kk[x_0..x_4];

i7 : I = pfaffians(4, sub(A, apply(10, i -> S_i => random(1, R))));

o7 : Ideal of R

i8 : assert(dim I == 2)

i9 : assert(genus I == 1)

i10 : assert(degree I == 5)

i11 : assertEmptyTerracini(2, I)
 -- 32.506 seconds elapsed
 ```

### del Pezzo surfaces

```m2
i1 : delPezzoSurface = t -> (
    kk := ZZ/32003;
    d := 9 - t;
    (x, y) := (symbol x, symbol y);
    R := kk[y_0..y_2];
    S := kk[x_0..x_d];
    P := intersect \\ ideal \ {
	{y_0, y_1}, {y_0, y_2}, {y_1, y_2}, {y_0 - y_1, y_0 - y_2}
	}_{0..t - 1};
    map(R, S, super basis(3, P)))

o1 = delPezzoSurface

o1 : FunctionClosure

i2 : -- Corollary 5.5
assertCorollary55 = t -> (
    I := elapsedTime terraciniLocus(2, delPezzoSurface t);
    comps := primaryDecomposition I;
    assert(#comps == (if t == 4 then 5 else t) and
	all(comps, J -> dim J - 2 == 3)))

o2 = assertCorollary55

o2 : FunctionClosure

i3 : assertCorollary55 1
 -- 0.557642 seconds elapsed

i4 : assertCorollary55 2
 -- 0.443816 seconds elapsed

i5 : assertCorollary55 3
 -- 0.512567 seconds elapsed

i6 : assertCorollary55 4
 -- 0.711542 seconds elapsed

i7 : -- Corollary 5.7
I = elapsedTime terraciniLocus(3, delPezzoSurface 1);
 -- 230.215 seconds elapsed

                ZZ
o7 : Ideal of -----[z   ..z   ]
              32003  0,0   2,2

i8 : comps = primaryDecomposition I;

i9 : assert(#comps == 4 and all(comps, J -> dim J - 3 == 5))
```

### Veronese varieties

```m2
i2 : needsPackage "Resultants"
 -- warning: symbol "resultant" in Elimination.Dictionary is shadowed by a symbol in Resultants.Dictionary
 --   use the synonym Elimination$resultant
 -- warning: symbol "discriminant" in Elimination.Dictionary is shadowed by a symbol in Resultants.Dictionary
 --   use the synonym Elimination$discriminant

o2 = Resultants

o2 : Package

i3 : assertEmptyTerracini(2, veronese(2, 3))
 -- 0.730559 seconds elapsed

i4 : I = elapsedTime terraciniLocus(3, veronese(2, 3));
 -- 617.377 seconds elapsed

o4 : Ideal of QQ[z   ..z   ]
                  0,0   2,2

i5 : assert(#primaryDecomposition I == 1 and dim I - 3 == 2 * 2 + 3 - 2)
```
