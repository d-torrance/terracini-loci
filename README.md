# TerraciniLoci package for Macaulay2

This package implements the algorithms from Section 8 of the paper
*Geometry of first nonempty Terracini loci* by F. Galuppi,
P. Santarsiero, D. Torrance, and E. Turatti.

## Examples

### Curves

#### Rational normal curves (Example 4.2)

Rational normal curves have empty Terracini loci.

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

When possible, it is better to use the `terraciniLocus` method with a ring
map than an ideal.  In addition to only working for 2nd Terracini loci, the
ideal method is also much slower.

```m2
i4 : assertEmptyTerracini(2, ker rationalNormalCurve 3)
 -- 1.49868 seconds elapsed

i5 : assertEmptyTerracini(2, ker rationalNormalCurve 4)
 -- 14.4332 seconds elapsed
```

#### Elliptic normal quintic (Example 4.2 continued)

Elliptic normal curves also have empty Terracini loci.  For this
example, we must use the ideal method.

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

####  Rational octic in $\mathbb P^7$ (Corollary 4.3)

The bound from this result is not sharp.  Indeed, it only guarantees that the 2nd Terracini locus of this curve is empty, but we see below that the 3rd Terracini locus is empty as well.

```m2
i3 : kk = ZZ/32003;

i4 : R = kk[x, y];

i5 : S = kk[z_0..z_7];

i6 : f = map(R, S, {x^8, x^7*y, x^6*y^2, x^5*y^3, x^4*y^4, x^3*y^5, x*y^7, y^8});

o6 : RingMap R <--- S

i7 : assertEmptyTerracini(2, f)
 -- 0.205354 seconds elapsed

i8 : assertEmptyTerracini(3, f)
 -- 10.1608 seconds elapsed
```

#### Rational quintic in $\mathbb P^4$ (Example 4.5)

A counterexample to Proposition 4.4 in even-dimensional space, i.e., a curve with an empty last Terracini locus that is not the rational normal curve.

```m2
i3 : kk = ZZ/32003;

i4 : R = kk[x, y];

i5 : S = kk[z_0..z_4];

i6 : f = map(R, S, {x^5, x^4*y, x^3*y^2, x*y^4, y^5});

o6 : RingMap R <--- S

i7 : assertEmptyTerracini(2, f)
 -- 0.167679 seconds elapsed
 ```


### del Pezzo surfaces

#### 2nd Terracini loci of del Pezzo surfaces (Corollary 5.5)

The 2nd Terracini locus of the blowup of $t$ points in $\mathbb P^2$ has $t$ 3-dimensional irreducible components when $t\in\{1,2,3\}$ and five 3-dimensional irreducible components when $t = 4$.

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
```

#### 3nd Terracini locus of a del Pezzo surface (Corollary 5.7)

The 3rd Terracini locus of the blowup of a single point in $\mathbb P^2$ has four 5-dimensional irreducible components.

```m2
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

#### 2nd Terracini locus of the Veronese cubic surface (Proposition 6.1)

When $2r < d + 2$, the $`r`$th Terracini locus of $V_n^d$ is empty.

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
```

#### 3rd Terracini locus of the Veronese cubic surface (Theorem 6.4)

When $r = \left\lceil\frac{d+2}{2}\right\rceil$, the $`r`$th Terracini locus of $V_n^d$ is irreducible of dimension $2n + r - 2$.

```m2

i4 : I = elapsedTime terraciniLocus(3, veronese(2, 3));
 -- 617.377 seconds elapsed

o4 : Ideal of QQ[z   ..z   ]
                  0,0   2,2

i5 : assert(#primaryDecomposition I == 1 and dim I - 3 == 2 * 2 + 3 - 2)
```
