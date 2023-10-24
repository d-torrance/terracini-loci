# TerraciniLoci package for Macaulay2

This package implements the algorithms from Section 8 of the paper
*Geometry of first nonempty Terracini loci* by F. Galuppi,
P. Santarsiero, D. Torrance, and E. Turatti.

## Description

This package contains a method `terraciniLocus` with two installed method functions.

### Parametrization

If $X\subset\mathbb P^m$ is parametrized by the rational map $f:\mathbb P^n\dashrightarrow\mathbb P^m$ (represented in Macaulay2 by a ring map from the coordinate ring of $\mathbb P^m$ to the coordinate ring of $\mathbb P^n$) then `terraciniLocus(r, f)` will give the ideal of the preimage of the $`r`$th Terracini locus of $X$ in $(\mathbb P^n)^r$.

### Ideal

If $X\subset\mathbb P^n$ has ideal $I$, then `terraciniLocus(2, I)` will give the ideal of the preimage of the 2nd Terracini locus of $X$ in $\mathbb P^n\times\mathbb P^n$.

## Installation

To install this package, download the file [`TerraciniLoci.m2`](https://github.com/d-torrance/terracini-loci/raw/master/TerraciniLoci.m2) to some local directory, say `~/Downloads`.

Then run the following in Macaulay2:

```m2
installPackage("TerraciniLoci", FileName => "~/Downloads/TerraciniLoci.m2")
```

## Examples

### Curves

#### Rational normal curves (Example 4.2)

Rational normal curves have empty Terracini loci.

We load the [`Resultants`](https://doi.org/10.2140/jsag.2018.8.21) package so that we can call the [`veronese`](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/Resultants/html/_veronese.html) function.  This returns the ring map corresponding to the Veronese embedding.

```m2
i1 : needsPackage "TerraciniLoci"

o1 = TerraciniLoci

o1 : Package

i2 : assertEmptyTerracini = (r, f) -> elapsedTime assert(terraciniLocus(r, f) == 1)

o2 = assertEmptyTerracini

o2 : FunctionClosure

i3 : needsPackage "Resultants"
 -- warning: symbol "resultant" in Elimination.Dictionary is shadowed by a symbol in Resultants.Dictionary
 --   use the synonym Elimination$resultant
 -- warning: symbol "discriminant" in Elimination.Dictionary is shadowed by a symbol in Resultants.Dictionary
 --   use the synonym Elimination$discriminant

o3 = Resultants

o3 : Package

i4 : -- twisted cubic
     assertEmptyTerracini(2, veronese(1, 3))
 -- 0.168295 seconds elapsed

i5 : -- rational normal quartic
     assertEmptyTerracini(2, veronese(1, 4))
 -- 0.174667 seconds elapsed

i6 : -- rational normal quintic
     assertEmptyTerracini(2, veronese(1, 5))
 -- 0.179793 seconds elapsed

i7 : assertEmptyTerracini(3, veronese(1, 5))
 -- 0.503009 seconds elapsed

i8 : -- rational normal sextic
     assertEmptyTerracini(2, veronese(1, 6))
 -- 0.17569 seconds elapsed

i9 : assertEmptyTerracini(3, veronese(1, 6))
 -- 0.916899 seconds elapsed

i10 : -- rational normal septic
      assertEmptyTerracini(2, veronese(1, 7))
 -- 0.249929 seconds elapsed

i11 : assertEmptyTerracini(3, veronese(1, 7))
 -- 2.11239 seconds elapsed

i12 : assertEmptyTerracini(4, veronese(1, 7))
 -- 877.63 seconds elapsed
 ```

When possible, it is better to use the `terraciniLocus` method with a ring
map than an ideal.  In addition to only working for 2nd Terracini loci, the
ideal method is also much slower.

```m2
i13 : assertEmptyTerracini(2, ker veronese(1, 3))
 -- 1.48714 seconds elapsed

i14 : assertEmptyTerracini(2, ker veronese(1, 4))
 -- 12.9082 seconds elapsed
```

#### Elliptic normal quintic (Example 4.2 continued)

Elliptic normal curves in even-dimensional space also have empty
Terracini loci.  For this example, we must use the ideal method.

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

####  Rational octic in projective 7-space (Corollary 4.3)

The bound from this result is not sharp.  Indeed, since $3r < 7 + 2$ implies $r\leq 2$, it only guarantees that the 2nd Terracini locus of this curve is empty, but we see below that the 3rd Terracini locus is empty as well.

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

#### Rational quintic in projective 4-space (Example 4.5)

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

#### Rational quartic in projective 3-space

So far, we have only looked at curves with empty Terracini loci.  Now consider
a rational quartic in $\mathbb P^3$.  Since it's not a rational normal curve,
we know by Proposition 4.4 that its 2nd Terracini locus is nonempty.

```m2
i1 : kk = ZZ/32003;

i2 : R = kk[x, y];

i3 : S = kk[z_0..z_3];

i4 : f = map(R, S, {x^4, x^3*y, x*y^3, y^4});

o4 : RingMap R <--- S

i5 : T1 = elapsedTime terraciniLocus(2, f)
 -- 0.0860385 seconds elapsed

            2   2                          2   2
o5 = ideal(z   z    + 4z   z   z   z    + z   z   )
            0,1 1,0     0,0 0,1 1,0 1,1    0,0 1,1

o5 : Ideal of kk[z   ..z   ]
                  0,0   1,1
```

Let's confirm that this also works using the ideal method.

```m2
i6 : T2 = elapsedTime terraciniLocus(2, ker f);
 -- 2.84349 seconds elapsed

o6 : Ideal of kk[z   ..z   ]
                  0,0   1,3

i7 : betti T2

            0  1
o7 = total: 1 20
         0: 1  .
         1: .  4
         2: . 14
         3: .  2

o7 : BettiTally
```

At first glance, these results seem to disagree.  But $T_1$ is the ideal of
a variety in $\mathbb P^1\times\mathbb P^1$ and $T_2$ is the ideal of a variety
in $\mathbb P^3\times\mathbb P^3$.  If we look at the image of $T_2$ under
our map and saturate out with respect to the irrelevant ideals, we see that
we indeed get the same thing.

```m2
i8 : use ring T1;

i9 : g = map(ring T1, ring T2, sub(matrix f, {x => z_(0,0), y => z_(0, 1)}) |
         sub(matrix f, {x => z_(1,0), y => z_(1, 1)}));

o9 : RingMap kk[z   ..z   ] <--- kk[z   ..z   ]
                 0,0   1,1           0,0   1,3

i10 : saturate(g T2, intersect(ideal(z_(0, 0), z_(0, 1)), ideal(z_(1, 0), z_(1, 1))))

             2   2                          2   2
o10 = ideal(z   z    + 4z   z   z   z    + z   z   )
             0,1 1,0     0,0 0,1 1,0 1,1    0,0 1,1

o10 : Ideal of kk[z   ..z   ]
                   0,0   1,1
```

### del Pezzo surfaces

Most del Pezzo surfaces are blowups of $t$ points in $\mathbb P^2$.  We embed them in $\mathbb P^{9-t}$ via the linear system of cubics that vanish on these $t$ points.  Since we are only interested in $t\leq 4$, we will assume the set of blowup points is a subset of $`\{[0:0:1], [0:1:0], [1:0:0], [1:1:1]\}`$.

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
```

#### 2nd Terracini loci of del Pezzo surfaces (Corollary 5.5)

The 2nd Terracini locus of the blowup of $t$ points in $\mathbb P^2$ has $t$ 3-dimensional irreducible components when $`t\in\{1,2,3\}`$ and five 3-dimensional irreducible components when $t = 4$.

The blowup of a single point has an additional 2-dimensional component, the set of pairs of points on the exceptional divisor, but this component disappears when pulling back to the plane, so it is missing when using the parametrization code.

```m2
i3 : elapsedTime apply(primaryDecomposition terraciniLocus(2, delPezzoSurface 1),
         I -> dim I - 2)
 -- 0.299992 seconds elapsed

o3 = {3}

o3 : List

i4 : elapsedTime apply(primaryDecomposition terraciniLocus(2, delPezzoSurface 2),
         I -> dim I - 2)
 -- 0.303849 seconds elapsed

o4 = {3, 3}

o4 : List

i5 : elapsedTime apply(primaryDecomposition terraciniLocus(2, delPezzoSurface 3),
         I -> dim I - 2)
 -- 0.330507 seconds elapsed

o5 = {3, 3, 3}

o5 : List

i6 : elapsedTime apply(primaryDecomposition terraciniLocus(2, delPezzoSurface 4),
         I -> dim I - 2)
 -- 0.649615 seconds elapsed

o6 = {3, 3, 3, 3, 3}

o6 : List
```

We check our work using the ideal method.  It takes considerably longer, and the computer may run out of memory computing the cases with $`t\in\{1,2\}`$.

```m2
i2 : elapsedTime apply(primaryDecomposition terraciniLocus(2, ker delPezzoSurface 3),
         I -> dim I - 2)
 -- 4604.65 seconds elapsed

o2 = {3, 3, 3}

o2 : List

i3 : elapsedTime apply(primaryDecomposition terraciniLocus(2, ker delPezzoSurface 4),
         I -> dim I - 2)
 -- 48960. seconds elapsed

o3 = {3, 3, 3, 3, 3}

o3 : List
```

#### 3rd Terracini locus of a del Pezzo surface (Corollary 5.7)

The 3rd Terracini locus of the blowup of a single point in $\mathbb P^2$ has four 5-dimensional irreducible components.

```m2
i7 : elapsedTime apply(primaryDecomposition terraciniLocus(3, delPezzoSurface 1),
         I -> dim I - 3)
 -- 175.677 seconds elapsed

o7 = {5, 5, 5, 5}
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

## Segre-Veronese varieties

### Segre-Veronese surfaces (Theorem 7.6)

If $`r = \left\lceil\frac{d_1+2}{2}\right\rceil`$ and $`\hat\jmath=\max\left\{i\in\{1,\ldots,k\}:\left\lceil\frac{d_i+2}{2}\right\rceil=r\right\}`$, then the $`r`$th Terracini locus of the Segre-Veronese embedding of $`\mathbb P^{n_1}\times\cdots\times\mathbb P^{n_k}`$ via $`\mathcal O(d_1,\ldots,d_k)`$ has $`\hat\jmath`$ components of dimension $`n_1+\dots+n_k+n_i + r - 2`$.

We demonstrate this for some Segre-Veronese surfaces, and in particular, for $`\mathbb P^1\times\mathbb P^1`$ embedded via $`\mathcal O(1, 2)`$, $`\mathcal O(1, 3)`$, and $`\mathcal O(2, 2)`$.  The latter is the del Pezzo surface mentioned at the beginning of Section 5.

```m2
i2 : segreVeronese = (n, d) -> (
    x := symbol x;
    r := #n;
    R := QQ new Array from splice apply(r, i -> x_(i, 0)..x_(i, n#i));
    y := symbol y;
    S := QQ[y_0..y_(product(n, d, (ni, di) -> binomial(ni + di, di)) - 1)];
    map(R, S, flatten entries first tensor apply(r, i -> (
		vector apply(subsets(n#i + d#i, d#i), A -> product(d#i, j ->
			x_(i, A#j - j)))))))


o2 = segreVeronese

o2 : FunctionClosure

i3 : assertTheorem76 = (n, d) -> (
    r := ceiling((d#0 + 2)/2);
    I := elapsedTime terraciniLocus(r, segreVeronese(n, d));
    jHat := max select(#d, i -> ceiling((d#i + 2)/2) == r) + 1;
    comps := primaryDecomposition I;
    assert(#comps == jHat);
    assert(sort apply(comps, J -> dim J - 2 - r) ==
	sort apply (jHat, i -> sum n + n#i + r - 2)))


o3 = assertTheorem76

o3 : FunctionClosure

i4 : assertTheorem76({1, 1}, {1, 2})
 -- 0.590006 seconds elapsed

i5 : assertTheorem76({1, 1}, {1, 3})
 -- 3.76348 seconds elapsed

i6 : assertTheorem76({1, 1}, {2, 2})
 -- 11.7846 seconds elapsed
```
