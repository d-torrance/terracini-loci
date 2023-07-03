restart
needsPackage "TerraciniLoci"

----------------------------
-- rational normal curves --
----------------------------

rationalNormalCurve = d -> (
    kk := ZZ/32003;
    (s, t, x) := (symbol s, symbol t, symbol x);
    R := kk[s, t];
    S := kk[x_0..x_d];
    map(R, S, apply(d + 1, i -> s^(d - i) * t^i)))

assertEmptyTerracini = (r, f) -> elapsedTime assert(terraciniLocus(r, f) == 1)

-- twisted cubic
assertEmptyTerracini(2, rationalNormalCurve 3)

-- rational normal quartic
assertEmptyTerracini(2, rationalNormalCurve 4)

-- rational normal quintic
assertEmptyTerracini(2, rationalNormalCurve 5)
assertEmptyTerracini(3, rationalNormalCurve 5)

-- rational normal sextic
assertEmptyTerracini(2, rationalNormalCurve 6)
assertEmptyTerracini(3, rationalNormalCurve 6)

-- rational normal septic
assertEmptyTerracini(2, rationalNormalCurve 7)
assertEmptyTerracini(3, rationalNormalCurve 7)
assertEmptyTerracini(4, rationalNormalCurve 7)

----------------------------
-- elliptic normal curves --
----------------------------

-- elliptic normal quintic
kk = ZZ/32003;
S = kk[y_0..y_9];
A = genericSkewMatrix(S, 5);
R = kk[x_0..x_4];
I = pfaffians(4, sub(A, apply(10, i -> S_i => random(1, R))));
assert(dim I == 2)
assert(genus I == 1)
assert(degree I == 5)
assertEmptyTerracini(2, I)

------------------------
-- del Pezzo surfaces --
------------------------

delPezzoSurface = t -> (
    kk := ZZ/32003;
    d := 9 - t;
    (x, y) := (symbol x, symbol y);
    R := kk[y_0..y_2];
    S := kk[x_0..x_d];
    P := intersect \\ ideal \ {
	{y_0, y_1}, {y_0, y_2}, {y_1, y_2}, {y_0 - y_1, y_0 - y_2}
	}_{0..t - 1};
    map(R, S, super basis(3, P)))

-- Corollary 5.5
assertCorollary55 = t -> (
    I := elapsedTime terraciniLocus(2, delPezzoSurface t);
    comps := primaryDecomposition I;
    assert(#comps == (if t == 4 then 5 else t) and
	all(comps, J -> dim J - 2 == 3)))

assertCorollary55 1
assertCorollary55 2
assertCorollary55 3
assertCorollary55 4

-- Corollary 5.7
I = elapsedTime terraciniLocus(3, delPezzoSurface 1);
comps = primaryDecomposition I;
assert(#comps == 4 and all(comps, J -> dim J - 3 == 5))
