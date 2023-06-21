restart
needsPackage "TerraciniLoci"

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
