-- this file contains the code used for the examples in README.md

needsPackage "TerraciniLoci"

----------------------------
-- rational normal curves --
----------------------------

needsPackage "Resultants"

assertEmptyTerracini = (r, f) -> elapsedTime assert(terraciniLocus(r, f) == 1)

-- twisted cubic
assertEmptyTerracini(2, veronese(1, 3))

-- rational normal quartic
assertEmptyTerracini(2, veronese(1, 4))

-- rational normal quintic
assertEmptyTerracini(2, veronese(1, 5))
assertEmptyTerracini(3, veronese(1, 5))

-- rational normal sextic
assertEmptyTerracini(2, veronese(1, 6))
assertEmptyTerracini(3, veronese(1, 6))

-- rational normal septic
assertEmptyTerracini(2, veronese(1, 7))
assertEmptyTerracini(3, veronese(1, 7))
assertEmptyTerracini(4, veronese(1, 7))

-- using the ideal method
assertEmptyTerracini(2, ker veronese(1, 3))
assertEmptyTerracini(2, ker veronese(1, 4))

---------------------------
-- other rational curves --
---------------------------

-- corollary 4.3 not sharp
kk = ZZ/32003;
R = kk[x, y];
S = kk[z_0..z_7];
f = map(R, S, {x^8, x^7*y, x^6*y^2, x^5*y^3, x^4*y^4, x^3*y^5, x*y^7, y^8});
assertEmptyTerracini(2, f)
assertEmptyTerracini(3, f)

-- example 4.5
kk = ZZ/32003;
R = kk[x, y];
S = kk[z_0..z_4];
f = map(R, S, {x^5, x^4*y, x^3*y^2, x*y^4, y^5});
assertEmptyTerracini(2, f)

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

------------------------
-- Veronese varieties --
------------------------

needsPackage "Resultants"

assertEmptyTerracini(2, veronese(2, 3))
I = elapsedTime terraciniLocus(3, veronese(2, 3));
assert(#primaryDecomposition I == 1 and dim I - 3 == 2 * 2 + 3 - 2)

------------------------------
-- Segre-Veronese varieties --
------------------------------

segreVeronese = (n, d) -> (
    x := symbol x;
    r := #n;
    R := QQ new Array from splice apply(r, i -> x_(i, 0)..x_(i, n#i));
    y := symbol y;
    S := QQ[y_0..y_(product(n, d, (ni, di) -> binomial(ni + di, di)) - 1)];
    map(R, S, flatten entries first tensor apply(r, i -> (
		vector apply(subsets(n#i + d#i, d#i), A -> product(d#i, j ->
			x_(i, A#j - j)))))))

assertTheorem76 = (n, d) -> (
    r := ceiling((d#0 + 2)/2);
    I := elapsedTime terraciniLocus(r, segreVeronese(n, d));
    jHat := max select(#d, i -> ceiling((d#i + 2)/2) == r) + 1;
    comps := primaryDecomposition I;
    assert(#comps == jHat);
    assert(sort apply(comps, J -> dim J - 2 - r) ==
	sort apply (jHat, i -> sum n + n#i + r - 2)))

assertTheorem76({1, 1}, {1, 2})
assertTheorem76({1, 1}, {1, 3})
assertTheorem76({1, 1}, {2, 2})
