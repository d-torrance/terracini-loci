-- Copyright (c) 2023 Francesco Galuppi, Pierpaola Santarsiero, Doug Torrance,
-- and Ettore Teixeira Turatti

-- Permission is hereby granted, free of charge, to any person obtaining
-- a copy of this software and associated documentation files (the
-- "Software"), to deal in the Software without restriction, including
-- without limitation the rights to use, copy, modify, merge, publish,
-- distribute, sublicense, and/or sell copies of the Software, and to
-- permit persons to whom the Software is furnished to do so, subject to
-- the following conditions:

-- The above copyright notice and this permission notice (including the
-- next paragraph) shall be included in all copies or substantial
-- portions of the Software.

-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
-- EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
-- MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
-- NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
-- LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
-- OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
-- WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

newPackage("TerraciniLoci",
    Headline => "Terracini loci of projective varieties",
    Version => "0.1",
    Date => "TBD",
    Authors => {
	{
	    Name => "Francesco Galuppi",
	    Email => "francesco.galuppi@impan.pl"},
	{
	    Name => "Pierpaola Santarsiero",
	    Email => "pierpaola.santarsiero@uni-osnabrueck.de"},
	{
	    Name => "Doug Torrance",
	    Email => "dtorrance@piedmont.edu",
	    HomePage => "https://webwork.piedmont.edu/~dtorrance"},
	{
	    Name => "Ettore Teixeira Turatti",
	    Email => "ettore.t.turatti@uit.no"}},
    HomePage => "https://github.com/d-torrance/terracini-loci",
    Keywords => {"Projective Algebraic Geometry"},
    PackageImports => {
	"CorrespondenceScrolls",
	"FastMinors",
	"MinimalPrimes"})

export {
    "terraciniLocus"
    }

importFrom("Core", {"concatRows"})

terraciniLocus = method()

terraciniLocus(ZZ, Matrix, Ideal) := (r, A, I) -> (
    if ring A =!= ring I then error "expected rings to agree";
    R := ring A;
    s := numRows A;
    t := numColumns A;
    rk := rank A;
    n := numgens R - 1;
    Q := productOfProjectiveSpaces(toList(r : n),
	CoefficientField => coefficientRing R,
	VariableName => "z");
    opts := apply(r, i -> apply(n + 1, j -> R_j => Q_((n + 1) * i + j)));
    Az := concatRows apply(r, i -> sub(A, opts#i));
    Ir := ideal apply(r, i -> sub(I, opts#i));
    result := recursiveMinors(min(r * rk, t), Az, Threads => 4) + Ir;
    Z := genericMatrix(Q, n + 1, r);
    duplicate := intersect apply(subsets(r, 2), ij ->
	recursiveMinors(2, Z_ij, Threads => 4));
    result = saturate(result, duplicate);
    blocksingular := recursiveMinors(rk, A);
    singular := intersect apply(r, i -> sub(blocksingular, opts#i));
    radical result : radical singular)

terraciniLocus(ZZ, RingMap) := (r, f) -> (
    terraciniLocus(r, jacobian matrix f, ideal 0_(target f)))

terraciniLocus(ZZ, Ideal) := (r, I) -> (
    terraciniLocus(r, transpose jacobian I, I))

beginDocumentation()

doc ///
  Key
    TerraciniLoci
  Headline
    package for computing Terracini loci
  Description
    Text
      This package implements the algorithms from Section 8 of the
      paper "Geometry of first nonempty Terracini loci" by F. Galuppi,
      P. Santarsiero, D. Torrance, and E. Turatti.

      The Terracini locus of projective variety $X$ is a subvariety of
      the symmetric power $X^{(r)}$ containing all tclosure of all
      sets $\{p_1,\ldots,p_r\}$ of points in $X$ for which the space
      $\langle T_{p_1}X,\ldots,T_{p_r}X\rangle$ has less than the expected
      dimension.

      This package exports one method, @TO terraciniLocus@, for computing the
      ideals of these varieties.
///

doc ///
  Key
    terraciniLocus
    (terraciniLocus, ZZ, Matrix, Ideal)
    (terraciniLocus, ZZ, RingMap)
    (terraciniLocus, ZZ, Ideal)
  Headline
    compute the Terracini locus of a projective variety
  Usage
    terraciniLocus(r, X)
  Inputs
    r:ZZ
    X:{RingMap,Ideal}
  Outputs
    :Ideal
  Description
    Text
      TODO
///

-----------
-- tests --
-----------
-- just the faster (< 1s) examples

TEST ///
-- rational normal curves
needsPackage "Resultants"
assertEmptyTerracini = (r, f) -> assert(terraciniLocus(r, f) == 1)

-- ring map
assertEmptyTerracini(2, veronese(1, 3))
assertEmptyTerracini(2, veronese(1, 4))

-- ideal (slower)
assertEmptyTerracini(2, ker veronese(1, 3))
///

TEST ///
-- del pezzo surfaces
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

assertCorollary55 = t -> (
    I := terraciniLocus(2, delPezzoSurface t);
    comps := primaryDecomposition I;
    assert(#comps == (if t == 4 then 5 else t) and
	all(comps, J -> dim J - 2 == 3)))

assertCorollary55 1
assertCorollary55 2
assertCorollary55 3
assertCorollary55 4
///

TEST ///
-- veronese
needsPackage "Resultants"
assertEmptyTerracini = (r, f) -> assert(terraciniLocus(r, f) == 1)

assertEmptyTerracini(2, veronese(2, 3))
///

TEST ///
-- segre-veronese

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
    I := terraciniLocus(r, segreVeronese(n, d));
    jHat := max select(#d, i -> ceiling((d#i + 2)/2) == r) + 1;
    comps := primaryDecomposition I;
    assert(#comps == jHat);
    assert(sort apply(comps, J -> dim J - 2 - r) ==
	sort apply (jHat, i -> sum n + n#i + r - 2)))

assertTheorem76({1, 1}, {1, 2})
///
