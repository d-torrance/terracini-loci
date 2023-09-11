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
