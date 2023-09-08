-- Copyright (c) 2023 Francesco Galuppi, Pierpaola Starsiero, Doug Torrance,
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
    verboseLog := if debugLevel > 0 then printerr else identity;
    if ring A =!= ring I then error "expected rings to agree";
    R := ring A;
    s := numRows A;
    t := numColumns A;
    rk := rank A;
    verboseLog("stacked matrix is ", toString(r * s), " x ", toString t);
    n := numgens R - 1;
    Q := productOfProjectiveSpaces(toList(r : n),
	CoefficientField => coefficientRing R,
	VariableName => "z");
    opts := apply(r, i -> apply(n + 1, j -> R_j => Q_((n + 1) * i + j)));
    Az := concatRows apply(r, i -> sub(A, opts#i));
    Ir := ideal apply(r, i -> sub(I, opts#i));
    verboseLog("computing ",
	toString binomial(max(r * rk, t), min(r * rk, t)),
	" minors of stacked matrix ...");
    result := recursiveMinors(min(r * rk, t), Az, Threads => 4) + Ir;
    Z := genericMatrix(Q, n + 1, r);
    verboseLog("computing ",
	toString(binomial(r, 2) * binomial(n + 1, 2)),
	" minors for ideal of duplicate points");
    duplicate := intersect apply(subsets(r, 2), ij ->
	recursiveMinors(2, Z_ij, Threads => 4));
    result = saturate(result, duplicate);
    verboseLog("computing ", toString binomial(s, rk),
	" minors for ideal of singular points");
    blocksingular := recursiveMinors(rk, A);
    singular := intersect apply(r, i -> sub(blocksingular, opts#i));
    radical result : radical singular)

terraciniLocus(ZZ, RingMap) := (r, f) -> (
    terraciniLocus(r, jacobian matrix f, ideal 0_(target f)))

terraciniLocus(ZZ, Ideal) := (r, I) -> (
    terraciniLocus(r, transpose jacobian I, I))
