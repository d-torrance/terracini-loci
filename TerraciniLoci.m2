newPackage("TerraciniLoci",
    PackageImports => {"CorrespondenceScrolls", "FastMinors"})

export {
    "terraciniLocus"
    }

importFrom("Core", {"concatRows"})

terraciniLocus = method()

terraciniLocus(ZZ, Matrix, Ideal) := (r, A, I) -> (
    verboseLog := if debugLevel > 0 then printerr else identity;
    if ring A =!= ring I then error "expected rings to agree";
    R := ring A;
    nrows := numRows A;
    ncols := numColumns A;
    verboseLog("stacked Jacobian is ", toString(r * nrows), " x ",
	toString ncols);
    n := numgens R - 1;
    S := productOfProjectiveSpaces(toList(r : n),
	CoefficientField => coefficientRing R,
	VariableName => "z");
    opts := apply(r, i -> apply(n + 1, j -> R_j => S_((n + 1) * i + j)));
    B := concatRows apply(r, i -> sub(A, opts#i));
    J := ideal apply(r, i -> sub(I, opts#i));
    verboseLog("computing ",
	toString binomial(max(r * nrows, ncols), min(r * nrows, ncols)),
	" minors of stacked Jacobian ...");
    result := recursiveMinors(min(r * nrows, ncols), B, Threads => 4) + J;
    G := genericMatrix(S, n + 1, r);
    verboseLog("computing ",
	toString(binomial(r, 2) * binomial(n + 1, 2)),
	" minors for ideal of duplicate points");
    duplicate := intersect apply(subsets(r, 2), ij ->
	recursiveMinors(2, G_ij, Threads => 4));
    result = saturate(result, duplicate);
    blockrank := if zero codim I then nrows else codim I;
    verboseLog("computing ",
	toString(r * binomial(nrows, blockrank)),
	" minors for ideal of singular points");
    singular := intersect apply(r, i ->
	recursiveMinors(blockrank, B^(toList(i * nrows..(i + 1) * nrows - 1)),
	    Threads => 4));
    saturate(result, singular))

terraciniLocus(ZZ, RingMap) := (r, f) -> (
    terraciniLocus(r, jacobian matrix f, ideal 0_(target f)))

terraciniLocus(ZZ, Ideal) := (r, I) -> (
    terraciniLocus(r, transpose jacobian I, I))
