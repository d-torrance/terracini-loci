newPackage("TerraciniLoci",
    PackageImports => {"CorrespondenceScrolls", "FastMinors"})

export {
    "terraciniLocus"
    }

importFrom("Core", {"concatRows"})

terraciniLocus = method()

terraciniLocus(ZZ, Matrix, Ideal) := (r, A, I) -> (
    if ring A =!= ring I then error "expected rings to agree";
    R := ring A;
    nrows := numRows A;
    ncols := numColumns A;
    n := numgens R - 1;
    S := productOfProjectiveSpaces(toList(r : n),
	CoefficientField => coefficientRing R,
	VariableName => "z");
    opts := apply(r, i -> apply(n + 1, j -> R_j => S_((n + 1) * i + j)));
    B := concatRows apply(r, i -> sub(A, opts#i));
    J := ideal apply(r, i -> sub(I, opts#i));
    result := recursiveMinors(min(r * nrows, ncols), B, Threads => 4) + J;
    G := genericMatrix(S, n + 1, r);
    blockrank := if zero codim I then nrows else codim I;
    duplicate := intersect apply(subsets(r, 2), ij ->
	recursiveMinors(2, G_ij, Threads => 4));
    result = saturate(result, duplicate);
    singular := intersect apply(r, i ->
	recursiveMinors(blockrank, B^(toList(i * nrows..(i + 1) * nrows - 1)),
	    Threads => 4));
    saturate(result, singular))

terraciniLocus(ZZ, RingMap) := (r, f) -> (
    terraciniLocus(r, jacobian matrix f, ideal 0_(target f)))

terraciniLocus(ZZ, Ideal) := (r, I) -> (
    terraciniLocus(r, transpose jacobian I, I))
