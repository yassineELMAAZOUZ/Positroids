using Positroids

# Projection analysis uses the target/anti-excedance convention by default, so
# this is accepted directly as a rank-two decorated permutation.
p = [4, 3, 1, 5, 2]

# A map Gr(2,5) -> Gr(2,4) requires a 5-by-4 matrix. This is the natural
# completion of the matrix in the question: Z = [I_4; -1 1 -1 1].
Z = [ 1  0  0  0;
      0  1  0  0;
      0  0  1  0;
      0  0  0  1;
     -1  1 -1  1]

parent = projection_jacobian_report(p, Z)
println("Parent: dimension $(parent.dimension), Jacobian rank $(parent.jacobian_rank)")

println("\nFacets")
facets = boundary_projection_jacobian_report(p, Z; strata=:facets)
for result in facets
    println(rpad(string(result.permutation), 22),
            " dimension=", result.dimension,
            " rank=", result.jacobian_rank,
            " certificate=", result.certificate)
end

println("\nFull proper boundary")
boundary = boundary_projection_jacobian_report(p, Z; strata=:boundary)
for dimension in 0:(parent.dimension-1)
    level = filter(result -> result.dimension == dimension, boundary)
    full = count(result -> result.full_rank, level)
    println("dimension $dimension: $full / $(length(level)) have full generic rank")
end

poset = projection_boundary_poset(p, Z)
println("\nClosed boundary poset: $(length(poset.nodes)) nodes and " *
        "$(length(poset.covers)) cover relations")
println("f-vector: ", poset.f_vector)

# The unique rank-deficient facet has bridge chart
#
# [1  0   0   0   0]
# [0  1  alpha3  alpha2  alpha1].
#
# Its projected plane has Plucker coordinates
# [1+alpha1 : alpha3-alpha1 : alpha1+alpha2 : 0 : 0 : 0].
# Thus its image lies in P^2 and the three-dimensional facet has
# one-dimensional generic fibers. For lambda sufficiently close to 1,
#
# alpha1' = lambda*(1+alpha1)-1
# alpha2' = 1+lambda*(alpha2-1)
# alpha3' = lambda*(1+alpha3)-1
#
# stays positive and gives the same projective Plucker point.
bad_facet = only(filter(result -> !result.full_rank, facets))
println("\nRank-deficient facet: ", bad_facet.permutation)
println("Projected symbolic matrix:")
display(bridge_parametrization(bad_facet.source_permutation) * Z)
