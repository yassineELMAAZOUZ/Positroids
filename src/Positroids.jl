module Positroids

using LinearAlgebra
import Plots
import HTTP
import JSON
import UUIDs

const PACKAGE_BUILD = v"1.3.0"

export decorated_permutations, decorated_excedances,
       shiftedOrder, cyclicInterval, iCompare, iCompareLists,
       minGrassmannNecklace, maxGrassmannNecklace, reverseGrassmannNecklace,
       fromNecklaceToPositroid, fromNecklaceToDecoratedPermutation,
       fromDecoratedPermToNecklace, fromDecoratedPermToMinNecklace,
       fromDecoratedPermToMaxNecklace, fromDecoratedPermToPositroid,
       fromPositroidToPermutation, countExceedences,
       alignmentNumber, dimensionOfPositroid, dimensionOfPermutation,
       le_diagram, le_diagram_matrix, draw_chords,
       PlabicGraph, plabic_graph, layout_plabic_graph,
       draw_plabic_graph, trip_permutation, is_bipartite,
       is_reduced, validate_plabic_graph, vertex_degrees, plabic_graph_from_drawing, PACKAGE_BUILD,
       square_cycles, square_move, square_move_by_label,
       PlabicEmbedding, plabic_embedding, graph_faces, graph_trips,
       face_labels, compare_plabic_graphs, InteractivePlabicGraph,
       interactive_plabic_graph, interactive_session, serve_positroids_web,
       is_orthogonal, orthogonal_dimension, orthogonal_dimensionOfPermutation,
       is_Omega_symplectic, is_Omega_symplecticTwo,
       symplectic_Omega_dimension, is_E_symplectic, symplectic_E_dimension,
       right_twist, left_twist, twist_map,
       plucker_coordinates, is_totally_nonnegative, is_totally_positive,
       BridgeParametrization, BoundaryMeasurementParametrization,
       bridge_parametrization, boundary_measurement_parametrization,
       parametrization_matrix, symbolic_matrix, symbolic_plucker_coordinates,
       symbolic_product, parameter_names,
       twistPerm, untwistPerm, twistNecklace,
       dimensionOfLGperm, is_cyclo_isotropic,
       permutationFromCycles, inversionsOfPerm, pminversionsOfPerm,
       is_codim1_face, is_child, randomMatching, up, iota, iota2, eta
export immediate_children

function _load_plotting()
    return nothing
end

"""
    PlabicGraph

A bicolored graph embedded in a disk. Boundary vertices are `1:n`; internal
vertices have color `:black` or `:white`. `positions` supplies a deterministic
disk layout, and `bridges` records the bridge decomposition used to construct
the graph.
"""
struct PlabicGraph
    n::Int
    colors::Vector{Symbol}
    edges::Vector{Tuple{Int,Int}}
    positions::Vector{NTuple{2,Float64}}
    permutation::Vector{Int}
    bridges::Vector{Tuple{Int,Int}}
end

"""A validated disk embedding with rotation, faces, and trips."""
struct PlabicEmbedding
    graph::PlabicGraph
    rotation::Vector{Vector{Int}}
    faces::Vector{Vector{Int}}
    halfedge_face::Dict{Tuple{Int,Int},Int}
    trips::Vector{Vector{Int}}
end

"""A live local-browser plabic graph session returned by `interactive_plabic_graph`."""
mutable struct InteractivePlabicGraph
    graph::Base.RefValue{PlabicGraph}
    server::Any
    url::String
end

Base.show(io::IO, session::InteractivePlabicGraph) =
    print(io,"InteractivePlabicGraph(\"",session.url,"\")")

abstract type AbstractCellParametrization end

"""A positive bridge-coordinate chart with one parameter per bridge."""
struct BridgeParametrization <: AbstractCellParametrization
    permutation::Vector{Int}
    gauge::Vector{Int}
    bridges::Vector{Tuple{Int,Int}}
    signs::Vector{Int}
    parameters::Vector{Symbol}
end

"""Boundary measurement of the normalized directed bridge network."""
struct BoundaryMeasurementParametrization <: AbstractCellParametrization
    permutation::Vector{Int}
    gauge::Vector{Int}
    network_bridges::Vector{Tuple{Int,Int}}
    signs::Vector{Int}
    edge_weights::Vector{Symbol}
end

function _inverse_decorated_permutation(p::Vector{Int})
    n = length(p)
    inverse = zeros(Int, n)
    for i in 1:n
        inverse[abs(p[i])] = i
    end
    for i in 1:n
        p[i] == -i && (inverse[i] = -i)
    end
    return inverse
end

function _parametrization_data(input)
    p = Int.(input)
    _validate_decorated_permutation(p)
    inverse = _inverse_decorated_permutation(p)
    base, removal = _bridge_decomposition(inverse)
    n = length(p)
    gauge = sort([base[i]-n for i in 1:n if base[i] > n])
    bridges = reverse(removal)
    signs = [(-1)^count(x -> a < x < b, gauge) for (a, b) in bridges]
    length(gauge) == decorated_excedances(p) || error("bridge gauge has incorrect rank")
    return p, gauge, bridges, signs
end

"""
    bridge_parametrization(p)

Construct the positive bridge chart of the positroid cell indexed by `p`.
There is one coordinate `αᵢ` per bridge and hence per cell dimension.
"""
function bridge_parametrization(p::AbstractVector{<:Integer})
    p, gauge, bridges, signs = _parametrization_data(p)
    parameters = [Symbol("α", i) for i in eachindex(bridges)]
    return BridgeParametrization(p, gauge, bridges, signs, parameters)
end

"""
    boundary_measurement_parametrization(p)

Construct the boundary-measurement chart of the normalized directed bridge
network for `p`. Its variables `wᵢ` are the nontrivial bridge-edge weights;
all gauge-tree edge weights are normalized to one.
"""
function boundary_measurement_parametrization(p::AbstractVector{<:Integer})
    p, gauge, bridges, signs = _parametrization_data(p)
    weights = [Symbol("w", i) for i in eachindex(bridges)]
    return BoundaryMeasurementParametrization(p, gauge, bridges, signs, weights)
end

parameter_names(P::BridgeParametrization) = copy(P.parameters)
parameter_names(P::BoundaryMeasurementParametrization) = copy(P.edge_weights)
_chart_bridges(P::BridgeParametrization) = P.bridges
_chart_bridges(P::BoundaryMeasurementParametrization) = P.network_bridges

function parametrization_matrix(P::AbstractCellParametrization,
                                values::AbstractVector=ones(length(parameter_names(P))))
    length(values) == length(parameter_names(P)) ||
        throw(ArgumentError("expected $(length(parameter_names(P))) parameter values"))
    T = isempty(values) ? Float64 : promote_type(map(typeof, values)...)
    k, n = length(P.gauge), length(P.permutation)
    A = zeros(T, k, n)
    for (row, column) in enumerate(P.gauge)
        A[row, column] = one(T)
    end
    for (j, (a, b)) in enumerate(_chart_bridges(P))
        A[:, b] .+= P.signs[j] * values[j] .* A[:, a]
    end
    return A
end

(P::AbstractCellParametrization)(values::AbstractVector) = parametrization_matrix(P, values)

function _symbolic_add(a::String, b::String)
    a == "0" && return b
    b == "0" && return a
    return "(" * a * " + " * b * ")"
end

function _symbolic_term(sign, variable, expression)
    expression == "0" && return "0"
    factor = string(variable) * (expression == "1" ? "" : "*" * expression)
    return sign == 1 ? factor : "-" * factor
end

"""Return a readable symbolic matrix for a cell parametrization."""
function symbolic_matrix(P::AbstractCellParametrization)
    k, n = length(P.gauge), length(P.permutation)
    A = fill("0", k, n)
    for (row, column) in enumerate(P.gauge)
        A[row, column] = "1"
    end
    names = parameter_names(P)
    for (j, (a, b)) in enumerate(_chart_bridges(P)), row in 1:k
        A[row, b] = _symbolic_add(A[row, b],
                                  _symbolic_term(P.signs[j], names[j], A[row, a]))
    end
    return A
end

const _ParamPoly = Dict{Tuple,Int}
_poly_zero() = _ParamPoly()
_poly_constant(c::Int) = c == 0 ? _poly_zero() : _ParamPoly(() => c)
_poly_variable(i::Int) = _ParamPoly((i,) => 1)

function _poly_add(a::_ParamPoly, b::_ParamPoly)
    result = copy(a)
    for (monomial, coefficient) in b
        result[monomial] = get(result, monomial, 0) + coefficient
        result[monomial] == 0 && delete!(result, monomial)
    end
    return result
end

_poly_negate(a::_ParamPoly) = _ParamPoly(m => -c for (m, c) in a)

function _poly_multiply(a::_ParamPoly, b::_ParamPoly)
    result = _poly_zero()
    for (ma, ca) in a, (mb, cb) in b
        monomial = Tuple(sort([ma...; mb...]))
        result[monomial] = get(result, monomial, 0) + ca*cb
        result[monomial] == 0 && delete!(result, monomial)
    end
    return result
end

function _polynomial_matrix(P::AbstractCellParametrization)
    k, n = length(P.gauge), length(P.permutation)
    A = [_poly_zero() for _ in 1:k, _ in 1:n]
    for (row, column) in enumerate(P.gauge)
        A[row, column] = _poly_constant(1)
    end
    for (j, (a, b)) in enumerate(_chart_bridges(P)), row in 1:k
        term = _poly_multiply(_poly_constant(P.signs[j]), _poly_variable(j))
        term = _poly_multiply(term, A[row, a])
        A[row, b] = _poly_add(A[row, b], term)
    end
    return A
end

function _polynomial_determinant(A)
    n = size(A, 1)
    n == size(A, 2) || throw(ArgumentError("determinant requires a square matrix"))
    n == 0 && return _poly_constant(1)
    n == 1 && return copy(A[1, 1])
    result = _poly_zero()
    for column in 1:n
        rows = 2:n
        columns = [j for j in 1:n if j != column]
        term = _poly_multiply(A[1, column], _polynomial_determinant(A[rows, columns]))
        isodd(1 + column) && (term = _poly_negate(term))
        result = _poly_add(result, term)
    end
    return result
end

function _format_polynomial(poly::_ParamPoly, names)
    isempty(poly) && return "0"
    monomials = sort(collect(keys(poly)); by=m -> (length(m), m))
    pieces = String[]
    for monomial in monomials
        coefficient = poly[monomial]
        variables = isempty(monomial) ? "" : join(string.(names[collect(monomial)]), "*")
        magnitude = abs(coefficient)
        body = isempty(variables) ? string(magnitude) :
               magnitude == 1 ? variables : string(magnitude, "*", variables)
        if isempty(pieces)
            push!(pieces, coefficient < 0 ? "-" * body : body)
        else
            push!(pieces, coefficient < 0 ? " - " * body : " + " * body)
        end
    end
    return join(pieces)
end

"""
    symbolic_plucker_coordinates(P)
    plucker_coordinates(P)

Return all maximal minors of a cell parametrization as simplified multivariate
polynomial strings indexed by increasing column tuples.
"""
function symbolic_plucker_coordinates(P::AbstractCellParametrization)
    A = _polynomial_matrix(P)
    k, n = size(A)
    names = parameter_names(P)
    return Dict(Tuple(I) => _format_polynomial(_polynomial_determinant(A[:, I]), names)
                for I in _subsets(Vector(1:n), k))
end

plucker_coordinates(P::AbstractCellParametrization) = symbolic_plucker_coordinates(P)

"""
    symbolic_product(P, Z)
    P * Z

Multiply a symbolic cell-parametrization matrix by an integer matrix `Z` and
return the simplified polynomial entries as strings.
"""
function symbolic_product(P::AbstractCellParametrization,
                          Z::AbstractMatrix{<:Integer})
    A = _polynomial_matrix(P)
    size(A, 2) == size(Z, 1) ||
        throw(DimensionMismatch("parametrization has $(size(A,2)) columns but Z has $(size(Z,1)) rows"))
    result = [_poly_zero() for _ in 1:size(A,1), _ in 1:size(Z,2)]
    for i in axes(A,1), j in axes(Z,2), l in axes(A,2)
        term = _poly_multiply(A[i,l], _poly_constant(Int(Z[l,j])))
        result[i,j] = _poly_add(result[i,j], term)
    end
    names = parameter_names(P)
    return [_format_polynomial(result[i,j], names)
            for i in axes(result,1), j in axes(result,2)]
end

Base.:*(P::AbstractCellParametrization, Z::AbstractMatrix{<:Integer}) =
    symbolic_product(P, Z)

function _bounded_affine_permutation(p::Vector{Int})
    n = length(p)
    [_affine_value_from_decorated_perm(p, i) for i in 1:n]
end

function _bridge_decomposition(p::Vector{Int})
    q = _bounded_affine_permutation(p)
    n = length(q)
    bridges = Tuple{Int,Int}[]

    # This is the default bridge chain used by the historical Mathematica
    # package: repeatedly swap the first adjacent active ascent.
    limit = n * n + sum(q .- collect(1:n))
    steps = 0
    while true
        fixed = [mod(q[i] - 1, n) + 1 == i for i in 1:n]
        active = findall(!, fixed)
        pair_index = findfirst(t -> q[t[1]] < q[t[2]],
                               [(active[j], active[j + 1]) for j in 1:max(0, length(active)-1)])
        isnothing(pair_index) && break
        a, b = active[pair_index], active[pair_index + 1]
        push!(bridges, (a, b))
        q[a], q[b] = q[b], q[a]
        steps += 1
        steps <= limit || error("bridge decomposition did not terminate")
    end
    return q, bridges
end

function _replace_boundary_endpoint!(edges, boundary::Int, replacement::Int)
    index = findfirst(e -> e[1] == boundary || e[2] == boundary, edges)
    isnothing(index) && error("boundary vertex $boundary has no incident edge")
    u, v = edges[index]
    edges[index] = u == boundary ? (replacement, v) : (u, replacement)
end

function _contract_same_color_edges(colors, edges, positions, n)
    total = length(colors)
    parent = collect(1:total)
    function root(v)
        while parent[v] != v
            parent[v] = parent[parent[v]]
            v = parent[v]
        end
        return v
    end
    function unite(u, v)
        ru, rv = root(u), root(v)
        ru == rv && return
        parent[max(ru, rv)] = min(ru, rv)
    end
    for (u, v) in edges
        if u > n && v > n && colors[u] == colors[v]
            unite(u, v)
        end
    end

    components = Dict{Int,Vector{Int}}()
    for v in (n+1):total
        push!(get!(components, root(v), Int[]), v)
    end
    roots = sort(collect(keys(components)))
    mapping = collect(1:total)
    new_colors = copy(colors[1:n])
    new_positions = copy(positions[1:n])
    for r in roots
        members = components[r]
        new_id = length(new_colors) + 1
        for v in members
            mapping[v] = new_id
        end
        push!(new_colors, colors[members[1]])
        push!(new_positions, (sum(positions[v][1] for v in members)/length(members),
                              sum(positions[v][2] for v in members)/length(members)))
    end
    new_edges = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    for (u, v) in edges
        a, b = mapping[u], mapping[v]
        a == b && continue
        edge = a < b ? (a, b) : (b, a)
        edge in seen && continue
        push!(seen, edge); push!(new_edges, edge)
    end
    return new_colors, new_edges, new_positions
end

function _contract_one_same_color_edge(colors, edges, positions, n)
    edge_index = findfirst(edges) do (u, v)
        u > n && v > n && colors[u] == colors[v]
    end
    isnothing(edge_index) && return colors, edges, positions, false, nothing
    u, v = edges[edge_index]
    kept, removed = minmax(u, v)
    mapping = zeros(Int, length(colors))
    new_colors = Symbol[]
    new_positions = NTuple{2,Float64}[]
    midpoint = ((positions[kept][1] + positions[removed][1]) / 2,
                (positions[kept][2] + positions[removed][2]) / 2)
    for old in eachindex(colors)
        old == removed && continue
        mapping[old] = length(new_colors) + 1
        push!(new_colors, colors[old])
        push!(new_positions, old == kept ? midpoint : positions[old])
    end
    mapping[removed] = mapping[kept]
    new_edges = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    for (a, b) in edges
        x, y = mapping[a], mapping[b]
        x == y && continue
        edge = x < y ? (x, y) : (y, x)
        edge in seen || (push!(seen, edge); push!(new_edges, edge))
    end
    return new_colors, new_edges, new_positions, true, (u, v)
end

function _delete_internal_vertex(colors, edges, positions, n, removed)
    mapping = zeros(Int, length(colors))
    new_colors = Symbol[]
    new_positions = NTuple{2,Float64}[]
    for old in eachindex(colors)
        old == removed && continue
        mapping[old] = length(new_colors) + 1
        push!(new_colors, colors[old])
        push!(new_positions, positions[old])
    end
    new_edges = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    for (u, v) in edges
        (u == removed || v == removed) && continue
        a, b = mapping[u], mapping[v]
        a == b && continue
        edge = a < b ? (a, b) : (b, a)
        edge in seen || (push!(seen, edge); push!(new_edges, edge))
    end
    return new_colors, new_edges, new_positions, mapping
end

function _suppress_one_bivalent(colors, edges, positions, n)
    degrees = zeros(Int, length(colors))
    neighbors = [Int[] for _ in colors]
    for (u, v) in edges
        degrees[u] += 1; degrees[v] += 1
        push!(neighbors[u], v); push!(neighbors[v], u)
    end
    vertex = findfirst(v -> degrees[v] == 2, (n+1):length(colors))
    isnothing(vertex) && return colors, edges, positions, false
    vertex = n + vertex
    a, b = neighbors[vertex]
    colors, edges, positions, mapping =
        _delete_internal_vertex(colors, edges, positions, n, vertex)
    a, b = mapping[a], mapping[b]
    if a != b
        edge = a < b ? (a, b) : (b, a)
        edge in edges || push!(edges, edge)
    end
    return colors, edges, positions, true
end

function _reduce_bipartite_graph(colors, edges, positions, n)
    while true
        old_vertices, old_edges = length(colors), length(edges)
        colors, edges, positions = _contract_same_color_edges(colors, edges, positions, n)
        while true
            colors, edges, positions, changed =
                _suppress_one_bivalent(colors, edges, positions, n)
            changed || break
        end
        colors, edges, positions = _contract_same_color_edges(colors, edges, positions, n)
        length(colors) == old_vertices && length(edges) == old_edges && break
    end
    return colors, edges, positions
end

function _drawing_reduction_animation(n, internal_colors, input_edges,
                                      internal_positions)
    # Compute and validate the mathematical result first.  The copies below are
    # then reduced independently so that every local operation can be retained
    # as a frame for the browser instead of disappearing inside the reducer.
    result = plabic_graph_from_drawing(n, internal_colors, input_edges,
                                       internal_positions; reduce=true)
    n = Int(n)
    colors = vcat(fill(:boundary, n), Symbol.(internal_colors))
    positions = NTuple{2,Float64}[]
    for i in 1:n
        theta = -2pi * (i - 1) / n - pi / 2
        push!(positions, (cos(theta), sin(theta)))
    end
    append!(positions, [(Float64(p[1]), Float64(p[2])) for p in internal_positions])
    edges = unique([(min(Int(e[1]), Int(e[2])), max(Int(e[1]), Int(e[2])))
                    for e in input_edges])
    stages = Tuple{String,PlabicGraph}[]
    snapshot(caption) = push!(stages, (caption,
        PlabicGraph(n, copy(colors), copy(edges), copy(positions),
                    copy(result.permutation), Tuple{Int,Int}[])))

    step = 0
    while true
        new_colors, new_edges, new_positions, changed, contracted =
            _contract_one_same_color_edge(colors, edges, positions, n)
        if changed
            colors, edges, positions = new_colors, new_edges, new_positions
            step += 1
            u, v = contracted
            snapshot("Contract same-color edge $(u)–$(v) · step $step")
            continue
        end
        new_colors, new_edges, new_positions, changed =
            _suppress_one_bivalent(colors, edges, positions, n)
        if changed
            colors, edges, positions = new_colors, new_edges, new_positions
            step += 1
            snapshot("Suppress one bivalent vertex · step $step")
            continue
        end
        break
    end
    isempty(stages) && snapshot("No local contractions or bivalent vertices")
    local_graph = PlabicGraph(n, colors, edges, positions,
                              collect(1:n), Tuple{Int,Int}[])
    local_graph = _with_recomputed_trip_permutation(local_graph)
    reduced_graph, move_stages = _reduce_by_explicit_moves(local_graph; animate=true)
    append!(stages, move_stages)
    return reduced_graph, stages
end

"""Return `true` when every internal edge joins opposite vertex colors."""
function is_bipartite(G::PlabicGraph)
    all(G.colors[u] == :boundary || G.colors[v] == :boundary ||
        G.colors[u] != G.colors[v] for (u, v) in G.edges)
end

"""
    is_reduced(G)

Check the reduced bridge certificate for a graph produced by `plabic_graph`:
bipartiteness, univalent boundary vertices, no bivalent internal vertices, equality of
bridge count with positroid dimension, and equality of graph and cell dimension.
"""
function is_reduced(G::PlabicGraph)
    is_bipartite(G) || return false
    degrees = _plabic_degrees(G)
    all(degrees[1:G.n] .== 1) || return false
    all(d -> d > 0 && d != 2, degrees[(G.n+1):end]) || return false
    k = decorated_excedances(G.permutation)
    cell_dimension = dimensionOfPermutation(k, G.n, G.permutation)
    graph_dimension = length(G.edges) - (length(G.colors) - G.n)
    bridge_certificate = isempty(G.bridges) || length(G.bridges) == cell_dimension
    return bridge_certificate && cell_dimension == graph_dimension
end

"""Validate structural plabic-graph invariants, throwing an error on failure."""
function validate_plabic_graph(G::PlabicGraph; require_reduced=true)
    is_bipartite(G) || error("plabic graph contains an internal same-color edge")
    degrees = _plabic_degrees(G)
    all(degrees[1:G.n] .== 1) || error("every boundary vertex must be univalent")
    if require_reduced && !is_reduced(G)
        error("graph does not satisfy the reduced bridge/dimension certificate")
    end
    return true
end

function _relax_disk_layout!(positions, edges, n; iterations=1400,
                             repulsion=0.0014, attraction=0.10,
                             ideal_length=0.21, lollipop_length=0.12,
                             gravity=0.0025)
    count = length(positions)
    x = [p[1] for p in positions]
    y = [p[2] for p in positions]
    fx = zeros(Float64, count)
    fy = zeros(Float64, count)
    internal = (n + 1):count
    isempty(internal) && return positions
    degrees = zeros(Int, count)
    for (u, v) in edges
        degrees[u] += 1; degrees[v] += 1
    end
    charge(v) = v > n && degrees[v] == 1 ? 0.22 : 1.0

    for iteration in 1:iterations
        fill!(fx, 0); fill!(fy, 0)
        # Coulomb repulsion prevents vertices and unrelated edges from piling up.
        for u in internal, v in (u + 1):count
            dx, dy = x[u] - x[v], y[u] - y[v]
            d2 = dx*dx + dy*dy + 1e-4
            force = repulsion * charge(u) * charge(v) / d2
            fx[u] += force*dx; fy[u] += force*dy
            fx[v] -= force*dx; fy[v] -= force*dy
        end
        for u in internal, v in 1:n
            dx, dy = x[u] - x[v], y[u] - y[v]
            d2 = dx*dx + dy*dy + 1e-4
            force = 0.20repulsion * charge(u) / d2
            fx[u] += force*dx; fy[u] += force*dy
        end
        # Edge springs and a weak pull toward the disk center.
        for (u, v) in edges
            dx, dy = x[v] - x[u], y[v] - y[u]
            distance = sqrt(dx*dx + dy*dy) + 1e-8
            target = (degrees[u] == 1 || degrees[v] == 1) ? lollipop_length : ideal_length
            force = attraction * (distance - target)
            if u > n
                fx[u] += force*dx/distance; fy[u] += force*dy/distance
            end
            if v > n
                fx[v] -= force*dx/distance; fy[v] -= force*dy/distance
            end
        end
        for u in internal
            fx[u] -= gravity*x[u]; fy[u] -= gravity*y[u]
            step = 0.22 * (1 - 0.88iteration/iterations)
            x[u] += step * clamp(fx[u], -0.08, 0.08)
            y[u] += step * clamp(fy[u], -0.08, 0.08)
            r = hypot(x[u], y[u])
            if r > 0.94
                x[u] *= 0.94/r; y[u] *= 0.94/r
            end
        end
    end
    for u in internal
        positions[u] = (x[u], y[u])
    end
    return positions
end

function _proper_edge_crossing(p1, q1, p2, q2)
    o1 = orientation(p1, q1, p2); o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1); o4 = orientation(p2, q2, q1)
    return o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4
end

function _layout_score(positions, edges, n)
    crossings = 0
    for i in 1:(length(edges)-1), j in (i+1):length(edges)
        u, v = edges[i]; s, t = edges[j]
        length(Set((u, v, s, t))) < 4 && continue
        crossings += _proper_edge_crossing(positions[u], positions[v], positions[s], positions[t])
    end
    lengths = [hypot(positions[u][1]-positions[v][1], positions[u][2]-positions[v][2])
               for (u, v) in edges]
    crowding = 0.0
    for u in (n+1):length(positions), v in (u+1):length(positions)
        d = hypot(positions[u][1]-positions[v][1], positions[u][2]-positions[v][2])
        crowding += max(0.0, 0.11-d)^2
    end
    mean_length = isempty(lengths) ? 0.0 : sum(lengths)/length(lengths)
    variance = sum((lengths .- mean_length).^2)
    return 1000crossings + 300crowding + variance
end

function _layout_start(G::PlabicGraph, restart::Int)
    restart == 1 && return copy(G.positions)
    positions = copy(G.positions)
    phi = (sqrt(5)-1)/2
    total = max(1, length(positions)-G.n)
    for (j, u) in enumerate((G.n+1):length(positions))
        theta = 2pi * mod(j*phi + 0.173restart, 1)
        radius = 0.18 + 0.68sqrt(j/(total+1))
        positions[u] = (radius*cos(theta), radius*sin(theta))
    end
    return positions
end

function _harmonic_disk_layout(G::PlabicGraph;lollipop_length=0.14)
    positions=copy(G.positions)
    degrees=_plabic_degrees(G)
    adjacency=[Int[] for _ in G.colors]
    for (u,v) in G.edges
        push!(adjacency[u],v);push!(adjacency[v],u)
    end

    # Boundary vertices and genuine lollipops are fixed anchors.  Without the
    # second class, a harmonic average places a lollipop exactly on top of its
    # boundary vertex.
    fixed=Set(1:G.n)
    for v in (G.n+1):length(G.colors)
        degrees[v]==1 || continue
        push!(fixed,v)
        neighbor=only(adjacency[v])
        if neighbor<=G.n
            x,y=positions[neighbor]
            radius=max(0.0,1-lollipop_length)
            norm=hypot(x,y)
            positions[v]=norm>0 ? (radius*x/norm,radius*y/norm) : (radius,0.0)
        end
    end

    unknown=[v for v in (G.n+1):length(G.colors) if !(v in fixed)]
    isempty(unknown) && return positions
    index=Dict(v=>i for (i,v) in enumerate(unknown))
    laplacian=zeros(Float64,length(unknown),length(unknown))
    bx=zeros(Float64,length(unknown));by=zeros(Float64,length(unknown))
    for v in unknown
        row=index[v]
        laplacian[row,row]=length(adjacency[v])
        for neighbor in adjacency[v]
            if haskey(index,neighbor)
                laplacian[row,index[neighbor]]-=1
            else
                bx[row]+=positions[neighbor][1]
                by[row]+=positions[neighbor][2]
            end
        end
    end
    xs=laplacian\bx;ys=laplacian\by
    for (j,v) in enumerate(unknown)
        positions[v]=(xs[j],ys[j])
    end
    return positions
end

"""
    layout_plabic_graph(G; method=:harmonic, iterations=1400, restarts=6, ...)

Improve the disk embedding.  The default `:harmonic` method places every
non-lollipop internal vertex at the barycenter of its neighbors, spreading the
graph through the disk while preserving the fixed circular boundary.  Use
`method=:force` for the earlier constrained charge model.
"""
function layout_plabic_graph(G::PlabicGraph; method=:harmonic,
                             iterations=1400, restarts=6,
                             repulsion=0.0014, attraction=0.10,
                             ideal_length=0.21, lollipop_length=0.12,
                             gravity=0.0025)
    iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
    restarts >= 1 || throw(ArgumentError("restarts must be positive"))
    method in (:harmonic,:force) ||
        throw(ArgumentError("method must be :harmonic or :force"))
    if method==:harmonic
        positions=_harmonic_disk_layout(G;lollipop_length=lollipop_length)
        _edge_crossing_count(positions,G.edges)==0 ||
            return layout_plabic_graph(G;method=:force,iterations=iterations,
                restarts=restarts,repulsion=repulsion,attraction=attraction,
                ideal_length=ideal_length,lollipop_length=lollipop_length,
                gravity=gravity)
        return PlabicGraph(G.n,copy(G.colors),copy(G.edges),positions,
                           copy(G.permutation),copy(G.bridges))
    end
    best = copy(G.positions)
    best_score = Inf
    for restart in 1:restarts
        candidate = _layout_start(G, restart)
        _relax_disk_layout!(candidate, G.edges, G.n; iterations=iterations,
                            repulsion=repulsion, attraction=attraction,
                            ideal_length=ideal_length,
                            lollipop_length=lollipop_length, gravity=gravity)
        score = _layout_score(candidate, G.edges, G.n)
        if score < best_score
            best, best_score = candidate, score
        end
        best_score < 1 && break
    end
    return PlabicGraph(G.n, copy(G.colors), copy(G.edges), best,
                       copy(G.permutation), copy(G.bridges))
end

function _edge_crossing_count(positions, edges)
    crossings = 0
    for i in 1:(length(edges)-1), j in (i+1):length(edges)
        u, v = edges[i]; s, t = edges[j]
        length(Set((u,v,s,t))) < 4 && continue
        crossings += _proper_edge_crossing(positions[u], positions[v],
                                           positions[s], positions[t])
    end
    return crossings
end

function _rotation_from_positions(G::PlabicGraph, augmented_edges)
    adjacency = [Int[] for _ in G.colors]
    for (u,v) in augmented_edges
        push!(adjacency[u],v); push!(adjacency[v],u)
    end
    for v in eachindex(adjacency)
        unique!(adjacency[v])
        sort!(adjacency[v]; by=u -> atan(G.positions[u][2]-G.positions[v][2],
                                         G.positions[u][1]-G.positions[v][1]))
    end
    # Boundary labels run clockwise. In the disk rotation system the unique
    # graph edge lies between the preceding and following boundary arcs,
    # irrespective of small tangential distortions in the force layout.
    for i in 1:G.n
        previous,next = mod1(i-1,G.n),mod1(i+1,G.n)
        # Read the graph half-edge from G.edges itself.  An actual graph edge
        # may join adjacent boundary vertices (for example the reduced graph
        # for [2,1]), in which case filtering the augmented adjacency by the
        # endpoint number would incorrectly discard it as a boundary arc.
        graph_neighbors = Int[]
        for (u,v) in G.edges
            u==i && push!(graph_neighbors,v)
            v==i && push!(graph_neighbors,u)
        end
        unique!(graph_neighbors)
        length(graph_neighbors)==1 || error("boundary vertex $i must have one graph neighbor")
        adjacency[i] = [previous,only(graph_neighbors),next]
    end
    return adjacency
end

function _left_successor(rotation, u, v)
    neighbors = rotation[v]
    index = findfirst(==(u), neighbors)
    isnothing(index) && error("rotation system is missing half-edge ($u,$v)")
    return neighbors[mod1(index-1, length(neighbors))]
end

function _enumerate_faces(rotation, edges)
    halfedge_face = Dict{Tuple{Int,Int},Int}()
    faces = Vector{Vector{Int}}()
    for (a,b) in edges, start in ((a,b),(b,a))
        haskey(halfedge_face,start) && continue
        face = Int[]
        u,v = start
        face_id = length(faces)+1
        while true
            haskey(halfedge_face,(u,v)) && (u,v) != start &&
                error("invalid rotation system while tracing faces")
            halfedge_face[(u,v)] = face_id
            push!(face,u)
            w = _left_successor(rotation,u,v)
            u,v = v,w
            (u,v) == start && break
            length(face) <= 2length(edges)+2 || error("face traversal did not close")
        end
        push!(faces,face)
    end
    return faces,halfedge_face
end

function _signed_polygon_area(face, positions)
    points = positions[face]
    return 0.5sum(points[j][1]*points[mod1(j+1,length(points))][2] -
                  points[j][2]*points[mod1(j+1,length(points))][1]
                  for j in eachindex(points))
end

function _trace_trip(G, rotation, start)
    boundary_neighbor = only([v for (u,v) in G.edges if u==start] ∪
                             [u for (u,v) in G.edges if v==start])
    path = [start,boundary_neighbor]
    previous,current = start,boundary_neighbor
    limit = 2length(G.edges)+2
    while current > G.n
        neighbors = rotation[current]
        index = findfirst(==(previous),neighbors)
        isnothing(index) && error("trip entered a missing half-edge")
        step = G.colors[current] == :white ? -1 : 1
        next = neighbors[mod1(index+step,length(neighbors))]
        push!(path,next)
        previous,current = current,next
        length(path) <= limit || error("trip from boundary $start did not terminate")
    end
    return path
end

function _decorated_trip_permutation(G::PlabicGraph)
    augmented=unique(vcat(G.edges,[(i,mod1(i+1,G.n)) for i in 1:G.n]))
    rotation=_rotation_from_positions(G,augmented)
    trips=[_trace_trip(G,rotation,i) for i in 1:G.n]
    endpoints=last.(trips)
    sort(endpoints)==collect(1:G.n) ||
        throw(ArgumentError("the trips do not define a permutation"))
    permutation=copy(endpoints)
    for i in 1:G.n
        endpoints[i]==i || continue
        neighbor=only([v for (u,v) in G.edges if u==i] ∪
                      [u for (u,v) in G.edges if v==i])
        permutation[i]=G.colors[neighbor]==:white ? -i : i
    end
    _validate_decorated_permutation(permutation)
    return permutation
end

function _with_recomputed_trip_permutation(G::PlabicGraph)
    permutation=_decorated_trip_permutation(G)
    return PlabicGraph(G.n,copy(G.colors),copy(G.edges),copy(G.positions),
                       permutation,copy(G.bridges))
end

"""
    plabic_embedding(G; iterations=2500, restarts=24)

Construct and validate a planar disk embedding, its counterclockwise rotation
system, disk faces, and trips. The selected layout must be crossing-free and the
computed trip endpoints must match the graph's decorated permutation.
"""
function plabic_embedding(G::PlabicGraph; iterations=2500,restarts=24,
                          method=:harmonic)
    H = layout_plabic_graph(G;method=method,iterations=iterations,restarts=restarts)
    _edge_crossing_count(H.positions,H.edges)==0 ||
        error("could not find a crossing-free embedding; increase restarts")
    boundary_edges = [(i,mod1(i+1,H.n)) for i in 1:H.n]
    augmented = unique(vcat(H.edges,boundary_edges))
    rotation = _rotation_from_positions(H,augmented)
    all_faces,all_halfedges = _enumerate_faces(rotation,augmented)
    exterior = argmin([_signed_polygon_area(f,H.positions) for f in all_faces])
    kept = [i for i in eachindex(all_faces) if i!=exterior]
    remap = Dict(old=>new for (new,old) in enumerate(kept))
    faces = all_faces[kept]
    halfedges = Dict(edge=>remap[id] for (edge,id) in all_halfedges if haskey(remap,id))
    trips = [_trace_trip(H,rotation,i) for i in 1:H.n]
    endpoints = [last(path) for path in trips]
    endpoints == abs.(H.permutation) ||
        error("embedding trip endpoints $endpoints do not match $(abs.(H.permutation))")
    return PlabicEmbedding(H,rotation,faces,halfedges,trips)
end

graph_faces(E::PlabicEmbedding)=copy.(E.faces)
graph_trips(E::PlabicEmbedding)=copy.(E.trips)
graph_faces(G::PlabicGraph;kwargs...)=graph_faces(plabic_embedding(G;kwargs...))
graph_trips(G::PlabicGraph;kwargs...)=graph_trips(plabic_embedding(G;kwargs...))

function _face_label_cardinality(permutation,side)
    left_count=count(eachindex(permutation)) do i
        permutation[i] < 0 || (permutation[i] > 0 && permutation[i] < i)
    end
    return side==:left ? left_count : length(permutation)-left_count
end

"""
    face_labels(E; convention=(:source, :left))

Label every disk face using oriented trips.  The first convention entry chooses
the number attached to a trip (`:source` or `:target`); the second chooses the
side of the oriented trip (`:left` or `:right`).
"""
function face_labels(E::PlabicEmbedding; convention=(:source,:left))
    convention isa Tuple && length(convention)==2 ||
        throw(ArgumentError("convention must be (:source, :left), (:target, :left), (:source, :right), or (:target, :right)"))
    label_kind,side=convention
    label_kind in (:source,:target) ||
        throw(ArgumentError("the first convention entry must be :source or :target"))
    side in (:left,:right) ||
        throw(ArgumentError("the second convention entry must be :left or :right"))
    labels = [Set{Int}() for _ in E.faces]
    undirected_edges = Set((min(u,v),max(u,v)) for (u,v) in keys(E.halfedge_face))
    for i in 1:E.graph.n
        path = E.trips[i]
        # A decorated fixed point is represented by a lollipop trip i→v→i.
        # Its strand is a small oriented loop.  A white lollipop is clockwise,
        # so the whole disk lies to its left; a black lollipop is
        # counterclockwise, so no disk face lies to its left.  Right-side
        # conventions exchange these two cases.
        if first(path)==last(path)
            color=E.graph.colors[path[2]]
            contributes=(side==:left && color==:white) ||
                        (side==:right && color==:black)
            if contributes
                trip_label=label_kind==:source ? i : last(path)
                for label in labels
                    push!(label,trip_label)
                end
            end
            continue
        end
        barriers = Set((min(u,v),max(u,v)) for (u,v) in zip(path[1:end-1],path[2:end]))
        dual = [Set{Int}() for _ in E.faces]
        for (u,v) in undirected_edges
            (u,v) in barriers && continue
            left = get(E.halfedge_face,(u,v),0)
            right = get(E.halfedge_face,(v,u),0)
            if left!=0 && right!=0 && left!=right
                push!(dual[left],right); push!(dual[right],left)
            end
        end
        reached = Set{Int}()
        queue = Int[]
        for (u,v) in zip(path[1:end-1],path[2:end])
            # `_left_successor` traces the geometric left face of an oriented
            # half-edge, independently of where boundary label 1 is placed.
            directed_edge = side==:left ? (u,v) : (v,u)
            face = get(E.halfedge_face,directed_edge,0)
            if face!=0 && !(face in reached)
                push!(reached,face); push!(queue,face)
            end
        end
        while !isempty(queue)
            face = popfirst!(queue)
            for neighbor in dual[face]
                if !(neighbor in reached)
                    push!(reached,neighbor); push!(queue,neighbor)
                end
            end
        end
        trip_label=label_kind==:source ? i : last(path)
        for face in reached
            push!(labels[face],trip_label)
        end
    end
    expected_size=_face_label_cardinality(E.graph.permutation,side)
    all(length(label)==expected_size for label in labels) ||
        error("face-label construction with convention $convention did not produce $expected_size-subsets: $labels")
    return [sort(collect(label)) for label in labels]
end

function face_labels(G::PlabicGraph;convention=(:source,:left),kwargs...)
    return face_labels(plabic_embedding(G;kwargs...);convention=convention)
end

function _point_in_polygon(x,y,polygon)
    inside=false
    j=length(polygon)
    for i in eachindex(polygon)
        xi,yi=polygon[i]; xj,yj=polygon[j]
        if (yi>y)!=(yj>y) && x < (xj-xi)*(y-yi)/(yj-yi)+xi
            inside=!inside
        end
        j=i
    end
    return inside
end

function _segment_distance(x,y,a,b)
    ax,ay=a; bx,by=b
    dx,dy=bx-ax,by-ay
    t=(dx==0 && dy==0) ? 0.0 : clamp(((x-ax)*dx+(y-ay)*dy)/(dx*dx+dy*dy),0,1)
    return hypot(x-(ax+t*dx),y-(ay+t*dy))
end

function _face_visual_center(face,positions)
    polygon=positions[face]
    xmin,xmax=extrema(first.(polygon)); ymin,ymax=extrema(last.(polygon))
    best=((xmin+xmax)/2,(ymin+ymax)/2); best_distance=-Inf
    cx,cy=best; half=max(xmax-xmin,ymax-ymin)/2
    for refinement in 1:4
        grid=range(-half,half;length=25)
        for dx in grid,dy in grid
            x,y=cx+dx,cy+dy
            _point_in_polygon(x,y,polygon) || continue
            distance=minimum(_segment_distance(x,y,polygon[j],polygon[mod1(j+1,length(polygon))])
                             for j in eachindex(polygon))
            if distance>best_distance
                best=(x,y);best_distance=distance
            end
        end
        cx,cy=best;half/=5
    end
    return best
end

function _face_barycenter(face,positions)
    polygon=positions[face]
    twice_area=0.0;cx=0.0;cy=0.0
    for j in eachindex(polygon)
        x1,y1=polygon[j];x2,y2=polygon[mod1(j+1,length(polygon))]
        cross=x1*y2-x2*y1
        twice_area+=cross
        cx+=(x1+x2)*cross;cy+=(y1+y2)*cross
    end
    if abs(twice_area)<1e-12
        return (sum(first.(polygon))/length(polygon),
                sum(last.(polygon))/length(polygon))
    end
    return (cx/(3twice_area),cy/(3twice_area))
end

function _face_label_layout(face,positions,text;max_size=10.0,
                            pixels_per_unit=250.0,font_unit_scale=1.0)
    polygon=positions[face]
    center=_face_barycenter(face,positions)
    _point_in_polygon(center...,polygon) ||
        (center=_face_visual_center(face,positions))

    # Approximate the rendered label by a centered rectangle.  Sampling its
    # interior as well as its boundary handles concave faces more safely than a
    # single inradius estimate.  The character-width factor is deliberately a
    # little conservative for braces, commas, and multi-digit labels.
    function fits(font_size)
        half_width=(0.65length(text)*font_size*font_unit_scale+3)/pixels_per_unit
        half_height=(0.90font_size*font_unit_scale+2)/pixels_per_unit
        for sx in range(-1,1;length=7),sy in range(-1,1;length=5)
            _point_in_polygon(center[1]+sx*half_width,
                              center[2]+sy*half_height,polygon) || return false
        end
        return true
    end
    low=0.5;high=Float64(max_size)
    for _ in 1:18
        middle=(low+high)/2
        fits(middle) ? (low=middle) : (high=middle)
    end
    return center,min(Float64(max_size),low)
end

function _strand_lane_data(trips)
    users=Dict{Tuple{Int,Int},Vector{Int}}()
    for (i,path) in enumerate(trips), (u,v) in zip(path[1:end-1],path[2:end])
        edge=(min(u,v),max(u,v))
        push!(get!(users,edge,Int[]),i)
    end
    for list in values(users)
        sort!(unique!(list))
    end
    return users
end

function _strand_points(path,trip_index,positions,users,spacing)
    segment_offsets=NTuple{2,Float64}[]
    for (u,v) in zip(path[1:end-1],path[2:end])
        x1,y1=positions[u];x2,y2=positions[v]
        dx,dy=x2-x1,y2-y1;segment_length=hypot(dx,dy)
        lanes=users[(min(u,v),max(u,v))]
        lane=findfirst(==(trip_index),lanes)-(length(lanes)+1)/2
        direction=u<v ? 1.0 : -1.0
        push!(segment_offsets,(direction*lane*spacing*(-dy/segment_length),
                               direction*lane*spacing*(dx/segment_length)))
    end
    points=NTuple{2,Float64}[]
    for j in eachindex(path)
        offset = if j==1
            segment_offsets[1]
        elseif j==length(path)
            segment_offsets[end]
        else
            ((segment_offsets[j-1][1]+segment_offsets[j][1])/2,
             (segment_offsets[j-1][2]+segment_offsets[j][2])/2)
        end
        push!(points,(positions[path[j]][1]+offset[1],
                      positions[path[j]][2]+offset[2]))
    end
    return points
end

function _medial_strand_points(path,positions,colors,n;
                               samples=14,wiggle=0.10,boundary_offset=0.075)
    midpoint(u,v)=((positions[u][1]+positions[v][1])/2,
                   (positions[u][2]+positions[v][2])/2)
    # A lollipop trip is an oriented closed loop based at one boundary vertex.
    # It crosses the lollipop edge at its midpoint, circles the internal vertex
    # clockwise for white and counterclockwise for black, and returns through
    # the same midpoint.
    if length(path)==3 && first(path)==last(path)
        boundary,vertex=path[1],path[2]
        bx,by=positions[boundary];vx,vy=positions[vertex]
        color=colors[vertex]
        side=color==:white ? 1.0 : -1.0
        angle=atan(by,bx)
        midpoint_point=midpoint(boundary,vertex)
        radius=hypot(midpoint_point[1]-vx,midpoint_point[2]-vy)
        start_angle=atan(midpoint_point[2]-vy,midpoint_point[1]-vx)
        direction=color==:white ? -1.0 : 1.0
        points=NTuple{2,Float64}[(cos(angle+side*boundary_offset),
                                  sin(angle+side*boundary_offset)),
                                 midpoint_point]
        for step in 1:(2samples)
            theta=start_angle+direction*2pi*step/(2samples)
            push!(points,(vx+radius*cos(theta),vy+radius*sin(theta)))
        end
        push!(points,(cos(angle-side*boundary_offset),
                      sin(angle-side*boundary_offset)))
        return points
    end

    # At a boundary edge the two medial strands meet the boundary on opposite
    # sides.  The outgoing strand starts on the right when the adjacent
    # internal vertex is white, and on the left when it is black.  The incoming
    # endpoint therefore uses the opposite side.
    start_color=colors[path[2]]
    end_color=colors[path[end-1]]
    start_side=start_color==:white ? 1.0 : -1.0
    end_side=end_color==:white ? -1.0 : 1.0
    start_angle=atan(positions[path[1]][2],positions[path[1]][1])+
                start_side*boundary_offset
    end_angle=atan(positions[path[end]][2],positions[path[end]][1])+
              end_side*boundary_offset
    points=NTuple{2,Float64}[(cos(start_angle),sin(start_angle)),midpoint(path[1],path[2])]
    for j in 2:(length(path)-1)
        v=path[j]
        v<=n && continue
        incoming=midpoint(path[j-1],v);outgoing=midpoint(v,path[j+1])
        vx,vy=positions[v]
        a1=atan(incoming[2]-vy,incoming[1]-vx)
        a2=atan(outgoing[2]-vy,outgoing[1]-vx)
        delta=colors[v]==:white ? -mod(a1-a2,2pi) : mod(a2-a1,2pi)
        r1=hypot(incoming[1]-vx,incoming[2]-vy)
        r2=hypot(outgoing[1]-vx,outgoing[2]-vy)
        middle_angle=a1+delta/2
        offset=(wiggle*min(r1,r2)*cos(middle_angle),
                wiggle*min(r1,r2)*sin(middle_angle))
        c1=(vx+0.16*(incoming[1]-vx)+offset[1],
            vy+0.16*(incoming[2]-vy)+offset[2])
        c2=(vx+0.16*(outgoing[1]-vx)+offset[1],
            vy+0.16*(outgoing[2]-vy)+offset[2])
        for step in 1:samples
            t=step/samples
            s=1-t
            x=s^3*incoming[1]+3s^2*t*c1[1]+3s*t^2*c2[1]+t^3*outgoing[1]
            y=s^3*incoming[2]+3s^2*t*c1[2]+3s*t^2*c2[2]+t^3*outgoing[2]
            push!(points,(x,y))
        end
    end
    push!(points,(cos(end_angle),sin(end_angle)))
    return points
end

"""
    plabic_graph(p)

Construct a canonical bridge plabic graph for the decorated permutation `p`.
The construction follows the bridge recursion in the reference Mathematica
package. Negative fixed points become white lollipops and positive fixed points
become black lollipops.
"""
function plabic_graph(input::AbstractVector{<:Integer})
    p = Int.(input)
    _validate_decorated_permutation(p)
    n = length(p)
    base, removal_bridges = _bridge_decomposition(p)

    colors = fill(:boundary, n)
    positions = NTuple{2,Float64}[]
    for i in 1:n
        theta = -2pi * (i - 1) / max(n, 1) - pi / 2
        push!(positions, (cos(theta), sin(theta)))
    end

    edges = Tuple{Int,Int}[]
    radial_depth = zeros(Int, n)
    for i in 1:n
        push!(colors, base[i] == i ? :black : :white)
        theta = -2pi * (i - 1) / max(n, 1) - pi / 2
        push!(positions, (0.82cos(theta), 0.82sin(theta)))
        push!(edges, (i, length(colors)))
        radial_depth[i] = 1
    end

    construction_bridges = reverse(removal_bridges)
    for (a, b) in construction_bridges
        radial_depth[a] += 1
        radial_depth[b] += 1
        radius_a = max(0.28, 0.88 - 0.105radial_depth[a])
        radius_b = max(0.28, 0.88 - 0.105radial_depth[b])
        theta_a = -2pi * (a - 1) / n - pi / 2
        theta_b = -2pi * (b - 1) / n - pi / 2

        white = length(colors) + 1
        push!(colors, :white)
        push!(positions, (radius_a*cos(theta_a), radius_a*sin(theta_a)))
        black = length(colors) + 1
        push!(colors, :black)
        push!(positions, (radius_b*cos(theta_b), radius_b*sin(theta_b)))

        _replace_boundary_endpoint!(edges, a, white)
        _replace_boundary_endpoint!(edges, b, black)
        append!(edges, [(a, white), (white, black), (black, b)])
    end

    colors, edges, positions = _reduce_bipartite_graph(colors, edges, positions, n)
    G = PlabicGraph(n, colors, edges, positions, p, construction_bridges)
    validate_plabic_graph(G)
    return G
end

"""
    plabic_graph_from_drawing(n, internal_colors, edges, internal_positions; reduce=true)

Create a plabic graph from a disk drawing. Boundary vertices `1:n` are placed
clockwise on the unit circle; internal vertices are numbered consecutively and
must have color `:black` or `:white`. The straight-line drawing must be planar
and every boundary vertex must be univalent. When `reduce=true`, same-color
internal edges are contracted and bivalent internal vertices are suppressed
before the decorated trip permutation is computed.
"""
function plabic_graph_from_drawing(n::Integer,
                                   internal_colors::AbstractVector,
                                   input_edges::AbstractVector,
                                   internal_positions::AbstractVector;
                                   reduce=true)
    n = Int(n)
    n >= 1 || throw(ArgumentError("the boundary count must be positive"))
    length(internal_colors) == length(internal_positions) ||
        throw(ArgumentError("every internal vertex needs one position"))
    colors = vcat(fill(:boundary,n), Symbol.(internal_colors))
    all(c -> c in (:black,:white), colors[(n+1):end]) ||
        throw(ArgumentError("internal colors must be :black or :white"))
    positions = NTuple{2,Float64}[]
    for i in 1:n
        theta = -2pi*(i-1)/n - pi/2
        push!(positions,(cos(theta),sin(theta)))
    end
    for p in internal_positions
        length(p) == 2 || throw(ArgumentError("positions must be planar pairs"))
        x,y = Float64(p[1]),Float64(p[2])
        isfinite(x) && isfinite(y) && x*x+y*y < 0.99^2 ||
            throw(ArgumentError("internal vertices must lie strictly inside the disk"))
        push!(positions,(x,y))
    end
    edges = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    for raw_edge in input_edges
        length(raw_edge) == 2 || throw(ArgumentError("edges must be vertex pairs"))
        u,v = Int(raw_edge[1]),Int(raw_edge[2])
        1 <= u <= length(colors) && 1 <= v <= length(colors) ||
            throw(ArgumentError("edge ($u,$v) uses a missing vertex"))
        u != v || throw(ArgumentError("loops are not allowed"))
        edge=(min(u,v),max(u,v))
        edge in seen || (push!(seen,edge);push!(edges,edge))
    end
    degrees=zeros(Int,length(colors))
    for (u,v) in edges
        degrees[u]+=1;degrees[v]+=1
    end
    all(degrees[1:n].==1) ||
        throw(ArgumentError("every boundary vertex must have exactly one edge"))
    reduce && ((colors,edges,positions)=_reduce_bipartite_graph(colors,edges,positions,n))
    _edge_crossing_count(positions,edges)==0 ||
        throw(ArgumentError("the drawing has crossing edges; move vertices or remove an edge"))
    provisional=PlabicGraph(n,colors,edges,positions,collect(1:n),Tuple{Int,Int}[])
    rotation=_rotation_from_positions(provisional,
                                      unique(vcat(edges,[(i,mod1(i+1,n)) for i in 1:n])))
    trips=[_trace_trip(provisional,rotation,i) for i in 1:n]
    endpoints=[last(path) for path in trips]
    sort(endpoints)==collect(1:n) ||
        throw(ArgumentError("the trips do not define a permutation; check connectivity and the embedding"))
    permutation=copy(endpoints)
    for i in 1:n
        endpoints[i]==i || continue
        neighbor=only([v for (u,v) in edges if u==i] ∪ [u for (u,v) in edges if v==i])
        permutation[i]=colors[neighbor]==:white ? -i : i
    end
    _validate_decorated_permutation(permutation)
    result=PlabicGraph(n,colors,edges,positions,permutation,Tuple{Int,Int}[])
    validate_plabic_graph(result;require_reduced=false)
    if reduce && !is_reduced(result)
        reduced, _ = _reduce_by_explicit_moves(result; animate=false)
        return reduced
    end
    return result
end

"""
    trip_permutation(G)

Return the decorated trip permutation represented by a bridge-constructed
plabic graph. The value is part of the graph's construction certificate.
"""
trip_permutation(G::PlabicGraph) = copy(G.permutation)

function _plabic_degrees(G::PlabicGraph)
    degrees = zeros(Int, length(G.colors))
    for (u, v) in G.edges
        degrees[u] += 1
        degrees[v] += 1
    end
    return degrees
end

"""Return the valence of every vertex of a plabic graph."""
vertex_degrees(G::PlabicGraph) = _plabic_degrees(G)

function _canonical_cycle4(cycle)
    rotations = [circshift(cycle, -(j-1)) for j in 1:4]
    reversed = reverse(cycle)
    append!(rotations, [circshift(reversed, -(j-1)) for j in 1:4])
    return minimum(rotations)
end

"""
    square_cycles(G)

List induced alternating internal four-cycles that are candidates for square
moves. Because `PlabicGraph` does not yet store a rotation system, a separating
four-cycle can also appear; the caller should select a cycle that is a face in
the displayed embedding.
"""
function square_cycles(G::PlabicGraph)
    adjacency = [Set{Int}() for _ in G.colors]
    edge_set = Set{Tuple{Int,Int}}()
    for (u, v) in G.edges
        push!(adjacency[u], v); push!(adjacency[v], u)
        push!(edge_set, u < v ? (u,v) : (v,u))
    end
    found = Dict{NTuple{4,Int},Vector{Int}}()
    internal = (G.n+1):length(G.colors)
    for u in internal, v in (u+1):length(G.colors)
        ((u,v) in edge_set) && continue
        common = sort(collect(intersect(adjacency[u], adjacency[v])))
        common = filter(>(G.n), common)
        for aidx in 1:(length(common)-1), bidx in (aidx+1):length(common)
            a, b = common[aidx], common[bidx]
            edge = a < b ? (a,b) : (b,a)
            edge in edge_set && continue
            cycle = _canonical_cycle4([u,a,v,b])
            colors = G.colors[cycle]
            all(colors[j] != colors[mod1(j+1,4)] for j in 1:4) || continue
            found[Tuple(cycle)] = cycle
        end
    end
    return sort(collect(values(found)))
end

function _validate_square_cycle(G, cycle)
    length(cycle) == 4 || throw(ArgumentError("a square cycle must contain four vertices"))
    length(unique(cycle)) == 4 || throw(ArgumentError("square vertices must be distinct"))
    all(v -> G.n < v <= length(G.colors), cycle) ||
        throw(ArgumentError("a square move requires four internal vertices"))
    edge_set = Set(u < v ? (u,v) : (v,u) for (u,v) in G.edges)
    for j in 1:4
        u, v = cycle[j], cycle[mod1(j+1,4)]
        (min(u,v),max(u,v)) in edge_set || throw(ArgumentError("the selected vertices do not form a cycle"))
        G.colors[u] != G.colors[v] || throw(ArgumentError("square colors must alternate"))
    end
    for (u,v) in ((cycle[1],cycle[3]), (cycle[2],cycle[4]))
        (min(u,v),max(u,v)) in edge_set && throw(ArgumentError("the selected four-cycle has a diagonal"))
    end
    return true
end

"""
    square_move(G, cycle)
    square_move(G)

Perform a square move on the selected facial four-cycle. Higher-valence square
vertices are split so that the square vertex is trivalent, the four colors are
switched, and the graph is reduced. The no-cycle form is available when exactly
one candidate square exists.
"""
function square_move(G::PlabicGraph, input_cycle::AbstractVector{<:Integer})
    cycle = Int.(input_cycle)
    _validate_square_cycle(G, cycle)
    colors, edges, positions = copy(G.colors), copy(G.edges), copy(G.positions)
    cycle_set = Set(cycle)

    # When a square vertex has several external branches, collect those branches
    # at a new hub so that the vertex on the square itself has valence three.
    for v in cycle
        external = [u == v ? w : u for (u,w) in edges
                    if (u == v || w == v) && !((u == v ? w : u) in cycle_set)]
        length(external) <= 1 && continue
        old_color = colors[v]
        edges = [(u,w) for (u,w) in edges if !((u == v && w in external) ||
                                                (w == v && u in external))]
        hub = length(colors) + 1
        push!(colors, old_color)
        px = sum(positions[u][1] for u in external)/length(external)
        py = sum(positions[u][2] for u in external)/length(external)
        push!(positions, ((positions[v][1]+px)/2, (positions[v][2]+py)/2))
        push!(edges, (v,hub))
        append!(edges, [(hub,u) for u in external])
    end

    for v in cycle
        colors[v] = colors[v] == :black ? :white : :black
    end
    colors, edges, positions = _reduce_bipartite_graph(colors, edges, positions, G.n)
    result = PlabicGraph(G.n, colors, edges, positions, copy(G.permutation), copy(G.bridges))
    validate_plabic_graph(result)
    return result
end

function _square_move_animation_stages(G::PlabicGraph,
                                       input_cycle::AbstractVector{<:Integer};
                                       require_reduced=true)
    cycle=Int.(input_cycle)
    _validate_square_cycle(G,cycle)
    colors,edges,positions=copy(G.colors),copy(G.edges),copy(G.positions)
    cycle_set=Set(cycle)
    stages=Tuple{String,PlabicGraph}[]
    snapshot(caption)=push!(stages,(caption,PlabicGraph(G.n,copy(colors),copy(edges),
                                                       copy(positions),copy(G.permutation),
                                                       copy(G.bridges))))
    split_count=0
    for v in cycle
        external=[u==v ? w : u for (u,w) in edges
                  if (u==v || w==v) && !((u==v ? w : u) in cycle_set)]
        length(external)<=1 && continue
        old_color=colors[v]
        edges=[(u,w) for (u,w) in edges if !((u==v && w in external) ||
                                              (w==v && u in external))]
        hub=length(colors)+1
        push!(colors,old_color)
        px=sum(positions[u][1] for u in external)/length(external)
        py=sum(positions[u][2] for u in external)/length(external)
        push!(positions,((positions[v][1]+px)/2,(positions[v][2]+py)/2))
        push!(edges,(v,hub));append!(edges,[(hub,u) for u in external])
        split_count+=1
        snapshot("Make square vertex $split_count trivalent")
    end
    split_count==0 && snapshot("The square vertices are already trivalent")
    for v in cycle colors[v]=colors[v]==:black ? :white : :black end
    snapshot("Switch the four square colors")
    reduction_step=0
    while true
        new_colors,new_edges,new_positions,changed,contracted=
            _contract_one_same_color_edge(colors,edges,positions,G.n)
        if changed
            colors,edges,positions=new_colors,new_edges,new_positions
            reduction_step+=1
            u,v=contracted
            snapshot("Contract same-color edge $(u)–$(v) · step $reduction_step")
            continue
        end
        new_colors,new_edges,new_positions,changed=
            _suppress_one_bivalent(colors,edges,positions,G.n)
        if changed
            colors,edges,positions=new_colors,new_edges,new_positions
            reduction_step+=1
            snapshot("Suppress one bivalent vertex · step $reduction_step")
            continue
        end
        break
    end
    result=PlabicGraph(G.n,colors,edges,positions,copy(G.permutation),copy(G.bridges))
    validate_plabic_graph(result;require_reduced=require_reduced)
    snapshot(require_reduced ? "Reduced square-move graph" :
             "Complete one facial square move")
    return result,stages
end

function _facial_square_cycles(G::PlabicGraph)
    _edge_crossing_count(G.positions,G.edges)==0 || return Vector{Vector{Int}}()
    boundary_edges=[(i,mod1(i+1,G.n)) for i in 1:G.n]
    augmented=unique(vcat(G.edges,boundary_edges))
    rotation=try
        _rotation_from_positions(G,augmented)
    catch
        return Vector{Vector{Int}}()
    end
    faces=try
        first(_enumerate_faces(rotation,augmented))
    catch
        return Vector{Vector{Int}}()
    end
    cycles=Vector{Vector{Int}}()
    seen=Set{NTuple{4,Int}}()
    for face in faces
        cycle=Int[]
        for v in face
            (isempty(cycle) || cycle[end]!=v) && push!(cycle,v)
        end
        length(cycle)>1 && cycle[1]==cycle[end] && pop!(cycle)
        length(cycle)==4 || continue
        length(unique(cycle))==4 || continue
        all(>(G.n),cycle) || continue
        all(G.colors[cycle[j]]!=G.colors[cycle[mod1(j+1,4)]] for j in 1:4) || continue
        key=Tuple(_canonical_cycle4(cycle))
        key in seen && continue
        push!(seen,key);push!(cycles,collect(key))
    end
    return sort(cycles)
end

function _drop_one_boundaryless_component(G::PlabicGraph)
    adjacency=[Int[] for _ in G.colors]
    for (u,v) in G.edges
        push!(adjacency[u],v);push!(adjacency[v],u)
    end
    visited=falses(length(G.colors))
    for seed in (G.n+1):length(G.colors)
        visited[seed] && continue
        component=Int[];stack=[seed];visited[seed]=true
        while !isempty(stack)
            v=pop!(stack);push!(component,v)
            for w in adjacency[v]
                visited[w] && continue
                visited[w]=true;push!(stack,w)
            end
        end
        any(<=(G.n),component) && continue
        removed=Set(component)
        mapping=zeros(Int,length(G.colors))
        colors=Symbol[];positions=NTuple{2,Float64}[]
        for old in eachindex(G.colors)
            old in removed && continue
            mapping[old]=length(colors)+1
            push!(colors,G.colors[old]);push!(positions,G.positions[old])
        end
        edges=Tuple{Int,Int}[]
        for (u,v) in G.edges
            (u in removed || v in removed) && continue
            push!(edges,(mapping[u],mapping[v]))
        end
        return PlabicGraph(G.n,colors,edges,positions,copy(G.permutation),
                           Tuple{Int,Int}[]),true,length(component)
    end
    return G,false,0
end

_explicit_reduction_key(G::PlabicGraph)=
    (Tuple(G.colors),
     Tuple(sort([(min(u,v),max(u,v)) for (u,v) in G.edges])))

function _strand_reduction_metric(G::PlabicGraph)
    augmented=unique(vcat(G.edges,[(i,mod1(i+1,G.n)) for i in 1:G.n]))
    rotation=try
        _rotation_from_positions(G,augmented)
    catch
        return (typemax(Int)÷4,typemax(Int)÷4)
    end
    trips=try
        [_trace_trip(G,rotation,i) for i in 1:G.n]
    catch
        return (typemax(Int)÷4,typemax(Int)÷4)
    end
    bad_crossings=0
    lens_length=0
    edge_orders=Vector{Dict{Tuple{Int,Int},Int}}(undef,G.n)
    for i in 1:G.n
        edge_orders[i]=Dict((min(u,v),max(u,v))=>j
            for (j,(u,v)) in enumerate(zip(trips[i][1:end-1],trips[i][2:end])))
    end
    for i in 1:(G.n-1),j in (i+1):G.n
        shared=collect(intersect(keys(edge_orders[i]),keys(edge_orders[j])))
        for a in 1:(length(shared)-1),b in (a+1):length(shared)
            first_delta=edge_orders[i][shared[a]]-edge_orders[i][shared[b]]
            second_delta=edge_orders[j][shared[a]]-edge_orders[j][shared[b]]
            sign(first_delta)==sign(second_delta) || continue
            bad_crossings+=1
            lens_length+=abs(first_delta)+abs(second_delta)
        end
    end
    return bad_crossings,lens_length
end

function _heap_push!(heap,item)
    push!(heap,item)
    index=length(heap)
    while index>1
        parent=index÷2
        heap[parent][1]<=heap[index][1] && break
        heap[parent],heap[index]=heap[index],heap[parent]
        index=parent
    end
    return heap
end

function _heap_pop!(heap)
    first_item=heap[1]
    last_item=pop!(heap)
    isempty(heap) && return first_item
    heap[1]=last_item
    index=1
    while true
        left=2index
        left>length(heap) && break
        right=left+1
        child=right<=length(heap) && heap[right][1]<heap[left][1] ? right : left
        heap[index][1]<=heap[child][1] && break
        heap[index],heap[child]=heap[child],heap[index]
        index=child
    end
    return first_item
end

function _explicit_reduction_priority(G,depth,serial)
    graph_dimension=length(G.edges)-(length(G.colors)-G.n)
    cell_dimension=dimensionOfPermutation(decorated_excedances(G.permutation),
                                          G.n,G.permutation)
    bad,lens=_strand_reduction_metric(G)
    return (abs(graph_dimension-cell_dimension),bad,lens,depth,length(G.colors),serial)
end

function _reduce_by_explicit_moves(G::PlabicGraph;animate=false,max_depth=18,max_states=30000)
    validate_plabic_graph(G;require_reduced=false)
    initial_stages=Tuple{String,PlabicGraph}[]
    current=G
    while true
        next,changed,count=_drop_one_boundaryless_component(current)
        changed || break
        current=next
        animate && push!(initial_stages,
            (count==2 ? "Remove one isolated dipole" :
                        "Remove one boundaryless component ($count vertices)",current))
    end
    is_reduced(current) && return current,initial_stages

    # Each node stores only its parent and the square used to reach it.  This
    # avoids retaining thousands of complete animation paths during search.
    nodes=Tuple{PlabicGraph,Int,Vector{Int},Int}[(current,0,Int[],0)]
    # Search states can gain vertices during trivalentization and lose them
    # during contractions.  Do not let Julia infer a fixed NTuple length from
    # the initial state's color tuple.
    visited=Set{Any}([_explicit_reduction_key(current)])
    serial=1
    heap=Tuple{NTuple{6,Int},Int}[]
    _heap_push!(heap,(_explicit_reduction_priority(current,0,serial),1))
    goal=0
    while !isempty(heap) && length(visited)<=max_states
        _,node_index=_heap_pop!(heap)
        graph,_,_,depth=nodes[node_index]
        depth>=max_depth && continue
        for cycle in _facial_square_cycles(graph)
            moved,_=try
                _square_move_animation_stages(graph,cycle;require_reduced=false)
            catch
                continue
            end
            moved=try
                _with_recomputed_trip_permutation(moved)
            catch
                continue
            end
            while true
                next,changed,count=_drop_one_boundaryless_component(moved)
                changed || break
                moved=next
            end
            key=_explicit_reduction_key(moved)
            key in visited && continue
            push!(visited,key)
            push!(nodes,(moved,node_index,copy(cycle),depth+1))
            child_index=length(nodes)
            if is_reduced(moved)
                goal=child_index
                break
            end
            serial+=1
            priority=_explicit_reduction_priority(moved,depth+1,serial)
            _heap_push!(heap,(priority,child_index))
            length(visited)>=max_states && break
        end
        goal>0 && break
    end
    if goal>0
        reduced=nodes[goal][1]
        animate || return reduced,Tuple{String,PlabicGraph}[]
        chain=Int[]
        node_index=goal
        while nodes[node_index][2]!=0
            push!(chain,node_index)
            node_index=nodes[node_index][2]
        end
        reverse!(chain)
        stages=copy(initial_stages)
        replay=current
        for index in chain
            cycle=nodes[index][3]
            replay,move_stages=
                _square_move_animation_stages(replay,cycle;require_reduced=false)
            append!(stages,move_stages)
            replay=_with_recomputed_trip_permutation(replay)
            while true
                next,changed,count=_drop_one_boundaryless_component(replay)
                changed || break
                replay=next
                push!(stages,(count==2 ? "Remove one isolated dipole" :
                                          "Remove one boundaryless component ($count vertices)",replay))
            end
        end
        return replay,stages
    end
    graph_dimension=length(current.edges)-(length(current.colors)-current.n)
    cell_dimension=dimensionOfPermutation(decorated_excedances(current.permutation),
                                          current.n,current.permutation)
    throw(ArgumentError("explicit local reduction search did not reach a reduced graph " *
                        "(graph dimension $graph_dimension, cell dimension $cell_dimension); " *
                        "the graph was left unchanged rather than replaced by a canonical graph"))
end

function square_move(G::PlabicGraph;face_label=nothing,iterations=2500,restarts=24)
    if !isnothing(face_label)
        return square_move_by_label(G,face_label;iterations=iterations,restarts=restarts)
    end
    candidates = square_cycles(G)
    length(candidates) == 1 ||
        throw(ArgumentError("expected exactly one square candidate, found $(length(candidates)); pass one from square_cycles(G)"))
    return square_move(G, only(candidates))
end

"""
    square_move_by_label(G, label; iterations=2500, restarts=24)
    square_move(G; face_label=label, ...)

Find the unique disk face carrying `label`, verify that its boundary is an
internal alternating square, and perform the square move on that face.
"""
function _square_cycle_by_label(G::PlabicGraph,
                                input_label::Union{AbstractVector{<:Integer},AbstractSet{<:Integer}};
                                iterations=2500,restarts=24)
    label=sort(Int.(collect(input_label)))
    expected_size=_face_label_cardinality(G.permutation,:left)
    length(label)==expected_size ||
        throw(ArgumentError("source-left face label must contain $expected_size elements"))
    E=plabic_embedding(G;iterations=iterations,restarts=restarts)
    labels=face_labels(E)
    matches=findall(==(label),labels)
    isempty(matches) && throw(ArgumentError("no face has label $label"))
    length(matches)==1 || throw(ArgumentError("face label $label is not unique"))
    boundary=E.faces[only(matches)]
    cycle=Int[]
    for v in boundary
        (isempty(cycle) || cycle[end]!=v) && push!(cycle,v)
    end
    length(cycle)>1 && cycle[1]==cycle[end] && pop!(cycle)
    length(cycle)==4 ||
        throw(ArgumentError("face $label has boundary $cycle, not four vertices"))
    all(>(G.n),cycle) ||
        throw(ArgumentError("face $label is a boundary face, not an internal square"))
    _validate_square_cycle(E.graph,cycle)
    return E.graph,cycle
end

function square_move_by_label(G::PlabicGraph,
                              input_label::Union{AbstractVector{<:Integer},AbstractSet{<:Integer}};
                              iterations=2500,restarts=24)
    embedded,cycle=_square_cycle_by_label(G,input_label;
                                          iterations=iterations,restarts=restarts)
    return square_move(embedded,cycle)
end


"""
    draw_plabic_graph(G; kwargs...)
    draw_plabic_graph(p; kwargs...)

Draw a plabic graph in the disk. Supported keywords include `size`, `radius`,
`vertex_size`, `line_width`, `labels`, `black_color`, and `white_color`.
Returns a Plots.jl plot object.
"""
function draw_plabic_graph(G::PlabicGraph; size=600, radius=1.0,
                           vertex_size=8, line_width=2.5, labels=true,
                           black_color=:black, white_color=:white,
                           optimize_layout=true, layout_method=:harmonic,
                           layout_iterations=1400,
                           layout_restarts=6, face_labels=false,
                           face_label_convention=(:source,:left),
                           face_label_color=:red, face_label_size=10,
                           adaptive_face_label_size=true,
                           face_label_position=:barycenter,
                           trips=false, strand_spacing=0.022,
                           strand_width=3.0, strand_alpha=0.9,
                           strand_arrows=true, strand_wiggle=0.10,
                           strand_samples=14, strand_boundary_offset=0.075)
    _load_plotting()
    face_data = nothing
    embedded_faces = nothing
    needs_embedding=face_labels || trips!==false
    trip_data=nothing
    if needs_embedding
        E = plabic_embedding(G;method=layout_method,
                             iterations=max(layout_iterations,2500),
                             restarts=max(layout_restarts,12))
        G = E.graph
        if face_labels
            face_data = getfield(@__MODULE__,:face_labels)(E;
                                      convention=face_label_convention)
            embedded_faces = E.faces
        end
        trip_data=E.trips
    elseif optimize_layout
        G = layout_plabic_graph(G; method=layout_method,
                                iterations=layout_iterations,
                                restarts=layout_restarts)
    end
    theta = range(0, 2pi; length=300)
    fig = Plots.plot(radius*cos.(theta), radius*sin.(theta);
                     color=:black, linewidth=1.5, label=false, aspect_ratio=:equal,
                     axis=false, ticks=false, legend=false, size=(size, size))
    for (u, v) in G.edges
        x1, y1 = G.positions[u]
        x2, y2 = G.positions[v]
        Plots.plot!(fig, radius .* [x1, x2], radius .* [y1, y2];
                    color=:black, linewidth=line_width, label=false)
    end

    if trips!==false
        selected = (trips===true || trips===:all) ? collect(1:G.n) : Int.(collect(trips))
        all(i->1<=i<=G.n,selected) || throw(ArgumentError("trip indices must lie in 1:$(G.n)"))
        colors=Plots.palette(:tab10,G.n)
        for i in selected
            points=_medial_strand_points(trip_data[i],G.positions,G.colors,G.n;
                                         samples=strand_samples,wiggle=strand_wiggle,
                                         boundary_offset=strand_boundary_offset)
            xs=radius .* first.(points);ys=radius .* last.(points)
            Plots.plot!(fig,xs,ys;color=colors[i],linewidth=strand_width,
                        alpha=strand_alpha,label=false,
                        arrow=strand_arrows ? :closed : false)
        end
    end

    for color in (:black, :white)
        ids = findall(==(color), G.colors)
        isempty(ids) && continue
        xs = radius .* first.(G.positions[ids])
        ys = radius .* last.(G.positions[ids])
        fillcolor = color == :black ? black_color : white_color
        Plots.scatter!(fig, xs, ys; markersize=vertex_size,
                       markercolor=fillcolor, markerstrokecolor=:black,
                       markerstrokewidth=1.5, label=false)
    end

    boundary_ids = 1:G.n
    Plots.scatter!(fig, radius .* first.(G.positions[boundary_ids]),
                   radius .* last.(G.positions[boundary_ids]); markersize=4,
                   markercolor=:black, label=false)
    if labels
        for i in boundary_ids
            x, y = G.positions[i]
            Plots.annotate!(fig, 1.12radius*x, 1.12radius*y,
                            Plots.text(string(i), 11, :black))
        end
    end
    if face_labels
        for (face,label) in zip(embedded_faces,face_data)
            text = "{" * join(label,",") * "}"
            label_size=Float64(face_label_size)
            x,y = if face_label_position==:barycenter
                if adaptive_face_label_size
                    center,label_size=_face_label_layout(face,G.positions,text;
                        max_size=face_label_size,pixels_per_unit=size/2.44,
                        font_unit_scale=4/3)
                    center
                else
                    _face_barycenter(face,G.positions)
                end
            elseif face_label_position==:visual_center
                _face_visual_center(face,G.positions)
            else
                throw(ArgumentError("face_label_position must be :barycenter or :visual_center"))
            end
            Plots.annotate!(fig,radius*x,radius*y,
                            Plots.text(text,label_size,face_label_color))
        end
    end
    Plots.plot!(fig; xlims=(-1.22radius, 1.22radius),
                ylims=(-1.22radius, 1.22radius))
    return fig
end

draw_plabic_graph(p::AbstractVector{<:Integer}; kwargs...) =
    draw_plabic_graph(plabic_graph(p); kwargs...)

"""
    compare_plabic_graphs(graphs...; titles=nothing, panel_size=600, kwargs...)

Draw two or more plabic graphs in independent side-by-side panels. Keyword
arguments are forwarded to `draw_plabic_graph`.
"""
function compare_plabic_graphs(graphs::PlabicGraph...;
                               titles=nothing,panel_size=600,kwargs...)
    length(graphs)>=2 || throw(ArgumentError("provide at least two graphs"))
    if !isnothing(titles)
        length(titles)==length(graphs) ||
            throw(ArgumentError("provide one title per graph"))
    end
    figures=Any[]
    for (i,G) in enumerate(graphs)
        fig=draw_plabic_graph(G;size=panel_size,kwargs...)
        isnothing(titles) || Plots.plot!(fig;title=string(titles[i]))
        push!(figures,fig)
    end
    return Plots.plot(figures...;layout=(1,length(figures)),
                      size=(panel_size*length(figures),panel_size))
end

function _face_cycle(face)
    cycle=Int[]
    for v in face
        (isempty(cycle) || cycle[end]!=v) && push!(cycle,v)
    end
    length(cycle)>1 && cycle[1]==cycle[end] && pop!(cycle)
    return cycle
end

function _interactive_graph_state(G;iterations=2500,restarts=12)
    E=plabic_embedding(G;iterations=iterations,restarts=restarts)
    labels=face_labels(E;convention=(:source,:left))
    counts=Dict{Tuple{Vararg{Int}},Int}()
    for label in labels
        key=Tuple(label); counts[key]=get(counts,key,0)+1
    end
    movable=Bool[]
    for (face,label) in zip(E.faces,labels)
        cycle=_face_cycle(face)
        valid=length(cycle)==4 && all(>(G.n),cycle) && counts[Tuple(label)]==1
        if valid
            valid=try
                _validate_square_cycle(E.graph,cycle)
            catch
                false
            end
        end
        push!(movable,valid)
    end
    return E,labels,movable
end

_json_number(x::Real)=isfinite(x) ? string(Float64(x)) : "null"
_json_ints(xs)=string("[",join(Int.(xs),","),"]")

function _boundary_dual_positions(E::PlabicEmbedding)
    boundary_arc_points=[NTuple{2,Float64}[] for _ in E.faces]
    for i in 1:E.graph.n
        j=mod1(i+1,E.graph.n)
        face_id=get(E.halfedge_face,(i,j),0)
        face_id==0 && (face_id=get(E.halfedge_face,(j,i),0))
        face_id==0 && continue
        x1,y1=E.graph.positions[i];x2,y2=E.graph.positions[j]
        mx,my=x1+x2,y1+y2
        norm=hypot(mx,my)
        midpoint=norm>1e-10 ? (mx/norm,my/norm) : E.graph.positions[i]
        push!(boundary_arc_points[face_id],midpoint)
    end
    result=Dict{Int,NTuple{2,Float64}}()
    for face_id in eachindex(boundary_arc_points)
        points=boundary_arc_points[face_id]
        isempty(points) && continue
        x=sum(first.(points))/length(points);y=sum(last.(points))/length(points)
        norm=hypot(x,y)
        result[face_id]=norm>1e-10 ? (x/norm,y/norm) : first(points)
    end
    return result
end

function _black_dual_faces(E::PlabicEmbedding)
    polygons=Vector{Vector{Int}}()
    G=E.graph
    for v in (G.n+1):length(G.colors)
        G.colors[v]==:black || continue
        face_ids=Int[]
        for u in E.rotation[v]
            face_id=get(E.halfedge_face,(v,u),0)
            face_id==0 && continue
            (isempty(face_ids) || face_ids[end]!=face_id) && push!(face_ids,face_id)
        end
        length(face_ids)>1 && face_ids[1]==face_ids[end] && pop!(face_ids)
        length(unique(face_ids))>=3 && push!(polygons,face_ids)
    end
    return polygons
end

_edge_key(u::Int,v::Int)=u<v ? (u,v) : (v,u)

function _perfect_orientation_matching(G::PlabicGraph,requested_sources=nothing)
    default_sources=Int.(minGrassmannNecklace(G.permutation)[1])
    chosen=requested_sources===nothing ? default_sources : Int.(requested_sources)
    length(chosen)==decorated_excedances(G.permutation) ||
        throw(ArgumentError("source set must contain exactly $(decorated_excedances(G.permutation)) boundary vertices"))
    length(unique(chosen))==length(chosen) || throw(ArgumentError("source vertices must be distinct"))
    all(i->1<=i<=G.n,chosen) || throw(ArgumentError("source vertices must lie in 1:$(G.n)"))
    sources=Set(chosen)
    matched=Set{Tuple{Int,Int}}()
    occupied=falses(length(G.colors))
    for i in 1:G.n
        neighbors=[u==i ? v : u for (u,v) in G.edges if u==i || v==i]
        length(neighbors)==1 || error("boundary vertex $i does not have degree one")
        v=only(neighbors)
        required=G.colors[v]==:white ? i in sources : !(i in sources)
        if required
            occupied[v] && error("boundary data do not admit an almost-perfect matching")
            push!(matched,_edge_key(i,v));occupied[v]=true
        end
    end
    internal_edges=[(u,v) for (u,v) in G.edges if u>G.n && v>G.n]
    function extend!()
        remaining=[v for v in (G.n+1):length(G.colors) if !occupied[v]]
        isempty(remaining) && return true
        v=remaining[argmin([count(e->(e[1]==u && !occupied[e[2]]) ||
                                         (e[2]==u && !occupied[e[1]]),internal_edges)
                             for u in remaining])]
        for (a,b) in internal_edges
            a==v ? (u=b) : b==v ? (u=a) : continue
            occupied[u] && continue
            occupied[v]=occupied[u]=true;push!(matched,_edge_key(u,v))
            extend!() && return true
            delete!(matched,_edge_key(u,v));occupied[v]=occupied[u]=false
        end
        return false
    end
    extend!() || error("could not find an almost-perfect matching for this graph")
    return sort!(collect(sources)),matched
end

function _directed_plabic_edges(G::PlabicGraph,matching)
    directed=Tuple{Int,Int,Int}[]
    for (edge_id,(a,b)) in enumerate(G.edges)
        ismatched=_edge_key(a,b) in matching
        if a<=G.n || b<=G.n
            boundary,internal=a<=G.n ? (a,b) : (b,a)
            c=G.colors[internal]
            from,to = c==:white ? (ismatched ? (boundary,internal) : (internal,boundary)) :
                                  (ismatched ? (internal,boundary) : (boundary,internal))
        else
            black,white=G.colors[a]==:black ? (a,b) : (b,a)
            from,to=ismatched ? (black,white) : (white,black)
        end
        push!(directed,(from,to,edge_id))
    end
    return directed
end

function _sneg(a::String)
    a=="0" && return "0";a=="1" && return "-1";a=="-1" && return "1"
    startswith(a,"-(") && endswith(a,")") && return a[3:end-1]
    return "-("*a*")"
end
function _sadd(a::String,b::String)
    a=="0" && return b;b=="0" && return a;a==_sneg(b) && return "0"
    return "("*a*" + "*b*")"
end
_ssub(a::String,b::String)=_sadd(a,_sneg(b))
function _smul(a::String,b::String)
    (a=="0" || b=="0") && return "0";a=="1" && return b;b=="1" && return a
    a=="-1" && return _sneg(b);b=="-1" && return _sneg(a)
    startswith(a,"-(") && return _sneg(_smul(_sneg(a),b))
    startswith(b,"-(") && return _sneg(_smul(a,_sneg(b)))
    return "("*a*")*("*b*")"
end
function _sdiv(a::String,b::String)
    a=="0" && return "0";b=="1" && return a;a==b && return "1"
    return "("*a*")/("*b*")"
end

function _symbolic_inverse(A::Matrix{String})
    n=size(A,1);size(A,2)==n || throw(DimensionMismatch("matrix must be square"))
    B=[A [i==j ? "1" : "0" for i in 1:n,j in 1:n]]
    for col in 1:n
        pivot=findfirst(r->B[r,col]!="0",col:n)
        pivot===nothing && error("singular symbolic path matrix")
        pivot=col+pivot-1
        if pivot!=col
            saved=copy(B[col,:]);B[col,:]=B[pivot,:];B[pivot,:]=saved
        end
        p=B[col,col];B[col,:]=[_sdiv(x,p) for x in B[col,:]]
        for row in 1:n
            row==col && continue
            factor=B[row,col];factor=="0" && continue
            B[row,:]=[_ssub(B[row,j],_smul(factor,B[col,j])) for j in axes(B,2)]
        end
    end
    return B[:,n+1:2n]
end

const _RPoly=Dict{Tuple,Int}
struct _RExpr
    num::_RPoly
    den::_RPoly
end
_rpconst(c::Int)=c==0 ? _RPoly() : _RPoly(()=>c)
_rpone()=_rpconst(1)
_rzero()=_RExpr(_rpconst(0),_rpone())
_rone()=_RExpr(_rpone(),_rpone())
_rvariable(name::String)=strip(name) in ("","1") ? _rone() :
    _RExpr(_RPoly((strip(name)=>1,)=>1),_rpone())

function _rpadd(a::_RPoly,b::_RPoly)
    c=copy(a)
    for (m,v) in b
        c[m]=get(c,m,0)+v;c[m]==0 && delete!(c,m)
    end
    c
end
_rpneg(a::_RPoly)=_RPoly(m=>-v for (m,v) in a)
function _rpmul(a::_RPoly,b::_RPoly)
    c=_RPoly()
    for (ma,ca) in a,(mb,cb) in b
        exponents=Dict{String,Int}()
        for p in ma exponents[p.first]=get(exponents,p.first,0)+p.second end
        for p in mb exponents[p.first]=get(exponents,p.first,0)+p.second end
        m=Tuple(sort(collect(exponents);by=first))
        c[m]=get(c,m,0)+ca*cb;c[m]==0 && delete!(c,m)
    end
    c
end

function _rp_content(p::_RPoly)
    isempty(p) && return 1,Dict{String,Int}()
    coefficient=foldl(gcd,abs.(collect(values(p)));init=0)
    variables=Set{String}(q.first for m in keys(p) for q in m)
    powers=Dict(v=>minimum(get(Dict(m),v,0) for m in keys(p)) for v in variables)
    return max(coefficient,1),powers
end
function _rpdivide_content(p::_RPoly,c::Int,powers::Dict{String,Int})
    result=_RPoly()
    for (m,value) in p
        exponents=Dict(m)
        for (v,e) in powers
            exponents[v]=get(exponents,v,0)-e;exponents[v]==0 && delete!(exponents,v)
        end
        result[Tuple(sort(collect(exponents);by=first))]=div(value,c)
    end
    result
end
function _rnormalize(num::_RPoly,den::_RPoly)
    isempty(den) && error("zero denominator in symbolic boundary measurement")
    isempty(num) && return _rzero()
    num==den && return _rone()
    cn,pn=_rp_content(num);cd,pd=_rp_content(den)
    common_c=gcd(cn,cd);common_p=Dict{String,Int}(
        v=>min(get(pn,v,0),get(pd,v,0)) for v in union(keys(pn),keys(pd)))
    num=_rpdivide_content(num,common_c,common_p);den=_rpdivide_content(den,common_c,common_p)
    first_den=first(values(den))
    if first_den<0
        num=_rpneg(num);den=_rpneg(den)
    end
    return _RExpr(num,den)
end
_rneg(a::_RExpr)=_RExpr(_rpneg(a.num),a.den)
_radd(a::_RExpr,b::_RExpr)=_rnormalize(_rpadd(_rpmul(a.num,b.den),_rpmul(b.num,a.den)),_rpmul(a.den,b.den))
_rsub(a::_RExpr,b::_RExpr)=_radd(a,_rneg(b))
_rmul(a::_RExpr,b::_RExpr)=_rnormalize(_rpmul(a.num,b.num),_rpmul(a.den,b.den))
_rdiv(a::_RExpr,b::_RExpr)=_rnormalize(_rpmul(a.num,b.den),_rpmul(a.den,b.num))
_riszero(a::_RExpr)=isempty(a.num)

function _rp_string(p::_RPoly;latex=false)
    isempty(p) && return "0"
    terms=sort(collect(p);by=x->string(first(x)))
    pieces=String[]
    for (m,c) in terms
        factors=String[]
        for (v,e) in m
            shown=latex ? replace(v,r"^([A-Za-z]+)_([A-Za-z0-9]+)$"=>s"\1_{\2}") : v
            push!(factors,e==1 ? shown : latex ? shown*"^{"*string(e)*"}" : shown*"^"*string(e))
        end
        monomial=join(factors,latex ? raw"\," : "*")
        magnitude=abs(c)
        body=isempty(monomial) ? string(magnitude) : magnitude==1 ? monomial :
             string(magnitude)*(latex ? raw"\," : "*")*monomial
        if isempty(pieces)
            push!(pieces,c<0 ? "-"*body : body)
        else
            push!(pieces,(c<0 ? " - " : " + ")*body)
        end
    end
    join(pieces)
end
function _rstring(a::_RExpr;latex=false)
    numerator=_rp_string(a.num;latex=latex);denominator=_rp_string(a.den;latex=latex)
    denominator=="1" && return numerator
    latex && return raw"\frac{"*numerator*"}{"*denominator*"}"
    num_complex=length(a.num)>1
    den_monomial=length(a.den)==1 ? only(keys(a.den)) : ()
    den_complex=length(a.den)>1 || length(den_monomial)>1 ||
        (!isempty(den_monomial) && first(den_monomial).second>1)
    return (num_complex ? "("*numerator*")" : numerator)*"/"*
           (den_complex ? "("*denominator*")" : denominator)
end

function _rational_inverse(A::Matrix{_RExpr})
    n=size(A,1);B=[A [i==j ? _rone() : _rzero() for i in 1:n,j in 1:n]]
    for col in 1:n
        pivot=findfirst(r->!_riszero(B[r,col]),col:n)
        pivot===nothing && error("singular symbolic path matrix")
        pivot=col+pivot-1
        if pivot!=col saved=copy(B[col,:]);B[col,:]=B[pivot,:];B[pivot,:]=saved end
        p=B[col,col];B[col,:]=[_rdiv(x,p) for x in B[col,:]]
        for row in 1:n
            row==col && continue
            factor=B[row,col];_riszero(factor) && continue
            B[row,:]=[_rsub(B[row,j],_rmul(factor,B[col,j])) for j in axes(B,2)]
        end
    end
    B[:,n+1:2n]
end

function _rdet(A::Matrix{_RExpr})
    n=size(A,1);n==0 && return _rone();n==1 && return A[1,1]
    result=_rzero()
    for j in 1:n
        _riszero(A[1,j]) && continue
        minor=A[2:end,[q for q in 1:n if q!=j]]
        term=_rmul(A[1,j],_rdet(minor))
        result=isodd(j) ? _radd(result,term) : _rsub(result,term)
    end
    result
end

function _tree_face_edge_weights(E::PlabicEmbedding,names::Vector{String},directed)
    F=length(E.faces);length(names)==F || throw(ArgumentError("expected $F face weights"))
    F<=1 && return fill(_rone(),length(E.graph.edges))
    candidates=Tuple{Int,Int,Int,Int}[]
    for (from,to,edge_id) in directed
        left=get(E.halfedge_face,(from,to),0);right=get(E.halfedge_face,(to,from),0)
        left!=0 && right!=0 && left!=right && push!(candidates,(left,right,edge_id,1))
    end
    parent=zeros(Int,F);parent_edge=zeros(Int,F);parent_sign=zeros(Int,F);parent[1]=-1
    queue=[1]
    while !isempty(queue)
        a=popfirst!(queue)
        for (left,right,e,_) in candidates
            if left==a && parent[right]==0
                parent[right]=a;parent_edge[right]=e;parent_sign[right]=1;push!(queue,right)
            elseif right==a && parent[left]==0
                parent[left]=a;parent_edge[left]=e;parent_sign[left]=-1;push!(queue,left)
            end
        end
    end
    all(parent[2:end].!=0) || error("face-adjacency graph is disconnected")
    weights=fill(_rone(),length(E.graph.edges))
    # A tree edge carries the product of variables in the child-side subtree.
    children=[Int[] for _ in 1:F]
    for f in 2:F push!(children[parent[f]],f) end
    function subtree(f)
        result=[f]
        for c in children[f] append!(result,subtree(c)) end
        result
    end
    for f in 2:F
        factors=subtree(f)
        product=foldl(_rmul,(_rvariable(names[j]) for j in factors);init=_rone())
        weights[parent_edge[f]]=parent_sign[f]==1 ? product : _rdiv(_rone(),product)
    end
    return weights
end

function _symbolic_det(A::Matrix{String})
    n=size(A,1);n==0 && return "1";n==1 && return A[1,1]
    result="0"
    for j in 1:n
        A[1,j]=="0" && continue
        minor=A[2:end,[q for q in 1:n if q!=j]]
        term=_smul(A[1,j],_symbolic_det(minor))
        result=j%2==1 ? _sadd(result,term) : _ssub(result,term)
    end
    return result
end

function _plabic_boundary_measurement(G::PlabicGraph;mode=:edge,weights=Dict{Int,String}(),
                                      sources=nothing,
                                      iterations=2500,restarts=12)
    E=plabic_embedding(G;iterations=iterations,restarts=restarts)
    sources,matching=_perfect_orientation_matching(E.graph,sources)
    directed=_directed_plabic_edges(E.graph,matching)
    edge_weights=fill(_rone(),length(E.graph.edges))
    if mode==:edge
        for (i,value) in weights edge_weights[i]=_rvariable(value) end
    elseif mode==:face
        names=[get(weights,i,"1") for i in 1:length(E.faces)]
        edge_weights=_tree_face_edge_weights(E,names,directed)
    else
        throw(ArgumentError("weight mode must be face or edge"))
    end
    internal=collect((G.n+1):length(E.graph.colors));index=Dict(v=>i for (i,v) in enumerate(internal))
    T=fill(_rzero(),length(internal),length(internal))
    for (u,v,e) in directed
        u>G.n && v>G.n && (T[index[u],index[v]]=_radd(T[index[u],index[v]],edge_weights[e]))
    end
    resolvent=_rational_inverse([_rsub(i==j ? _rone() : _rzero(),T[i,j])
                                  for i in eachindex(internal),j in eachindex(internal)])
    A=fill(_rzero(),length(sources),G.n)
    for (row,source) in enumerate(sources)
        A[row,source]=_rone()
        q=fill(_rzero(),length(internal));out=Dict{Int,Vector{Tuple{Int,_RExpr}}}()
        for (u,v,e) in directed
            u==source && v>G.n && (q[index[v]]=_radd(q[index[v]],edge_weights[e]))
            u>G.n && v<=G.n && push!(get!(out,v,Tuple{Int,_RExpr}[]),(index[u],edge_weights[e]))
        end
        for target in 1:G.n
            target in sources && continue
            value=_rzero()
            for (a,qa) in enumerate(q), (b,w) in get(out,target,Tuple{Int,_RExpr}[])
                value=_radd(value,_rmul(_rmul(qa,resolvent[a,b]),w))
            end
            lo,hi=minmax(source,target)
            between=count(s->lo<s<hi,sources)
            A[row,target]=isodd(between) ? _rneg(value) : value
        end
    end
    pluckers=Dict{String,_RExpr}()
    for I in _subsets(collect(1:G.n),length(sources))
        pluckers[join(I,",")]=_rdet(A[:,I])
    end
    return (matrix=_rstring.(A),matrix_latex=_rstring.(A;latex=true),
            pluckers=Dict(k=>_rstring(v) for (k,v) in pluckers),
            pluckers_latex=Dict(k=>_rstring(v;latex=true) for (k,v) in pluckers),
            sources=sources,edge_weights=_rstring.(edge_weights))
end

function _measurement_payload(G,request_body;iterations=2500,restarts=12)
    input=JSON.parse(String(request_body))
    mode=Symbol(get(input,"mode","face"))
    sources=haskey(input,"sources") ? Int.(input["sources"]) : nothing
    raw=get(input,"weights",Dict{String,Any}())
    weights=Dict(parse(Int,string(k))=>string(v) for (k,v) in raw)
    R=_plabic_boundary_measurement(G;mode=mode,weights=weights,sources=sources,
                                   iterations=iterations,restarts=restarts)
    matrix_julia="["*join((join(("("*entry*")" for entry in R.matrix[i,:])," ")
                             for i in axes(R.matrix,1)),"; ")*"]"
    matrix_m2="matrix{"*join(("{"*join(R.matrix[i,:],",")*"}"
                               for i in axes(R.matrix,1)),",")*"}"
    latex_matrix=raw"\begin{pmatrix}"*
                 join((join(R.matrix_latex[i,:]," & ")
                       for i in axes(R.matrix,1)),raw"\\\\")*raw"\end{pmatrix}"
    ordered=sort(collect(R.pluckers);by=x->parse.(Int,split(first(x),',')))
    plucker_julia="Dict("*join(("("*I*") => "*v for (I,v) in ordered),", ")*")"
    plucker_m2="{"*join(("p_"*replace(I,","=>"_")*" => "*v for (I,v) in ordered),",")*"}"
    latex_plucker=raw"\begin{aligned}"*
                    join((raw"\Delta_{"*replace(I,","=>"")*"} &= "*
                          R.pluckers_latex[I] for (I,v) in ordered),raw"\\\\")*
                    raw"\end{aligned}"
    return JSON.json(Dict("sources"=>R.sources,"matrix"=>R.matrix,
                          "pluckers"=>R.pluckers,"matrix_latex"=>latex_matrix,
                          "plucker_latex"=>latex_plucker,
                          "matrix_julia"=>matrix_julia,"matrix_m2"=>matrix_m2,
                          "plucker_julia"=>plucker_julia,"plucker_m2"=>plucker_m2))
end

function _interactive_state_json(G;iterations=2500,restarts=12)
    E,labels,movable=_interactive_graph_state(G;iterations=iterations,restarts=restarts)
    default_sources=Int.(minGrassmannNecklace(G.permutation)[1])
    vertices=join((string("{\"id\":",i,",\"x\":",_json_number(x),
                          ",\"y\":",_json_number(y),",\"color\":\"",
                          E.graph.colors[i],"\"}")
                     for (i,(x,y)) in enumerate(E.graph.positions)),",")
    edges=join((string("[",u,",",v,"]") for (u,v) in E.graph.edges),",")
    boundary_positions=_boundary_dual_positions(E)
    faces=String[]
    for (face_id,(face,label,can_move)) in enumerate(zip(E.faces,labels,movable))
        text="{"*join(label,",")*"}"
        (x,y),font_size=_face_label_layout(face,E.graph.positions,text;
                                             max_size=16,pixels_per_unit=300)
        polygon=join((string("[",_json_number(E.graph.positions[v][1]),",",
                              _json_number(E.graph.positions[v][2]),"]")
                      for v in face),",")
        boundary_face=haskey(boundary_positions,face_id)
        dual_x,dual_y=boundary_face ? boundary_positions[face_id] : (x,y)
        push!(faces,string("{\"id\":",face_id,",\"label\":",_json_ints(label),
                           ",\"x\":",_json_number(x),",\"y\":",_json_number(y),
                           ",\"dual_x\":",_json_number(dual_x),
                           ",\"dual_y\":",_json_number(dual_y),
                           ",\"boundary\":",boundary_face,
                           ",\"font_size\":",_json_number(font_size),
                           ",\"polygon\":[",polygon,"]",
                           ",\"movable\":",can_move,"}"))
    end
    dual_edge_set=Set{Tuple{Int,Int}}()
    for (u,v) in E.graph.edges
        a=get(E.halfedge_face,(u,v),0)
        b=get(E.halfedge_face,(v,u),0)
        a!=0 && b!=0 && a!=b && push!(dual_edge_set,(min(a,b),max(a,b)))
    end
    dual_edges=join((string("[",a,",",b,"]")
                     for (a,b) in sort(collect(dual_edge_set))),",")
    black_dual_faces=join((string("[",join(face_ids,","),"]")
                           for face_ids in _black_dual_faces(E)),",")
    strands=String[]
    for i in 1:G.n
        points=_medial_strand_points(E.trips[i],E.graph.positions,
                                     E.graph.colors,E.graph.n)
        encoded=join((string("[",_json_number(x),",",_json_number(y),"]")
                      for (x,y) in points),",")
        push!(strands,string("{\"source\":",i,",\"points\":[",encoded,"]}"))
    end
    return string("{\"n\":",G.n,",\"k\":",length(default_sources),
                  ",\"default_sources\":",_json_ints(default_sources),
                  ",\"permutation\":",_json_ints(G.permutation),
                  ",\"vertices\":[",vertices,
                  "],\"edges\":[",edges,"],\"faces\":[",join(faces,","),
                  "],\"dual_edges\":[",dual_edges,
                  "],\"dual_black_faces\":[",black_dual_faces,
                  "],\"strands\":[",join(strands,","),"]}")
end

_interactive_blank_state_json() =
    "{\"blank\":true,\"n\":0,\"k\":0,\"default_sources\":[],\"permutation\":[],\"vertices\":[],\"edges\":[],\"faces\":[],\"dual_edges\":[],\"dual_black_faces\":[],\"strands\":[]}"

function _interactive_animation_frame_json(caption::AbstractString,G::PlabicGraph)
    vertices=join((string("{\"id\":",i,",\"x\":",_json_number(x),
                          ",\"y\":",_json_number(y),",\"color\":\"",
                          G.colors[i],"\"}")
                     for (i,(x,y)) in enumerate(G.positions)),",")
    edges=join((string("[",u,",",v,"]") for (u,v) in G.edges),",")
    return string("{\"caption\":",JSON.json(String(caption)),
                  ",\"n\":",G.n,",\"vertices\":[",vertices,
                  "],\"edges\":[",edges,"]}")
end

function _interactive_animation_json(stages,state)
    frames=join((_interactive_animation_frame_json(caption,G)
                 for (caption,G) in stages),",")
    return string("{\"frames\":[",frames,"],\"state\":",state,"}")
end

function _interactive_history_json(state::AbstractString,can_undo::Bool,can_redo::Bool=false)
    endswith(state,"}") || error("invalid interactive state")
    return string(state[begin:prevind(state,lastindex(state))],
                  ",\"can_undo\":",can_undo,
                  ",\"can_redo\":",can_redo,"}")
end

const _INTERACTIVE_PLABIC_HTML = raw"""<!doctype html>
<html><head><meta charset="utf-8"><title>Interactive plabic graph</title>
<script>window.MathJax={tex:{inlineMath:[['\\(','\\)']],displayMath:[['\\[','\\]']]},svg:{fontCache:'global'}};</script><script async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
<style>
html,body{margin:0;height:100%;background:#fafafa;font-family:system-ui,sans-serif}
#wrap{height:100%;display:grid;grid-template-rows:auto 1fr}.bar{padding:9px 14px;background:white;border-bottom:1px solid #ddd;display:flex;align-items:center;gap:8px;flex-wrap:wrap}
#main{--side-width:350px;min-height:0;display:grid;grid-template-columns:minmax(0,1fr) var(--side-width)}#workspace{min-width:0;min-height:0;overflow:auto;scrollbar-gutter:stable}#canvas{min-width:0;height:max(620px,calc(100vh - 58px));position:relative}#side{overflow:auto;background:white;border-left:1px solid #ddd;padding:14px}#side.collapsed{padding:8px 5px;overflow:hidden}#side-content[hidden]{display:none}#side h3{margin:14px 0 8px}#side h3:first-child{margin-top:0}.panel-size-controls{position:sticky;top:0;z-index:4;background:white;border-bottom:1px solid #ddd;padding-bottom:9px;margin-bottom:9px}.panel-size-controls label{display:grid;grid-template-columns:auto minmax(80px,1fr) 42px;gap:6px;align-items:center;font-size:12px}.panel-size-controls input{width:auto;min-width:0}.face-item{font-family:ui-monospace,monospace;padding:3px 5px;user-select:text}.square-item{color:#d7191c;font-weight:700}.side-note{font-size:12px;color:#666;margin:5px 0 10px}.controls,.measurement-controls,.copy-controls{display:flex;gap:6px;align-items:center;flex-wrap:wrap}input{width:230px;padding:5px 7px;font-family:ui-monospace,monospace}button,select{padding:5px 9px;cursor:pointer}#source-set{width:90px}#status{color:#555;margin-left:4px}.weight-row{display:grid;grid-template-columns:minmax(0,1fr) 120px;gap:6px;align-items:center;margin:4px 0;font-family:ui-monospace,monospace;font-size:12px}.weight-row input{width:auto;min-width:0}#measurement-panel{background:white;border-top:1px solid #ccc;padding:12px 18px;min-height:120px;overflow:auto}.measurement-output{overflow:auto;margin-top:8px;padding:14px;background:#fafafa;border:1px solid #ddd;min-height:54px;font-size:11px;text-align:center}svg{width:100%;height:100%;touch-action:none}.edge{stroke:#222;stroke-width:2.5}.dual-black-face{fill:#30343b;fill-opacity:.72;stroke:none;pointer-events:none}.dual-edge{stroke:#357a8a;stroke-width:2.5;stroke-dasharray:8 6;opacity:.9;pointer-events:none}.dual-label{fill:#d7191c;text-anchor:middle;dominant-baseline:central;font-family:ui-monospace,monospace;font-weight:650}.dual-movable{pointer-events:all;cursor:pointer;font-weight:700}
.internal{stroke:#111;stroke-width:1.5}.boundary{fill:#111;cursor:pointer}.boundary:hover{fill:#0066cc}.boundary-label{font-size:16px;cursor:pointer;user-select:none}.face-label{fill:#d7191c;text-anchor:middle;dominant-baseline:central;user-select:none}
.square-face{fill:#cfeeff;stroke:none}.movable{fill:#d7191c;cursor:pointer;font-weight:700}.movable:hover{text-decoration:underline}.strand{fill:none;stroke-width:3;stroke-linecap:round;stroke-linejoin:round;opacity:.9;pointer-events:none}.hint{font-size:14px;color:#333;margin-left:12px}
.measurement-output mjx-container svg{width:auto!important;height:auto!important;max-width:none}
.edge-id{fill:#24506b;font-size:12px;font-weight:700;text-anchor:middle;dominant-baseline:central;paint-order:stroke;stroke:white;stroke-width:4px;stroke-linejoin:round;pointer-events:none}
#side.collapsed .panel-size-controls{border:0;margin:0;padding:0}#side.collapsed .panel-size-controls label{display:block;font-size:0}#side.collapsed #side-size{display:block;width:28px;height:160px;padding:0;writing-mode:vertical-lr;direction:rtl}#side.collapsed #side-size-value{display:none}
#main{grid-template-columns:minmax(0,1fr) 7px var(--side-width)}.panel-size-controls{display:none}#side-dragger{background:#e3e3e3;cursor:col-resize;touch-action:none;position:relative;z-index:6}#side-dragger:hover,#side-dragger.dragging{background:#8bb7cc}#side-dragger:after{content:"";position:absolute;inset:0 -4px}#side.collapsed{padding:0;border-left:0}
html,body{overflow:hidden}.bar>b{white-space:nowrap}.bar>.controls{min-width:0;flex:1}.bar button{padding:5px 7px;font-size:12px}.bar input{width:min(230px,22vw)}
.builder-panel{margin-bottom:18px;padding-bottom:16px;border-bottom:1px solid var(--line)}.builder-panel[hidden]{display:none}.builder-grid{display:grid;grid-template-columns:1fr 1fr;gap:6px;margin:8px 0}.builder-grid button{font-size:12px}.builder-grid button.active{background:var(--accent);border-color:var(--accent);color:white}.builder-panel input{width:65px}.builder-actions{display:flex;gap:6px;flex-wrap:wrap}.builder-vertex{cursor:pointer;stroke-width:2}.builder-selected{stroke:#ef9b27!important;stroke-width:5!important}.builder-edge{stroke:#44545b;stroke-width:3;cursor:pointer}.builder-edge:hover{stroke:var(--danger);stroke-width:5}.builder-boundary{fill:#34434a;cursor:pointer}.builder-number{font-size:16px;font-weight:600;pointer-events:none}.builder-hint{fill:#708087;font-size:18px;text-anchor:middle}.builder-permutation{font-family:ui-monospace,monospace;color:var(--accent-dark);font-weight:650;margin-top:8px}
.animation-panel{padding-bottom:14px;margin-bottom:10px;border-bottom:1px solid var(--line)}.animation-panel label{display:flex;align-items:center;gap:8px;color:var(--muted-ink);font-size:12px;font-weight:600}.animation-panel select{flex:1}.animation-caption{position:absolute;left:50%;bottom:24px;transform:translateX(-50%);z-index:7;padding:7px 13px;border-radius:18px;background:rgba(38,50,56,.88);color:white;font-size:13px;box-shadow:0 4px 14px rgba(0,0,0,.16);pointer-events:none}.animation-caption[hidden]{display:none}.animation-frame .edge{stroke:#4b5f67}.animation-frame .internal{filter:drop-shadow(0 2px 3px rgba(0,0,0,.16))}
.face-hit{fill:transparent;stroke:none;pointer-events:all}.edge-hit{stroke:transparent;stroke-width:16;pointer-events:stroke;cursor:context-menu}.parameter-label{fill:#176f86;font-family:ui-monospace,monospace;font-size:15px;font-weight:700;text-anchor:middle;dominant-baseline:central;paint-order:stroke;stroke:white;stroke-width:4px;stroke-linejoin:round;pointer-events:none}.parameter-menu{position:fixed;z-index:1000;width:220px;padding:10px;background:white;border:1px solid #c9d5d9;border-radius:9px;box-shadow:0 10px 30px rgba(25,43,50,.2)}.parameter-menu[hidden]{display:none}.parameter-menu b{display:block;font-size:12px;margin-bottom:7px}.parameter-menu input{box-sizing:border-box;width:100%;margin-bottom:7px}.parameter-menu-actions{display:flex;gap:6px}.parameter-menu-note{font-size:11px;color:var(--muted-ink);margin-bottom:7px}
:root{--ink:#263238;--muted-ink:#65747b;--line:#dce3e6;--soft:#f4f7f8;--paper:#fff;--accent:#287b91;--accent-dark:#1f6477;--danger:#c43d4b}html,body{color:var(--ink);background:var(--soft)}#wrap{background:linear-gradient(135deg,#f7fafb 0%,#eef3f4 100%)}.bar{padding:11px 16px;gap:10px;background:rgba(255,255,255,.94);border-color:var(--line);box-shadow:0 2px 12px rgba(28,50,58,.06);z-index:8}.bar>b{font-size:15px;letter-spacing:.01em;margin-right:5px}.controls{gap:7px}.controls label,.measurement-controls label{color:var(--muted-ink);font-size:12px;font-weight:600}button,select,input{border:1px solid #cbd6da;border-radius:7px;background:var(--paper);color:var(--ink);transition:border-color .15s,box-shadow .15s,background .15s,transform .08s}button:hover,select:hover,input:hover{border-color:#91aeb8}button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid rgba(40,123,145,.28);outline-offset:1px;border-color:var(--accent)}button:active{transform:translateY(1px)}button:disabled{opacity:.45;cursor:not-allowed}#draw{background:var(--accent);border-color:var(--accent);color:white;font-weight:650}#draw:hover{background:var(--accent-dark)}#back,#forward{background:#f8fafb}#status{font-size:12px;color:var(--muted-ink)}#workspace{background:radial-gradient(circle at 50% 38%,#fff 0%,#f7f9fa 68%,#eef2f3 100%)}#canvas{padding:10px;box-sizing:border-box}#side{background:rgba(255,255,255,.97);border-color:var(--line);padding:18px 16px}#side h3{font-size:13px;letter-spacing:.035em;text-transform:uppercase;color:#49616a;margin:22px 0 10px;padding-bottom:7px;border-bottom:1px solid #e8edef}#side h3:first-child{margin-top:2px}.side-note{color:var(--muted-ink);line-height:1.45}.face-item{border-radius:5px;padding:4px 7px}.face-item:nth-child(odd){background:#f7f9fa}.square-item{color:var(--danger)}.weight-row{padding:3px 0;gap:8px}.weight-row span{overflow:hidden;text-overflow:ellipsis}.weight-row input{background:#fbfcfc}#side-dragger{background:#dfe6e8}#side-dragger:before{content:"";position:absolute;top:50%;left:2px;width:3px;height:42px;transform:translateY(-50%);border-radius:3px;background:#a8b7bc}#side-dragger:hover,#side-dragger.dragging{background:#c9e1e8}#side-dragger:hover:before,#side-dragger.dragging:before{background:var(--accent)}#measurement-panel{box-sizing:border-box;background:rgba(255,255,255,.96);border-color:var(--line);padding:18px 22px}#measurement-panel b{font-size:14px}.measurement-output{background:linear-gradient(180deg,#fff,#fbfcfc);border-color:var(--line);border-radius:10px;box-shadow:inset 0 1px 2px rgba(29,50,57,.035);padding:20px;color:#202b30}.disk{fill:#fff;stroke:#aab8bd;stroke-width:1.8;filter:drop-shadow(0 8px 14px rgba(35,55,62,.10))}.edge{stroke:#34434a;stroke-width:2.35}.internal{stroke:#26343a;stroke-width:1.4}.boundary{fill:#34434a}.boundary:hover{fill:var(--accent)}.boundary-label{fill:#34434a;font-weight:600}.face-label{fill:#c53645}.square-face{fill:#d9eef5}.dual-edge{stroke:#31869b;opacity:.8}.edge-id{fill:#24687b}.strand{filter:drop-shadow(0 1px 1px rgba(0,0,0,.12))}@media(max-width:760px){.bar{padding:9px}.bar>b{width:100%}#main{--side-width:280px}.controls input{width:170px}#canvas{height:620px}#side{padding:12px 10px}}
</style></head><body><div id="wrap"><div class="bar"><b>Interactive plabic graph</b><div class="controls"><label for="permutation">Permutation</label><input id="permutation" spellcheck="false"><button id="draw">Draw / replace graph</button><button id="back" disabled>Back</button><button id="forward" disabled>Forward</button><button id="all-strands">Draw all strands</button><button id="primal-toggle">Hide plabic graph</button><button id="dual-toggle">Show dual graph</button></div><span id="status"></span></div><div id="main"><div id="workspace"><div id="canvas"><svg id="graph" viewBox="0 0 800 800"></svg></div><section id="measurement-panel"><div class="measurement-controls"><b>Boundary measurement result</b><button id="copy-julia" disabled>Copy Julia</button><button id="copy-m2" disabled>Copy Macaulay2</button></div><div id="measurement-result" class="measurement-output">Assign variables, then press Compute.</div></section></div><aside id="side"><div class="panel-size-controls"><label>Panel <input id="side-size" type="range" min="0" max="600" step="10" value="350"><span id="side-size-value">350</span></label></div><div id="side-content"><h3>Face labels</h3><button id="copy-labels">Copy all</button><div class="side-note">Square-face labels are red. Each line can also be selected and copied normally.</div><div id="face-list"></div><h3>Boundary measurement</h3><div class="measurement-controls"><select id="weight-mode"><option value="face">Face variables</option><option value="edge">Edge variables</option></select><label>Sources <input id="source-set" spellcheck="false"></label><button id="assign-edge-variables" disabled>Assign t₁,t₂,… to all edges</button><button id="compute-measurement">Compute</button></div><div class="side-note">Choose a source (k\)-subset such as 1,2. Blank edge weights equal 1. Face variables use Postnikov's clockwise convention; one face is dependent through the product-one relation.</div><div id="weight-list"></div><h3>Result display</h3><div class="measurement-controls"><label>Show <select id="result-mode"><option value="matrix">Matrix</option><option value="pluckers">Plücker coordinates</option></select></label><label>Text size <input id="matrix-size" type="range" min="7" max="22" step="1" value="11"><span id="matrix-size-value">11</span></label><label>Result height <input id="result-height" type="range" min="120" max="700" step="20" value="300"><span id="result-height-value">300</span></label></div></div></aside></div></div>
<script>
const svg=document.getElementById('graph'),status=document.getElementById('status'),faceList=document.getElementById('face-list'),permInput=document.getElementById('permutation'),backButton=document.getElementById('back'),forwardButton=document.getElementById('forward'),allStrandsButton=document.getElementById('all-strands'),primalButton=document.getElementById('primal-toggle'),dualButton=document.getElementById('dual-toggle'),weightMode=document.getElementById('weight-mode'),sourceInput=document.getElementById('source-set'),assignEdgeVariablesButton=document.getElementById('assign-edge-variables'),resultMode=document.getElementById('result-mode'),weightList=document.getElementById('weight-list'),measurementResult=document.getElementById('measurement-result'),copyJulia=document.getElementById('copy-julia'),copyM2=document.getElementById('copy-m2');
const NS='http://www.w3.org/2000/svg',X=x=>400+300*x,Y=y=>400-300*y;
const selectedStrands=new Set(),palette=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#17a2b8','#6a3d9a','#1b9e77'];
let currentState=null,primalVisible=true,dualVisible=false,labelsVisible=true,weightSignature='',measurementData=null;const faceWeights=new Map(),edgeWeights=new Map();
const main=document.getElementById('main'),side=document.getElementById('side'),sideContent=document.getElementById('side-content'),sideSize=document.getElementById('side-size'),sideSizeValue=document.getElementById('side-size-value'),matrixSize=document.getElementById('matrix-size'),matrixSizeValue=document.getElementById('matrix-size-value'),resultHeight=document.getElementById('result-height'),resultHeightValue=document.getElementById('result-height-value'),measurementPanel=document.getElementById('measurement-panel');
function setSideSize(){const value=Number(sideSize.value);main.style.setProperty('--side-width',value+'px');side.classList.toggle('collapsed',value===0);sideContent.hidden=value===0;sideSizeValue.textContent=value===0?'hidden':value}
function setResultHeight(){measurementPanel.style.height=resultHeight.value+'px';resultHeightValue.textContent=resultHeight.value}
const sideDragger=document.createElement('div');sideDragger.id='side-dragger';sideDragger.title='Drag to resize the control panel';main.insertBefore(sideDragger,side);let resizingSide=false;sideDragger.addEventListener('pointerdown',e=>{resizingSide=true;sideDragger.classList.add('dragging');sideDragger.setPointerCapture(e.pointerId);e.preventDefault()});sideDragger.addEventListener('pointermove',e=>{if(!resizingSide)return;sideSize.value=String(Math.max(0,Math.min(600,window.innerWidth-e.clientX)));setSideSize()});function stopSideResize(){resizingSide=false;sideDragger.classList.remove('dragging')}sideDragger.addEventListener('pointerup',stopSideResize);sideDragger.addEventListener('pointercancel',stopSideResize);sideDragger.addEventListener('dblclick',()=>{sideSize.value=side.classList.contains('collapsed')?'350':'0';setSideSize()});
measurementResult.style.fontSize=matrixSize.value+'px';matrixSize.addEventListener('input',()=>{measurementResult.style.fontSize=matrixSize.value+'px';matrixSizeValue.textContent=matrixSize.value;if(measurementData)renderMeasurementResult()});resultHeight.addEventListener('input',setResultHeight);setSideSize();setResultHeight();
function el(name,attrs,text){const e=document.createElementNS(NS,name);for(const [k,v] of Object.entries(attrs||{}))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;return e}
const toolbarControls=document.querySelector('.bar .controls'),manualButton=document.createElement('button');manualButton.id='manual-draw';manualButton.textContent='Draw graph manually';toolbarControls.append(manualButton);const editGraphButton=document.createElement('button');editGraphButton.id='edit-current-graph';editGraphButton.textContent='Edit current graph';editGraphButton.disabled=true;toolbarControls.append(editGraphButton);
const faceLabelsButton=document.createElement('button');faceLabelsButton.id='face-label-toggle';faceLabelsButton.textContent='Hide face labels';toolbarControls.append(faceLabelsButton);
const parameterMenu=document.createElement('div');parameterMenu.className='parameter-menu';parameterMenu.hidden=true;parameterMenu.innerHTML='<b id="parameter-menu-title">Assign parameter</b><div id="parameter-menu-note" class="parameter-menu-note"></div><input id="parameter-name" spellcheck="false" placeholder="a or x_1"><div class="parameter-menu-actions"><button id="parameter-assign">Assign</button><button id="parameter-clear">Clear</button><button id="parameter-cancel">Cancel</button></div>';document.body.append(parameterMenu);const parameterName=document.getElementById('parameter-name'),parameterTitle=document.getElementById('parameter-menu-title'),parameterNote=document.getElementById('parameter-menu-note'),parameterAssign=document.getElementById('parameter-assign');let parameterTarget=null;
function closeParameterMenu(){parameterMenu.hidden=true;parameterTarget=null}
function openParameterMenu(kind,id,event){event.preventDefault();event.stopPropagation();parameterTarget={kind,id};const map=kind==='face'?faceWeights:edgeWeights,dependent=kind==='face'&&id===1;parameterTitle.textContent='Assign parameter to '+kind+' '+id;parameterNote.textContent=dependent?'This reference face is dependent and has no independent parameter.':'Unassigned weights are 1.';parameterName.value=map.get(id)||'';parameterName.disabled=dependent;parameterAssign.disabled=dependent;parameterMenu.hidden=false;const width=240,height=145;parameterMenu.style.left=Math.max(8,Math.min(window.innerWidth-width,event.clientX))+'px';parameterMenu.style.top=Math.max(8,Math.min(window.innerHeight-height,event.clientY))+'px';dependent||setTimeout(()=>parameterName.focus(),0)}
function saveParameter(clear=false){if(!parameterTarget)return;const map=parameterTarget.kind==='face'?faceWeights:edgeWeights,value=clear?'':parameterName.value.trim();value&&value!=='1'?map.set(parameterTarget.id,value):map.delete(parameterTarget.id);measurementData=null;const description=parameterTarget.kind+' '+parameterTarget.id;closeParameterMenu();if(currentState)render(currentState);status.textContent=value?'Assigned '+value+' to '+description+'.':'Cleared '+description+' parameter.'}
parameterAssign.addEventListener('click',()=>saveParameter(false));document.getElementById('parameter-clear').addEventListener('click',()=>saveParameter(true));document.getElementById('parameter-cancel').addEventListener('click',closeParameterMenu);parameterName.addEventListener('keydown',e=>{if(e.key==='Enter'&&!parameterAssign.disabled)saveParameter(false);if(e.key==='Escape')closeParameterMenu()});document.addEventListener('pointerdown',e=>{if(!parameterMenu.hidden&&!parameterMenu.contains(e.target))closeParameterMenu()});faceLabelsButton.addEventListener('click',()=>{labelsVisible=!labelsVisible;faceLabelsButton.textContent=labelsVisible?'Hide face labels':'Show face labels';if(currentState)render(currentState)});
const builderPanel=document.createElement('section');builderPanel.className='builder-panel';builderPanel.hidden=true;builderPanel.innerHTML='<h3>Graph drawing</h3><label>Boundary vertices <input id="builder-n" type="number" min="1" max="30" value="6"></label><div class="builder-grid"><button data-builder-tool="black">Black vertex</button><button data-builder-tool="white">White vertex</button><button data-builder-tool="edge">Add edges</button><button data-builder-tool="delete">Delete</button></div><div class="builder-actions"><button id="builder-undo">Undo</button><button id="builder-clear">Clear</button><button id="builder-cancel">Cancel</button><button id="builder-finish">Reduce & use graph</button></div><div class="side-note">Place internal vertices, then drag them whenever you want to adjust the drawing. In Add edges mode, click two endpoints. Delete removes a selected vertex or edge. Boundary vertices must each have one edge.</div><div id="builder-result" class="builder-permutation"></div>';
sideContent.insertBefore(builderPanel,sideContent.firstChild);const builderN=document.getElementById('builder-n'),builderResult=document.getElementById('builder-result');
const animationPanel=document.createElement('section');animationPanel.className='animation-panel';animationPanel.innerHTML='<h3>Animation</h3><label>Transformation speed <select id="animation-speed"><option value="0">Instant</option><option value="350">Fast</option><option value="900" selected>Slow</option><option value="1600">Very slow</option><option value="2800">Movie mode</option></select></label><div class="side-note">Square moves and hand-drawn graph reductions continuously morph through their individual operations.</div>';builderPanel.after(animationPanel);const animationSpeed=document.getElementById('animation-speed'),animationCaption=document.createElement('div');animationCaption.className='animation-caption';animationCaption.hidden=true;document.getElementById('canvas').append(animationCaption);let mutationPlaying=false;
function renderAnimationLabels(labels){const layer=el('g',{class:'animation-labels'});for(const f of labels||[])layer.append(el('text',{x:X(f.x),y:Y(f.y),'font-size':f.font_size,class:'face-label'},labelText(f.label)));svg.append(layer)}
function renderAnimationFrame(frame,labels=[]){svg.replaceChildren();svg.append(el('circle',{cx:400,cy:400,r:300,class:'disk'}));const layer=el('g',{class:'animation-frame'}),by=new Map(frame.vertices.map(v=>[v.id,v]));svg.append(layer);for(const [a,b] of frame.edges){const u=by.get(a),v=by.get(b);if(u&&v)layer.append(el('line',{x1:X(u.x),y1:Y(u.y),x2:X(v.x),y2:Y(v.y),class:'edge'}))}for(const v of frame.vertices){if(v.id<=frame.n)continue;layer.append(el('circle',{cx:X(v.x),cy:Y(v.y),r:8,fill:v.color==='black'?'black':'white',class:'internal'}))}for(let i=1;i<=frame.n;i++){const v=by.get(i);if(!v)continue;layer.append(el('circle',{cx:X(v.x),cy:Y(v.y),r:5,class:'boundary'}));layer.append(el('text',{x:X(1.12*v.x),y:Y(1.12*v.y),class:'boundary-label','text-anchor':'middle','dominant-baseline':'central'},String(i)))}renderAnimationLabels(labels)}
function stateAnimationFrame(state,caption=''){return {caption,n:state.n,vertices:state.vertices.map(v=>({id:v.id,x:v.x,y:v.y,color:v.color})),edges:state.edges.map(e=>[e[0],e[1]])}}
function nearestVertex(vertex,candidates){let best=candidates[0],score=Infinity;for(const candidate of candidates){const dx=vertex.x-candidate.x,dy=vertex.y-candidate.y,colorPenalty=vertex.color===candidate.color?0:.025,value=dx*dx+dy*dy+colorPenalty;if(value<score){score=value;best=candidate}}return best}
function morphCorrespondence(from,to){
 const targets=new Map(to.vertices.map(v=>[v.id,v])),sources=new Map(from.vertices.map(v=>[v.id,v])),forward=new Map(),reverse=new Map(),samePoint=(a,b)=>Math.hypot(a.x-b.x,a.y-b.y)<1e-8;
 const contraction=String(to.caption||'').match(/^Contract same-color edge (\d+)[–-](\d+)/);
 if(contraction){const kept=Math.min(Number(contraction[1]),Number(contraction[2])),removed=Math.max(Number(contraction[1]),Number(contraction[2]));for(const vertex of from.vertices){const targetId=vertex.id===removed?kept:(vertex.id>removed?vertex.id-1:vertex.id);forward.set(vertex.id,targets.get(targetId)||vertex)}for(const target of to.vertices){const sourceId=target.id===kept?kept:(target.id>=removed?target.id+1:target.id);reverse.set(target.id,sources.get(sourceId)||target)}return {forward,reverse}}
 if(String(to.caption||'').startsWith('Suppress one bivalent vertex')&&from.vertices.length===to.vertices.length+1){let removed=null;for(let candidate=from.n+1;candidate<=from.vertices.length;candidate++){const remaining=from.vertices.filter(v=>v.id!==candidate);if(remaining.every((v,j)=>to.vertices[j]&&v.color===to.vertices[j].color&&samePoint(v,to.vertices[j]))){removed=candidate;break}}if(removed!==null){const remap=id=>id>removed?id-1:id,neighbors=from.edges.flatMap(([a,b])=>a===removed?[b]:b===removed?[a]:[]),internal=neighbors.filter(id=>id>from.n),anchor=(internal.length?internal:neighbors).sort((a,b)=>{const u=sources.get(a),v=sources.get(b),r=sources.get(removed);return Math.hypot(u.x-r.x,u.y-r.y)-Math.hypot(v.x-r.x,v.y-r.y)})[0],anchorTarget=targets.get(remap(anchor));for(const vertex of from.vertices)forward.set(vertex.id,vertex.id===removed?anchorTarget:(targets.get(remap(vertex.id))||vertex));for(const target of to.vertices){const sourceId=target.id>=removed?target.id+1:target.id;reverse.set(target.id,sources.get(sourceId)||target)}return {forward,reverse}}}
 const sameSize=from.vertices.length===to.vertices.length;
 for(const vertex of from.vertices){let target;if(vertex.id<=from.n)target=targets.get(vertex.id);else if(sameSize&&targets.has(vertex.id))target=targets.get(vertex.id);else target=nearestVertex(vertex,to.vertices.filter(v=>v.id>to.n));forward.set(vertex.id,target||vertex)}
 for(const target of to.vertices){let source;if(target.id<=to.n)source=sources.get(target.id);else source=nearestVertex(target,from.vertices.filter(v=>v.id>from.n));reverse.set(target.id,source||target)}return {forward,reverse}
}
function drawMorph(from,to,t,labels=[]){
 const ease=t*t*(3-2*t),maps=morphCorrespondence(from,to),fromBy=new Map(from.vertices.map(v=>[v.id,v])),toBy=new Map(to.vertices.map(v=>[v.id,v])),edgeKey=(a,b)=>a<b?a+'-'+b:b+'-'+a,point=(a,b,s)=>({x:a.x+(b.x-a.x)*s,y:a.y+(b.y-a.y)*s});
 svg.replaceChildren();svg.append(el('circle',{cx:400,cy:400,r:300,class:'disk'}));const edgesLayer=el('g',{class:'animation-frame'}),nodes=el('g');svg.append(edgesLayer);svg.append(nodes);
 const movedOld=new Map(),movedNew=new Map(),claimedTargets=new Set();
 for(const v of from.vertices){const q=maps.forward.get(v.id)||v;movedOld.set(v.id,point(v,q,ease));if(q.id<=to.vertices.length)claimedTargets.add(q.id)}
 for(const q of to.vertices){const p=maps.reverse.get(q.id)||q;movedNew.set(q.id,point(p,q,ease))}
 const sourceEdges=new Map(from.edges.map(([a,b],index)=>[edgeKey(a,b),index])),targetEdges=new Map(to.edges.map(([a,b],index)=>[edgeKey(a,b),index])),usedSources=new Set(),usedTargets=new Set(),neededTargetNodes=new Set();
 const drawEdge=(u,v)=>edgesLayer.append(el('line',{x1:X(u.x),y1:Y(u.y),x2:X(v.x),y2:Y(v.y),class:'edge'}));
 for(const [sourceIndex,[a,b]] of from.edges.entries()){const targetA=maps.forward.get(a)||fromBy.get(a),targetB=maps.forward.get(b)||fromBy.get(b);if(targetA.id===targetB.id)continue;const targetIndex=targetEdges.get(edgeKey(targetA.id,targetB.id));if(targetIndex===undefined||usedTargets.has(targetIndex))continue;drawEdge(movedOld.get(a),movedOld.get(b));usedSources.add(sourceIndex);usedTargets.add(targetIndex)}
 for(const [targetIndex,[a,b]] of to.edges.entries()){if(usedTargets.has(targetIndex))continue;const sourceA=maps.reverse.get(a)||toBy.get(a),sourceB=maps.reverse.get(b)||toBy.get(b),sourceIndex=sourceEdges.get(edgeKey(sourceA.id,sourceB.id));if(sourceIndex===undefined||usedSources.has(sourceIndex))continue;drawEdge(movedNew.get(a),movedNew.get(b));neededTargetNodes.add(a);neededTargetNodes.add(b);usedSources.add(sourceIndex);usedTargets.add(targetIndex)}
 for(const [sourceIndex,[a,b]] of from.edges.entries()){if(usedSources.has(sourceIndex))continue;drawEdge(movedOld.get(a),movedOld.get(b))}
 for(const [targetIndex,[a,b]] of to.edges.entries()){if(usedTargets.has(targetIndex))continue;drawEdge(movedNew.get(a),movedNew.get(b));neededTargetNodes.add(a);neededTargetNodes.add(b)}
 for(const v of from.vertices){if(v.id<=from.n)continue;const q=maps.forward.get(v.id)||v,p=movedOld.get(v.id),color=ease<.5?v.color:q.color;nodes.append(el('circle',{cx:X(p.x),cy:Y(p.y),r:8,fill:color==='black'?'black':'white',class:'internal'}))}
 for(const q of to.vertices){if(q.id<=to.n||claimedTargets.has(q.id)&&!neededTargetNodes.has(q.id))continue;const v=movedNew.get(q.id);nodes.append(el('circle',{cx:X(v.x),cy:Y(v.y),r:8,fill:q.color==='black'?'black':'white',class:'internal'}))}
 for(let i=1;i<=to.n;i++){const q=movedNew.get(i)||toBy.get(i);nodes.append(el('circle',{cx:X(q.x),cy:Y(q.y),r:5,class:'boundary'}));nodes.append(el('text',{x:X(1.12*q.x),y:Y(1.12*q.y),class:'boundary-label','text-anchor':'middle','dominant-baseline':'central'},String(i)))}renderAnimationLabels(labels)
}
function tweenFrames(from,to,duration,labels=[]){if(duration<=0){renderAnimationFrame(to,labels);return Promise.resolve()}return new Promise(resolve=>{let started=null;function tick(timestamp){if(started===null)started=timestamp;const t=Math.min(1,(timestamp-started)/duration);drawMorph(from,to,t,labels);if(t<1)requestAnimationFrame(tick);else{renderAnimationFrame(to,labels);resolve()}}requestAnimationFrame(tick)})}
async function playMutation(payload,options={}){const duration=Number(animationSpeed.value),persistentLabels=options.labels===undefined?(labelsVisible?(currentState.faces||[]).map(f=>({label:[...f.label],x:f.x,y:f.y,font_size:f.font_size})):[]):options.labels,transform=options.transform||((frame)=>frame);mutationPlaying=true;animationCaption.hidden=false;try{let from=options.startFrame||stateAnimationFrame(currentState);for(const rawFrame of payload.frames||[]){const frame=transform(rawFrame);animationCaption.textContent=frame.caption;status.textContent=frame.caption;await tweenFrames(from,frame,duration,persistentLabels);from=frame}const finalFrame=stateAnimationFrame(payload.state,options.finalCaption||'Recompute layout and face labels');animationCaption.textContent=finalFrame.caption;status.textContent=finalFrame.caption;await tweenFrames(from,finalFrame,duration,persistentLabels);render(payload.state);status.textContent=options.completeMessage||'Square move complete.'}finally{mutationPlaying=false;animationCaption.hidden=true}}
async function transitionState(nextState,caption){const duration=Number(animationSpeed.value);if(!duration||!currentState||currentState.blank||nextState.blank||currentState.n!==nextState.n){render(nextState);return}const labels=labelsVisible?(currentState.faces||[]).map(f=>({label:[...f.label],x:f.x,y:f.y,font_size:f.font_size})):[];mutationPlaying=true;animationCaption.textContent=caption;animationCaption.hidden=false;status.textContent=caption;try{await tweenFrames(stateAnimationFrame(currentState),stateAnimationFrame(nextState,caption),duration,labels);render(nextState)}finally{animationCaption.hidden=true;mutationPlaying=false}}
const builder={active:false,n:6,vertices:[],edges:[],tool:'black',selected:null,history:[],drag:null,suppressClickUntil:0};
function builderSnapshot(){return JSON.stringify({vertices:builder.vertices,edges:builder.edges})}
function saveBuilder(){builder.history.push(builderSnapshot());if(builder.history.length>100)builder.history.shift()}
function setBuilderTool(tool){builder.tool=tool;builder.selected=null;for(const b of builderPanel.querySelectorAll('[data-builder-tool]'))b.classList.toggle('active',b.dataset.builderTool===tool);renderBuilder()}
function builderPoint(event){const point=svg.createSVGPoint();point.x=event.clientX;point.y=event.clientY;const matrix=svg.getScreenCTM();if(matrix){const local=point.matrixTransform(matrix.inverse());return [local.x,local.y]}const r=svg.getBoundingClientRect(),scale=Math.min(r.width/800,r.height/800),left=r.left+(r.width-800*scale)/2,top=r.top+(r.height-800*scale)/2;return [(event.clientX-left)/scale,(event.clientY-top)/scale]}
function boundaryBuilderPosition(i){const theta=-Math.PI/2-2*Math.PI*(i-1)/builder.n;return {id:i,x:Math.cos(theta),y:Math.sin(theta)}}
function builderVertex(id){if(id<=builder.n)return boundaryBuilderPosition(id);const v=builder.vertices[id-builder.n-1];return v&&{id,x:v.x,y:v.y,color:v.color}}
function builderAnimationFrame(){const vertices=[];for(let i=1;i<=builder.n;i++)vertices.push({...boundaryBuilderPosition(i),color:'boundary'});builder.vertices.forEach((v,j)=>vertices.push({id:builder.n+j+1,x:v.x,y:v.y,color:v.color}));return {caption:'Edited graph',n:builder.n,vertices,edges:builder.edges.map(e=>[e[0],e[1]])}}
function selectBuilderVertex(id,event){event.stopPropagation();if(performance.now()<builder.suppressClickUntil)return;if(builder.tool==='delete'){if(id<=builder.n){status.textContent='Boundary vertices cannot be deleted.';return}saveBuilder();builder.vertices.splice(id-builder.n-1,1);builder.edges=builder.edges.filter(e=>!e.includes(id)).map(e=>e.map(v=>v>id?v-1:v));builder.selected=null;renderBuilder();return}if(builder.tool!=='edge')return;if(builder.selected===null){builder.selected=id;renderBuilder();return}if(builder.selected===id){builder.selected=null;renderBuilder();return}const edge=[Math.min(builder.selected,id),Math.max(builder.selected,id)];if(!builder.edges.some(e=>e[0]===edge[0]&&e[1]===edge[1])){saveBuilder();builder.edges.push(edge)}builder.selected=null;renderBuilder()}
function renderBuilder(){svg.replaceChildren();svg.append(el('circle',{cx:400,cy:400,r:300,class:'disk'}));const layer=el('g');svg.append(layer);for(const [j,e] of builder.edges.entries()){const u=builderVertex(e[0]),v=builderVertex(e[1]);if(!u||!v)continue;const line=el('line',{x1:X(u.x),y1:Y(u.y),x2:X(v.x),y2:Y(v.y),class:'builder-edge'});line.addEventListener('click',event=>{if(performance.now()<builder.suppressClickUntil||builder.tool!=='delete')return;event.stopPropagation();saveBuilder();builder.edges.splice(j,1);renderBuilder()});layer.append(line)}for(let i=1;i<=builder.n;i++){const v=boundaryBuilderPosition(i),dot=el('circle',{cx:X(v.x),cy:Y(v.y),r:6,class:'builder-boundary'+(builder.selected===i?' builder-selected':'')});dot.addEventListener('click',e=>selectBuilderVertex(i,e));layer.append(dot);layer.append(el('text',{x:X(1.12*v.x),y:Y(1.12*v.y),'text-anchor':'middle','dominant-baseline':'central',class:'builder-number'},String(i)))}builder.vertices.forEach((v,j)=>{const id=builder.n+j+1,dot=el('circle',{cx:X(v.x),cy:Y(v.y),r:10,fill:v.color==='black'?'black':'white',stroke:'#26343a',class:'builder-vertex'+(builder.selected===id?' builder-selected':'')});dot.addEventListener('pointerdown',e=>beginBuilderDrag(id,e));dot.addEventListener('click',e=>selectBuilderVertex(id,e));layer.append(dot)});if(builder.vertices.length===0)layer.append(el('text',{x:400,y:400,class:'builder-hint'},'Click inside the disk to place a '+builder.tool+' vertex'));status.textContent='Manual drawing mode · '+builder.vertices.length+' internal vertices · '+builder.edges.length+' edges · drag a vertex to reposition it'}
function beginBuilderDrag(id,event){if(event.button!==0)return;const [px,py]=builderPoint(event);builder.drag={id,pointerId:event.pointerId,startX:px,startY:py,moved:false,snapshot:builderSnapshot()}}
svg.addEventListener('pointermove',event=>{const drag=builder.drag;if(!builder.active||!drag||event.pointerId!==drag.pointerId)return;const [px,py]=builderPoint(event);if(!drag.moved&&Math.hypot(px-drag.startX,py-drag.startY)<4)return;if(!drag.moved){drag.moved=true;svg.setPointerCapture(event.pointerId);builder.history.push(drag.snapshot);if(builder.history.length>100)builder.history.shift()}const v=builder.vertices[drag.id-builder.n-1];if(!v)return;let x=(px-400)/300,y=(400-py)/300;const radius=Math.hypot(x,y),limit=.94;if(radius>limit){x*=limit/radius;y*=limit/radius}v.x=x;v.y=y;renderBuilder()});
function endBuilderDrag(event){if(!builder.drag)return;if(builder.drag.moved)builder.suppressClickUntil=performance.now()+150;builder.drag=null;if(svg.hasPointerCapture(event.pointerId))svg.releasePointerCapture(event.pointerId)}
svg.addEventListener('pointerup',endBuilderDrag);svg.addEventListener('pointercancel',endBuilderDrag);
function startBuilder(){const n=Number(builderN.value);if(!Number.isInteger(n)||n<1||n>30){status.textContent='Choose between 1 and 30 boundary vertices.';return}builder.active=true;builder.n=n;builder.vertices=[];builder.edges=[];builder.selected=null;builder.history=[];builderResult.textContent='';builderPanel.hidden=false;setBuilderTool('black')}
function editCurrentGraph(){if(!currentState||currentState.blank||mutationPlaying)return;builder.active=true;builder.n=currentState.n;builderN.value=String(builder.n);builder.vertices=currentState.vertices.filter(v=>v.id>builder.n).sort((a,b)=>a.id-b.id).map(v=>({x:v.x,y:v.y,color:v.color}));builder.edges=currentState.edges.map(e=>[e[0],e[1]]);builder.selected=null;builder.history=[];builder.drag=null;builderResult.textContent='Editing permutation: ['+currentState.permutation.join(', ')+']';builderPanel.hidden=false;selectedStrands.clear();setBuilderTool('black')}
function leaveBuilder(){builder.active=false;builderPanel.hidden=true;if(currentState)render(currentState);else load()}
async function finishBuilder(){if(!builder.active||mutationPlaying)return;status.textContent='Preparing the reduction movie and decorated trip permutation…';try{const startFrame=builderAnimationFrame(),backendVertices=builder.vertices.map(v=>({x:v.x,y:v.y,color:v.color}));const r=await fetch('/custom-graph',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({n:builder.n,vertices:backendVertices,edges:builder.edges})});if(!r.ok)throw Error(await r.text());const payload=await r.json(),state=payload.state;builderResult.textContent='Permutation: ['+state.permutation.join(', ')+']';builder.active=false;builderPanel.hidden=true;selectedStrands.clear();await playMutation(payload,{startFrame,labels:[],finalCaption:'Arrange the reduced representative and compute its face labels',completeMessage:'Reduction complete. Decorated trip permutation: ['+state.permutation.join(', ')+'].'})}catch(e){status.textContent=e.message}}
svg.addEventListener('click',event=>{if(performance.now()<builder.suppressClickUntil||!builder.active||!['black','white'].includes(builder.tool)||event.target!==svg&&event.target.getAttribute('class')!=='disk')return;const [px,py]=builderPoint(event),x=(px-400)/300,y=(400-py)/300;if(x*x+y*y>=.94*.94)return;saveBuilder();builder.vertices.push({x,y,color:builder.tool});renderBuilder()});
manualButton.addEventListener('click',startBuilder);editGraphButton.addEventListener('click',editCurrentGraph);builderN.addEventListener('input',()=>{if(builder.active)startBuilder()});for(const b of builderPanel.querySelectorAll('[data-builder-tool]'))b.addEventListener('click',()=>setBuilderTool(b.dataset.builderTool));document.getElementById('builder-undo').addEventListener('click',()=>{if(!builder.history.length)return;const old=JSON.parse(builder.history.pop());builder.vertices=old.vertices;builder.edges=old.edges;builder.selected=null;renderBuilder()});document.getElementById('builder-clear').addEventListener('click',()=>{saveBuilder();builder.vertices=[];builder.edges=[];builder.selected=null;renderBuilder()});document.getElementById('builder-cancel').addEventListener('click',leaveBuilder);document.getElementById('builder-finish').addEventListener('click',finishBuilder);
async function load(){const r=await fetch('/state');if(!r.ok)throw Error(await r.text());render(await r.json())}
function labelText(label){return '{'+label.join(',')+'}'}
function defaultVariable(i){return i<=26?String.fromCharCode(96+i):'x_'+i}
function graphSignature(s){return JSON.stringify([s.permutation,s.edges,s.faces.map(f=>f.label)])}
function collectWeights(){const map=weightMode.value==='face'?faceWeights:edgeWeights;for(const input of weightList.querySelectorAll('input'))map.set(Number(input.dataset.id),input.value)}
function renderWeightEditor(s){const signature=graphSignature(s);if(signature!==weightSignature){faceWeights.clear();edgeWeights.clear();weightSignature=signature;measurementData=null;sourceInput.value=(s.default_sources||[]).join(',')}weightList.replaceChildren();const faceMode=weightMode.value==='face',items=faceMode?s.faces:s.edges;items.forEach((item,j)=>{const id=faceMode?item.id:j+1,map=faceMode?faceWeights:edgeWeights;const row=document.createElement('label');row.className='weight-row';const description=document.createElement('span');description.textContent=faceMode?'face '+labelText(item.label)+(id===1?' (reference)':''):'e'+id+' ('+item[0]+'—'+item[1]+')';const input=document.createElement('input');input.dataset.id=id;input.disabled=faceMode&&id===1;input.value=faceMode&&id===1?'dependent':(map.get(id)||'');input.placeholder='1';input.addEventListener('change',()=>{const value=input.value.trim();value&&value!=='1'?map.set(id,value):map.delete(id);measurementData=null;if(currentState)render(currentState)});row.append(description,input);weightList.append(row)});renderMeasurementResult()}
function assignAllEdgeVariables(){if(!currentState||currentState.blank)return;edgeWeights.clear();currentState.edges.forEach((_,j)=>edgeWeights.set(j+1,'t_'+(j+1)));weightMode.value='edge';measurementData=null;render(currentState);status.textContent='Assigned t_1 through t_'+currentState.edges.length+' to all '+currentState.edges.length+' edges.'}
function renderMeasurementResult(){if(!measurementData){measurementResult.textContent='Assign variables, then press Compute.';copyJulia.disabled=copyM2.disabled=true;return}const pl=resultMode.value==='pluckers';measurementResult.innerHTML='\\['+(pl?measurementData.plucker_latex:measurementData.matrix_latex)+'\\]';copyJulia.disabled=copyM2.disabled=false;if(window.MathJax?.typesetPromise)MathJax.typesetPromise([measurementResult])}
async function computeMeasurement(){if(!currentState||currentState.blank)return;collectWeights();const sources=sourceInput.value.split(/[ ,]+/).filter(Boolean).map(Number);if(sources.length!==currentState.k||sources.some(x=>!Number.isInteger(x))){status.textContent='Enter exactly '+currentState.k+' distinct integer sources.';return}status.textContent='Computing boundary measurement…';const map=weightMode.value==='face'?faceWeights:edgeWeights;try{const r=await fetch('/measurement',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({mode:weightMode.value,sources,weights:Object.fromEntries(map)})});if(!r.ok)throw Error(await r.text());measurementData=await r.json();renderMeasurementResult();status.textContent='Boundary measurement computed from sources '+measurementData.sources.join(', ')+'.'}catch(e){status.textContent=e.message}}
function renderDual(s,layer){const by=new Map(s.faces.map(f=>[f.id,f]));for(const ids of (s.dual_black_faces||[])){const faces=ids.map(id=>by.get(id)).filter(Boolean);if(faces.length>=3)layer.append(el('polygon',{points:faces.map(f=>X(f.dual_x)+','+Y(f.dual_y)).join(' '),class:'dual-black-face'}))}for(const [a,b] of s.dual_edges){const u=by.get(a),v=by.get(b);layer.append(el('line',{x1:X(u.dual_x),y1:Y(u.dual_y),x2:X(v.dual_x),y2:Y(v.dual_y),class:'dual-edge'}))}}
function render(s){builder.active=false;builderPanel.hidden=true;currentState=s;permInput.value=s.blank?'':'['+s.permutation.join(', ')+']';svg.replaceChildren();svg.append(el('circle',{cx:400,cy:400,r:300,class:'disk'}));const primal=el('g',{display:primalVisible?'inline':'none'}),dual=el('g',{display:dualVisible?'inline':'none'}),labels=el('g',{display:(primalVisible||dualVisible)?'inline':'none'});svg.append(primal);svg.append(dual);svg.append(labels);
 for(const f of s.faces){if(!f.movable)continue;primal.append(el('polygon',{points:f.polygon.map(p=>X(p[0])+','+Y(p[1])).join(' '),class:'square-face'}))}
 for(const f of s.faces){const hit=el('polygon',{points:f.polygon.map(p=>X(p[0])+','+Y(p[1])).join(' '),class:'face-hit'});hit.addEventListener('contextmenu',e=>openParameterMenu('face',f.id,e));if(f.movable)hit.addEventListener('dblclick',()=>move(f.label));primal.append(hit)}
 const by=new Map(s.vertices.map(v=>[v.id,v]));for(const [j,[a,b]] of s.edges.entries()){const u=by.get(a),v=by.get(b),line=el('line',{x1:X(u.x),y1:Y(u.y),x2:X(v.x),y2:Y(v.y),class:'edge'});line.addEventListener('contextmenu',e=>openParameterMenu('edge',j+1,e));primal.append(line);const hit=el('line',{x1:X(u.x),y1:Y(u.y),x2:X(v.x),y2:Y(v.y),class:'edge-hit'});hit.addEventListener('contextmenu',e=>openParameterMenu('edge',j+1,e));primal.append(hit)}
 const defs=el('defs');primal.append(defs);for(const strand of s.strands){if(!selectedStrands.has(strand.source))continue;const color=palette[(strand.source-1)%palette.length],marker=el('marker',{id:'arrow-'+strand.source,viewBox:'0 0 10 10',refX:8,refY:5,markerWidth:6,markerHeight:6,orient:'auto-start-reverse'});marker.append(el('path',{d:'M 0 0 L 10 5 L 0 10 z',fill:color}));defs.append(marker);const d=strand.points.map((p,j)=>(j?'L':'M')+X(p[0])+' '+Y(p[1])).join(' ');primal.append(el('path',{d,class:'strand',stroke:color,'marker-end':'url(#arrow-'+strand.source+')'}))}
 for(const v of s.vertices){if(v.id<=s.n)continue;primal.append(el('circle',{cx:X(v.x),cy:Y(v.y),r:8,fill:v.color==='black'?'black':'white',class:'internal'}))}
 for(let i=1;i<=s.n;i++){const v=by.get(i),toggle=()=>{selectedStrands.has(i)?selectedStrands.delete(i):selectedStrands.add(i);render(s)};const dot=el('circle',{cx:X(v.x),cy:Y(v.y),r:5,class:'boundary'});dot.addEventListener('dblclick',toggle);dot.append(el('title',{},'Double-click to toggle strand '+i));primal.append(dot);const number=el('text',{x:X(1.12*v.x),y:Y(1.12*v.y),class:'boundary-label','text-anchor':'middle','dominant-baseline':'central'},String(i));number.addEventListener('dblclick',toggle);primal.append(number)}
 renderDual(s,dual);
 if(labelsVisible)for(const f of s.faces){const labelX=dualVisible?f.dual_x:f.x,labelY=dualVisible?f.dual_y:f.y,t=el('text',{x:X(labelX),y:Y(labelY),'font-size':f.font_size,class:'face-label'+(f.movable?' movable':'')},labelText(f.label));t.addEventListener('contextmenu',e=>openParameterMenu('face',f.id,e));if(f.movable){t.setAttribute('title','Double-click to apply square move');t.addEventListener('dblclick',()=>move(f.label))}labels.append(t)}
 for(const f of s.faces){const value=faceWeights.get(f.id);if(value)labels.append(el('text',{x:X(f.x),y:Y(f.y)+(labelsVisible?16:0),class:'parameter-label'},value))}
 s.edges.forEach(([a,b],j)=>{const value=edgeWeights.get(j+1);if(!value)return;const u=by.get(a),v=by.get(b),dx=X(v.x)-X(u.x),dy=Y(v.y)-Y(u.y),length=Math.hypot(dx,dy)||1,offset=12;labels.append(el('text',{x:(X(u.x)+X(v.x))/2-offset*dy/length,y:(Y(u.y)+Y(v.y))/2+offset*dx/length,class:'parameter-label'},value))});
 faceList.replaceChildren();for(const f of s.faces){const item=document.createElement('div');item.className='face-item'+(f.movable?' square-item':'');item.textContent=labelText(f.label);item.title=f.movable?'Square face':'';faceList.append(item)}
 renderWeightEditor(s);
 primalButton.textContent=primalVisible?'Hide plabic graph':'Show plabic graph';dualButton.textContent=dualVisible?'Hide dual graph':'Show dual graph';
 backButton.disabled=!s.can_undo;forwardButton.disabled=!s.can_redo;editGraphButton.disabled=!!s.blank||s.n===0;allStrandsButton.disabled=s.n===0;assignEdgeVariablesButton.disabled=!!s.blank||s.edges.length===0;allStrandsButton.textContent=s.n>0&&selectedStrands.size===s.n?'Undraw all strands':'Draw all strands';
 status.textContent=s.blank?'Enter a decorated permutation and draw its graph.':'';}
async function move(label){if(mutationPlaying)return;status.textContent='Preparing square-move animation…';try{const r=await fetch('/move/'+label.join(','),{method:'POST'});if(!r.ok)throw Error(await r.text());await playMutation(await r.json())}catch(e){mutationPlaying=false;animationCaption.hidden=true;status.textContent=e.message}}
async function drawPermutation(){if(mutationPlaying)return;status.textContent='Replacing graph…';try{const r=await fetch('/permutation',{method:'POST',body:permInput.value});if(!r.ok)throw Error(await r.text());selectedStrands.clear();await transitionState(await r.json(),'Construct the new permutation graph');status.textContent='Graph replaced.'}catch(e){status.textContent=e.message}}
async function goBack(){if(mutationPlaying)return;status.textContent='Restoring previous graph…';try{const r=await fetch('/undo',{method:'POST'});if(!r.ok)throw Error(await r.text());selectedStrands.clear();await transitionState(await r.json(),'Restore the previous graph');status.textContent='Previous graph restored.'}catch(e){status.textContent=e.message}}
async function goForward(){if(mutationPlaying)return;status.textContent='Restoring next graph…';try{const r=await fetch('/redo',{method:'POST'});if(!r.ok)throw Error(await r.text());selectedStrands.clear();await transitionState(await r.json(),'Restore the next graph');status.textContent='Next graph restored.'}catch(e){status.textContent=e.message}}
document.getElementById('draw').addEventListener('click',drawPermutation);permInput.addEventListener('keydown',e=>{if(e.key==='Enter')drawPermutation()});
backButton.addEventListener('click',goBack);
forwardButton.addEventListener('click',goForward);
primalButton.addEventListener('click',()=>{primalVisible=!primalVisible;if(currentState)render(currentState)});dualButton.addEventListener('click',()=>{dualVisible=!dualVisible;if(currentState)render(currentState)});
allStrandsButton.addEventListener('click',()=>{if(!currentState)return;if(selectedStrands.size===currentState.n)selectedStrands.clear();else{selectedStrands.clear();for(let i=1;i<=currentState.n;i++)selectedStrands.add(i)}render(currentState)});
document.getElementById('copy-labels').addEventListener('click',async()=>{if(!currentState)return;const text=currentState.faces.map(f=>labelText(f.label)).join('\n');try{await navigator.clipboard.writeText(text);status.textContent='Face labels copied.'}catch(e){status.textContent='Select the labels in the list and copy them manually.'}});
weightMode.addEventListener('change',()=>{if(currentState)render(currentState)});assignEdgeVariablesButton.addEventListener('click',assignAllEdgeVariables);resultMode.addEventListener('change',renderMeasurementResult);document.getElementById('compute-measurement').addEventListener('click',computeMeasurement);
copyJulia.addEventListener('click',async()=>{if(measurementData)await navigator.clipboard.writeText(resultMode.value==='pluckers'?measurementData.plucker_julia:measurementData.matrix_julia)});copyM2.addEventListener('click',async()=>{if(measurementData)await navigator.clipboard.writeText(resultMode.value==='pluckers'?measurementData.plucker_m2:measurementData.matrix_m2)});
load().catch(e=>status.textContent=e.message);
</script></body></html>"""

function _parse_interactive_permutation(raw)
    cleaned=strip(replace(String(raw),r"[\[\]\(\)]"=>" "))
    parts=filter(!isempty,split(cleaned,r"[,\s]+"))
    isempty(parts) && throw(ArgumentError("enter a decorated permutation"))
    return parse.(Int,parts)
end

function _parse_interactive_drawing(raw)
    input=JSON.parse(String(raw))
    n=Int(get(input,"n",0))
    vertices=get(input,"vertices",Any[])
    colors=Symbol[Symbol(lowercase(string(get(v,"color","")))) for v in vertices]
    positions=[(Float64(v["x"]),Float64(v["y"])) for v in vertices]
    edges=[(Int(e[1]),Int(e[2])) for e in get(input,"edges",Any[])]
    return plabic_graph_from_drawing(n,colors,edges,positions;reduce=true)
end

function _parse_interactive_drawing_animation(raw)
    input=JSON.parse(String(raw))
    n=Int(get(input,"n",0))
    vertices=get(input,"vertices",Any[])
    colors=Symbol[Symbol(lowercase(string(get(v,"color","")))) for v in vertices]
    positions=[(Float64(v["x"]),Float64(v["y"])) for v in vertices]
    edges=[(Int(e[1]),Int(e[2])) for e in get(input,"edges",Any[])]
    return _drawing_reduction_animation(n,colors,edges,positions)
end

"""
    interactive_plabic_graph(G; port=8765, open_browser=true, iterations=2500, restarts=12)

Start a local interactive face-labelled drawing.  Valid internal square faces
have a light-blue background; double-clicking their red label applies the
square move and redraws the graph.
Double-clicking boundary vertex `i` toggles the medial strand originating at
`i`; selected strands persist across square moves.  The page also provides a
permutation input, an all-strands toggle, and a copyable face-label sidebar in
which square labels are red.  The default face-label convention is
`(:source, :left)`. Back/Forward navigate graph history. The boundary-
measurement panel accepts Postnikov face variables or individual edge weights
and displays the matrix or Plücker coordinates with Julia and Macaulay2 copy
formats. Keep the returned session alive and call `close(session)` when finished.
"""
function interactive_plabic_graph(G::PlabicGraph;port=8765,open_browser=true,
                                  iterations=2500,restarts=12,
                                  initial_blank=false)
    current=Ref(G)
    has_graph=Ref(!initial_blank)
    history=PlabicGraph[]
    future=PlabicGraph[]
    guard=ReentrantLock()
    function handler(request)
        target=String(request.target)
        if request.method=="GET" && target=="/"
            return HTTP.Response(200,["Content-Type"=>"text/html; charset=utf-8"],
                                 _INTERACTIVE_PLABIC_HTML)
        elseif request.method=="GET" && target=="/state"
            body=lock(guard) do
                state=has_graph[] ?
                    _interactive_state_json(current[];iterations=iterations,restarts=restarts) :
                    _interactive_blank_state_json()
                _interactive_history_json(state,!isempty(history),!isempty(future))
            end
            return HTTP.Response(200,["Content-Type"=>"application/json"],body)
        elseif request.method=="POST" && target=="/permutation"
            try
                body=lock(guard) do
                    permutation=_parse_interactive_permutation(String(request.body))
                    replacement=plabic_graph(permutation)
                    has_graph[] && push!(history,current[])
                    empty!(future)
                    current[]=replacement
                    has_graph[]=true
                    _interactive_history_json(
                        _interactive_state_json(current[];iterations=iterations,restarts=restarts),
                        !isempty(history),false)
                end
                return HTTP.Response(200,["Content-Type"=>"application/json"],body)
            catch err
                return HTTP.Response(400,sprint(showerror,err))
            end
        elseif request.method=="POST" && target=="/custom-graph"
            try
                body=lock(guard) do
                    replacement,stages=_parse_interactive_drawing_animation(request.body)
                    has_graph[] && push!(history,current[])
                    empty!(future)
                    current[]=replacement
                    has_graph[]=true
                    state=_interactive_history_json(
                        _interactive_state_json(current[];iterations=iterations,restarts=restarts),
                        !isempty(history),false)
                    _interactive_animation_json(stages,state)
                end
                return HTTP.Response(200,["Content-Type"=>"application/json"],body)
            catch err
                return HTTP.Response(400,sprint(showerror,err))
            end
        elseif request.method=="POST" && startswith(target,"/move/")
            has_graph[] || return HTTP.Response(400,"enter and draw a permutation first")
            raw=target[length("/move/")+1:end]
            label=try
                isempty(raw) ? Int[] : parse.(Int,split(raw,','))
            catch
                return HTTP.Response(400,"invalid face label")
            end
            try
                body=lock(guard) do
                    embedded,cycle=_square_cycle_by_label(current[],label;
                                                           iterations=iterations,restarts=restarts)
                    moved,stages=_square_move_animation_stages(embedded,cycle)
                    push!(history,current[])
                    empty!(future)
                    current[]=moved
                    state=_interactive_history_json(
                        _interactive_state_json(current[];iterations=iterations,restarts=restarts),
                        true,false)
                    _interactive_animation_json(stages,state)
                end
                return HTTP.Response(200,["Content-Type"=>"application/json"],body)
            catch err
                return HTTP.Response(400,sprint(showerror,err))
            end
        elseif request.method=="POST" && target=="/measurement"
            has_graph[] || return HTTP.Response(400,"enter and draw a permutation first")
            try
                body=lock(guard) do
                    _measurement_payload(current[],request.body;
                                         iterations=iterations,restarts=restarts)
                end
                return HTTP.Response(200,["Content-Type"=>"application/json"],body)
            catch err
                return HTTP.Response(400,sprint(showerror,err))
            end
        elseif request.method=="POST" && target=="/undo"
            try
                body=lock(guard) do
                    isempty(history) && throw(ArgumentError("there is no previous graph"))
                    push!(future,current[])
                    current[]=pop!(history)
                    has_graph[]=true
                    _interactive_history_json(
                        _interactive_state_json(current[];iterations=iterations,restarts=restarts),
                        !isempty(history),true)
                end
                return HTTP.Response(200,["Content-Type"=>"application/json"],body)
            catch err
                return HTTP.Response(400,sprint(showerror,err))
            end
        elseif request.method=="POST" && target=="/redo"
            try
                body=lock(guard) do
                    isempty(future) && throw(ArgumentError("there is no next graph"))
                    push!(history,current[])
                    current[]=pop!(future)
                    has_graph[]=true
                    _interactive_history_json(
                        _interactive_state_json(current[];iterations=iterations,restarts=restarts),
                        true,!isempty(future))
                end
                return HTTP.Response(200,["Content-Type"=>"application/json"],body)
            catch err
                return HTTP.Response(400,sprint(showerror,err))
            end
        end
        return HTTP.Response(404,"not found")
    end
    server=nothing
    chosen_port=Int(port)
    for candidate in chosen_port:(chosen_port+20)
        try
            server=HTTP.serve!(handler,"127.0.0.1",candidate;verbose=false)
            chosen_port=candidate
            break
        catch err
            candidate==chosen_port+20 && rethrow(err)
        end
    end
    url="http://127.0.0.1:$chosen_port/"
    session=InteractivePlabicGraph(current,server,url)
    if open_browser
        try
            Sys.isapple() ? run(`open $url`) :
            Sys.iswindows() ? run(`cmd /c start $url`) : run(`xdg-open $url`)
        catch err
            @warn "Could not open a browser automatically; open the session URL manually" url exception=err
        end
    end
    return session
end

interactive_plabic_graph(p::AbstractVector{<:Integer};kwargs...) =
    interactive_plabic_graph(plabic_graph(p);kwargs...)

"""
    interactive_session(; kwargs...)

Launch a blank interactive plabic-graph workspace. Enter a decorated
permutation in the browser and press **Draw / replace graph** to begin. Keywords
such as `port`, `open_browser`, `iterations`, and `restarts` are forwarded to
`interactive_plabic_graph`.
"""
interactive_session(;kwargs...) =
    interactive_plabic_graph(plabic_graph([1]);initial_blank=true,kwargs...)

function Base.close(session::InteractivePlabicGraph)
    close(session.server)
    return nothing
end

"""Internal, dependency-free materialization of the `k`-subsets of `items`."""
function _subsets(items, k::Integer)
    k < 0 && return Vector{Vector{eltype(items)}}()
    k == 0 && return [eltype(items)[]]
    k > length(items) && return Vector{Vector{eltype(items)}}()
    result = Vector{Vector{eltype(items)}}()
    chosen = Vector{eltype(items)}()
    function visit(first)
        length(chosen) == k && (push!(result, copy(chosen)); return)
        needed = k - length(chosen)
        for j in first:(length(items) - needed + 1)
            push!(chosen, items[j]); visit(j + 1); pop!(chosen)
        end
    end
    visit(1)
    return result
end

function _permutation_cycles(p::AbstractVector{<:Integer})
    q = abs.(p)
    seen = falses(length(q))
    result = Vector{Vector{Int}}()
    for start in eachindex(q)
        seen[start] && continue
        cycle = Int[]
        i = start
        while !seen[i]
            push!(cycle, i)
            seen[i] = true
            i = q[i]
        end
        push!(result, cycle)
    end
    return result
end

"""
    plucker_coordinates(A)

Return the maximal minors of a `k × n` matrix, indexed by increasing tuples of
column indices. These are the Plücker coordinates in the orientation of `A`.
"""
function plucker_coordinates(A::AbstractMatrix)
    k, n = size(A)
    k <= n || throw(ArgumentError("expected a k × n matrix with k ≤ n"))
    return Dict(Tuple(I) => det(A[:, I]) for I in _subsets(Vector(1:n), k))
end

function _minor_nonnegative(x, atol)
    x isa Real || throw(ArgumentError("positivity is defined here only for real matrices"))
    return x isa AbstractFloat ? x >= -atol : x >= zero(x)
end


function _minor_positive(x, atol)
    x isa Real || throw(ArgumentError("positivity is defined here only for real matrices"))
    return x isa AbstractFloat ? x > atol : x > zero(x)
end

"""Return whether every maximal minor of `A` is nonnegative."""
function is_totally_nonnegative(A::AbstractMatrix; atol=1e-10)
    all(x -> _minor_nonnegative(x, atol), values(plucker_coordinates(A)))
end

"""Return whether every maximal minor of `A` is strictly positive."""
function is_totally_positive(A::AbstractMatrix; atol=1e-10)
    all(x -> _minor_positive(x, atol), values(plucker_coordinates(A)))
end


function shiftedOrder(i::Int64,n::Int64)
    return vcat(i:n, 1:i-1)
end

# returns the cyclic interval a ..... b or a .... n 1 .... b-1 depending on a>b or b<a.
function cyclicInterval(n::Int64, i::Int64, j::Int64)
    if i<=j 
        return Vector(i:j)
    end

    if i>j 
        return reduce(vcat, (Vector(i:n), Vector(1:j-1)))
    end
end

# Compares two elements in the Gale order of i
function iCompare(i::Int64, n::Int64, x::Int64, y::Int64)
    
    L = shiftedOrder(i,n);
    
    posx = (findall(t->t==x, L))[1];
    posy = (findall(t->t==y, L))[1];

    return posx >= posy;
end


# Compares two Lists in the induced Gale order of i
function iCompareLists(i::Int64,n::Int64, X::Vector{Int64}, Y::Vector{Int64})
    length(X) == length(Y) || throw(ArgumentError("Gale order requires subsets of equal size"))
    isempty(X) && return true
    
    L = shiftedOrder(i,n);
    
    Xind = [[(findall(t->t==x, L))[1], x] for x in X];
    Yind = [[(findall(t->t==y, L))[1], y] for y in Y];

    Xind = sort(Xind)
    Yind = sort(Yind)

    Xsor = [x[2] for x in Xind]
    Ysor = [y[2] for y in Yind]

    comparisons = unique([iCompare(i,n, Xsor[j], Ysor[j]) for j in 1:size(X)[1] ])

    if(size(comparisons)[1] > 1) 
        return false
    end
    
    if(comparisons[1] == false)
        return false
    end

    return true;
end


# Gets the minimal element in the ith Gale order of sequences.
function getMinList(i::Int64,n::Int64, B::Set{Vector{Int64}})

    X = collect(B)[1];

    for Y in B

        if iCompareLists(i,n, X, Y)
            X = Y;
        end

    end

    return X

end



# Gets the min Necklace given a positroid
function getNeckLace(n::Int64, B::Set{Vector{Int64}})

    return [getMinList(i,n,B) for i in 1:n]
    
end


function minGrassmannNecklace(n::Int64, B::Set{Vector{Int64}})
    return getNeckLace(n, B)
end


# Gets the minimal element in the ith Gale order of sequences.
function getMaxList(i::Int64,n::Int64, B::Set{Vector{Int64}})

    X = collect(B)[1];

    for Y in B

        if iCompareLists(i,n, Y, X)
            X = Y;
        end

    end

    return X

end


# Gets the max Necklace given a positroid
function getMaxNeckLace(n::Int64, B::Set{Vector{Int64}})

    return [getMaxList(i,n,B) for i in 1:n]
    
end


function maxGrassmannNecklace(n::Int64, B::Set{Vector{Int64}})
    return getMaxNeckLace(n, B)
end


function fromNecklaceToPositroid(k::Int64,n::Int64, N::Vector{Vector{Int64}})
    
    L = _subsets(Vector(1:n), k);
    Intervals = [];
    Positroid = [];

    interval = [];

    for j in 1:n
        I = N[j];
        interval = [];

        for J in L
            
            if iCompareLists(j,n, J, I) 
                append!(interval, [J]);
            end
        end

        append!(Intervals, [interval]);

    end

    S = Set{Vector{Int64}}(Intervals[1]);
    
    for j in 2:n
        S = intersect(S, Set{Vector{Int64}}(Intervals[j]));
    end

    return S;
end



function fromNecklaceToDecoratedPermutation(N::Vector{Vector{Int64}})
    
    n = size(N)[1];
    p = Vector{Int64}(1:n);

    for i in 1:n
        
        j = i+1;

        if j == n+1
            j = 1;
        end

        I1 = Set(N[i])
        I2 = Set(N[j])

        if I1 != I2
            p_i = collect(setdiff(I2,I1))[1]
            p[p_i] = i
        end
        
        if I1 == I2

            if i in I1
                p[i] = -i
            end
            if ! (i in I1)
                p[i] = i
            end

        end
    end
    return p
end


# Computes the Grassmann Necklace associated to the decorated permutation
function fromDecoratedPermToNecklace(p::Vector{Int64})

    n = size(p)[1];

    N = Vector{Vector{Int64}}([]);
    
    for k in 1:n

        I = Vector{Int64}([]);

        for i in 1:n
            
            if abs(p[i]) == i
                if abs(p[i]) == - p[i]
                    append!(I, [i])
                end
            else 
                if iCompare(k,n,   abs(p[i]) , i)
                    append!(I, [i])
                end
            end
        end

        append!(N, [I]);
    end


    return N;
    
end


function fromDecoratedPermToMinNecklace(p::Vector{Int64})
    return fromDecoratedPermToNecklace(p)
end


function minGrassmannNecklace(p::Vector{Int64})
    return fromDecoratedPermToMinNecklace(p)
end


function fromDecoratedPermToMaxNecklace(p::Vector{Int64})

    _validate_decorated_permutation(p)

    n = length(p)
    k = countExceedences(p)

    if k == 0
        return [Int64[] for i in 1:n]
    end

    if k == n
        return [Vector(1:n) for i in 1:n]
    end

    return getMaxNeckLace(n, fromDecoratedPermToPositroid(k, n, p))
end


function maxGrassmannNecklace(p::Vector{Int64})
    return fromDecoratedPermToMaxNecklace(p)
end


function _affine_value_from_decorated_perm(p::Vector{Int64}, i::Int64)

    n = length(p)
    pi = abs(p[i])

    if pi == i
        return p[i] < 0 ? i + n : i
    end

    return pi > i ? pi : pi + n
end


function reverseGrassmannNecklace(p::Vector{Int64})

    _validate_decorated_permutation(p)

    n = length(p)
    shifts = [_affine_value_from_decorated_perm(p, i) - i for i in 1:n]
    N = Vector{Vector{Int64}}([])

    for q in 1:n
        I = Vector{Int64}([])

        for r in 1:n
            representative = r <= q ? r : r - n

            if representative <= q && q < representative + shifts[r]
                push!(I, modN(representative, n))
            end
        end

        push!(N, sort(I))
    end

    return N
end


# This is either wrong or the next function is wrong
function fromDecoratedPermToPositroid(k::Int64,n::Int64, p::Vector{Int64})
    _validate_decorated_permutation(p)
    length(p) == n || throw(ArgumentError("n=$n does not match the permutation length $(length(p))"))
    countExceedences(p) == k || throw(ArgumentError("k=$k does not match the permutation's $(countExceedences(p)) excedances"))
    N = fromDecoratedPermToNecklace(p);

    Pos = fromNecklaceToPositroid(k,n, N); 
    return Pos;
end


function fromPositroidToPermutation(n::Int64, M::Set{Vector{Int64}})
    
    N = getNeckLace(n,M);

    return fromNecklaceToDecoratedPermutation(N);
end


# Function to generate vertices of a regular n-gon
function generate_ngon_vertices(n::Int, r::Float64 = 1.0)
    angles = range(0, 2π, length=n+1)[1:end-1]  # n angles evenly spaced
    return [(r*cos(-θ-π/2), r*sin(-θ-π/2)) for θ in angles]
end


# Function to draw the n-gon and chords
function draw_chords(p::Vector{Int64})

    _load_plotting()

    r = 1.0;

    n = size(p)[1]

    vertices = generate_ngon_vertices(n)

    x_coords = [v[1] for v in vertices]
    y_coords = [v[2] for v in vertices]


    # Plot the disk (circle)
    θ = range(0, 2π, length=100)
    Plots.plot(r*cos.(θ), r*sin.(θ), seriestype = :path, lw=2, label="Disk", legend=false, ratio=1, fill=:blue, fillalpha=0.1)

    # Plot the n-gon
    #plot!(x_coords, y_coords, seriestype = :shape, lw=3, label="n-gon", legend=false, ratio=1)

    # Draw all possible chords
    for i in 1:n
        
        
        if abs(p[i])!= i
            Plots.GR.setarrowsize(2)
            j = abs(p[i]);
            Plots.plot!([x_coords[i], x_coords[j]], [y_coords[i], y_coords[j]], lw=4, linecolor=:black, label=false, arrow=(:closed, 2.0))
        end


        if  p[i] < 0
            Plots.GR.setarrowsize(2)
            x_text = x_coords[i] * 1.2   # Slightly offset the text position
            y_text = y_coords[i] * 1.2
            Plots.annotate!(x_text, y_text, Plots.text("-", "red", 26))
        end

    end

    # Number the vertices clockwise from 1 to n
    for i in 1:n
       
        label = string(i)
        x_text = x_coords[i] * 1.05  # Slightly offset the text position
        y_text = y_coords[i] * 1.05
        Plots.annotate!(x_text, y_text, label)
       
       
        #if i%2 == 1
        #    x_text = x_coords[i] * 1.2  # Slightly offset the text position
        #    y_text = y_coords[i] * 1.2
        #    annotate!(x_text, y_text, text("+", "blue", 20))
        #end
        
        #if i%2 == 0
        #    x_text = x_coords[i] * 1.2  # Slightly offset the text position
        #    y_text = y_coords[i] * 1.2
        #    annotate!(x_text, y_text, text("-", "red", 20))
        #end
    end


    # Display the plot
    return Plots.plot!(showaxis = false, ticks = false, ylimits=(-1.5,1.5),  xlimits=(-1.5,1.5))
end


















function signOfList(I::Vector{Int64})

    indexedI = [[I[i], i] for i in 1:size(I)[1]]

    sortedI = sort(indexedI);

    p = [x[2] for x in sortedI];

    inversions = count(i > j for (a, i) in enumerate(I) for j in I[(a + 1):end])
    return iseven(inversions) ? 1 : -1
end

function is_orthogonal(k::Int64, n::Int64, M::Set{Vector{Int64}}, verbose=false)

    W = _subsets(Vector(1:n), k-1);

    N = binomial(n,k-1);

    for i in 1:N
        
        for j in i:N

            I = W[i];
            J = W[j];

            IJcomp = sort(collect(setdiff( Set(Vector(1:n)) ,union(Set(I), Set(J)) )));

            L1 = filter( l-> signOfList(reduce(vcat, (I, [l])) ) * signOfList(reduce(vcat, (J, [l])) ) * ( sort(reduce(vcat, (I, [l])) ) in M) * ( sort(reduce(vcat, (J, [l])) ) in M) * (-1)^l ==  1 , IJcomp) 
            L2 = filter( l-> signOfList(reduce(vcat, (I, [l])) ) * signOfList(reduce(vcat, (J, [l])) ) * ( sort(reduce(vcat, (I, [l])) ) in M) * ( sort(reduce(vcat, (J, [l])) ) in M) * (-1)^l == -1 , IJcomp) 

            if  ( isempty(L1) &&  ! isempty(L2) ) || (!isempty(L1) && isempty(L2)) 
                
                if verbose
                    println("Problem at: I= ", I," and J= ",J);
                    println("+ side indices : ",L1, " and - side indices:", L2 )
                end

                return false;

            end

        end

    end


    return true
end


function is_orthogonal(p::Vector{Int64}, verbose=false)
    k = countExceedences(p);
    n = size(p)[1];
    M = fromDecoratedPermToPositroid(k,n,p)
    return is_orthogonal(k,n,M);
end

function alignmentNumber(p::Vector{Int64})
    
    n = size(p)[1]

    counter = 0


    for i in 1:n
        for j in 1:n

            if i != j


                pi = abs(p[i])
                pj = abs(p[j])

                

                i_pj_interval = cyclicInterval(n, i,pj)
                pj_i_interval = cyclicInterval(n, pj,i)
                
                v = (pi in i_pj_interval) &&  (j in pj_i_interval)

                if (pi == i)
                    v = v && (p[i] < 0)
                end

                if (pj == j)
                    v = v && (p[j] > 0)
                end
                
                if v 
                    
                    counter = counter + 1;
                end

            end

        end
    end
    
    return counter
    
end


function dimensionOfPositroid(k::Int64, n::Int64, M::Set{Vector{Int64}})

    p = fromPositroidToPermutation(n, M);

    return k*(n-k) - alignmentNumber(p);
    
end

function dimensionOfPermutation(k::Int64, n::Int64, p::Vector{Int64})
    
    return k*(n-k) - alignmentNumber(p);
    
end


function numVertices(n)

    if n in 0:3
        return 0
    end

    s = 0

    for l in 1:trunc(Int64, n/2)

        if n%2 ==0

            s = s+ (l-1)^2 +   trunc(Int64, ((n-2*l)/2)^2);

        end

        if n%2 == 1

            s = s+ (l-1)^2 +   trunc(Int64, ((n+1-2*l)/2) * ((n-1-2*l)/2) );

        end

    end
    
    return s + numVertices(n-1)
    
end

function countExceedences(p::Vector{Int64})
    return sum([( i < abs(p[i])) + (i == - p[i]) for i in 1:size(p)[1]]);
end



function permutationFromCycles(n::Int64, C::Vector{Vector{Int64}})

    p = [0 for i in 1:n]

    for i in 1:n
        
        for I in C
            if i in I
                j = findall(t->t==i, I)[1];

                if j < size(I)[1]
                    p[i] = I[j+1]
                else
                    p[i] = I[1]
                end
            end
        end


        if p[i] == 0
            p[i] = i
        end
    end

    return p;


end




# Define a function to check the orientation of three points
function orientation(p, q, r)
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if val > 0
        return 1  # Clockwise
    elseif val < 0
        return -1 # Counter-clockwise
    else
        return 0  # Collinear
    end
end

# Check if a point is on the segment defined by two other points
function on_segment(p, q, r)
    return min(p[1], r[1]) <= q[1] <= max(p[1], r[1]) && min(p[2], r[2]) <= q[2] <= max(p[2], r[2])
end

# Check if two line segments (p1, q1) and (p2, q2) intersect
function do_intersect(p1, q1, p2, q2)
    # Find the orientations of the ordered triplets (p1, q1, p2), (p1, q1, q2), (p2, q2, p1), (p2, q2, q1)
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case: segments intersect if orientations are different
    if o1 != o2 && o3 != o4
        return true
    end

    # Special cases
    if o1 == 0 && on_segment(p1, p2, q1)
        return true
    end
    if o2 == 0 && on_segment(p1, q2, q1)
        return true
    end
    if o3 == 0 && on_segment(p2, p1, q2)
        return true
    end
    if o4 == 0 && on_segment(p2, q1, q2)
        return true
    end

    return false
end

# Function to check if two polygons intersect
function polygons_intersect(poly1, poly2)
    n1 = length(poly1)
    n2 = length(poly2)

    # Check every edge of poly1 against every edge of poly2
    for i in 1:n1
        p1 = poly1[i]
        q1 = poly1[(i % n1) + 1]
        for j in 1:n2
            p2 = poly2[j]
            q2 = poly2[(j % n2) + 1]
            if do_intersect(p1, q1, p2, q2)
                return true  # An intersection was found
            end
        end
    end
    return false  # No intersection found
end

# Function to group only the data based on polygon connected components
function group_data_by_polygon_components(polygons_with_data)
    n = length(polygons_with_data)
    visited = falses(n)  # Array to keep track of visited polygons
    components = []      # Array to store connected components of data

    # Define a DFS function to explore connected polygons and collect data
    function dfs(polygon_idx, component_data)
        visited[polygon_idx] = true
        # Collect the data associated with the polygon
        push!(component_data, polygons_with_data[polygon_idx][2])

        # Check all other polygons for intersection
        for i in 1:n
            if !visited[i] && polygons_intersect(polygons_with_data[polygon_idx][1], polygons_with_data[i][1])
                dfs(i, component_data)
            end
        end
    end

    # Loop over all polygons and perform DFS for unvisited polygons
    for i in 1:n
        if !visited[i]
            component_data = []
            dfs(i, component_data)
            push!(components, component_data)
        end
    end

    return components
end


function orthogonal_dimension(k::Int64, n::Int64, M::Set{Vector{Int64}})

 

    p = fromPositroidToPermutation(n,M);
    k = countExceedences(p);
    vertices = generate_ngon_vertices(n);


    cyc = filter(c -> length(c) > 1, _permutation_cycles(p));
    polygons_with_cycles = [ [[vertices[i] for i in c], c] for c in cyc];




    N = 0;

    for i in 1:size(cyc)[1]
    
        ki = countExceedences(permutationFromCycles(n,[polygons_with_cycles[i][2]]));

        N = N + binomial(ki+1,2);

        for j in i+1:size(cyc)[1]
            if polygons_intersect(polygons_with_cycles[i][1], polygons_with_cycles[j][1])
                
                kj = countExceedences(permutationFromCycles(n,[polygons_with_cycles[j][2]]))
                N = N + ki*kj
            end
        end
    end


    return dimensionOfPositroid(k,n,M) - N;

end





function orthogonal_dimensionOfPermutation(p::Vector{Int64})
    n = size(p)[1];
    k = countExceedences(p);

    M = fromDecoratedPermToPositroid(k,n,p)
    
    return orthogonal_dimension(k,n,M);
    
end







function inversionsOfPerm(p::Vector{Int64})

    n = size(p)[1];

    Invs = []

    for i in 1:n-1
        for j in i+1:n

            if p[i]> p[j]
                append!(Invs, [(i,j)]);
            end
        end
    end
    
    return Invs
end


function pminversionsOfPerm(p::Vector{Int64})

    n = size(p)[1];

    Invs = []

    for i in 1:n-1
        for j in i+1:n

            if p[i]> p[j] && (-1)^(i+j) == -1
                append!(Invs, [(i,j)]);
            end
        end
    end
    
    return size(Invs)[1]
end



# The order on the positroids is not containement!
function is_codim1_face(p::Vector{Int64}, q::Vector{Int64})
    
    n = size(p)[1]
    k = countExceedences(p);
    
    if countExceedences(q) !=k
        return "The two permutations should have the same k";
    end

    d = orthogonal_dimensionOfPermutation(p)

    if d != orthogonal_dimensionOfPermutation(q) -1
        return false;
    end

    Mp = fromDecoratedPermToPositroid(k,n,p)
    Mq = fromDecoratedPermToPositroid(k,n,q)

    return isempty(setdiff(Mp, Mq));
end


function is_child(p::Vector{Int64}, q::Vector{Int64})
    
    n = size(p)[1]
    k = countExceedences(p);
    
    if countExceedences(q) !=k
        return "The two permutations should have the same k";
    end

    Mp = fromDecoratedPermToPositroid(k,n,p)
    Mq = fromDecoratedPermToPositroid(k,n,q)

    return isempty(setdiff(Mp, Mq));
end

"""
    immediate_children(p)

Return the codimension-one boundary cells covered by the positroid cell `p`.
Each result has the same rank, dimension `dimension(p)-1`, and a set of bases
contained in that of `p`. The current exhaustive algorithm is intended for
moderate ground-set sizes.
"""
function immediate_children(input::AbstractVector{<:Integer})
    p = Int.(input)
    _validate_decorated_permutation(p)
    n = length(p)
    k = decorated_excedances(p)
    dimension = dimensionOfPermutation(k, n, p)
    dimension == 0 && return Vector{Vector{Int}}()
    parent_bases = fromDecoratedPermToPositroid(k, n, p)
    children = Vector{Vector{Int}}()
    for candidate in decorated_permutations(n, k)
        dimensionOfPermutation(k, n, candidate) == dimension - 1 || continue
        child_bases = fromDecoratedPermToPositroid(k, n, candidate)
        isempty(setdiff(child_bases, parent_bases)) || continue
        push!(children, candidate)
    end
    sort!(children)
    return children
end

function _validate_decorated_permutation(p::Vector{Int64})

    n = length(p)

    if any(x -> x == 0 || abs(x) > n, p)
        error("A decorated permutation of [1,$n] can only contain signed values in 1:$n.")
    end

    if sort(abs.(p)) != Vector(1:n)
        error("The absolute values of a decorated permutation must be a permutation of 1:$n.")
    end

    for i in 1:n
        if abs(p[i]) != i && p[i] < 0
            error("Only fixed points may carry a negative decoration; p[$i] = $(p[i]).")
        end
    end

    return true
end


function _le_shape_data(I::Vector{Int64}, n::Int64)

    I = sort(I)
    k = length(I)
    Ic = sort(collect(setdiff(Set(Vector(1:n)), Set(I))), rev=true)

    row_lengths = [n - k + r - I[r] for r in 1:k]

    if any(l -> l < 0 || l > n-k, row_lengths)
        error("The first necklace term $I does not determine a Young diagram inside a $k x $(n-k) rectangle.")
    end

    if !issorted(row_lengths, rev=true)
        error("The row lengths $row_lengths are not a partition shape.")
    end

    return I, Ic, row_lengths
end


function le_diagram_matrix(p::Vector{Int64})

    _validate_decorated_permutation(p)

    n = length(p)
    N = fromDecoratedPermToNecklace(p)
    I, Ic, row_lengths = _le_shape_data(N[1], n)

    T = [zeros(Int, row_lengths[r]) for r in 1:length(row_lengths)]

    row_index = Dict(I[r] => r for r in 1:length(I))
    column_index = Dict(Ic[c] => c for c in 1:length(Ic))

    for j in 2:n

        J = sort(N[j])

        A = sort(collect(setdiff(Set(I), Set(J))), rev=true)
        B = sort(collect(setdiff(Set(J), Set(I))))

        if length(A) != length(B)
            error("Invalid Grassmann necklace step at index $j: removed $A but added $B.")
        end

        for l in 1:length(A)
            a, b = A[l], B[l]

            rrow = row_index[a]
            ccolumn = column_index[b]

            if ccolumn > length(T[rrow])
                error("The box labeled by row $a and column $b lies outside the Young diagram.")
            end

            T[rrow][ccolumn] = 1
        end
    end

    return T
end


function le_diagram(p::Vector{Int64}; as_tableau::Bool=true)

    T = le_diagram_matrix(p)

    # `as_tableau` is retained for source compatibility. Earlier prototypes
    # returned an Oscar tableau here; the package now returns its dependency-free
    # row-vector representation in both modes.
    return T
end



function is_Omega_symplectic(k::Int64, n::Int64, M::Set{Vector{Int64}})

    m = n//2;

    V = _subsets(Vector(1:n), k-2);

    for S in V 

        A_odd  = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 1   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );
        A_even = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 0   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );

        if  ( isempty(A_odd) &&  ! isempty(A_even) ) || (!isempty(A_odd) && isempty(A_even)) 
          return false
        end

    end

    return true
end




function is_Omega_symplecticTwo(k::Int64, n::Int64, M::Set{Vector{Int64}})

    m = n//2;

    V = _subsets(Vector(1:n), k-2);

    for S in V 

        A_odd  = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 1   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );
        A_even = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 0   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );


        if (size(A_odd)[1] in [0,1]) && (size(A_even)[1] in [0,1])

            if  (isempty(A_odd) &&  ! isempty(A_even)) || (!isempty(A_odd) && isempty(A_even)) 
                return false
            end

        end

    end

    return true
end




function symplectic_Omega_dimension(k::Int64, n::Int64, M::Set{Vector{Int64}})

    if k > 2 return "Not implemented for k>2 yet!" end 

    m = n//2;

    A_odd  = filter(t-> (t%2==1 && [t,m+t] in M) == true, Vector(1:m));
    A_even = filter(t-> (t%2==0 && [t,m+t] in M) == true, Vector(1:m));

    if !isempty(A_odd) && !isempty(A_even)
        return dimensionOfPositroid(k,n,M)-1;
    end

    return dimensionOfPositroid(k,n,M);
end




function modN(i,N)

    if i > N
        return i - N
    elseif i < 1
        return i + N
    else
        return i
    end
    
end




function is_E_symplectic(k::Int64, n::Int64, M::Set{Vector{Int64}})

    m = n//2;

    V = _subsets(Vector(1:n), k-2);

    for S in V 

        A_odd  = filter(t-> (t%2==1 && sort(unique(vcat(S,[t,2*m-t+1]))) in M) == true, Vector(1:m));
        A_even = filter(t-> (t%2==0 && sort(unique(vcat(S,[t,2*m-t+1]))) in M) == true, Vector(1:m));

        if  ( isempty(A_odd) &&  ! isempty(A_even) ) || (!isempty(A_odd) && isempty(A_even)) 
            return false
        end

    end

    return true
end



function symplectic_E_dimension(k::Int64, n::Int64, M::Set{Vector{Int64}})

    if k > 2 return "Not implemented for k>2 yet!" end 

    m = n//2;

    A_odd  = filter(t-> (t%2==1 && [t,2*m-t+1] in M) == true, Vector(1:m));
    A_even = filter(t-> (t%2==0 && [t,2*m-t+1] in M) == true, Vector(1:m));

    if !isempty(A_odd) && !isempty(A_even)
        return dimensionOfPositroid(k,n,M)-1;
    end

    return dimensionOfPositroid(k,n,M);
end



function up(k::Int64, n::Int64,  M::Set{Vector{Int64}})

    if n!=2k+1 
        error("n != 2k + 1");
    end

    newM = Set{Vector{Int64}}([]);

    for I in M 
        J = reduce(vcat, [I, [n+1]])
        push!(newM, J);

        Jc = filter( t-> ! (t in J) == true,  1:n+1);
        push!(newM, Jc);
    end

    return newM;
end


function randomMatching(n::Int64)

    if n%2 != 0
        error("n is not even!!");
    end

    V = Vector(1:n);

    C = Vector{Vector{Int64}}([]);
    
    for i in 1:n//2
    
        a = rand(V);
        V = filter(j-> !(j==a)==true, V);
        b = rand(V);
        V = filter(j-> !(j==b)==true, V);
    
        C = append!(C, [[a,b]]);
    end
    

    return permutationFromCycles(n, C);

end












function twistPerm(p)

    m = size(p)[1];
    n = trunc(Int64, m/2)
    q = Vector{Int64}([]);

    pp = abs.(p);

    for i in 1:m

        e = 1;

        if pp[i] == -p[i]
            e = -1;
        end

        push!(q, e* modN( pp[i]+n,m));

    end


    return q

end





function untwistPerm(p)

    m = size(p)[1];
    n = trunc(Int64, m/2)
    q = Vector{Int64}([]);

    pp = abs.(p);


    for i in 1:m
        
        e = 1

        if p[i] < 0
            e = -1;
        end

        push!(q,  e * modN( pp[i] - n, m));
    end

    return q

end

function twistNecklace(I:: Vector{Vector{Int64}})

    m = size(I)[1];
    n = trunc(Int64,m/2);
    
    Res = Vector{Vector{Int64}}([]);

    J = [I[ modN(i+n,m) ] for i in 1:m];

    for i in 1:m

        S = J[i];

        rotS = [modN(x+n,m) for x in S]

        push!(Res, filter( x-> !(x in rotS) , 1:m))

    end

    return Res

end


function _validate_twist_necklace(N::Vector{Vector{Int64}}, k::Int64, n::Int64)

    if length(N) != n
        error("Expected a Grassmann necklace with $n terms, but got $(length(N)).")
    end

    for i in 1:n
        I = sort(N[i])

        if length(I) != k
            error("Necklace term $i has length $(length(I)); expected $k.")
        end

        if length(unique(I)) != k || any(j -> j < 1 || j > n, I)
            error("Necklace term $i is not a $k-subset of 1:$n: $(N[i]).")
        end
    end

    return [sort(I) for I in N]
end


function _is_zero_for_twist(x, atol)
    if x isa AbstractFloat
        return abs(x) <= atol
    end

    return iszero(x)
end


function _is_nonzero_for_twist(x, atol)
    if x isa AbstractFloat
        return abs(x) > atol
    end

    return !iszero(x)
end


function _check_twist_basis(A::AbstractMatrix, I::Vector{Int64}, i::Int64, atol)

    d = det(A[:, I])

    if !_is_nonzero_for_twist(d, atol)
        error("The columns indexed by necklace term I_$i = $I are not a basis for the given matrix.")
    end

    return nothing
end


function _check_twist_loop_column(A::AbstractMatrix, i::Int64, atol)

    if any(x -> !_is_zero_for_twist(x, atol), A[:, i])
        error("Necklace term I_$i does not contain $i, so the cell expects column $i to be a loop/zero column.")
    end

    return nothing
end


function _solve_twist_column(A::AbstractMatrix, I::Vector{Int64}, position::Int64)

    rhs = zeros(eltype(A), size(A, 1))
    rhs[position] = one(eltype(A))

    return transpose(A[:, I]) \ rhs
end


# Right twist of a point A in the positroid cell with Grassmann necklace N.
# This is symbolic-friendly: pass a matrix over SymPy.Sym, rationals, or another
# exact field/ring and the returned entries are computed by exact linear solves.
function right_twist(A::AbstractMatrix, N::Vector{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)
    N = _validate_twist_necklace(N, k, n)

    twisted_columns = Vector{Any}(undef, n)

    for i in 1:n
        I = N[i]
        position = findfirst(==(i), I)

        if isnothing(position)
            if check_cell
                _check_twist_loop_column(A, i, atol)
            end

            twisted_columns[i] = nothing
            continue
        end

        if check_cell
            _check_twist_basis(A, I, i, atol)
        end

        twisted_columns[i] = _solve_twist_column(A, I, position)
    end

    nonloop_column = findfirst(c -> !isnothing(c), twisted_columns)

    if isnothing(nonloop_column)
        return zeros(eltype(A), k, n)
    end

    zero_column = zero(twisted_columns[nonloop_column])

    for i in 1:n
        if isnothing(twisted_columns[i])
            twisted_columns[i] = zero_column
        end
    end

    return hcat(twisted_columns...)
end


function right_twist(A::AbstractMatrix, p::Vector{Int64}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if length(p) != n
        error("The decorated permutation has length $(length(p)), but the matrix has $n columns.")
    end

    if countExceedences(p) != k
        error("The decorated permutation has k = $(countExceedences(p)), but the matrix has $k rows.")
    end

    return right_twist(A, fromDecoratedPermToNecklace(p); check_cell=check_cell, atol=atol)
end


function right_twist(A::AbstractMatrix, M::Set{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if !isempty(M) && length(collect(M)[1]) != k
        error("The positroid consists of $(length(collect(M)[1]))-subsets, but the matrix has $k rows.")
    end

    return right_twist(A, getNeckLace(n, M); check_cell=check_cell, atol=atol)
end


function left_twist(A::AbstractMatrix, N::Vector{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    return right_twist(A, N; check_cell=check_cell, atol=atol)
end


function left_twist(A::AbstractMatrix, p::Vector{Int64}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if length(p) != n
        error("The decorated permutation has length $(length(p)), but the matrix has $n columns.")
    end

    if countExceedences(p) != k
        error("The decorated permutation has k = $(countExceedences(p)), but the matrix has $k rows.")
    end

    return left_twist(A, reverseGrassmannNecklace(p); check_cell=check_cell, atol=atol)
end


function left_twist(A::AbstractMatrix, M::Set{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if !isempty(M) && length(collect(M)[1]) != k
        error("The positroid consists of $(length(collect(M)[1]))-subsets, but the matrix has $k rows.")
    end

    return left_twist(A, reverseGrassmannNecklace(fromPositroidToPermutation(n, M)); check_cell=check_cell, atol=atol)
end


function twist_map(A::AbstractMatrix, cell; side::Symbol=:right, check_cell::Bool=true, atol=1e-10)

    if side == :right
        return right_twist(A, cell; check_cell=check_cell, atol=atol)
    end

    if side == :left
        return left_twist(A, cell; check_cell=check_cell, atol=atol)
    end

    error("Unknown twist side $side. Use side=:right or side=:left.")
end





## Function that computes dimension of a decorated involution


function dimensionOfLGperm(p::Vector{Int64})

    q = abs.(twistPerm(p));

    m = size(p)[1];
    n = trunc(Int64, m/2);

    chords = Vector{Vector{Int64}}([]);
    
    visited = Vector{Int64}([]);

    for i in 1:m

        if !(i in visited)

            j = q[i];

            if j != i

                push!(chords, [i,j])
                append!(visited, [i,j]);   
            end

            append!(visited, [i,j]);   
        end

    end


    D = 0;

    for c in chords
        
        D = D + min( modN(c[1] - c[2],m), modN(c[2] - c[1],m)  );

    end

    if size(chords)[1] > 1   

            for a in 1:size(chords)[1]-1

                for b in a+1:size(chords)[1]

                    c1 = chords[a];
                    c2 = chords[b];

                    i,j = 1, modN(c1[2] - c1[1]+1,m)
                    k,l = modN(c2[1] - c1[1]+1,m) , modN(c2[2] - c1[1]+1,m)    
                    
                    if min(k,l) < j && j < max(k,l) 

                        D = D-1;

                    end

                end
            end
    end


    return trunc(Int64, n*(n+1)/2) - D

end










function iota(I::Vector{Int64})

    n = size(I)[1];
    J = [modN(i+n,2*n) for i in I];    
    return sort(filter( x-> !(x in J) , 1:2*n))

end





function eta(i::Int64,n::Int64)

    modN(i+1,2*n)
end


function eta(I::Vector{Int64},n::Int64)
    return sort([eta(i,n) for i in I])
end




function iota2(I::Vector{Int64})

    n = size(I)[1];
    J = eta(I,n)
    return sort(filter( x-> !(x in J) , 1:2*n))

end





function is_cyclo_isotropic(k::Int64, n::Int64, M::Set{Vector{Int64}})

    for I in M

        if !(iota2(I) in M)
            return false
        end
    end
    return true
end

"""
    decorated_permutations(n, k)

Return all signed-vector encodings of decorated permutations on `1:n` with
exactly `k` excedances. A negative fixed point contributes one excedance;
only fixed points may be negative.
"""
function decorated_permutations(n::Integer, k::Integer)
    n >= 0 || throw(ArgumentError("n must be nonnegative"))
    0 <= k <= n || return Vector{Vector{Int}}()

    results = Vector{Vector{Int}}()
    decorated = Vector{Int}(undef, n)
    used = falses(n)

    function generate(position::Int, excedances::Int)
        excedances > k && return
        excedances + (n - position + 1) < k && return
        if position > n
            excedances == k && push!(results, copy(decorated))
            return
        end
        for value in 1:n
            used[value] && continue
            used[value] = true
            if value == position
                decorated[position] = value
                generate(position + 1, excedances)
                decorated[position] = -value
                generate(position + 1, excedances + 1)
            else
                decorated[position] = value
                generate(position + 1, excedances + (position < value))
            end
            used[value] = false
        end
    end

    generate(1, 0)
    return results
end

"""Count excedances in the signed-vector encoding of a decorated permutation."""
function decorated_excedances(p::AbstractVector{<:Integer})
    count(i -> i < abs(p[i]) || p[i] == -i, eachindex(p))
end

include("web_server.jl")

end # module Positroids
