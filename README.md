# Positroids.jl

Julia tools for positroids and plabic graphs.

## Install and load locally

In Julia:

using Pkg;

Pkg.add(url="https://github.com/yassineELMAAZOUZ/Positroids.git");

<<<<<<< HEAD
using Positroids;

interactive_session();
=======
Decorated permutations use signed one-line notation. Only fixed points may be
negative. The convention is uniform throughout the package: `p[i] = -i` is a
**coloop**, counts toward the rank, and is drawn as a **white lollipop**;
`p[i] = i` is a **loop** and is drawn as a **black lollipop**.

For the target/trip convention used by plabic graphs and face labels, use
`positroid_rank(p)`. For example,

```julia
positroid_rank([4, 3, 1, 5, 2])  # 2
```

`decorated_excedances` and the legacy necklace conversion routines retain the
historical source convention for compatibility.

```julia
p = [-1, 3, 4, 2]
k = decorated_excedances(p)                 # 3
N = minGrassmannNecklace(p)                 # Grassmann necklace
M = fromDecoratedPermToPositroid(k, 4, p)   # set of bases

fromNecklaceToDecoratedPermutation(N) == p  # true
fromPositroidToPermutation(4, M) == p       # true
dimensionOfPermutation(k, 4, p)
le_diagram_matrix(p)
```

Enumerate all rank-`k` decorated permutations on `[n]`:

```julia
cells = decorated_permutations(4, 2)
```

For a matrix representing a point of a positroid cell:

```julia
A = [1.0 0.0 1.0 1.0;
     0.0 1.0 1.0 2.0]
p = [3, 4, 1, 2]
T = right_twist(A, p; check_cell=false)

coords = plucker_coordinates(A)
is_totally_nonnegative(A)
```

Construct and draw a canonical bridge plabic graph:

```julia
p = [3, 4, 1, 2]
G = plabic_graph(p)
trip_permutation(G) == p
is_bipartite(G)  # true
is_reduced(G)    # true
draw_plabic_graph(G)
```

The drawing uses a constrained force-directed charge model by default. For
harder examples, increase the layout search:

```julia
draw_plabic_graph(G; layout_iterations=2500, layout_restarts=12)
```

Construct either of the two cell charts:

```julia
p = [3, 4, 1, 2]

B = bridge_parametrization(p)
parameter_names(B)
symbolic_matrix(B)
A = parametrization_matrix(B, [1.0, 2.0, 3.0, 4.0])

N = boundary_measurement_parametrization(p)
parameter_names(N)
plucker_coordinates(N)  # symbolic Plücker polynomials in w₁,…,w₄
C = parametrization_matrix(N, [1.0, 2.0, 3.0, 4.0])
plucker_coordinates(C)
```

Enumerate codimension-one boundary cells in the positroid poset:

```julia
children = immediate_children([3, 4, 1, 2])
```

The boundary helpers use the target/trip convention by default, matching the
interactive graph and its face labels. Thus every child of
`[4,3,1,5,2]` also has `positroid_rank(child) == 2`; callers do not need to
apply `invperm`. Pass `permutation_convention=:source` only for historical
source-form input.

Enumerate the entire proper boundary or count its cells by dimension:

```julia
boundary = boundary_cells([3, 4, 1, 2])
f = boundary_f_vector([3, 4, 1, 2])
# [6, 12, 10, 4], counting dimensions 0, 1, 2, 3

closed_f = boundary_f_vector([3, 4, 1, 2]; include_cell=true)
# [6, 12, 10, 4, 1]
```

The traversal deduplicates cells reached along different chains. Both functions
accept `max_cells` to place an explicit limit on large boundary computations.

Test the generic differential rank of a linear projection
`C -> rowspace(C*Z)` on facets or on the full positroid boundary:

```julia
p = [4, 3, 1, 5, 2]  # rank 2; no invperm conversion is needed
Z = [ 1  0  0  0;
      0  1  0  0;
      0  0  1  0;
      0  0  0  1;
     -1  1 -1  1]

projection_jacobian_report(p, Z)
facet_ranks = boundary_projection_jacobian_report(p, Z; strata=:facets)
all_boundary_ranks = boundary_projection_jacobian_report(p, Z; strata=:boundary)
poset = projection_boundary_poset(p, Z)
# poset.nodes, poset.covers, poset.levels, poset.f_vector
```

These projection functions use the target/anti-excedance permutation
convention by default. Pass `permutation_convention=:source` only when supplying
a permutation in the package's historical source-trip convention.

The calculation uses affine Plucker ratios and exact rational arithmetic.
`certificate == :full_rank` proves generic local injectivity on the chart;
`certificate == :rank_drop` proves a generic rank loss using a symbolic
projected-coordinate upper bound. Jacobian rank alone does not prove global
one-to-one behavior. The complete worked example is in
[`examples/projection_boundary_analysis.jl`](examples/projection_boundary_analysis.jl).

In an interactive session, press **Show facets** to place these codimension-one
permutations below the main graph. Each facet appears as copyable permutation
text beside a compact, static plabic-graph thumbnail. The list is part of the
main scrollable workspace and can be hidden with the same button. Facet
thumbnails are deliberately read-only: square moves, strands,
weights, and boundary measurements continue to apply only to the original
large graph. **Compute f-vector** traverses the entire proper boundary and shows
only its copyable f-vector and total cell count; it never constructs descendant
graphs.

Perform a square move on a displayed square face:

```julia
G = plabic_graph([3, 4, 1, 2])
squares = square_cycles(G)
H = square_move(G, squares[1])
draw_plabic_graph(H)
```

Or select the square by its face label:

```julia
H = square_move_by_label(G,[1,3])
H = square_move(G;face_label=[1,3])
```

Compare the graph before and after the move without replacing either panel:

```julia
comparison = compare_plabic_graphs(
    G,H;
    titles=("Before","After"),
    face_labels=true,
)
display(comparison)
```

Compute and display trip target labels on all disk faces:

```julia
face_labels(G)
draw_plabic_graph(G;face_labels=true)
```

The default `(:source,:left)` convention uses the literal geometric left side
of each oriented trip from `i` to `π(i)`; rotating boundary vertex `1` does not
reverse this convention.

Select a different face-label convention when desired:

```julia
face_labels(G;convention=(:target,:right))
draw_plabic_graph(G;face_labels=true,
                  face_label_convention=(:target,:right))
```

Overlay all trips, or only selected strands:

```julia
draw_plabic_graph(G;trips=:all)
draw_plabic_graph(G;trips=[1,3,5])
```

Open a live face-labelled drawing and perform square moves by double-clicking a
valid square-face label:

```julia
session = interactive_plabic_graph(G)
# later:
close(session)
```

Or launch a blank workspace and enter the permutation entirely in the browser:

```julia
session = interactive_session()
```

For a hosted, multi-user service use `serve_positroids_web`. Deployment and
Docker instructions are in [`DEPLOYMENT.md`](DEPLOYMENT.md).

The blank workspace also has **Draw graph manually**. Choose the number of
boundary vertices, place black and white internal vertices by clicking in the
disk, and drag any internal vertex to reposition it while its incident edges
follow live. Boundary vertex `1` starts at the bottom and the remaining labels
run clockwise. Then use **Add edges** to join pairs of vertices. **Reduce & use graph**
contracts same-color edges, suppresses bivalent internal vertices, checks
planarity, and computes the decorated trip permutation. If local reductions do
not suffice, it replaces the drawing by the canonical reduced graph for the
same trip permutation, then loads it into the usual interactive tools. The drawing editor includes undo,
clear, delete, and cancel controls.

After drawing a graph from a permutation—or after applying square moves—press
**Edit current graph** to transfer that exact graph and layout into the drawing
editor. You can move vertices, add or delete vertices and edges, then press
**Reduce & use graph** to compute its new decorated trip permutation and return
to all of the usual interactive tools.

The interactive view uses `(:source,:left)` labels, shades movable square faces
light blue while keeping every label red, and updates in place after each move. Double-click a
boundary vertex to toggle its originating medial strand. The current graph is
available as `session.graph[]`. Face labels are centered at polygon barycenters
and automatically shrink to fit their individual faces. Boundary-face labels
move to their adjacent boundary arcs only while the dual graph is displayed.

The interactive page also has a permutation input, a draw/undraw-all-strands
button, and a copyable face-label sidebar whose square-face entries are red.
Independent show/hide buttons overlay the plabic graph and its face-adjacency
dual on the same disk. The face labels themselves are the dual vertices (there
are no extra circles or duplicate labels), and dual faces corresponding to
black plabic vertices are shaded dark gray. Without the dual graph, boundary-
face labels stay at their face barycenters; they move to their boundary arcs
only when the dual graph is shown. Square labels can still be
double-clicked to perform moves. The **Back** button restores the graph before
the most recent square move or permutation replacement and can be used
repeatedly to walk backward through the session history; **Forward** walks
through drawings that were undone.

The **Graph transforms** controls rotate the whole drawing—including boundary
vertices and their indices—without changing the underlying graph. **Indices
−1** and **Indices +1** cyclically relabel the boundary at its fixed geometric
locations, recompute the decorated permutation and face data, and enter the
Back/Forward history. These index changes are immediate and are not animated.
**Swap black ↔ white** interchanges every colored internal
vertex; boundary vertices remain uncolored, and the reversed trip permutation
and fixed-point decorations are updated automatically. Rotation also remains
active in the hand-drawing editor, where pointer positions are converted back
to disk coordinates so vertices still land exactly where they are clicked.

Graph changes are animated continuously. During a square move, vertices travel
smoothly to their new positions, contracted edges visibly shrink to zero,
created edges grow from their source locations, colors blend, and removed
vertices collapse away. The motion covers trivalentization, color switching,
same-color contractions, bivalent removal, and final layout recomputation. The
original face labels remain visible above the graph throughout the motion and
are replaced by the recomputed labels at completion. The **Mutation speed**
selector ranges from Instant through Movie mode.

The boundary-measurement panel supports Postnikov face weights and arbitrary
edge weights. **Hide face labels** independently hides the combinatorial
`k`-subset labels. Right-click a face or edge to assign, rename, or clear a
parameter directly on the drawing; it appears in blue at that face or edge and
stays synchronized with the sidebar. Unassigned weights equal one. One face is
retained as the dependent reference face for face-weight coordinates. One
context-sensitive button assigns `s_1,s_2,\ldots` to all independent faces in
face mode or `t_1,t_2,\ldots` to all edges in edge mode. Only the selected
mode's parameters are drawn; the other assignments are hidden but retained.
Individual values can still be changed with the sidebar or right-click menu. The
resulting `k × n` matrix or its Plücker coordinates are typeset with
LaTeX below the graph in the same scrollable workspace and can be copied in
Julia or Macaulay2 syntax. A source-set field chooses the `k` pivot columns; admissible
choices are rendered as an identity submatrix. Symbolic products and fractions
are simplified before display. The right control panel contains the matrix /
Plücker selector, result text-size and result-height controls. Dragging the
divider along its left edge resizes the panel; dragging it fully right hides
the panel, and dragging left reopens it. Double-clicking the divider also
toggles between hidden and the default width.

Graph drawings use a deterministic harmonic (Tutte-style) layout by default,
which spreads internal vertices through the disk without repulsive charges.
The former layout remains available as
`layout_plabic_graph(G;method=:force)` or
`draw_plabic_graph(G;layout_method=:force)`.

For decorated fixed points, white lollipops are drawn with clockwise closed
strands and contribute their indices to every source-left face label; black
lollipops are counterclockwise and contribute to none.

## Scope and maturity

- **Core:** decorated permutations, cyclic/Gale order, necklace/positroid
  conversions, enumeration, cell dimensions, and Le-diagrams.
- **Experimental:** matrix twists and some legacy necklace helpers.

The current basis representation is `Set{Vector{Int64}}`; because vectors are
mutable keys, callers must not mutate a basis after inserting it into a set.
Future releases should replace this with an immutable typed representation.

See [the API and capability guide](docs/src/index.md) for conventions, function
groups, limitations, and migration status.

## Test

```julia
using Pkg
Pkg.test()
```
>>>>>>> f7b41bf (Fix permutation conventions and add interactive graph transforms)
