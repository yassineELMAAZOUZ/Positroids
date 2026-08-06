# Positroids.jl

Experimental Julia tools for positroids, decorated permutations, Grassmann
necklaces, total positivity, and isotropic variants. The combinatorial core is
dependency-free; Plots.jl is used only by `draw_chords`.

This repository consolidates the Julia prototypes in `Positroids_code.jl` and
`decorated_permutations.jl`. The historical Mathematica implementation
(`positroids.m`) and demonstration notebook are retained as reference sources;
see [the migration inventory](docs/src/mathematica-migration.md).

## Install and load locally

From Julia in this directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Positroids
```

## Quick start

Decorated permutations use signed one-line notation. Only fixed points may be
negative, and a negative fixed point counts as an excedance.

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
retained as the dependent reference face for face-weight coordinates. The
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
- **Experimental:** twists; orthogonal, symplectic, and cyclic-isotropy tests;
  specialized dimension formulas and face predicates.
- **Reference only:** Mathematica boundary-poset traversal, symbolic
  matrix/residue machinery, and scattering-amplitude code.
  These have not yet been ported and are not silently advertised as Julia APIs.

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

The unrelated Java hyperbolic-coloring README previously at the repository root
is preserved as `HYPERBOLIC_README.md`.
