# Positroids.jl capability guide

## Mathematical conventions

A decorated permutation is stored in signed one-line notation. The vector entry
`p[i] = -i` denotes a negatively decorated fixed point; `p[i] = i` denotes a
positively decorated fixed point. Nonfixed entries must be positive. This is the
convention used by the original Julia prototype. It differs from the bounded
affine-permutation notation used by the Mathematica package, where values may be
shifted by `n`.

A Grassmann necklace is a `Vector{Vector{Int64}}` with one `k`-subset per cyclic
position. A positroid is currently a `Set{Vector{Int64}}` of sorted bases.

## Core combinatorics

| Purpose | Julia API | Status |
|---|---|---|
| Enumerate decorated permutations | `decorated_permutations(n, k)` | Tested core |
| Count excedances/rank | `decorated_excedances`, `countExceedences` | Tested core |
| Cyclic order and intervals | `shiftedOrder`, `cyclicInterval` | Core, legacy names |
| Gale comparisons | `iCompare`, `iCompareLists` | Core, legacy names |
| Minimum necklace | `minGrassmannNecklace(p)` or `(n, bases)` | Tested core |
| Maximum/reverse necklaces | `maxGrassmannNecklace`, `reverseGrassmannNecklace` | Experimental |
| Necklace to positroid bases | `fromNecklaceToPositroid(k, n, N)` | Tested core |
| Necklace to permutation | `fromNecklaceToDecoratedPermutation(N)` | Tested core |
| Permutation to positroid | `fromDecoratedPermToPositroid(k, n, p)` | Tested core |
| Positroid to permutation | `fromPositroidToPermutation(n, M)` | Tested core |
| Dimension | `dimensionOfPermutation`, `dimensionOfPositroid` | Core |
| Le-diagram | `le_diagram_matrix`, `le_diagram` | Experimental; row-vector output |

Conversion routines enumerate all `k`-subsets when reconstructing bases, so
their cost grows as `binomial(n,k)`. They are intended for moderate examples.

## Matrix and total-positivity tools

`plucker_coordinates(A)` computes all maximal minors of a `k × n` matrix and
indexes them by increasing column tuples. `is_totally_nonnegative(A)` tests that
all of these coordinates are nonnegative, while `is_totally_positive(A)` uses
strict positivity. Both accept an `atol` keyword for floating-point matrices and
use exact comparisons for exact real element types. The signs refer to the
orientation of the supplied representative matrix.

### Cell parametrizations

`bridge_parametrization(p)` returns a `BridgeParametrization` implementing the
Mathematica-style bridge chart. It applies the bridge recursion to the inverse
decorated permutation, which translates Mathematica's bounded-affine convention
to this package's excedance convention. There is one positive coordinate `αᵢ`
per bridge, and the parameter count equals the cell dimension.

`boundary_measurement_parametrization(p)` returns a separate
`BoundaryMeasurementParametrization`. This is the boundary measurement of the
normalized directed bridge network: `wᵢ` is the nontrivial bridge-edge weight,
while gauge-tree weights are fixed to one. Matrix entries are the corresponding
path sums.

```julia
p = [3, 4, 1, 2]
B = bridge_parametrization(p)
parameter_names(B)
symbolic_matrix(B)
A = parametrization_matrix(B, [1.0, 2.0, 3.0, 4.0])

N = boundary_measurement_parametrization(p)
symbolic_minors = plucker_coordinates(N)
C = N([1.0, 2.0, 3.0, 4.0])
plucker_coordinates(C)
```

Calling `plucker_coordinates` directly on either parametrization returns
symbolic polynomial strings. Like terms and determinant cancellations are
combined internally; calling it on an evaluated matrix continues to return
numeric or exact scalar minors.

For an integer matrix `Z`, `symbolic_product(P, Z)` (or simply `P * Z`)
multiplies the symbolic parametrization matrix by `Z` and returns simplified
polynomial strings. For numerical parameters, evaluate first and use ordinary
matrix multiplication: `P(weights) * Z`.

Positive coordinates give a totally nonnegative representative whose nonzero
maximal minors are exactly the bases of the requested positroid.

`right_twist`, `left_twist`, and `twist_map` accept a `k × n` matrix and either
a necklace, decorated permutation, or set of bases. With `check_cell=true`, the
code verifies required necklace bases and loop columns. Floating-point checks use
the `atol` keyword; exact element types use exact zero tests.

These routines are experimental. In particular, the right and left twists share
a linear-solve kernel and still require broader comparison against the
Mathematica implementation and published convention choices.

`draw_chords(p)` visualizes the chord diagram through Plots/GR.

## Isotropic variants

The package contains research-stage predicates and formulas:

- orthogonal membership and dimension: `is_orthogonal`,
  `orthogonal_dimension`, `orthogonal_dimensionOfPermutation`;
- symplectic membership: `is_Omega_symplectic`,
  `is_Omega_symplecticTwo`, `is_E_symplectic`;
- specialized symplectic dimensions: `symplectic_Omega_dimension`,
  `symplectic_E_dimension` (currently only implemented for `k ≤ 2`);
- twists/involutions and Lagrangian dimension: `twistPerm`, `untwistPerm`,
  `twistNecklace`, `dimensionOfLGperm`;
- cyclic isotropy: `iota`, `iota2`, `eta`, `is_cyclo_isotropic`.

These functions preserve the prototype algorithms and should be treated as
experimental until examples and theorem-level regression tests are added.

## Boundary helpers

`is_child` compares basis containment. `is_codim1_face` combines containment
with the prototype orthogonal-dimension calculation. Despite their names, these
do not yet reproduce the Mathematica package's complete `boundary` and
`inverseBoundary` machinery.

`immediate_children(p)` enumerates the ordinary positroid-poset covers below
`p`: same-rank cells of dimension one less whose bases are contained in those of
`p`. It currently checks all decorated permutations of the same `(k,n)`, so it
is reliable for moderate `n` but scales factorially.

## Plabic graphs

`plabic_graph(p)` constructs a deterministic bridge graph from a decorated
permutation. It ports the default recursive bridge construction used by the
Mathematica package. Positive and negative decorated fixed points become black
and white lollipops respectively. The result is a `PlabicGraph` containing the
bicolored vertices, edges, disk positions, bridge decomposition, and input trip
permutation certificate.

`draw_plabic_graph(G)` draws this graph with Plots.jl, and the convenience form
`draw_plabic_graph(p)` constructs and draws it in one call. Use
`trip_permutation(G)` to retrieve the graph's decorated permutation.

Raw bridge reconstruction can create same-color edges and internal bivalent
vertices. Julia repeatedly suppresses each bivalent vertex and contracts every
same-color component until neither remains. Thus a chain
`black—white—black` loses its bivalent white vertex and the two resulting
adjacent black vertices are merged. `is_bipartite(G)` checks the resulting
coloring. `is_reduced(G)` verifies the bridge construction certificate: valid
degrees including the absence of valence two, univalent boundary vertices,
bipartiteness, and equality among bridge count, graph dimension, and
positroid-cell dimension. `plabic_graph(p)` calls
`validate_plabic_graph(G)` and throws rather than returning a graph if these
certificates fail.

`square_cycles(G)` lists induced alternating internal four-cycles eligible as
square-move candidates. Apply `square_move(G, cycle)` to a candidate that is a
face in the displayed embedding. The operation trivalentizes higher-valence
square vertices, switches their colors, reduces same-color and bivalent
configurations, and revalidates the result. `square_move(G)` is shorthand when
there is exactly one candidate. A rotation system is not yet stored, so candidate
enumeration cannot always distinguish a facial square from a separating
four-cycle; explicit selection keeps that choice visible.

Once faces are labeled, `square_move_by_label(G,label)` selects a square by its
target `k`-subset. The equivalent keyword form is
`square_move(G;face_label=label)`. The label must identify a unique internal
face whose boundary consists of four alternating vertices; boundary or
nonsquare faces are rejected.

```julia
H = square_move_by_label(G,[1,3])
# equivalently:
H = square_move(G;face_label=Set([1,3]))
```

### Trips, faces, and target labels

`plabic_embedding(G)` finds a crossing-free disk layout, adds boundary-circle
arcs, constructs the counterclockwise rotation system, and enumerates faces by
directed half-edge traversal. It independently traces trips by turning maximally
left at white vertices and maximally right at black vertices, and rejects the
embedding unless their endpoints reproduce the decorated permutation.

`graph_faces(G)`, `graph_trips(G)`, and `face_labels(G)` expose these data. For
each trip from `i` to `π(i)`, the default convention adds `i` to every face on
the left of the oriented trip. The trip edges form a barrier in the dual graph; a flood fill
from the locally left faces captures the entire left-hand region, not merely the
faces touching the trip. Here “left” and “right” are literal geometric sides of
the oriented trip and are unaffected by rotating boundary vertex `1`.

```julia
G = plabic_graph([3,4,1,2])
face_labels(G)
draw_plabic_graph(G;face_labels=true)
```

Choose any of the four source/target and left/right conventions with a tuple of
symbols. The same option is available in the drawing function:

```julia
face_labels(G;convention=(:target,:right))
draw_plabic_graph(G;face_labels=true,
                  face_label_convention=(:target,:right))
```

Source-left labels have size equal to the number of strict anti-exceedances
`π(i)<i` plus negatively decorated fixed points. Right-side labels have the
complementary size.

Decorated fixed points use the lollipop-loop convention. A white lollipop
(negative decoration) carries a clockwise closed strand based at its boundary
vertex, so its index belongs to every source-left face label. A black lollipop
(positive decoration) carries a counterclockwise loop, so its index belongs to
no source-left face label. Right-side labels make the complementary choice.
Thus `face_labels(plabic_graph([-1,-2,3,4]))` is `[[1,2]]`, even though the
graph has only one disk face.

Labels are red by default and placed at the area-weighted polygon barycenter of
each face using the shoelace formula. Each label is independently shrunk from
`face_label_size` until a conservative estimate of its rendered text box fits
inside the face. For a concave face whose barycenter lies outside its polygon,
the interior visual center is used instead. Disable fitting with
`adaptive_face_label_size=false`, or explicitly select the farthest-from-edges
placement with `face_label_position=:visual_center`. Customize appearance with
`face_label_color` and `face_label_size`. Interactive drawings use the same
barycenter-and-fit calculation and recompute sizes after every square move.

`compare_plabic_graphs(G,H;...)` preserves both drawings in independent
side-by-side panels, which is useful before and after a square move. Supply
`titles=("Before","After")` and any ordinary drawing keywords such as
`face_labels=true`.

The embedding routine rejects layouts with crossings; increasing its `restarts`
keyword gives difficult graphs more candidate layouts.

Set `trips=:all` (or `trips=true`) in `draw_plabic_graph` to overlay every
oriented medial strand. Pass `trips=[1,3,5]` to show selected strands. A strand
starts beside its boundary vertex, crosses every traversed graph edge at its
midpoint, and follows a curved arc around each internal vertex: clockwise at
white and counterclockwise at black. Appearance is controlled by
`strand_wiggle`, `strand_samples`, `strand_boundary_offset`, `strand_width`,
`strand_alpha`, and `strand_arrows`.

```julia
draw_plabic_graph(G;trips=:all)
draw_plabic_graph(G;trips=[1,3],face_labels=true)
```

### Interactive square moves

Start a live local-browser drawing whose face labels use the default
`(:source,:left)` convention. Valid internal square faces have a light-blue
background and all labels remain red; double-click a square's label to apply
the move, reduce the graph, recompute the embedding and labels, and update the
same drawing. Double-click boundary vertex
`i` to show or hide the medial strand originating at `i`. Multiple strands can
be shown at once, and selected strands remain visible after square moves.

The page includes a permutation input and **Draw / replace graph** button. It
replaces the current graph, face sidebar, and strand data with the graph of the
new decorated permutation without returning to the Julia prompt. **Draw all
strands** toggles between displaying and hiding every strand while individual
boundary-vertex toggles remain available. A sidebar lists every face label as
selectable text; **Copy all** copies the newline-separated list, and square-face
labels are red in the sidebar.

Press **Show dual graph** to overlay the planar dual on the same disk as the
plabic graph. Every disk face becomes a labelled dual vertex, with the face
label itself serving as the vertex rather than a separate circle, and two dual
vertices are joined by a dashed edge when the corresponding faces share a
plabic edge. **Show/Hide plabic graph** and **Show/Hide dual graph** control the
two layers independently, so either graph or both can be displayed. Square moves,
permutation replacement, strand controls, and the face sidebar continue to work
while the dual panel is visible; it is recomputed after every graph update.

The **Back** button undoes the most recent square move or permutation
replacement. Each successful graph-changing operation is saved, so the button
can be pressed repeatedly to return through earlier drawings in the session.
**Forward** restores drawings that were undone; a new move or replacement
starts a new branch and clears the forward history.

### Interactive boundary measurements

The sidebar includes a symbolic boundary-measurement editor. In **Face
variables** mode it uses Postnikov's clockwise face-weight convention. Since
the product of all face weights is one, the first face is the dependent
product-inverse coordinate and every other face variable can be renamed. In
**Edge variables** mode each graph edge has its own input; an empty input means
weight `1`. Right-click an edge to assign its parameter directly on the graph.

The graph's source set and an almost-perfect matching determine a perfect
orientation. Directed path sums, including the rational resolvent when the
orientation has directed cycles, produce the `k × n` boundary-measurement
matrix. The result selector switches between the matrix and all maximal
Plücker coordinates. Both are typeset as LaTeX, with separate copy buttons for
Julia and Macaulay2 text. The result sits below the graph in the same scrollable
workspace. The **Sources** field accepts a comma- or space-separated `k`-subset and
uses those columns as the identity pivot submatrix; a clear error is shown if
the chosen subset does not admit an almost-perfect matching. Rational
expressions are combined, common monomial factors are cancelled, products are
printed without redundant parentheses, and LaTeX uses ordinary stacked
fractions. Controls in the right panel select matrix versus Plücker output and
adjust both the result text size and result-panel height. Drag the divider at
the left edge of the right panel to resize it or hide it completely; the
divider remains available to reopen it. Double-clicking it toggles the panel.

The combinatorial face labels can be hidden independently with **Hide face
labels**. Right-clicking a face or edge opens its parameter menu; assigned
values are drawn directly on the corresponding face or edge and synchronized
with the weight editor. Clearing a parameter restores weight one. Face 1 is the
dependent reference face and has no independent face parameter.

Boundary-face labels stay at their face barycenters in the ordinary plabic
view. When the dual graph is shown, they move to the midpoint of the adjacent
boundary arc and are reused as the corresponding dual vertices, so a boundary
face is never labelled twice. Internal labels retain their interior positions.
A dual face is shaded dark gray exactly when it corresponds to a
black internal plabic vertex. Square labels are double-clickable, so square
moves remain available when the dual is shown or even when the plabic layer is
hidden.

Square moves play as continuous geometric morphs: vertices move smoothly,
contracted edges shorten until their endpoints coincide, new edges grow into
place, colors blend, and removed vertices collapse away. The motion passes
through trivalentization, color switches, same-color contractions, bivalent
removal, and final layout recomputation. Choose Instant, Fast, Slow, Very slow,
or Movie mode from **Mutation speed**. Face labels stay in a top overlay for the
entire morph and update only when the final graph is ready. Undo/redo stores only
completed graphs.

```julia
session = interactive_plabic_graph(G)
# The browser opens automatically. When finished:
close(session)
```

To start without constructing a graph or supplying a permutation in Julia, use
the no-argument workspace launcher:

```julia
session = interactive_session()
```

This opens a blank disk. Enter the decorated permutation in the browser and
press **Draw / replace graph**. Afterward, the current graph is available as
`session.graph[]` just as it is with `interactive_plabic_graph`.

Alternatively, press **Draw graph manually**, choose `n`, and place black and
white vertices inside the disk. In **Add edges** mode, click two vertices to
join them. Internal vertices can be dragged at any time and their incident edges
follow continuously. Boundary vertex `1` is at the bottom and the labels proceed
clockwise. The editor also supports undo, deletion, clearing, and cancellation.
Press **Reduce & use graph** to contract same-color edges, suppress bivalent
internal vertices, validate the planar drawing, and compute its decorated trip
permutation. If the remaining graph is still nonreduced, the package constructs
the canonical reduced graph for that same permutation. The result is transferred
to all the existing strand, face-label, dual, square-move, and boundary-measurement tools.

For a graph already constructed from a permutation, press **Edit current graph**
to reopen its current representative in the same editor. This also works after
square moves. Edit the layout or combinatorics, then press **Reduce & use graph**
to replace the interactive graph with the reduced result.

If a browser cannot be opened automatically, visit `session.url`. Choose a port
or suppress automatic opening with
`interactive_plabic_graph(G;port=9000,open_browser=false)`. The Julia process and
`session` must remain alive while the view is in use.

By default, drawing calls `layout_plabic_graph` with `method=:harmonic`. Boundary
vertices remain fixed on the circle, genuine lollipops receive compact inward
anchors, and every other internal vertex is placed at the barycenter of its
neighbors. This Tutte-style layout uses the middle of the disk more evenly and
does not depend on repulsive charges. A crossing check is always performed; an
unusual graph for which harmonic placement crosses falls back automatically to
the force layout.

The previous optimizer remains available with `method=:force` or
`layout_method=:force` in `draw_plabic_graph`. Its controls are `iterations`,
`restarts`, `repulsion`, `attraction`, `ideal_length`, `lollipop_length`, and
`gravity`.

```julia
G = plabic_graph([3, 4, 1, 2])
Gnice = layout_plabic_graph(G; method=:harmonic)
draw_plabic_graph(Gnice; optimize_layout=false)
```

The harmonic layout is deterministic. The optional force-directed optimization
is heuristic, so increasing its restart count searches more candidate minima.

## Known limitations

- Most legacy methods require `Vector{Int64}` exactly rather than generic integer
  vectors.
- Bases use mutable vectors as `Set` keys.
- Several specialized algorithms have only prototype-level tests.
- Face-labeled plabic graphs and scattering-amplitude operations remain
  Mathematica-only.
- The historical spelling `countExceedences` is retained for compatibility;
  `decorated_excedances` is the preferred public spelling.
