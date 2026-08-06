# Mathematica migration inventory

The historical `positroids.m` package combines several domains in one large
source file. A literal line-by-line translation would couple the Julia core to
notebook formatting and amplitude-specific global state, so migration is grouped
by mathematical subsystem.

| Mathematica area / representative functions | Julia status |
|---|---|
| Legal/decorated permutations (`legalPermQ`, `decorate`, `permK`, `rotate`) | Partial: signed decorated-permutation convention, validation, rank, and cyclic helpers |
| Grassmann necklaces (`necklace`, `necklaceR`) | Implemented as minimum, maximum, and reverse necklace routines |
| Positroid dimension and incidence (`dimension`, `boundary`, `inverseBoundary`) | Dimension implemented; full boundary traversal not ported |
| Matrix/cell conversion (`matToPerm`, `permToMatrix`, `preferredGauge`) | Positive bridge matrix parametrization ported; inverse matrix recognition remains |
| Positivity (`positiveQ`) | Ported and clarified as `is_totally_nonnegative` and `is_totally_positive` |
| Bridge coordinates (`transpositionChain`, `bridgeToMinors`) | Bridge and normalized boundary-measurement charts ported; symbolic minor formulas remain |
| Plabic graph data and drawing (`plabicGraphData`, `plabicGraph`) | Core bridge construction and disk drawing ported; face labels, trip overlays, and optimized layout remain |
| Twists | Experimental Julia matrix twists are present; convention parity needs regression fixtures |
| Residues and scattering amplitudes (`permToResidue`, BCFW and contour functions) | Out of the initial Julia-core scope |
| Notebook presentation/export helpers | Intentionally not part of the mathematical API |

## Recommended next migration slices

1. Add an immutable `DecoratedPermutation` type and bounded-affine conversion,
   making the Mathematica and Julia conventions explicit.
2. Port `boundary`/`inverseBoundary` with cross-language fixtures from the demo
   notebook.
3. Add maximal-minor positivity and matrix-to-cell recognition.
4. Introduce typed plabic graphs, then port bridge decompositions and drawing.
5. Keep scattering-amplitude routines in an optional extension package so the
   combinatorics core remains lightweight.

The original `.m` and `.nb` files are retained unchanged as migration fixtures.
