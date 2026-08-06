# Hyperbolic nearest-point coloring

This is a dependency-free Java Swing app for the hyperbolic upper half-plane.

The folder also contains `CompleteGraphApplet.java`, a dependency-free Swing
applet for placing `n` points in the plane and drawing the complete graph on
those points. Points remain draggable after the lines are drawn.

```sh
javac CompleteGraphApplet.java
java CompleteGraphApplet
```

## Run it

```sh
javac HyperbolicColoring.java
java HyperbolicColoring
```

Choose `k`, click **Restart**, and then:

1. Click `k` points. They receive distinct colors.
2. Every later click receives the color of its nearest previously chosen point
   in the hyperbolic metric.
3. Each new point then joins the set used by later clicks.
4. The colored background shows the current hyperbolic Voronoi regions: it
   displays which color a new point would receive at every visible location.
   The regions update after each new point.

Use the mouse wheel to zoom, drag horizontally to pan, **Undo** to remove the
last point, and **Reset view** to restore the initial view.

## Distance and search

For points `(x1,y1)` and `(x2,y2)`, where both `y` coordinates are positive,
the app uses

```text
d = acosh(1 + ((x1-x2)^2 + (y1-y2)^2)/(2*y1*y2)).
```

Nearest-neighbor queries use an exact VP-tree. New points first enter a small
buffer, and the tree is rebuilt in batches. This avoids rebuilding after every
click while preserving exact answers. The displayed Voronoi background is
rendered at full screen-pixel resolution; color assignment for clicked points
is exact.
