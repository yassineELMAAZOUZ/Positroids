import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Interactive nearest-neighbor coloring in the hyperbolic upper half-plane.
 *
 * Compile: javac HyperbolicColoring.java
 * Run:     java HyperbolicColoring
 */
public final class HyperbolicColoring {
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Hyperbolic Nearest-Point Coloring");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(new AppPanel());
            frame.setMinimumSize(new Dimension(760, 560));
            frame.setSize(1050, 760);
            frame.setLocationByPlatform(true);
            frame.setVisible(true);
        });
    }

    private static final class AppPanel extends JPanel {
        private final PlaneCanvas canvas = new PlaneCanvas();
        private final JSpinner kSpinner = new JSpinner(new SpinnerNumberModel(6, 1, 100, 1));
        private final JLabel status = new JLabel();

        AppPanel() {
            super(new BorderLayout(8, 8));
            setBorder(new EmptyBorder(10, 10, 10, 10));

            JPanel controls = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 0));
            controls.add(new JLabel("Initial colored points k:"));
            controls.add(kSpinner);

            JButton restart = new JButton("Restart");
            restart.addActionListener(e -> {
                canvas.restart((Integer) kSpinner.getValue());
                updateStatus();
            });
            controls.add(restart);

            JButton undo = new JButton("Undo");
            undo.addActionListener(e -> {
                canvas.undo();
                updateStatus();
            });
            controls.add(undo);

            JButton fit = new JButton("Reset view");
            fit.addActionListener(e -> canvas.resetView());
            controls.add(fit);

            status.setBorder(new EmptyBorder(0, 8, 0, 0));
            controls.add(status);

            canvas.setStateListener(this::updateStatus);
            add(controls, BorderLayout.NORTH);
            add(canvas, BorderLayout.CENTER);

            JLabel help = new JLabel(
                    "Click to add points. Mouse wheel zooms; drag horizontally to pan. "
                            + "The first k points get distinct colors; every later point inherits "
                            + "the color of its nearest earlier point.");
            help.setBorder(new EmptyBorder(2, 4, 0, 4));
            add(help, BorderLayout.SOUTH);

            canvas.restart((Integer) kSpinner.getValue());
            updateStatus();
        }

        private void updateStatus() {
            int n = canvas.pointCount();
            int k = canvas.initialCount();
            if (n < k) {
                status.setText("Choose " + (k - n) + " initial point" + (k - n == 1 ? "" : "s") + ".");
            } else {
                status.setText(n + " points — new clicks inherit a nearest-neighbor color.");
            }
        }
    }

    private static final class ColoredPoint {
        final double x;
        final double y;
        final Color color;
        final int insertionOrder;

        ColoredPoint(double x, double y, Color color, int insertionOrder) {
            this.x = x;
            this.y = y;
            this.color = color;
            this.insertionOrder = insertionOrder;
        }
    }

    /**
     * Exact dynamic nearest-neighbor index:
     * a VP-tree for the stable points plus a short linear insertion buffer.
     */
    private static final class DynamicVPTree {
        private final List<ColoredPoint> stable = new ArrayList<>();
        private final List<ColoredPoint> buffer = new ArrayList<>();
        private VPNode root;

        void clear() {
            stable.clear();
            buffer.clear();
            root = null;
        }

        void add(ColoredPoint point) {
            buffer.add(point);
            int rebuildAt = Math.max(32, Math.max(1, stable.size() / 4));
            if (buffer.size() >= rebuildAt) rebuild();
        }

        ColoredPoint nearest(double x, double y) {
            Best best = new Best();
            search(root, x, y, best);
            for (ColoredPoint point : buffer) best.consider(point, distance(x, y, point.x, point.y));
            return best.point;
        }

        private void rebuild() {
            stable.addAll(buffer);
            buffer.clear();
            root = build(new ArrayList<>(stable));
        }

        private static VPNode build(List<ColoredPoint> points) {
            if (points.isEmpty()) return null;
            ColoredPoint vantage = points.remove(points.size() - 1);
            if (points.isEmpty()) return new VPNode(vantage, 0.0, null, null);

            points.sort(Comparator.comparingDouble(p ->
                    distance(vantage.x, vantage.y, p.x, p.y)));
            int middle = points.size() / 2;
            double threshold = distance(
                    vantage.x, vantage.y, points.get(middle).x, points.get(middle).y);
            List<ColoredPoint> inside = new ArrayList<>(points.subList(0, middle));
            List<ColoredPoint> outside = new ArrayList<>(points.subList(middle, points.size()));
            return new VPNode(vantage, threshold, build(inside), build(outside));
        }

        private static void search(VPNode node, double x, double y, Best best) {
            if (node == null) return;
            double d = distance(x, y, node.point.x, node.point.y);
            best.consider(node.point, d);

            if (node.inside == null && node.outside == null) return;
            if (d < node.threshold) {
                if (d - best.distance <= node.threshold) search(node.inside, x, y, best);
                if (d + best.distance >= node.threshold) search(node.outside, x, y, best);
            } else {
                if (d + best.distance >= node.threshold) search(node.outside, x, y, best);
                if (d - best.distance <= node.threshold) search(node.inside, x, y, best);
            }
        }

        private static final class VPNode {
            final ColoredPoint point;
            final double threshold;
            final VPNode inside;
            final VPNode outside;

            VPNode(ColoredPoint point, double threshold, VPNode inside, VPNode outside) {
                this.point = point;
                this.threshold = threshold;
                this.inside = inside;
                this.outside = outside;
            }
        }

        private static final class Best {
            ColoredPoint point;
            double distance = Double.POSITIVE_INFINITY;

            void consider(ColoredPoint candidate, double candidateDistance) {
                if (candidateDistance < distance
                        || (candidateDistance == distance
                        && (point == null || candidate.insertionOrder < point.insertionOrder))) {
                    point = candidate;
                    distance = candidateDistance;
                }
            }
        }
    }

    private static final class PlaneCanvas extends JPanel {
        private static final int BOTTOM_MARGIN = 28;
        private static final int POINT_RADIUS = 7;
        private static final Color BACKGROUND = new Color(248, 249, 252);

        private final List<ColoredPoint> points = new ArrayList<>();
        private final DynamicVPTree index = new DynamicVPTree();
        private int initialCount = 6;
        private double pixelsPerUnit = 90.0;
        private double centerX = 0.0;
        private int dragStartX;
        private double dragStartCenterX;
        private Runnable stateListener = () -> {};
        private BufferedImage voronoiImage;
        private double voronoiCenterX;
        private double voronoiScale;

        PlaneCanvas() {
            setOpaque(true);
            setBackground(BACKGROUND);
            setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));

            MouseAdapter mouse = new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    if (SwingUtilities.isLeftMouseButton(e)) {
                        dragStartX = e.getX();
                        dragStartCenterX = centerX;
                    }
                }

                @Override
                public void mouseDragged(MouseEvent e) {
                    if ((e.getModifiersEx() & MouseEvent.BUTTON1_DOWN_MASK) != 0) {
                        int delta = e.getX() - dragStartX;
                        if (Math.abs(delta) > 4) {
                            centerX = dragStartCenterX - delta / pixelsPerUnit;
                            invalidateVoronoi();
                            repaint();
                        }
                    }
                }

                @Override
                public void mouseReleased(MouseEvent e) {
                    if (!SwingUtilities.isLeftMouseButton(e)) return;
                    if (Math.abs(e.getX() - dragStartX) <= 4) addAtScreen(e.getX(), e.getY());
                }

                @Override
                public void mouseWheelMoved(MouseWheelEvent e) {
                    double oldScale = pixelsPerUnit;
                    double worldXAtMouse = screenToWorldX(e.getX());
                    pixelsPerUnit *= Math.pow(1.14, -e.getPreciseWheelRotation());
                    pixelsPerUnit = Math.max(18.0, Math.min(600.0, pixelsPerUnit));
                    centerX = worldXAtMouse - (e.getX() - getWidth() / 2.0) / pixelsPerUnit;
                    if (oldScale != pixelsPerUnit) {
                        invalidateVoronoi();
                        repaint();
                    }
                }
            };
            addMouseListener(mouse);
            addMouseMotionListener(mouse);
            addMouseWheelListener(mouse);
        }

        void setStateListener(Runnable listener) {
            stateListener = listener;
        }

        int pointCount() {
            return points.size();
        }

        int initialCount() {
            return initialCount;
        }

        void restart(int k) {
            initialCount = k;
            points.clear();
            index.clear();
            invalidateVoronoi();
            repaint();
            stateListener.run();
        }

        void undo() {
            if (points.isEmpty()) return;
            points.remove(points.size() - 1);
            rebuildIndex();
            invalidateVoronoi();
            repaint();
            stateListener.run();
        }

        void resetView() {
            pixelsPerUnit = 90.0;
            centerX = 0.0;
            invalidateVoronoi();
            repaint();
        }

        private void rebuildIndex() {
            index.clear();
            for (ColoredPoint point : points) index.add(point);
        }

        private void addAtScreen(int sx, int sy) {
            double x = screenToWorldX(sx);
            double y = screenToWorldY(sy);
            if (!(y > 0.0)) return;

            Color color;
            if (points.size() < initialCount) {
                color = paletteColor(points.size(), initialCount);
            } else {
                ColoredPoint nearest = index.nearest(x, y);
                if (nearest == null) return;
                color = nearest.color;
            }

            ColoredPoint point = new ColoredPoint(x, y, color, points.size());
            points.add(point);
            index.add(point);
            invalidateVoronoi();
            repaint();
            stateListener.run();
        }

        private void invalidateVoronoi() {
            voronoiImage = null;
        }

        private double screenToWorldX(double sx) {
            return centerX + (sx - getWidth() / 2.0) / pixelsPerUnit;
        }

        private double screenToWorldY(double sy) {
            return (getHeight() - BOTTOM_MARGIN - sy) / pixelsPerUnit;
        }

        private double worldToScreenX(double x) {
            return getWidth() / 2.0 + (x - centerX) * pixelsPerUnit;
        }

        private double worldToScreenY(double y) {
            return getHeight() - BOTTOM_MARGIN - y * pixelsPerUnit;
        }

        @Override
        protected void paintComponent(Graphics rawGraphics) {
            super.paintComponent(rawGraphics);
            Graphics2D g = (Graphics2D) rawGraphics.create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            drawVoronoiRegions(g);
            drawBoundary(g);
            for (int i = 0; i < points.size(); i++) drawPoint(g, points.get(i), i);
            g.dispose();
        }

        private void drawVoronoiRegions(Graphics2D g) {
            int width = getWidth();
            int planeHeight = getHeight() - BOTTOM_MARGIN;
            if (points.isEmpty() || width <= 0 || planeHeight <= 0) return;

            if (voronoiImage == null
                    || voronoiImage.getWidth() != width
                    || voronoiImage.getHeight() != planeHeight
                    || voronoiCenterX != centerX
                    || voronoiScale != pixelsPerUnit) {
                voronoiImage = renderVoronoi(width, planeHeight);
                voronoiCenterX = centerX;
                voronoiScale = pixelsPerUnit;
            }
            g.drawImage(voronoiImage, 0, 0, null);
        }

        private BufferedImage renderVoronoi(int width, int height) {
            BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            int[] tintedColors = new int[points.size()];
            for (int i = 0; i < points.size(); i++) {
                tintedColors[i] = tint(points.get(i).color, BACKGROUND, 0.34).getRGB();
            }

            /*
             * For a fixed query (x,y), minimizing hyperbolic distance is
             * equivalent to minimizing
             *
             *     ((x-x_i)^2 + (y-y_i)^2) / y_i.
             *
             * This avoids acosh and square roots while producing exactly the
             * same Voronoi regions.
             */
            for (int sy = 0; sy < height; sy++) {
                double y = screenToWorldY(sy + 0.5);
                for (int sx = 0; sx < width; sx++) {
                    double x = screenToWorldX(sx + 0.5);
                    int nearestIndex = 0;
                    double bestScore = Double.POSITIVE_INFINITY;

                    for (int i = 0; i < points.size(); i++) {
                        ColoredPoint point = points.get(i);
                        double dx = x - point.x;
                        double dy = y - point.y;
                        double score = (dx * dx + dy * dy) / point.y;
                        if (score < bestScore) {
                            bestScore = score;
                            nearestIndex = i;
                        }
                    }
                    image.setRGB(sx, sy, tintedColors[nearestIndex]);
                }
            }
            return image;
        }

        private void drawBoundary(Graphics2D g) {
            int y = getHeight() - BOTTOM_MARGIN;
            g.setColor(new Color(52, 58, 72));
            g.setStroke(new BasicStroke(2f));
            g.drawLine(0, y, getWidth(), y);
            g.setFont(g.getFont().deriveFont(Font.ITALIC, 12f));
            g.drawString("boundary y = 0", 8, y + 18);
        }

        private void drawPoint(Graphics2D g, ColoredPoint point, int number) {
            double sx = worldToScreenX(point.x);
            double sy = worldToScreenY(point.y);
            if (sx < -20 || sx > getWidth() + 20 || sy < -20 || sy > getHeight() + 20) return;

            Ellipse2D circle = new Ellipse2D.Double(
                    sx - POINT_RADIUS, sy - POINT_RADIUS,
                    2.0 * POINT_RADIUS, 2.0 * POINT_RADIUS);
            g.setColor(point.color);
            g.fill(circle);
            g.setColor(new Color(30, 33, 40));
            g.setStroke(new BasicStroke(number < initialCount ? 2.2f : 1.2f));
            g.draw(circle);

            if (number < initialCount) {
                g.setFont(g.getFont().deriveFont(Font.BOLD, 11f));
                g.drawString(Integer.toString(number + 1),
                        (float) sx + POINT_RADIUS + 3, (float) sy - POINT_RADIUS);
            }
        }
    }

    /**
     * Hyperbolic distance in the upper-half-plane model:
     * acosh(1 + ((x1-x2)^2 + (y1-y2)^2)/(2*y1*y2)).
     */
    static double distance(double x1, double y1, double x2, double y2) {
        double dx = x1 - x2;
        double dy = y1 - y2;
        double coshDistance = 1.0 + (dx * dx + dy * dy) / (2.0 * y1 * y2);
        coshDistance = Math.max(1.0, coshDistance);
        return Math.log(coshDistance + Math.sqrt(coshDistance * coshDistance - 1.0));
    }

    private static Color paletteColor(int index, int total) {
        float hue = (float) ((index * 0.6180339887498949) % 1.0);
        float saturation = total <= 2 ? 0.78f : 0.68f;
        return Color.getHSBColor(hue, saturation, 0.92f);
    }

    private static Color tint(Color foreground, Color background, double amount) {
        double inverse = 1.0 - amount;
        int red = (int) Math.round(amount * foreground.getRed() + inverse * background.getRed());
        int green = (int) Math.round(amount * foreground.getGreen() + inverse * background.getGreen());
        int blue = (int) Math.round(amount * foreground.getBlue() + inverse * background.getBlue());
        return new Color(red, green, blue);
    }

}
