import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

/**
 * Interactive complete graph drawing.
 *
 * Compile: javac CompleteGraphApplet.java
 * Run:     java CompleteGraphApplet
 */
public final class CompleteGraphApplet {
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Complete Graph Point Applet");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(new AppPanel());
            frame.setMinimumSize(new Dimension(680, 480));
            frame.setSize(980, 700);
            frame.setLocationByPlatform(true);
            frame.setVisible(true);
        });
    }

    private static final class AppPanel extends JPanel {
        private final GraphCanvas canvas = new GraphCanvas();
        private final JSpinner nSpinner = new JSpinner(new SpinnerNumberModel(5, 1, 200, 1));
        private final JLabel status = new JLabel();

        AppPanel() {
            super(new BorderLayout(8, 8));
            setBorder(new EmptyBorder(10, 10, 10, 10));

            JPanel controls = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 0));
            controls.add(new JLabel("n:"));
            controls.add(nSpinner);

            JButton reset = new JButton("New");
            reset.addActionListener(e -> {
                canvas.setTargetPointCount((Integer) nSpinner.getValue());
                updateStatus();
            });
            controls.add(reset);

            JButton removeLast = new JButton("Undo");
            removeLast.addActionListener(e -> {
                canvas.removeLastPoint();
                updateStatus();
            });
            controls.add(removeLast);

            JButton randomize = new JButton("Randomize");
            randomize.addActionListener(e -> {
                canvas.randomizePoints((Integer) nSpinner.getValue());
                updateStatus();
            });
            controls.add(randomize);

            status.setBorder(new EmptyBorder(0, 8, 0, 0));
            controls.add(status);

            nSpinner.addChangeListener(e -> {
                canvas.setTargetPointCount((Integer) nSpinner.getValue());
                updateStatus();
            });
            canvas.setStateListener(this::updateStatus);

            add(controls, BorderLayout.NORTH);
            add(canvas, BorderLayout.CENTER);

            canvas.setTargetPointCount((Integer) nSpinner.getValue());
            updateStatus();
        }

        private void updateStatus() {
            int placed = canvas.pointCount();
            int target = canvas.targetPointCount();
            int lines = placed * (placed - 1) / 2;
            if (placed < target) {
                int remaining = target - placed;
                status.setText(placed + "/" + target + " points, " + remaining + " to place");
            } else {
                status.setText(placed + " points, " + lines + " lines");
            }
        }
    }

    private static final class GraphCanvas extends JPanel {
        private static final int POINT_RADIUS = 8;
        private static final int HIT_RADIUS = 13;

        private final List<Point2D.Double> points = new ArrayList<>();
        private int targetPointCount = 5;
        private int draggedPoint = -1;
        private Runnable stateListener = () -> {};

        GraphCanvas() {
            setOpaque(true);
            setBackground(new Color(248, 249, 252));
            setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));

            MouseAdapter mouse = new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    requestFocusInWindow();
                    draggedPoint = findPointAt(e.getX(), e.getY());
                    if (draggedPoint < 0 && SwingUtilities.isLeftMouseButton(e)
                            && points.size() < targetPointCount) {
                        points.add(new Point2D.Double(e.getX(), e.getY()));
                        draggedPoint = points.size() - 1;
                        stateListener.run();
                        repaint();
                    }
                }

                @Override
                public void mouseDragged(MouseEvent e) {
                    if (draggedPoint < 0) return;
                    Point2D.Double point = points.get(draggedPoint);
                    point.x = clamp(e.getX(), POINT_RADIUS, Math.max(POINT_RADIUS, getWidth() - POINT_RADIUS));
                    point.y = clamp(e.getY(), POINT_RADIUS, Math.max(POINT_RADIUS, getHeight() - POINT_RADIUS));
                    repaint();
                }

                @Override
                public void mouseReleased(MouseEvent e) {
                    draggedPoint = -1;
                }
            };
            addMouseListener(mouse);
            addMouseMotionListener(mouse);
        }

        void setStateListener(Runnable listener) {
            stateListener = listener;
        }

        int pointCount() {
            return points.size();
        }

        int targetPointCount() {
            return targetPointCount;
        }

        void setTargetPointCount(int n) {
            targetPointCount = n;
            points.clear();
            draggedPoint = -1;
            repaint();
            stateListener.run();
        }

        void removeLastPoint() {
            if (points.isEmpty()) return;
            points.remove(points.size() - 1);
            draggedPoint = -1;
            repaint();
            stateListener.run();
        }

        void randomizePoints(int n) {
            targetPointCount = n;
            points.clear();
            int width = Math.max(getWidth(), 320);
            int height = Math.max(getHeight(), 240);
            for (int i = 0; i < n; i++) {
                double x = POINT_RADIUS + Math.random() * Math.max(1, width - 2.0 * POINT_RADIUS);
                double y = POINT_RADIUS + Math.random() * Math.max(1, height - 2.0 * POINT_RADIUS);
                points.add(new Point2D.Double(x, y));
            }
            draggedPoint = -1;
            repaint();
            stateListener.run();
        }

        @Override
        protected void paintComponent(Graphics rawGraphics) {
            super.paintComponent(rawGraphics);
            Graphics2D g = (Graphics2D) rawGraphics.create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            drawGrid(g);
            drawLines(g);
            drawPoints(g);
            g.dispose();
        }

        private void drawGrid(Graphics2D g) {
            g.setColor(new Color(225, 229, 236));
            for (int x = 0; x < getWidth(); x += 40) g.drawLine(x, 0, x, getHeight());
            for (int y = 0; y < getHeight(); y += 40) g.drawLine(0, y, getWidth(), y);
        }

        private void drawLines(Graphics2D g) {
            g.setStroke(new BasicStroke(1.25f));
            g.setColor(new Color(52, 70, 95, 95));
            for (int i = 0; i < points.size(); i++) {
                Point2D.Double a = points.get(i);
                for (int j = i + 1; j < points.size(); j++) {
                    Point2D.Double b = points.get(j);
                    g.draw(new Line2D.Double(a, b));
                }
            }
        }

        private void drawPoints(Graphics2D g) {
            for (int i = 0; i < points.size(); i++) {
                Point2D.Double point = points.get(i);
                Ellipse2D circle = new Ellipse2D.Double(
                        point.x - POINT_RADIUS, point.y - POINT_RADIUS,
                        2.0 * POINT_RADIUS, 2.0 * POINT_RADIUS);
                g.setColor(i == draggedPoint ? new Color(238, 117, 71) : new Color(43, 126, 211));
                g.fill(circle);
                g.setColor(new Color(20, 30, 45));
                g.setStroke(new BasicStroke(1.8f));
                g.draw(circle);

                g.setFont(g.getFont().deriveFont(Font.BOLD, 11f));
                g.drawString(Integer.toString(i + 1),
                        (float) point.x + POINT_RADIUS + 3,
                        (float) point.y - POINT_RADIUS);
            }
        }

        private int findPointAt(double x, double y) {
            for (int i = points.size() - 1; i >= 0; i--) {
                if (points.get(i).distance(x, y) <= HIT_RADIUS) return i;
            }
            return -1;
        }

        private static double clamp(double value, double min, double max) {
            return Math.max(min, Math.min(max, value));
        }
    }
}
