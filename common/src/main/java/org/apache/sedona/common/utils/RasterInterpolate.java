package org.apache.sedona.common.utils;

import org.geotools.coverage.GridSampleDimension;
import org.geotools.coverage.grid.GridCoverage2D;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.index.strtree.STRtree;
import org.apache.commons.math3.linear.*;

import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.util.*;

public class RasterInterpolate {
    private RasterInterpolate() {}

    public static STRtree generateSTRtree(GridCoverage2D inputRaster, int band) {
        Raster rasterData = inputRaster.getRenderedImage().getData();
        int width = rasterData.getWidth();
        int height = rasterData.getHeight();
        Double noDataValue = RasterUtils.getNoDataValue(inputRaster.getSampleDimension(band));
        GeometryFactory geometryFactory = new GeometryFactory();
        STRtree rtree = new STRtree();

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double value = rasterData.getSampleDouble(x, y, band);
                if (!Double.isNaN(value) && value != noDataValue) {
                    Point point = geometryFactory.createPoint(new Coordinate(x, y));
                    RasterPoint rasterPoint = new RasterPoint(point, value, 0.0);
                    rtree.insert(new Envelope(point.getCoordinate()), rasterPoint);
                }
            }
        }
        rtree.build();
        return rtree;
    }

    public static double interpolateIDW(int x, int y, STRtree strtree, int width, int height, double power, String mode, Double numPointsOrRadius, Double maxRadiusOrMinPoints) {
        GeometryFactory geometryFactory = new GeometryFactory();
        PriorityQueue<RasterPoint> minHeap = new PriorityQueue<>(Comparator.comparingDouble(RasterPoint::getDistance));

        if (mode.equalsIgnoreCase("variable")) {
            Double numPoints = (numPointsOrRadius==null) ? 12:numPointsOrRadius; // Default no. of points -> 12
            Double maxRadius = (maxRadiusOrMinPoints==null) ? Math.sqrt((width*width)+(height*height)):maxRadiusOrMinPoints; // Default max radius -> diagonal of raster
            List<RasterPoint> queryResult = strtree.query(new Envelope(x - maxRadius, x + maxRadius, y - maxRadius, y + maxRadius));
            if (mode.equalsIgnoreCase("variable") && strtree.size() < numPointsOrRadius) {
                throw new IllegalArgumentException("Parameter 'numPoints' defaulted to 12 which is larger than no. of valid pixels within the max search radius. Please choose an appropriate value");
            }
            for (RasterPoint rasterPoint : queryResult) {
                if (numPoints<=0) {
                    break;
                }
                Point point = rasterPoint.getPoint();
                double distance = point.distance(geometryFactory.createPoint(new Coordinate(x, y)));
                rasterPoint.setDistance(distance);
                minHeap.add(rasterPoint);
                numPoints --;
            }
        } else if (mode.equalsIgnoreCase("fixed")) {
            Double radius = (numPointsOrRadius==null) ? Math.sqrt((width*width)+(height*height)):numPointsOrRadius; // Default radius -> diagonal of raster
            Double minPoints = (maxRadiusOrMinPoints==null) ? 0:maxRadiusOrMinPoints; // Default min no. of points -> 0
            List<RasterPoint> queryResult = new ArrayList<>();
            do {
                queryResult.clear();
                Envelope searchEnvelope = new Envelope(x - radius, x + radius, y - radius, y + radius);
                queryResult = strtree.query(searchEnvelope);
                // If minimum points requirement met, break the loop
                if (queryResult.size() >= minPoints) {
                    break;
                }
                radius *= 1.5; // Increase radius by 50%
            } while (true);
//            System.out.println("\nFor ("+x+" "+y+") -");
//            printQueryResult(queryResult);
            for (RasterPoint rasterPoint : queryResult) {
                Point point = rasterPoint.getPoint();
                double distance = point.distance(geometryFactory.createPoint(new Coordinate(x, y)));
                if (distance <= 0 || distance > radius) {
                    continue;
                }
                rasterPoint.setDistance(distance);
                minHeap.add(rasterPoint);
            }
        }

        double numerator = 0.0;
        double denominator = 0.0;

        while (!minHeap.isEmpty()) {
            RasterPoint rasterPoint = minHeap.poll();
            double value = rasterPoint.getValue();
            double distance = rasterPoint.getDistance();
            double weight = 1.0 / Math.pow(distance, power);
            numerator += weight * value;
            denominator += weight;
        }

        double interpolatedValue = (denominator > 0 ? numerator / denominator : Double.NaN);
        return interpolatedValue;
    }

    public static void printQueryResult(List<RasterPoint> list) {
        String res = "";
        for (RasterPoint rasterPoint : list) {
            res += rasterPoint.getPoint().toString() + ", ";
        }
        System.out.println(res);
    }

    public static GridCoverage2D interpolateSpline(GridCoverage2D inputRaster, Integer band) throws IllegalArgumentException {
        // Extract valid points from raster data
        List<RasterPoint> validPoints = extractValidPoints(inputRaster, band);
        System.out.println("\nvalidPoints size: "+validPoints.size());
        printQueryResult(validPoints);

        // Initialize the spline interpolator with the valid points
        SplineInterpolation splineInterpolation = new SplineInterpolation(validPoints, 0.5);

        // Get the raster data and prepare for writing interpolated values
        Raster rasterData = inputRaster.getRenderedImage().getData();
        WritableRaster writableRaster = rasterData.createCompatibleWritableRaster(rasterData.getWidth(), rasterData.getHeight());
        Double noDataValue = RasterUtils.getNoDataValue(inputRaster.getSampleDimension(band));

        // Interpolate values for each point in the raster
        for (int y = 0; y < writableRaster.getHeight(); y++) {
            for (int x = 0; x < writableRaster.getWidth(); x++) {
                // Check if the current raster value is a noDataValue
                double value = rasterData.getSampleDouble(x, y, band);
                if (Double.isNaN(value) || value == noDataValue) {
                    // Use the spline interpolator to get the interpolated value
                    double interpolatedValue = splineInterpolation.interpolate(x, y);
                    // Set the interpolated value in the writable raster
                    writableRaster.setSample(x, y, band, interpolatedValue);
                } else {
                    writableRaster.setSample(x, y, band, value);
                }
            }
        }

        // Create and return the new GridCoverage2D object with the interpolated raster
        GridSampleDimension[] gridSampleDimensions = inputRaster.getSampleDimensions();
        return RasterUtils.clone(writableRaster, inputRaster.getGridGeometry(), gridSampleDimensions, inputRaster, null, true);
    }

    public static List<RasterPoint> extractValidPoints(GridCoverage2D inputRaster, Integer band) {
        List<RasterPoint> validPoints = new ArrayList<>();
        Raster rasterData = inputRaster.getRenderedImage().getData();
        int width = rasterData.getWidth();
        int height = rasterData.getHeight();
        GeometryFactory geometryFactory = new GeometryFactory();

        Double noDataValue = RasterUtils.getNoDataValue(inputRaster.getSampleDimension(band));

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double value = rasterData.getSampleDouble(x, y, band);
                if (!Double.isNaN(value) && value != noDataValue) {
                    Point point = geometryFactory.createPoint(new Coordinate(x, y));
                    validPoints.add(new RasterPoint(point, value));
                }
            }
        }
//        System.out.println("validPoints size: "+validPoints.size());
        return validPoints;
    }

    public static class RasterPoint {
        private Point point; // JTS Point
        private double value; // The associated value
        private double distance; // Distance measure

        public RasterPoint(Point point, double value, double distance) {
            this.point = point;
            this.value = value;
            this.distance = distance;
        }
        public RasterPoint(Point point, double value) {
            this.point = point;
            this.value = value;
            this.distance = 0.0;
        }

        public Point getPoint() {
            return point;
        }

        public double getX() {
            return point.getX();
        }

        public double getY() {
            return point.getY();
        }

        public double getValue() {
            return value;
        }

        public double getDistance() {
            return distance;
        }

        public void setDistance(double distance) {
            this.distance = distance;
        }
    }
}
