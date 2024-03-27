package org.apache.sedona.common.utils;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.special.BesselJ;
import org.apache.sedona.common.utils.RasterInterpolate.RasterPoint;
import java.util.List;

public class SplineInterpolation {
    private RealVector coefficients; // Coefficients for the spline
    private List<RasterPoint> points;
    private RealVector polynomialCoefficients;
    private double tau; // The weight parameter for regularized spline

    public SplineInterpolation(List<RasterPoint> points) {
        this(points, 0.0);
    }
    public SplineInterpolation(List<RasterPoint> points, double tau) {
        this.points = points;
        this.tau = tau;
        this.coefficients = calculateCoefficients(points); // Make sure this is implemented
        // Extract polynomial coefficients, assuming they are the last 3 elements
        int n = coefficients.getDimension();
        this.polynomialCoefficients = coefficients.getSubVector(n - 3, 3);
    }

    private RealVector calculateCoefficients(List<RasterPoint> points) {
        int n = points.size();
        RealMatrix K = new Array2DRowRealMatrix(n, n);
        RealMatrix P = new Array2DRowRealMatrix(n, 3);
        RealVector V = new ArrayRealVector(n);

        // Construct matrices K, P, and vector V
        for (int i = 0; i < n; i++) {
            RasterPoint pi = points.get(i);
            V.setEntry(i, pi.getValue());
            P.setEntry(i, 0, 1);
            P.setEntry(i, 1, pi.getX());
            P.setEntry(i, 2, pi.getY());

            for (int j = 0; j < n; j++) {
                RasterPoint pj = points.get(j);
                double r = Math.hypot(pi.getX() - pj.getX(), pi.getY() - pj.getY());
                K.setEntry(i, j, r == 0 ? 0 : r * r * Math.log(r));
            }
        }

        // Combine K and P to form system matrix A and vector b
        RealMatrix A = MatrixUtils.createRealMatrix(n + 3, n + 3);
        A.setSubMatrix(K.getData(), 0, 0);
        A.setSubMatrix(P.getData(), 0, n);
        A.setSubMatrix(P.transpose().getData(), n, 0);
        RealVector b = new ArrayRealVector(n + 3);
        b.setSubVector(0, V);

        // Apply the constraints for the polynomial part (last three rows)
        for (int i = n; i < n + 3; i++) {
            b.setEntry(i, 0);
        }

        // Perform SVD on matrix A
        SingularValueDecomposition svd = new SingularValueDecomposition(A);

        // Solve Ax = b using the pseudo-inverse of A
        RealVector x = svd.getSolver().solve(b);

        // The first n entries are the spline coefficients, the remaining are the polynomial coefficients
        return x.getSubVector(0, n);
    }

    // Method to interpolate a value at a given (x, y)
    // Implement the actual interpolation formula
    public double interpolate(double x, double y) {
        // Calculate T(x, y) using polynomial coefficients
        double Txy = polynomialCoefficients.getEntry(0) // a1
                + polynomialCoefficients.getEntry(1) * x // a2 * x
                + polynomialCoefficients.getEntry(2) * y; // a3 * y

        // Initialize the sum for the radial basis function
        double sumR = 0;

        // Loop over all points to calculate the R(r) function part
        for (int i = 0; i < points.size(); i++) {
            RasterPoint p = points.get(i);
            double r = Math.hypot(x - p.getX(), y - p.getY()); // distance from (x,y) to the point p

            // Avoid log of 0 if r is 0
            if (r != 0) {
                // R(r) for regularized option with the modified Bessel function of the second kind
                // and assuming r is already squared in the following formula
                double rr = r * r;
                double logTerm = rr * Math.log(r);
                double besselTerm = rr * (Math.log(rr) - 1) * BesselJ.value(0, Math.sqrt(rr * tau));
                double Rr = logTerm - tau * besselTerm; // using r^2 instead of r in the formula

                sumR += coefficients.getEntry(i) * Rr; // add weighted R(r)
            }
        }

        // The interpolated value is the sum of T(x,y) and the sum of weighted R(r) values
        return Txy + sumR;
    }
}
