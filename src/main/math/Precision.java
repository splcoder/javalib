package main.math;

public class Precision {
	public static final double PRECISION = 1e-12;
	public static boolean equals( double _1, double _2 ){
		double dif = _1 - _2;
		return dif < 0. ? -dif <= PRECISION : dif <= PRECISION;
	}
	public static boolean isZero( double v ){
		return v < 0. ? -v <= PRECISION : v <= PRECISION;
	}
}
