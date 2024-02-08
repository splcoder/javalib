package main.math;

public class Complex implements Comparable<Complex> {
	public static final Complex ONE	= new Complex( 1., 0. );
	public static final Complex I	= new Complex( 0., 1. );
	public static final Complex _I	= new Complex( 0., -1. );
	
	private final double real, imag;

	// Helpers for a high intensive use
	private Double abs = null, phase = null;

	public Complex( double real, double imag ) {
		this.real = real;
		this.imag = imag;
	}

	public Complex( double real ) {
		this.real = real;
		this.imag = 0.;
	}
	private Complex( double real, double imag, Double abs, Double phase ){
		this.real = real;
		this.imag = imag;
		this.abs = (abs == null ? null : (double)abs);
		this.phase = (phase == null ? null : (double)phase);
	}
	public Complex copy(){
		return new Complex( real, imag, abs, phase );
	}
	public static Complex fromComplex( Complex a ){	// Copy only the basic variables
		return new Complex( a.real, a.imag );
	}
	// Conversion ------------------------------------------------------------------------------------------------------
	public static Complex fromDouble( double real ){
		return new Complex( real, 0. );
	}
	public static Complex fromDouble( double real, double imag ){
		return new Complex( real, imag );
	}
	public static Complex fromPolar( double magnitude, double phase ){
		return new Complex( magnitude * Math.cos( phase ), magnitude * Math.sin( phase ) );
	}
	/**
	 * Format: "real,imag"
	 */
	public static Complex parse( String str ) throws NumberFormatException {
		if( str == null )	throw new NumberFormatException( "Complex.parse -> null argument" );
		double rp = 0, ip = 0;
		String[] arr = str.split( "," );
		switch( arr.length ){
			case 0:	break;	// ","
			case 2: {
				String aux = arr[1];
				ip = aux.isEmpty() ? 0. : Double.parseDouble( aux );
			}
			case 1: {
				String aux = arr[0];
				rp = aux.isEmpty() ? 0. : Double.parseDouble( aux );
				break;
			}
			default:	throw new NumberFormatException( "Complex.parse -> there must be at max 2 components (max 1 comma)" );
		}
		return new Complex( rp, ip );
	}

	public String toString() {
		if( imag == 0 )	return real + "";
		if( real == 0 )	return imag + "i";
		return real + (imag < 0. ? (" - " + (-imag) ) : (" + " + imag)) + "i";
	}

	// Basic fast access values ----------------------------------------------------------------------------------------

	public double real(){ return real; }
	public double imag(){ return imag; }
	
	public boolean isReal(){ return Precision.isZero( imag ); }
	public boolean isImaginary(){ return ! Precision.isZero( imag ); }
	public boolean isZero(){ return Precision.isZero( real ) && Precision.isZero( imag ); }
	public boolean isINF(){ return Double.isInfinite( real ) || Double.isInfinite( imag ); }
	public boolean isNaN(){ return Double.isNaN( real ) || Double.isNaN( imag ); }
	/**
	 * Returns the absolute/modulus/magnitude value of this complex number.
	 */
	public double abs() {
		//return Math.hypot( real, imag );
		if( abs == null )	abs = Math.hypot( real, imag );
		return abs;
	}
	/**
	 * abs^2
	 */
	public double norm() {
		return real * real + imag * imag;
	}
	/**
	 * Returns the phase/angle/argument of this complex number, in (-PI, PI]
	 */
	public double phase() {
		//return Math.atan2( imag, real );
		if( phase == null )	phase = Math.atan2( imag, real );
		return phase;
	}
	public Complex conjugate() {
		return new Complex( real, -imag );
	}
	public Complex neg() {
		return new Complex( -real, -imag );
	}
	/**
	 * Returns the reciprocal of this complex number: (1 / this)
	 */
	public Complex reciprocal() {
		double scale = real * real + imag * imag;
		return new Complex(real / scale, -imag / scale );
	}
	public Complex square() {
		return new Complex( this.real * this.real - this.imag * this.imag, this.real * this.imag + this.imag * this.real );
	}

	// Fast access for basic operations (+ - * /) ----------------------------------------------------------------------

	public Complex add( Complex a ){ return new Complex( real + a.real, imag + a.imag ); }
	public Complex add( double a ){ return new Complex( real + a, imag ); }
	public Complex sub( Complex a ){ return new Complex( real - a.real, imag - a.imag ); }
	public Complex sub( double a ){ return new Complex( real - a, imag ); }
	public Complex mul( Complex a ){ return new Complex( real * a.real - imag * a.imag, real * a.imag + imag * a.real ); }
	public Complex mul( double a ){ return new Complex( real * a, imag * a ); }
	public Complex mulI( double a ){ return new Complex( -imag * a, real * a ); }		// this * (i*a)
	public Complex div( Complex a ) {
		double n = a.norm();
		return new Complex( (real * a.real + imag * a.imag) / n, (imag * a.real - real * a.imag) / n );
	}
	public Complex div( double a ){ return new Complex( real / a, imag / a ); }

	// Addition and Subtraction ----------------------------------------------------------------------------------------

	public static Complex add( Complex a, Complex b ) {
		return new Complex( a.real + b.real, a.imag + b.imag );
	}
	public static Complex add( Complex a, double b ) {
		return new Complex( a.real + b, a.imag );
	}
	public static Complex add( double a, Complex b ) {
		return new Complex( a + b.real, b.imag );
	}
	public static Complex add( double a, double b ) {
		return new Complex( a + b, 0. );
	}
	public static Complex sub( Complex a, Complex b ) {
		return new Complex( a.real - b.real, a.imag - b.imag );
	}
	public static Complex sub( Complex a, double b ) {
		return new Complex( a.real - b, a.imag );
	}
	public static Complex sub( double a, Complex b ) {
		return new Complex( a - b.real, -b.imag );
	}
	public static Complex sub( double a, double b ) {
		return new Complex( a - b, 0. );
	}

	// Multiplication and Division -------------------------------------------------------------------------------------

	public static Complex mul( Complex a, Complex b ) {
		return new Complex( a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real );
	}
	public static Complex mul( Complex a, double b ) {
		return new Complex( a.real * b, a.imag * b );
	}
	public static Complex mul( double a, Complex b ) {
		return new Complex( a * b.real, a * b.imag );
	}
	public static Complex mul( double a, double b ) {
		return new Complex( a * b, 0. );
	}
	public static Complex mulI( Complex a, Complex b ) {
		return new Complex( -(a.real * b.imag + a.imag * b.real), a.real * b.real - a.imag * b.imag );
	}
	public static Complex mulI( Complex a, double b ) {	// a * (i*b)
		return new Complex( -a.imag * b, a.real * b );
	}
	public static Complex mulI( Complex a ) {	// a * (i*b)
		return new Complex( -a.imag, a.real );
	}
	public static Complex div( Complex a, Complex b ) {
		//return Complex.mul( a, b.reciprocal() );
		double n = b.norm();
		return new Complex( (a.real * b.real + a.imag * b.imag) / n, (a.imag * b.real - a.real * b.imag) / n );
	}
	public static Complex div( Complex a, double b ) {
		return new Complex( a.real / b, a.imag / b );
	}
	public static Complex div( double a, Complex b ) {
		double n = b.norm();
		return new Complex( a * b.real / n,  -a * b.imag / n );
	}
	public static Complex div( double a, double b ) {
		return new Complex( a / b, 0. );
	}

	// Power and Root --------------------------------------------------------------------------------------------------

	public static Complex pow( Complex base, int exp ){
		//double modulus = base.abs(), angle = base.phase();
		//return Complex.fromPolar( Math.pow( modulus, exp ), angle * exp );
		// Better
		double n = base.norm(), angle = base.phase();	// Integer exponent: it doesn't affect to sin/cos
		return Complex.fromPolar( Math.pow( n, 0.5 * exp ), angle * exp );
	}
	/**
	 * https://en.wikipedia.org/wiki/Complex_number	>>> Integer and fractional exponents
	 *		The nth root is a n-valued function of z
	 *		i.e. (z^n)^(1/n) != z, since the left-hand side consists of n values, and the right-hand side is a single value
	 * @param k		branch of the polar complex number: (k*2*PI)
	 */
	public static Complex pow( Complex base, Complex exp, int k ){
		// Set base in polar form, and exp in cartesian
		//double modulus = base.abs(), angle = base.phase() + Constants.M_2PI * k;
		//return Complex.fromPolar( Math.pow( modulus, exp.real ) * Math.exp( -angle * exp.imag ), angle * exp.real + exp.imag * Math.log( modulus ) );
		// Better
		double n = base.norm(), angle = base.phase() + Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( n, 0.5 * exp.real ) * Math.exp( -angle * exp.imag ), angle * exp.real + 0.5 * exp.imag * Math.log( n ) );
	}
	public static Complex pow( Complex base, Complex exp ){
		return Complex.pow( base, exp, 0 );
	}
	public static Complex pow( Complex base, double exp, int k ){
		//double modulus = base.abs(), angle = base.phase();
		//return Complex.fromPolar( Math.pow( modulus, exp ), angle * exp );
		// Better
		double n = base.norm(), angle = base.phase() + Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( n, 0.5 * exp ), angle * exp );
	}
	public static Complex pow( Complex base, double exp ){
		return Complex.pow( base, exp, 0 );
	}
	public static Complex pow( double base, Complex exp, int k ){
		if( base < 0. ){
			double angle = Math.PI + Constant.M_2PI * k;
			return Complex.fromPolar( Math.pow( -base, exp.real ) * Math.exp( -angle * exp.imag ), angle * exp.real + exp.imag * Math.log( -base ) );
		}
		double angle = Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( base, exp.real ) * Math.exp( -angle * exp.imag ), angle * exp.real + exp.imag * Math.log( base ) );
	}
	public static Complex pow( double base, Complex exp ){
		return Complex.pow( base, exp, 0 );
	}
	public static Complex pow( double base, double exp, int k ){
		if( base < 0. ){
			double angle = Math.PI + Constant.M_2PI * k;
			return Complex.fromPolar( Math.pow( -base, exp ), angle * exp );
		}
		double angle = Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( base, exp ), angle * exp );
	}
	public static Complex pow( double base, double exp ){
		return Complex.pow( base, exp, 0 );
	}
	public static Complex sqrt( Complex a, int k ){
		//double modulus = base.abs(), angle = base.phase();
		//return Complex.fromPolar( Math.pow( modulus, 0.5 ), angle * 0.5 );
		// Better
		double n = a.norm(), angle = a.phase() + Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( n, 0.25 ), angle * 0.5 );
	}
	public static Complex sqrt( Complex a ){
		return Complex.sqrt( a, 0 );
	}
	public static Complex sqrt( double a, int k ){
		if( a < 0. ){
			double angle = Math.PI + Constant.M_2PI * k;
			return Complex.fromPolar( Math.sqrt( -a ), angle * 0.5 );
		}
		double angle = Constant.M_2PI * k;
		return Complex.fromPolar( Math.sqrt( a ), angle * 0.5 );
	}
	public static Complex sqrt( double a ){
		return Complex.sqrt( a, 0 );
	}
	public static void sqrt( Complex a, Complex[] res ){	// res = new Complex[2];
		double n = a.norm(), angle = a.phase();
		double modulus = Math.pow( n, 0.25 );
		res[0] = Complex.fromPolar( modulus, angle * 0.5 );
		res[1] = Complex.fromPolar( modulus, (angle + Constant.M_2PI) * 0.5 );
	}
	public static Complex[] sqrt( Complex a, boolean allSolutions ){
		Complex[] res = new Complex[2];
		Complex.sqrt( a, res );
		return res;
	}
	public static void sqrt( double a, Complex[] res ){	// res = new Complex[2];
		if( a < 0. ){
			double modulus = Math.sqrt( -a ), angle = Math.PI;
			res[0] = Complex.fromPolar( modulus, angle * 0.5 );
			res[1] = Complex.fromPolar( modulus, (angle + Constant.M_2PI) * 0.5 );
		}
		double modulus = Math.sqrt( a );
		res[0] = Complex.fromPolar( modulus, 0. );
		res[1] = Complex.fromPolar( modulus, Constant.M_2PI * 0.5 );
	}
	public static Complex[] sqrt( double a, boolean allSolutions ){
		Complex[] res = new Complex[2];
		Complex.sqrt( a, res );
		return res;
	}
	public static Complex root( Complex base, Complex exp, int k ){
		// Set base in polar form, and exp in cartesian
		double n = base.norm(), angle = base.phase() + Constant.M_2PI * k;
		Complex rec = exp.reciprocal();
		return Complex.fromPolar( Math.pow( n, 0.5 * rec.real ) * Math.exp( -angle * rec.imag ), angle * rec.real + 0.5 * rec.imag * Math.log( n ) );
	}
	public static Complex root( Complex base, Complex exp ){
		return Complex.root( base, exp, 0 );
	}
	public static Complex root( Complex base, int exp, int k ){
		double n = base.norm(), angle = base.phase() + Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( n, 0.5 / exp ), angle / exp );
	}
	public static Complex[] root( Complex base, int exp, boolean allSolutions ){
		double n = base.norm(), angle = base.phase();	// + 2PI * k
		double modulus = Math.pow( n, 0.5 / exp );
		// Get all the solutions
		int max = exp < 0 ? -exp : exp;
		Complex[] res = new Complex[ max ];
		for( int i = 0; i < max; i++ ){
			res[ i ] = Complex.fromPolar( modulus, (angle + Constant.M_2PI * i)/ exp );
		}
		return res;
	}
	public static Complex root( Complex base, double exp, int k ){
		double n = base.norm(), angle = base.phase() + Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( n, 0.5 / exp ), angle / exp );
	}
	public static Complex root( Complex base, double exp ){
		return Complex.root( base, exp, 0 );
	}
	public static Complex root( double base, Complex exp, int k ){
		Complex rec = exp.reciprocal();
		if( base < 0. ){
			double angle = Math.PI + Constant.M_2PI * k;
			return Complex.fromPolar( Math.pow( -base, rec.real ) * Math.exp( -angle * rec.imag ), angle * rec.real + rec.imag * Math.log( -base ) );
		}
		double angle = Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( base, rec.real ) * Math.exp( -angle * rec.imag ), angle * rec.real + rec.imag * Math.log( base ) );
	}
	public static Complex root( double base, Complex exp ){
		return Complex.root( base, exp, 0 );
	}
	public static Complex root( double base, double exp, int k ){
		if( base < 0. ){
			double angle = Math.PI + Constant.M_2PI * k;
			return Complex.fromPolar( Math.pow( -base, 1./exp ), angle / exp );
		}
		double angle = Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( base, 1./exp ), angle / exp );
	}
	public static Complex root( double base, double exp ){
		return Complex.root( base, exp, 0 );
	}
	public static Complex root( double base, int exp, int k ){
		if( base < 0. ){
			double angle = Math.PI + Constant.M_2PI * k;
			return Complex.fromPolar( Math.pow( -base, 1./exp ), angle / exp );
		}
		double angle = Constant.M_2PI * k;
		return Complex.fromPolar( Math.pow( base, 1./exp ), angle / exp );
	}
	public static Complex[] root( double base, int exp, boolean allSolutions ){
		double n = base, angle = 0.;	// + 2PI * k
		if( base < 0. ){
			n = -base;
			angle = Math.PI;
		}
		double modulus = Math.pow( n, 1./exp );
		// Get all the solutions
		int max = exp < 0 ? -exp : exp;
		Complex[] res = new Complex[ max ];
		for( int i = 0; i < max; i++ ){
			res[ i ] = Complex.fromPolar( modulus, (angle + Constant.M_2PI * i)/ exp );
		}
		return res;
	}

	// Logarithm and Exponential ---------------------------------------------------------------------------------------

	public static Complex log( Complex a, int k ) {
		//return new Complex( Math.log( a.abs() ), a.phase() );
		return new Complex( 0.5 * Math.log( a.norm() ), a.phase() + Constant.M_2PI * k );
	}
	public static Complex log( Complex a ){
		return Complex.log( a, 0 );
	}
	public static Complex log( double a, int k ) {
		//return Complex.log( Complex.fromDouble( a ) );
		if( a < 0. )	return new Complex( Math.log( -a ), Math.PI + Constant.M_2PI * k );
		return new Complex( Math.log( a ), Constant.M_2PI * k );
	}
	public static Complex log( double a ){
		return Complex.log( a, 0 );
	}
	public static Complex log( Complex a, Complex base, int ka, int kb ){
		return Complex.div( Complex.log( a, ka ), Complex.log( base, kb ) );
	}
	public static Complex log( Complex a, Complex base ){
		return Complex.log( a, base, 0, 0 );
	}
	public static Complex log( Complex a, double base, int ka, int kb ){
		return Complex.div( Complex.log( a, ka ), Complex.log( base, kb ) );
	}
	public static Complex log( Complex a, double base ){
		return Complex.log( a, base, 0, 0 );
	}
	public static Complex log( double a, Complex base, int ka, int kb ){
		return Complex.div( Complex.log( a, ka ), Complex.log( base, kb ) );
	}
	public static Complex log( double a, Complex base ){
		return Complex.log( a, base, 0, 0 );
	}
	public static Complex log( double a, double base, int ka, int kb ){
		return Complex.div( Complex.log( a, ka ), Complex.log( base, kb ) );
	}
	public static Complex log( double a, double base ){
		return Complex.log( a, base, 0, 0 );
	}
	public static Complex exp( Complex a ) {
		double expPart = Math.exp( a.real );
		return new Complex(expPart * Math.cos( a.imag ), expPart * Math.sin( a.imag ) );
	}
	public static Complex exp( double a ){
		return new Complex( Math.exp( a ), 0. );
	}
	/**
	 * https://en.wikipedia.org/wiki/Lambert_W_function
	 * z*(e^(z)) = e^(x) * (x*Cos(y) - y*Sin(y)) + i * e^(x) * (x*Sin(y) + y*Cos(y))
	 */
	public static Complex expLW( Complex a ) {
		double e = Math.exp( a.real ), s = Math.sin( a.imag ), c = Math.cos( a.imag );
		return new Complex( e * (a.real * c - a.imag * s), e * (a.real * s + a.imag * c) );
	}
	public static Complex expLW( double a ){
		return new Complex( a * Math.exp( a ), 0. );
	}

	// Trigonometric ---------------------------------------------------------------------------------------------------

	public static Complex sin( Complex a ) {
		return new Complex(Math.sin( a.real ) * Math.cosh( a.imag ), Math.cos( a.real ) * Math.sinh( a.imag ) );
	}
	public static Complex sin( double a ) {
		return new Complex( Math.sin( a ), 0. );
	}
	public static Complex cos( Complex a ) {
		return new Complex(Math.cos( a.real ) * Math.cosh( a.imag ), -Math.sin( a.real ) * Math.sinh( a.imag ) );
	}
	public static Complex cos( double a ) {
		return new Complex( Math.cos( a ), 0. );
	}
	public static void sinCos( Complex a, Complex[] res ){	// Complex[] res = new Complex[2];
		double sc = Math.sin( a.real ), cs = Math.cosh( a.imag ), cc = Math.cos( a.real ), ss = Math.sinh( a.imag );
		res[0] = new Complex(sc * cs, cc * ss );
		res[1] = new Complex(cc * cs, -sc * ss );
	}
	public static Complex[] sinCos( Complex a ){
		Complex[] res = new Complex[2];
		Complex.sinCos( a, res );
		return res;
	}
	// fsincos ???
	public static void sinCos( double a, Complex[] res ){	// Complex[] res = new Complex[2];
		double sc = Math.sin( a ), cc = Math.cos( a );
		res[0] = new Complex( sc, 0. );
		res[1] = new Complex( cc, 0. );
	}
	public static Complex[] sinCos( double a ){
		Complex[] res = new Complex[2];
		Complex.sinCos( a, res );
		return res;
	}
	public static Complex tan( Complex a ){
		//if( a.imag == 0. )	return new Complex( Math.tan( a.real ), 0. );
		//if( Precision.isZero( a.imag ) )	return new Complex( Math.tan( a.real ), 0. );
		double sc = Math.sin( a.real ), cs = Math.cosh( a.imag ), cc = Math.cos( a.real ), ss = Math.sinh( a.imag );
		Complex _sin = new Complex(sc * cs, cc * ss );
		Complex _cos = new Complex(cc * cs, -sc * ss );
		return Complex.div( _sin, _cos );
	}
	public static Complex tan( double a ){
		return new Complex( Math.tan( a ), 0. );
	}
	public static Complex asin( Complex a, int ks, int kl ){	// -I*ln( I*x + sqrt( 1 - x^2 ) )
		//Complex _1_xx = new Complex( 1. - a.real * a.real + a.imag * a.imag, -2. * a.real * a.imag );
		Complex _1_xx = new Complex( (1. - a.real)*(1. + a.real) + a.imag * a.imag, -2. * a.real * a.imag );
		Complex r = Complex.sqrt( _1_xx, ks );
		//Complex ix = new Complex( -a.imag, a.real );
		Complex ix = Complex.mulI( a );
		Complex l = Complex.log( Complex.add( ix, r ), kl );
		//return Complex.mul( _I, l );
		return Complex.mulI( l, -1. );
	}
	public static Complex asin( Complex a ){
		return Complex.asin( a, 0, 0 );
	}
	/*public static Complex asin( double a ){
		if( a < -1. || 1. < a )	return Complex.asin( Complex.fromDouble( a ) );
		return Complex.fromDouble( Math.asin( a ) );
	}*/
	public static Complex asin( double a, int ks, int kl ){
		return Complex.asin( Complex.fromDouble( a ), ks, kl );
	}
	public static Complex asin( double a ){
		return Complex.asin( Complex.fromDouble( a ), 0, 0 );
	}
	public static Complex acos( Complex a, int ks, int kl ){	// -I*ln( x + sqrt( x^2 - 1 ) )
		//Complex xx_1 = new Complex(  a.real * a.real - a.imag * a.imag - 1., 2. * a.real * a.imag );
		Complex xx_1 = new Complex(  (a.real - 1.) * (a.real + 1.) - a.imag * a.imag, 2. * a.real * a.imag );
		Complex r = Complex.sqrt( xx_1, ks );
		Complex l = Complex.log( Complex.add( a, r ), kl );
		//return Complex.mul( _I, l );
		return Complex.mulI( l, -1. );
	}
	public static Complex acos( Complex a ){
		return Complex.acos( a, 0, 0 );
	}
	/*public static Complex acos( double a ){
		if( a < -1. || 1. < a )	return Complex.acos( Complex.fromDouble( a ) );
		return Complex.fromDouble( Math.acos( a ) );
	}*/
	public static Complex acos( double a, int ks, int kl ){
		return Complex.acos( Complex.fromDouble( a ), ks, kl );
	}
	public static Complex acos( double a ){
		return Complex.acos( Complex.fromDouble( a ), 0, 0 );
	}
	public static Complex atan( Complex a, int k ){	// 0.5*I*ln( (I + x)/(I - x) )
		double den = a.real * a.real + (a.imag - 1.) * (a.imag - 1.);
		Complex d = new Complex( -(a.real * a.real + (a.imag - 1.)*(a.imag + 1.)) / den, -2. * a.real / den );
		Complex l = Complex.log( d, k );
		return Complex.mulI( l, 0.5 );
	}
	public static Complex atan( Complex a ){
		return Complex.atan( a, 0 );
	}
	public static Complex atan( double a, int k ){
		//return Complex.fromDouble( Math.atan( a ) );
		return Complex.atan( Complex.fromDouble( a ), k );
	}
	public static Complex atan( double a ){
		return Complex.atan( Complex.fromDouble( a ), 0 );
	}

	// Hyperbolic ------------------------------------------------------------------------------------------------------

	public static Complex sinh( Complex a ) {
		return new Complex(Math.sinh( a.real ) * Math.cos( a.imag ), Math.cosh( a.real ) * Math.sin( a.imag ) );
	}
	public static Complex sinh( double a ) {
		return new Complex( Math.sinh( a ), 0. );
	}
	public static Complex cosh( Complex a ) {
		return new Complex(Math.cosh( a.real ) * Math.cos( a.imag ), Math.sinh( a.real ) * Math.sin( a.imag ) );
	}
	public static Complex cosh( double a ) {
		return new Complex( Math.cosh( a ), 0. );
	}
	public static void sinhCosh( Complex a, Complex[] res ){	// Complex[] res = new Complex[2];
		double sc = Math.sinh( a.real ), cs = Math.cos( a.imag ), cc = Math.cosh( a.real ), ss = Math.sin( a.imag );
		res[0] = new Complex(sc * cs, cc * ss );
		res[1] = new Complex(cc * cs, sc * ss );
	}
	public static Complex[] sinhCosh( Complex a ){
		Complex[] res = new Complex[2];
		Complex.sinhCosh( a, res );
		return res;
	}
	public static void sinhCosh( double a, Complex[] res ){	// Complex[] res = new Complex[2];
		double sc = Math.sinh( a ), cc = Math.cosh( a );
		res[0] = new Complex( sc, 0. );
		res[1] = new Complex( cc, 0. );
	}
	public static Complex[] sinhCosh( double a ){
		Complex[] res = new Complex[2];
		Complex.sinhCosh( a, res );
		return res;
	}
	public static Complex tanh( Complex a ){
		//if( a.imag == 0. )	return new Complex( Math.tanh( a.real ), 0. );
		//if( Precision.isZero( a.imag ) )	return new Complex( Math.tanh( a.real ), 0. );
		double sc = Math.sinh( a.real ), cs = Math.cos( a.imag ), cc = Math.cosh( a.real ), ss = Math.sin( a.imag );
		Complex _sin = new Complex(sc * cs, cc * ss );
		Complex _cos = new Complex(cc * cs, sc * ss );
		return Complex.div( _sin, _cos );
	}
	public static Complex tanh( double a ){
		return new Complex( Math.tanh( a ), 0. );
	}
	public static Complex asinh( Complex a, int ks, int kl ){	// ln( x + sqrt( x^2 + 1 ) )
		Complex xxp1 = new Complex( 1. + a.real * a.real - a.imag * a.imag, 2. * a.real * a.imag );
		Complex r = Complex.sqrt( xxp1, ks );
		return Complex.log( Complex.add( a, r ), kl );
	}
	public static Complex asinh( Complex a ){
		return Complex.asinh( a, 0, 0 );
	}
	public static Complex asinh( double a, int ks, int kl ){
		//return new Complex( Math.log( a + Math.sqrt( a * a + 1. ) ), 0. );
		return Complex.asinh( Complex.fromDouble( a ), ks, kl );
	}
	public static Complex asinh( double a ){
		return Complex.asinh( Complex.fromDouble( a ), 0, 0 );
	}
	public static Complex acosh( Complex a, int ks, int kl ){	// ln( x + sqrt( x^2 - 1 ) )
		//Complex xx_1 = new Complex( a.real * a.real - 1. - a.imag * a.imag, 2. * a.real * a.imag );
		Complex xx_1 = new Complex( (a.real - 1.)*(a.real + 1.) - a.imag * a.imag, 2. * a.real * a.imag );
		Complex r = Complex.sqrt( xx_1, ks );
		return Complex.log( Complex.add( a, r ), kl );
	}
	public static Complex acosh( Complex a ){
		return Complex.acosh( a, 0, 0 );
	}
	/*public static Complex acosh( double a ){
		if( a < 1. )	return Complex.acosh( Complex.fromDouble( a ) );
		return new Complex( Math.log( a + Math.sqrt( (a - 1.)*(a + 1.) ) ), 0. );
	}*/
	public static Complex acosh( double a, int ks, int kl ){
		return Complex.acosh( Complex.fromDouble( a ), ks, kl );
	}
	public static Complex acosh( double a ){
		return Complex.acosh( Complex.fromDouble( a ), 0, 0 );
	}
	public static Complex atanh( Complex a, int k ){	// 0.5*ln( (1 + x)/(1 - x) )
		double den = (1. - a.real)*(1. - a.real) + a.imag * a.imag;
		Complex d = new Complex( ((1. - a.real)*(1. + a.real) - a.imag * a.imag) / den, 2. * a.imag / den );
		return Complex.mul( 0.5, Complex.log( d, k ) );
	}
	public static Complex atanh( Complex a ){
		return Complex.atanh( a, 0 );
	}
	/*public static Complex atanh( double a ){
		if( a <= -1. || 1. <= a )	return Complex.atanh( Complex.fromDouble( a ) );
		return new Complex( 0.5 * Math.log( (1. + a)/(1. - a) ), 0. );
	}*/
	public static Complex atanh( double a, int k ){
		return Complex.atanh( Complex.fromDouble( a ), k );
	}
	public static Complex atanh( double a ){
		return Complex.atanh( Complex.fromDouble( a ), 0 );
	}

	// Other useful functions ------------------------------------------------------------------------------------------

	public static double distance( Complex a, Complex b ){
		double difR = a.real - b.real, difI = a.imag - b.imag;
		return Math.sqrt( difR * difR + difI * difI );
	}
	public static double distance2( Complex a, Complex b ){
		double difR = a.real - b.real, difI = a.imag - b.imag;
		return ( difR * difR + difI * difI );
	}
	
	// https://en.wikipedia.org/wiki/Lambert_W_function
	// TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//public static Complex lambertW( Complex a )		Numerical evaluation
	//public static Complex lambertW( double a )
	// https://github.com/DarkoVeberic/LambertW
	// https://ingalidakis.com/math/LWCalculating2.html


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	public boolean equals( Object x ) {
		if( x == null || this.getClass() != x.getClass() )	return false;
		if( this == x )	return true;
		Complex that = (Complex) x;
		//return (this.real == that.real) && (this.imag == that.imag);
		return Precision.equals( this.real, that.real ) && Precision.equals( this.imag, that.imag );
	}
	@Override
	public int compareTo( Complex x ){		// Sort by real part
		if( x == null )	return 1;
		//double diff = this.real - x.real;
		//return diff < 0. ? -1 : (diff > 0. ? 1 : 0);
		int res = Double.compare( this.real, x.real );
		if( res == 0 )	return Double.compare( this.imag, x.imag );
		return res;
	}
}
