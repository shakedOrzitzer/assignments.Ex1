package Ex1;

import static java.lang.Double.NaN;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	public static final double[] MINUS_ONE = {-1};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p,x1);
        double x12 = (x1+x2)/2;
        double f12 = f(p,x12);
        if (Math.abs(f12)<eps) {return x12;}
        if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
        else {return root_rec(p, x12, x2, eps);}
    }
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            if (lx==3&&ly==3){
                double x1=xx[0],y1=yy[0],x2=xx[1],y2=yy[1],x3=xx[2],y3=yy[2];
                double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                double A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                double B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
                double C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                double xv = -B / (2*A);
                double yv = C - B*B / (4*A);
                ans = new double[]{C,B,A};}
            else {
                if(lx==2&&ly==2){
                    double x1=xx[0],y1=yy[0],x2=xx[1],y2=yy[1];
                    double leny= Math.abs(y1-y2);
                    double lenx= Math.abs(x1-x2);
                    double m=leny/lenx;
                    double b=y1-m-x1;
                    ans = new double[]{b,m};
                    ans = new double[]{b,m};
                }
            }
		}
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true if p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        int bigPow =0;
        if(p1.length>=p2.length) {
            bigPow=p1.length;}
        else {bigPow=p2.length;}
        for(int i=0;i<=bigPow;i++) {
            if(Math.abs((f(p1,i)-f(p2,i))) > EPS) {return false;}
        }
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
     *
     * int len=poly.length
     * for(int i=len;i>=0;i-=1){
     *      ans+= string(poly[i])+"X^"+string(i)
     * }
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            poly=compact(poly);
            int len=poly.length;
            for(int i=len-1;i>=0;i-=1){
                if (i==len-1) {
                    ans += (poly[i]);
                }
                else{
                        if (poly[i] == 0) {
                        }
                        else {
                            if (i==0){
                                if(poly[i]<0) {ans += " " + (poly[i]);}
                                if(poly[i]>0) {ans += " +" + (poly[i]);}
                            }
                            else {
                                if (poly[i] < 0) {
                                    ans += " " + (poly[i]) + "x^" + (i);
                                }
                                if (poly[i] > 0) {
                                    ans += " +" + (poly[i]) + "x^" + (i);
                                }
                            }
                        }
                    }
                }
            }
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = -1;
        for(double i=x1;i<=x2;i+=0.00001) {
            if(Math.abs((f(p1,i)-f(p2,i))) < eps) {return i;}
        }
		return x1;
	}

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param poly - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] poly, double x1, double x2, int numberOfSegments) {
		double ans = 0;
        double p= (x2-x1) / (numberOfSegments+1);
        double x3=0,x4=0, f3=0,f4=0,lenX=0,lenY=0;
        for (double i=x1;i<=x2-p;i+=p) {
            x3=i;
            x4=i+p;
            f3=f(poly,x3);
            f4=f(poly,x4);
            lenX=Math.pow((x4-x3),2);
            lenY=Math.pow((f4-f3),2);
            ans+=Math.sqrt(lenX+lenY);
        }
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        double len= (x2-x1)/ numberOfTrapezoid;
        double base1=0,base2=0,temp=0,xi=0,nextX=0,trapezoidArea=0;
        for(double i=x1;i<=x2-len;i+=len) {
            xi=i;
            nextX=xi+len;
            base1 = Math.abs(f(p1,xi) - f(p2,xi));
            base2 = Math.abs(f(p1,nextX) - f(p2,nextX));
            double p= sameValue(p1,p2,i,i+len,EPS);
            if(i==p){
                trapezoidArea=((Math.abs(base1)+Math.abs(base2)) * len)/2;
                ans+=trapezoidArea;
            }
            else{
                double tri1=(Math.abs(base1* (p-i)))/2;
                double tri2=(Math.abs(base2* (p+i)))/2;
                ans+= tri1+tri2;
            }
        }
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0, "-1.2x^3 +3.1x^2 +2.0"
        String[] p2= p.split(" ");
        int len=p2.length;
        int bigPow= getBFromMonom(p2[0]);
        String[] revP2= new String[bigPow+1];
        for(int i=p2.length-1;i>=0;i--){
            if(getBFromMonom(p2[i]) ==i){}
                revP2[i]=p2[i];
        }
        double[] poly = new double[revP2.length];
        for (int i=revP2.length-1;i>=0;i-=1) {
            if(getBFromMonom(revP2[i])==i){
                poly[i]=getAFromMonom(revP2[i]);
            }
            int run=p2.length-1;
            for(int i=0;i<revP2.length;i+=1){
                if( getBFromMonom(p2[run])==i){
                    poly[i]=getAFromMonom(p2[run]);
                }
                else{poly[i]=0}
            }
        }
		return poly;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     * if l2>l1
     *      replace l1,l2
     * for (int i=0;i<l2;i++){
     *      p1[i]+=p2[i]
     * }
     *
	 * @param p1 longer
	 * @param p2 shorter
	 * @return p2
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        if (p1==null)
            return p2;
        if (p2==null)
            return p1;
        p1=compact(p1);
        p2=compact(p2);
        int l1=p1.length;
        int l2=p2.length;
        int maxl=Math.max(l1,l2);
        int minl=Math.min(l1,l2);
        double [] res=new double[maxl];
        for(int i=0;i<minl;i++){
            res[i]=p1[i]+p2[i];
        }

        if (p2.length>p1.length){
            for(int i=minl;i<maxl;i++){
                res[i]=p2[i];
            }
        }
        else{
            if(p1.length>p2.length){
            for(int i=minl;i<maxl;i++){
                res[i]=p1[i];}
            }
        }
		return compact(res);
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     * p1,p2 compact
     * double[] res= new double[]
     * for (int i=0;i<p1.len;i++){
     *      double[] c= mul(p2,p1[i],i]
     *      and=add (ans,c)
     *
     *     }
     * }
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        compact(p1);
        compact(p2);
        double[] res= new double[p1.length+p2.length];
        for (int i=0;i<p1.length-1;i++){
            double[] c= mul1(p2,p1[i],i);
            res=add (res,c);}
        compact(res);
		return res;

    }
	/**
	 * This function computes the derivative of the p0 polynomial function.
     * input (poly)
     * output ans= ZERO
     * if ( poly != null && poly.length>1)
     * int len = poly.length
     * ans= new double [len-1]
     * for(int i=0;i<ans.length
     *  ans[i]= poly[i+1] * (i+1)
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;
        if ( po!= null && po.length>1) {
            int len = po.length;
            ans= new double [len-1] ;
        for(int i=0;i<ans.length;i++)
            ans[i]= po[i+1] * (i+1);
        }
		return ans;
	}
    /**
     * פונקציית עזר לmul
     * function gets poly, m= (a*x**b)
     * returns poly*m
     * poly= {1,2,3} m=2*x^3
     *
     * mul (poly,a,b)
     * res= new double[poly.len+b]
     * for(int i=0;i<res.len;i++)
     *      res[i]=0
     * for(int i=b;b<res.len;i++)
     *      res[i]=a*poly[i]
     * @param poly
     * @param b
     * @return
     */
    public static double[] mul1(double[] poly,double a, int b) {
        double[] res1 = new double[poly.length+b];
        for(int i = b; i< res1.length; i+=1)
            res1[i]=a*poly[i-b];
        res1 =compact(res1);
        return res1;
    }
    /**
     *this function gets a polynome, and returns polynome with no 0's in the biggest pow
     * for {1,2,4,3,0,4,0,0}
     * returns {1,2,4,3,0,4}
     *
     * run on poly and find biggest index for which poly[i] !=0
     * double[] res=new double[found index]
     * for (int i=0;i<res.length;i++)
     *      res[i]=poly[i]
     * return res
     * @param p1
     * @return
     */
    public static double[] compact(double[] p1) {
        int max=0;
        for (int i=0;i<p1.length;i++) {
            if (p1[i]!=0) {max=i;}
        }
        double[] res=new double[max+1];
        for (int i=0;i<res.length;i++) {res[i]=p1[i];}
        return res;
    }
    /**
     * function checks if input in a number
     */
    public static Double stringToNumber (String s) {
        int ans= -1;
        try{
            return parseDouble(s);}
        catch (NumberFormatException e) {
            return -1.0;
        }
    }

    /**
     * function gets a monom a*x^b and returns a
     * @param monom
     * @return
     */
    public static double getAFromMonom (String monom) {
        double ans= -1;
        if (stringToNumber(monom)>=0) return stringToNumber(monom);
        String[] n= monom.split("x");
        return stringToNumber(n[0]);

        }

    /**
     * this function gets a monom a*x^b and returns b
     */
    public static int getBFromMonom (String monom) {
        monom= monom.substring(1);
        double b= stringToNumber(monom);
        if (b>=0) return 0;
        else { // monom ax^b
            String[]n=monom.split("\\^");
                if(stringToNumber(n[1])!=-1){return parseInt(n[1]);}

        }

        return parseInt(monom);
    }

    /**
     * this function computes the difference of 2 polynomes
     */
    public static double[] diff(double[] p1, double[] p2) {
        if (p1==null)
            return p2;
        if (p2==null)
            return p1;
        p1=compact(p1);
        p2=compact(p2);
        int l1=p1.length;
        int l2=p2.length;
        double[]ans;
        if (l2>l1){
            ans=new double[l2];
            for (int i = 0; i <l1 ; i++) {
                ans[i] = p1[i]-p2[i];
            }
            for (int i = l1; i < l2 ; i++) {
                ans[i]=p2[i];
            }
        }
        else {
            ans=new double[l1];
            for (int i = 0; i <l2 ; i++) {
                ans[i] = p1[i]-p2[i];
            }
            for (int i = l2; i < l1 ; i++) {
                ans[i]=p1[i];
            }
        }
        return ans;
    }
}
