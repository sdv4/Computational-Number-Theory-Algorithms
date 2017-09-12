import java.util.Scanner;
import java.math.BigInteger;

/* 
 * @author Shane Sims
 * @version  24 July 2017
*/
public class DetermineQuadRes {
	

	/*
 	 * This algorithm determines the Jacobi symbol (a/m). When
	 * m is an prime number, this is identical to the Legendre
	 * symbol.
	 * Based on the pseudo-code for this purpose in the Crandall 
 	 * & Pomerance text.
	 * @param a  is an integer for which this algorithm will determine
	 * its quadratic residuocity modulo m.
	 * @param m is an odd prime moduli.
	 * @return t is -1 if a is not a quadratic residue mod m, 1 if it is
	 * and 0 if a is congruent to 0 mod m
	 */
	public static int LegJacSym(int a, int m) {
		a = a%m;
		int t = 1;
		while(a != 0) {
			while(a%2 == 0) {
				a = a/2;
				if(m%8 == 3 || m%8 == 5)
					t = t * -1;
			}
			int temp = a;
			a = m;
			m = temp;				//(a,m) = (m,a)
			if(a % 4 == 3 && m % 4 == 3)
				t = t * -1;
			a = a%m;
		}
		if (m == 1)
			return t;
		return 0;
	}


	// BigInteger a and m version of the above method
	public static int LegJacSym(BigInteger biga, BigInteger bigm) {
		BigInteger a = biga;
		BigInteger m = bigm;
		a = a.mod(m);

		int t = 1;
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger three = BigInteger.valueOf(3L);
		BigInteger four = BigInteger.valueOf(4L);		
		BigInteger five = BigInteger.valueOf(5L);		
		BigInteger eight = BigInteger.valueOf(8L);		

		while(a.compareTo(BigInteger.ZERO) != 0){
			while((a.mod(two)).compareTo(BigInteger.ZERO) == 0){
				a = a.shiftRight(1);
				if((m.mod(eight)).compareTo(three) == 0 || (m.mod(eight)).compareTo(five) == 0)
					t = t * -1;
			}
			BigInteger bigtemp = a;
			a = m;
			m = bigtemp;
			if((a.mod(four)).compareTo(three) == 0 && (m.mod(four)).compareTo(three) == 0)
				t = t * (-1);
			a = a.mod(m);		
		}
		if(m.compareTo(BigInteger.ONE) == 0)
			return t;
		return 0;
	}// end imethod



	// BigInteger version of the above method
	public static int LegJacSym(BigInteger a, int m) {
		BigInteger bigm = BigInteger.valueOf((long) m);
		a = a.remainder(bigm);
		int t = 1;
		BigInteger zero = BigInteger.valueOf(0L);
		BigInteger one = BigInteger.valueOf(1L);
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger three = BigInteger.valueOf(3L);
		BigInteger four = BigInteger.valueOf(4L);		
		BigInteger five = BigInteger.valueOf(5L);		
		BigInteger eight = BigInteger.valueOf(8L);		
		


		while(a.compareTo(zero) != 0){
			while((a.remainder(two)).compareTo(zero) == 0){
				a = a.shiftRight(1);
				if((bigm.remainder(eight)).compareTo(three) == 0 || (bigm.remainder(eight)).compareTo(five) == 0)
					t = t * (-1);
			}
			BigInteger bigtemp = a;
			a = bigm;
			bigm = bigtemp;
			if((a.remainder(four)).compareTo(three) == 0 && (bigm.remainder(four)).compareTo(three) == 0)
				t = t * -1;
			a = a.remainder(bigm);		
		}
		if(bigm.compareTo(one) == 0)
			return t;
		return 0;
	}// end imethod





	/*
	* Main to allow calculation of Legendre/Jacobi symbols.
	* NOTE: no exception handelling provided. Use for testing purposes only.
	*/
	public static void main(String[] args) {
		System.out.println("To calculate the Legendre/Jacobi symbol (a/m), ");
		System.out.print("enter a: ");
		Scanner sc = new Scanner(System.in);
     		int a = sc.nextInt();
		System.out.print("enter m: ");
		int m = sc.nextInt();
		System.out.println("You entered:");
		System.out.println("a = " + a);
		System.out.println("m = " + m);
		System.out.println();
		System.out.print("(" + a + "/" + m + ") = ");
		System.out.println(LegJacSym(a, m)); 
		System.out.println();
	}





}
