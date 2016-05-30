import java.util.Scanner;

/* 
 * @author Shane Sims
 * @version  30 May 2016
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
	 * @param m is an odd integer moduli.
	 * @return is -1 if a is not a quadratic residue mod m, 1 if it is
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
			if(a % 4 == m && m % 4 == 3)
				t = t * -1;
			a = a%m;
		}
		if (m == 1)
			return t;
		return 0;
	}

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
