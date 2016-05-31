import java.lang.Math;

/*
 * @author  Shane Sims
 * @version  31 May 2016
 */

public class TonelliSqrtModP{
	
	/* This method computes the square root of an integer
	 * modulo a prime moduli, based on the pseduo-code for A. Tonelli's
	 * algorithm for this purpose. Pseudo-code found 
	 * in Crandall and Pomerance's Prime Numbers: a Computational 
	 * perspective.
	 * @int a  a quadratic residue of p, whoes sqrt will be calculated
	 * @ int p  an odd prime
	 */
	public static int SqrtModP(int a, int p) {
		int x;					//sqrt(a) to be returned.


	/* Check cases p \equiv 3, 5 or 7 (mod 8) */
		a = a % p;					
		if(p % 8 == 3 || p % 8 == 7){
			x = (int)Math.pow(a, ((p+1)/4)) % p;
			return x;
		}
		if(p % 8 == 5) {
			x = (int)Math.pow(a, ((p+3)/8)) % p;
			int c = (x*x) % p;
			if (c % p != a)
				x = x * (int)Math.pow(2, ((p-1)/4)) % p;
			return x;
		}

	/* Then we have the case: p \equiv 1 (mod 8) */
		
		/* Find a d such that d is a Quadratic Non-Residue of p*/
		int d;						//will hold value of QNR of p
		d = 2;						//initialised to 2 to begin searching
		while(DetermineQuadRes.LegJacSym(d, p) != -1)
			d++;					//inc d until QNR of p
		System.out.println("d: " +d);
		System.out.println(DetermineQuadRes.LegJacSym(d, p));	
		/* Find values s, t such that p - 1 = (2^s)t*/
		int pmo = p - 1;				//pmo will represent P Minus One (p-1)
		int s = 0;
		int t;
		while(pmo % 2 == 0){				//extract factors of 2 from pmo
			s++;
			pmo /= 2;
		}
		System.out.println("s: " + s);
		t = (p-1) / ((int)Math.pow(2, s));
		System.out.println("t: " + t);	
		int A = (int)Math.pow(a, t) % p;
		int D = (int)Math.pow(d, t) % p;
		int m = 0;					//m is 2\mu from pseudocode
		for(int i = 0; i < s; i++) {
			if((int)Math.pow(A*(Math.pow(D, m)), Math.pow(2, s-1-i)) % p == p-1)
				m = m + (int)Math.pow(2, i);
		}
		System.out.println("AD^m = " + (A * Math.pow(D, m) % p)   );
		x = (int)((Math.pow(a,((t+1)/2)) * Math.pow(D, (m/2))) % p);
		return x;
	}//end method

	public static void main(String[] args) {
		System.out.println(SqrtModP(7, 113));  //should be 32
		System.out.println(SqrtModP(2, 41));  //should be 17
		System.out.println(SqrtModP(2, 113));//should be 62
	}


}//end class
