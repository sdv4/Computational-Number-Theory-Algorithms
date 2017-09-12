import java.lang.Math;
import java.util.*;
import java.math.BigInteger;
import java.lang.IllegalArgumentException;

/*
 * @author  Shane Sims
 * @version  23 July 2017
 */

public class TonelliSqrtModP{
	
	/* This method computes the square root of an integer
	 * modulo a prime moduli, based on the pseduo-code for A. Tonelli's
	 * algorithm for this purpose. Pseudo-code found 
	 * in Crandall and Pomerance's Prime Numbers: a Computational 
	 * perspective.
	 * @param a  quadratic residue of p (where Lagendre symbol (a/p)=1), whoes sqrt will be calculated
	 * @param p  an odd prime
	 */
	// TODO: need to check if a has a square root mod p by checking Lagendre symbol. 
	// Currently we assume taht 
	public static int SqrtModP(int a, int p) {
		int x;					//sqrt(a) to be returned.

		//TODO: convert following series of if statements to if-else if- else statements

	/* Check cases p \equiv 3, 5 or 7 (mod 8) */
		a = a % p;					
		if(p % 8 == 3 || p % 8 == 7){
			x = (int)Math.pow(a, ((p+1)/4)) % p;
			return x;
		}
		if(p % 8 == 5) {
			x = (int)Math.pow(a, ((p+3)/8)) % p;
			int c = (x*x) % p;
			if (c != (a % p))
				x = x * (int)Math.pow(2, ((p-1)/4)) % p;
			return x;
		}

	/* Then we have the case: p \equiv 1 (mod 8) */
		
		/* Find a d such that d is a Quadratic Non-Residue of p*/
		int d;						//will hold value of QNR of p
		d = 2;						//initialised to 2 to begin searching
		while(DetermineQuadRes.LegJacSym(d, p) != -1)
			d++;					//inc d until QNR of p
	
		/* Find values s, t such that p - 1 = (2^s)t*/
		int pmo = p - 1;				//pmo will represent P Minus One (p-1)
		int s = 0;
		int t;
		while(pmo % 2 == 0){				//extract factors of 2 from pmo
			s++;
			pmo /= 2;
		}
		t = (p-1) / ((int)Math.pow(2, s));
		int A = (int)Math.pow(a, t) % p;
		int D = (int)Math.pow(d, t) % p;
		int m = 0;					//m is 2\mu from pseudocode
		for(int i = 0; i < s; i++) {
			if((int)Math.pow(A*(Math.pow(D, m)), Math.pow(2, s-1-i)) % p == p-1)
				m = m + (int)Math.pow(2, i);
		}
		x = (int)((Math.pow(a,((t+1)/2)) * Math.pow(D, (m/2))) % p);
		return x;
	}//end method

	//Big Integer a and int p version of above method
	public static int SqrtModP(BigInteger a, int p) {

		BigInteger bigp = BigInteger.valueOf((long) p);

		int quadres = DetermineQuadRes.LegJacSym(a,bigp);
		if(DetermineQuadRes.LegJacSym(a, bigp) != 1){
			throw new IllegalArgumentException();
		}
		int x;					//sqrt(a) to be returned.
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger three = BigInteger.valueOf(3L);
		BigInteger four = BigInteger.valueOf(4L);
		BigInteger five = BigInteger.valueOf(5L);
		BigInteger seven = BigInteger.valueOf(7L);
		BigInteger eight = BigInteger.valueOf(8L);

		
		/* Check cases p \equiv 3, 7 or 5 (mod 8) */
		a = a.mod(bigp);
		if((bigp.mod(eight)).compareTo(three) == 0 || (bigp.mod(eight)).compareTo(seven) == 0){
			BigInteger exponent = (bigp.add(BigInteger.ONE)).divide(four);
			BigInteger bigx = (a.modPow(exponent, bigp));
			return bigx.intValue();

		}
		else if((bigp.mod(eight)).compareTo(five) == 0){
			BigInteger exponent = (bigp.add(three)).divide(eight);
			BigInteger bigx = (a.modPow(exponent, bigp));
			BigInteger bigc = (bigx.pow(2)).mod(bigp);
			if (bigc.compareTo(a) != 0){
				exponent = (bigp.subtract(BigInteger.ONE)).divide(four);
				bigx = (bigx.multiply(two.modPow(exponent, bigp))).mod(bigp);
			}
			return bigx.intValue();

		}
		/* Then we have the case: p \equiv 1 (mod 8) */
		else{


			/* Find a d such that d is a Quadratic Non-Residue of p*/
			BigInteger d = two;						//will hold value of QNR of p. Initialised to 2 to begin searching
			while(DetermineQuadRes.LegJacSym(d, bigp) != (-1)){
				d = d.add(BigInteger.ONE);
			}

			/* Find values s, t such that p - 1 = (2^s)t*/
			BigInteger pmo = bigp.subtract(BigInteger.ONE);
			int s = 0;
			BigInteger t;
			while((pmo.remainder(two)).compareTo(BigInteger.ZERO) == 0){
				s++;
				pmo = pmo.divide(two);
			}
			t = (bigp.subtract(BigInteger.ONE)).divide(BigInteger.valueOf((long)Math.pow(2, s)));

			BigInteger bigA = a.modPow(t, bigp);
			BigInteger bigD = d.modPow(t, bigp);
			int m = 0;
			for(int i = 0; i < s; i++){
				if((((bigA.multiply(bigD.pow(m))).pow(two.pow(s-1-i).intValue())).mod(bigp)).compareTo(bigp.subtract(BigInteger.ONE)) == 0)
					m = m + (int)Math.pow(2, i);
			}
			BigInteger a1 = (a.modPow((t.add(BigInteger.ONE)).divide(two), bigp));
			BigInteger a2 = bigD.pow(m/2);
			BigInteger bigx = (a1.multiply(a2)).mod(bigp);
			return bigx.intValue();

		}//end else
	}//end method


	//Big Integer a and BigInteger p version of above method
	public static BigInteger SqrtModP(BigInteger a, BigInteger bigp) {


		int quadres = DetermineQuadRes.LegJacSym(a,bigp);
		if(DetermineQuadRes.LegJacSym(a, bigp) != 1){
			throw new IllegalArgumentException();
		}
		int x;					//sqrt(a) to be returned.
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger three = BigInteger.valueOf(3L);
		BigInteger four = BigInteger.valueOf(4L);
		BigInteger five = BigInteger.valueOf(5L);
		BigInteger seven = BigInteger.valueOf(7L);
		BigInteger eight = BigInteger.valueOf(8L);

		
		/* Check cases p \equiv 3, 7 or 5 (mod 8) */
		a = a.mod(bigp);
		if((bigp.mod(eight)).compareTo(three) == 0 || (bigp.mod(eight)).compareTo(seven) == 0){
			BigInteger exponent = (bigp.add(BigInteger.ONE)).divide(four);
			BigInteger bigx = (a.modPow(exponent, bigp));
			return bigx;

		}
		else if((bigp.mod(eight)).compareTo(five) == 0){
			BigInteger exponent = (bigp.add(three)).divide(eight);
			BigInteger bigx = (a.modPow(exponent, bigp));
			BigInteger bigc = (bigx.pow(2)).mod(bigp);
			if (bigc.compareTo(a) != 0){
				exponent = (bigp.subtract(BigInteger.ONE)).divide(four);
				bigx = (bigx.multiply(two.modPow(exponent, bigp))).mod(bigp);
			}
			return bigx;

		}
		/* Then we have the case: p \equiv 1 (mod 8) */
		else{


			/* Find a d such that d is a Quadratic Non-Residue of p*/
			BigInteger d = two;						//will hold value of QNR of p. Initialised to 2 to begin searching
			while(DetermineQuadRes.LegJacSym(d, bigp) != (-1)){
				d = d.add(BigInteger.ONE);
			}

			/* Find values s, t such that p - 1 = (2^s)t*/
			BigInteger pmo = bigp.subtract(BigInteger.ONE);
			int s = 0;
			BigInteger t;
			while((pmo.remainder(two)).compareTo(BigInteger.ZERO) == 0){
				s++;
				pmo = pmo.divide(two);
			}
			t = (bigp.subtract(BigInteger.ONE)).divide(BigInteger.valueOf((long)Math.pow(2, s)));

			BigInteger bigA = a.modPow(t, bigp);
			BigInteger bigD = d.modPow(t, bigp);
			int m = 0;
			for(int i = 0; i < s; i++){
				if((((bigA.multiply(bigD.pow(m))).pow(two.pow(s-1-i).intValue())).mod(bigp)).compareTo(bigp.subtract(BigInteger.ONE)) == 0)
					m = m + (int)Math.pow(2, i);
			}
			BigInteger a1 = (a.modPow((t.add(BigInteger.ONE)).divide(two), bigp));
			BigInteger a2 = bigD.pow(m/2);
			BigInteger bigx = (a1.multiply(a2)).mod(bigp);
			return bigx;

		}//end else
	}//end method




	public static void main(String[] args) {
		Scanner reader = new Scanner(System.in);
		System.out.println("Enter odd prime: ");
		int p = reader.nextInt();
		System.out.println("Enter quadratic residue (mod p) whose square root will be calculated: ");
		BigInteger a = BigInteger.valueOf((long) reader.nextInt());
		System.out.println();
		System.out.print("Square root of " + a + " (mod " + p + ") is : ");
		System.out.println(SqrtModP(a, p));
	}


}//end class
