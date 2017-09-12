import java.lang.Math.*;
import java.util.*;
import java.math.BigInteger;
import java.io.*;
import java.util.Scanner;

/** This is an implementation of the Sieving portion of the Multiple Polynomial 
* Quadratic Sieve factoring algorithm. Built on code written in BasicQS and based 
* on the pseduocode provided in 
* the thesis 'Factoring Integers with the Self-Initializing Quadratic Sieve' by
* Scott P. Contini was consulted in the design of this code. The notation matches
* that found in Contini's thesis as far as possible so as to increase the didactic 
* value of this class. 
*
* This version is for numbers N of approximately 60 digits.
*
* TODO: Split sieving array into blocks for sieving
*
* @author Shane Sims <shane.sims.ss@gmail.com>
* @version 31 July 2017 
*
*/
public class MPQSsixty {	

	public static final int ERROR_TERM = 9;
	public static final int M = 600000;
	public static final int DESIRED_NUM_RELATIONS = 10;



	/**
	* This method determine the smooth relations required to factor N.
	* 
	* @param N	The integer to be factored. Must be an odd composite
	* free of prime powers. For this versions N is approx 60 digits in length.   
	* @param F	The smoothness bound to which the returned relations will
	* be bound by.
	* @return 	An array with each element containing an array with length 2, 
	* representing an ordered pair of the form (x, (ax+b)^2 - N), where the value 
	* x leads (ax + b)^2-n to be B-smooth.
	*
	*/
	public static BigInteger[][] getSmoothRelations(BigInteger N, int F){

		BigInteger[][] S;								// Array that will hold the smooth relations when returned.
		ArrayList<Integer> factorBase;							// Primes which will be used in the sieving stage. 
		int[] modSqrts;									// Will hold t in t^2 \equiv N (mod p) for each p in Factor Base
		int[] solution1;								// Will hold a^-1(t - b) \equiv N (mod p) for each prime in factor base
		int[] solution2;								// Will hold a^-1(-t - b) \equiv N (mod p) for each prime in factor base
		int[] logp;									// Approximate logarithm for each prime in the factor base
		int K;										// Method will return when K+1 smooth relations have been found.
		int offset = 0;									// Used if/when sieving stage returned to
		BigInteger a;									// Term 'a' in g_a,b(x) = (ax+b)^2 - N 
		BigInteger b;									// Term 'b' in g_a,b(x)
		BigInteger c;									// Term 'c' in g_a,b(x) = (ax+b)^2 - N = a(ax^2 + 2bx + c)
	
	//TODO: check preconditions of N are met here
	
	/******************** Compute Startup Data **********************************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/


		System.out.println("Computing Startup Data...");


		BigInteger squareRootN = BigIntegerSqrt.bigIntegerSqrtCeiling(N);		
		factorBase = determineFactorBase(F, N);						// Initialize factorBase.
		K = factorBase.size();								// Save factor base size. Need K+1 smooth valued poly to ensure factorization
		S = new BigInteger[K+1][4];							// Elements are pairs (x,(x+b)^2-N)				
		logp = new int[K];								
		modSqrts = new int[K];								// Initialize modSqrts to have space for each element of factorBase
		modSqrts[0] = 1;								// TonelliSqrtModP requires odd prime, so only even prime hard coded
		solution1 = new int[K];
		solution2 = new int[K];

		for(int i = 1; i < modSqrts.length; i++) {					// Solve t in t^2 \equiv N (mod p) for each p in factorBase
			modSqrts[i] = TonelliSqrtModP.SqrtModP(N, factorBase.get(i));
		}

		for(int i = 0; i < logp.length; i++) {
			logp[i] = (int) Math.round(Math.log(factorBase.get(i))/Math.log(2));	// Get approximate log_2(p) using change of base for each factor base prime
		}

		System.out.println("Compute Startup Data complete...");

	/******************* Initialization Stage ***********************************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/


		System.out.println("Starting Initialization Stage...");

		int smoothRelFound = 0;								// Tracks the number of smooth relations verified by trial division.
		BigInteger bigq = BigInteger.valueOf(0L);
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger bigM = BigInteger.valueOf((long) M);
		
//TODO: fix how q is chosen. Below method is not optimal and result is 'a' that is too big


	while(smoothRelFound < DESIRED_NUM_RELATIONS){

		/* Choose paramaters a,b,c */

		// Find q and set a
		if(bigq.compareTo(BigInteger.ZERO) == 0){
			bigq = N.multiply(two);						
			bigq = BigIntegerSqrt.bigIntegerSqrtCeiling(bigq);
			bigq = bigq.divide(bigM);
			bigq = BigIntegerSqrt.bigIntegerSqrtCeiling(bigq);
			bigq = bigq.subtract(BigInteger.ONE);	
		}
		else
			bigq = bigq.add(BigInteger.ONE);

		do {
			bigq = bigq.nextProbablePrime();
		} while ((DetermineQuadRes.LegJacSym(N, bigq) != 1));

		a = bigq.multiply(bigq);							// Set a = q^2
	
		// Compute b	
		BigInteger bigbPrime = TonelliSqrtModP.SqrtModP(N, bigq);			// Get b'=b in b^2 \equiv N (mod q)
		b = henselLift(bigbPrime, bigq, N);						// Set b for this polynomial
		c = ((b.pow(2)).subtract(N)).divide(a);						// Set c for this polynomial
	
	
		// Compute soln1_p and soln2_p for each p in factor base (the roots of g_a,b(x) mod p
		
		BigInteger sol1, sol2;								// Temp. variables to hold BigInteger roots before storing as int
		for(int i = 0; i < solution1.length; i++){
			if(bigq.compareTo(BigInteger.valueOf((long) factorBase.get(i))) != 0){
				BigInteger bigp = new BigInteger(Integer.toString(factorBase.get(i)));
				BigInteger bigaModInv = a.modInverse(bigp);
				sol1 = (bigaModInv.multiply((BigInteger.valueOf((long) modSqrts[i])).subtract(b))).mod(BigInteger.valueOf((long) factorBase.get(i)));
				sol2 = (bigaModInv.multiply((BigInteger.valueOf((long) modSqrts[i]*(-1L))).subtract(b))).mod(BigInteger.valueOf((long) factorBase.get(i)));

				solution1[i] = sol1.intValue();
				solution2[i] = sol2.intValue();
			}
			else{
				solution1[i] = -1;
				solution2[i] = -1;						// Mark unusable solutions.
			}
		}

		System.out.println("Initialization Stage complete...");

	/**************************** Sieving Stage *********************************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/
	
		System.out.println("Starting Sieving Stage...");

		int index1 = 0;									// Will hold latest value of (solution1_p + ip) or (solution2_p + ip)
		System.out.println("Smooth relations found: " + smoothRelFound);
		int[] sieveArray = new int[2*M + 1];						// Initialize sieve array with 0's.
		for(int j = 0; j < factorBase.size(); j++) {					// For each prime in the factor base.
		int k = 0;	
			int p = factorBase.get(j);
			int s1 = solution1[j];
			int s2 = solution2[j];

			if(s1 != -1){
				int i = -1;
				//find lower bound value for i such that soln1 +i*p >= -M
				while((s1 + (i*p)) >= -M && (s1 + (i*p)) <= M){
					index1 = s1 + (i*p);	
//System.out.println(("soln1_" + p + " + " + i*p + " = " + index1));
					sieveArray[index1 + M] += logp[j];		// add M to offset sieve array starts at -M = sieveArray[0]
					i--;
				}
				i = 0;
				//find upper bound value for i such that soln1 +i*p <= M					
				while((s1 + (i*p)) >= -M && (s1 + (i*p)) <= M){						
					index1 = s1 + (i*p);
//System.out.println(("soln1_" + p + " + " + i*p + " = " + index1));
					sieveArray[index1 + M] += logp[j];
					i++;
				}
				if(p != 2){						// If p = 2, sieve only with solution1
					i = -1;
					while((s2 + (i*p)) >= -M && (s2 + (i*p)) <= M){
						index1 = s2 + (i*p);
//System.out.println(("soln2_" + p + " + " + i*p + " = " + index1));
						sieveArray[index1 + M] += logp[j];
						i--;
					}
					i = 0;
					while((s2 + (i*p)) >= -M && (s2 + (i*p)) <= M){
						index1 = s2 + i*p;
//System.out.println(("soln2_" + p + " + " + i*p + " = " + index1));
						sieveArray[index1 + M] += logp[j];
						i++;
					}
				}
			}//end if
		}//end for
		
		System.out.println("Sieving Stage complete...");

	/**************************** Trial Division Stage **************************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/


		System.out.println("Starting Trial Division Stage...");

		BigInteger negOne = BigInteger.valueOf(-1L);
		BigInteger xActual;

		int possiblySmooth = 0;									// Track number of elements in candidates array
		BigInteger[][] candidates = new BigInteger[M][2];

		for(int x = 0; x < sieveArray.length; x++){						// Scan sieve array for locations x indicating potential g_a,b(x) smooth
			xActual = BigInteger.valueOf((long) x - M);
			BigInteger testOperand = squareRootN.multiply(BigInteger.valueOf((long) M));// Test condition w/o error: 2x*sqrt(N)
			int testTerm = testOperand.bitLength();
			if(sieveArray[x] >= (testTerm-ERROR_TERM)){						
				candidates[possiblySmooth][0] = xActual;				// Mark x as possibly having g_a,b(x) as F-Smooth
				candidates[possiblySmooth][1] =						// Save corresponding polynomial value to array
					((a.multiply(xActual.pow(2))).add((two.multiply(b)).multiply(xActual))).add(c);
				possiblySmooth++;					
			}
		}

		for(int i = 0; i < possiblySmooth; i++){						// For each potentially smooth g(x), trial divide to check smoothness
			BigInteger checkMe = candidates[i][1];						// Number to trial divide for smoothness
			if(checkMe.compareTo(BigInteger.ZERO) < 0){					// Factor out -1 for trial division
				checkMe = checkMe.multiply(negOne);					// Make positive for trial division
			}
			for(int j = 0; j < factorBase.size(); j++){					// For each prime in the factor base
				int p = factorBase.get(j);
				while((checkMe.compareTo(BigInteger.ONE) > 0) && ((checkMe.mod(BigInteger.valueOf((long) p))).compareTo(BigInteger.ZERO) == 0))
					checkMe = checkMe.divide(BigInteger.valueOf((long) p));
			}
			if(checkMe.compareTo(BigInteger.ONE) == 0 || checkMe.compareTo(negOne) == 0){
				if(smoothRelFound == K+1)
					break;
				if(i >= S.length)
					break;
				S[smoothRelFound][0] = candidates[i][0];
				S[smoothRelFound][1] = candidates[i][1];
				S[smoothRelFound][2] = a;
				S[smoothRelFound][3] = b;
				smoothRelFound++;
			}
		} //end trial div outer for

	/* Print output for each sieving round - testing */

		System.out.println("  x --------- g_a,b(x)/a --- a " + "--- b ");
		for(int j = 0; j < S.length; j++){
			System.out.printf(" %-12d %-12d %-5d %d \n", S[j][0], S[j][1], S[j][2], S[j][3]);
			if(S[j][0] == null)
				break;
		}

	}//end outer while smooth rel < needed



		
	return S;

	}//end get smooth relations 



	/****************** Helper Methods *********************************************************************/


	/** Determine the factor base: the
	 * smooth primes such that N is a 
	 * quadratic residue (mod p). I.e. Legendre symbol
	 * (N/p) = 1
	 *
	 * @param F  the smoothness bound
	 * @param N  the composite being factored
	 * @return  a list containing the elements of the 
	 * factor base. 
	 */
	private static ArrayList<Integer> determineFactorBase(int F, BigInteger N) {
		// First get primes <= F into an array
		int[] primes = SieveOfEratosthenes.basicEratos(F);  
		ArrayList<Integer> factorBase = new ArrayList<Integer>();
		// Next determine which ones give leg (N/p)=1
		// by calling DetermineQuadRes
		for(int i = 0; i < primes.length; i++) {
			int p = primes[i];
			if(DetermineQuadRes.LegJacSym(N, p) == 1)
				factorBase.add(p);
		}
		

		return factorBase;
	}

	/** This algorithm lifts by one degree via Hensel's lemma; Based on the pseudocode 
	 * given in Algorithm 2.3.11 in Prime Numbers: Computational Perspective; 
	 * Given b and p such that b^2 \equiv N (mod p), this algorithm returns a new b such
	 * that b^2 \equiv N (mod p^2)
	 *
	 * @param b the square root of N (mod p)
	 * @param N a large composite number
	 * @param p the prime in the congruence
	 * @return bb the square root of N (mod p^2)
	 *
	 */
	private static BigInteger henselLift(BigInteger b, BigInteger p, BigInteger N){
		BigInteger x,y,z,bb;
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger negOne = BigInteger.valueOf(-1L);

		x = ((b.pow(2)).subtract(N)).divide(p);					// x = (b^2 - N)*p^-1
		z = (b.multiply(two)).modInverse(p);					// z = (2b)^-1 (mod p)
		y = ((x.multiply(negOne)).multiply(z)).mod(p);				// y = -xz (mod p)
		bb = b.add(y.multiply(p));
		return bb;

	}


	public static void main(String[] args){
		BigInteger b1 = new BigInteger("416064700201658306196320137931");		// 30 digit prime
     		BigInteger b2 = new BigInteger("513821217024129243948411056803");		// 30 digit prime
	//	BigInteger b1 = new BigInteger("1451730470513778492236629598992166035067");	// 40 digit prime
	//	BigInteger b2 = new BigInteger("2425967623052370772757633156976982469681");	// 40 digit prime
		
		BigInteger N = b1.multiply(b2);
		//N = 213782870618395542957440178483171817733685561437850435894593
	//	getSmoothRelations(BigInteger.valueOf(101868649L), 233);	
		getSmoothRelations(N, 70000);							// for 30 digit prime
	//	getSmoothRelations(N, 900000);							// for 40 digit number
	
	}//end main

}//end BasicQS

