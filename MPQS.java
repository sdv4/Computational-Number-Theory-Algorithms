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
* @author Shane Sims <shane.sims.ss@gmail.com>
* @version 28 June 2017 
*
*/

public class MPQS {	

	/**
	* This method determine the smooth relations required to factor N.
	* 
	* @param N	The integer to be factored. Must be an odd composite
	* free of prime powers.   
	* @param F	The smoothness bound to which the returned relations will
	* be bound by.
	* @return 	An array with each element containing an array with length 2, 
	* representing an ordered pair of the form (x, (ax+b)^2 - N), where the value 
	* x leads (ax + b)^2-n to be B-smooth.
	*
	*/
		public static int[][] getSmoothRelations(int N, int F){
	
		int[][] S;									// Array that will hold the smooth relations when returned.
		ArrayList<Integer> factorBase;							// Primes which will be used in the sieving stage. 
		int[] modSqrts;									// Will hold t in t^2 \equiv N (mod p) for each p in Factor Base
		int[] solution1;								// Will hold a^-1(t - b) \equiv N (mod p) for each prime in factor base
		int[] solution2;								// Will hold a^-1(-t - b) \equiv N (mod p) for each prime in factor base
		int[] logp;									// Approximate logarithm for each prime in the factor base
		int K;										// Method will return when K+1 smooth relations have been found.
		int offset = 0;									// Used if/when sieving stage returned to
		int a;										// Term 'a' in g_a,b(x) = (ax+b)^2 - N 
		int M;										// Term in bound [-M,M]
		int b;										// Term 'b' in g_a,b(x)
	//
	//check preconditions of N are met here
	//


	/* Compute Startup Data */
		M = F;		
		double squareRootN = Math.sqrt(N);		// can delete???
		b = (int) Math.ceil(squareRootN);						// Compute the constant term \ceil(\sqrt(N)). Used in the Sieving Stage
		factorBase = determineFactorBase(F, N);						// Initialize factorBase.
		K = factorBase.size();
		S = new int[K+1][4];								// Elements are pairs (x,(x+b)^2-N)				
		logp = new int[K];
		modSqrts = new int[K];								// Initialize modSqrts to have space for each element of factorBase
		modSqrts[0] = 1;								// TonelliSqrtModP requires odd prime, so only even prime hard coded
		solution1 = new int[K];
		solution2 = new int[K];
System.out.println("tmem_" + factorBase.get(0) + " = " + modSqrts[0]);

		for(int i = 1; i < modSqrts.length; i++) {					// Solve t in t^2 \equiv N (mod p) for each p in factorBase
			modSqrts[i] = TonelliSqrtModP.SqrtModP(N, factorBase.get(i));
System.out.println("tmem_" + factorBase.get(i) + " = " + modSqrts[i]);
			
		}
		for(int i = 0; i < logp.length; i++) {
			logp[i] = (int) Math.round(Math.log(factorBase.get(i))/Math.log(2));	// Get approximate log_2(p) using change of base for each factor base prime
		}





	/* Initialization Stage */
		int smoothRelFound = 0;								// Tracks the number of smooth relations verified by trial division.
		int q = 0;
	while(smoothRelFound <= K){


		// Find q and set a
		if(q == 0){
			q = (int) Math.round(Math.sqrt(Math.sqrt(2 * N)/M));			// q approx sqrt(sqrt(2N)/M)
		}
		else
			q++;
		while(!factorBase.contains(q))							// Get F-smooth a with Jacobi (q/N)=1 
			q++;

		a = q * q;									// Set a = q^2

		// Compute b
		int NmodA = N % a;								// Reduce N (mod a)		
    		int bPrime = TonelliSqrtModP.SqrtModP(N, q);					// Get b'=b in b^2 \equiv N (mod q)

		while(((bPrime * bPrime) % a) != NmodA)
			bPrime += q;
		b = bPrime;

		// Compute soln1_p and soln2_p for each p in factor base

		for(int i = 0; i < solution1.length; i++){
			if(q != factorBase.get(i)){
			BigInteger biga = new BigInteger(Integer.toString(a));
			BigInteger bigp = new BigInteger(Integer.toString(factorBase.get(i)));

			BigInteger bigaModInv = biga.modInverse(bigp);
			int aInverseModp = bigaModInv.intValue();
			solution1[i] = Math.floorMod((aInverseModp*(modSqrts[i] - b)), factorBase.get(i)); 
			solution2[i] = Math.floorMod(((aInverseModp*(-1*modSqrts[i])) - b), factorBase.get(i));
			}
			else{
				solution1[i] = -1;
				solution2[i] = -1;		// just marking unusable solutions.
			}
		}
		System.out.println("a: " + a + "       b: " + b);
		for(int i = 0; i< solution1.length; i++)
			System.out.println("solution1: " + solution1[i]);

System.out.println("Initialization Stage complete...");		
System.out.println("Press Enter to continue.");
Scanner readinput = new Scanner(System.in);	
readinput.nextLine();
	
	/* Sieving Stage */
	
		int index1 = 0;									// Will hold latest value of solution1_p + ip
		int index2 = 0;									// Will hold latest value of solution2_p + ip
		int[] ithMultipleP = new int[K];						// Will hold i'th multiple of p sieved with so fa
		
		
		System.out.println("Smooth relations found: " + smoothRelFound);
		System.out.println("Size of factor base: " + K);
		int[] sieveArray = new int[2*M + 1];						// Initialize sieve array with 0's.
		for(int j = 0; j < factorBase.size(); j++) {					// For each prime in the factor base.
		int k = 0;	
			while(true){				// TODO resume debugging here. Need to account for index 0 of sieve array being index 0 and size of array being 2M+1
				int p = factorBase.get(j);
	System.out.println("P: " + p);
				int s1 = solution1[j];
				int s2 = solution2[j];
				if(s1 != -1){
					int i = -1;
				//find lower bound value for i such that soln1 +i*p >= -M
					while(s1 + i*p >= -M && s1 + i*p <= M){
						index1 = s1 + i*p;
	System.out.println("solution1: " + s1);
	System.out.println("index1: " + index1);
						
	
						sieveArray[index1 + M] += logp[j];		// add M to offset sieve array starts at -M = sieveArray[0]
						i--;
					}
					i = 0;
					while(s1 + i*p <= M){						
						index1 = s1 + i*p;
	System.out.println("index1: " + index1);


						sieveArray[index1 + M] += logp[j];
						i++;
					}
					if(p != 2){						// If p = 2, sieve only with solution1
						i = -1;
						while(s2 + i*p >= -M && s2 + i*p <= M){
							index1 = s2 + i*p;
	System.out.println("index1*: " + index1);
	
							sieveArray[index1 + M] += logp[j];
							i--;
						}
						i = 0;
						while(s2 + i*p <= M){
							index1 = s2 + i*p;

	System.out.println("index1**: " + index1);
	
							sieveArray[index1 + M] += logp[j];
							i++;
						}
					}
				}//end if
				break;
			}//end while
		}//end for
	

	/* Trial Division Stage */

		int possiblySmooth = 0;
		int[][] candidates = new int[M][2];
		for(int x = 0; x < sieveArray.length; x++){						// Scan sieve array for locations x indicating potential g_a,b(x) smooth

			double testTerm = Math.log(M*(squareRootN))/Math.log(2);			// Test condition for values at location x: log_2(2*x*\sqrt(N))
			System.out.println("Sieve array index " + x + ": " + sieveArray[x] + "    testTerm:" + (testTerm-4.5));

			if(sieveArray[x] >= (testTerm-4.5)){						//TODO Needs experimentation to optimize
				
				candidates[possiblySmooth][0] = x - M;					// Mark x as possibly having g_a,b(x) as F-Smooth
	System.out.println("Possibly smooth: x = " + candidates[possiblySmooth][0]);
				candidates[possiblySmooth][1] = (a*(x-M) + b)*(a*(x-M) + b) - N;
				possiblySmooth++;					
			}
		}

		for(int i = 0; i < possiblySmooth; i++){						// For each potentially smooth g(x)
			int checkMe = candidates[i][1];							// Number to trial divide for smoothness
			if(checkMe < 0){								// Factor out -1 for trial division
				checkMe *= -1;								// Make positive for trial division
			}
			for(int j = 0; j < factorBase.size(); j++){					// For each prime in the factor base
				int p = factorBase.get(j);
				while((checkMe > 1) && (checkMe % p == 0))
					checkMe /= p;
			}
			if(checkMe == 1 || checkMe == -1){
				if(smoothRelFound == K+1)
					break;
				S[smoothRelFound][0] = candidates[i][0];
				S[smoothRelFound][1] = candidates[i][1];
				S[smoothRelFound][2] = a;
				S[smoothRelFound][3] = b;
				smoothRelFound++;
				//need to break if S is full
			}
		} //end trial div outer for

	/* Print output for each sieving round - testing */

			System.out.println("  x --------- g_a,b(x) --- a " + "--- b ");
		for(int j = 0; j < S.length; j++)
			System.out.printf(" %-12d %-12d %-5d %d \n", S[j][0], S[j][1], S[j][2], S[j][3]);

	}//end outer while smooth rel < needed



		
		return S;

	}//end get smooth relations 


	/**
	* This method is the same as above but calculates the smoothness bound 
	* B based on the heuristic arguments presented in Prime Numbers.
	* 
	* @param N	The integer to be factored. Must be an odd composite
	* free of prime powers.   
	*
	*/
//	public static int[] getSmoothRelations(int N){}


	// Helper Methods //


	/* Determine the factor base: the
	 * smooth primes such that N is a 
	 * quadratic residue (mod p). I.e. Legendre symbol
	 * (N/p) = 1
	 *
	 * @param F  the smoothness bound
	 * @param N  the composite begin factored
	 * @return  a list containing the elements of the 
	 * factor base. 
	 */
	private static ArrayList<Integer> determineFactorBase(int F, int N) {
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

	public static void main(String[] args){
		getSmoothRelations(62113, 37);	
	}//end main

}//end BasicQS
