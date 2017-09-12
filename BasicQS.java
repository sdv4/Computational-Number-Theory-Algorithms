import java.lang.Math.*;
import java.util.*;

/** This is a implementation of the Sieving portion of the basic Quadratic Sieve 
* factorization algorithm. 
* This version of the QS includes enhancements suggested in Prime Numbers, such as
* sieving with approximate logarithms of factor base primes and centreing the values 
* of x around the number to be factored.  (TODO)							 //TODO Decide if we will sieve with small primes//
*
* Based on the pseudocode found in Prime Numbers: A Computational Perspective by 
* Pomerance and Crandall.
* The thesis 'Factoring Integers with the Self-Initializing Quadratic Sieve' by
* Scott P. Contini was consulted in the design of this code. The notation matches
* that found in Contini's thesis as far as possible so as to increase the didactic 
* value of this class. 
*
* @author Shane Sims <shane.sims.ss@gmail.com>
* @version 2 May 2016 
*
*/

public class BasicQS {

	//	private static final int F				TODO		
	//	private static final int SIEVE_INTERVAL			TODO
	//	private static final int ERROR_TERM			TODO

	/**
	* This method determine the smooth relations required to factor N.
	* 
	* @param N	The integer to be factored. Must be an odd composite
	* free of prime powers.   
	* @param F	The smoothness bound to which the returned relations will
	* be bound by.
	* @return 	An array with each element containing an array with length 2, 
	* representing an ordered pair of the form (x, x^2 - N), where the value 
	* x leads x^2-n to be B-smooth.
	*
	*/
	public static int[][] getSmoothRelations(int N, int F){
	
		int[][] S;									// Array that will hold the smooth relations when returned.
		ArrayList<Integer> factorBase;							// Primes which will be used in the sieving stage. 
		int[] modSqrts;									// Will hold t in t^2 \equiv N (mod p) for each p in Factor Base
		int[] solution1;								// Will hold t - b \equiv N (mod p) for each prime in factor base
		int[] solution2;								// Will hold -t - b \equiv N (mod p) for each prime in factor base
		int[] logp;									// Approximate logarithm for each prime in the factor base
		int K;										// Method will return when K+1 smooth relations have been found.
		int offset = 0;									// Used if/when sieving stage returned to
	//
	//check preconditions of N are met here
	//


	/* Compute Startup Data */								// Note: the term 'Initialization' will be reserved for use in SIQS
		double squareRootN = Math.sqrt(N);
		int b = (int) Math.ceil(squareRootN);						// Compute the constant term \ceil(\sqrt(N)). Used in the Sieving Stage
		factorBase = determineFactorBase(F, N);						// Initialize factorBase.
		K = factorBase.size();
		S = new int[K+1][2];								// Elements are pairs (x,(x+b)^2-N)				
		logp = new int[K];
		modSqrts = new int[K];								// Initialize modSqrts to have space for each element of factorBase
		modSqrts[0] = 1;								// TonelliSqrtModP requires odd prime, so only even prime hard coded
		solution1 = new int[K];
		solution2 = new int[K];

		for(int i = 1; i < modSqrts.length; i++) {					// Solve t in t^2 \equiv N (mod p) for each p in factorBase
			modSqrts[i] = TonelliSqrtModP.SqrtModP(N, factorBase.get(i));
		}
		
		for(int i = 0; i < solution1.length; i++){					// Calc t-b (mod p) and -t-b (mod p) for each p in factor base and store
			solution1[i] = Math.floorMod((modSqrts[i] - b), factorBase.get(i));
			solution2[i] = Math.floorMod(((-1*modSqrts[i]) - b), factorBase.get(i));
		}

		for(int i = 0; i < logp.length; i++) {
			logp[i] = (int) Math.round(Math.log(factorBase.get(i))/Math.log(2));	// Get approximate log_2(p) using change of base for each factor base prime
		}

	/* Sieving Stage */
	
		int smoothRelFound = 0;								// Tracks the number of smooth relations verified by trial division.
		int size = F * F;								// Set interval size F^2 as per Pomerance suggestion of optimal
		int index1 = 0;									// Will hold latest value of solution1_p + ip
		int index2 = 0;									// Will hold latest value of solution2_p + ip
		int[] ithMultipleP = new int[K];						// Will hold i'th multiple of p sieved with so fa
	while(smoothRelFound <= K){
		System.out.println("Smooth relations found: " + smoothRelFound);
		System.out.println("Size of factor base: " + K);
		int[] sieveArray = new int[size];						// Initialize sieve array with 0's.
		for(int j = 0; j < factorBase.size(); j++) {					// For each prime in the factor base (!=2)
			int i = ithMultipleP[j];
			while(true){
				index1 = solution1[j] + i*factorBase.get(j) - offset;		// Increment location to add factor base prime
				index2 = solution2[j] + i*factorBase.get(j) - offset;		// Increment location to add factor base prime
				if((index1) >= size || index1 < 0 || index2 >= size || index2 < 0){// Stop if index will be out of array bounds
					ithMultipleP[j] = i;					// Save last multiple of p used in index 1 or 2	
					break;
				}
				else {
					sieveArray[index1] += logp[j];				// Otherwise add approx log of p to index in sieve array
					if(factorBase.get(j) != 2){				// Only sieve with solution1_p if p = 2
						sieveArray[index2] += logp[j];

					}
					i++;
				}
			}//end while
		}//end for
	

	/* Trial Division Stage */

		int possiblySmooth = 0;
		int[][] candidates = new int[size][2];
		for(int x = 0; x < sieveArray.length; x++){					// Scan sieve array for locations x indicating potential g(x) smooth

			double testTerm = Math.log((x + offset)*(squareRootN))/Math.log(2);	// Test condition for values at location x: log_2(2*x*\sqrt(N))
			System.out.println("Sieve array index " + (x+offset) + ": " + sieveArray[x] + "    testTerm:" + (testTerm-6.5));

			if(sieveArray[x] >= (testTerm-6.5)){					//TODO Needs experimentation to optimize
				
				candidates[possiblySmooth][0] = x + offset;			// Mark x as possibly having g(x) as F-Smooth
				candidates[possiblySmooth][1] = (x + offset + b)*(x + offset + b) - N;
				possiblySmooth++;					
			}
		}
		
		for(int i = 0; i < possiblySmooth; i++){					// For each potentially smooth g(x)
			int checkMe = candidates[i][1];						// Number to trial divide for smoothness
			for(int j = 0; j < factorBase.size(); j++){				// For each prime in the factor base
				int p = factorBase.get(j);
				while((checkMe > 1) && (checkMe % p == 0))
					checkMe /= p;
			}
			if(checkMe == 1){
				if(smoothRelFound == K+1)
					break;
				S[smoothRelFound][0] = candidates[i][0];
				S[smoothRelFound][1] = candidates[i][1];
				smoothRelFound++;
				//need to break if S is full
			}
		} //end trial div outer for
		offset += size;									// Adjust offset for next sieve interval

	/* Print output for each sieving round - testing */

			System.out.println("x---------g(x)");
		for(int j = 0; j < S.length; j++)
			System.out.println(S[j][0] + "         " + S[j][1]);

	}//end outer while



		
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
		getSmoothRelations(101868649, 233);	
	}//end main

}//end BasicQS
