import java.lang.Math;

/* @file SieveOfEratosthenes.java
 * @author Shane Sims, <shane.sims.ss@gmail.com>
 * 
 * Implementation of the Sieve of Eratosthenes based on 
 * the pseudo-code in Crandall and Pomerance's 'Prime Numbers'.
 *
 */

public class SieveOfEratosthenes{

	

	/**
	* Find primes in the interval (L,R) and print result to output.
	* @pre L and R are both even.
	* @pre R > L
	* @pre B|L-R (B is the number of odds in the interval)
	* @pre L > ceil(\sqrt(R))
	*  
	* @param L Left (non-inclusive) interval bound
	* @param R Right (non-inclusive) interval bound
	*/

	public static void practicalEratos(int L, int R) {
		// Get list of pi(P) primes
		int P = (int)Math.ceil(Math.sqrt(R));				// sqrt of R (ceiling).
		int[] piPprimes = basicEratos(P);				// Array of primes from 2 to P

		// Initialise the offsets
		int piP = piPprimes.length;					// Number of primes not exceeding P
		int[] Q = new int[piP];						// Array of offsets
		for(int k = 2; k <= piP; k++) {					// Start with second prime because we are sieving with odd numbers
			int negMod = (int)(-.5 * (L + 1 + piPprimes[k-1]));
			while(negMod < 0)					// Perform modular arithmetic
				negMod += piPprimes[k-1];
			Q[k-1] = negMod;					// Initialize offset q_k
		}


		// Process the blocks
		int T = L;
		int B = ((R - L)/2);						// Number of odd integers between L and R
		int[] blocks = new int[B]; 					// Array of primality bits; index i represents i+1'th odd num in (L,R)
		while(T < R) {
			for(int j = 0; j < blocks.length; j++)
				blocks[j] = 1;					// Init j+1'th odd num primality bit in interval to 1


			for(int k = 2; k <= piP; k++){				// For each multiple of prime from 3 to piP, change primality bit to 0 to rep composite 
				for(int j = Q[k-1]; j<B; j += piPprimes[k-1]){
					blocks[j] = 0;
				}
				int negMod = Q[k-1] - B;
				while(negMod < 0)
					negMod += piPprimes[k-1];
				Q[k-1] = negMod;
			}

			// Print out result
			for(int j = 0; j < B; j++) {
				if(blocks[j] == 1)
					System.out.println(T + 2*j + 1);
			}
			T = T + 2*B;
		}//end while

	}//end method





	//will return an array of the prime numbers from 2-N
	//such that 2 is in index 0, 3 in index 1 and so on.
	public static int[] basicEratos(int N) {
		int[] primes = new int[N+1];
		for(int i = 2; i < primes.length; i++) {
			primes[i] = 1;
		} //array initialised
		for(int j = 2; j< primes.length; j++) {
			if(primes[j] == 1) {
				int k = j + j;
				while(k<primes.length) {
					primes[k] = 0;
					k += j;	
				}//end while
			}//end if
		}
		int piP=0;
		for(int i = 0; i <primes.length; i++) {
			if(primes[i] == 1)
				piP++;
		}
		int[] finalPrimes = new int[piP];;
		int count = 0;//index to place next prime
		for(int i = 0; i < primes.length; i++){
			if(primes[i] == 1) {
				finalPrimes[count] = i;
				count++;
			}
		}
		return finalPrimes;
		
	}//end basicEratos





}// end class
