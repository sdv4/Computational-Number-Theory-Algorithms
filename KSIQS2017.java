import java.lang.Math.*;
import java.util.*;
import java.math.BigInteger;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.io.*;
import java.lang.StringBuilder.*;

/** Implementation of the Sieving portion of the Self Initializing 
* Quadratic Sieve factoring algorithm, with improved Initialization being done
* according to T. Kleinjung's method found in 'Quadratic Sieving' (2016).
* Built on code written for my implementation of SIQS and based on the pseduocode provided in 
* the thesis 'Factoring Integers with the Self-Initializing Quadratic Sieve' by
* Scott P. Contini.
* 
* This version is for numbers N of approximately 60 digits.
*
*
* @author Shane Sims <shane.sims.ss@gmail.com>
* @version 9 August 2017 
*
*/

public class KSIQS2017 {

	public static final int ERROR_TERM = 27;
	public static final int M = 360000;
	public static final int DESIRED_NUM_RELATIONS = 30;
	public static final int MIN_FACTOR_IN_A = 2000;						// 170th prime in FB as per Contini reccomendation
	public static final int n = 10;								// desired number of factors of term 'a'. Choose so that more than half of the 2^n-1 polys needed
	
	public static final int MAX_FACTOR_IN_A = 2353;						// Hard coded value of (sqrt(N)/M)^(1/10) - needs to change if N or n change. N 80 dig, n = 10

	/**
     	 * Hide utility class.
     	 */
    	private KSIQS2017() {}


	/**
       	 * Determine the smooth relations required to factor N by QS.
	 *
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
	public static void getSmoothRelations(BigInteger N, int F)
	{


//		Scanner scanner = new Scanner(System.in);
 // 		System.out.println("Press Enter to begin\t");
 //    	    	scanner.nextLine();

		// Method Variables

		BigInteger[][] S;								// Array that will hold the smooth relations when returned.
		ArrayList<Integer> factorBase;							// Primes which will be used in the sieving stage. 
		int[] tmem_p;									// Will hold t in t^2 \equiv N (mod p) for each p in Factor Base
		int[] solution1;								// Will hold a^-1(t - b) \equiv N (mod p) for each prime in factor base
		int[] solution2;								// Will hold a^-1(-t - b) \equiv N (mod p) for each prime in factor base
		int[] log_p;									// Approximate logarithm for each prime in the factor base
		int K;										// Method will return when K+1 smooth relations have been found.
		BigInteger a;									// Term 'a' in g_a,b(x) = (ax+b)^2 - N 
		BigInteger b;									// Term 'b' in g_a,b(x)
		BigInteger c;									// Term 'c' in g_a,b(x) = (ax+b)^2 - N = a(ax^2 + 2bx + c)
		ArrayList<Integer> aFactors;							// The prime factors, the product of which forms a
		BigInteger aInv_p[];								// a^-1 (mod p) for each p in FB
		BigInteger squareRootN;
		BigInteger[] B;									// B_j values where j is polynomial number
		int smoothRelFound;								// Tracks the number of smooth relations verified by trial division.
		int polynomialsInitialized;
		ArrayList<int[]> sievingEvents;							// All sieving events induced by p > 2M + n/2
		BigInteger[] allBValues;

		// Method Constants
		
		BigInteger negOne = BigInteger.valueOf(-1L);
		BigInteger one = BigInteger.valueOf(1L);
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger bigM = BigInteger.valueOf((long) M);	
		

/******************** Compute Startup Data *********************************************************************************************************************************************************/
/***************************************************************************************************************************************************************************************************/                                                                                     

		System.out.println("Computing Startup Data...");
		smoothRelFound = 0;								// Tracks the number of smooth relations verified by trial division.
		polynomialsInitialized = 0;							// Number of polynomials initialized so far
		squareRootN = BigIntegerSqrt.bigIntegerSqrtCeiling(N);		
		factorBase = determineFactorBase(F, N);						// Initialize factorBase.
		K = factorBase.size();								// Save factor base size. Need K+1 smooth valued poly to ensure factorization
		S = new BigInteger[K+1][4];							// Elements are pairs (x,(x+b)^2-N)				
		log_p = new int[K];								
		tmem_p = new int[K];								// Initialize modSqrts to have space for each element of factorBase
		tmem_p[0] = 1;									// TonelliSqrtModP requires odd prime, so only even prime hard coded
		solution1 = new int[K];
		solution2 = new int[K];
		aInv_p = new BigInteger[K];
		b = BigInteger.ZERO;
		int index = 0; 									// Will hold index of sieving event from SievingArray for last event inspected

		for(int i = 1; i < tmem_p.length; i++) {					// Solve t in t^2 \equiv N (mod p) for each p in factorBase
			tmem_p[i] = TonelliSqrtModP.SqrtModP(N, factorBase.get(i));
		}

		for(int i = 0; i < log_p.length; i++) {
			log_p[i] = (int) Math.round(Math.log(factorBase.get(i))/Math.log(2));	// Get approximate log_2(p) using change of base for each factor base prime
		}	

		System.out.println("Compute Startup Data complete...");

/******************** Improved Initialization ******************************************************************************************************************************************************/
/***************************************************************************************************************************************************************************************************/

		System.out.println("Starting Initialization Stage...");
			
		/* Compute and set 'a' term */
		
		a = one;
		BigInteger targetA = BigInteger.ZERO;	
		targetA = N.multiply(two);						
		targetA = BigIntegerSqrt.bigIntegerSqrtCeiling(targetA);
		targetA = targetA.divide(bigM);
		aFactors = new ArrayList<Integer>();		
		int factorIndex = 0;
		while(factorBase.get(factorIndex) <= MAX_FACTOR_IN_A){				// Find index of first factor base prime > MIN_FACTOR_IN_A
			factorIndex++;
		}

		BigInteger tempFactor = BigInteger.valueOf((long) factorBase.get(factorIndex));
		BigInteger tempA = tempFactor;							// Temp. value of a to ensure size is as close to target w/o going over.
		int aNumFactors = 0;								// Track number of factors in a, since .size() not necessarily accurate.


		while(aNumFactors < n && factorIndex < factorBase.size()){			
			a = tempA;
			aFactors.add(tempFactor.intValue());
			aNumFactors++;
			factorIndex--;
			tempFactor = BigInteger.valueOf((long) factorBase.get(factorIndex));	// Get value in factorbase at index 'factorIndex'
			tempA = tempA.multiply(tempFactor);							
		}

		/* Calculate and save all signs permutations to make 'b' terms */

		B = new BigInteger[n];
		// Get B_i for all i factors of a		
		BigInteger gamma;
		BigInteger bigTmemp;
		BigInteger tempTerm2;
		BigInteger bigAFactor;
		
		for(int i = 0; i < aNumFactors; i++)
		{
			tempTerm2 = BigInteger.ONE;
			bigAFactor = BigInteger.valueOf((long) aFactors.get(i));

			for(int j = 0; j < aFactors.size(); j++)
			{
				if(i != j)
					tempTerm2 = tempTerm2.multiply(BigInteger.valueOf((long) aFactors.get(j)));			// compute a/bigAFactor
			}					
			bigTmemp = BigInteger.valueOf((long) tmem_p[factorBase.indexOf(aFactors.get(i))]);				// get tmem_p, where p is current factor of 'a'
			tempTerm2 = tempTerm2.modInverse(bigAFactor);
			gamma = (bigTmemp.multiply(tempTerm2)).mod(bigAFactor);
			float gammaFloat = gamma.floatValue();
			float qover2 = bigAFactor.floatValue()/2;
			if(i == 0 && gammaFloat > bigAFactor.intValue())
				gamma = bigAFactor.subtract(gamma);
			else if(gammaFloat > qover2)
				gamma = bigAFactor.subtract(gamma);
			else{}
			B[i] = gamma.multiply(a.divide(bigAFactor));	
		}

		int numOfPossiblePolys = (int) Math.pow(2, n - 1);
	
		/* For each prime p > 2M + n/2, p in FB, detect all sieving events induced by (p, tmem_p) */				
		sievingEvents = new ArrayList<int[]>();
		int indexOfFirstBigPrime;
		if(F > 2.2*M)
		{	
			indexOfFirstBigPrime = 0;												// will hold index of first prime >2M in FB
			while(indexOfFirstBigPrime < factorBase.size() && factorBase.get(indexOfFirstBigPrime) <= (2.2*M))
				indexOfFirstBigPrime++;


			int indexOfBigPrime = indexOfFirstBigPrime;										// Save to mark max prime in FB for regular initialization.
			sievingEvents = discoverSievingEventsInduced(N, F, a, factorBase.get(indexOfBigPrime), tmem_p[indexOfBigPrime], B);
			sievingEvents.addAll(discoverSievingEventsInduced(N, F, a, factorBase.get(indexOfBigPrime), (-1*tmem_p[indexOfBigPrime]), B));
			indexOfBigPrime++;
		
		
		
			while(indexOfBigPrime < factorBase.size())
			{
				sievingEvents.addAll(discoverSievingEventsInduced(N, F, a, factorBase.get(indexOfBigPrime), tmem_p[indexOfBigPrime], B));
				sievingEvents.addAll(discoverSievingEventsInduced(N, F, a, factorBase.get(indexOfBigPrime), (-1*tmem_p[indexOfBigPrime]), B));
				
				indexOfBigPrime++;
			}
		}
		else
			indexOfFirstBigPrime = factorBase.size();

		sievingEvents = sortTupleBySecondElement(sievingEvents);
		
		while(polynomialsInitialized < numOfPossiblePolys && smoothRelFound < DESIRED_NUM_RELATIONS){

		/* First Polynoial Initialization */
			if(polynomialsInitialized == 0){
		
				
				// Set first value of b and first solution1 and solution2
				b = determineValueOfB(polynomialsInitialized, B, n);
				BigInteger bigaInv_p;
				BigInteger fb_p;
				BigInteger bigTemp;
				BigInteger bigtmem_p;

				for(int i = 0; i < indexOfFirstBigPrime; i++){									// For each p in FB > 2M + n/2 such that 
					int p = factorBase.get(i);
					if(!aFactors.contains(p)){										// NOT p|a
						fb_p = BigInteger.valueOf((long) p);			
						bigaInv_p = a.modInverse(fb_p);									// compute a^-1 (mod p).
						aInv_p[i] = bigaInv_p;
						bigtmem_p = BigInteger.valueOf((long) tmem_p[i]);
						solution1[i] = ((bigaInv_p.multiply(bigtmem_p.subtract(b))).mod(fb_p)).intValue();
						solution2[i] = ((bigaInv_p.multiply((BigInteger.ZERO.subtract(bigtmem_p)).subtract(b))).mod(fb_p)).intValue();

					}//end if
					else{
						solution1[i] = -1;
						solution2[i] = -1;										// Mark unusable solutions.
					}

				}//end for loop precomputing values for subsequent poly init
				c = calculateC(a, b, N);

				polynomialsInitialized++;				// Increment number of polynomials initialized
			}//end first poly initialization

		/* Subsequent Polynoial Initialization */

		else
		{
			// Switch b
			b = determineValueOfB(polynomialsInitialized, B, n);

			
			// Compute solution1 and solution2 for this polynomial
			int p;	
			BigInteger bigtmem_p;
			BigInteger fb_p;
			for(int i = 0; i < indexOfFirstBigPrime; i++){
				p = factorBase.get(i);
				if(!aFactors.contains(p)){				
					fb_p = BigInteger.valueOf((long) factorBase.get(i));
					bigtmem_p = BigInteger.valueOf((long) tmem_p[factorBase.indexOf(p)]);
					solution1[i] = ((aInv_p[i].multiply((bigtmem_p.subtract(b)))).mod(fb_p)).intValue();
					solution2[i] = ((aInv_p[i].multiply(((bigtmem_p.multiply(negOne)).subtract(b)))).mod(fb_p)).intValue();
				}
				else{
					solution1[i] = -1;
					solution2[i] = -1;							// Mark unusable solutions.
				}
			}
			
System.out.println("a: " + a + "       b: " + b);
System.out.println("If a and b are correct, this should be 0:  " + (b.multiply(b).mod(a)).subtract(N.mod(a)));

			c = calculateC(a, b, N);
			polynomialsInitialized++;

		}
		System.out.println("Initialization Stage complete...");

/******************** Sieving Stage ****************************************************************************************************************************************************************/
/***************************************************************************************************************************************************************************************************/

		System.out.println("Starting Sieving Stage...");
		int[] sieveArray = new int[2*M + 1];						// Initialize sieve array with 0's.
		
		for(int j = 0; j < indexOfFirstBigPrime; j++) {					// For each prime in the factor base.
			int p = factorBase.get(j);
			int s1 = solution1[j];
			int s2 = solution2[j];
			int y = 0;

			if(s1 != -1){
				int i = -1;
				//find lower bound value for i such that soln1 +i*p >= -M
				while((s1 + (i*p)) >= -M && (s1 + (i*p)) <= M){
					y = s1 + (i*p);
					sieveArray[y + M] += log_p[j];			// add M to offset sieve array starts at -M = sieveArray[0]	
					i--;
				}
				i = 0;
				//find upper bound value for i such that soln1 +i*p <= M					
				while((s1 + (i*p)) >= -M && (s1 + (i*p)) <= M){						
					y = s1 + (i*p);
					sieveArray[y + M] += log_p[j];			// add M to offset sieve array starts at -M = sieveArray[0]
					i++;
				}
				if(p != 2){							// If p = 2, sieve only with solution1
					i = -1;
					while((s2 + (i*p)) >= -M && (s2 + (i*p)) <= M){
						y = s2 + (i*p);
						sieveArray[y + M] += log_p[j];			// add M to offset sieve array starts at -M = sieveArray[0]
						i--;
					}
					i = 0;
					while((s2 + (i*p)) >= -M && (s2 + (i*p)) <= M){
						y = s2 + i*p;
						sieveArray[y + M] += log_p[j];			// add M to offset sieve array starts at -M = sieveArray[0]	
						i++;
					}
				}
			}//end if

		}//end for

		if(!sievingEvents.isEmpty())
		{
			int j = polynomialsInitialized - 1;
			int p;
			int y = 0;
			
			//Find index of first occurence of polynomial j
			while(index < sievingEvents.size() && sievingEvents.get(index)[1] == j)
			{
				y = sievingEvents.get(index)[2];
				p = sievingEvents.get(index)[0];
				if(y >= -M && y <= M)
				{	
					sieveArray[y + M] += log_p[factorBase.indexOf(p)];
					index++;					
				}
				else if ((y - p) >= -M && (y - p) <= M)
				{
					sieveArray[(y - p) + M] += log_p[factorBase.indexOf(p)];
					index++;				
				}
				else if ((y + p) >= -M && (y + p) <= M)
				{
					sieveArray[(y + p) + M] += log_p[factorBase.indexOf(p)];
					index++;				
				}
				else{
					index++;					
				}
			}

		}
		

		System.out.println("Sieving Stage complete...");


/********************************** Trial Division Stage *******************************************************************************************************************************************/
/***************************************************************************************************************************************************************************************************/

		System.out.println("Starting Trial Division Stage...");

		BigInteger xActual;

		int possiblySmooth = 0;									// Track number of elements in candidates array
		BigInteger[][] candidates = new BigInteger[M][2];
		BigInteger testOperand = squareRootN.multiply(BigInteger.valueOf((long) M));		// Test condition w/o error: 2x*sqrt(N)
		int testTerm = testOperand.bitLength();

		for(int x = 0; x < sieveArray.length; x++){						// Scan sieve array for locations x indicating potential g_a,b(x)a smooth
			xActual = BigInteger.valueOf((long) x - M);	

			if(sieveArray[x] >= (testTerm - ERROR_TERM)){						
				candidates[possiblySmooth][0] = xActual;				// Mark x as possibly having g_a,b(x)/a as F-Smooth
				candidates[possiblySmooth][1] =						// Save corresponding g_a,b(x)/a value to array
					((a.multiply(xActual.pow(2))).add((two.multiply(b)).multiply(xActual))).add(c);
				possiblySmooth++;					
			}
		}


		System.out.println("*******************************     number of values to trial divide: " + possiblySmooth);
		for(int i = 0; i < possiblySmooth; i++){						// For each potentially smooth g(x)/a, trial divide to check smoothness
			BigInteger checkMe = candidates[i][1];						// Number to trial divide for smoothness
			if(checkMe.compareTo(BigInteger.ZERO) < 0){					// Factor out -1 for trial division
				checkMe = checkMe.multiply(negOne);					// Make positive for trial division
			}
			for(int m = 0; m < factorBase.size(); m++){					// For each prime in the factor base
				int p = factorBase.get(m);
				while((checkMe.compareTo(BigInteger.ONE) > 0) && ((checkMe.mod(BigInteger.valueOf((long) p))).compareTo(BigInteger.ZERO) == 0))
					checkMe = checkMe.divide(BigInteger.valueOf((long) p));
			}
			if(checkMe.compareTo(BigInteger.ONE) == 0 || checkMe.compareTo(negOne) == 0){
				S[smoothRelFound][0] = candidates[i][0];
				S[smoothRelFound][1] = candidates[i][1];
				S[smoothRelFound][2] = a;
				S[smoothRelFound][3] = b;
				smoothRelFound++;
			}
		}

		/* Print output for each sieving round - testing */

		System.out.println("  x --------- g_a,b(x)/a --- a " + "--- b ");
		for(int n = 0; n < S.length; n++){
			System.out.printf(" %-12d %-12d %-5d %d \n", S[n][0], S[n][1], S[n][2], S[n][3]);
			if(S[n][0] == null)
				break;
		}



	}//end outer while smooth rel < needed
System.out.println("Polynomials initialized: " + polynomialsInitialized);
System.out.println("Size of FB: " + factorBase.size());
	}// End getSmoothRelations













////////////Helper Methods////////////////////////////



	/** Takes in ArrayList of sieving event tuples (p,j,y) and sorts them by j.
	 *
	 * ref:https://stackoverflow.com/questions/39323208/sorting-an-array-list-of-int-arrays-in-terms-of-the-first-and-the-second-element
	 * 
	 */
	private static ArrayList<int[]> sortTupleBySecondElement(ArrayList<int[]> SievingEvents){
		ArrayList<int[]> se = SievingEvents;
		Collections.sort(se, new Comparator<int[]>() {
         		public int compare(int[] a, int[] b) {
         	        	return a[1]-(b[1]);
         			
         		}
    		});
		return se;

	}

	/** Takes in ArrayList of sieving event tuples (p,j,y) and sorts them by j.
	 *
	 * ref:https://stackoverflow.com/questions/39323208/sorting-an-array-list-of-int-arrays-in-terms-of-the-first-and-the-second-element
	 * 
	 */
	private static ArrayList<int[]> sortTupleByFirstElement(ArrayList<int[]> SievingEvents){
		ArrayList<int[]> se = SievingEvents;
		Collections.sort(se, new Comparator<int[]>() {
         		public int compare(int[] a, int[] b) {
         	        	return a[0]-b[0];
         			
         		}
    		});
		return se;

	}
	/** Implementation of T.Kleinjung's algorithm 2.2 in referenced paper. Note that 
	 * part (1) done in getSmoothRelations, as the values are common to/needed for
	 * other parts of initialization and need only be done once. 
	 *
	 *
	 * @param N	number to be factored
	 * @param F	smoothness bound
	 * @param A	product of primes in FB, 'a' term in KSIQS
	 * @param p	prime in FB > 2M whose siving events will be returned
	 * @param s	a square root of N (mod p). As in tmem_p in KSIQS/SIQS
	 * @param B	the partial terms whose sum, up to value of sign, form 'b' term
	 *
	 * @return	array of 3-tuples (p,j,x) containing all sieving events induced by (p,s).
	 */
	private static ArrayList<int[]> discoverSievingEventsInduced(BigInteger N, int F, BigInteger A, int p,
		       	int s, BigInteger[] B)
	{
		
		BigInteger negOne = BigInteger.valueOf(-1L);		
		BigInteger two = BigInteger.valueOf(2L);		
		ArrayList<int[]> sievingEvents = new ArrayList<int[]>();
		double[] beta = new double[n];									// beta_k = (B_k*(2A)^(-1) (mod p)) + B_k/2A
		int n1, n2;
		ArrayList<int[]> Tau1, Tau2;
		

		/* Part (2) of Alg 2.2 - Note: beta values are not required for this simplified case */
		double modBOverA;
		BigInteger bigP = BigInteger.valueOf((long) p);
		BigInteger invAmodp = (A).modInverse(bigP);							// (2A)^(-1) mod p

		

		int[] modBOverACollection = new int[n];
		for(int k = 0; k < n; k++)
		{
			modBOverA = ((B[k].multiply(invAmodp)).mod(bigP)).doubleValue();			//first term
			modBOverACollection[k] = (int) modBOverA;						//save for use in part 4a
		}


		/* Part (3) of Alg 2.2 */
		n1 = (int) Math.ceil(n/2.0);
		n2 = (int) Math.floor(n/2.0);

		
		/* Part (4a) of Alg 2.2 - Implemented according to special case of M > n - See remark 2.7 in Kleinjung's paper */	
		BigInteger sumAllBOverA = BigInteger.ZERO;
		BigInteger bigS = BigInteger.valueOf((long) s);
		BigInteger modSOverA = (bigS.multiply(invAmodp)).mod(bigP);
		
	       	BigInteger rem_p;	
		int[] tupleForTau1;
		Tau1 = new ArrayList<int[]>();
		BigInteger x;
		int numOfBitPermutations1 = (int) Math.pow(2, n1 - 1);							// number of permutations of n1-1 bits
		String binaryK = "";
		for(int k = 0; k < numOfBitPermutations1; k++)
		{
			sumAllBOverA = BigInteger.valueOf((long) modBOverACollection[0]);
			tupleForTau1 = new int[2];
			binaryK = Integer.toBinaryString(k);
			while(binaryK.length() < n1)
				binaryK = "0" + binaryK;



			//calculate partial sum
			for(int i = 1; i < n1; i++)
			{
				if(binaryK.charAt(i) == '0')
					sumAllBOverA = sumAllBOverA.add((BigInteger.valueOf((long) modBOverACollection[i])));
				else
					sumAllBOverA = sumAllBOverA.subtract((BigInteger.valueOf((long) modBOverACollection[i])));
			}
			x = ((negOne.multiply(sumAllBOverA)).add(modSOverA));
			rem_p = x.mod(bigP);	
			tupleForTau1[0] = rem_p.intValue();
			tupleForTau1[1] = k;
			Tau1.add(tupleForTau1);	
		}
		Tau1 = sortTupleByFirstElement(Tau1);

		/* Part (4b) of Alg 2.2 - assuming M > n, as in part (4a) above. */

		int[] tupleForTau2;
		Tau2 = new ArrayList<int[]>();
		int numOfBitPermutations2 = (int) Math.pow(2, n2);							// number of permutations of n1-1 bits
		String binaryI;	

		for(int i = 0; i < numOfBitPermutations2; i++)
		{
			sumAllBOverA = BigInteger.ZERO;	
			tupleForTau2 = new int[2];
			tupleForTau2[1] = i;
			binaryI = Integer.toBinaryString(i);
			while(binaryI.length() < n2)
				binaryI = "0" + binaryI;

			// Calculate partial sum
			for(int m = n1; m < n; m++)
			{
				if(binaryI.charAt(m-n1) == '0')
					sumAllBOverA = sumAllBOverA.add((BigInteger.valueOf((long) modBOverACollection[m])));
				else
					sumAllBOverA = sumAllBOverA.subtract((BigInteger.valueOf((long) modBOverACollection[m])));
					
			}
			rem_p = ((sumAllBOverA).mod(bigP));
			tupleForTau2[0] = rem_p.intValue();
			tupleForTau2[1] = i;
			Tau2.add(tupleForTau2);
		}
		Tau2 = sortTupleByFirstElement(Tau2);

//NOTE: comparting sieving arrays with               sdiff -s SIQSsievearry.txt KSIQSsievearry.txt | cat -n


		/* Part (5a) of Alg 2.2 */
		
		int k0 = 0;
		int maxi1 = (int) Math.pow(2, n1 - 1);
		int maxK = (int) Math.pow(2, n2);						
		int t1;
		int k1;
		for(int i1 = 0; i1 < maxi1; i1++)
		{
			t1 = (Tau1.get(i1))[0];
			while(k0 < maxK && (Tau2.get(k0))[0] < (t1 - M))			
				k0++;
			
			if(k0 > 0)
				k1 = k0-1;
			else
		       		k1 = k0;
			while((k1) < maxK && Tau2.get(k1)[0] < (t1 + M))
				k1++;	

			for(int i2 = k0; i2 < k1; i2++) 
			{
				int y = Tau1.get(i1)[0] - Tau2.get(i2)[0];
				int[] tau1 = Tau1.get(i1);
				int[] tau2 = Tau2.get(i2);
				int[] event = constructSievingEvent(p, tau1, tau2, y);
				sievingEvents.add(event);
			}	
		}
		

		/* Part (5b) of Alg 2.2 */
		
		k0 = 0;
		for(int i1 = 0; i1 < maxi1; i1++)
		{
			t1 = (Tau1.get(i1))[0];
			while(k0 < maxK && (Tau2.get(k0))[0] < (t1 + p - M))			
				k0++;
						
			for(int i2 = k0; i2 < maxK; i2++) 
			{
				int y = Tau1.get(i1)[0] - Tau2.get(i2)[0];
				if(y < -M || y > M){
					int[] tau1 = Tau1.get(i1);
					int[] tau2 = Tau2.get(i2);
					int[] event = constructSievingEvent(p, tau1, tau2, y);
					sievingEvents.add(event);	
				} // prevent 5b from detecting events from 5a in the event that M is such that there is overlap in areas detected from - prevent double detection
			}	
		}


		/* Part (5c) of Alg 2.2 */

		k1 = 0;
		for(int i1 = 0; i1 < maxi1; i1++)
		{
			t1 = (Tau1.get(i1))[0];
			while(k1 < maxK && (Tau2.get(k1+1))[0] < (t1 + M - p))			
				k1++;
						
			for(int i2 = 0; i2 <= k1; i2++) 
			{
				int y = Tau1.get(i1)[0] - Tau2.get(i2)[0];
				if(y < -M || y > M){
					if(y < -p || y > (-p + M)){
						int[] tau1 = Tau1.get(i1);
						int[] tau2 = Tau2.get(i2);
						int[] event = constructSievingEvent(p, tau1, tau2, y);
						sievingEvents.add(event);						
					}
				} // prevent 5b from detecting events from 5a in the event that M is such that there is overlap in areas detected from - prevent double detection

			}	
		}
		return sievingEvents;
	}





	// Helper Methods //

	/** Determines the value of the B term for polynomial j
	 *
	 *
	 *
	 */ 
	private static BigInteger determineValueOfB(int polynomialNumber, BigInteger[] B, int n){
		BigInteger b = B[0];

		String binaryPolynomialNumber = Integer.toBinaryString(polynomialNumber);
		while(binaryPolynomialNumber.length() != n)
			binaryPolynomialNumber = "0" + binaryPolynomialNumber;		
		for(int i = 1; i < binaryPolynomialNumber.length(); i++)
		{
			if(binaryPolynomialNumber.charAt(i) == '0')
				b = b.add(B[i]);
			else
				b = b.subtract(B[i]);

		}
		return b;
	}

	/** Determines all possible subsets (of fixed length) of a set of integers. 
	 * Subsets are then ordered based on the distance of product of elements from
	 * some target. 
	 *
	 * @param	factorBase is a set of primes forming the FB from a KSIQS operation
	 * @param	targetA the ideal A value for KSIQS
	 * @param	minFactor is the minimum size of a factor that can be included in a product
	 * @param	maxFactor is the maximum size of a factor that can be included in a product
	 * @param	n is the number of factors for A.
	 * @return	an ArrayList of int[] (all possible n-tuples) where product of elements can form A. 
	 *
	 * source: https://stackoverflow.com/questions/12548312/find-all-subsets-of-length-k-in-an-array
	 */
	private static ArrayList<int[]> determineAllPossibleA(ArrayList<Integer> factorBase, BigInteger targetA, int minFactor, int maxFactor, int n){
		int[] suitableFactors;
		int factorIndex = 0;

		while(factorBase.get(factorIndex) <= minFactor){				// Find index of first factor base prime > MIN_FACTOR_IN_A
			factorIndex++;
		}
		int maxFactorIndex = 0;
		while(factorBase.get(maxFactorIndex) <= maxFactor)
			maxFactorIndex++;

		suitableFactors = new int[maxFactorIndex - factorIndex];
		for(int i = 0; i < suitableFactors.length; i++)
			suitableFactors[i] = factorBase.get(factorIndex + i).intValue();


		ArrayList<int[]> subsets = new ArrayList<>();

		int[] s = new int[n];                  // here we'll keep indices 
		                                       // pointing to elements in input array
		if (n <= suitableFactors.length) {						 
    		// first index sequence: 0, 1, 2, ...
    			for (int i = 0; (s[i] = i) < n - 1; i++); 	       
    			subsets.add(getSubset(suitableFactors, s));
    			for(;;) {
        			int i;
        		// find position of item that can be incremented
				i = n - 1;
				while(i >= 0 && s[i] == (suitableFactors.length - n + i))
					i--;
			        if (i < 0) {
            				break;
        			}
        			s[i]++;                    // increment this item
        			for (++i; i < n; i++) {    // fill up remaining items
            				s[i] = s[i - 1] + 1; 
        			}
        			subsets.add(getSubset(suitableFactors, s));
    			}
		}
		subsets = sortByDistanceFromTarget(subsets, targetA);
		return subsets;
	}


	private static ArrayList<int[]> sortByDistanceFromTarget(ArrayList<int[]> toSort, BigInteger tgt){
		boolean change;
		BigInteger prod1;
		BigInteger prod2;
		ArrayList<int[]> toReturnSorted = toSort;
		do {
			change = false;
			for(int i = toReturnSorted.size(); i > 0; i--){
				
				prod1 = BigInteger.ONE;
				prod2 = BigInteger.ONE;
				for(int j = 0; j < toReturnSorted.get(i).length; j++){
					prod1 = prod1.multiply(BigInteger.valueOf((long)toReturnSorted.get(i)[j]));
					prod2 = prod2.multiply(BigInteger.valueOf((long)toReturnSorted.get(i-1)[j]));

				}
				


				if((prod1.subtract(tgt)).abs().compareTo((prod2.subtract(tgt)).abs()) < 0){
					int[] temp = toReturnSorted.get(i);
					toReturnSorted.set(i, toReturnSorted.get(i-1));
					toReturnSorted.set((i-1), temp);
					change = true;
				}

			}


		} while(change);

		return toReturnSorted;
	}


	// generate actual subset by index sequence
	private static int[] getSubset(int[] input, int[] subset) {
    		int[] result = new int[subset.length]; 
   		 for (int i = 0; i < subset.length; i++) 
        		result[i] = input[subset[i]];
    		return result;
	}



	private static BigInteger calculateC(BigInteger a, BigInteger b, BigInteger N){
		BigInteger c = ((b.pow(2)).subtract(N)).divide(a);
		return c;
	}



	private static int[] constructSievingEvent(int p, int[] tau1, int[] tau2, int y)
	{
		int n2 = (int) Math.floor(n/2.0);
		int jpart1 = tau1[1];
		int jpart2 = tau2[1];
		int j = jpart1 << n2;
		j += jpart2;

		int[] sieveEvent = {p, j, y};
		return sieveEvent;

	}

	/* Determine the factor base: the
	 * smooth primes such that N is a 
	 * quadratic residue (mod p). I.e. Legendre symbol
	 * (N/p) = 1
	 *
	 * @param F  the smoothness bound
	 * @param N  the composite being factored
	 * @return  a list containing the elements of the 
	 * factor base. 
	 */
	private static ArrayList<Integer> determineFactorBase(int F, BigInteger N) 
	{
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




	// TO profile in VisualVM, run as: java -Xverify:none KSIQS2017
	public static void main(String[] args)
	{
	//	BigInteger b1 = new BigInteger("416064700201658306196320137931");		// 30 digit prime
	//	BigInteger b2 = new BigInteger("513821217024129243948411056803");		// 30 digit prime
		BigInteger b1 = new BigInteger("1451730470513778492236629598992166035067");	// 40 digit prime
		BigInteger b2 = new BigInteger("2425967623052370772757633156976982469681");	// 40 digit prime
	//	BigInteger b1 = new BigInteger("103582180924623748121674293193717277486911");	// 42 digit prime
	//	BigInteger b1 = new BigInteger("123123412345123456123456712345678123456789");	// 42 digit prime

		BigInteger N = b1.multiply(b2);
		//60 digit N = 213782870618395542957440178483171817733685561437850435894593
		//79 digit N = 3521851118865011044136429217528930691441965435121409905222808922963363310303627
		//82 digit N = 12753391573589630958480469096673850446715878700449135789848997684642473014021588779	
		//
	//	getSmoothRelations(BigInteger.valueOf(101868649L), 233);	
	//	getSmoothRelations(N, 60000);
		getSmoothRelations(N, 900000);
	}//end main

/*
1041718812002259943259738289892712+19292445969865979347728497259531+1944329959097339189612111413389784+1861424791128709788609778912700832
 +527010423870396795419677806295364+757141965508564557951901379713951 +1482242903743321043940274962148793 +385910191782383402682251485098186 
 +2414506586176976738192155008422589 +779244978421970508782963427277404
*/
} // End KSIQS2017 class

