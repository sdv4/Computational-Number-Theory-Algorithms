import java.lang.Math.*;
import java.util.*;
import java.math.BigInteger;
import java.io.*;
import java.util.Scanner;

/** This is an implementation of the Sieving portion of the Self Initializing 
* Quadratic Sieve factoring algorithm. Built on code written in authours
* MPQS implementation and based on the pseduocode provided in 
* the thesis 'Factoring Integers with the Self-Initializing Quadratic Sieve' by
* Scott P. Contini. The notation matches
* that found in this thesis as far as possible so as to increase the didactic 
* value of this class. 
*
* This version is for numbers N of approximately 60 digits.
*
* TODO: Split sieving array into blocks for sieving
*
* @author Shane Sims <shane.sims.ss@gmail.com>
* @version 1 August 2017 
*
*/
public class SIQSsixty {	

	public static final int ERROR_TERM = 27;
	public static final int M = 360000;
	public static final int DESIRED_NUM_RELATIONS = 30;
	public static final int MIN_FACTOR_IN_A = 2000;						// 170th prime in FB as per Contini reccomendation
	public static final int n = 10;								// desired number of factors of term 'a'
//	public static final int MAX_FACTOR_IN_A = 3965;						// Hard coded value of (sqrt(N)/M)^(1/7) - needs to change if M, N or n change. N 60 dig, n = 7
	public static final int MAX_FACTOR_IN_A = 2353;						// Hard coded value of (sqrt(N)/M)^(1/7) - needs to change if M, N or n change. N 60 dig, n = 7


	/**
     	 * Hide utility class.
     	 */
    	private SIQSsixty() {}	


	/**
	* This method determines the smooth relations required to factor N.
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
	public static BigInteger[][] getSmoothRelations(BigInteger N, int F){

//		Scanner scanner = new Scanner(System.in);
//    		System.out.println("Press Enter to begin\t");
//        	scanner.nextLine();
		
		BigInteger[][] S;								// Array that will hold the smooth relations when returned.
		ArrayList<Integer> factorBase;							// Primes which will be used in the sieving stage.
		int[] tmem_p;									// Will hold t in t^2 \equiv N (mod p) for each p in Factor Base
		int[] solution1;								// Will hold a^-1(t - b) \equiv N (mod p) for each prime in factor base
		int[] solution2;								// Will hold a^-1(-t - b) \equiv N (mod p) for each prime in factor base
		int[] logp;									// Approximate logarithm for each prime in the factor base
		int K;										// Method will return when K+1 smooth relations have been found.
		int offset = 0;									// Used if/when sieving stage returned to
		BigInteger a;									// Term 'a' in g_a,b(x) = (ax+b)^2 - N 
		BigInteger b;									// Term 'b' in g_a,b(x)
		BigInteger c;									// Term 'c' in g_a,b(x) = (ax+b)^2 - N = a(ax^2 + 2bx + c)
		ArrayList<Integer> aFactors;							// The prime factors, the product of which forms a
		int[][] Bainv2_jp;								// Precomputation term for time saving on initialization. j = polynomial id; p = FB prime
		BigInteger aInv_p[];								// a^-1 (mod p) for each p in FB
	
	//TODO: check preconditions of N are met here
	
	/******************** Compute Startup Data **********************************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/


		System.out.println("Computing Startup Data...");

		int polys = 0;									// Number of polynomials initialized so far
		BigInteger squareRootN = BigIntegerSqrt.bigIntegerSqrtCeiling(N);		
		factorBase = determineFactorBase(F, N);						// Initialize factorBase.
		K = factorBase.size();								// Save factor base size. Need K+1 smooth valued poly to ensure factorization
		S = new BigInteger[K+1][4];							// Elements are pairs (x,(x+b)^2-N)				
		logp = new int[K];								
		tmem_p = new int[K];								// Initialize modSqrts to have space for each element of factorBase
		tmem_p[0] = 1;									// TonelliSqrtModP requires odd prime, so only even prime hard coded
		solution1 = new int[K];
		solution2 = new int[K];
		aInv_p = new BigInteger[K];
		b = BigInteger.ZERO;

		for(int i = 1; i < tmem_p.length; i++) {					// Solve t in t^2 \equiv N (mod p) for each p in factorBase
			tmem_p[i] = TonelliSqrtModP.SqrtModP(N, factorBase.get(i));
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

	int smoothRelFound = 0;															// Tracks the number of smooth relations verified by trial division.

	//Assumption: will get required number of relations with one value of a
	// TODO: Implement loop if above assumption doesn't hold
	
		a = BigInteger.valueOf(1L);													// Will want a approx sqrt(2*N)/M
		BigInteger targetA = BigInteger.valueOf(0L);	
		BigInteger two = BigInteger.valueOf(2L);
		BigInteger negOne = BigInteger.valueOf(-1L);
		BigInteger bigM = BigInteger.valueOf((long) M);	
		targetA = N.multiply(two);						
		targetA = BigIntegerSqrt.bigIntegerSqrtCeiling(targetA);
		targetA = targetA.divide(bigM);												// Target value of a now set
		aFactors = new ArrayList<Integer>();
		int factorIndex = 0;
		//while(factorBase.get(factorIndex) <= MIN__FACTOR_IN_A){
		while(factorBase.get(factorIndex) <= MAX_FACTOR_IN_A){
			factorIndex++;
		}

		BigInteger tempFactor = BigInteger.valueOf((long) factorBase.get(factorIndex));
		BigInteger tempA = tempFactor;							// Temp. value of a to ensure size is as close to target w/o going over - to keep a as small as possible
		int aNumFactors = 0;												// Track number of factors in a, since .size() not necessarily accurate.

		// Set a	
		while(aNumFactors < n && factorIndex < factorBase.size()){			
			a = tempA;
			aFactors.add(tempFactor.intValue());
			aNumFactors++;
			//factorIndex++;
			factorIndex--;
			tempFactor = BigInteger.valueOf((long) factorBase.get(factorIndex));	// Get value in factorbase at index 'factorIndex'
			tempA = tempA.multiply(tempFactor);							
		}
		





System.out.println("number of factors in a: " + aFactors.size());
		int numOfPossiblePolys = (int) Math.pow(2, (aFactors.size() - 1));		// Number of polynomials that can be generated from a

System.out.println("a: " + a);

System.out.println("Number of possible b's for this a: " + numOfPossiblePolys);
		Bainv2_jp = new int[numOfPossiblePolys][factorBase.size()];			
		BigInteger[] B = new BigInteger[aNumFactors];
		
		while(polys < numOfPossiblePolys && smoothRelFound < DESIRED_NUM_RELATIONS){
			/* Initialization for first polynomial */
			if(polys == 0){
		
				// Get B_i for all i factors of a		
				BigInteger gamma;
				BigInteger bigTmemp;
				BigInteger tempTerm2;
				BigInteger bigAFactor;
		
				for(int i = 0; i < aNumFactors; i++){
					tempTerm2 = BigInteger.ONE;
					bigAFactor = BigInteger.valueOf((long) aFactors.get(i));
					for(int j = 0; j < aFactors.size(); j++){
						if(i != j)
							tempTerm2 = tempTerm2.multiply(BigInteger.valueOf((long) aFactors.get(j)));			// compute a/bigAFactor
					}					
					bigTmemp = BigInteger.valueOf((long) tmem_p[factorBase.indexOf(aFactors.get(i))]);				// get tmem_p, where p is current factor of 'a'
					tempTerm2 = tempTerm2.modInverse(bigAFactor);
					gamma = (bigTmemp.multiply(tempTerm2)).mod(bigAFactor);
					float gammaFloat = gamma.floatValue();
					float qover2 = bigAFactor.floatValue()/2;
						if(gammaFloat > qover2)
							gamma = bigAFactor.subtract(gamma);
					B[i] = gamma.multiply(a.divide(bigAFactor));
//System.out.println("B[" + i + "] = " + B[i]);					
				}
//		Scanner scanner = new Scanner(System.in);
//     		System.out.println("Press Enter to begin\t");
//        	scanner.nextLine();

				// Set first value of b and first solution1 and solution2
		//		for(int i = 0; i < B.length; i++)
		//			b = b.add(B[i]);

			b = determineValueOfB(polys, B, n);
				
//				System.out.println("b = " + b);
//				Scanner scanner = new Scanner(System.in);
//		     		System.out.println("Press Enter to begin\t");
//		        	scanner.nextLine();


				BigInteger bigaInv_p;
				BigInteger fb_p;
				BigInteger bigTemp;
				BigInteger bigtmem_p;

				for(int i = 0; i < factorBase.size(); i++){				// For each p in FB such that 
					int p = factorBase.get(i);
					if(!aFactors.contains(p)){					// NOT p|a
						fb_p = BigInteger.valueOf((long) p);			
						bigaInv_p = a.modInverse(fb_p);				// compute a^-1 (mod p).
						aInv_p[i] = bigaInv_p;

						for(int j = 0; j < aFactors.size(); j++){		// For each factor of 'a',
							Bainv2_jp[j][i] = (((two.multiply(B[j])).multiply(bigaInv_p)).mod(fb_p)).intValue();	// compute 2 * B_j * a^(-1) (mod p)
						}
						bigtmem_p = BigInteger.valueOf((long) tmem_p[i]);
						solution1[i] = ((bigaInv_p.multiply(bigtmem_p.subtract(b))).mod(fb_p)).intValue();
						solution2[i] = ((bigaInv_p.multiply((BigInteger.ZERO.subtract(bigtmem_p)).subtract(b))).mod(fb_p)).intValue();

					}//end if
					else{
						solution1[i] = -1;
						solution2[i] = -1;							// Mark unusable solutions.
					}

				}//end for loop precomputing values for subsequent poly init

			

System.out.println("a: " + a + "       b: " + b);
System.out.println("If a and b are correct, this should be 0:  " + (b.multiply(b).mod(a)).subtract(N.mod(a)));
				c = calculateC(a, b, N);
				polys++;										// Increment number of polynomials initialized
			}//end initialization of first poly
		

			/* Initialization for subsequent polynomials */
			// Consists of getting new 'b' value and computing new solution1 and solution2
			else{
				// Switch b										// See Contini thesis Appendix C for Grey code information
		//		int ii = polys;										
		//		int v = 1;										// Will hold position of right most set bit. Init to 1 to check if LSB set
		//		while((ii & 1) != 1){
		//			ii = ii >>>1;
		//			v++;
		//		}

		//		int greyCodeExp = (int) Math.ceil(polys / (Math.pow(2, v)));
		//		BigInteger negOneToPowGrey = negOne.pow(greyCodeExp);

		//		BigInteger tempppp = negOneToPowGrey.multiply(two);
			b = determineValueOfB(polys, B, n);
			//	b = b.add(tempppp.multiply(B[v]));	

				// Compute solution1 and solution2 for this polynomial

		
////////////////////////////////////////////////////////////////////////////////////////////////////// need to get this working to realize savings in changing solutions 1 and 2
			/*	BigInteger fb_p;
				BigInteger BaInv_jp;							
				for(int i = 0; i < factorBase.size(); i++){				// For each p in FB such that 
					int p = factorBase.get(i);
					if(!aFactors.contains(p)){					// NOT p|a
						fb_p = BigInteger.valueOf((long) p);
						BaInv_jp = BigInteger.valueOf((long) Bainv2_jp[v][i]);
						solution1[i] += ((negOneToPowGrey.multiply(BaInv_jp)).mod(fb_p)).intValue();
						solution2[i] += ((negOneToPowGrey.multiply(BaInv_jp)).mod(fb_p)).intValue();
					}

				}		
			*/
/////////////////////////////////////////////////////////////////////////////////////////////////////

	
				int p;	
				BigInteger bigtmem_p;
				BigInteger fb_p;
				for(int i = 0; i < factorBase.size(); i++){
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
				polys++;										// Increment number of polynomials initialized

			}//end else/init subsequent polys

		System.out.println("Initialization Stage complete...");

	/**************************** Sieving Stage *********************************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/

		System.out.println("Starting Sieving Stage...");

		int index1 = 0;									// Will hold latest value of (solution1_p + ip) or (solution2_p + ip)
		int[] sieveArray = new int[2*M + 1];						// Initialize sieve array with 0's.
		for(int j = 0; j < factorBase.size(); j++) {					// For each prime in the factor base.
			int p = factorBase.get(j);
			int s1 = solution1[j];
			int s2 = solution2[j];


			if(s1 != -1){
				int i = -1;
				//find lower bound value for i such that soln1 +i*p >= -M
				while((s1 + (i*p)) >= -M && (s1 + (i*p)) <= M){
					index1 = s1 + (i*p);	
					sieveArray[index1 + M] += logp[j];			// add M to offset sieve array starts at -M = sieveArray[0]
					i--;
				}
				i = 0;
				//find upper bound value for i such that soln1 +i*p <= M					
				while((s1 + (i*p)) >= -M && (s1 + (i*p)) <= M){						
					index1 = s1 + (i*p);
					sieveArray[index1 + M] += logp[j];
					i++;
				}
				if(p != 2){							// If p = 2, sieve only with solution1
					i = -1;
					while((s2 + (i*p)) >= -M && (s2 + (i*p)) <= M){
						index1 = s2 + (i*p);
						sieveArray[index1 + M] += logp[j];
						i--;
					}
					i = 0;
					while((s2 + (i*p)) >= -M && (s2 + (i*p)) <= M){
						index1 = s2 + i*p;
						sieveArray[index1 + M] += logp[j];
						i++;
					}
				}
			}//end if

		}//end for
		
System.out.println("Sieving Stage complete...");

	/********************************** Trial Division Stage ********************************/
	/*                                                                                      */
	/*                                                                                      */
	/*                                                                                      */
	/****************************************************************************************/
/*
if(b.compareTo(new BigInteger("6734053319198016744014228498473658")) == 0)
{
System.out.println("b: " + b);

	try{
    PrintWriter writer = new PrintWriter("SIQSsievearry.txt", "UTF-8");
    int line = 0;
    for(int elmnt : sieveArray){
    	writer.println("index" + line + ": " + elmnt);
	line++;
    }
    writer.close();
} catch (IOException e) {
   // do something
}
		Scanner scanner = new Scanner(System.in);
     		System.out.println("Press Enter to begin\t");
        	scanner.nextLine();



}
*/

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
			for(int j = 0; j < factorBase.size(); j++){					// For each prime in the factor base
				int p = factorBase.get(j);
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
		for(int j = 0; j < S.length; j++){
			System.out.printf(" %-12d %-12d %-5d %d \n", S[j][0], S[j][1], S[j][2], S[j][3]);
			if(S[j][0] == null)
				break;
		}

	}//end outer while smooth rel < needed	


System.out.println("Polynomials initialized: " + polys);
System.out.println("Size of FB: " + factorBase.size());



	return S;

}//end get smooth relations 




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
//		System.out.println("b_0 = " + b);
//		Scanner scanner = new Scanner(System.in);
//System.out.println("Press Enter to begin\t");
//scanner.nextLine();

		return b;



	}



private static BigInteger calculateC(BigInteger a, BigInteger b, BigInteger N){
	BigInteger c = ((b.pow(2)).subtract(N)).divide(a);
	return c;

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
private static BigInteger[] determineAllBValues(int maxNumPolys, BigInteger[] B, ArrayList<int[]> epsilon)
	{

		BigInteger[] allBValues = new BigInteger[maxNumPolys];		
	
		for(int i = 0; i < allBValues.length; i++){
			int[] delta = epsilon.get(i);
			BigInteger bValue = BigInteger.ZERO;
			BigInteger termForB; 
			for(int j = 0; j < delta.length; j++){
				if(delta[j] == 1)
					bValue = bValue.add(B[j]);
				else
					bValue = bValue.subtract(B[j]);

			}
			allBValues[i] = bValue;
		}
		return allBValues;
	}



	public static void main(String[] args){
	//	BigInteger b1 = new BigInteger("416064700201658306196320137931");		// 30 digit prime
	//	BigInteger b2 = new BigInteger("513821217024129243948411056803");		// 30 digit prime
		BigInteger b1 = new BigInteger("1451730470513778492236629598992166035067");	// 40 digit prime
		BigInteger b2 = new BigInteger("2425967623052370772757633156976982469681");	// 40 digit prime
	//	BigInteger b1 = new BigInteger("103582180924623748121674293193717277486911");	// 42 digit prime
	//	BigInteger b1 = new BigInteger("123123412345123456123456712345678123456789");	// 42 digit prime
		
		BigInteger N = b1.multiply(b2);
		//N = 213782870618395542957440178483171817733685561437850435894593				// 60 digit prime
		//N = 3521851118865011044136429217528930691441965435121409905222808922963363310303627		// 79 digit prime
		//N = 12753391573589630958480469096673850446715878700449135789848997684642473014021588779	// 82 digit prime
	//	getSmoothRelations(BigInteger.valueOf(101868649L), 233);	
	//	getSmoothRelations(N, 60000);
		getSmoothRelations(N, 900000);
	}//end main

}//end BasicQS

