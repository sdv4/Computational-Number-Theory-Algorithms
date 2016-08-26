import java.lang.Math;
/* 
 * @file ModularInverseCalc.java
 * @author Shane Sims, <shane.sims@ucalgary.ca>
 * @version  26 August 2016
 * 
 * This class contains method to calculate the modular inverse 
 * of an integer. Will be updated as additional methods are
 * needed. 
*/

public class ModularInverseCalc {
	
	/* Method to calculate the integer congruence a^-1 of a mod p
	 * when a and p are coprime (gcd(a, p) =1). This method uses 
	 * the Euler's Toitent Function technique of determining 
  	 * modular invverse. 
	 *
 	 * @param a  is an integer whose inverse will be determined
	 * @param p  is the modulus used in the calculation, coprime 
	 * with a. 
	 * @return aInv the integer congruent to a^-1 (mod p)
	 */
	public static double calcModInverse1(int a, double p) {
		double phiP = (1 - (1/p)) * p;
		double aInv = Math.pow(a, phiP-1);
		aInv = aInv % p;
		return aInv;
	}

	/* Method to calculate the integer congruence a^-1 of a mod p
	 * when a and p are coprime (gcd(a, p) =1). This method uses 
	 * the Euler's Toitent Function technique of determining 
  	 * modular invverse in the special case where p is prime.
	 *
 	 * @param a  is an integer whose inverse will be determined
	 * @param p  the prime modulus used in the calculation.
	 * @return aInv the integer congruent to a^-1 (mod p)
	 */
	public static double calcModInverse2(int a, double p) {
		double aInv = Math.pow(a, p-2);
		aInv = aInv % p;
		return aInv;
	}

}
