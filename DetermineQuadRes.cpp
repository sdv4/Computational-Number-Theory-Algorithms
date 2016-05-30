#include <iostream>
using namespace std;

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
	 * @author  Shane Sims
	 * @version 30 May 2016
	 */
	int LegJacSym(int a, int m) {
		a = a%m;
		int t = 1;
		while(a != 0) {
			while(a % 2 == 0) {
				a = a/2;
				if(m % 8 == 3 || m % 8 == 5)
					t = t * -1;
			}
			int temp = a;
			a = m;
			m = temp;
			if(a % 4 == m && m % 4 == 3)
				t = t * -1;
			a = a % m;
		}
		if (m == 1) 
			return t;
		return 0;

	}

	int main() {
		int a, m;
		cout <<"To calculate the Legendre/Jacobi symbol (a/m), \n";
		cout <<" enter a: ";
		cin >> a;
		cout <<"\n enter m: ";
		cin >> m;
		cout << "\n You entered: \n";
		cout << " a = " << a << "\n";
		cout << " m = " << m << "\n\n";

		cout << "(" << a << "/" << m <<") = ";
		cout << LegJacSym(a,m) << "\n";

	}
