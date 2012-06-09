// random.cpp. R. Loos, 16. 11. 2000
#ifndef INFO_UT_HPP
#define INFO_UT_HPP
#define random_double()  ( (double)(++rand % 100000001) * 1e-8) /*[0,1]*/

namespace randlib {

//  Revision History
//   2005  Random double and random subarray generation (Sebastian Wernicke)
//   2000  Revision by Rüdiger Loos, C++ with 
//   1999  boost interface (Beman Dawes).
//   1996  Revision by Roland Weiss, C++ version
//   1982  Initial C version for Aldes 

//  (C) Copyright Rüdiger Loos 2000. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

//  (C) Copyright Beman Dawes 1994-99. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// rand -----------------------------------------------------------------------//

// a rand object initializes a pseudo random number generator

  class rand {
  public: 
    explicit rand(long seed = 310952) 
      : ran_arr_sentinel(-1), ran_arr_ptr(&ran_arr_sentinel) {
      //assert( 0 <= seed && seed < MM );
      ran_start(seed); 
    }
    rand& operator=( long new_value ) {
      //assert( 0 <= new_value && new_value < MM );
      ran_start(new_value); 
      return *this;
    }
    operator long() const       { return *ran_arr_ptr; }
    double fvalue() const       { return double(*ran_arr_ptr) / MM; }
    long operator++() {
      return *ran_arr_ptr>=0? *ran_arr_ptr++: ran_arr_cycle();
    }
    long operator++(int) { 
      long temp = *ran_arr_ptr; 
      operator++(); 
      return temp; 
    }
    long dice(long n) { return operator++() % n + 1; }
    bool trueWithProb(double d) { return ((operator++() & 0x0FFFFFFFL) <= (d*0x0FFFFFFFL)); }
    //  satisfy std::RandomNumberGenerator and std::Generator requirements:
    typedef long argument_type;
    typedef long result_type;
    long operator()( long n )   { return operator++() % n; }
    long operator()()           { return operator++(); }
    
  private:
    enum constants { 
      KK = 100,            // the long lag
      LL = 37,             // the short lag
      MM = (1L<<30),       // the modulus 
      TT = 70,             // guaranteed separation between streams
      QUALITY = 1009};     // recommended quality level for high-res use 
    long ran_x[KK];        // the generator state 
    long ran_arr_buf[QUALITY];
    long ran_arr_sentinel;
    long *ran_arr_ptr;     // the next random number, or -1
    bool is_odd(const int x) { return x & 1; }
    int evenize(const int x) { return x & (MM-2); }
    int mod_diff(const int x, const int y) { // subtraction mod MM  
      return (x - y) & (MM - 1);
    }
    void ran_start(long seed);        // ctor. seed selects different streams
    void ran_array(long aa[], int n); // put n new random numbers in ran_arr_buf 
                                      // array length (must be at least KK)
    long ran_arr_cycle();
  
  }; // class rand
} // namespace random

//    This program by D E Knuth is in the public domain and freely copyable
//    AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
//    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
//    (or in the errata to the 2nd edition --- see
//        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
//    in the changes to pages 171 and following of Volume 2).              */

//    If you find any bugs, please report them immediately to
//                 taocp@cs.stanford.edu
//    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
// moved to C++ by R. Loos, 12. 11. 2000, with ideas from Beman Dawes.


inline void gen_selection (long start, long end, 
	                       double prob, unsigned long *array, 
						   randlib::rand& rand) {
	long len = end-start;
	if (len == 0)
	{
		array[0] = end;
		return;
	}
	bool *sel = new bool[len];
	long point = 0;
	double rem = double(len) * prob;
	long num = (long)rem;
	rem -= num;
	if ( random_double() <= rem )
		{++num;}
	for (int i = 0; i!= len; ++i)
		{sel[i] = true;}
	if (prob < 0.6) {   // CASE OF SMALL PROBABILITIES
	    //Use Floyds random selection algorithm
		for (long J = len-num; J != len; ++J) {
			long pos = ++rand % (J+1);
			if (sel[pos]) {
				sel[pos] = false;
				array[point] = start + pos;
			} else {
				sel[J] = false;
				array[point] = start + J;
			}
			++point;
		}
		//while (num != 0) {
		//	long pos = ++rand % len;
		//	if (sel[pos]) {
		//		sel[pos] = false;
		//		array[point] = start + pos;
		//		//--num;
		//		//++point;
		//	}
		//}
	} else { 	// CASE OF LARGE PROBABILITIES
		num = len-num;
		for (long J = len-num; J != len; ++J) {
			long pos = ++rand % (J+1);
			if (sel[pos]) {
				sel[pos] = false;
				array[point] = start + pos;
			} else {
				sel[J] = false;
				array[point] = start + J;
			}
			++point;
		}
		for (int i = 0; i != len; ++i) {
			if (sel[i]) {
				array[point] = start + i;
				++point;
			}
		}
	}
	array[point] = end;
	delete[] sel;
	return;
}

#endif  // INFO_UT_HPP

