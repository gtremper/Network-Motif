// random.cpp. For specification and copyright see random.hpp. R. Loos, 16. 11. 2000

#include "random.hpp"

namespace randlib {


//  rand class implementation  ----------------------------------------------//

  void randlib::rand::ran_array(long aa[],int n) 
  {
    //assert(n >= KK);
    int i,j;
    for (j=0;j<KK;j++) 
      aa[j]=ran_x[j];
    for (;j<n;j++) 
      aa[j]=mod_diff(aa[j-KK],aa[j-LL]);
    for (i=0;i<LL;i++,j++) 
      ran_x[i]=mod_diff(aa[j-KK],aa[j-LL]);
    for (;i<KK;i++,j++) 
     ran_x[i]=mod_diff(aa[j-KK],ran_x[i-LL]);
  }
  
  void randlib::rand::ran_start(long seed) 
  {
    int t,j;
    long x[KK+KK-1];              /* the preparation buffer */
    long ss=evenize(seed+2);
    for (j=0;j<KK;j++) {
      x[j]=ss;                      /* bootstrap the buffer */
      ss<<=1; if (ss>=MM) ss-=MM-2; /* cyclic shift 29 bits */
    }
    for (;j<KK+KK-1;j++) x[j]=0;
    x[1]++;              /* make x[1] (and only x[1]) odd */
    ss=seed&(MM-1);
    t=TT-1; while (t) {
      for (j=KK-1;j>0;j--) x[j+j]=x[j];  /* "square" */
      for (j=KK+KK-2;j>KK-LL;j-=2) x[KK+KK-1-j]=evenize(x[j]);
      for (j=KK+KK-2;j>=KK;j--) if(is_odd(x[j])) {
        x[j-(KK-LL)]=mod_diff(x[j-(KK-LL)],x[j]);
        x[j-KK]=mod_diff(x[j-KK],x[j]);
      }
      if (is_odd(ss)) {              /* "multiply by z" */
        for (j=KK;j>0;j--)  x[j]=x[j-1];
        x[0]=x[KK];            /* shift the buffer cyclically */
        if (is_odd(x[KK])) x[LL]=mod_diff(x[LL],x[KK]);
      }
      if (ss) ss>>=1; else t--;
    }
    for (j=0;j<LL;j++) ran_x[j+KK-LL]=x[j];
    for (;j<KK;j++) ran_x[j-LL]=x[j];
  }
  
  /* the following routines are from exercise 3.6--15 */
  
  long randlib::rand::ran_arr_cycle()
  {
    ran_array(ran_arr_buf,QUALITY);
    ran_arr_buf[100]=-1;
    ran_arr_ptr=ran_arr_buf+1;
    return ran_arr_buf[0];
  }
  
} // namespace randlib
