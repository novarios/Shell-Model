#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
#include <numeric>
#include <omp.h>
using namespace std;

double matrixe(double strength1, double strength2, vector<int> indvec, int bra, int ket, vector<double> onebody, vector<vector<int> > twobodybraket, vector<double> twobody)
{
  int length2 = int(twobody.size());
  int size = int(indvec.size());
  double temp = 0.0;
  vector<double> tempvec(length2);
  
  /*int m, n, l, k;
  int mbit, nbit, lbit, kbit;
  double phase;
  int tempket, bcount, comp, flag = 0;*/

  if (bra == ket)
    {
      for (int i = 0; i < size; ++i)
	{
	  if ((1 << i & ket) != 0){ temp = temp + onebody[indvec[i] - 1]; };
	};
    };

  std::cout << "<" << bra << "|" << ket << ">" << endl;

  /*for (int a = 0; a < length2; ++a)
    {
      m = twobodybraket[a][0];
      n = twobodybraket[a][1];
      l = twobodybraket[a][2];
      k = twobodybraket[a][3];
      mbit = 1 << (m - 1); nbit = 1 << (n - 1); lbit = 1 << (l - 1); kbit = 1 << (k - 1);
      phase = 1.0; tempket = ket;
      flag = 0;
      if ( ((((ket^lbit)^kbit)^nbit)^mbit) == bra )
	{ flag = 1; }
      else if( ((((ket^mbit)^nbit)^kbit)^lbit) == bra )
	{ flag = 1; swap(mbit,lbit); swap(m,l); swap(nbit,kbit); swap(n,k); }
      if (flag == 1)
	{
	  int ket2 = ((((ket^lbit)^kbit)^nbit)^mbit);
	  std::cout << ket2 << " " << bra << endl;
	  comp = tempket & ~(~0 << (l - 1));
	  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	  phase = phase*pow(-1, bcount); tempket = tempket^lbit;
	  comp = tempket & ~(~0 << (k - 1));
	  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	  phase = phase*pow(-1, bcount); tempket = tempket^kbit;
	  if ((tempket & nbit) == 0 && (tempket & mbit) == 0)
	    //{
	      comp = tempket & ~(~0 << (n - 1));
	      for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	      phase = phase*pow(-1, bcount); tempket = tempket^nbit;
	      comp = tempket & ~(~0 << (m - 1));
	      for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	      phase = phase*pow(-1, bcount); tempket = tempket^mbit;
	      //if (tempket == bra)
	      //{
		  temp = temp + phase*strength2*twobody[a];
		  //};
		  //};
	};
	};*/
  
  
  #pragma omp parallel for
  for (int a = 0; a < length2; ++a)
    {
      double temp2 = 0.0;
      int m = twobodybraket[a][0];
      int n = twobodybraket[a][1];
      int l = twobodybraket[a][2];
      int k = twobodybraket[a][3];
      int mbit = 1 << (m - 1); int nbit = 1 << (n - 1); int lbit = 1 << (l - 1); int kbit = 1 << (k - 1);
      double phase = 1.0; int tempket = ket;
      int flag = 0;
      if ( (lbit & ket) != 0 && (kbit & ket) != 0 && (nbit & bra) != 0 && (mbit & bra) != 0 &&
	   ((mbit & ket) == 0 || (mbit == lbit || mbit == kbit)) && ((nbit & ket) == 0 || (nbit == lbit || nbit == kbit)) )
	{ flag = 1; }
      else if( (mbit & ket) != 0 && (nbit & ket) != 0 && (kbit & bra) != 0 && (lbit & bra) != 0 &&
	       ((lbit & ket) == 0 || (lbit == mbit || lbit == nbit)) && ((kbit & ket) == 0 || (kbit == mbit || kbit == nbit)) )
	{ flag = 1; swap(mbit,lbit); swap(m,l); swap(nbit,kbit); swap(n,k); }
      if (flag == 1)
	{
	  int bcount;
	  int comp = tempket & ~(~0 << (l - 1));
	  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	  phase = phase*pow(-1, bcount); tempket = tempket^lbit;
	  comp = tempket & ~(~0 << (k - 1));
	  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	  phase = phase*pow(-1, bcount); tempket = tempket^kbit;
	  comp = tempket & ~(~0 << (n - 1));
	  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	  phase = phase*pow(-1, bcount); tempket = tempket^nbit;
	  comp = tempket & ~(~0 << (m - 1));
	  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
	  phase = phase*pow(-1, bcount); tempket = tempket^mbit;
	  if (tempket == bra)
	    {
	      temp2 = temp2 + phase*strength2*twobody[a];
	    };
	  std::cout << m << "," << n << "," << l << "," << k << "  " << phase << "," << twobody[a] << "," << temp2 << " " << endl;
	};
      tempvec[a] = temp2;
    };

  temp = std::accumulate(tempvec.begin(),tempvec.end(),temp);
  
  std::cout << temp << endl;
  
  return temp;
}
