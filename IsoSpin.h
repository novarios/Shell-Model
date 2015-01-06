#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
using namespace std;

double StateT(int count, vector<double> state, vector<int> configs, vector<double> angmomentum, vector<double> levelsm, vector<double> projections)
{
  double temp;
  int comp, bcount;
  int bra, ket;
  int mbit, nbit, lbit, kbit;
  int tempket; //temporary ket to act on
  double phase; //phase from applying operators
  double Tsquared = 0.0; //isospin squared
  double T; //return isospin
  
  for (int i = 0; i < int(state.size()); ++i)
    {
      bra = configs[i];
      
      for (int j = 0; j < int(state.size()); ++j)
	{
	  ket = configs[j];
	  
	  //T_z, T_z^2 part
	  if (bra == ket)
	    {
	      temp = 0.0;
	      for (int a = 0; a < int(projections.size()); ++a)
		{
		  if (((1 << a) & ket) != 0){ temp = temp + projections[a]; };
		};
	      
	      Tsquared = Tsquared + state[i] * state[j] * (temp*temp + temp);
	    };
	  
	  //T_-T_+ part
	  temp = 0.0;
	  for (int m = int(projections.size()) / 2; m < int(projections.size()); ++m) //don't include highest projection
	    {
	      for (int n = 0; n < int(projections.size()) / 2; ++n) //protons
		{
		  if (angmomentum[m] == angmomentum[n] && levelsm[m] == levelsm[n])
		    {
		      for (int l = 0; l < int(projections.size()) / 2; ++l) //don't include lowest projection
			{
			  for (int k = int(projections.size()) / 2; k < int(projections.size()); ++k) //don't include higher (or equal) projections
			    {
			      if (angmomentum[l] == angmomentum[k] && levelsm[l] == levelsm[k])
				{
				  mbit = 1 << m; nbit = 1 << n; lbit = 1 << l; kbit = 1 << k;
				  tempket = ket;
				  phase = 1.0;
				  if ((tempket & kbit) != 0)
				    {
				      comp = tempket & ~(~0 << k);
				      for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
				      phase = phase*pow(-1.0, bcount); tempket = tempket^kbit;
				      if ((tempket & lbit) == 0)
					{
					  comp = tempket & ~(~0 << l);
					  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
					  phase = phase*pow(-1.0, bcount); tempket = tempket^lbit;
					  if ((tempket & nbit) != 0)
					    {
					      comp = tempket & ~(~0 << n);
					      for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
					      phase = phase*pow(-1.0, bcount); tempket = tempket^nbit;
					      if ((tempket & mbit) == 0)
						{
						  comp = tempket & ~(~0 << m);
						  for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
						  phase = phase*pow(-1.0, bcount); tempket = tempket^mbit;
						  if (bra == tempket){
						    temp = temp + phase;
						  };
						};
					    };
					};
				    };
				};
			    };
			};
		    };
		};
	    };
	  
	  Tsquared = Tsquared + state[i] * state[j] * temp;
	  
	};
    };
  
  T = 0.5 * int(sqrt(1.0 + 4.0 * Tsquared) - 1.0 + 0.1);
  
  return T;
}
