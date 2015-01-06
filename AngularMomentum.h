#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
using namespace std;

double StateJ(int count, vector<double> state, vector<int> configs, vector<int> indvec, vector<double> angularmom, vector<double> projections, vector<vector<int> > shellprojections)
{
  double temp;
  int comp, bcount;
  int bra, ket;
  int m, n, l, k;
  int mbit, nbit, lbit, kbit;
  int tempket; //temporary ket to act on
  double phase; //phase from applying operators
  double Jsquared = 0.0; //isospin squared
  double J; //return isospin
  vector<int> indicies1, indicies2; //vectors of relevant shells
  
  for (int i = 0; i < int(state.size()); ++i)
    {
      indicies1.clear();
      bra = configs[i];
      for (int a = 0; a < int(shellprojections.size()); ++a)
	{
	  for (int b = 0; b < int(shellprojections[a].size()); ++b)
	    {
	      if (((1 << (shellprojections[a][b] - 1)) & bra) != 0){ indicies1.push_back(a); break; };
	    };
	};
      
      for (int j = 0; j < int(state.size()); ++j)
	{
	  indicies2.clear();
	  ket = configs[j];
	  for (int a = 0; a < int(shellprojections.size()); ++a)
	    {
	      for (int b = 0; b < int(shellprojections[a].size()); ++b)
		{
		  if (((1 << (shellprojections[a][b] - 1)) & ket) != 0)
		    {
		      for (int c = 0; c < int(indicies1.size()); ++c)
			{ 
			  if (indicies1[c] == a){ indicies2.push_back(a); break; }; 
			}; break;
		    };
		};
	    };
	  
	  //J_z, J_z^2 part
	  if (bra == ket)
	    {
	      temp = 0;
	      for (int a = 0; a < int(angularmom.size()); ++a)
		{
		  if ((1 << a & ket) != 0){ temp = temp + projections[a]; };
		};
	      Jsquared = Jsquared + state[i] * state[j] * (temp*temp + temp);
	    };

	  //J_-J_+ part
	  temp = 0;
	  for (int a = 0; a < int(indicies2.size()); ++a)
	    {
	      for (int a1 = 0; a1 < int(shellprojections[indicies2[a]].size()) - 1; ++a1) //don't include highest projection
		{
		  for (int a2 = a1 + 1; a2 < int(shellprojections[indicies2[a]].size()); ++a2) //don't include lower (or equal) projections
		    {
		      m = shellprojections[indicies2[a]][a1] - 1, n = shellprojections[indicies2[a]][a2] - 1;
		      if (projections[m] != projections[n] - 1.0){ continue; }
		      else{
			for (int b = 0; b < int(indicies2.size()); ++b)
			  {
			    for (int b1 = 1; b1 < int(shellprojections[indicies2[b]].size()); ++b1) //don't include lowest projection
			      {
				for (int b2 = 0; b2 < b1; ++b2) //don't include higher (or equal) projections
				  {
				    l = shellprojections[indicies2[b]][b1] - 1, k = shellprojections[indicies2[b]][b2] - 1;
				    if (projections[l] != projections[k] + 1.0){ continue; }
				    else{
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
						      if (bra == tempket)
							{
							  temp += phase*sqrt(angularmom[n]*(angularmom[n]+1.0)-projections[n]*(projections[n]-1.0))*sqrt(angularmom[k]*(angularmom[k]+1.0)-projections[k]*(projections[k]+1.0));
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
		};
	    };
	  Jsquared = Jsquared + state[i] * state[j] * temp;
	};
    };

  J = 0.5 * (sqrt(1.0 + 4.0 * Jsquared) - 1.0);
  
  return J;
}
