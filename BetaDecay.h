#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
#include "CGC.h"
using namespace std;

double BetaMDecay(vector<double> statei, vector<double> statef, vector<int> configsi, vector<int> configsf, int indp, int indn, vector<int> levelsn, vector<int> levelsl, vector<double> levelsj, vector<double> levelsm)
{
  double temp;
  double M = 0.0;
  int comp, bcount;
  int bra, ket;
  int abit, bbit;
  int tempket; //temporary ket to act on
  double phase; //phase from applying operators
  double factor1, factor2, factor3, factor4, factor5; //factors for GT SPME
  
  for (unsigned int f = 0; f < statef.size(); ++f)
    {
      bra = configsf[f];
      for (unsigned int i = 0; i < statei.size(); ++i)
	{
	  ket = configsi[i];

	  if (bra == ket)
	    {
	      continue;
	    }
	  
	  temp = 0;
	  for (int a = indp; a < indp + indn; ++a) //neutron indicies
	    {
	      for (int b = 0; b < indp; ++b) //proton indicies
		{
		  abit = 1 << a;
		  bbit = 1 << b;
		  if ((ket & abit) != 0 && (ket & bbit) == 0 &&  levelsl[b] == levelsl[a] && levelsn[b] == levelsn[a])
		    {
		      tempket = ket;
		      phase = 1.0;
		      comp = tempket & ~(~0 << a);
		      for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };               //<b|a>
		      phase = phase*pow(-1.0, bcount);
		      tempket = tempket^abit;
		      comp = tempket & ~(~0 << b);
		      for (bcount = 0; comp; ++bcount){ comp ^= comp & -comp; };
		      phase = phase*pow(-1.0, bcount);
		      tempket = tempket^bbit;
		      if (bra == tempket)
			{
			  factor1 = sqrt(6.0*(2*levelsj[b] + 1.0)*(2*levelsj[a] + 1.0));
			  factor2 = pow(-1.0, levelsl[b] + int(levelsj[b] + 1.5));
			  factor3 = CGC6(0.5, 0.5, 1.0, levelsj[a], levelsj[b], levelsl[b]);
			  factor4 = pow(-1.0, int(levelsj[b] - levelsm[b]));
			  factor5 = CGC3(levelsj[b], -levelsm[b], 1.0, levelsm[b]-levelsm[a], levelsj[a], levelsm[a]);
			  temp = temp + phase * factor1 * factor2 * factor3 * factor4 * factor5;
			  /*if (abs(phase * factor1 * factor2 * factor3 * factor4 * factor5) > 0.001 && statef[f] != 0.0 && statei[i] != 0.0)
			    {
			      std::cout << "<" << (bitset<16>) bra << "|" << (bitset<16>) ket << "> : " << phase << ", " << factor1 * factor2 * factor3 * factor4 * factor5 << endl;
			      }*/
			  //std::cout << "bbit = " << (bitset<16>) bbit << ", abit = " << (bitset<16>) abit;
			  //std::cout << levelsj[b] << ", " << levelsm[b] << ", " << levelsj[a] << ", " << levelsm[a] << endl;
			  //std::cout << factor1 << " " << factor2 << " " << factor3 << " " << factor4 << " " << factor5 << endl << endl;
			  //std::cout << phase << endl;
			};
		    };
		};
	    };
	  /*if (abs(temp) > 0.001 && statef[f] != 0.0 && statei[i] != 0.0)
	    {
	      std::cout << statei[i] << " " << statef[f] << " " << statef[f] * statei[i] * temp << endl;
	      }*/
	  M = M + statef[f] * statei[i] * temp;
	};
    };
  return M;
}
  
