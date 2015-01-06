#include <iostream>
#include <math.h>
#include <vector>
#include <bitset>
using namespace std;

double StateP(vector<double> state, vector<int> configs, vector<int> levelsl)
{
  double parity = 1.0;
  int bra;
  
  for (int i = 0; i < int(state.size()); ++i)
    {
      bra = configs[i];
      
      for (int a = 0; a < int(levelsl.size()); ++a)
	{
	  if (((1 << a) & bra) != 0){ parity = parity * pow(-1.0, levelsl[a]); };
	};
      
    };
  
  return parity;
}
