#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

vector<int> bitconfigsetup(vector<int> newconfigs, int N)
{
  int blength1 = int(newconfigs.size()) / N;
  vector<int> bconfigs(blength1);
  int btemp;
  int bocc;
  int shift;
  int bsd;
  
  for (int i = 0; i < blength1; i++)
    {
      bsd = 0;
      for (int j = 0; j < N; j++)
	{
	  bocc = newconfigs[N*i + j]; shift = bocc - 1; btemp = 1 << shift; bsd = bsd + btemp;
	};
      bconfigs[i] = bsd;
    };
  return bconfigs;
}
