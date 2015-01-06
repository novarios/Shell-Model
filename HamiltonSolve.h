#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

vector<double> solution(vector<double> H, size_t length)
{
  char jobz = 'V';
  char uplo = 'U';
  int n = int(length);
  vector<double> a = H;
  int lda = n;
  vector<double> w(n,0);
  int lwork = (3+2)*n;
  vector<double> work(lwork,0);
  int info;
  
  /*for (int i = 0; i < length; ++i)
    {
      for (int j = 0; j < length; ++j)
	{
	  std::cout << H[length * i + j] << " ";
	}
      std::cout << endl;
    }
    std::cout << endl;*/

  /*MatrixXd Ham(length, length);
  for (int i = 0; i < length; ++i)
    {
      for (int j = 0; j < length; ++j)
	{
	  Ham(i, j) = H[length*i + j];
	};
	};*/
  
  /*SelfAdjointEigenSolver<MatrixXd> eigensolver(Ham);
    if (eigensolver.info() != Success){ abort(); };*/
  
  dsyev_(&jobz, &uplo, &n, & *a.begin(), &lda, & *w.begin(), & *work.begin(), &lwork, &info);

  vector<double> solution(n * (n + 1));

  for (int p = 0; p < n; ++p)
    {
      solution[p * (n + 1)] = w[p];
      for (int i = 0; i < n; ++i)
	{
	  solution[p * (n + 1) + i + 1] = a[n * p + i];
	};
    };
  
  return solution;
}
