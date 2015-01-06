#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

vector<vector<double> > Lanczos(vector<double> Hamiltonian)
{
  vector<vector<double> > alphabeta;
  
  int veclength = int(sqrt(Hamiltonian.size()));
  vector<vector<double> > V;
  vector<double> v0(veclength, 0.0), v1;
  V.push_back(v0);
  
  vector<double> alpha, beta;
  alpha.push_back(0.0);
  beta.push_back(0.0);
  beta.push_back(0.0);
  
  double normal = 0.0;
  
  vector<vector<double> > Omega;
  Omega.push_back(v0);
  
  //build random vector v1
  for (int i = 0; i < veclength; ++i)
    {
      int random = rand();
      v1.push_back(random);
      normal = normal + random*random; //keep track of vector magnitude
    }
  
  //normalize vector
  for (int i = 0; i < int(v1.size()); ++i)
    {
      v1[i] = v1[i] / sqrt(normal);
    }
  
  V.push_back(v1);
  
  for (int i = 1; i < veclength; ++i)
    {
      vector<double> Hdot; //H(dot)vi
      for (int j = 0; j < veclength; ++j)
	{
	  double Hdotrow = 0.0;
	  for (int h = 0; h < veclength; ++h)
	    {
	      Hdotrow = Hdotrow + V[i][h] * Hamiltonian[veclength*j + h];
	    }
	  Hdot.push_back(Hdotrow);
	}
      
      Omega.push_back(Hdot);
      
      double alphatemp = 0.0;
      for (int j = 0; j < veclength; ++j)
	{
	  alphatemp = alphatemp + Omega[i][j] * V[i][j];
	}
      
      alpha.push_back(alphatemp);
      
      vector<double> newomega;
      for (int j = 0; j < veclength; ++j)
	{
	  newomega.push_back(Omega[i][j] - alpha[i] * V[i][j] - beta[i] * V[i - 1][j]);
	}
      
      Omega[i] = newomega;
      
      double tempnorm = 0.0;
      for (int j = 0; j < veclength; ++j)
	{
	  tempnorm = tempnorm + Omega[i][j] * Omega[i][j];
	}
      
      beta.push_back(sqrt(tempnorm));
      
      vector<double> tempv;
      for (int j = 0; j < veclength; ++j)
	{
	  tempv.push_back(Omega[i][j] / beta[i + 1]);
	}
      
      V.push_back(tempv);
    }
  
  vector<double> Hdot; //H(dot)vm
  for (int j = 0; j < veclength; ++j)
    {
      double Hdotrow = 0.0;
      for (int h = 0; h < veclength; ++h)
	{
	  Hdotrow = Hdotrow + V[veclength][h] * Hamiltonian[veclength*j + h];
	}
      Hdot.push_back(Hdotrow);
    }
  
  Omega.push_back(Hdot);
  
  double alphatemp = 0.0;
  for (int j = 0; j < veclength; ++j)
    {
      alphatemp = alphatemp + Omega[veclength][j] * V[veclength][j];
    }
  
  alpha.push_back(alphatemp);
  
  //get rid of zeros at beginning of vectors
  alpha.erase(alpha.begin());
  beta.erase(beta.begin());
  beta.erase(beta.begin());

  //erase v0 to get matrix of Lanczos vectors
  V.erase(V.begin());

  //get V into single vector form (columns)
  vector<double> Vm;
  for (int i = 0; i < int(V.size()); ++i)
    {for (int j = 0; j < int(V[i].size()); ++j)
	{
	  Vm.push_back(V[j][i]);
	}
    }

  alphabeta.push_back(alpha);
  alphabeta.push_back(beta);
  alphabeta.push_back(Vm);
  
  return alphabeta;
}



vector<vector<double> > Householder(vector<double> Hamiltonian)
{
  vector<vector<double> > alphabeta;
  vector<double> alpha, beta;
  vector<double> A = Hamiltonian;
  
  int veclength = int(sqrt(Hamiltonian.size()));
  
  vector<double> S;
  for (int i = 0; i < veclength; ++i)
    {
      for (int j = 0; j < veclength; ++j)
	{
	  if (i == j){ S.push_back(1.0); }
	  else{ S.push_back(0.0); }
	}
    }
  
  for (int n = 0; n < veclength - 2; ++n)
    {
      vector<double> V;
      double Vnorm2 = 0.0;
      for (int i = 0; i < veclength - n - 1; ++i)
	{
	  V.push_back(A[veclength*n + (n + 1) + i]);
	  Vnorm2 = Vnorm2 + V[i] * V[i];
	}
      
      double sign = 1.0;
      if (V[0] < 0){ sign = -1.0; }
      
      double Vnorm = -sign * sqrt(Vnorm2);
      double r = sqrt((Vnorm2 - V[0] * Vnorm) / 2.0);
      
      vector<double> U;
      for (int i = 0; i < int(V.size()); ++i)
	{
	  if (i == 0){ U.push_back((V[i] - Vnorm) / (2.0*r)); }
	  else{ U.push_back(V[i] / (2.0*r)); }
	}
      
      vector<double> Stemp;
      for (int i = 0; i < veclength; ++i)
	{
	  for (int j = 0; j < veclength; ++j)
	    {
	      if (i == j && i <= n){ Stemp.push_back(1.0); }
	      else if (i != j && (i <= n || j <= n)){ Stemp.push_back(0.0); }
	      else if (i == j){ Stemp.push_back(1.0 - 2.0*U[i - n - 1] * U[j - n - 1]); }
	      else{ Stemp.push_back(-2.0*U[i - n - 1] * U[j - n - 1]); }
	    }
	}
      
      vector<double> Atemp; // A*S
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + A[veclength*a + d] * Stemp[veclength*d + b];
		}
	      Atemp.push_back(tempel);
	    }
	}
      A = Atemp;
      Atemp.clear();
      // S^T dot A
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + Stemp[veclength*a + d] * A[veclength*d + b]; //Stemp transpose
		}
	      Atemp.push_back(tempel);
	    }
	}
      A = Atemp;

      //Get new S = S*Stemp for vectors
      vector<double> Stottemp;
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + S[veclength*a + d] * Stemp[veclength*d + b];
		}
	      Stottemp.push_back(tempel);
	    }
	}
      S = Stottemp;
    }
  
  for (int i = 0; i < veclength; ++i)
    {
      alpha.push_back(A[veclength*i + i]);
    }
  
  for (int i = 1; i < veclength; ++i)
    {
      beta.push_back(A[veclength*i + (i - 1)]);
    }
  
  alphabeta.push_back(alpha);
  alphabeta.push_back(beta);
  alphabeta.push_back(S);
  
  return alphabeta;
}


vector<vector<double> > QR(vector<double> Alpha, vector<double> Beta)
{
  vector<vector<double> > givens;
  
  int veclength = int(Alpha.size());
  
  vector<double> A;
  for (int i = 0; i < veclength; ++i)
    {
      for (int j = 0; j < veclength; ++j)
	{
	  if (i == j){ A.push_back(Alpha[i]); }
	  else if (i + 1 == j){ A.push_back(Beta[i]); }
	  else if (i - 1 == j){ A.push_back(Beta[i - 1]); }
	  else{ A.push_back(0.0); }
	}
    }
  
  vector<double> Q;
  for (int i = 0; i < veclength; ++i)
    {
      for (int j = 0; j < veclength; ++j)
	{
	  if (i == j){ Q.push_back(1.0); }
	  else{ Q.push_back(0.0); }
	}
    }
  
  for (int i = 0; i < veclength - 1; ++i)
    {
      double r = sqrt(A[veclength*i + i] * A[veclength*i + i] + A[veclength*(i + 1) + i] * A[veclength*(i + 1) + i]);
      double c = A[veclength*i + i] / r;
      double s = -A[veclength*(i + 1) + i] / r;
      
      //std::cout << endl << "alpha,beta = " << A[veclength*i + i] << "," << A[veclength*(i + 1) + i] << "   r,c,s = " << r << ", " << c << ", " << s << endl;
      
      vector<double> Q1;
      for (int a = 0; a < veclength; ++a)
	{
	  for (int b = 0; b < veclength; ++b)
	    {
	      if (a == b && (a == i || a == i + 1)){ Q1.push_back(c); }
	      else if (a + 1 == b && a == i){ Q1.push_back(-s); }
	      else if (a == b + 1 && a == i + 1){ Q1.push_back(s); }
	      else if (a == b){ Q1.push_back(1.0); }
	      else { Q1.push_back(0.0); }
	    }
	}
      
      vector<double> Atemp;
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + Q1[veclength*a + d] * A[veclength*d + b];
		}
	      Atemp.push_back(tempel);
	    }
	}
      A = Atemp;
      
      vector<double> Qtemp;
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + Q[veclength*a + d] * Q1[veclength*b + d]; //Q1 transpose
		}
	      Qtemp.push_back(tempel);
	    }
	}
      Q = Qtemp;
      
    }
  
  /*for (int a = 0; a < veclength; ++a)
    {
    for (int b = 0; b < veclength; ++b)
    {
    std::cout << A[veclength*a + b] << " ";
    }
    std::cout << endl;
    }
    std::cout << endl;
    
    for (int a = 0; a < veclength; ++a)
    {
    for (int b = 0; b < veclength; ++b)
    {
    std::cout << Q[veclength*a + b] << " ";
    }
    std::cout << endl;
    }
    std::cout << endl;*/
  
  givens.push_back(A);
  givens.push_back(Q);
  
  return givens;
}

vector<double> QR2(vector<double> Alpha, vector<double> Beta, vector<double> Vm)
{

  vector<vector<double> > QR2;
  int veclength = int(Alpha.size());

  vector<double> A;
  vector<double> Qtot;
  //initialize A = tridiagonal, Qtot = I
  for (int i = 0; i < veclength; ++i)
    {
      for (int j = 0; j < veclength; ++j)
	{
	  if (i == j){ A.push_back(Alpha[i]); Qtot.push_back(1.0); }
	  else if (i + 1 == j){ A.push_back(Beta[i]); Qtot.push_back(0.0); }
	  else if (i - 1 == j){ A.push_back(Beta[i - 1]); Qtot.push_back(0.0); }
	  else{ A.push_back(0.0); Qtot.push_back(0.0); }
	}
    }

  double err = 1000.0;
  
  while (err > 1e-6)
    {
      //Get diagonal (Alphas) and off-diagonal (Betas) of A
      vector<double> Alphas, Betas;
      for (int i = 0; i < veclength; ++i)
	{
	  for (int j = 0; j < veclength; ++j)
	    {
	      if (i == j){ Alphas.push_back(A[veclength*i + j]); }
	      else if (i + 1 == j){ Betas.push_back(A[veclength*i + j]); }
	    }
	}
      
      //Get QR decomposition of A
      vector<vector<double> > factorization = QR(Alphas, Betas);
      vector<double> R = factorization[0];
      vector<double> Q = factorization[1];
      
      vector<double> Atemp;
      vector<double> Qtottemp;
      err = 0.0;
      int count = 0;

      //Get new A = R*Q
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + R[veclength*a + d] * Q[veclength*d + b];
		}
	      Atemp.push_back(tempel);
	      if (a != b){ err = err + abs(tempel); ++count; }
	    }
	}

      //Get new Qtot = Qtot*Q
      for (int a = 0; a < veclength; ++a) //left matrix row index
	{
	  for (int b = 0; b < veclength; ++b) //right matrix column index
	    {
	      double tempel = 0.0;
	      for (int d = 0; d < veclength; ++d) //move across row/column
		{
		  tempel = tempel + Qtot[veclength*a + d] * Q[veclength*d + b];
		}
	      Qtottemp.push_back(tempel);
	    }
	}
      
      A = Atemp;
      Qtot = Qtottemp;
      err = err / count;
      
    }
  
  /*std::cout << endl;
  for (int a = 0; a < veclength; ++a)
    {
      for (int b = 0; b < veclength; ++b)
	{
	  std::cout << Qtot[veclength*a + b] << " ";
	}
      std::cout << endl;
    }
  std::cout << endl;

  std::cout << endl << "test" << endl;*/

  //Get eigenvector matrix = Vm*Qtot (eigenvectors are columns)
  vector<double> Vec;
  for (int a = 0; a < veclength; ++a) //left matrix row index
    {
      for (int b = 0; b < veclength; ++b) //right matrix column index
	{
	  double tempel = 0.0;
	  for (int d = 0; d < veclength; ++d) //move across row/column
	    {
	      tempel = tempel + Vm[veclength*a + d] * Qtot[veclength*d + b];
	    }
	  Vec.push_back(tempel);
	}
    }
  
  vector<double> Output;
  //Build output matrix (E1,v11,v12,...,v1n,E2,v21,v22,...,vnn)
  for (int a = 0; a < veclength; ++a)
    {
      Output.push_back(A[veclength*a + a]);
      for (int b = 0; b < veclength; ++b)
	{
	  Output.push_back(Vec[veclength*b + a]); //indices reversed because vectors are columns
	}
    }
	  

  /*for (int a = 0; a < veclength; ++a)
    {
      for (int b = 0; b < veclength; ++b)
	{
	  std::cout << A[veclength*a + b] << " ";
	}
      std::cout << endl;
    }
  std::cout << endl;


  for (int a = 0; a < veclength; ++a)
    {
      for (int b = 0; b < veclength; ++b)
	{
	  std::cout << Vec[veclength*a + b] << " ";
	}
      std::cout << endl;
    }
    std::cout << endl;*/
  
  QR2.push_back(A);
  QR2.push_back(Vec);
  return Output;
}
