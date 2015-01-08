#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <bitset>
#include <iomanip>
#include <omp.h>
#include "BetaDecay.h"
#include "BitMaxtrixElements.h"
#include "AngularMomentum.h"
#include "bitconfigsetup.h"
#include "HamiltonSolve.h"
#include "IsoSpin.h"
#include "Parity.h"
#include "lanczos.h"

const string PATH = "/home/sam/Documents/ShellModel/files/";

using namespace std;

vector<int> tempconfigs, combination;
//gives all n-combinations of levels, initialize level to 0
void slater(int n, int level, vector<int> inlevels)
{
  if (n == 0)
    { 
      tempconfigs.insert(tempconfigs.end(), combination.begin(), combination.end());
      return;
    };
  for (int i = level; i <= int(inlevels.size()) - n; ++i)
    {
      combination.push_back(inlevels[i]);
      slater(n - 1, i + 1, inlevels);
      combination.pop_back();
    };
};


struct Input_Parameters{
  int N; //number of valence neutrons
  int P; //number of valence protons
  int Mmax; //maximum angular momentum projection (x2)
  int COM; //flag for center-of-mass matrix elements
  int beta; //flag for beta decay
  string LevelScheme; //level scheme path
  string MatrixElements; //matrix elements path
  string COMMatrixElements; //com matrix elements path
};

struct Model_Space{
  int indp; //number of proton major orbits
  int indn; //number of neutron major orbits
  int indtot; //number of total orbits
  vector<int> levelsind; //list of single particle state indicies (1,2...)
  vector<int> levelsn; //list of single particle state principal quantum numbers
  vector<int> levelsl; //list of single particle state orbital angular momentum
  vector<double> levelsj; //list of single particle state total angular momentum
  vector<double> levelsm; //list of single particle state total angular momentum projection
  vector<double> levelst; //list of single particle state isospins
  vector<int> levelschemeind; //list of major shell indicies (1,2...)
  vector<double> shellsj; //list of major shell angular momentums
  vector<double> shellst; //list of major shell isospin projection
  vector<vector<int> > shellsm; //list of single particle state indicies for each major shell
  vector<string> shellsname;
};

struct Many_Body_States{
  vector<int> bit_configs; //list of bit representation multiparticle states as integers
  vector<int> H_Split; //list of submatrix sizes, distinguished by parity and j-projection
  vector<double> proj; //list of configs angular momentum projection
  vector<double> parit; //list of configs parity
};

struct Matrix_Elements{
  vector<vector<int> > Braket; //list of 4 J-scheme states that make up TBME
  vector<int> J; //list of coupled angular momentum (2x)
  vector<int> T; //list of coupled isospin (2x)
  vector<double> OBME; //list of one-body matrix element values
  vector<double> TBME; //list of two-body matrix element values
};

/*struct M_Matrix_Elements{
  vector<vector<int> > Braket; //list of 4 M-scheme states that make up TBME
  vector<double> OBME; //list of one-body matrix element values
  vector<double> TBME; //list of two-body matrix element values
  };*/

struct Hamiltonian{
  vector<vector<double> > Matrices; //list of vectors of length H_Split[i]^2 for each submatrix
  vector<int> H_Split; //list of submatrix sizes, distinguished by parity and j-projection
};

struct Eigen_System{
  vector<vector<double> > Vectors; //vector of each eigenvector for total Hamiltonian
  vector<vector<double> > GroundState; //vector of each ground-states eigenvector
  vector<double> GroundMs; //vector of each ground_state projections
  vector<double> Energies; //vector of eigenenergies corresponding to each eigenvector
  vector<double> proj; //list of configs angular momentum projection
  vector<double> parit; //list of configs parity
};

struct Level_Structure{
  vector<vector<double> > Vectors; //vector of each eigenvector for total Hamiltonian
  vector<vector<double> > GroundState; //vector of each ground-states eigenvector
  vector<double> GroundMs; //vector of each ground-state projection
  vector<double> Energies; //vector of eigenenergies
  vector<double> Js; //vector of final state angular momentum
  vector<double> Ms; //vector of final state angular momentum projections
  vector<double> Ps; //vector of final state parities
  vector<double> Ts; //vector of final state isospins
};






//Initialize program from input file
Input_Parameters Get_Input_Parameters(string infile)
{ 
  Input_Parameters Input;
  
  std::cout << "Reading Input File" << endl;
  std::cout << "------------------" << endl;

  /***Read Input File***/
  string line;    //initialize placeholder string
  ifstream filestream;
  string path = PATH + infile;
  filestream.open(path.c_str());
  if (!filestream.is_open())
    {
      cerr << "Input file, " << path << ", does not exist" << endl; exit(1);
    }

  //find lines that start with '\*'
  int index = 0; //keep track of input line
  while (getline(filestream, line))
    { if(line[0] == '\\' && line[1] == '*')
	{  
	  ++index;
	  size_t colon = line.find(':');
	  if( colon == line.size() - 1 ){ continue; };
	  string substr = line.substr(colon + 2, line.size());
	  switch(index)
	    {
	    case 1:
	      Input.P = atoi(substr.c_str());
	      break;
	    case 2:
	      Input.N = atoi(substr.c_str());
	      break;
	    case 3:
	      Input.Mmax = atoi(substr.c_str());
	      break;
	    case 4:
	      Input.LevelScheme = substr;
	      break;
	    case 5:
	      Input.MatrixElements = substr;
	      if (Input.MatrixElements[Input.MatrixElements.size() - 2] == '_' && Input.MatrixElements[Input.MatrixElements.size() - 1] == 'M')
		{
		  Input.MatrixElements.erase(Input.MatrixElements.end() - 2, Input.MatrixElements.end());
		}
	      break;
	    case 6:
	      Input.COM = atoi(substr.c_str());
	      break;
	    case 7:
	      if(Input.COM == 1)
		{ 
		  Input.COMMatrixElements = substr;
		  if (Input.COMMatrixElements[Input.COMMatrixElements.size() - 2] == '_' && Input.COMMatrixElements[Input.COMMatrixElements.size() - 1] == 'M')
		    {
		      Input.COMMatrixElements.erase(Input.COMMatrixElements.end() - 2, Input.COMMatrixElements.end());
		    }
		  break;
		}
	      else{ break; }
	    case 8:
	      Input.beta = atoi(substr.c_str());
	      break;
	    } 
	}
      else{ continue; };
    }

  std::cout << Input.P << "  valence protons and  " << Input.N << "  valence neutrons" << endl;
  std::cout << "Maximum Angular Momentum = " << Input.Mmax << "/2" << endl;
  std::cout << "Level Scheme = " << Input.LevelScheme << endl;
  std::cout << "Matrix Elements = " << Input.MatrixElements << " " << Input.COMMatrixElements << endl << endl;

  return Input;
}






Model_Space Build_Model_Space(Input_Parameters Parameters)
{
  Model_Space Space;

  std::cout << "Reading Model Space File" << endl;
  std::cout << "------------------------" << endl;

  //Unpacking the Level Scheme
  string fileformat;    //initialize placeholder string and level format
  int coreA, coreZ, TotOrbs, Shells, POrbs, NOrbs;    //initialize core A, core Z, # total orbits, # major shells, # p orbits, # n orbits
  int ind, n, l, j2;    //initialize level index, n, l, and 2j from file. m is added later.
                        //also initialize corresponding vectors. j from file is 2j

  Space.indtot = 0; //keep running counts of single-particle states
  Space.indp = 0; //number of proton single particle states
  Space.indn = 0; //number of neutron single particle states
  
  //open level scheme file named splevels

  ifstream splevels;
  string fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open())
    {
      cerr << "Level Scheme file does not exist" << endl; exit(1);
    };

  //skip lines that start with '!'
  string phline;
  getline(splevels, phline);
  while (phline[0] == '!'){ getline(splevels, phline); }

  //read rest of level scheme parameters
  fileformat = phline;
  getline(splevels, phline);
  istringstream(phline) >> coreA >> coreZ;
  getline(splevels, phline);
  istringstream(phline) >> TotOrbs;
  getline(splevels, phline);
  istringstream(phline) >> Shells >> POrbs >> NOrbs;

  vector<int> indvec(TotOrbs);
  vector<int> nvec(TotOrbs);
  vector<int> lvec(TotOrbs);
  vector<int> j2vec(TotOrbs);

  Space.shellsj.resize(TotOrbs);
  Space.shellst.resize(TotOrbs);
  Space.shellsm.resize(TotOrbs);
  Space.shellsname.resize(TotOrbs);

  std::cout << "Core: Z = " << coreZ << ", A = " << coreA << endl;
  std::cout << "Proton Valence Shells:" << endl;
  int sp_count = 0; //count single-particle states
  char shell;
  for(int i = 0; i < TotOrbs; ++i)
    {
      if(i == POrbs){ std::cout << "Neutron Valence Shells:" << endl; };
      getline(splevels, phline);
      istringstream(phline) >> ind >> n >> l >> j2;
      //std::cout << ind << " " << n << " " << l << " " << j2 << endl;
      indvec[i] = ind;
      nvec[i] = n;
      lvec[i] = l;
      j2vec[i] = j2;
      Space.shellsj[i] = 0.5*j2;
      if(i < POrbs){ Space.shellst[i] = 0.5; }
      else{ Space.shellst[i] = -0.5; };

      switch(l)
	{
	case 0:
	  shell = 's';
	  break;
	case 1:
	  shell = 'p';
	  break;
	case 2:
	  shell = 'd';
	  break;
	case 3:
	  shell = 'f';
	  break;
	case 4:
	  shell = 'g';
	  break;
	case 5:
	  shell = 'h';
	  break;
	default:
	  shell = 'M';
	  break;
	};

      stringstream a, b, c;
      a << n;
      b << shell;
      c << j2;
      string a1 = a.str(), b1 = b.str(), c1 = c.str();
      Space.shellsname[i] = a1 + b1 + c1 + "/2";
      std::cout << Space.shellsname[i] << endl;

      sp_count += j2 + 1;
      Space.shellsm[i].resize(j2 + 1);
    };
  std::cout << endl;

  Space.levelst.resize(sp_count);
  Space.levelsind.resize(sp_count);
  Space.levelsn.resize(sp_count);
  Space.levelsl.resize(sp_count);
  Space.levelsj.resize(sp_count);
  Space.levelsm.resize(sp_count);
  Space.levelschemeind.resize(sp_count);

  if (fileformat == "pn")
    {
      for(int i = 0; i < TotOrbs; ++i)
	{
	  int shelllength = j2vec[i] + 1;
	  for(int j = 0; j < shelllength; ++j)
	    {
	      Space.indtot++; if(indvec[i] <= POrbs){ Space.indp++; Space.levelst[Space.indtot - 1] = 0.5; }
	      else{ Space.indn++; Space.levelst[Space.indtot - 1] = -0.5; };	//protons T_z = 1/2, neutrons T_z = -1/2
	      Space.levelsind[Space.indtot - 1] = Space.indtot;
	      Space.levelsn[Space.indtot - 1] = nvec[i];
	      Space.levelsl[Space.indtot - 1] = lvec[i];
	      Space.levelsj[Space.indtot - 1] = 0.5*j2vec[i];
	      Space.levelsm[Space.indtot - 1] = -0.5*j2vec[i] + j;
	      Space.levelschemeind[Space.indtot - 1] = indvec[i];
	      Space.shellsm[i][j] = Space.indtot;
	    };
	};
    }
  else
    {
      cerr << "Level Scheme is not formatted as 'pn'" << endl; splevels.close(); exit(1);
    }

  //Print Single Particle State Details
  /*for (int i = 0; i < Space.indtot; ++i)
    {
      std::cout << Space.levelst[i] << " " << Space.levelsind[i] << " " << Space.levelsn[i] << " " << Space.levelsl[i] << " " << Space.levelsj[i] << " " << Space.levelsm[i] << " " << Space.levelschemeind[i] << endl;
      }*/

  splevels.close();

  return Space;
}





Many_Body_States Build_Many_Body_States(Input_Parameters Parameters, Model_Space Space)
{
  Many_Body_States States;

  std::cout << "Building Many-Body States" << endl;
  std::cout << "-------------------------" << endl;

  //get all Slater determinants using the recurse function
  long long Num_Pconfigs = factorial(Space.indp)/(factorial(Parameters.P) * factorial(Space.indp - Parameters.P));
  vector<int> Plevels(Space.indp);
  vector<int> Pconfigs(Num_Pconfigs * Parameters.P);
  for (int i = 1; i <= Space.indp; ++i){ Plevels[i-1] = i; }
  tempconfigs.clear();
  slater(Parameters.P, 0, Plevels);
  Pconfigs = tempconfigs;
  tempconfigs.clear();

  long long Num_Nconfigs = factorial(Space.indn)/(factorial(Parameters.N) * factorial(Space.indn - Parameters.N));
  vector<int> Nlevels(Space.indn);
  vector<int> Nconfigs(Num_Nconfigs * Parameters.N);
  for (int i = Space.indp + 1; i <= Space.indtot; ++i){ Nlevels[i-Space.indp-1] = i; } // start index at indp+1 for neutrons
  tempconfigs.clear();
  slater(Parameters.N, 0, Nlevels);
  Nconfigs = tempconfigs;
  tempconfigs.clear();

  vector<int> Totalconfigs(Num_Pconfigs * Num_Nconfigs * (Parameters.P + Parameters.N));

  for (int i = 0; i < Num_Pconfigs; i++)
    {
      for (int j = 0; j < Num_Nconfigs; j++)
	{
	  for (int pi = 0; pi < Parameters.P; pi++)
	    {
	      Totalconfigs[(Num_Nconfigs*(Parameters.P+Parameters.N))*i + (Parameters.P+Parameters.N)*j + pi] = 
		Pconfigs[(Parameters.P)*i + pi];
	    };
	  for (int ni = 0; ni < Parameters.N; ni++)
	    {
	      Totalconfigs[(Num_Nconfigs*(Parameters.P+Parameters.N))*i + (Parameters.P+Parameters.N)*j + Parameters.P + ni] =
		Nconfigs[(Parameters.N)*j + ni];
	    };
	};
    };

  
  // get total angular momentum projection and parity for PN-particle state
  int PN = Parameters.P + Parameters.N;
  int Npn = int(Totalconfigs.size()) / PN;
  States.proj.resize(Npn);
  States.parit.resize(Npn);
  for (int p = 0; p < Npn; p++)
    {
      double tempproj = 0.0;
      double tempparit = 1.0;
      for (int q = 0; q < PN; q++)
	{
	  int pos1 = Totalconfigs[PN * p + q] - 1; // get level index and shift down by one for vector index
	  tempproj = tempproj + Space.levelsm[pos1];
	  tempparit = tempparit * pow(-1.0, Space.levelsl[pos1]);
	};
      States.proj[p] = tempproj;
      States.parit[p] = tempparit;
    };

  // order states by angular momentum projection and parity
  int testint = 0; // count states inside maximum projection
  for(int i = 0; i < int(States.proj.size()); ++i)
    {
      if(abs(States.proj[i]) <= Parameters.Mmax)
	{
	  ++testint;
	};
    };
  
  States.bit_configs = bitconfigsetup(Totalconfigs, PN);

  int tempstate;
  double tempproj;
  double tempparit;

  int HamSplitTemp;

  for(int i = 0; i < Npn; ++i)
    { 
      HamSplitTemp = 1;
      for(int j = i + 1; j < Npn; ++j)
	{
	  if(States.proj[j] == States.proj[i] && States.parit[j] == States.parit[i])
	    {
	      ++i;
	      ++HamSplitTemp;
	      tempstate = States.bit_configs[j];
	      tempproj = States.proj[j];
	      tempparit = States.parit[j];
	      States.bit_configs[j] = States.bit_configs[i];
	      States.proj[j] = States.proj[i];
	      States.parit[j] = States.parit[i];
	      States.bit_configs[i] = tempstate;
	      States.proj[i] = tempproj;
	      States.parit[i] = tempparit;
	    }
	}
      States.H_Split.push_back(HamSplitTemp);
    }
  
  //Print State Details
  /*std::cout << endl;
  for(int i = 0; i < int(States.H_Split.size()); ++i)
    {
      std::cout << States.H_Split[i] << " ";
    }
  std::cout << endl;
  for(int i = 0; i < int(States.bit_configs.size()); ++i)
    {
      std::cout << States.bit_configs[i] << "\t" << (std::bitset<16>) States.bit_configs[i] << " " << States.proj[i] << " " << States.parit[i] << endl;
    }
    std::cout << endl;*/
  
  std::cout << "Total Number of Slater Determinants = " << States.bit_configs.size() << endl << endl;

  return States;

}
  




Matrix_Elements Convert_To_M_Matrix_Elements(string MatrixElements, Model_Space Space, Matrix_Elements J_ME)
{
  Matrix_Elements M_ME;

  std::cout << "Converting Matrix Elements from J-Scheme to M-Scheme" << endl;
  std::cout << "----------------------------------------------------" << endl;

  M_ME.OBME = J_ME.OBME;

  /*for(int i = 0; i < int(J_ME.TBME.size()); ++i)
    {
      int m3, n3, l3, k3, j3, t3;
      double tbme3;
      m3 = J_ME.Braket[i][0];
      n3 = J_ME.Braket[i][1];
      l3 = J_ME.Braket[i][2];
      k3 = J_ME.Braket[i][3];
      j3 = J_ME.J[i];
      t3 = J_ME.T[i];
      tbme3 = J_ME.TBME[i];
      std::cout << m3 << " " << n3 << " " << l3 << " " << k3 << "   " << j3 << " " << t3 << " " << tbme3 << endl;
      };*/

  //Change J-Scheme matrix elements to M-Scheme
  int length1 = int(Space.levelsind.size());
  int length2 = int(J_ME.TBME.size());
  int l2, k2, m2, n2;

  int size = int(length1*length1*(length1 - 1)*(length1 - 1)/8);
  int ind = 0;
  vector<int> braket(4);
  M_ME.Braket.resize(size);
  for(int i = 0; i < size; ++i){ M_ME.Braket[i].resize(4); };
  M_ME.TBME.resize(size);

  double tbodyphase;

  int testint1 = 0;
  
  for (int m1 = 0; m1 < length1 - 1; ++m1)
    {
      for (int n1 = m1 + 1; n1 < length1; ++n1)
	{
	  for (int l1 = m1; l1 < length1 - 1; ++l1)
	    {
	      for (int k1 = l1 + 1; k1 < length1; ++k1)
		{
		  if(m1 == l1 && n1 > k1){ continue; };
		  double angproj1 = Space.levelsm[m1] + Space.levelsm[n1];
		  double angproj2 = Space.levelsm[l1] + Space.levelsm[k1];
		  double isoproj1 = Space.levelst[m1] + Space.levelst[n1];
		  double isoproj2 = Space.levelst[l1] + Space.levelst[k1];
		  if (angproj1 == angproj2 && isoproj1 == isoproj2)
		    {
		      ++testint1;
		      
		      double tempmom1 = abs(Space.levelsj[m1] + Space.levelsj[n1]);
		      int ind1 = 0;
		      vector<double> angmom1(int(tempmom1 - abs(Space.levelsj[m1] - Space.levelsj[n1]) + 1.1)); 
		      while (tempmom1 >= abs(Space.levelsj[m1] - Space.levelsj[n1]) && tempmom1 >= abs(angproj1))
			{ 
			  angmom1[ind1] = tempmom1; ++ind1;
			  tempmom1 = tempmom1 - 1.0;
			};
		      angmom1.resize(ind1);

		      double tempmom2 = abs(Space.levelsj[l1] + Space.levelsj[k1]);
		      int ind2 = 0;
		      vector<double> angmom2(int(tempmom2 - abs(Space.levelsj[l1] - Space.levelsj[k1]) + 1.1)); 
		      while (tempmom2 >= abs(Space.levelsj[l1] - Space.levelsj[k1]) && tempmom2 >= abs(angproj2))
			{
			  angmom2[ind2] = tempmom2; ++ind2;
			  tempmom2 = tempmom2 - 1.0;
			};
		      angmom2.resize(ind2);

		      int size3 = (ind1 >= ind2 ? ind1 : ind2);
		      int ind3 = 0;
		      vector<double> angmom3(size3);
		      for (int i = 0; i < int(angmom1.size()); ++i)
			{
			  for (int j = 0; j < int(angmom2.size()); ++j)
			    {
			      if (angmom1[i] == angmom2[j]){ angmom3[ind3] = angmom1[i]; ++ind3; break; };
			    };
			};
		      angmom3.resize(ind3);

		      vector<double> iso3(2);
		      iso3[0] = 1.0;
		      if (abs(isoproj1) == 0.0)
			{ iso3[1] = 0.0; }
		      else{ iso3.resize(1); };

		      // get interaction file indices that correspond to the level scheme indices
		      m2 = Space.levelschemeind[m1];
		      n2 = Space.levelschemeind[n1];
		      l2 = Space.levelschemeind[l1];
		      k2 = Space.levelschemeind[k1];
		      double tempmatel = 0.0;

		      for (int b = 0; b < int(angmom3.size()); ++b)	//	loop over relevant angular momentum
			{
			  for (int c = 0; c < int(iso3.size()); ++c)	//	loop over relevant isospin
			    {
			      for (int a = 0; a < length2; ++a)	//	loop over interaction file lines
				{
				  tbodyphase = 1.0;
				  double tbody1, tbody2, tbody3, tbody4;
				  tbody1 = J_ME.Braket[a][0]; tbody2 = J_ME.Braket[a][1]; tbody3 = J_ME.Braket[a][2]; tbody4 = J_ME.Braket[a][3]; 
				  if (((tbody1 == m2 && tbody2 == n2 && tbody3 == l2 && tbody4 == k2) ||
				       (tbody1 == l2 && tbody2 == k2 && tbody3 == m2 && tbody4 == n2)) && J_ME.J[a] == angmom3[b] && J_ME.T[a] == iso3[c])
				    {
				      double CGC1 = CGC(Space.levelsj[m1], Space.levelsm[m1], Space.levelsj[n1], Space.levelsm[n1], angmom3[b], angproj1);
				      double CGC2 = CGC(Space.levelsj[l1], Space.levelsm[l1], Space.levelsj[k1], Space.levelsm[k1], angmom3[b], angproj2);
				      double CGC3 = CGC(0.5, Space.levelst[m1], 0.5, Space.levelst[n1], iso3[c], isoproj1);
				      double CGC4 = CGC(0.5, Space.levelst[l1], 0.5, Space.levelst[k1], iso3[c], isoproj2);
				      tempmatel = tempmatel + tbodyphase*CGC1*CGC2*CGC3*CGC4*J_ME.TBME[a];
				    };
				};
			    };
			};

		      double factor1 = 1.0, factor2 = 1.0, factor3;
		      if (Space.levelsn[m1] == Space.levelsn[n1] && Space.levelsj[m1] == Space.levelsj[n1] && Space.levelsl[m1] == Space.levelsl[n1])
			{ factor1 = 2.0; };
		      if (Space.levelsn[l1] == Space.levelsn[k1] && Space.levelsj[l1] == Space.levelsj[k1] && Space.levelsl[l1] == Space.levelsl[k1])
			{ factor2 = 2.0; };
		      factor3 = sqrt(factor1*factor2);
		      
		      braket[0] = m1 + 1; 
		      braket[1] = n1 + 1;
		      braket[2] = l1 + 1;
		      braket[3] = k1 + 1;
		      M_ME.Braket[ind] = braket;
		      M_ME.TBME[ind] = factor3*tempmatel;
		      ++ind;
		    };
		};
	    };
	};
    };
  
  std::cout << endl;

  M_ME.Braket.resize(ind);
  M_ME.TBME.resize(ind);


  ofstream mschemefile;
  mschemefile.open((PATH + MatrixElements + "_M.int").c_str());
  for (int i = 0; i < int(Space.shellsname.size()); ++i)
    {
      if(i < int(Space.shellsname.size())/2)
	{ mschemefile << "! P - " << Space.shellsname[i] << " = "; }
      else{ mschemefile << "! N - " << Space.shellsname[i] << " = "; }
      for (int j = 0; j < int(Space.shellsm[i].size()); ++j)
	{
	  if(j != int(Space.shellsm[i].size()) - 1)
	    { mschemefile << Space.shellsm[i][j] << ", "; }
	  else{ mschemefile << Space.shellsm[i][j] << "\n"; }
	}
    };
  mschemefile << ind << "\t";
  for (int i = 0; i < int(M_ME.OBME.size()); ++i)
    {
      mschemefile << M_ME.OBME[i] << "\t";
    }
  mschemefile << "\n";
  for (int i = 0; i < int(M_ME.TBME.size()); ++i)
    {
      mschemefile << M_ME.Braket[i][0] << "\t" << M_ME.Braket[i][1] << "\t" << M_ME.Braket[i][2] << "\t" << M_ME.Braket[i][3] << "\t";
      mschemefile << M_ME.TBME[i] << "\n";
    }
  mschemefile.close();

  
  std::cout << "Number of M-Scheme Two-Body Matrix Elements = " << ind << endl << endl;
  
  return M_ME;
}




Matrix_Elements Get_Matrix_Elements(Input_Parameters Parameters, Model_Space Space)
{
  Matrix_Elements ME;
  Matrix_Elements COM_ME;

  std::cout << "Reading Matrix Elements" << endl;
  std::cout << "-----------------------" << endl;

  //Get Matrix Elements
  int NumElements;
  double OBME, TBME;
  int shell1, shell2, shell3, shell4, coupJ, coupT; // J-scheme interaction file contents
  ifstream interaction;	// interaction file
  string interactionline; // interaction file line
  char type, comtype; // j or m scheme
  
  string fullpath2 = PATH + Parameters.MatrixElements + ".int";
  string fullpath3 = PATH + Parameters.MatrixElements + "_M.int";
  interaction.open(fullpath3.c_str()); // try m-scheme first
  size_t intsize = Parameters.MatrixElements.size();

  //open interaction file
  if (interaction.is_open())
    { type = 'm'; }
  else if (Parameters.MatrixElements[intsize-2] == '_' && Parameters.MatrixElements[intsize-1] == 'M')
    { type = 'm'; interaction.open(fullpath2.c_str()); }
  else
    { type = 'j'; interaction.open(fullpath2.c_str()); }

  if (!interaction.is_open())
    {
      cerr << "Matrix Element file does not exist" << endl; exit(1);
    }
  
  //skip lines that start with '!'
  getline(interaction, interactionline);
  while (interactionline[0] == '!'){ getline(interaction, interactionline); }

  //read matrix element parameters and one-body matrix elements
  istringstream filestring(interactionline);
  filestring >> NumElements;

  //get one-body matrix elements that correspond to the ordered proton/neutron shells
  while (filestring >> OBME)
    {
      ME.OBME.push_back(OBME);
    }

  if(ME.OBME.size() != Space.shellsname.size())
    { cerr << "Space/Interaction Mismatch" << endl; exit(1); }

  vector<int> braket(4);
  ME.Braket.resize(NumElements);
  for(int i = 0; i < NumElements; ++i){ ME.Braket[i].resize(4); };
  ME.J.resize(NumElements);
  ME.T.resize(NumElements);
  ME.TBME.resize(NumElements);
  
  //read two-body parameters and two-body matrix elements
  double tempS1 = 0, tempS2 = 0, tempS3 = 0, tempS4 = 0, tempJ = -1, tempT = -1;
  for(int i = 0; i < NumElements; ++i)
    {
      getline(interaction, interactionline);
      if(type == 'j')
	{ istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> coupJ >> coupT >> TBME;
	  if (shell1 == tempS3 && shell2 == tempS4 && shell3 == tempS1 && shell4 == tempS2 && coupJ == tempJ && coupT == tempT)
	    { --NumElements; --i; continue; };
	  tempS1 = shell1; tempS2 = shell2; tempS3 = shell3; tempS4 = shell4; tempJ = coupJ; tempT = coupT;
	  if(shell2 < shell1)
	    {
	      swap(shell1, shell2);
	      TBME = TBME * pow(-1.0, int(Space.shellsj[shell1 - 1] + Space.shellsj[shell2 - 1] - coupJ - coupT));
	    }
	  if(shell4 < shell3)
	    {
	      swap(shell3, shell4);
	      TBME = TBME * pow(-1.0, int(Space.shellsj[shell3 - 1] + Space.shellsj[shell4 - 1] - coupJ - coupT));
	    }
	  if((shell3 < shell1) || (shell3 == shell1 && shell4 < shell2))
	    {
	      swap(shell1, shell3);
	      swap(shell2, shell4);
	    }
	  braket[0] = shell1;
	  braket[1] = shell2;
	  braket[2] = shell3;
	  braket[3] = shell4;
	  ME.Braket[i] = braket;
	  ME.J[i] = coupJ;
	  ME.T[i] = coupT;
	  ME.TBME[i] = TBME;
	}
      else if(type == 'm')
	{ istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
	  if (shell1 == tempS3 && shell2 == tempS4 && shell3 == tempS1 && shell4 == tempS2)
	    { --NumElements; --i; continue; };
	  tempS1 = shell1; tempS2 = shell2; tempS3 = shell3; tempS4 = shell4;
	  if(shell2 < shell1)
	    {
	      swap(shell1, shell2);
	      TBME = TBME * -1.0;
	    }
	  if(shell4 < shell3)
	    {
	      swap(shell3, shell4);
	      TBME = TBME * -1.0;
	    }
	  if((shell3 < shell1) || (shell3 == shell1 && shell4 < shell2))
	    {
	      swap(shell1, shell3);
	      swap(shell2, shell4);
	    }
	  braket[0] = shell1;
	  braket[1] = shell2;
	  braket[2] = shell3;
	  braket[3] = shell4;
	  ME.Braket[i] = braket;
	  ME.TBME[i] = TBME;
	}
    }
  interaction.close();

  ME.Braket.resize(NumElements);
  ME.J.resize(NumElements);
  ME.T.resize(NumElements);
  ME.TBME.resize(NumElements);

  //COM matrix elements
  if(Parameters.COM == 1)
    {
      fullpath2 = PATH + Parameters.COMMatrixElements + ".int";
      fullpath3 = PATH + Parameters.COMMatrixElements + "_M.int";
      interaction.open(fullpath3.c_str()); // try m-scheme first
      intsize = Parameters.COMMatrixElements.size();
      
      //open interaction file
      if (interaction.is_open())
	{ comtype = 'm'; }
      else if (Parameters.COMMatrixElements[intsize-2] == '_' && Parameters.COMMatrixElements[intsize-1] == 'M')
	{ comtype = 'm'; interaction.open(fullpath2.c_str()); }
      else
	{ comtype = 'j'; interaction.open(fullpath2.c_str()); }
      
      if (!interaction.is_open())
	{
	  cerr << "COM Matrix Element file does not exist" << endl; exit(1);
	}
      
      //skip lines that start with '!'
      getline(interaction, interactionline);
      while (interactionline[0] == '!'){ getline(interaction, interactionline); }
      
      //read matrix element parameters and one-body matrix elements
      istringstream filestring(interactionline);
      filestring >> NumElements;
      
      //get one-body matrix elements that correspond to the ordered proton/neutron shells
      while (filestring >> OBME)
	{
	  COM_ME.OBME.push_back(OBME);
	}
      
      if(ME.OBME.size() != Space.shellsname.size())
	{ cerr << "Space/COM Interaction Mismatch" << endl; exit(1); }

      vector<int> braket(4);
      COM_ME.Braket.resize(NumElements);
      for(int i = 0; i < NumElements; ++i){ COM_ME.Braket[i].resize(4); };
      COM_ME.J.resize(NumElements);
      COM_ME.T.resize(NumElements);
      COM_ME.TBME.resize(NumElements);
      
      //read two-body parameters and two-body matrix elements
      tempS1 = 0; tempS2 = 0; tempS3 = 0; tempS4 = 0; tempJ = -1; tempT = -1;
      for(int i = 0; i < NumElements; ++i)
	{
	  getline(interaction, interactionline);
	  if(comtype == 'j')
	    { istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> coupJ >> coupT >> TBME;
	      if ((shell1 == tempS3 && shell2 == tempS4 && shell3 == tempS1 && shell4 == tempS2 && coupJ == tempJ && coupT == tempT) || abs(TBME) < 0.000001)
		{ --NumElements; --i; continue; };
	      tempS1 = shell1; tempS2 = shell2; tempS3 = shell3; tempS4 = shell4; tempJ = coupJ; tempT = coupT;
	      if(shell2 < shell1)
		{
		  swap(shell1, shell2);
		  TBME = TBME * pow(-1.0, int(Space.shellsj[shell1 - 1] + Space.shellsj[shell2 - 1] - coupJ - coupT));
		}
	      if(shell4 < shell3)
		{
		  swap(shell3, shell4);
		  TBME = TBME * pow(-1.0, int(Space.shellsj[shell3 - 1] + Space.shellsj[shell4 - 1] - coupJ - coupT));
		}
	      if((shell3 < shell1) || (shell3 == shell1 && shell4 < shell2))
		{
		  swap(shell1, shell3);
		  swap(shell2, shell4);
		}
	      braket[0] = shell1;
	      braket[1] = shell2;
	      braket[2] = shell3;
	      braket[3] = shell4;
	      COM_ME.Braket[i] = braket;
	      COM_ME.J[i] = coupJ;
	      COM_ME.T[i] = coupT;
	      COM_ME.TBME[i] = TBME;
	    }
	  else if(comtype == 'm')
	    { istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
	      if ((shell1 == tempS3 && shell2 == tempS4 && shell3 == tempS1 && shell4 == tempS2) || abs(TBME) < 0.0001)
		{ --NumElements; --i; continue; };
	      tempS1 = shell1; tempS2 = shell2; tempS3 = shell3; tempS4 = shell4;
	      if(shell2 < shell1)
		{
		  swap(shell1, shell2);
		  TBME = TBME * -1.0;
		}
	      if(shell4 < shell3)
		{
		  swap(shell3, shell4);
		  TBME = TBME * -1.0;
		}
	      if((shell3 < shell1) || (shell3 == shell1 && shell4 < shell2))
		{
		  swap(shell1, shell3);
		  swap(shell2, shell4);
		}
	      braket[0] = shell1;
	      braket[1] = shell2;
	      braket[2] = shell3;
	      braket[3] = shell4;
	      COM_ME.Braket[i] = braket;
	      COM_ME.TBME[i] = TBME;
	    }
	}
      interaction.close();
      
      COM_ME.Braket.resize(NumElements);
      COM_ME.J.resize(NumElements);
      COM_ME.T.resize(NumElements);
      COM_ME.TBME.resize(NumElements);
    }

  std::cout << "Number of One-Body Matrix Elements = " << ME.OBME.size() << endl;
  if(type == 'j')
    {std::cout << "Number of J-Scheme Two-Body Matrix Elements = " << ME.TBME.size() << endl << endl;}
  else if(type == 'm')
    {std::cout << "Number of M-Scheme Two-Body Matrix Elements = " << ME.TBME.size() << endl << endl;}

  if(Parameters.COM == 1)
    {
      std::cout << "Number of COM One-Body Matrix Elements = " << COM_ME.OBME.size() << endl;
      if(comtype == 'j')
	{std::cout << "Number of COM J-Scheme Two-Body Matrix Elements = " << COM_ME.TBME.size() << endl << endl;}
      else if(comtype == 'm')
	{std::cout << "Number of COM M-Scheme Two-Body Matrix Elements = " << COM_ME.TBME.size() << endl << endl;}
    }


  /*for(int i = 0; i < int(ME.TBME.size()); ++i)
    {
      int m3, n3, l3, k3, j3, t3;
      double tbme3;
      m3 = ME.Braket[i][0];
      n3 = ME.Braket[i][1];
      l3 = ME.Braket[i][2];
      k3 = ME.Braket[i][3];
      j3 = ME.J[i];
      t3 = ME.T[i];
      tbme3 = ME.TBME[i];
      std::cout << m3 << " " << n3 << " " << l3 << " " << k3 << "   " << j3 << " " << t3 << " " << tbme3 << endl;
      };*/

					  
  if(type == 'j')
    { 
      ME = Convert_To_M_Matrix_Elements(Parameters.MatrixElements, Space, ME);
    }
  if(Parameters.COM == 1)
    {
      if(comtype == 'j')
	{ 
	  COM_ME = Convert_To_M_Matrix_Elements(Parameters.COMMatrixElements, Space, COM_ME);
	}
      //add COM
      int m1, n1, l1, k1, m2, n2, l2, k2;
      for(int i = 0; i < int(ME.OBME.size()); ++i)
	{ ME.OBME[i] += 50.0 * COM_ME.OBME[i]; }
      for(int i = 0; i < int(ME.TBME.size()); ++i)
	{
	  m1 = ME.Braket[i][0], n1 = ME.Braket[i][1], l1 = ME.Braket[i][2], k1 = ME.Braket[i][3];
	  for(int j = 0; j < int(COM_ME.TBME.size()); ++j)
	    {
	      m2 = COM_ME.Braket[j][0], n2 = COM_ME.Braket[j][1], l2 = COM_ME.Braket[j][2], k2 = COM_ME.Braket[j][3];
	      if(m1 == m2 && n1 == n2 && l1 == l2 && k1 == k2)
		{ ME.TBME[i] += 50.0 * COM_ME.TBME[j]; break; }
	    }
	} 
    }

  int m1, n1, l1, k1, j1, t1, m2, n2, l2, k2, j2, t2;
  for(int i = 0; i < int(ME.OBME.size()); ++i)
    { ME.OBME[i] += 50.0 * COM_ME.OBME[i]; }
  for(int i = 0; i < int(ME.TBME.size()); ++i)
    {
      m1 = ME.Braket[i][0], n1 = ME.Braket[i][1], l1 = ME.Braket[i][2], k1 = ME.Braket[i][3];
      if(type == 'j'){ j1 = ME.J[i]; t1 = ME.T[i]; };
      for(int j = 0; j < int(COM_ME.TBME.size()); ++j)
	{
	  m2 = COM_ME.Braket[j][0], n2 = COM_ME.Braket[j][1], l2 = COM_ME.Braket[j][2], k2 = COM_ME.Braket[j][3];
	  if(comtype == 'j'){ j2 = COM_ME.J[j]; t2 = COM_ME.T[j]; };
	  if((m1 == m2 && n1 == n2 && l1 == l2 && k1 == k2) && ((type == 'j' && comtype == 'j' && j1 == j2 && t1 == t2) || (type == 'm' && comtype == 'm')))
	    { ME.TBME[i] += 50.0 * COM_ME.TBME[j]; break; }
	}
    }
  
  return ME;

};





Hamiltonian Build_Hamiltonian(Model_Space Space, Many_Body_States States, Matrix_Elements M_ME)
{
   Hamiltonian Ham;

   std::cout << "Building Hamiltonian Matrix" << endl;
   std::cout << "---------------------------" << endl;

   Ham.H_Split = States.H_Split;
  
   /*for (int n = 0; n < int(States.bit_configs.size()); ++n)
     {
     std::cout << States.bit_configs[n] << endl;
     }*/
   
  Ham.Matrices.resize(Ham.H_Split.size());
  int totallength = 0;
  for (int n = 0; n < int(Ham.H_Split.size()); ++n)
    {
      totallength += (Ham.H_Split[n] * (Ham.H_Split[n] + 1)) / 2;
      Ham.Matrices[n].resize(Ham.H_Split[n]*Ham.H_Split[n]);
    }

  //multiplicative strengths for all one and two-body matrix elements, respectively
  double strength1 = 1.0;
  double strength2 = 1.0;  //pow(18.0/(16.0 + PN),0.3);

  /*int HamiltonianCount = 0;
  double perc = 0.0;
  int indexoffset = 0;
  int length, bra, ket;
  double ME;

  for (int n = 0; n < Ham.H_Split.size(); ++n)
    {
      length = Ham.H_Split[n];
      for (int i = 0; i < length; ++i)
	{
	  bra = States.bit_configs[indexoffset + i];
	  for (int j = i; j < length; ++j)
	    {
	      ++HamiltonianCount; ket = States.bit_configs[indexoffset + j];
	      if (((HamiltonianCount * 100.0 / totallength) - perc) >= 1.0)
		{
		  perc = HamiltonianCount * 100.0 / totallength;
		  std::cerr << "Building Hamiltonian: " << int(HamiltonianCount * 100.0 / totallength) << " % \r";
		};
	      ME = matrixe(strength1, strength2, Space.levelschemeind, bra, ket, M_ME.OBME, M_ME.Braket, M_ME.TBME);
	      Ham.Matrices[n][i * length + j] = ME;
	      if (i != j){ Ham.Matrices[n][j * length + i] = ME; };
	    };
	};
      indexoffset += length;
      };*/

  //vector<int> HamiltonianCount(Ham.H_Split.size(), 0);
  //double perc = 0.0;

  #pragma omp parallel for
  for(int n = 0; n < int(Ham.H_Split.size()); ++n)
    {
      int indexoffset = 0;
      for (int i = 0; i < n; ++i){ indexoffset += Ham.H_Split[i]; }
      int length = Ham.H_Split[n];
      //#pragma omp parallel for
      for (int i = 0; i < length; ++i)
	{
	  int bra = States.bit_configs[indexoffset + i];
	  //#pragma omp parallel for
	  for (int j = i; j < length; ++j)
	    {
	      /*++HamiltonianCount[n];
	      double tempperc = (std::accumulate(HamiltonianCount.begin(), HamiltonianCount.end(), 0) / totallength) * 100.0;
	      std::cout << tempperc << endl;
	      if(tempperc - perc > 1.0)
		{
		  std::cerr << "Building Hamiltonian: " << int(tempperc) << " % \r";
		  perc = tempperc;
		}*/
	      int ket = States.bit_configs[indexoffset + j];
	      double ME = matrixe(strength1, strength2, Space.levelschemeind, bra, ket, M_ME.OBME, M_ME.Braket, M_ME.TBME);
	      //std::cout << bra << " " << ket << "    " << ME << endl;
	      Ham.Matrices[n][i * length + j] = ME;
	      if (i != j){ Ham.Matrices[n][j * length + i] = ME; };
	    }
	}
    }

  std::cout << endl;
  std::cout.precision(5);
  for(int i = 0; i < int(Ham.H_Split.size()); ++i)
    {
      for(int t = 0; t < sqrt(int(Ham.Matrices[i].size())); ++t)
	{
	  for(int h = 0; h < sqrt(int(Ham.Matrices[i].size())); ++h)
	    {
	      std::cout << Ham.Matrices[i][t*sqrt(Ham.Matrices[i].size()) + h] << " ";
	    }
	  std::cout << endl;
	}
      std::cout << endl;
    }
  
  return Ham;
}
  



Eigen_System Solve_Hamiltonian(Hamiltonian Ham, Many_Body_States States)
{
  Eigen_System Eigen;

  Eigen.proj = States.proj;
  Eigen.parit = States.parit;

  int Hamtot = 0;
  for (int i = 0; i < int(Ham.H_Split.size()); ++i)
    { 
      Hamtot = Hamtot + Ham.H_Split[i];
    };

  Eigen.Energies.resize(Hamtot);
  Eigen.Vectors.resize(Hamtot);
  for(int i = 0; i < Hamtot; ++i)
    {
      Eigen.Vectors[i].resize(Hamtot);
    };
  

  int testcount1 = 0;
  int Hamtemp = 0;
  for (int i = 0; i < int(Ham.Matrices.size()); ++i)  //!!! Add Option to choose Lanczos or Householder !!!
    {
      vector<double> solutionstemp = solution(Ham.Matrices[i], Ham.H_Split[i]);
      //vector<vector<double> > tridiagonal1 = Lanczos(Ham.Matrices[i]);
      //vector<vector<double> > tridiagonal2 = Householder(Ham.Matrices[i]);
      //vector<double> solutionstemp = QR2(tridiagonal1[0], tridiagonal1[1], tridiagonal1[2]);
      //vector<double> solutionstemp = QR2(tridiagonal2[0], tridiagonal2[1], tridiagonal2[2]);
      ++testcount1;
      std::cerr << "Solving Hamiltonian: " << int(testcount1 * 100.0 / Ham.Matrices.size()) << " % \r";
      for (int j = 0; j < int(solutionstemp.size())/(Ham.H_Split[i] + 1); ++j)
	{ 
	  //std::cout << solutionstemp[j * (Ham.H_Split[i] + 1)] << endl;
	  Eigen.Energies[Hamtemp + j] = solutionstemp[j*(Ham.H_Split[i] + 1)];
	  for (int k = 0; k < Hamtot; ++k)
	    { 
	      if(k >= Hamtemp && k < (Hamtemp + Ham.H_Split[i]))
		{ Eigen.Vectors[Hamtemp + j][k] = solutionstemp[j*(Ham.H_Split[i] + 1) + k + 1 - Hamtemp]; }
	      else{ Eigen.Vectors[Hamtemp + j][k] = 0.0; };
	    }
	}
      Hamtemp = Hamtemp + Ham.H_Split[i];
    }  
  std::cout << endl;

  //Order States by Energy
  for (int i = 0; i < int(Eigen.Energies.size()); ++i)
    {
      double temp1 = Eigen.Energies[i];
      vector<double> tempvec = Eigen.Vectors[i];
      double tempproj = Eigen.proj[i];
      double tempparit = Eigen.parit[i];
      int index = i; //initialize
      for (int j = i + 1; j < int(Eigen.Energies.size()); ++j)
	{
	  double temp2 = Eigen.Energies[j];
	  if(temp1 - temp2 <= 0.05){ continue; }
	  else{
	    temp1 = Eigen.Energies[j];
	    tempvec = Eigen.Vectors[j];
	    tempproj = Eigen.proj[j];
	    tempparit = Eigen.parit[j];
	    index = j;
	  };
	};
      Eigen.Energies[index] = Eigen.Energies[i];
      Eigen.Vectors[index] = Eigen.Vectors[i];
      Eigen.proj[index] = Eigen.proj[i];
      Eigen.parit[index] = Eigen.parit[i];
      Eigen.Energies[i] = temp1;
      Eigen.Vectors[i] = tempvec;
      Eigen.proj[i] = tempproj;
      Eigen.parit[i] = tempparit;
    };
  
  // Get Ground State
  Eigen.GroundState.push_back(Eigen.Vectors[0]);
  Eigen.GroundMs.push_back(Eigen.proj[0]);
  for (int i = 1; i < int(Eigen.Energies.size()); ++i)
    {
      if (abs(Eigen.Energies[i] - Eigen.Energies[i-1]) < 0.05)
	{
	  Eigen.GroundState.push_back(Eigen.Vectors[i]);
	  Eigen.GroundMs.push_back(Eigen.proj[i]);
	}
      else
	{
	  break;
	}
    };

  //Remove Degenerate States
  vector<double> EnergiesTemp(Hamtot);
  vector<int> flagvec(Hamtot);
  int flagcount = 0;
  vector<vector<double> > VectorsTemp(Hamtot);
  vector<double> ProjTemp(Hamtot);
  vector<double> ParitTemp(Hamtot);
  int ind = 0;
  for (int i = 0; i < int(Eigen.Energies.size()); ++i)
    {
      double temp1 = Eigen.Energies[i];
      for (int k = 0; k < flagcount; ++k)
	{
	  int flag2 = flagvec[k];
	  if(i == flag2){ goto contflag; };
	};
      EnergiesTemp[ind] = Eigen.Energies[i];
      VectorsTemp[ind] = Eigen.Vectors[i];
      ProjTemp[ind] = Eigen.proj[i];
      ParitTemp[ind] = Eigen.parit[i];
      ++ind;
      for (int j = i + 1; j < int(Eigen.Energies.size()); ++j)
	{
	  double temp2 = Eigen.Energies[j];
	  if(abs(temp1 - temp2) < 0.05)
	    {
	      flagvec[flagcount] = j;
	      ++flagcount;
	    };
	};
    contflag:;
    };

  Eigen.Energies = EnergiesTemp;
  Eigen.Vectors = VectorsTemp;
  Eigen.proj = ProjTemp;
  Eigen.parit = ParitTemp;
  vector<double>().swap(EnergiesTemp);
  vector<vector<double> >().swap(VectorsTemp);
  vector<double>().swap(ProjTemp);
  vector<double>().swap(ParitTemp);
  vector<int>().swap(flagvec);
  Eigen.Energies.resize(ind);
  Eigen.Vectors.resize(ind);
  Eigen.proj.resize(ind);
  Eigen.parit.resize(ind);

  return Eigen;
}



Level_Structure Get_Level_Structure(Model_Space Space, Many_Body_States States, Eigen_System Eigen)
{
  Level_Structure Structure;
  Structure.Energies = Eigen.Energies;
  Structure.Vectors = Eigen.Vectors;
  Structure.GroundState = Eigen.GroundState;
  Structure.GroundMs = Eigen.GroundMs;

  int length3 = int(Eigen.Energies.size());
  Structure.Js.resize(length3);
  Structure.Ts.resize(length3);

  int testcount1 = 0;
  for (int i = 0; i < length3; ++i)
    {
      ++testcount1;
      vector<double> statebra = Eigen.Vectors[i];
      std::cerr << "Calculating Angular Momentum: " << int(testcount1 * 100.0 / length3) << " % \r";
      double stateJ = StateJ(testcount1, statebra, States.bit_configs, Space.levelschemeind, Space.levelsj, Space.levelsm, Space.shellsm);
      Structure.Js[i] = stateJ;
    };
  std::cout << endl;
  
  Structure.Ms = Eigen.proj;
  Structure.Ps = Eigen.parit;

  return Structure;
}



void Results(Input_Parameters Parameters, Level_Structure Structure)
{ 
  stringstream a, b;
  a << Parameters.P;
  b << Parameters.N;
  string filenum1 = a.str(), filenum2 = b.str();
  
  ofstream myfile;
  myfile.open((PATH + Parameters.LevelScheme + "." + Parameters.MatrixElements + "." + 
	       filenum1 + "." + filenum2 + ".solution.txt").c_str());
  for (int i = 0; i < int(Structure.Energies.size()); ++i)
    {
      myfile << Structure.Energies[i] << "     \t" << Structure.Js[i] << "\t" << Structure.Ps[i] << "\n";
    };
    myfile.close();
}


int main(int argc, char* argv[])
{ 
  clock_t t1, t2;    	//initialize program clock
  t1 = clock();

  string inputfile = "inputbeta.dat";
  Input_Parameters Parameters = Get_Input_Parameters(inputfile);
  Model_Space Space = Build_Model_Space(Parameters);
  Many_Body_States States = Build_Many_Body_States(Parameters, Space);
  Matrix_Elements ME = Get_Matrix_Elements(Parameters, Space);
  //M_Matrix_Elements M_ME = Convert_To_M_Matrix_Elements(Parameters, Space, J_ME);

  Hamiltonian Ham = Build_Hamiltonian(Space, States, ME);
  Eigen_System Eigen = Solve_Hamiltonian(Ham, States);
  Level_Structure Structure = Get_Level_Structure(Space, States, Eigen);
  Results(Parameters, Structure);
  std::cout << endl << endl;

  if (Parameters.beta == 1)
    {
      Input_Parameters Parameters2 = Parameters;
      Parameters2.P = Parameters2.P + 1;
      Parameters2.N = Parameters2.N - 1;
      Many_Body_States States2 = Build_Many_Body_States(Parameters2, Space);
      Hamiltonian Ham2 = Build_Hamiltonian(Space, States2, ME);
      Eigen_System Eigen2 = Solve_Hamiltonian(Ham2, States2);
      Level_Structure Structure2 = Get_Level_Structure(Space, States2, Eigen2);
      Results(Parameters2, Structure2);

      double betamtot = 0.0;

      for (int i = 0; i < int(Structure.GroundState.size()); ++i)
	{
	  for (int j = 0; j < int(Structure2.GroundState.size()); ++j)
	    {
	      double betam = BetaMDecay(Structure.GroundState[i], Structure2.GroundState[j], States.bit_configs, States2.bit_configs, Space.indp, Space.indn, Space.levelsn, Space.levelsl, Space.levelsj, Space.levelsm);
	      betamtot += betam*betam;
	    }
	}
      std::cout << endl << sqrt(betamtot) << endl << endl;
    }

  t2 = clock();
  float diff((float)t2 - (float)t1);
  std::cout << diff / CLOCKS_PER_SEC << "sec" << endl;

  int a;
  std::cin >> a;

  return 0;

}
