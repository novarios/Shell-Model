#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <bitset>
#include "BetaDecay.h"
#include "BitMaxtrixElements.h"
#include "AngularMomentum.h"
#include "CGC.h"
#include "bitconfigsetup.h"
#include "HamiltonSolve.h"
#include "IsoSpin.h"
#include "Parity.h"
#include "lanczos.h"

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
  for (int i = level; i <= inlevels.size() - n; ++i)
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
  vector<vector<int> > shellsm; //list of single particle state indicies for each major shell
};

struct Many_Body_States{
  vector<int> bit_configs; //list of bit representation multiparticle states as integers
  vector<int> H_Split; //list of submatrix sizes, distinguished by parity and j-projection
};

struct J_Matrix_Elements{
  vector<vector<int> > Braket; //list of 4 J-scheme states that make up TBME
  vector<int> J; //list of coupled angular momentum (2x)
  vector<int> T; //list of coupled isospin (2x)
  vector<double> OBME; //list of one-body matrix element values
  vector<double> TBME; //list of two-body matrix element values
};

struct M_Matrix_Elements{
  vector<vector<int> > Braket; //list of 4 M-scheme states that make up TBME
  vector<double> OBME; //list of one-body matrix element values
  vector<double> TBME; //list of two-body matrix element values
};

struct Hamiltonian{
  vector<vector<double> > Matrices; //list of vectors of length H_Split[i]^2 for each submatrix
  vector<int> H_Split; //list of submatrix sizes, distinguished by parity and j-projection
};

struct Eigen_System{
  vector<vector<double> > Vectors; //vector of each eigenvector for total Hamiltonian
  vector<double> Energies; //vector of eigenenergies corresponding to each eigenvector
};

struct Level_Structure{
  vector<vector<double> > Vectors; //vector of each eigenvector for total Hamiltonian
  vector<double> Energies; //vector of eigenenergies
  vector<int> Js; //vector of final state angular momentum (x2)
  vector<int> Ps; //vector of final state parities
  vector<int> Ts; //vector of final state isospins
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
  string path = "C:\\Users\\Sam\\Documents\\Shell Model\\files\\" + infile;
  filestream.open(path);
  if (!filestream.is_open())
    {
      cerr << "Input file does not exist" << endl; exit(1);
    }

  //find lines that start with '\*'
  int index = 0; //keep track of input line
  while (getline(filestream, line))
    { if(line[0] == '\\' && line[1] == '*')
	{  
	  ++index;
	  size_t colon = line.find(':');
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
	      break;
	    case 6:
	      Input.COM = atoi(substr.c_str());
	      break;
	    case 7:
	      if(Input.COM == 1){ Input.COMMatrixElements = substr; }
	      else{ break; }
	    case 8:
	      Input.beta = substr;
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
  string fullpath = "C:\\Users\\Sam\\Documents\\Shell Model\\files\\" + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath);
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
  Space.shellsm.resize(TotOrbs);

  std::cout << "Core: Z = " << coreZ << ", A = " << coreA << endl;
  std::cout << "Proton Valence Shells:" << endl;
  int sp_count = 0; //count single-particle states
  char shell;
  for(int i = 0; i < TotOrbs; ++i)
    {
      if(i == POrbs){ std::cout << "Neutron Valence Shells:" << endl; };
      getline(splevels, phline);
      istringstream(phline) >> ind >> n >> l >> j2;
      indvec[i] = ind;
      nvec[i] = n;
      lvec[i] = l;
      j2vec[i] = j2;
      Space.shellsj[i] = 0.5*j2;

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

      std::cout << n << shell << " " << j2 << "/2" << endl;
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
      splevels.close();
    }
  else
    {
      cerr << "Level Scheme is not formatted as 'pn'" << endl; exit(1);
    }

  return Space;
}





Many_Body_States Build_Many_Body_States(Input_Parameters Parameters, Model_Space Space)
{
  Many_Body_States States;

  std::cout << "Building Many-Body States" << endl;
  std::cout << "-------------------------" << endl;

  //get all Slater determinants using the recurse function
  int Num_Pconfigs = intfac(Space.indp)/(intfac(Parameters.P)*intfac(Space.indp - Parameters.P));
  vector<int> Plevels(Space.indp);
  vector<int> Pconfigs(Num_Pconfigs * Parameters.P);
  for (int i = 1; i <= Space.indp; ++i){ Plevels[i-1] = i; }
  tempconfigs.clear();
  slater(Parameters.P, 0, Plevels);
  Pconfigs = tempconfigs;
  tempconfigs.clear();

  int Num_Nconfigs = intfac(Space.indn)/(intfac(Parameters.N)*intfac(Space.indn - Parameters.N));
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
  size_t Npn = Totalconfigs.size() / PN;
  vector<double> proj(Npn);
  vector<double> parit(Npn);
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
      proj[p] = tempproj;
      parit[p] = tempparit;
    };

  // order states by angular momentum projection and parity
  int testint = 0; // count states inside maximum projection
  for(int i = 0; i < proj.size(); ++i)
    {
      if(abs(proj[i]) <= Parameters.Mmax)
	{
	  ++testint;
	};
    };

  States.bit_configs = bitconfigsetup(Totalconfigs, PN);

  int tempstate;
  double tempproj;
  double tempparit;

  int HamSplitTemp;
  //std::cout << endl;
  for(int i = 0; i < Npn; ++i)
    { 
      /*for(int j = 0; j < Npn; ++j)
	{ std::cout << proj[j] << " "; }
	std::cout << endl;*/
      if(abs(proj[i]) <= Parameters.Mmax)
	{
	  HamSplitTemp = 1;
	  for(int j = i + 1; j < Npn; ++j)
	    {
	      if(proj[j] == proj[i] && parit[j] == parit[i])
		{
		  ++i;
		  ++HamSplitTemp;
		  tempstate = States.bit_configs[j];
		  tempproj = proj[j];
		  tempparit = parit[j];
		  States.bit_configs[j] = States.bit_configs[i];
		  proj[j] = proj[i];
		  parit[j] = parit[i];
		  States.bit_configs[i] = tempstate;
		  proj[i] = tempproj;
		  parit[i] = tempparit;
		}
	      States.H_Split.push_back(HamSplitTemp);
	    }
	}
    }

  //std::cout << States.bit_configs[5] << " " << States.bit_configs[6] << " " << States.bit_configs[7] << endl;

  std::cout << "Total Number of Slater Determinants = " << States.bit_configs.size() << endl << endl;

  return States;

}
  




J_Matrix_Elements Get_J_Matrix_Elements(Input_Parameters Parameters)
{
  J_Matrix_Elements J_ME;

  std::cout << "Reading J-Scheme Matrix Elements" << endl;
  std::cout << "--------------------------------" << endl;

  //Get Matrix Elements
  int NumElements;
  double OBME, TBME, shell1, shell2, shell3, shell4, coupJ, coupT; // J-scheme interaction file contents
  ifstream interaction;	// interaction file
  string interactionline; // interaction file line
  
  //open interaction file
  string fullpath2 = "C:\\Users\\Sam\\Documents\\Shell Model\\files\\" + Parameters.MatrixElements + ".int";
  interaction.open(fullpath2);
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
      J_ME.OBME.push_back(OBME);
    }

  vector<int> braket(4);
  J_ME.Braket.resize(NumElements);
  for(int i = 0; i < NumElements; ++i){ J_ME.Braket[i].resize(4); };
  J_ME.J.resize(NumElements);
  J_ME.T.resize(NumElements);
  J_ME.TBME.resize(NumElements);
  
  //read two-body parameters and two-body matrix elements
  double tempS1 = 0, tempS2 = 0, tempS3 = 0, tempS4 = 0, tempJ = -1, tempT = -1;
  int count = 0;
  for(int i = 0; i < NumElements; ++i)
    {
      getline(interaction, interactionline);
      istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> coupJ >> coupT >> TBME;	//	get two-body m.e.
      if (shell1 == tempS3 && shell2 == tempS4 && shell3 == tempS1 && shell4 == tempS2 && coupJ == tempJ && coupT == tempT)
	{ --NumElements; --i; continue; };
      tempS1 = shell1; tempS2 = shell2; tempS3 = shell3; tempS4 = shell4; tempJ = coupJ; tempT = coupT;
      braket[0] = shell1;										//	P1-4 = ordered proton/neutron shell
      braket[1] = shell2;
      braket[2] = shell3;
      braket[3] = shell4;
      J_ME.Braket[i] = braket;
      J_ME.J[i] = coupJ;								     //	J = ang momentum of coupled pair
      J_ME.T[i] = coupT;								     //	T = isospin of coupled pair
      J_ME.TBME[i] = TBME;
      //std::cout << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << "  " << coupJ << " " << coupT << " " << TBME << endl;
    }
  interaction.close();

  std::cout << "Number of One-Body Matrix Elements = " << J_ME.OBME.size() << endl;
  std::cout << "Number of J-Scheme Two-Body Matrix Elements = " << J_ME.TBME.size() << endl << endl;

  return J_ME;

};





M_Matrix_Elements Convert_To_M_Matrix_Elements(Model_Space Space, J_Matrix_Elements J_ME)
{
  M_Matrix_Elements M_ME;

  std::cout << "Converting Matrix Elements from J-Scheme to M-Scheme" << endl;
  std::cout << "----------------------------------------------------" << endl;

  M_ME.OBME = J_ME.OBME;

  //Change J-Scheme matrix elements to M-Scheme
  size_t length1 = Space.levelsind.size();
  size_t length2 = J_ME.TBME.size();
  int l2, k2, m2, n2;

  int size = length1*length1*(length1 - 1)*(length1 - 1)/8;
  int ind = 0;
  vector<int> braket(4);
  M_ME.Braket.resize(size);
  for(int i = 0; i < size; ++i){ M_ME.Braket[i].resize(4); };
  M_ME.TBME.resize(size);

  for (int m1 = 0; m1 < length1; ++m1)
    {
      for (int n1 = m1 + 1; n1 < length1; ++n1)
	{
	  for (int l1 = 0; l1 < length1; ++l1)
	    {
	      for (int k1 = l1 + 1; k1 < length1; ++k1)
		{
		  double angproj1 = Space.levelsm[m1] + Space.levelsm[n1];
		  double angproj2 = Space.levelsm[l1] + Space.levelsm[k1];
		  double isoproj1 = Space.levelst[m1] + Space.levelst[n1];
		  double isoproj2 = Space.levelst[l1] + Space.levelst[k1];
		  if (angproj1 == angproj2 && isoproj1 == isoproj2)
		    {
		      double tempmom1 = abs(Space.levelsj[m1] + Space.levelsj[n1]);
		      int size1 = int(tempmom1 - abs(Space.levelsj[m1] - Space.levelsj[n1]) + 1.1);
		      int ind1 = 0;
		      vector<double> angmom1(int(tempmom1 - abs(Space.levelsj[m1] - Space.levelsj[n1]) + 1.1)); 
		      while (tempmom1 >= abs(Space.levelsj[m1] - Space.levelsj[n1]) && tempmom1 >= abs(angproj1))
			{ 
			  angmom1[ind1] = tempmom1; ++ind1;
			  tempmom1 = tempmom1 - 1.0;
			};
		      angmom1.resize(ind1);

		      double tempmom2 = abs(Space.levelsj[l1] + Space.levelsj[k1]);
		      int size2 = int(tempmom2 - abs(Space.levelsj[l1] - Space.levelsj[k1]) + 1.1);
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
		      for (int i = 0; i < angmom1.size(); ++i)
			{
			  for (int j = 0; j < angmom2.size(); ++j)
			    {
			      if (angmom1[i] == angmom2[j]){ angmom3[ind3] = angmom1[i]; ++ind3; };
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

		      for (int b = 0; b < angmom3.size(); ++b)	//	loop over interaction file lines
			{
			  for (int c = 0; c < iso3.size(); ++c)	//	loop over relevant angular momentum
			    {
			      for (int a = 0; a < length2; ++a)	//	loop over relevant isospin
				{
				  double tbodyphase = 1.0;
				  double tbody1, tbody2, tbody3, tbody4;
				  if (J_ME.Braket[a][0] <= J_ME.Braket[a][1])
				    { tbody1 = J_ME.Braket[a][0]; tbody2 = J_ME.Braket[a][1]; }
				  else
				    { 
				      tbody1 = J_ME.Braket[a][1]; tbody2 = J_ME.Braket[a][0];
				      int power = int(Space.shellsj[tbody1 - 1] + Space.shellsj[tbody2 - 1] - angmom3[b] + 1.0 - iso3[c] + 0.01);
				      tbodyphase = -1.0 * tbodyphase * pow(-1.0, power); 
				    };
				  if (J_ME.Braket[a][2] <= J_ME.Braket[a][3])
				    { tbody3 = J_ME.Braket[a][2]; tbody4 = J_ME.Braket[a][3]; }
				  else
				    {
				      tbody3 = J_ME.Braket[a][3]; tbody4 = J_ME.Braket[a][2];
				      int power = int(Space.shellsj[tbody3 - 1] + Space.shellsj[tbody4 - 1] - angmom3[b] + 1.0 - iso3[c] + 0.01);
				      tbodyphase = -1.0 * tbodyphase * pow(-1.0, power);
				    };
				  if (((tbody1 == m2 && tbody2 == n2 && tbody3 == l2 && tbody4 == k2) || 
				       (tbody1 == l2 && tbody2 == k2 && tbody3 == m2 && tbody4 == n2)) && 
				       J_ME.J[a] == angmom3[b] && 
				       J_ME.T[a] == iso3[c])
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
		      if (Space.levelsj[m1] == Space.levelsj[n1]){ factor1 = 2.0; };
		      if (Space.levelsj[l1] == Space.levelsj[k1]){ factor2 = 2.0; };
		      factor3 = sqrt(factor1*factor2);
		      int flag = 0;
		      //don't include duplicates
		      for (int i = 0; i < M_ME.TBME.size(); ++i)
			{ 
			  if (M_ME.Braket[i][0] == l1 + 1 && M_ME.Braket[i][1] == k1 + 1 && 
			      M_ME.Braket[i][2] == m1 + 1 && M_ME.Braket[i][3] == n1 + 1 )
			    { 
			      flag = 1;
			    }; 
			};
		      if (abs(tempmatel) > 0.001 && flag == 0)
			{
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
    };

  M_ME.Braket.resize(ind);
  M_ME.TBME.resize(ind);

  std::cout << "Number of M-Scheme Two-Body Matrix Elements = " << M_ME.TBME.size() << endl << endl;

  return M_ME;
}





Hamiltonian Build_Hamiltonian(Model_Space Space, Many_Body_States States, M_Matrix_Elements M_ME)
{
  Hamiltonian Ham;

  std::cout << "Building Hamiltonian Matrix" << endl;
  std::cout << "---------------------------" << endl;

  Ham.H_Split = States.H_Split;

  Ham.Matrices.resize(Ham.H_Split.size());
  int totallength = 0;
  for (int n = 0; n < Ham.H_Split.size(); ++n)
    {
      totallength += (Ham.H_Split[n] * (Ham.H_Split[n] + 1) / 2);
      Ham.Matrices[n].resize(Ham.H_Split[n]*Ham.H_Split[n]);
    }
  int HamiltonianCount = 0;
  int indexoffset = 0;
  double perc = 0.0;

  //multiplicative strengths for all one and two-body matrix elements, respectively
  double strength1 = 1.0;
  double strength2 = 1.0;  //pow(18.0/(16.0 + PN),0.3);

  for (int n = 0; n < Ham.H_Split.size(); ++n)
    {
      int length = Ham.H_Split[n];
      for (int i = 0; i < length; ++i)
	{
	  int bra = States.bit_configs[indexoffset + i];
	  for (int j = i; j < length; ++j)
	    {
	      ++HamiltonianCount; int ket = States.bit_configs[indexoffset + j];
	      if (((HamiltonianCount * 100.0 / totallength) - perc) >= 1.0)
		{
		  perc = HamiltonianCount * 100.0 / totallength;
		  std::cerr << "Building Hamiltonian: " << int(HamiltonianCount * 100.0 / totallength) << " % \r";
		};
	      double c = matrixe(strength1, strength2, Space.levelschemeind, bra, ket, M_ME.OBME, M_ME.Braket, M_ME.TBME);
	      Ham.Matrices[n][i * length + j] = c;
	      if (i != j){ Ham.Matrices[n][j * length + i] = c; };
	    };
	};
      indexoffset += length;
    };

  std::cout << endl;

  return Ham;
}
  
  /*for (int n = 0; n < H.size(); ++n)
    {for(int i = 0; i < HamSplit[n]; ++i)
	{for(int j = 0; j < HamSplit[n]; ++j)
	    {std::cout << H[n][HamSplit[n]*i + j] << " ";};
	  std::cout << endl;
	};
      std::cout << endl;
    };
  
  double Htest1[] = { 4.0, 1.0, -2.0, 2.0,
		      1.0, 2.0, 0.0, 1.0,
		      -2.0, 0.0, 3.0, -2.0,
		      2.0, 1.0, -2.0, -1.0 };
  vector<double> Htest(Htest1, Htest1 + sizeof(Htest1) / sizeof(double));
  
  vector<vector<double>> tridiagonal1 = Lanczos(Htest);
  vector<vector<double>> tridiagonal2 = Householder(Htest);
  
  std::cout << endl << "alpha = ";
  for (int i = 0; i < tridiagonal1[0].size(); ++i)
    {
      std::cout << tridiagonal1[0][i] << " ";
    }
  std::cout << endl << "beta = ";
  for (int i = 0; i < tridiagonal1[1].size(); ++i)
    {
      std::cout << tridiagonal1[1][i] << " ";
    }
  std::cout << endl;
  std::cout << endl << "vector = " << endl;
  for (int i = 0; i < sqrt(tridiagonal1[2].size()); ++i)
    {for (int j = 0; j < sqrt(tridiagonal1[2].size()); ++j)
	{
	  std::cout << tridiagonal1[2][i*sqrt(tridiagonal1[2].size()) + j] << " ";
	}
      std::cout << endl;
    }
    
  vector<vector<double>> diagonal1 = QR2(tridiagonal1[0], tridiagonal1[1], tridiagonal1[2]);


  std::cout << endl << "alpha = ";
  for (int i = 0; i < tridiagonal2[0].size(); ++i)
    {
      std::cout << tridiagonal2[0][i] << " ";
    }
  std::cout << endl << "beta = ";
  for (int i = 0; i < tridiagonal2[1].size(); ++i)
    {
      std::cout << tridiagonal2[1][i] << " ";
    }
  std::cout << endl;
  std::cout << endl << "vector = " << endl;
  for (int i = 0; i < sqrt(tridiagonal2[2].size()); ++i)
    {for (int j = 0; j < sqrt(tridiagonal2[2].size()); ++j)
	{
	  std::cout << tridiagonal2[2][i*sqrt(tridiagonal2[2].size()) + j] << " ";
	}
      std::cout << endl;
    }
    
    vector<vector<double>> diagonal2 = QR2(tridiagonal2[0], tridiagonal2[1], tridiagonal2[2]);*/


Eigen_System Solve_Hamiltonian(Hamiltonian Ham)
{
  Eigen_System Eigen;

  int Hamtot = 0;
  for (int i = 0; i < Ham.H_Split.size(); ++i)
    { 
      Hamtot = Hamtot + Ham.H_Split[i];
    };

  Eigen.Energies.resize(Hamtot);
  Eigen.Vectors.resize(Hamtot);
  for(int i = 0; i < Hamtot; ++i)
    {
      Eigen.Vectors[i].resize(Hamtot);
    };

  /*for (int i = 0; i < Ham.H_Split.size(); ++i)
    {
      for (int j = 0; j < Ham.H_Split[i]; ++j)
	{
	  for (int k = 0; k < Ham.H_Split[i]; ++k)
	    {
	      std::cout << Ham.Matrices[i][Ham.H_Split[i]*j + k] << " ";
	    }
	  std::cout << endl;
	}
      std::cout << endl;
      }*/

  double perc1 = 0.0;
  int testcount1 = 0;
  int Hamtemp = 0;
  for (int i = 0; i < Ham.Matrices.size(); ++i)  //!!! Add Option to choose Lanczos or Householder !!!
    {
      vector<double> solutionstemp = solution(Ham.Matrices[i], Ham.H_Split[i]);
      //vector<vector<double> > tridiagonal1 = Lanczos(Ham.Matrices[i]);
      //vector<vector<double> > tridiagonal2 = Householder(Ham.Matrices[i]);
      //vector<double> solutionstemp = QR2(tridiagonal1[0], tridiagonal1[1], tridiagonal1[2]);
      //vector<double> solutionstemp = QR2(tridiagonal2[0], tridiagonal2[1], tridiagonal2[2]);
      ++testcount1;
      if (((testcount1 * 100.0 / Ham.Matrices.size()) - perc1) >= 1.0)
	{
	  perc1 = testcount1 * 100.0 / Ham.Matrices.size();
	  std::cerr << "Solving Hamiltonian: " << int(testcount1 * 100.0 / Ham.Matrices.size()) << " % \r";
	};
      for (int j = 0; j < solutionstemp.size()/(Ham.H_Split[i] + 1); ++j)
	{ 
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


  //Remove Degenerate States
  vector<double> EnergiesTemp(Hamtot);
  vector<int> flagvec(Hamtot);
  int flagcount = 0;
  vector<vector<double> > VectorsTemp(Hamtot);
  int ind = 0;
  for (int i = 0; i < Eigen.Energies.size(); ++i)
    {
      double temp1 = Eigen.Energies[i];
      for (int k = 0; k < flagcount; ++k)
	{
	  int flag2 = flagvec[k];
	  if(i == flag2){ goto contflag; };
	};
      EnergiesTemp[ind] = Eigen.Energies[i];
      VectorsTemp[ind] = Eigen.Vectors[i];
      ++ind;
      for (int j = i + 1; j < Eigen.Energies.size(); ++j)
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
  vector<double>().swap(EnergiesTemp);
  vector<vector<double> >().swap(VectorsTemp);
  vector<int>().swap(flagvec);
  Eigen.Energies.resize(ind);
  Eigen.Vectors.resize(ind);


  //Order States by Energy
  for (int i = 0; i < Eigen.Energies.size(); ++i)
    {
      double temp1 = Eigen.Energies[i];
      vector<double> temp2 = Eigen.Vectors[i]; 
      int index = i; //initialize
      for (int j = i + 1; j < Eigen.Energies.size(); ++j)
	{
	  double temp2 = Eigen.Energies[j];
	  if(temp1 < temp2){ continue; }
	  else{ temp1 = Eigen.Energies[j]; index = j; };
	};
      Eigen.Energies[index] = Eigen.Energies[i];
      Eigen.Vectors[index] = Eigen.Vectors[i];
      Eigen.Energies[i] = temp1;
      Eigen.Vectors[i] = temp2;
    };

  std::cout << endl;

  return Eigen;
}



Level_Structure Get_Level_Structure(Model_Space Space, Many_Body_States States, Eigen_System Eigen)
{
  Level_Structure Structure;
  Structure.Energies = Eigen.Energies;
  Structure.Vectors = Eigen.Vectors;

  size_t length3 = Eigen.Energies.size();
  double perc1 = 0.0;
  int testcount1 = 0;
  for (int i = 0; i < length3; ++i)
    {
      ++testcount1;
      vector<double> statebra = Eigen.Vectors[i];
      if (((testcount1 * 100.0 / length3) - perc1) >= 1.0)
	{
	  perc1 = testcount1 * 100.0 / length3;
	  std::cerr << "Calculating Angular Momentum: " << int(testcount1 * 100.0 / length3) << " % \r";
	};
      double stateJ = StateJ(testcount1, statebra, States.bit_configs, Space.levelschemeind, Space.levelsj, Space.levelsm, Space.shellsm);
      Structure.Js.push_back(stateJ);
    };
  std::cout << endl;
  
  double perc2 = 0.0;
  int testcount2 = 0;
  for (int i = 0; i < length3; ++i)
    {
      ++testcount2;
      vector<double> statebra = Eigen.Vectors[i];
      if (((testcount2 * 100.0 / length3) - perc2) >= 1.0)
	{
	  perc2 = testcount2 * 100.0 / length3;
	  std::cerr << "Calculating Isospin: " << int(testcount2 * 100.0 / length3) << " % \r";
	};
      double stateT = StateT(testcount2, statebra, States.bit_configs, Space.levelsj, Space.levelsm, Space.levelst);
      Structure.Ts.push_back(stateT);
    };
  std::cout << endl;
  
  double perc3 = 0.0;
  int testcount3 = 0;
  for (int i = 0; i < length3; ++i)
    {
      ++testcount3;
      vector<double> statebra = Eigen.Vectors[i];
      if (((testcount3 * 100.0 / length3) - perc3) >= 1.0)
	{ 
	  perc3 = testcount3 * 100.0 / length3;
	  std::cerr << "Calculating Parity: " << int(testcount3 * 100.0 / length3) << " % \r";
	};
      double stateP = StateP(statebra, States.bit_configs, Space.levelsl);
      Structure.Ps.push_back(stateP);
    };

  std::cout << endl;

  return Structure;
}



void Results(Input_Parameters Parameters, Level_Structure Structure)
{ 
  stringstream a, b;
  a << Parameters.P;
  b << Parameters.N;
  string filenum1 = a.str(), filenum2 = b.str();
  
  ofstream myfile;
  myfile.open("C:\\Users\\Sam\\Documents\\Shell Model\\files\\" + 
	      Parameters.LevelScheme + "." + Parameters.MatrixElements + "." + 
	      filenum1 + "." + filenum2 + ".solution.txt");
  for (int i = 0; i < Structure.Energies.size(); ++i)
    {
      myfile << Structure.Energies[i] << "     \t" << Structure.Js[i] << "\t" << Structure.Ps[i] << "\t" << Structure.Ts[i] << "\n";
    };
    myfile.close();
}


int main(int argc, char* argv[])
{ 
  clock_t t1, t2;    	//initialize program clock
  t1 = clock();


  string inputfile = "input.dat";
  Input_Parameters Parameters = Get_Input_Parameters(inputfile);
  Model_Space Space = Build_Model_Space(Parameters);
  Many_Body_States States = Build_Many_Body_States(Parameters, Space);
  J_Matrix_Elements J_ME = Get_J_Matrix_Elements(Parameters);
  M_Matrix_Elements M_ME = Convert_To_M_Matrix_Elements(Space, J_ME);

  if (Parameters.beta == 0)
    {
      Hamiltonian Ham = Build_Hamiltonian(Space, States, M_ME);
      Eigen_System Eigen = Solve_Hamiltonian(Ham);
      Level_Structure Structure = Get_Level_Structure(Space, States, Eigen);
      Results(Parameters, Structure);
    }
  else
    {
      Input_Parameters Parameters2 = Parameters;
      Parameters2.P = Parameters2.P + 1;
      Parameters2.N = Parameters2.N - 1;
      Many_Body_States States2 = Build_Many_Body_States(Parameters2, Space);
      Hamiltonian Ham2 = Build_Hamiltonian(Space, States, M_ME);
      Eigen_System Eigen2 = Solve_Hamiltonian(Ham);
      Level_Structure Structure2 = Get_Level_Structure(Space, States, Eigen);
    }

  double betam = BetaMDecay(Structure.Vectors[0], Structure1.Vectors[0], States.bit_configs, States2.bit_configs, Space.indp, Space.indn, Space.levelsl, Space.levelsj, Space.levelsm);
  
  std::cout << "Beta M = " << betam;

  t2 = clock();
  float diff((float)t2 - (float)t1);
  std::cout << diff / CLOCKS_PER_SEC << "sec" << endl;

  int a;
  std::cin >> a;

  return 0;

}
