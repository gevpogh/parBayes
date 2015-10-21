#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "BA_EoS_F.h"
#include <mpi.h> 

using namespace std;

//char* EOSNAME = "EOS\\HD_H_-HDD_s.dat";
char* EOSNAME = "EOS/EoS-APR_Vex_r0=0.4.dat";
char statfile[255];
int dbg_infos=1;	//Debugging switch to show details about MPI communication

ofstream Kakerlake;
ofstream NSCDB; // Neutron star configuration database

EoS_F eosf; 

const double N1 = 11;
const double N2 = 11;
const double N3 = 11;

/*****************************************/
const double C1 = 1.11269e-5; //[4Pi/c^2]
const double C2 = 1.4766; //[]
const double m_N = 938.272046;// Proton mass [MeV] //931.494; //[MeV]
//const double n_B = 1.; //
/*****************************************/

/*****************************************/
const double gama_min = 0.0;
const double gama_max = 1.0;

const double Cs_sqr_min = 0.3;
const double Cs_sqr_max = 1.0;

const double eps_crit_min = 400.;
const double eps_crit_max = 1000.;

const double eps_low = 172.2075;// = exp(5.1487)
const double eps_high = 201.5424;// = exp(5.306)

/*****************************************/

/*****************************************/
double Cs_sqr;// = 0.9; // Cs_sqr is in [0.3 ... 1]
double eps_crit;// = 800.; // eps_crit is in [100 ... 1000]
double gama;// = 0.2; // gama is in [0 ... 1]
double eps_delta;// = gama * eps_crit;
double eps_quark;// = eps_crit + eps_delta;
/*****************************************/

vector< double > m, m_B, p, r, e;
vector< double > M, M_B, R, Eps_0;

vector< vector< double > > M_vec, M_B_vec, R_vec;
//vector< double > Eps_0_vec;

double Log_P_low(double e_)
{
	double x = log(e_);
	return -6.3135 + 0.88656*x - 0.014317*pow(x,2) + 0.016476*pow(x,3) + 0.0022707*pow(x,4) + 0.000010627*pow(x,5) - 0.000014646*pow(x,6) - 9.0534e-7*pow(x,7) - 1.6578e-8*pow(x,8);
}

double Log_P_high(double e_)
{
	double x = log(e_);
	return 152.21 - 97.574 * x + 22.436 * x*x - 2.1846 * x*x*x + 0.077895 * x*x*x*x;
}

double Log_P_middle(double e_)
{
/*
	const double v1 = -0.00108343;
	const double v2 = 0.0187079;
	const double v3 = -0.101814;
	const double v4 = 0.185342;
	const double v5 = 0.136305;
	const double v6 = 0.648336;
	const double v7 = -6.73771;
	double x = log(e_);
	return v1*pow(x, 6) + v2*pow(x, 5) + v3*pow(x, 4) +
	       v4*pow(x, 3) + v5*pow(x, 2) + v6*x + v7;
*/
	double tga = Log_P_high(eps_high);
	tga-= Log_P_low(eps_low);
	tga/= eps_high - eps_low;
	return tga*(e_ - eps_high) + Log_P_high(eps_high);
}

double P_APR_cat(double e_)
{
//	exp(5.1487) = 172.2075
//	exp(5.306)  = 201.5424
	double result;
	if(e_<eps_low) result = exp(Log_P_low(e_));
	if(e_>=eps_low && (e_<=eps_high)) result = exp(Log_P_middle(e_));
	if(e_>eps_high) result = exp(Log_P_high(e_));
	return result;
}

double P_hadronic(double e_)
{
#ifdef APR
	return P_APR_cat(e_);
#else
	return eosf.get_P(e_);
#endif
}

double P_quark(double e_)
{
	return Cs_sqr*(e_ - eps_quark) + P_hadronic(eps_crit);
//	Cs_sqr*e_  - Cs_sqr*eps_bar + Pcrit
//	Cs_sqr*e_  + B*(Cs_sqr + 1)
//	B = 
}

double P(double e_)
{
	double result;
	
	eps_delta = gama * eps_crit;
	eps_quark = eps_crit + eps_delta;
	
	if(e_<=eps_crit) result = P_hadronic(e_);
	if((e_>eps_crit) && (e_<eps_quark)) result = P_hadronic(eps_crit);
	if(e_>=eps_quark) result = P_quark(e_);
	
	return result;
//	return P_quark(e_);
}

double dPde(double e_)
{
	double result;
	double h_e = 1e-5;
	if(e_ - h_e < 0.)
	{
		result = P(e_+h_e);
		result-= P(e_);
		result/= h_e;
    }
	else
	{
		result = P(e_+h_e);
		result-= P(e_-h_e);
		result/= h_e;
		result/= 2.;
	}
	return result;
}

double Log_n_B_APR_low(double e_)
{
	double v1,v2,v3,v4,v5,v6,v7;
	double x = log(e_);
	v1 =-3.28711e-10;
	v2 = -4.99547e-8;
	v3 = -2.37048e-6;
	v4 =-0.0000512534;
	v5 = -0.000589898;
	v6 = 0.996108;
	v7 = -6.84956;
	return v1*pow(x, 6) + v2*pow(x, 5) + v3*pow(x, 4) +
		   v4*pow(x, 3) + v5*pow(x, 2) + v6*x + v7;
}

double Log_n_B_APR_middle(double e_)
{
	double v1,v2,v3,v4,v5,v6,v7;
	double x = log(e_);
	v1 = -0.0000215335;
	v2 = 0.00024168;
	v3 = -0.000863071;
	v4 = 0.000754409;
	v5 = 0.000945294;
	v6 = 0.9978;
	v7 = -6.84706;
	return v1*pow(x, 6) + v2*pow(x, 5) + v3*pow(x, 4) +
		   v4*pow(x, 3) + v5*pow(x, 2) + v6*x + v7;
}

double Log_n_B_APR_high(double e_)
{
	double v1,v2,v3,v4,v5,v6,v7;
	double x = log(e_);
	v1 = -0.0019048;
	v2 = 0.0764515;
	v3 = -1.25998;
	v4 = 10.8996;
	v5 = -52.2216;
	v6 = 132.511;
	v7 = -142.987;
	return v1*pow(x, 6) + v2*pow(x, 5) + v3*pow(x, 4) +
		   v4*pow(x, 3) + v5*pow(x, 2) + v6*x + v7;
}

double n_B_APR(double e_)
{
#ifdef APR
	double x = log(e_);
	double result;
	if(x < -1.19122) result = exp(Log_n_B_APR_low(e_)); // APR_Cat
	if((x >= -1.19122) && (x < 5.306)) result = exp(Log_n_B_APR_middle(e_)); // APR_Cat
	if(x >= 5.306) result = exp(Log_n_B_APR_high(e_)); // APR_Cat
	return result;
#else
return eosf.get_n_B(e_);
#endif
}

double n_B_trans(double e_)
{
	double P_crit = P(eps_crit);
	double n_B_APR_crit = n_B_APR(eps_crit);
	double result = P_crit + e_;
	result /= P_crit + eps_crit;
	result *= n_B_APR_crit;
	return result;
}

double n_B_quark(double e_)
{
	eps_delta = gama * eps_crit;
	eps_quark = eps_crit + eps_delta;
	double P_crit = P(eps_crit);
	double B  = Cs_sqr*eps_quark - P_crit;
		   B /= 1. + Cs_sqr;
	double n_B_APR_crit = n_B_APR(eps_crit);

	double result  = e_ - B;
		   result /= eps_quark - B;
		   result  = pow(result, 1./(1.+Cs_sqr));
		   result *= n_B_APR_crit;
		   result *= P_crit + eps_quark;
		   result /= P_crit + eps_crit;

	return result;
}

double n_B(double e_)
{
	eps_delta = gama * eps_crit;
	eps_quark = eps_crit + eps_delta;
	double result;
	if(e_ < eps_crit) result = n_B_APR(e_); // APR_Cat
	if((e_ >= eps_crit) && (e_ <= eps_quark)) result = n_B_trans(e_); // Transition
	if(e_ > eps_quark) result = n_B_quark(e_); // APR_Cat
	return result;
}

double mu_B(double e_)
{
	return (P(e_)+e_)/n_B(e_);
}

//m'=f1(r,e,m,m_B,p)
double f1(double e_, double r_, double m_, double m_B_, double p_)
{
	double result = C1;
	result *= e_;
	result *= r_;
	result *= r_;
	return result;
//	return C1_*e_*r_*r_;
}

//mb'=f2(r,eps,m,m_B,p)
double f2(double e_, double r_, double m_, double m_B_, double p_)
{
//	return 0.;
	double result = C1;
	result *= n_B(e_);
	result *= m_N;
	result *= r_;
	result *= r_;
	if(1. <= 2. * C2 * m_/r_)
	{
		cout <<"Error 1!"<<endl;
		cout <<"m = "<< m_ << endl;
		cout <<"r = "<< r_ << endl;
		cout <<"e = "<< e_ << endl;
		cout <<"m_B = "<< m_B_ << endl;
		cout <<"P(e) = "<< p_ << endl;
		cout <<"2. * C2 * m/r  = "<< 2. * C2 * m_/r_ << endl;
		exit (EXIT_FAILURE);
	}
	result /= sqrt(1. - 2. * C2 * m_/r_);
	return result;
//	return C1 * n_B * m_N * r_ * r_ / sqrt(1. - 2. * C2 * m_/r_);
}

//p'=f3(r,e,m,m_B,p)
double f3(double e_, double r_, double m_, double m_B_, double p_)
{
	double result = C1;
	result *= r_;
	result *= r_;
	result *= r_;
	result *= p_;
	result += m_;
	result *= p_ + e_;
	result *= C2;
	if(r_ == 2. * C2 * m_)
	{
		cout <<"Error 2!"<<endl;
		exit (EXIT_FAILURE);
	}
	result /= r_ - 2. * C2 * m_;
	result /= r_;
	result *= -1.;
	return result;
//	return -C2 * (p_ + e_) * (m_ + C1 * p_*r_*r_*r_) / (r_*(r_ - 2. * C2 * m_));
}

//dr/deps = g3(eps,r,m,m_B,p)
double g3(double e_, double m_, double m_B_, double r_)
{
	double result = dPde(e_);
	double tmp = f3(e_, r_, m_, m_B_, P(e_));
	if(tmp == 0.)
	{
		cout <<"Error 3!"<<endl;
		exit (EXIT_FAILURE);
	}
	result /= tmp;
	return result;
//	return P'(e_) / f3(r_, m_, m_B_, p_);
}

//dm_B/deps = g2(eps,r,m,m_B,p)
double g2(double e_, double m_, double m_B_, double r_)
{
//	return 0.;
	double result = f2(e_, r_, m_, m_B_, P(e_));
	result *= g3(e_, m_, m_B_, r_);
	return result;
//	return f2(e_, r_, m_, m_B_, p_) * g3(e_, m_, m_B_, r_);
}

//dm/deps = g1(r,m,m_B,p)
double g1(double e_, double m_, double m_B_, double r_)
{
	double result = f1(e_, r_, m_, m_B_, P(e_));
	result *= g3(e_, m_, m_B_, r_);
	return result;
//	return f1(e_, r_, m_, m_B_, p_) * g3(e_, m_, m_B_, r_);
}

void Runge_Kutta(double t_0, double dt, double x_0, double y_0, double z_0)
{
	double k11,k12,k13,k14;
	double k21,k22,k23,k24;
	double k31,k32,k33,k34;

	double X = x_0;
	double Y = y_0;
	double Z = z_0;

	m.push_back(X);
	m_B.push_back(Y);
	r.push_back(Z);
	e.push_back(t_0);
	p.push_back(P(t_0));

	double t = t_0;
	do{
		k11 = g1(t, X, Y, Z)*dt;
		k21 = g2(t, X, Y, Z)*dt;
		k31 = g3(t, X, Y, Z)*dt;

		k12 = g1(t+dt/2., X+k11/2., Y+k21/2., Z+k31/2.)*dt;
		k22 = g2(t+dt/2., X+k11/2., Y+k21/2., Z+k31/2.)*dt;
		k32 = g3(t+dt/2., X+k11/2., Y+k21/2., Z+k31/2.)*dt;

		k13 = g1(t+dt/2., X+k12/2., Y+k22/2., Z+k32/2.)*dt;
		k23 = g2(t+dt/2., X+k12/2., Y+k22/2., Z+k32/2.)*dt;
		k33 = g3(t+dt/2., X+k12/2., Y+k22/2., Z+k32/2.)*dt;

		k14 = g1(t+dt, X+k13, Y+k23, Z+k33)*dt;
		k24 = g2(t+dt, X+k13, Y+k23, Z+k33)*dt;
		k34 = g3(t+dt, X+k13, Y+k23, Z+k33)*dt;

		X = X + (k11 + 2.*k12 + 2.*k13 + k14)/6.;
		Y = Y + (k21 + 2.*k22 + 2.*k23 + k24)/6.;
		Z = Z + (k31 + 2.*k32 + 2.*k33 + k34)/6.;

		t += dt;

		m.push_back(X);
		m_B.push_back(Y);
		r.push_back(Z);
		p.push_back(P(t));
		e.push_back(t);
		
/*		if (dbg_infos){	
			FILE* file =fopen(statfile,"a");
			fprintf(file,"\n[%f] RungeKutta with m#%d=%g mB#%d=%g r#%d=%g p#%d=%g e#%d=%g",MPI_Wtime(),m.size(),m[m.size()-1],m_B.size(),m_B[m_B.size()-1],r.size(),r[r.size()-1],p.size(),p[p.size()-1],e.size(),e[e.size()-1]);
			fprintf(file,"\n[%f] RungeKutta with m#X=%g mB#Y=%g r#Z=%g p#P(t)=%g e#t=%g",MPI_Wtime(),m[m.size()-1],m_B[m_B.size()-1],r[r.size()-1],p[p.size()-1],e[e.size()-1]);
			fclose(file);
		}
*/
	}while(t-fabs(dt) > 0.2);
//	cout << t << endl;
//	cout << endl;
}

void TOV_solver(double eps_0)
{
	double r_0 = 0.05;
	double m_0 = C1 * eps_0 * r_0 * r_0 * r_0 / 3.;
	double m_B_0 = m_0;
//	double p_0 = P(eps_0);
	double h_e = 0.05;
//	double h_r = 0.05;
	Runge_Kutta(eps_0,-h_e, m_0, m_B_0, r_0);
}

void DeleteMemorySpace()
{
	m.clear();
	m_B.clear();
	p.clear();
	r.clear();
	e.clear();
}

//---------------------------------------------------------------------------
#define FINISHED 0
#define ACONFIG 1
	static double eps_0_min = 100.;
	static double eps_0_max = 2000.;
	double h_eps_0 = 10.;
	int rank, numtasks;
	
void calcFamilyConfig()
{	double eps_0 = eps_0_min;
	
	int nStars = int((eps_0_max - eps_0_min)/h_eps_0) + 10;

	double datablock[5];
	eps_0 += rank*h_eps_0;
	for(int iStar = rank; iStar < nStars; iStar+=(numtasks-1))
	{
		if (dbg_infos){	
			char numstr[9];
			strcpy(statfile,"mpiid-");
			sprintf(numstr, "%d", rank);
			strcat(statfile,numstr);
			strcat(statfile,"-runprotocol.txt");
			FILE* file =fopen(statfile,"a");
			fprintf(file,"\n[%f] TOV:Starts with #%d from %d with eps_0 =%g",MPI_Wtime(),iStar,nStars,eps_0);
			fclose(file);
		}
		
		TOV_solver(eps_0);

		datablock[0] = (double)iStar;
		datablock[1] = eps_0;
		datablock[2] = m[m.size()-1];
		datablock[3] = m_B[m_B.size()-1];
		datablock[4] = r[r.size()-1];
		DeleteMemorySpace();

		if (dbg_infos){	
			FILE* file =fopen(statfile,"a");
			fprintf(file,"\n[%f] TOV:Ended with #%d values for iStar eps_0 m mB r",MPI_Wtime(),iStar);
			for (int i=0;i<5;i++) fprintf(file," %g",datablock[i]);
			fclose(file);
		}
		
		MPI_Send(datablock,5,MPI_DOUBLE,0,ACONFIG,MPI_COMM_WORLD);
		
		if(eps_0 >= eps_0_max) break;
		eps_0 += (numtasks-1)*h_eps_0;
	}
//Send FINISHED tag when all confings of particular slave are are finisched
	MPI_Send(datablock,5,MPI_DOUBLE,0,FINISHED,MPI_COMM_WORLD);

/*	M.clear();
	M_B.clear();
	R.clear();
	Eps_0.clear();
*/
}
//---------------------------------------------------------------------------

int calcConfigs()
{
double parblock[5];

if (rank == 0) {
	
	double hCs_sqr = (Cs_sqr_max - Cs_sqr_min) / (N1-1);
	double hGama = (gama_max - gama_min) / (N2-1);
	double hEps_crit = (eps_crit_max - eps_crit_min) / (N3-1);

	for(int iPar1 = 0; iPar1 < N1; iPar1++)
		for(int iPar2 = 0; iPar2 < N2; iPar2++)
			for(int iPar3 = 0; iPar3 < N3; iPar3++)
			{int slavecount = 0;
				parblock[0] = Cs_sqr = Cs_sqr_min + iPar1*hCs_sqr;
				parblock[1] = gama = gama_min + iPar2*hGama;
				parblock[2] = eps_crit = eps_crit_min + iPar3*hEps_crit;
				parblock[3] = dbg_infos;

				eps_delta = gama * eps_crit;
				eps_quark = eps_crit + eps_delta;

				cout <<eps_crit<<" ";
				cout <<gama<<" ";
				cout <<Cs_sqr<<endl;
				
				MPI_Bcast(parblock,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
				if (dbg_infos){
					cout <<"["<<MPI_Wtime()<<"]"<<" Broadcasted eps_crit="<<parblock[2]<<" gamma="<<parblock[1]<<" Cs_sqr="<<parblock[0]<< endl;
				}

				NSCDB <<"#MOD" <<"	";
				NSCDB << eps_crit <<"	";
				NSCDB << gama <<"	";
				NSCDB << Cs_sqr << endl;
				
				do{
				MPI_Status status;
				MPI_Recv(parblock,5,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if ( status.MPI_TAG == ACONFIG ) {
					int iStar = (int)parblock[0];
					if ( iStar > Eps_0.size() ){
						Eps_0.resize(iStar);
						M.resize(iStar);
						M_B.resize(iStar);
						R.resize(iStar);
					}
					Eps_0[iStar]=parblock[1];
					M[iStar]=parblock[2];
					M_B[iStar]=parblock[3];
					R[iStar]=parblock[4];
				if (dbg_infos){
					cout <<"["<<MPI_Wtime()<<"]"<<" Slave#"<<status.MPI_SOURCE<<" finished his Config eps0"<<Eps_0[iStar]<<" M="<<M[iStar]<<" R="<<R[iStar]<< endl;
				}
				}
				if ( status.MPI_TAG == FINISHED ) {
					slavecount++;
					if (dbg_infos){
					cout <<"["<<MPI_Wtime()<<"]"<<" Slave#"<<status.MPI_SOURCE<<" finished his part for iPar1="<<iPar1<<" iPar2="<<iPar2<<" iPar3="<<iPar3<<" "<<(numtasks-slavecount)<<"remaining"<<endl;
					}
				}
				} while(slavecount < numtasks - 1);

				//////  Write down configuration to the file using standart format                 
				int nStars = int((eps_0_max - eps_0_min)/h_eps_0) + 10;
				double delta_e = 1.;
				double delta_m = 1.;
				for (int iStar = 0 ; iStar < nStars; iStar ++) //writeout matrix into file
				{
				if(iStar != 0)
				{
					delta_e = Eps_0[iStar] - Eps_0[iStar-1];
					delta_m = M[iStar] - M[iStar-1];
					cout << delta_m << endl;
				}
				if((delta_m / delta_e) >= 0.) // Zel'dovich stability condition
				{
					NSCDB << Eps_0[iStar] <<"	";
					NSCDB << R[iStar] <<"	";
					NSCDB << M[iStar] <<"	";
					NSCDB << M_B[iStar] << endl;
					if(iStar != 0)
						if(M[iStar] - M[iStar-1] < 0.) cout << (M[iStar] - M[iStar-1]) << endl;
				}
				}
			}
		if (dbg_infos){
			cout <<"["<<MPI_Wtime()<<"]"<<" Master Finished!"<<endl;
		}
		parblock[1] = 12345678.9 ;
		MPI_Bcast(parblock,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
		NSCDB.close();

	}
	else {
		do {
			MPI_Bcast(parblock,5,MPI_DOUBLE,0,MPI_COMM_WORLD);

			Cs_sqr = parblock[0];
			gama =	parblock[1];
			eps_crit = parblock[2];
			dbg_infos = parblock[3];
			eps_delta = gama * eps_crit;
			eps_quark = eps_crit + eps_delta;

			if (dbg_infos){
				cout <<"["<<MPI_Wtime()<<"]"<<" Slave #"<<rank<<"Got Broadcast with eps_crit="<<parblock[2]<<" gamma="<<parblock[1]<<" Cs_sqr="<<parblock[0]<< endl;
			}

			if ( gama == 12345678.9 ) break;
		
			calcFamilyConfig();
			if (dbg_infos){
				cout <<"["<<MPI_Wtime()<<"]"<<" Slave #"<<rank<<" Finished his part of TOV!"<<endl;
			}
		} while (eps_crit > 0);
	}
	return 0;
}

void PrintTable(int N, vector<double> A, vector<double>  B, vector<double>  C, vector<double>  D)
{
	for(int i = 0; i < A.size(); i++)
		NSCDB << A[i] <<"	"<< B[i] <<"	"<< C[i] <<"	"<< D[i] << endl;
}


int main(int argc, char *argv[])
{	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 
	eosf.Read_EoS_FILE_(EOSNAME);
	if (rank == 0)
	{
	//Reading the Command line Arguments                  
		if(argc<1)
		{
			//Give default values to the parameters
			printf("Use D as command line argument for debugging mode, if you need to generating detailed information about parallel computation.\n");
			dbg_infos=0;  
		}
		else
		{
			int i;
			printf("Input: arguments read from the command line are ");
			for(i=1;i<argc;i++) printf("%s ",argv[i]);

			if(!strcmp(argv[1],"D")||!strcmp(argv[1],"d")) 
				{ 
					dbg_infos=1; 
					printf("Debugging mode will be used to protocol steps done by MPI-slaves. Check the \"mpiid-$ID$-runprotocol.txt\" , where $ID$ is the unique identification of parallel running slave rank\n");
				}
		}
		ofstream EoS;
		EoS.open("EoS.dat");
		for(double iE = 10; iE <= 2000.; iE+=10.)
		{
			EoS << iE <<"	"<< P_hadronic(iE);
			EoS << "	"<< n_B_APR(iE);
			EoS << "	"<< eosf.get_mu(iE)
			<< endl;
		}
	//		EoS << iE <<"	"<< P(iE) <<"	"<< dPde(iE) << endl;
		EoS.close();

		NSCDB.open("NSCDB_EoS_APR_exv.dat");
		NSCDB << "#EoS hybrid (based on APR_exv)" << endl;
		
	}
	calcConfigs();
	MPI_Finalize();
	return 0;
}