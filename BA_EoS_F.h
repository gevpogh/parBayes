
#pragma once

#define peos 1510

class EoS_F 
{
public:
	EoS_F(void);
	~EoS_F(void);

	int EOSTABFILEsize;
	double mu_crit
		//,	mu_min
		;
double	get_n_B(double e_);
double	get_P(double e_);
double 	get_mu(double e_);
double epsf[peos],Pf[peos],nf[peos],muf[peos];
bool m_file;
char name[100];

void Read_EoS_FILE_(char*);
};
