// Code for pentavalent model
// Created by Andromachi Karachaliou Prasinou
// Date 28/05/2022



#include <Rcpp.h>
using namespace Rcpp;



// Helper functions for 2,3 matrix indexing
inline int index2(int i, int a, int y)
{
    return y*i + a;
}


// [[Rcpp::export]]


List pentavalent(double t, NumericVector x, NumericVector params, NumericMatrix Cm, NumericVector mu, double max_age, NumericVector B, NumericVector ir_A, NumericVector ir_R, NumericVector stoc, NumericVector epi_A, NumericVector epi_R)
{
    
    NumericVector dxdt(x.length());
    
    double phi = params[0]; //loss of immunity
    double rho = params[1]; // recovery of disease
    double alpha_A = params[2]; // duration of carriage menA
	double alpha_R = params[3]; // duration of carriage menCWYX
    double epsbeta = params[4]; // seasonal amplitude
    double beta_A = params[5]; // contact parameter
	double beta_R = params[6];
	double ksi = params[7];
	double delta = params[8];
	double ws = params[9];
	double wl = params[10];
	double delta_R = params[11];
	double ksi_R = params[12];
    
    double pi = M_PI;
    
    double age = 1.0/365.0; // Daily Ageing
    
    int comp=68;
    
    enum  state_variables {S,C1A,C1R,I1A,I1R,RA,RR,C2A,C2R,I2A,I2R,R,SVlong,C1AVlong,C1RVlong,I1AVlong,I1RVlong,RAVlong,RRVlong,C2AVlong,C2RVlong,I2AVlong,I2RVlong,RVlong,SVshort,C1AVshort,C1RVshort,I1AVshort,I1RVshort,RAVshort,RRVshort,C2AVshort,C2RVshort,I2AVshort,I2RVshort,RVshort,SPlong,C1APlong,C1RPlong,I1APlong,I1RPlong,RAPlong,RRPlong,C2APlong,C2RPlong,I2APlong,I2RPlong,RPlong,SPshort,C1APshort,C1RPshort,I1APshort,I1RPshort,RAPshort,RRPshort,C2APshort,C2RPshort,I2APshort,I2RPshort,RPshort,IncA,IncAV,IncAP,IncA_all,IncR,IncRV,IncRP,IncR_all};
    
    NumericVector FOI_A(max_age);
	NumericVector FOI_R(max_age);
    
    double stoc1 = stoc[t/365];
    
    double birth = B[t/365];
    
    double N = 0.0;
    for( int j=0; j<(comp-8)*max_age; j++)
    {
        N += x[j];
    }
    
    // Calculate force of infection
    for(int a=0 ; a<max_age ; a++ )
    {
        FOI_A[a] = 0.0;
		FOI_R[a] = 0.0;
		
        for( int a2=0; a2<max_age; a2++ )
        {
            
            double Z = (1+epsbeta*cos(2*pi*t/365.0));
			
            FOI_A[a] += beta_A*Z*Cm(a,a2)*(x[index2(C1A,a2,max_age)] + x[index2(C2A,a2,max_age)] + x[index2(I1A,a2,max_age)] + x[index2(I2A,a2,max_age)] 
			+ x[index2(C1AVlong,a2,max_age)] + x[index2(C2AVlong,a2,max_age)] + x[index2(I1AVlong,a2,max_age)] + x[index2(I2AVlong,a2,max_age)]
			+ x[index2(C1AVshort,a2,max_age)] + x[index2(C2AVshort,a2,max_age)] + x[index2(I1AVshort,a2,max_age)] + x[index2(I2AVshort,a2,max_age)]
			+ x[index2(C1APlong,a2,max_age)] + x[index2(C2APlong,a2,max_age)] + x[index2(I1APlong,a2,max_age)] + x[index2(I2APlong,a2,max_age)]
			+ x[index2(C1APshort,a2,max_age)] + x[index2(C2APshort,a2,max_age)] + x[index2(I1APshort,a2,max_age)] + x[index2(I2APshort,a2,max_age)])/N;
			
			FOI_R[a] += beta_R*Z*Cm(a,a2)*(x[index2(C1R,a2,max_age)] + x[index2(C2R,a2,max_age)] + x[index2(I1R,a2,max_age)] + x[index2(I2R,a2,max_age)] 
			+ x[index2(C1RVlong,a2,max_age)] + x[index2(C2RVlong,a2,max_age)] + x[index2(I1RVlong,a2,max_age)] + x[index2(I2RVlong,a2,max_age)]
			+ x[index2(C1RVshort,a2,max_age)] + x[index2(C2RVshort,a2,max_age)] + x[index2(I1RVshort,a2,max_age)] + x[index2(I2RVshort,a2,max_age)]
			+ x[index2(C1RPlong,a2,max_age)] + x[index2(C2RPlong,a2,max_age)] + x[index2(I1RPlong,a2,max_age)] + x[index2(I2RPlong,a2,max_age)]
			+ x[index2(C1RPshort,a2,max_age)] + x[index2(C2RPshort,a2,max_age)] + x[index2(I1RPshort,a2,max_age)] + x[index2(I2RPshort,a2,max_age)])/N;
        }
    }
    
    
    
    dxdt[ index2(S , 0, max_age )] = birth*N + phi*x[index2(RA,0,max_age)] + phi*x[index2(RR,0,max_age)] + phi*x[index2(R,0,max_age)] - stoc1*FOI_A[0]*x[index2(S,0,max_age)] - stoc1*FOI_R[0]*x[index2(S,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(S,0,max_age)] + ws*x[index2(SVshort,0,max_age)] + ws*x[index2(SPshort,0,max_age)] + wl*x[index2(SVlong,0,max_age)] + wl*x[index2(SPlong,0,max_age)];
    
    dxdt[ index2(C1A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1A,0,max_age)] + stoc1*FOI_A[0]*x[index2(S,0,max_age)]
	+ ws*x[index2(C1AVshort,0,max_age)] + ws*x[index2(C1APshort,0,max_age)] + wl*x[index2(C1AVlong,0,max_age)] + wl*x[index2(C1APlong,0,max_age)];
	
	dxdt[ index2(C1R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1R,0,max_age)] + stoc1*FOI_R[0]*x[index2(S,0,max_age)]
	+ ws*x[index2(C1RVshort,0,max_age)] + ws*x[index2(C1RPshort,0,max_age)] + wl*x[index2(C1RVlong,0,max_age)] + wl*x[index2(C1RPlong,0,max_age)];
    
    dxdt[ index2(I1A , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1A,0,max_age)];
	
	dxdt[ index2(I1R , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1R,0,max_age)];
    
    dxdt[ index2(RA , 0, max_age )] = rho*x[index2(I1A,0,max_age)]+alpha_A*x[index2(C1A,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RA,0,max_age)] - stoc1*FOI_R[0]*x[index2(RA,0,max_age)]
	+ ws*x[index2(RAVshort,0,max_age)] + ws*x[index2(RAPshort,0,max_age)] + wl*x[index2(RAVlong,0,max_age)] + wl*x[index2(RAPlong,0,max_age)];
    
	dxdt[ index2(RR , 0, max_age )] = rho*x[index2(I1R,0,max_age)]+alpha_R*x[index2(C1R,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RR,0,max_age)] - stoc1*FOI_A[0]*x[index2(RR,0,max_age)]
	+ ws*x[index2(RRVshort,0,max_age)] + ws*x[index2(RRPshort,0,max_age)] + wl*x[index2(RRVlong,0,max_age)] + wl*x[index2(RRPlong,0,max_age)];
	
	dxdt[ index2(C2A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2A,0,max_age)] + stoc1*FOI_A[0]*x[index2(RR,0,max_age)]
	+ ws*x[index2(C2AVshort,0,max_age)] + ws*x[index2(C2APshort,0,max_age)] + wl*x[index2(C2AVlong,0,max_age)] + wl*x[index2(C2APlong,0,max_age)];
	
	dxdt[ index2(C2R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2R,0,max_age)] + stoc1*FOI_R[0]*x[index2(RA,0,max_age)]
	+ ws*x[index2(C2RVshort,0,max_age)] + ws*x[index2(C2RPshort,0,max_age)] + wl*x[index2(C2RVlong,0,max_age)] + wl*x[index2(C2RPlong,0,max_age)];
	
	dxdt[ index2(I2A , 0, max_age )] = ir_A[0]*x[index2(C2A,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2A,0,max_age)];
	
	dxdt[ index2(I2R , 0, max_age )] = ir_R[0]*x[index2(C2R,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2R,0,max_age)];
	
	dxdt[ index2(R, 0, max_age )] = alpha_A*x[index2(C2A,0,max_age)] + alpha_R*x[index2(C2R,0,max_age)] + rho*x[index2(I2A,0,max_age)] + rho*x[index2(I2R,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(R,0,max_age)]+ ws*x[index2(RVshort,0,max_age)] + ws*x[index2(RPshort,0,max_age)] + wl*x[index2(RVlong,0,max_age)] + wl*x[index2(RPlong,0,max_age)]; 
	
	// MenAfriVac
	 
	dxdt[ index2(SVlong , 0, max_age )] = phi*x[index2(RAVlong,0,max_age)] + phi*x[index2(RRVlong,0,max_age)] + phi*x[index2(RVlong,0,max_age)] - stoc1*(1-delta)*FOI_A[0]*x[index2(SVlong,0,max_age)] - stoc1*FOI_R[0]*x[index2(SVlong,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(SVlong,0,max_age)] - wl*x[index2(SVlong,0,max_age)];
    
    dxdt[ index2(C1AVlong , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1AVlong,0,max_age)] 
	+ stoc1*(1-delta)*FOI_A[0]*x[index2(SVlong,0,max_age)] - wl*x[index2(C1AVlong,0,max_age)];
	
	dxdt[ index2(C1RVlong , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1RVlong,0,max_age)] + stoc1*FOI_R[0]*x[index2(SVlong,0,max_age)]
	- wl*x[index2(C1RVlong,0,max_age)];
    
    dxdt[ index2(I1AVlong , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1AVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1AVlong,0,max_age)];
	
	dxdt[ index2(I1RVlong , 0, max_age )] = ir_R[0]*x[index2(C1RVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1RVlong,0,max_age)];
    
    dxdt[ index2(RAVlong , 0, max_age )] = rho*x[index2(I1AVlong,0,max_age)]+alpha_A*x[index2(C1AVlong,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RAVlong,0,max_age)] 
	- stoc1*FOI_R[0]*x[index2(RAVlong,0,max_age)] - wl*x[index2(RAVlong,0,max_age)];
    
	dxdt[ index2(RRVlong , 0, max_age )] = rho*x[index2(I1RVlong,0,max_age)]+alpha_R*x[index2(C1RVlong,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RRVlong,0,max_age)]
	 - stoc1*(1-delta)*FOI_A[0]*x[index2(RRVlong,0,max_age)] - wl*x[index2(RRVlong,0,max_age)];
	
	dxdt[ index2(C2AVlong , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2AVlong,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(RRVlong,0,max_age)]
	- wl*x[index2(C2AVlong,0,max_age)];
	
	dxdt[ index2(C2RVlong , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2RVlong,0,max_age)] + stoc1*FOI_R[0]*x[index2(RAVlong,0,max_age)]
	- wl*x[index2(C2RVlong,0,max_age)];
	
	dxdt[ index2(I2AVlong , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C2AVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2AVlong,0,max_age)];
	
	dxdt[ index2(I2RVlong , 0, max_age )] = ir_R[0]*x[index2(C2RVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2RVlong,0,max_age)];
	
	dxdt[ index2(RVlong, 0, max_age )] = alpha_A*x[index2(C2AVlong,0,max_age)] + alpha_R*x[index2(C2RVlong,0,max_age)] + rho*x[index2(I2AVlong,0,max_age)] + rho*x[index2(I2RVlong,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(RVlong,0,max_age)] - wl*x[index2(RVlong,0,max_age)]; 
   
   
   	dxdt[ index2(SVshort , 0, max_age )] = phi*x[index2(RAVshort,0,max_age)] + phi*x[index2(RRVshort,0,max_age)] + phi*x[index2(RVshort,0,max_age)] - stoc1*(1-delta)*FOI_A[0]*x[index2(SVshort,0,max_age)] - stoc1*FOI_R[0]*x[index2(SVshort,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(SVshort,0,max_age)] - ws*x[index2(SVshort,0,max_age)];
    
    dxdt[ index2(C1AVshort , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1AVshort,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(SVshort,0,max_age)]
	- ws*x[index2(C1AVshort,0,max_age)];
	
	dxdt[ index2(C1RVshort , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1RVshort,0,max_age)] + stoc1*FOI_R[0]*x[index2(SVshort,0,max_age)]
	- ws*x[index2(C1RVshort,0,max_age)];
    
    dxdt[ index2(I1AVshort , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1AVshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1AVshort,0,max_age)];
	
	dxdt[ index2(I1RVshort , 0, max_age )] = ir_R[0]*x[index2(C1RVshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1RVshort,0,max_age)];
    
    dxdt[ index2(RAVshort , 0, max_age )] = rho*x[index2(I1AVshort,0,max_age)]+alpha_A*x[index2(C1AVshort,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RAVshort,0,max_age)]
	 - stoc1*FOI_R[0]*x[index2(RAVshort,0,max_age)] - ws*x[index2(RAVshort,0,max_age)];
    
	dxdt[ index2(RRVshort , 0, max_age )] = rho*x[index2(I1RVshort,0,max_age)]+alpha_R*x[index2(C1RVshort,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RRVshort,0,max_age)]
	 - stoc1*(1-delta)*FOI_A[0]*x[index2(RRVshort,0,max_age)] - ws*x[index2(RRVshort,0,max_age)];
	
	dxdt[ index2(C2AVshort , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2AVshort,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(RRVshort,0,max_age)]
	- ws*x[index2(C2AVshort,0,max_age)];
	
	dxdt[ index2(C2RVshort , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2RVshort,0,max_age)] + stoc1*FOI_R[0]*x[index2(RAVshort,0,max_age)]
	- ws*x[index2(C2RVshort,0,max_age)];
	
	dxdt[ index2(I2AVshort , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C2AVshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2AVshort,0,max_age)];
	
	dxdt[ index2(I2RVshort , 0, max_age )] = ir_R[0]*x[index2(C2RVshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2RVshort,0,max_age)];
	
	dxdt[ index2(RVshort, 0, max_age )] = alpha_A*x[index2(C2AVshort,0,max_age)] + alpha_R*x[index2(C2RVshort,0,max_age)] + rho*x[index2(I2AVshort,0,max_age)] + rho*x[index2(I2RVshort,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(RVshort,0,max_age)] - ws*x[index2(RVshort,0,max_age)]; 
	 
	 
	 // Pentavalent
	 
	 
	dxdt[ index2(SPlong , 0, max_age )] = phi*x[index2(RAPlong,0,max_age)] + phi*x[index2(RRPlong,0,max_age)] + phi*x[index2(RPlong,0,max_age)] - stoc1*(1-delta)*FOI_A[0]*x[index2(SPlong,0,max_age)] - stoc1*(1-delta_R)*FOI_R[0]*x[index2(SPlong,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(SPlong,0,max_age)] - wl*x[index2(SPlong,0,max_age)];
    
    dxdt[ index2(C1APlong , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1APlong,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(SPlong,0,max_age)]
	- wl*x[index2(C1APlong,0,max_age)];
	
	dxdt[ index2(C1RPlong , 0, max_age )] = -((1-ksi_R)*ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1RPlong,0,max_age)] + stoc1*(1-delta_R)*FOI_R[0]*x[index2(SPlong,0,max_age)]
	- wl*x[index2(C1RPlong,0,max_age)];
    
    dxdt[ index2(I1APlong , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1APlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1APlong,0,max_age)];
	
	dxdt[ index2(I1RPlong , 0, max_age )] = (1-ksi_R)*ir_R[0]*x[index2(C1RPlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1RPlong,0,max_age)];
    
    dxdt[ index2(RAPlong , 0, max_age )] = rho*x[index2(I1APlong,0,max_age)]+alpha_A*x[index2(C1APlong,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RAPlong,0,max_age)]
	 - stoc1*(1-delta_R)*FOI_R[0]*x[index2(RAPlong,0,max_age)] - wl*x[index2(RAPlong,0,max_age)];
    
	dxdt[ index2(RRPlong , 0, max_age )] = rho*x[index2(I1RPlong,0,max_age)]+alpha_R*x[index2(C1RPlong,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RRPlong,0,max_age)]
	 - stoc1*(1-delta)*FOI_A[0]*x[index2(RRPlong,0,max_age)] - wl*x[index2(RRPlong,0,max_age)];
	
	dxdt[ index2(C2APlong , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2APlong,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(RRPlong,0,max_age)]
	- wl*x[index2(C2APlong,0,max_age)];
	
	dxdt[ index2(C2RPlong , 0, max_age )] = -((1-ksi_R)*ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2RPlong,0,max_age)] + stoc1*(1-delta_R)*FOI_R[0]*x[index2(RAPlong,0,max_age)]
	- wl*x[index2(C2RPlong,0,max_age)];
	
	dxdt[ index2(I2APlong , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C2APlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2APlong,0,max_age)];
	
	dxdt[ index2(I2RPlong , 0, max_age )] = (1-ksi_R)*ir_R[0]*x[index2(C2RPlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2RPlong,0,max_age)];
	
	dxdt[ index2(RPlong, 0, max_age )] = alpha_A*x[index2(C2APlong,0,max_age)] + alpha_R*x[index2(C2RPlong,0,max_age)] + rho*x[index2(I2APlong,0,max_age)] + rho*x[index2(I2RPlong,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(RPlong,0,max_age)] - wl*x[index2(RPlong,0,max_age)]; 
   
   
   	dxdt[ index2(SPshort , 0, max_age )] = phi*x[index2(RAPshort,0,max_age)] + phi*x[index2(RRPshort,0,max_age)] + phi*x[index2(RPshort,0,max_age)] - stoc1*(1-delta)*FOI_A[0]*x[index2(SPshort,0,max_age)] - stoc1*(1-delta_R)*FOI_R[0]*x[index2(SPshort,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(SPshort,0,max_age)] - ws*x[index2(SPshort,0,max_age)];
    
    dxdt[ index2(C1APshort , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1APshort,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(SPshort,0,max_age)]
	- ws*x[index2(C1APshort,0,max_age)];
	
	dxdt[ index2(C1RPshort , 0, max_age )] = -((1-ksi_R)*ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1RPshort,0,max_age)] + stoc1*(1-delta_R)*FOI_R[0]*x[index2(SPshort,0,max_age)]
	- ws*x[index2(C1RPshort,0,max_age)];
    
    dxdt[ index2(I1APshort , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1APshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1APshort,0,max_age)];
	
	dxdt[ index2(I1RPshort , 0, max_age )] = (1-ksi_R)*ir_R[0]*x[index2(C1RPshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1RPshort,0,max_age)];
    
    dxdt[ index2(RAPshort , 0, max_age )] = rho*x[index2(I1APshort,0,max_age)]+alpha_A*x[index2(C1APshort,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RAPshort,0,max_age)]
	 - stoc1*(1-delta_R)*FOI_R[0]*x[index2(RAPshort,0,max_age)] - ws*x[index2(RAPshort,0,max_age)];
    
	dxdt[ index2(RRPshort , 0, max_age )] = rho*x[index2(I1RPshort,0,max_age)]+alpha_R*x[index2(C1RPshort,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RRPshort,0,max_age)]
	 - stoc1*(1-delta)*FOI_A[0]*x[index2(RRPshort,0,max_age)] - ws*x[index2(RRPshort,0,max_age)];
	
	dxdt[ index2(C2APshort , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2APshort,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(RRPshort,0,max_age)]
	- ws*x[index2(C2APshort,0,max_age)];
	
	dxdt[ index2(C2RPshort , 0, max_age )] = -((1-ksi_R)*ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2RPshort,0,max_age)] + stoc1*(1-delta)*FOI_R[0]*x[index2(RAPshort,0,max_age)]
	- ws*x[index2(C2RPshort,0,max_age)];

	dxdt[ index2(I2APshort , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C2APshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2APshort,0,max_age)];
	
	dxdt[ index2(I2RPshort , 0, max_age )] = (1-ksi_R)*ir_R[0]*x[index2(C2RPshort,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2RPshort,0,max_age)];
	
	dxdt[ index2(RPshort, 0, max_age )] = alpha_A*x[index2(C2APshort,0,max_age)] + alpha_R*x[index2(C2RPshort,0,max_age)] + rho*x[index2(I2APshort,0,max_age)] + rho*x[index2(I2RPshort,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(RPshort,0,max_age)] - ws*x[index2(RPshort,0,max_age)]; 
	 
   
   // Incidence
    
    dxdt[ index2(IncA , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)] + ir_A[0]*x[index2(C2A,0,max_age)];
	dxdt[ index2(IncAV , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1AVlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2AVlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C1AVshort,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2AVshort,0,max_age)];
	dxdt[ index2(IncAP , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1APlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2APlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C1APshort,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2APshort,0,max_age)];
	dxdt[ index2(IncA_all , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)] + ir_A[0]*x[index2(C2A,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C1AVlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2AVlong,0,max_age)] 
		+ (1-ksi)*ir_A[0]*x[index2(C1AVshort,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2AVshort,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C1APlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2APlong,0,max_age)] + 
		(1-ksi)*ir_A[0]*x[index2(C1APshort,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2APshort,0,max_age)];
	
	dxdt[ index2(IncR , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)] + ir_R[0]*x[index2(C2R,0,max_age)];
	dxdt[ index2(IncRV , 0, max_age )] = ir_R[0]*x[index2(C1RVlong,0,max_age)] + ir_R[0]*x[index2(C2RVlong,0,max_age)] + ir_R[0]*x[index2(C1RVshort,0,max_age)] + ir_R[0]*x[index2(C2RVshort,0,max_age)];
	dxdt[ index2(IncRP , 0, max_age )] = (1-ksi_R)*ir_R[0]*x[index2(C1RPlong,0,max_age)] + (1-ksi_R)*ir_R[0]*x[index2(C2RPlong,0,max_age)] + (1-ksi_R)*ir_R[0]*x[index2(C1RPshort,0,max_age)] + (1-ksi_R)*ir_R[0]*x[index2(C2RPshort,0,max_age)];
	dxdt[ index2(IncR_all , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)] + ir_R[0]*x[index2(C2R,0,max_age)] + ir_R[0]*x[index2(C1RVlong,0,max_age)] + ir_R[0]*x[index2(C2RVlong,0,max_age)] 
		+ ir_R[0]*x[index2(C1RVshort,0,max_age)] + ir_R[0]*x[index2(C2RVshort,0,max_age)] + (1-ksi_R)*ir_R[0]*x[index2(C1RPlong,0,max_age)] + (1-ksi_R)*ir_R[0]*x[index2(C2RPlong,0,max_age)] + 
		(1-ksi_R)*ir_R[0]*x[index2(C1RPshort,0,max_age)] + (1-ksi_R)*ir_R[0]*x[index2(C2RPshort,0,max_age)];
	
    
	
	//Second age group
	
	    
    dxdt[ index2(S , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(S,0,max_age)]+phi*x[index2(RA,1,max_age)] + phi*x[index2(RR,1,max_age)] + phi*x[index2(R,1,max_age)] - stoc1*FOI_A[1]*x[index2(S,1,max_age)] - stoc1*FOI_R[1]*x[index2(S,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(S,1,max_age)] + ws*x[index2(SVshort,1,max_age)] + wl*x[index2(SVlong,1,max_age)] + ws*x[index2(SPshort,1,max_age)] + wl*x[index2(SPlong,1,max_age)];
    
    dxdt[ index2(C1A , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C1A,0,max_age)]-(ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1A,1,max_age)] + stoc1*FOI_A[1]*x[index2(S,1,max_age)]
	+ ws*x[index2(C1AVshort,1,max_age)] + wl*x[index2(C1AVlong,1,max_age)] + ws*x[index2(C1APshort,1,max_age)] + wl*x[index2(C1APlong,1,max_age)];
	
	dxdt[ index2(C1R , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C1R,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1R,1,max_age)] + stoc1*FOI_R[1]*x[index2(S,1,max_age)]
	+ ws*x[index2(C1RVshort,1,max_age)] + wl*x[index2(C1RVlong,1,max_age)] + ws*x[index2(C1RPshort,1,max_age)] + wl*x[index2(C1RPlong,1,max_age)];
    
    dxdt[ index2(I1A , 1, max_age )] = age*x[index2(I1A,0,max_age)]+ir_A[1]*x[index2(C1A,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1A,1,max_age)];
	
	dxdt[ index2(I1R , 1, max_age )] = age*x[index2(I1R,0,max_age)]+ir_R[1]*x[index2(C1R,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1R,1,max_age)];
    
    dxdt[ index2(RA , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(RA,0,max_age)]+rho*x[index2(I1A,1,max_age)]+alpha_A*x[index2(C1A,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RA,1,max_age)] - stoc1*FOI_R[1]*x[index2(RA,1,max_age)]
	+ ws*x[index2(RAVshort,1,max_age)] + wl*x[index2(RAVlong,1,max_age)] + ws*x[index2(RAPshort,1,max_age)] + wl*x[index2(RAPlong,1,max_age)];
	
	dxdt[ index2(RR , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(RR,0,max_age)]+rho*x[index2(I1R,1,max_age)]+alpha_R*x[index2(C1R,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RR,1,max_age)] - stoc1*FOI_A[1]*x[index2(RR,1,max_age)]
	+ ws*x[index2(RRVshort,1,max_age)] + wl*x[index2(RRVlong,1,max_age)] + ws*x[index2(RRPshort,1,max_age)] + wl*x[index2(RRPlong,1,max_age)];
	
	dxdt[ index2(C2A , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C2A,0,max_age)]-(ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2A,1,max_age)] + stoc1*FOI_A[1]*x[index2(RR,1,max_age)]
	+ ws*x[index2(C2AVshort,1,max_age)] + wl*x[index2(C2AVlong,1,max_age)] + ws*x[index2(C2APshort,1,max_age)] + wl*x[index2(C2APlong,1,max_age)];
	
	dxdt[ index2(C2R , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C2R,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2R,1,max_age)] + stoc1*FOI_R[1]*x[index2(RA,1,max_age)]
	+ ws*x[index2(C2RVshort,1,max_age)] + wl*x[index2(C2RVlong,1,max_age)] + ws*x[index2(C2RPshort,1,max_age)] + wl*x[index2(C2RPlong,1,max_age)];
	
	dxdt[ index2(I2A , 1, max_age )] = age*x[index2(I2A,0,max_age)]+ir_A[1]*x[index2(C2A,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2A,1,max_age)];
	
	dxdt[ index2(I2R , 1, max_age )] = age*x[index2(I2R,0,max_age)]+ir_R[1]*x[index2(C2R,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2R,1,max_age)];
	
	dxdt[ index2(R, 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(R,0,max_age)]+alpha_A*x[index2(C2A,1,max_age)] + alpha_R*x[index2(C2R,1,max_age)] + rho*x[index2(I2A,1,max_age)] + rho*x[index2(I2R,1,max_age)]
	 - (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(R,1,max_age)] + ws*x[index2(RVshort,1,max_age)] + wl*x[index2(RVlong,1,max_age)] + ws*x[index2(RPshort,1,max_age)] + wl*x[index2(RPlong,1,max_age)]; 
	
	// MenAfriVac
	 
	dxdt[ index2(SVlong , 1, max_age )] = age*x[index2(SVlong,0,max_age)]+phi*x[index2(RAVlong,1,max_age)] + phi*x[index2(RRVlong,1,max_age)] + phi*x[index2(RVlong,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(SVlong,1,max_age)] - stoc1*FOI_R[1]*x[index2(SVlong,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(SVlong,1,max_age)] - wl*x[index2(SVlong,1,max_age)];
    
    dxdt[ index2(C1AVlong , 1, max_age )] = age*x[index2(C1AVlong,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1AVlong,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(SVlong,1,max_age)]
	-wl*x[index2(C1AVlong,1,max_age)];
	
	dxdt[ index2(C1RVlong , 1, max_age )] = age*x[index2(C1RVlong,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1RVlong,1,max_age)] + stoc1*FOI_R[1]*x[index2(SVlong,1,max_age)]
	-wl*x[index2(C1RVlong,1,max_age)];
    
    dxdt[ index2(I1AVlong , 1, max_age )] = age*x[index2(I1AVlong,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C1AVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1AVlong,1,max_age)];
	
	dxdt[ index2(I1RVlong , 1, max_age )] = age*x[index2(I1RVlong,0,max_age)]+ir_R[1]*x[index2(C1RVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1RVlong,1,max_age)];
    
    dxdt[ index2(RAVlong , 1, max_age )] = age*x[index2(RAVlong,0,max_age)]+rho*x[index2(I1AVlong,1,max_age)]+alpha_A*x[index2(C1AVlong,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RAVlong,1,max_age)] - stoc1*FOI_R[1]*x[index2(RAVlong,1,max_age)]
	-wl*x[index2(RAVlong,1,max_age)];
    
	dxdt[ index2(RRVlong , 1, max_age )] = age*x[index2(RRVlong,0,max_age)]+rho*x[index2(I1RVlong,1,max_age)]+alpha_R*x[index2(C1RVlong,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RRVlong,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(RRVlong,1,max_age)]
	-wl*x[index2(RRVlong,1,max_age)];
	
	dxdt[ index2(C2AVlong , 1, max_age )] = age*x[index2(C2AVlong,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2AVlong,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(RRVlong,1,max_age)]
	-wl*x[index2(C2AVlong,1,max_age)];
	
	dxdt[ index2(C2RVlong , 1, max_age )] = age*x[index2(C2RVlong,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2RVlong,1,max_age)] + stoc1*FOI_R[1]*x[index2(RAVlong,1,max_age)]
	-wl*x[index2(C2RVlong,1,max_age)];
	
	dxdt[ index2(I2AVlong , 1, max_age )] = age*x[index2(I2AVlong,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C2AVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2AVlong,1,max_age)];
	
	dxdt[ index2(I2RVlong , 1, max_age )] = age*x[index2(I2RVlong,0,max_age)]+ir_R[1]*x[index2(C2RVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2RVlong,1,max_age)];
	
	dxdt[ index2(RVlong, 1, max_age )] = age*x[index2(RVlong,0,max_age)]+alpha_A*x[index2(C2AVlong,1,max_age)] + alpha_R*x[index2(C2RVlong,1,max_age)] + rho*x[index2(I2AVlong,1,max_age)] + rho*x[index2(I2RVlong,1,max_age)]
	 	- (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(RVlong,1,max_age)] - wl*x[index2(RVlong,1,max_age)]; 
   
   
   	dxdt[ index2(SVshort , 1, max_age )] = age*x[index2(SVshort,0,max_age)]+phi*x[index2(RAVshort,1,max_age)] + phi*x[index2(RRVshort,1,max_age)] + phi*x[index2(RVshort,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(SVshort,1,max_age)] - stoc1*FOI_R[1]*x[index2(SVshort,1,max_age)]
	 	- (age+mu[index2(floor(t/365),1,max_age)])*x[index2(SVshort,1,max_age)] + epi_A[t/365]*age*x[index2(S,0,max_age)] - ws*x[index2(SVshort,1,max_age)];
    
    dxdt[ index2(C1AVshort , 1, max_age )] = age*x[index2(C1AVshort,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1AVshort,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(SVshort,1,max_age)]
		+ epi_A[t/365]*age*x[index2(C1A,0,max_age)] - ws*x[index2(C1AVshort,1,max_age)];
	
	dxdt[ index2(C1RVshort , 1, max_age )] = age*x[index2(C1RVshort,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1RVshort,1,max_age)] + stoc1*FOI_R[1]*x[index2(SVshort,1,max_age)]
		+ epi_A[t/365]*age*x[index2(C1R,0,max_age)] - ws*x[index2(C1RVshort,1,max_age)];
    
    dxdt[ index2(I1AVshort , 1, max_age )] = age*x[index2(I1AVshort,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C1AVshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1AVshort,1,max_age)];
	
	dxdt[ index2(I1RVshort , 1, max_age )] = age*x[index2(I1RVshort,0,max_age)]+ir_R[1]*x[index2(C1RVshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1RVshort,1,max_age)];
    
    dxdt[ index2(RAVshort , 1, max_age )] = age*x[index2(RAVshort,0,max_age)]+rho*x[index2(I1AVshort,1,max_age)]+alpha_A*x[index2(C1AVshort,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RAVshort,1,max_age)] - stoc1*FOI_R[1]*x[index2(RAVshort,1,max_age)]
		+ epi_A[t/365]*age*x[index2(RA,0,max_age)] - ws*x[index2(RAVshort,1,max_age)];
    
	dxdt[ index2(RRVshort , 1, max_age )] = age*x[index2(RRVshort,0,max_age)]+rho*x[index2(I1RVshort,1,max_age)]+alpha_R*x[index2(C1RVshort,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RRVshort,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(RRVshort,1,max_age)]
		+ epi_A[t/365]*age*x[index2(RR,0,max_age)] - ws*x[index2(RRVshort,1,max_age)];
	
	dxdt[ index2(C2AVshort , 1, max_age )] = age*x[index2(C2AVshort,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2AVshort,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(RRVshort,1,max_age)]
		+ epi_A[t/365]*age*x[index2(C2A,0,max_age)] - ws*x[index2(C2AVshort,1,max_age)];
	
	dxdt[ index2(C2RVshort , 1, max_age )] = age*x[index2(C2RVshort,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2RVshort,1,max_age)] + stoc1*FOI_R[1]*x[index2(RAVshort,1,max_age)]
		+ epi_A[t/365]*age*x[index2(C2R,0,max_age)] - ws*x[index2(C2RVshort,1,max_age)];
	
	dxdt[ index2(I2AVshort , 1, max_age )] = age*x[index2(I2AVshort,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C2AVshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2AVshort,1,max_age)];
	
	dxdt[ index2(I2RVshort , 1, max_age )] = age*x[index2(I2RVshort,0,max_age)]+ir_R[1]*x[index2(C2RVshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2RVshort,1,max_age)];
	
	dxdt[ index2(RVshort, 1, max_age )] = age*x[index2(RVshort,0,max_age)]+alpha_A*x[index2(C2AVshort,1,max_age)] + alpha_R*x[index2(C2RVshort,1,max_age)] + rho*x[index2(I2AVshort,1,max_age)] + rho*x[index2(I2RVshort,1,max_age)]
		- (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(RVshort,1,max_age)] + epi_A[t/365]*age*x[index2(R,0,max_age)] - ws*x[index2(RVshort,1,max_age)]; 
	 
	 
	 // Pentavalent
	 
	 
	dxdt[ index2(SPlong , 1, max_age )] = age*x[index2(SPlong,0,max_age)]+phi*x[index2(RAPlong,1,max_age)] + phi*x[index2(RRPlong,1,max_age)] + phi*x[index2(RPlong,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(SPlong,1,max_age)] - stoc1*(1-delta_R)*FOI_R[1]*x[index2(SPlong,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(SPlong,1,max_age)] - wl*x[index2(SPlong,1,max_age)];
    
    dxdt[ index2(C1APlong , 1, max_age )] = age*x[index2(C1APlong,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1APlong,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(SPlong,1,max_age)]
	- wl*x[index2(C1APlong,1,max_age)];
	
	dxdt[ index2(C1RPlong , 1, max_age )] = age*x[index2(C1RPlong,0,max_age)]-((1-ksi_R)*ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1RPlong,1,max_age)] + stoc1*(1-delta_R)*FOI_R[1]*x[index2(SPlong,1,max_age)]
	- wl*x[index2(C1RPlong,1,max_age)];
    
    dxdt[ index2(I1APlong , 1, max_age )] = age*x[index2(I1APlong,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C1APlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1APlong,1,max_age)];
	
	dxdt[ index2(I1RPlong , 1, max_age )] = age*x[index2(I1RPlong,0,max_age)]+(1-ksi_R)*ir_R[1]*x[index2(C1RPlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1RPlong,1,max_age)];
    
    dxdt[ index2(RAPlong , 1, max_age )] = age*x[index2(RAPlong,0,max_age)]+rho*x[index2(I1APlong,1,max_age)]+alpha_A*x[index2(C1APlong,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RAPlong,1,max_age)] - stoc1*(1-delta)*FOI_R[1]*x[index2(RAPlong,1,max_age)]
	- wl*x[index2(RAPlong,1,max_age)];
    
	dxdt[ index2(RRPlong , 1, max_age )] = age*x[index2(RRPlong,0,max_age)]+rho*x[index2(I1RPlong,1,max_age)]+alpha_R*x[index2(C1RPlong,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RRPlong,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(RRPlong,1,max_age)]
	- wl*x[index2(RRPlong,1,max_age)];
	
	dxdt[ index2(C2APlong , 1, max_age )] = age*x[index2(C2APlong,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2APlong,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(RRPlong,1,max_age)]
	- wl*x[index2(C2APlong,1,max_age)];
	
	dxdt[ index2(C2RPlong , 1, max_age )] = age*x[index2(C2RPlong,0,max_age)]-((1-ksi_R)*ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2RPlong,1,max_age)] + stoc1*(1-delta_R)*FOI_R[1]*x[index2(RAPlong,1,max_age)]
	- wl*x[index2(C2RPlong,1,max_age)];
	
	dxdt[ index2(I2APlong , 1, max_age )] = age*x[index2(I2APlong,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C2APlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2APlong,1,max_age)];
	
	dxdt[ index2(I2RPlong , 1, max_age )] = age*x[index2(I2RPlong,0,max_age)]+(1-ksi_R)*ir_R[1]*x[index2(C2RPlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2RPlong,1,max_age)];
	
	dxdt[ index2(RPlong, 1, max_age )] = age*x[index2(RPlong,0,max_age)]+alpha_A*x[index2(C2APlong,1,max_age)] + alpha_R*x[index2(C2RPlong,1,max_age)] + rho*x[index2(I2APlong,1,max_age)] + rho*x[index2(I2RPlong,1,max_age)]
	 - (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(RPlong,1,max_age)] - wl*x[index2(RPlong,1,max_age)]; 
   
   
   	dxdt[ index2(SPshort , 1, max_age )] = age*x[index2(SPshort,0,max_age)]+phi*x[index2(RAPshort,1,max_age)] + phi*x[index2(RRPshort,1,max_age)] + phi*x[index2(RPshort,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(SPshort,1,max_age)] - stoc1*(1-delta_R)*FOI_R[1]*x[index2(SPshort,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(SPshort,1,max_age)] + epi_R[t/365]*age*x[index2(S,0,max_age)] - ws*x[index2(SPshort,1,max_age)];
    
    dxdt[ index2(C1APshort , 1, max_age )] = age*x[index2(C1APshort,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1APshort,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(SPshort,1,max_age)]
		+ epi_R[t/365]*age*x[index2(C1A,0,max_age)] - ws*x[index2(C1APshort,1,max_age)];
	
	dxdt[ index2(C1RPshort , 1, max_age )] = age*x[index2(C1RPshort,0,max_age)]-((1-ksi_R)*ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1RPshort,1,max_age)] + stoc1*(1-delta_R)*FOI_R[1]*x[index2(SPshort,1,max_age)]
		+ epi_R[t/365]*age*x[index2(C1R,0,max_age)] - ws*x[index2(C1RPshort,1,max_age)];
    
    dxdt[ index2(I1APshort , 1, max_age )] = age*x[index2(I1APshort,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C1APshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1APshort,1,max_age)];
	
	dxdt[ index2(I1RPshort , 1, max_age )] = age*x[index2(I1RPshort,0,max_age)]+(1-ksi_R)*ir_R[1]*x[index2(C1RPshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1RPshort,1,max_age)];
    
    dxdt[ index2(RAPshort , 1, max_age )] = age*x[index2(RAPshort,0,max_age)]+rho*x[index2(I1APshort,1,max_age)]+alpha_A*x[index2(C1APshort,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RAPshort,1,max_age)] - stoc1*(1-delta)*FOI_R[1]*x[index2(RAPshort,1,max_age)]
		+ epi_R[t/365]*age*x[index2(RA,0,max_age)] - ws*x[index2(RAPshort,1,max_age)];
    
	dxdt[ index2(RRPshort , 1, max_age )] = age*x[index2(RRPshort,0,max_age)]+rho*x[index2(I1RPshort,1,max_age)]+alpha_R*x[index2(C1RPshort,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RRPshort,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(RRPshort,1,max_age)]
		+ epi_R[t/365]*age*x[index2(RR,0,max_age)] - ws*x[index2(RRPshort,1,max_age)];
	
	dxdt[ index2(C2APshort , 1, max_age )] = age*x[index2(C2APshort,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2APshort,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(RRPshort,1,max_age)]
		+ epi_R[t/365]*age*x[index2(C2A,0,max_age)] - ws*x[index2(C2APshort,1,max_age)];
	
	dxdt[ index2(C2RPshort , 1, max_age )] = age*x[index2(C2RPshort,0,max_age)]-((1-ksi_R)*ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2RPshort,1,max_age)] + stoc1*(1-delta_R)*FOI_R[1]*x[index2(RAPshort,1,max_age)]
		+ epi_R[t/365]*age*x[index2(C2R,0,max_age)] - ws*x[index2(C2RPshort,1,max_age)];
	
	dxdt[ index2(I2APshort , 1, max_age )] = age*x[index2(I2APshort,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C2APshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2APshort,1,max_age)];
	
	dxdt[ index2(I2RPshort , 1, max_age )] = age*x[index2(I2RPshort,0,max_age)]+(1-ksi_R)*ir_R[1]*x[index2(C2RPshort,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2RPshort,1,max_age)];
	
	dxdt[ index2(RPshort, 1, max_age )] = age*x[index2(RPshort,0,max_age)]+alpha_A*x[index2(C2APshort,1,max_age)] + alpha_R*x[index2(C2RPshort,1,max_age)] + rho*x[index2(I2APshort,1,max_age)] + rho*x[index2(I2RPshort,1,max_age)]
	 - (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(RPshort,1,max_age)] + epi_R[t/365]*age*x[index2(R,0,max_age)] - ws*x[index2(RPshort,1,max_age)]; 
	 
	 // Incidence
	   
	 dxdt[ index2(IncA , 1, max_age )] = ir_A[1]*x[index2(C1A,1,max_age)] + ir_A[1]*x[index2(C2A,1,max_age)];
	dxdt[ index2(IncAV , 1, max_age )] = (1-ksi)*ir_A[1]*x[index2(C1AVlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2AVlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C1AVshort,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2AVshort,1,max_age)];
	dxdt[ index2(IncAP , 1, max_age )] = (1-ksi)*ir_A[1]*x[index2(C1APlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2APlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C1APshort,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2APshort,1,max_age)];
	dxdt[ index2(IncA_all , 1, max_age )] = ir_A[1]*x[index2(C1A,1,max_age)] + ir_A[1]*x[index2(C2A,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C1AVlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2AVlong,1,max_age)] 
		+ (1-ksi)*ir_A[1]*x[index2(C1AVshort,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2AVshort,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C1APlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2APlong,1,max_age)] + 
		(1-ksi)*ir_A[1]*x[index2(C1APshort,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2APshort,1,max_age)];
	
	dxdt[ index2(IncR , 1, max_age )] = ir_R[1]*x[index2(C1R,1,max_age)] + ir_R[1]*x[index2(C2R,1,max_age)];
	dxdt[ index2(IncRV , 1, max_age )] = ir_R[1]*x[index2(C1RVlong,1,max_age)] + ir_R[1]*x[index2(C2RVlong,1,max_age)] + ir_R[1]*x[index2(C1RVshort,1,max_age)] + ir_R[1]*x[index2(C2RVshort,1,max_age)];
	dxdt[ index2(IncRP , 1, max_age )] = (1-ksi_R)*ir_R[1]*x[index2(C1RPlong,1,max_age)] + (1-ksi_R)*ir_R[1]*x[index2(C2RPlong,1,max_age)] + (1-ksi_R)*ir_R[1]*x[index2(C1RPshort,1,max_age)] + (1-ksi_R)*ir_R[1]*x[index2(C2RPshort,1,max_age)];
	dxdt[ index2(IncR_all , 1, max_age )] = ir_R[1]*x[index2(C1R,1,max_age)] + ir_R[1]*x[index2(C2R,1,max_age)] + ir_R[1]*x[index2(C1RVlong,1,max_age)] + ir_R[1]*x[index2(C2RVlong,1,max_age)] 
		+ ir_R[1]*x[index2(C1RVshort,1,max_age)] + ir_R[1]*x[index2(C2RVshort,1,max_age)] + (1-ksi_R)*ir_R[1]*x[index2(C1RPlong,1,max_age)] + (1-ksi_R)*ir_R[1]*x[index2(C2RPlong,1,max_age)] + 
		(1-ksi_R)*ir_R[1]*x[index2(C1RPshort,1,max_age)] + (1-ksi_R)*ir_R[1]*x[index2(C2RPshort,1,max_age)];
	
	 
    
    // Annual age groups
    for( size_t a=2 ; a<(max_age-1) ; ++a )
    {
        
            dxdt[ index2(S , a, max_age )] = age*(x[index2(S,a-1,max_age)]-x[index2(S,a,max_age)]) + phi*x[index2(RA,a,max_age)] + phi*x[index2(RR,a,max_age)] + phi*x[index2(R,a,max_age)] - stoc1*FOI_A[a]*x[index2(S,a,max_age)]
			 - stoc1*FOI_R[a]*x[index2(S,a,max_age)] - (mu[index2(floor(t/365),a,max_age)])*x[index2(S,a,max_age)] + ws*x[index2(SVshort,a,max_age)] + wl*x[index2(SVlong,a,max_age)] + ws*x[index2(SPshort,a,max_age)] + wl*x[index2(SPlong,a,max_age)];
    
    dxdt[ index2(C1A , a, max_age )] = age*(x[index2(C1A,a-1,max_age)]-x[index2(C1A,a,max_age)]) -(ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)])*x[index2(C1A,a,max_age)] + stoc1*FOI_A[a]*x[index2(S,a,max_age)]
	+ ws*x[index2(C1AVshort,a,max_age)] + wl*x[index2(C1AVlong,a,max_age)] + ws*x[index2(C1APshort,a,max_age)] + wl*x[index2(C1APlong,a,max_age)];
	
	dxdt[ index2(C1R , a, max_age )] = age*(x[index2(C1R,a-1,max_age)]-x[index2(C1R,a,max_age)]) -(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)])*x[index2(C1R,a,max_age)] + stoc1*FOI_R[a]*x[index2(S,a,max_age)]
	+ ws*x[index2(C1RVshort,a,max_age)] + wl*x[index2(C1RVlong,a,max_age)] + ws*x[index2(C1RPshort,a,max_age)] + wl*x[index2(C1RPlong,a,max_age)];
    
    dxdt[ index2(I1A , a, max_age )] = age*(x[index2(I1A,a-1,max_age)]-x[index2(I1A,a,max_age)]) +ir_A[a]*x[index2(C1A,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I1A,a,max_age)];
	
	dxdt[ index2(I1R , a, max_age )] = age*(x[index2(I1R,a-1,max_age)]-x[index2(I1R,a,max_age)]) +ir_R[a]*x[index2(C1R,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I1R,a,max_age)];
    
    dxdt[ index2(RA , a, max_age )] = age*(x[index2(RA,a-1,max_age)]-x[index2(RA,a,max_age)]) +rho*x[index2(I1A,a,max_age)]+alpha_A*x[index2(C1A,a,max_age)]-(phi+mu[index2(floor(t/365),a,max_age)])*x[index2(RA,a,max_age)]
	 - stoc1*FOI_R[a]*x[index2(RA,a,max_age)] + ws*x[index2(RAVshort,a,max_age)] + wl*x[index2(RAVlong,a,max_age)] + ws*x[index2(RAPshort,a,max_age)] + wl*x[index2(RAPlong,a,max_age)];
    
	dxdt[ index2(RR , a, max_age )] = age*(x[index2(RR,a-1,max_age)]-x[index2(RR,a,max_age)]) +rho*x[index2(I1R,a,max_age)]+alpha_R*x[index2(C1R,a,max_age)]-(phi+mu[index2(floor(t/365),a,max_age)])*x[index2(RR,a,max_age)]
	 - stoc1*FOI_A[a]*x[index2(RR,a,max_age)] + ws*x[index2(RRVshort,a,max_age)] + wl*x[index2(RRVlong,a,max_age)] + ws*x[index2(RRPshort,a,max_age)] + wl*x[index2(RRPlong,a,max_age)];
	
	dxdt[ index2(C2A , a, max_age )] = age*(x[index2(C2A,a-1,max_age)]-x[index2(C2A,a,max_age)]) -(ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)])*x[index2(C2A,a,max_age)] + stoc1*FOI_A[a]*x[index2(RR,a,max_age)]
	+ ws*x[index2(C2AVshort,a,max_age)] + wl*x[index2(C2AVlong,a,max_age)] + ws*x[index2(C2APshort,a,max_age)] + wl*x[index2(C2APlong,a,max_age)];
	
	dxdt[ index2(C2R , a, max_age )] = age*(x[index2(C2R,a-1,max_age)]-x[index2(C2R,a,max_age)]) -(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)])*x[index2(C2R,a,max_age)] + stoc1*FOI_R[a]*x[index2(RA,a,max_age)]
	+ ws*x[index2(C2RVshort,a,max_age)] + wl*x[index2(C2RVlong,a,max_age)] + ws*x[index2(C2RPshort,a,max_age)] + wl*x[index2(C2RPlong,a,max_age)];
	
	dxdt[ index2(I2A , a, max_age )] = age*(x[index2(I2A,a-1,max_age)]-x[index2(I2A,a,max_age)]) +ir_A[a]*x[index2(C2A,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I2A,a,max_age)];
	
	dxdt[ index2(I2R , a, max_age )] = age*(x[index2(I2R,a-1,max_age)]-x[index2(I2R,a,max_age)]) +ir_R[a]*x[index2(C2R,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I2R,a,max_age)];
	
	dxdt[ index2(R, a, max_age )] = age*(x[index2(R,a-1,max_age)]-x[index2(R,a,max_age)]) +alpha_A*x[index2(C2A,a,max_age)] + alpha_R*x[index2(C2R,a,max_age)] + rho*x[index2(I2A,a,max_age)] + rho*x[index2(I2R,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)])*x[index2(R,a,max_age)] + ws*x[index2(RVshort,a,max_age)] + wl*x[index2(RVlong,a,max_age)] + ws*x[index2(RPshort,a,max_age)] + wl*x[index2(RPlong,a,max_age)]; 
	 
	 
	 	// MenAfriVac
	 
	dxdt[ index2(SVlong , a, max_age )] = age*x[index2(SVlong,a-1,max_age)]+phi*x[index2(RAVlong,a,max_age)] + phi*x[index2(RRVlong,a,max_age)] + phi*x[index2(RVlong,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(SVlong,a,max_age)] - stoc1*FOI_R[a]*x[index2(SVlong,a,max_age)]
	 - (age+mu[index2(floor(t/365),a,max_age)])*x[index2(SVlong,a,max_age)] - wl*x[index2(SVlong,a,max_age)];
    
    dxdt[ index2(C1AVlong , a, max_age )] = age*x[index2(C1AVlong,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1AVlong,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(SVlong,a,max_age)]
	- wl*x[index2(C1AVlong,a,max_age)];
	
	dxdt[ index2(C1RVlong , a, max_age )] = age*x[index2(C1RVlong,a-1,max_age)]-(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1RVlong,a,max_age)] + stoc1*FOI_R[a]*x[index2(SVlong,a,max_age)]
	- wl*x[index2(C1RVlong,a,max_age)];
    
    dxdt[ index2(I1AVlong , a, max_age )] = age*x[index2(I1AVlong,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C1AVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1AVlong,a,max_age)];
	
	dxdt[ index2(I1RVlong , a, max_age )] = age*x[index2(I1RVlong,a-1,max_age)]+ir_R[a]*x[index2(C1RVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1RVlong,a,max_age)];
    
    dxdt[ index2(RAVlong , a, max_age )] = age*x[index2(RAVlong,a-1,max_age)]+rho*x[index2(I1AVlong,a,max_age)]+alpha_A*x[index2(C1AVlong,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RAVlong,a,max_age)] - stoc1*FOI_R[a]*x[index2(RAVlong,a,max_age)]
	- wl*x[index2(RAVlong,a,max_age)];
    
	dxdt[ index2(RRVlong , a, max_age )] = age*x[index2(RRVlong,a-1,max_age)]+rho*x[index2(I1RVlong,a,max_age)]+alpha_R*x[index2(C1RVlong,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RRVlong,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(RRVlong,a,max_age)]
	- wl*x[index2(RRVlong,a,max_age)];
	
	dxdt[ index2(C2AVlong , a, max_age )] = age*x[index2(C2AVlong,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2AVlong,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(RRVlong,a,max_age)]
	- wl*x[index2(C2AVlong,a,max_age)];
	
	dxdt[ index2(C2RVlong , a, max_age )] = age*x[index2(C2RVlong,a-1,max_age)]-(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2RVlong,a,max_age)] + stoc1*FOI_R[a]*x[index2(RAVlong,a,max_age)]
	- wl*x[index2(C2RVlong,a,max_age)];
	
	dxdt[ index2(I2AVlong , a, max_age )] = age*x[index2(I2AVlong,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C2AVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2AVlong,a,max_age)];
	
	dxdt[ index2(I2RVlong , a, max_age )] = age*x[index2(I2RVlong,a-1,max_age)]+ir_R[a]*x[index2(C2RVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2RVlong,a,max_age)];
	
	dxdt[ index2(RVlong, a, max_age )] = age*x[index2(RVlong,a-1,max_age)]+alpha_A*x[index2(C2AVlong,a,max_age)] + alpha_R*x[index2(C2RVlong,a,max_age)] + rho*x[index2(I2AVlong,a,max_age)] + rho*x[index2(I2RVlong,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(RVlong,a,max_age)] - wl*x[index2(RVlong,a,max_age)]; 
   
   
   	dxdt[ index2(SVshort , a, max_age )] = age*x[index2(SVshort,a-1,max_age)]+phi*x[index2(RAVshort,a,max_age)] + phi*x[index2(RRVshort,a,max_age)] + phi*x[index2(RVshort,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(SVshort,a,max_age)] - stoc1*FOI_R[a]*x[index2(SVshort,a,max_age)]
	 - (age+mu[index2(floor(t/365),a,max_age)])*x[index2(SVshort,a,max_age)] - ws*x[index2(SVshort,a,max_age)];
    
    dxdt[ index2(C1AVshort , a, max_age )] = age*x[index2(C1AVshort,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1AVshort,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(SVshort,a,max_age)]
	- ws*x[index2(C1AVshort,a,max_age)];
	
	dxdt[ index2(C1RVshort , a, max_age )] = age*x[index2(C1RVshort,a-1,max_age)]-(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1RVshort,a,max_age)] + stoc1*FOI_R[a]*x[index2(SVshort,a,max_age)]
	- ws*x[index2(C1RVshort,a,max_age)];
    
    dxdt[ index2(I1AVshort , a, max_age )] = age*x[index2(I1AVshort,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C1AVshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1AVshort,a,max_age)];
	
	dxdt[ index2(I1RVshort , a, max_age )] = age*x[index2(I1RVshort,a-1,max_age)]+ir_R[a]*x[index2(C1RVshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1RVshort,a,max_age)];
    
    dxdt[ index2(RAVshort , a, max_age )] = age*x[index2(RAVshort,a-1,max_age)]+rho*x[index2(I1AVshort,a,max_age)]+alpha_A*x[index2(C1AVshort,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RAVshort,a,max_age)] - stoc1*FOI_R[a]*x[index2(RAVshort,a,max_age)]
	- ws*x[index2(RAVshort,a,max_age)];
    
	dxdt[ index2(RRVshort , a, max_age )] = age*x[index2(RRVshort,a-1,max_age)]+rho*x[index2(I1RVshort,a,max_age)]+alpha_R*x[index2(C1RVshort,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RRVshort,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(RRVshort,a,max_age)]
	- ws*x[index2(RRVshort,a,max_age)];
	
	dxdt[ index2(C2AVshort , a, max_age )] = age*x[index2(C2AVshort,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2AVshort,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(RRVshort,a,max_age)]
	- ws*x[index2(C2AVshort,a,max_age)];
	
	dxdt[ index2(C2RVshort , a, max_age )] = age*x[index2(C2RVshort,a-1,max_age)]-(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2RVshort,a,max_age)] + stoc1*FOI_R[a]*x[index2(RAVshort,a,max_age)]
	- ws*x[index2(C2RVshort,a,max_age)];
	
	dxdt[ index2(I2AVshort , a, max_age )] = age*x[index2(I2AVshort,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C2AVshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2AVshort,a,max_age)];
	
	dxdt[ index2(I2RVshort , a, max_age )] = age*x[index2(I2RVshort,a-1,max_age)]+ir_R[a]*x[index2(C2RVshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2RVshort,a,max_age)];
	
	dxdt[ index2(RVshort, a, max_age )] = age*x[index2(RVshort,a-1,max_age)]+alpha_A*x[index2(C2AVshort,a,max_age)] + alpha_R*x[index2(C2RVshort,a,max_age)] + rho*x[index2(I2AVshort,a,max_age)] + rho*x[index2(I2RVshort,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(RVshort,a,max_age)] - ws*x[index2(RVshort,a,max_age)]; 
	 
	 
	 // Pentavalent
	 
	 
	dxdt[ index2(SPlong , a, max_age )] = age*x[index2(SPlong,a-1,max_age)]+phi*x[index2(RAPlong,a,max_age)] + phi*x[index2(RRPlong,a,max_age)] + phi*x[index2(RPlong,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(SPlong,a,max_age)] - stoc1*(1-delta_R)*FOI_R[a]*x[index2(SPlong,a,max_age)]
	 - (age+mu[index2(floor(t/365),a,max_age)])*x[index2(SPlong,a,max_age)] - wl*x[index2(SPlong,a,max_age)];
    
    dxdt[ index2(C1APlong , a, max_age )] = age*x[index2(C1APlong,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1APlong,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(SPlong,a,max_age)]
	- wl*x[index2(C1APlong,a,max_age)];
	
	dxdt[ index2(C1RPlong , a, max_age )] = age*x[index2(C1RPlong,a-1,max_age)]-((1-ksi_R)*ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1RPlong,a,max_age)] + stoc1*(1-delta_R)*FOI_R[a]*x[index2(SPlong,a,max_age)]
	- wl*x[index2(C1RPlong,a,max_age)];
    
    dxdt[ index2(I1APlong , a, max_age )] = age*x[index2(I1APlong,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C1APlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1APlong,a,max_age)];
	
	dxdt[ index2(I1RPlong , a, max_age )] = age*x[index2(I1RPlong,a-1,max_age)]+(1-ksi_R)*ir_R[a]*x[index2(C1RPlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1RPlong,a,max_age)];
    
    dxdt[ index2(RAPlong , a, max_age )] = age*x[index2(RAPlong,a-1,max_age)]+rho*x[index2(I1APlong,a,max_age)]+alpha_A*x[index2(C1APlong,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RAPlong,a,max_age)] - stoc1*(1-delta)*FOI_R[a]*x[index2(RAPlong,a,max_age)]
	- wl*x[index2(RAPlong,a,max_age)];
    
	dxdt[ index2(RRPlong , a, max_age )] = age*x[index2(RRPlong,a-1,max_age)]+rho*x[index2(I1RPlong,a,max_age)]+alpha_R*x[index2(C1RPlong,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RRPlong,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(RRPlong,a,max_age)]
	- wl*x[index2(RRPlong,a,max_age)];
	
	dxdt[ index2(C2APlong , a, max_age )] = age*x[index2(C2APlong,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2APlong,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(RRPlong,a,max_age)]
	- wl*x[index2(C2APlong,a,max_age)];
	
	dxdt[ index2(C2RPlong , a, max_age )] = age*x[index2(C2RPlong,a-1,max_age)]-((1-ksi_R)*ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2RPlong,a,max_age)] + stoc1*(1-delta_R)*FOI_R[a]*x[index2(RAPlong,a,max_age)]
	- wl*x[index2(C2RPlong,a,max_age)];
	
	dxdt[ index2(I2APlong , a, max_age )] = age*x[index2(I2APlong,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C2APlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2APlong,a,max_age)];
	
	dxdt[ index2(I2RPlong , a, max_age )] = age*x[index2(I2RPlong,a-1,max_age)]+(1-ksi_R)*ir_R[a]*x[index2(C2RPlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2RPlong,a,max_age)];
	
	dxdt[ index2(RPlong, a, max_age )] = age*x[index2(RPlong,a-1,max_age)]+alpha_A*x[index2(C2APlong,a,max_age)] + alpha_R*x[index2(C2RPlong,a,max_age)] + rho*x[index2(I2APlong,a,max_age)] + rho*x[index2(I2RPlong,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(RPlong,a,max_age)] - wl*x[index2(RPlong,a,max_age)]; 
   
   
   	dxdt[ index2(SPshort , a, max_age )] = age*x[index2(SPshort,a-1,max_age)]+phi*x[index2(RAPshort,a,max_age)] + phi*x[index2(RRPshort,a,max_age)] + phi*x[index2(RPshort,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(SPshort,a,max_age)] - stoc1*(1-delta_R)*FOI_R[a]*x[index2(SPshort,a,max_age)]
	 - (age+mu[index2(floor(t/365),a,max_age)])*x[index2(SPshort,a,max_age)] - ws*x[index2(SPshort,a,max_age)];
    
    dxdt[ index2(C1APshort , a, max_age )] = age*x[index2(C1APshort,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1APshort,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(SPshort,a,max_age)]
	- ws*x[index2(C1APshort,a,max_age)];
	
	dxdt[ index2(C1RPshort , a, max_age )] = age*x[index2(C1RPshort,a-1,max_age)]-((1-ksi_R)*ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1RPshort,a,max_age)] + stoc1*(1-delta_R)*FOI_R[a]*x[index2(SPshort,a,max_age)]
	- ws*x[index2(C1RPshort,a,max_age)];
    
    dxdt[ index2(I1APshort , a, max_age )] = age*x[index2(I1APshort,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C1APshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1APshort,a,max_age)];
	
	dxdt[ index2(I1RPshort , a, max_age )] = age*x[index2(I1RPshort,a-1,max_age)]+(1-ksi_R)*ir_R[a]*x[index2(C1RPshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1RPshort,a,max_age)];
    
    dxdt[ index2(RAPshort , a, max_age )] = age*x[index2(RAPshort,a-1,max_age)]+rho*x[index2(I1APshort,a,max_age)]+alpha_A*x[index2(C1APshort,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RAPshort,a,max_age)] - stoc1*(1-delta)*FOI_R[a]*x[index2(RAPshort,a,max_age)]
	- ws*x[index2(RAPshort,a,max_age)];
    
	dxdt[ index2(RRPshort , a, max_age )] = age*x[index2(RRPshort,a-1,max_age)]+rho*x[index2(I1RPshort,a,max_age)]+alpha_R*x[index2(C1RPshort,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RRPshort,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(RRPshort,a,max_age)]
	- ws*x[index2(RRPshort,a,max_age)];
	
	dxdt[ index2(C2APshort , a, max_age )] = age*x[index2(C2APshort,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2APshort,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(RRPshort,a,max_age)]
	- ws*x[index2(C2APshort,a,max_age)];
	
	dxdt[ index2(C2RPshort , a, max_age )] = age*x[index2(C2RPshort,a-1,max_age)]-((1-ksi_R)*ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2RPshort,a,max_age)] + stoc1*(1-delta_R)*FOI_R[a]*x[index2(RAPshort,a,max_age)]
	- ws*x[index2(C2RPshort,a,max_age)];
	
	dxdt[ index2(I2APshort , a, max_age )] = age*x[index2(I2APshort,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C2APshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2APshort,a,max_age)];
	
	dxdt[ index2(I2RPshort , a, max_age )] = age*x[index2(I2RPshort,a-1,max_age)]+(1-ksi_R)*ir_R[a]*x[index2(C2RPshort,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2RPshort,a,max_age)];
	
	dxdt[ index2(RPshort, a, max_age )] = age*x[index2(RPshort,a-1,max_age)]+alpha_A*x[index2(C2APshort,a,max_age)] + alpha_R*x[index2(C2RPshort,a,max_age)] + rho*x[index2(I2APshort,a,max_age)] + rho*x[index2(I2RPshort,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(RPshort,a,max_age)] - ws*x[index2(RPshort,a,max_age)]; 
	 
	 
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , a, max_age )] = ir_A[a]*x[index2(C1A,a,max_age)] + ir_A[a]*x[index2(C2A,a,max_age)];
	dxdt[ index2(IncAV , a, max_age )] = (1-ksi)*ir_A[a]*x[index2(C1AVlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2AVlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C1AVshort,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2AVshort,a,max_age)];
	dxdt[ index2(IncAP , a, max_age )] = (1-ksi)*ir_A[a]*x[index2(C1APlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2APlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C1APshort,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2APshort,a,max_age)];
	dxdt[ index2(IncA_all , a, max_age )] = ir_A[a]*x[index2(C1A,a,max_age)] + ir_A[a]*x[index2(C2A,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C1AVlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2AVlong,a,max_age)] 
		+ (1-ksi)*ir_A[a]*x[index2(C1AVshort,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2AVshort,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C1APlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2APlong,a,max_age)] + 
		(1-ksi)*ir_A[a]*x[index2(C1APshort,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2APshort,a,max_age)];
	
	dxdt[ index2(IncR , a, max_age )] = ir_R[a]*x[index2(C1R,a,max_age)] + ir_R[a]*x[index2(C2R,a,max_age)];
	dxdt[ index2(IncRV , a, max_age )] = ir_R[a]*x[index2(C1RVlong,a,max_age)] + ir_R[a]*x[index2(C2RVlong,a,max_age)] + ir_R[a]*x[index2(C1RVshort,a,max_age)] + ir_R[a]*x[index2(C2RVshort,a,max_age)];
	dxdt[ index2(IncRP , a, max_age )] = (1-ksi_R)*ir_R[a]*x[index2(C1RPlong,a,max_age)] + (1-ksi_R)*ir_R[a]*x[index2(C2RPlong,a,max_age)] + (1-ksi_R)*ir_R[a]*x[index2(C1RPshort,a,max_age)] + (1-ksi_R)*ir_R[a]*x[index2(C2RPshort,a,max_age)];
	dxdt[ index2(IncR_all , a, max_age )] = ir_R[a]*x[index2(C1R,a,max_age)] + ir_R[a]*x[index2(C2R,a,max_age)] + ir_R[a]*x[index2(C1RVlong,a,max_age)] + ir_R[a]*x[index2(C2RVlong,a,max_age)] 
		+ ir_R[a]*x[index2(C1RVshort,a,max_age)] + ir_R[a]*x[index2(C2RVshort,a,max_age)] + (1-ksi_R)*ir_R[a]*x[index2(C1RPlong,a,max_age)] + (1-ksi_R)*ir_R[a]*x[index2(C2RPlong,a,max_age)] + 
		(1-ksi_R)*ir_R[a]*x[index2(C1RPshort,a,max_age)] + (1-ksi_R)*ir_R[a]*x[index2(C2RPshort,a,max_age)];
        
        
        
        
    }
	
    
    // Final Age Group
    dxdt[ index2(S , max_age-1, max_age )] = age*(x[index2(S,(max_age-2),max_age)])+phi*x[index2(RA,max_age-1,max_age)] + phi*x[index2(RR,max_age-1,max_age)] + phi*x[index2(R,max_age-1,max_age)] - stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)]
			 - stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)] - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(S,max_age-1,max_age)]
			 + ws*x[index2(SVshort,max_age-1,max_age)] + wl*x[index2(SVlong,max_age-1,max_age)] + ws*x[index2(SPshort,max_age-1,max_age)] + wl*x[index2(SPlong,max_age-1,max_age)];
    
    dxdt[ index2(C1A , max_age-1, max_age )] = age*(x[index2(C1A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)]
	+ ws*x[index2(C1AVshort,max_age-1,max_age)] + wl*x[index2(C1AVlong,max_age-1,max_age)] + ws*x[index2(C1APshort,max_age-1,max_age)] + wl*x[index2(C1APlong,max_age-1,max_age)];
	
	dxdt[ index2(C1R , max_age-1, max_age )] = age*(x[index2(C1R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)]
	+ ws*x[index2(C1RVshort,max_age-1,max_age)] + wl*x[index2(C1RVlong,max_age-1,max_age)] + ws*x[index2(C1RPshort,max_age-1,max_age)] + wl*x[index2(C1RPlong,max_age-1,max_age)];
    
    dxdt[ index2(I1A , max_age-1, max_age )] = age*(x[index2(I1A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1A,max_age-1,max_age)];
	
	dxdt[ index2(I1R , max_age-1, max_age )] = age*(x[index2(I1R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1R,max_age-1,max_age)];
    
    dxdt[ index2(RA , max_age-1, max_age )] = age*(x[index2(RA,(max_age-2),max_age)])+rho*x[index2(I1A,max_age-1,max_age)]+alpha_A*x[index2(C1A,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RA,max_age-1,max_age)]
	 - stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)] + ws*x[index2(RAVshort,max_age-1,max_age)] + wl*x[index2(RAVlong,max_age-1,max_age)] + ws*x[index2(RAPshort,max_age-1,max_age)] + wl*x[index2(RAPlong,max_age-1,max_age)];
    
	dxdt[ index2(RR , max_age-1, max_age )] = age*(x[index2(RR,(max_age-2),max_age)])+rho*x[index2(I1R,max_age-1,max_age)]+alpha_R*x[index2(C1R,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RR,max_age-1,max_age)]
	 - stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)] + ws*x[index2(RAVshort,max_age-1,max_age)] + wl*x[index2(RAVlong,max_age-1,max_age)] + ws*x[index2(RAPshort,max_age-1,max_age)] + wl*x[index2(RAPlong,max_age-1,max_age)];
	
	dxdt[ index2(C2A , max_age-1, max_age )] = age*(x[index2(C2A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)]
	+ ws*x[index2(C2AVshort,max_age-1,max_age)] + wl*x[index2(C2AVlong,max_age-1,max_age)] + ws*x[index2(C2APshort,max_age-1,max_age)] + wl*x[index2(C2APlong,max_age-1,max_age)];
	
	dxdt[ index2(C2R , max_age-1, max_age )] = age*(x[index2(C2R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)]
	+ ws*x[index2(C2RVshort,max_age-1,max_age)] + wl*x[index2(C2RVlong,max_age-1,max_age)] + ws*x[index2(C2RPshort,max_age-1,max_age)] + wl*x[index2(C2RPlong,max_age-1,max_age)];
	
	dxdt[ index2(I2A , max_age-1, max_age )] = age*(x[index2(I2A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2A,max_age-1,max_age)];
	
	dxdt[ index2(I2R , max_age-1, max_age )] = age*(x[index2(I2R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2R,max_age-1,max_age)];
	
	dxdt[ index2(R, max_age-1, max_age )] = age*(x[index2(R,(max_age-2),max_age)])+alpha_A*x[index2(C2A,max_age-1,max_age)] + alpha_R*x[index2(C2R,max_age-1,max_age)] + rho*x[index2(I2A,max_age-1,max_age)] + rho*x[index2(I2R,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(R,max_age-1,max_age)] + ws*x[index2(RVshort,max_age-1,max_age)] + wl*x[index2(RVlong,max_age-1,max_age)] + ws*x[index2(RPshort,max_age-1,max_age)] + wl*x[index2(RPlong,max_age-1,max_age)]; 
	 
	 	 	// MenAfriVac
	 
	dxdt[ index2(SVlong , max_age-1, max_age )] = age*(x[index2(SVlong,(max_age-2),max_age)])+phi*x[index2(RAVlong,max_age-1,max_age)] + phi*x[index2(RRVlong,max_age-1,max_age)] + phi*x[index2(RVlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SVlong,max_age-1,max_age)] - stoc1*FOI_R[max_age-1]*x[index2(SVlong,max_age-1,max_age)]
	 - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(SVlong,max_age-1,max_age)] - wl*x[index2(SVlong,max_age-1,max_age)];
    
    dxdt[ index2(C1AVlong , max_age-1, max_age )] = age*(x[index2(C1AVlong,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1AVlong,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SVlong,max_age-1,max_age)]
	- wl*x[index2(C1AVlong,max_age-1,max_age)];
	
	dxdt[ index2(C1RVlong , max_age-1, max_age )] = age*(x[index2(C1RVlong,(max_age-2),max_age)])+-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1RVlong,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(SVlong,max_age-1,max_age)]
	- wl*x[index2(C1RVlong,max_age-1,max_age)];
    
    dxdt[ index2(I1AVlong , max_age-1, max_age )] = age*(x[index2(I1AVlong,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C1AVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1AVlong,max_age-1,max_age)];
	
	dxdt[ index2(I1RVlong , max_age-1, max_age )] = age*(x[index2(I1RVlong,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1RVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1RVlong,max_age-1,max_age)];
    
    dxdt[ index2(RAVlong , max_age-1, max_age )] = age*(x[index2(RAVlong,(max_age-2),max_age)])+rho*x[index2(I1AVlong,max_age-1,max_age)]+alpha_A*x[index2(C1AVlong,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RAVlong,max_age-1,max_age)] - stoc1*FOI_R[max_age-1]*x[index2(RAVlong,max_age-1,max_age)]
	- wl*x[index2(RAVlong,max_age-1,max_age)];
    
	dxdt[ index2(RRVlong , max_age-1, max_age )] = age*(x[index2(RRVlong,(max_age-2),max_age)])+rho*x[index2(I1RVlong,max_age-1,max_age)]+alpha_R*x[index2(C1RVlong,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RRVlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRVlong,max_age-1,max_age)]
	- wl*x[index2(RRVlong,max_age-1,max_age)];
	
	dxdt[ index2(C2AVlong , max_age-1, max_age )] = age*(x[index2(C2AVlong,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2AVlong,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRVlong,max_age-1,max_age)]
	- wl*x[index2(C2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(C2RVlong , max_age-1, max_age )] = age*(x[index2(C2RVlong,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2RVlong,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RAVlong,max_age-1,max_age)]
	- wl*x[index2(C2RVlong,max_age-1,max_age)];
	
	dxdt[ index2(I2AVlong , max_age-1, max_age )] = age*(x[index2(I2AVlong,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C2AVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(I2RVlong , max_age-1, max_age )] = age*(x[index2(I2RVlong,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2RVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2RVlong,max_age-1,max_age)];
	
	dxdt[ index2(RVlong, max_age-1, max_age )] = age*(x[index2(RVlong,(max_age-2),max_age)])+alpha_A*x[index2(C2AVlong,max_age-1,max_age)] + alpha_R*x[index2(C2RVlong,max_age-1,max_age)] + rho*x[index2(I2AVlong,max_age-1,max_age)] + rho*x[index2(I2RVlong,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RVlong,max_age-1,max_age)] - wl*x[index2(RVlong,max_age-1,max_age)]; 
   
   
   	dxdt[ index2(SVshort , max_age-1, max_age )] = age*(x[index2(SVshort,(max_age-2),max_age)])+phi*x[index2(RAVshort,max_age-1,max_age)] + phi*x[index2(RRVshort,max_age-1,max_age)] + phi*x[index2(RVshort,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SVshort,max_age-1,max_age)] - stoc1*FOI_R[max_age-1]*x[index2(SVshort,max_age-1,max_age)]
	 - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(SVshort,max_age-1,max_age)] - ws*x[index2(SVshort,max_age-1,max_age)];
    
    dxdt[ index2(C1AVshort , max_age-1, max_age )] = age*(x[index2(C1AVshort,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1AVshort,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SVshort,max_age-1,max_age)]
	- ws*x[index2(C1AVshort,max_age-1,max_age)];
	
	dxdt[ index2(C1RVshort , max_age-1, max_age )] = age*(x[index2(C1RVshort,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1RVshort,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(SVshort,max_age-1,max_age)]
	- ws*x[index2(C1RVshort,max_age-1,max_age)];
    
    dxdt[ index2(I1AVshort , max_age-1, max_age )] = age*(x[index2(I1AVshort,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C1AVshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1AVshort,max_age-1,max_age)];
	
	dxdt[ index2(I1RVshort , max_age-1, max_age )] = age*(x[index2(I1RVshort,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1RVshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1RVshort,max_age-1,max_age)];
    
    dxdt[ index2(RAVshort , max_age-1, max_age )] = age*(x[index2(RAVshort,(max_age-2),max_age)])+rho*x[index2(I1AVshort,max_age-1,max_age)]+alpha_A*x[index2(C1AVshort,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RAVshort,max_age-1,max_age)] - stoc1*FOI_R[max_age-1]*x[index2(RAVshort,max_age-1,max_age)]
	- ws*x[index2(RAVshort,max_age-1,max_age)];
    
	dxdt[ index2(RRVshort , max_age-1, max_age )] = age*(x[index2(RRVshort,(max_age-2),max_age)])+rho*x[index2(I1RVshort,max_age-1,max_age)]+alpha_R*x[index2(C1RVshort,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RRVshort,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRVshort,max_age-1,max_age)]
	- ws*x[index2(RRVshort,max_age-1,max_age)];
	
	dxdt[ index2(C2AVshort , max_age-1, max_age )] = age*(x[index2(C2AVshort,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2AVshort,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRVshort,max_age-1,max_age)]
	- ws*x[index2(C2AVshort,max_age-1,max_age)];
	
	dxdt[ index2(C2RVshort , max_age-1, max_age )] = age*(x[index2(C2RVshort,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2RVshort,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RAVshort,max_age-1,max_age)]
	- ws*x[index2(C2RVshort,max_age-1,max_age)];
	
	dxdt[ index2(I2AVshort , max_age-1, max_age )] = age*(x[index2(I2AVshort,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C2AVshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2AVshort,max_age-1,max_age)];
	
	dxdt[ index2(I2RVshort , max_age-1, max_age )] = age*(x[index2(I2RVshort,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2RVshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2RVshort,max_age-1,max_age)];
	
	dxdt[ index2(RVshort, max_age-1, max_age )] = age*(x[index2(RVshort,(max_age-2),max_age)])+alpha_A*x[index2(C2AVshort,max_age-1,max_age)] + alpha_R*x[index2(C2RVshort,max_age-1,max_age)] + rho*x[index2(I2AVshort,max_age-1,max_age)] + rho*x[index2(I2RVshort,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RVshort,max_age-1,max_age)] - ws*x[index2(RVshort,max_age-1,max_age)]; 
	 
	 
	 // Pentavalent
	 
	 
	dxdt[ index2(SPlong , max_age-1, max_age )] = age*(x[index2(SPlong,(max_age-2),max_age)])+phi*x[index2(RAPlong,max_age-1,max_age)] + phi*x[index2(RRPlong,max_age-1,max_age)] + phi*x[index2(RPlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SPlong,max_age-1,max_age)] - stoc1*(1-delta_R)*FOI_R[max_age-1]*x[index2(SPlong,max_age-1,max_age)]
	 - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(SPlong,max_age-1,max_age)] - wl*x[index2(SPlong,max_age-1,max_age)];
    
    dxdt[ index2(C1APlong , max_age-1, max_age )] = age*(x[index2(C1APlong,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1APlong,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SPlong,max_age-1,max_age)]
	- wl*x[index2(C1APlong,max_age-1,max_age)];
	
	dxdt[ index2(C1RPlong , max_age-1, max_age )] = age*(x[index2(C1RPlong,(max_age-2),max_age)])-((1-ksi_R)*ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1RPlong,max_age-1,max_age)] + stoc1*(1-delta_R)*FOI_R[max_age-1]*x[index2(SPlong,max_age-1,max_age)]
	- wl*x[index2(C1RPlong,max_age-1,max_age)];
    
    dxdt[ index2(I1APlong , max_age-1, max_age )] = age*(x[index2(I1APlong,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C1APlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1APlong,max_age-1,max_age)];
	
	dxdt[ index2(I1RPlong , max_age-1, max_age )] = age*(x[index2(I1RPlong,(max_age-2),max_age)])+(1-ksi_R)*ir_R[max_age-1]*x[index2(C1RPlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1RPlong,max_age-1,max_age)];
    
    dxdt[ index2(RAPlong , max_age-1, max_age )] = age*(x[index2(RAPlong,(max_age-2),max_age)])+rho*x[index2(I1APlong,max_age-1,max_age)]+alpha_A*x[index2(C1APlong,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RAPlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_R[max_age-1]*x[index2(RAPlong,max_age-1,max_age)]
	- wl*x[index2(RAPlong,max_age-1,max_age)];
    
	dxdt[ index2(RRPlong , max_age-1, max_age )] = age*(x[index2(RRPlong,(max_age-2),max_age)])+rho*x[index2(I1RPlong,max_age-1,max_age)]+alpha_R*x[index2(C1RPlong,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RRPlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRPlong,max_age-1,max_age)]
	- wl*x[index2(RRPlong,max_age-1,max_age)];
	
	dxdt[ index2(C2APlong , max_age-1, max_age )] = age*(x[index2(C2APlong,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2APlong,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRPlong,max_age-1,max_age)]
	- wl*x[index2(C2APlong,max_age-1,max_age)];
	
	dxdt[ index2(C2RPlong , max_age-1, max_age )] = age*(x[index2(C2RPlong,(max_age-2),max_age)])-((1-ksi_R)*ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2RPlong,max_age-1,max_age)] + stoc1*(1-delta_R)*FOI_R[max_age-1]*x[index2(RAPlong,max_age-1,max_age)]
	- wl*x[index2(C2RPlong,max_age-1,max_age)];
	
	dxdt[ index2(I2APlong , max_age-1, max_age )] = age*(x[index2(I2APlong,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C2APlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2APlong,max_age-1,max_age)];
	
	dxdt[ index2(I2RPlong , max_age-1, max_age )] = age*(x[index2(I2RPlong,(max_age-2),max_age)])+(1-ksi_R)*ir_R[max_age-1]*x[index2(C2RPlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2RPlong,max_age-1,max_age)];
	
	dxdt[ index2(RPlong, max_age-1, max_age )] = age*(x[index2(RPlong,(max_age-2),max_age)])+alpha_A*x[index2(C2APlong,max_age-1,max_age)] + alpha_R*x[index2(C2RPlong,max_age-1,max_age)] + rho*x[index2(I2APlong,max_age-1,max_age)] + rho*x[index2(I2RPlong,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RPlong,max_age-1,max_age)] - wl*x[index2(RPlong,max_age-1,max_age)]; 
   
   
   	dxdt[ index2(SPshort , max_age-1, max_age )] = age*(x[index2(SPshort,(max_age-2),max_age)])+phi*x[index2(RAPshort,max_age-1,max_age)] + phi*x[index2(RRPshort,max_age-1,max_age)] + phi*x[index2(RPshort,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SPshort,max_age-1,max_age)] - stoc1*(1-delta_R)*FOI_R[max_age-1]*x[index2(SPshort,max_age-1,max_age)]
	 - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(SPshort,max_age-1,max_age)] - ws*x[index2(SPshort,max_age-1,max_age)];
    
    dxdt[ index2(C1APshort , max_age-1, max_age )] = age*(x[index2(C1APshort,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1APshort,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SPshort,max_age-1,max_age)]
	- ws*x[index2(C1APshort,max_age-1,max_age)];
	
	dxdt[ index2(C1RPshort , max_age-1, max_age )] = age*(x[index2(C1RPshort,(max_age-2),max_age)])-((1-ksi_R)*ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1RPshort,max_age-1,max_age)] + stoc1*(1-delta_R)*FOI_R[max_age-1]*x[index2(SPshort,max_age-1,max_age)]
	- ws*x[index2(C1RPshort,max_age-1,max_age)];
    
    dxdt[ index2(I1APshort , max_age-1, max_age )] = age*(x[index2(I1APshort,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C1APshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1APshort,max_age-1,max_age)];
	
	dxdt[ index2(I1RPshort , max_age-1, max_age )] = age*(x[index2(I1RPshort,(max_age-2),max_age)])+(1-ksi_R)*ir_R[max_age-1]*x[index2(C1RPshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1RPshort,max_age-1,max_age)];
    
    dxdt[ index2(RAPshort , max_age-1, max_age )] = age*(x[index2(RAPshort,(max_age-2),max_age)])+rho*x[index2(I1APshort,max_age-1,max_age)]+alpha_A*x[index2(C1APshort,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RAPshort,max_age-1,max_age)] - stoc1*(1-delta)*FOI_R[max_age-1]*x[index2(RAPshort,max_age-1,max_age)]
	- ws*x[index2(RAPshort,max_age-1,max_age)];
    
	dxdt[ index2(RRPshort , max_age-1, max_age )] = age*(x[index2(RRPshort,(max_age-2),max_age)])+rho*x[index2(I1RPshort,max_age-1,max_age)]+alpha_R*x[index2(C1RPshort,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RRPshort,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRPshort,max_age-1,max_age)]
	- ws*x[index2(RRPshort,max_age-1,max_age)];
	
	dxdt[ index2(C2APshort , max_age-1, max_age )] = age*(x[index2(C2APshort,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2APshort,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRPshort,max_age-1,max_age)]
	- ws*x[index2(C2APshort,max_age-1,max_age)];
	
	dxdt[ index2(C2RPshort , max_age-1, max_age )] = age*(x[index2(C2RPshort,(max_age-2),max_age)])-((1-ksi_R)*ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2RPshort,max_age-1,max_age)] + stoc1*(1-delta_R)*FOI_R[max_age-1]*x[index2(RAPshort,max_age-1,max_age)]
	- ws*x[index2(C2RPshort,max_age-1,max_age)];
	
	dxdt[ index2(I2APshort , max_age-1, max_age )] = age*(x[index2(I2APshort,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C2APshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2APshort,max_age-1,max_age)];
	
	dxdt[ index2(I2RPshort , max_age-1, max_age )] = age*(x[index2(I2RPshort,(max_age-2),max_age)])+(1-ksi_R)*ir_R[max_age-1]*x[index2(C2RPshort,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2RPshort,max_age-1,max_age)];
	
	dxdt[ index2(RPshort, max_age-1, max_age )] = age*(x[index2(RPshort,(max_age-2),max_age)])+alpha_A*x[index2(C2APshort,max_age-1,max_age)] + alpha_R*x[index2(C2RPshort,max_age-1,max_age)] + rho*x[index2(I2APshort,max_age-1,max_age)] + rho*x[index2(I2RPshort,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RPshort,max_age-1,max_age)] - ws*x[index2(RPshort,max_age-1,max_age)]; 
    
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , max_age-1, max_age )] = ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)] + ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)];
	dxdt[ index2(IncAV , max_age-1, max_age )] = (1-ksi)*ir_A[max_age-1]*x[index2(C1AVlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2AVlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C1AVshort,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2AVshort,max_age-1,max_age)];
	dxdt[ index2(IncAP , max_age-1, max_age )] = (1-ksi)*ir_A[max_age-1]*x[index2(C1APlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2APlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C1APshort,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2APshort,max_age-1,max_age)];
	dxdt[ index2(IncA_all , max_age-1, max_age )] = ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)] + ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C1AVlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2AVlong,max_age-1,max_age)] 
		+ (1-ksi)*ir_A[max_age-1]*x[index2(C1AVshort,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2AVshort,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C1APlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2APlong,max_age-1,max_age)] + 
		(1-ksi)*ir_A[max_age-1]*x[index2(C1APshort,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2APshort,max_age-1,max_age)];
	
	dxdt[ index2(IncR , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)];
	dxdt[ index2(IncRV , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1RVlong,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2RVlong,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C1RVshort,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2RVshort,max_age-1,max_age)];
	dxdt[ index2(IncRP , max_age-1, max_age )] = (1-ksi_R)*ir_R[max_age-1]*x[index2(C1RPlong,max_age-1,max_age)] + (1-ksi_R)*ir_R[max_age-1]*x[index2(C2RPlong,max_age-1,max_age)] + (1-ksi_R)*ir_R[max_age-1]*x[index2(C1RPshort,max_age-1,max_age)] + (1-ksi_R)*ir_R[max_age-1]*x[index2(C2RPshort,max_age-1,max_age)];
	dxdt[ index2(IncR_all , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C1RVlong,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2RVlong,max_age-1,max_age)] 
		+ ir_R[max_age-1]*x[index2(C1RVshort,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2RVshort,max_age-1,max_age)] + (1-ksi_R)*ir_R[max_age-1]*x[index2(C1RPlong,max_age-1,max_age)] + (1-ksi_R)*ir_R[max_age-1]*x[index2(C2RPlong,max_age-1,max_age)] + 
		(1-ksi_R)*ir_R[max_age-1]*x[index2(C1RPshort,max_age-1,max_age)] + (1-ksi_R)*ir_R[max_age-1]*x[index2(C2RPshort,max_age-1,max_age)];
        
    
    return List::create(dxdt, stoc1,birth, N);
}




// [[Rcpp::export]]


List PentaSt(double t, NumericVector x, NumericVector params, NumericMatrix Cm, NumericVector mu, double max_age, double B, NumericVector ir_A, NumericVector ir_R, NumericVector stoc)
{
    
    NumericVector dxdt(x.length());
    
    double phi = params[0]; //loss of immunity
    double rho = params[1]; // recovery of disease
    double alpha_A = params[2]; // duration of carriage menA
	double alpha_R = params[3]; // duration of carriage menCWYX
    double epsbeta = params[4]; // seasonal amplitude
    double beta_A = params[5]; // contact parameter
    double beta_R = params[6];
	
    double pi = M_PI;
    
    double age = 1.0/365.0; // Daily Ageing
    
    int comp=14;
    
    enum  state_variables {S,C1A,C1R,I1A,I1R,RA,RR,C2A,C2R,I2A,I2R,R,IncA,IncR};
    
    NumericVector FOI_A(max_age);
	NumericVector FOI_R(max_age);
    
    double stoc1 = stoc[t/365];
    
  //  double birth = B[t/365];
    
    double N = 0.0;
    for( int j=0; j<(comp-1)*max_age; j++)
    {
        N += x[j];
    }
    
    // Calculate force of infection
    for(int a=0 ; a<max_age ; a++ )
    {
        FOI_A[a] = 0.0;
		FOI_R[a] = 0.0;
		
        for( int a2=0; a2<max_age; a2++ )
        {
            
            double Z = (1+epsbeta*cos(2*pi*t/365.0));
			
            FOI_A[a] += beta_A*Z*Cm(a,a2)*(x[index2(C1A,a2,max_age)] + x[index2(C2A,a2,max_age)] + x[index2(I1A,a2,max_age)] + x[index2(I2A,a2,max_age)])/N;
			
			FOI_R[a] += beta_R*Z*Cm(a,a2)*(x[index2(C1R,a2,max_age)] + x[index2(C2R,a2,max_age)] + x[index2(I1R,a2,max_age)] + x[index2(I2R,a2,max_age)])/N;
        }
    }
    
    
    
    dxdt[ index2(S , 0, max_age )] = B + phi*x[index2(RA,0,max_age)] + phi*x[index2(RR,0,max_age)] + phi*x[index2(R,0,max_age)] - stoc1*FOI_A[0]*x[index2(S,0,max_age)] - stoc1*FOI_R[0]*x[index2(S,0,max_age)]
	 - (age+mu[0])*x[index2(S,0,max_age)];
    
    dxdt[ index2(C1A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[0]+age)*x[index2(C1A,0,max_age)] + stoc1*FOI_A[0]*x[index2(S,0,max_age)];
	
	dxdt[ index2(C1R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[0]+age)*x[index2(C1R,0,max_age)] + stoc1*FOI_R[0]*x[index2(S,0,max_age)];
    
    dxdt[ index2(I1A , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)]-(rho+age+mu[0])*x[index2(I1A,0,max_age)];
	
	dxdt[ index2(I1R , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)]-(rho+age+mu[0])*x[index2(I1R,0,max_age)];
    
    dxdt[ index2(RA , 0, max_age )] = rho*x[index2(I1A,0,max_age)]+alpha_A*x[index2(C1A,0,max_age)]-(phi+age+mu[0])*x[index2(RA,0,max_age)] - stoc1*FOI_R[0]*x[index2(RA,0,max_age)];
    
	dxdt[ index2(RR , 0, max_age )] = rho*x[index2(I1R,0,max_age)]+alpha_R*x[index2(C1R,0,max_age)]-(phi+age+mu[0])*x[index2(RR,0,max_age)] - stoc1*FOI_A[0]*x[index2(RR,0,max_age)];
	
	dxdt[ index2(C2A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[0]+age)*x[index2(C2A,0,max_age)] + stoc1*FOI_A[0]*x[index2(RR,0,max_age)];
	
	dxdt[ index2(C2R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[0]+age)*x[index2(C2R,0,max_age)] + stoc1*FOI_R[0]*x[index2(RA,0,max_age)];
	
	dxdt[ index2(I2A , 0, max_age )] = ir_A[0]*x[index2(C2A,0,max_age)]-(rho+age+mu[0])*x[index2(I2A,0,max_age)];
	
	dxdt[ index2(I2R , 0, max_age )] = ir_R[0]*x[index2(C2R,0,max_age)]-(rho+age+mu[0])*x[index2(I2R,0,max_age)];
	
	dxdt[ index2(R, 0, max_age )] = alpha_A*x[index2(C2A,0,max_age)] + alpha_R*x[index2(C2R,0,max_age)] + rho*x[index2(I2A,0,max_age)] + rho*x[index2(I2R,0,max_age)]
	 - (phi+mu[0]+age)*x[index2(R,0,max_age)]; 
	
	 
   
   // Incidence
    
    dxdt[ index2(IncA , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)] + ir_A[0]*x[index2(C2A,0,max_age)];
	
	dxdt[ index2(IncR , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)] + ir_R[0]*x[index2(C2R,0,max_age)];
	
    
	
	//Second age group
	
	    
    dxdt[ index2(S , 1, max_age )] = age*x[index2(S,0,max_age)]+phi*x[index2(RA,1,max_age)] + phi*x[index2(RR,1,max_age)] + phi*x[index2(R,1,max_age)] - stoc1*FOI_A[1]*x[index2(S,1,max_age)] - stoc1*FOI_R[1]*x[index2(S,1,max_age)]
	 - (age+mu[1])*x[index2(S,1,max_age)];
    
    dxdt[ index2(C1A , 1, max_age )] = age*x[index2(C1A,0,max_age)]-(ir_A[1]+alpha_A+mu[1]+age)*x[index2(C1A,1,max_age)] + stoc1*FOI_A[1]*x[index2(S,1,max_age)];
	
	dxdt[ index2(C1R , 1, max_age )] = age*x[index2(C1R,0,max_age)]-(ir_R[1]+alpha_R+mu[1]+age)*x[index2(C1R,1,max_age)] + stoc1*FOI_R[1]*x[index2(S,1,max_age)];
    
    dxdt[ index2(I1A , 1, max_age )] = age*x[index2(I1A,0,max_age)]+ir_A[1]*x[index2(C1A,1,max_age)]-(rho+age+mu[1])*x[index2(I1A,1,max_age)];
	
	dxdt[ index2(I1R , 1, max_age )] = age*x[index2(I1R,0,max_age)]+ir_R[1]*x[index2(C1R,1,max_age)]-(rho+age+mu[1])*x[index2(I1R,1,max_age)];
    
    dxdt[ index2(RA , 1, max_age )] = age*x[index2(RA,0,max_age)]+rho*x[index2(I1A,1,max_age)]+alpha_A*x[index2(C1A,1,max_age)]-(phi+age+mu[1])*x[index2(RA,1,max_age)] - stoc1*FOI_R[1]*x[index2(RA,1,max_age)];
    
	dxdt[ index2(RR , 1, max_age )] = age*x[index2(RR,0,max_age)]+rho*x[index2(I1R,1,max_age)]+alpha_R*x[index2(C1R,1,max_age)]-(phi+age+mu[1])*x[index2(RR,1,max_age)] - stoc1*FOI_A[1]*x[index2(RR,1,max_age)];
	
	dxdt[ index2(C2A , 1, max_age )] = age*x[index2(C2A,0,max_age)]-(ir_A[1]+alpha_A+mu[1]+age)*x[index2(C2A,1,max_age)] + stoc1*FOI_A[1]*x[index2(RR,1,max_age)];
	
	dxdt[ index2(C2R , 1, max_age )] = age*x[index2(C2R,0,max_age)]-(ir_R[1]+alpha_R+mu[1]+age)*x[index2(C2R,1,max_age)] + stoc1*FOI_R[1]*x[index2(RA,1,max_age)];
	
	dxdt[ index2(I2A , 1, max_age )] = age*x[index2(I2A,0,max_age)]+ir_A[1]*x[index2(C2A,1,max_age)]-(rho+age+mu[1])*x[index2(I2A,1,max_age)];
	
	dxdt[ index2(I2R , 1, max_age )] = age*x[index2(I2R,0,max_age)]+ir_R[1]*x[index2(C2R,1,max_age)]-(rho+age+mu[1])*x[index2(I2R,1,max_age)];
	
	dxdt[ index2(R, 1, max_age )] = age*x[index2(R,0,max_age)]+alpha_A*x[index2(C2A,1,max_age)] + alpha_R*x[index2(C2R,1,max_age)] + rho*x[index2(I2A,1,max_age)] + rho*x[index2(I2R,1,max_age)]
	 - (phi+mu[1]+age)*x[index2(R,1,max_age)]; 
	
	 // Incidence
	   
	 dxdt[ index2(IncA , 1, max_age )] = ir_A[1]*x[index2(C1A,1,max_age)] + ir_A[1]*x[index2(C2A,1,max_age)];
	
	dxdt[ index2(IncR , 1, max_age )] = ir_R[1]*x[index2(C1R,1,max_age)] + ir_R[1]*x[index2(C2R,1,max_age)];
	 
    
    // Annual age groups
    for( size_t a=2 ; a<(max_age-1) ; ++a )
    {
        
            dxdt[ index2(S , a, max_age )] = age*(x[index2(S,a-1,max_age)]-x[index2(S,a,max_age)]) + phi*x[index2(RA,a,max_age)] + phi*x[index2(RR,a,max_age)] + phi*x[index2(R,a,max_age)] - stoc1*FOI_A[a]*x[index2(S,a,max_age)]
			 - stoc1*FOI_R[a]*x[index2(S,a,max_age)] - (mu[a])*x[index2(S,a,max_age)];
    
    dxdt[ index2(C1A , a, max_age )] = age*(x[index2(C1A,a-1,max_age)]-x[index2(C1A,a,max_age)]) -(ir_A[a]+alpha_A+mu[a])*x[index2(C1A,a,max_age)] + stoc1*FOI_A[a]*x[index2(S,a,max_age)];
	
	dxdt[ index2(C1R , a, max_age )] = age*(x[index2(C1R,a-1,max_age)]-x[index2(C1R,a,max_age)]) -(ir_R[a]+alpha_R+mu[a])*x[index2(C1R,a,max_age)] + stoc1*FOI_R[a]*x[index2(S,a,max_age)];
    
    dxdt[ index2(I1A , a, max_age )] = age*(x[index2(I1A,a-1,max_age)]-x[index2(I1A,a,max_age)]) +ir_A[a]*x[index2(C1A,a,max_age)]-(rho+mu[a])*x[index2(I1A,a,max_age)];
	
	dxdt[ index2(I1R , a, max_age )] = age*(x[index2(I1R,a-1,max_age)]-x[index2(I1R,a,max_age)]) +ir_R[a]*x[index2(C1R,a,max_age)]-(rho+mu[a])*x[index2(I1R,a,max_age)];
    
    dxdt[ index2(RA , a, max_age )] = age*(x[index2(RA,a-1,max_age)]-x[index2(RA,a,max_age)]) +rho*x[index2(I1A,a,max_age)]+alpha_A*x[index2(C1A,a,max_age)]-(phi+mu[a])*x[index2(RA,a,max_age)]
	 - stoc1*FOI_R[a]*x[index2(RA,a,max_age)];
    
	dxdt[ index2(RR , a, max_age )] = age*(x[index2(RR,a-1,max_age)]-x[index2(RR,a,max_age)]) +rho*x[index2(I1R,a,max_age)]+alpha_R*x[index2(C1R,a,max_age)]-(phi+mu[a])*x[index2(RR,a,max_age)]
	 - stoc1*FOI_A[a]*x[index2(RR,a,max_age)];
	
	dxdt[ index2(C2A , a, max_age )] = age*(x[index2(C2A,a-1,max_age)]-x[index2(C2A,a,max_age)]) -(ir_A[a]+alpha_A+mu[a])*x[index2(C2A,a,max_age)] + stoc1*FOI_A[a]*x[index2(RR,a,max_age)];
	
	dxdt[ index2(C2R , a, max_age )] = age*(x[index2(C2R,a-1,max_age)]-x[index2(C2R,a,max_age)]) -(ir_R[a]+alpha_R+mu[a])*x[index2(C2R,a,max_age)] + stoc1*FOI_R[a]*x[index2(RA,a,max_age)];
	
	dxdt[ index2(I2A , a, max_age )] = age*(x[index2(I2A,a-1,max_age)]-x[index2(I2A,a,max_age)]) +ir_A[a]*x[index2(C2A,a,max_age)]-(rho+mu[a])*x[index2(I2A,a,max_age)];
	
	dxdt[ index2(I2R , a, max_age )] = age*(x[index2(I2R,a-1,max_age)]-x[index2(I2R,a,max_age)]) +ir_R[a]*x[index2(C2R,a,max_age)]-(rho+mu[a])*x[index2(I2R,a,max_age)];
	
	dxdt[ index2(R, a, max_age )] = age*(x[index2(R,a-1,max_age)]-x[index2(R,a,max_age)]) +alpha_A*x[index2(C2A,a,max_age)] + alpha_R*x[index2(C2R,a,max_age)] + rho*x[index2(I2A,a,max_age)] + rho*x[index2(I2R,a,max_age)]
	 - (phi+mu[a])*x[index2(R,a,max_age)]; 
	 
	 
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , a, max_age )] = ir_A[a]*x[index2(C1A,a,max_age)] + ir_A[a]*x[index2(C2A,a,max_age)];
	
	dxdt[ index2(IncR , a, max_age )] = ir_R[a]*x[index2(C1R,a,max_age)] + ir_R[a]*x[index2(C2R,a,max_age)];
	    
        
    }
	
    
    // Final Age Group
    dxdt[ index2(S , max_age-1, max_age )] = age*(x[index2(S,(max_age-2),max_age)])+phi*x[index2(RA,max_age-1,max_age)] + phi*x[index2(RR,max_age-1,max_age)] + phi*x[index2(R,max_age-1,max_age)] - stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)]
			 - stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)] - (mu[max_age-1])*x[index2(S,max_age-1,max_age)];
    
    dxdt[ index2(C1A , max_age-1, max_age )] = age*(x[index2(C1A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[max_age-1])*x[index2(C1A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)];
	
	dxdt[ index2(C1R , max_age-1, max_age )] = age*(x[index2(C1R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[max_age-1])*x[index2(C1R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)];
    
    dxdt[ index2(I1A , max_age-1, max_age )] = age*(x[index2(I1A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)]-(rho+mu[max_age-1])*x[index2(I1A,max_age-1,max_age)];
	
	dxdt[ index2(I1R , max_age-1, max_age )] = age*(x[index2(I1R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)]-(rho+mu[max_age-1])*x[index2(I1R,max_age-1,max_age)];
    
    dxdt[ index2(RA , max_age-1, max_age )] = age*(x[index2(RA,(max_age-2),max_age)])+rho*x[index2(I1A,max_age-1,max_age)]+alpha_A*x[index2(C1A,max_age-1,max_age)]-(phi+mu[max_age-1])*x[index2(RA,max_age-1,max_age)]
	 - stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)];
    
	dxdt[ index2(RR , max_age-1, max_age )] = age*(x[index2(RR,(max_age-2),max_age)])+rho*x[index2(I1R,max_age-1,max_age)]+alpha_R*x[index2(C1R,max_age-1,max_age)]-(phi+mu[max_age-1])*x[index2(RR,max_age-1,max_age)]
	 - stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)];
	
	dxdt[ index2(C2A , max_age-1, max_age )] = age*(x[index2(C2A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[max_age-1])*x[index2(C2A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)];
	
	dxdt[ index2(C2R , max_age-1, max_age )] = age*(x[index2(C2R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[max_age-1])*x[index2(C2R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)];
	
	dxdt[ index2(I2A , max_age-1, max_age )] = age*(x[index2(I2A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)]-(rho+mu[max_age-1])*x[index2(I2A,max_age-1,max_age)];
	
	dxdt[ index2(I2R , max_age-1, max_age )] = age*(x[index2(I2R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)]-(rho+mu[max_age-1])*x[index2(I2R,max_age-1,max_age)];
	
	dxdt[ index2(R, max_age-1, max_age )] = age*(x[index2(R,(max_age-2),max_age)])+alpha_A*x[index2(C2A,max_age-1,max_age)] + alpha_R*x[index2(C2R,max_age-1,max_age)] + rho*x[index2(I2A,max_age-1,max_age)] + rho*x[index2(I2R,max_age-1,max_age)]
	 - (phi+mu[max_age-1])*x[index2(R,max_age-1,max_age)]; 
	 
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , max_age-1, max_age )] = ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)] + ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)];

	dxdt[ index2(IncR , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)];
        
    
    return List::create(dxdt, stoc1,B, N);
}





// [[Rcpp::export]]


List PentaDyn(double t, NumericVector x, NumericVector params, NumericMatrix Cm, NumericVector mu, double max_age, NumericVector B, NumericVector ir_A, NumericVector ir_R, NumericVector stoc)
{
    
    NumericVector dxdt(x.length());
    
    double phi = params[0]; //loss of immunity
    double rho = params[1]; // recovery of disease
    double alpha_A = params[2]; // duration of carriage menA
	double alpha_R = params[3]; // duration of carriage menCWYX
    double epsbeta = params[4]; // seasonal amplitude
    double beta_A = params[5]; // contact parameter
    double beta_R = params[6];
	
    double pi = M_PI;
    
    double age = 1.0/365.0; // Daily Ageing
    
    int comp=14;
    
    enum  state_variables {S,C1A,C1R,I1A,I1R,RA,RR,C2A,C2R,I2A,I2R,R,IncA,IncR};
    
    NumericVector FOI_A(max_age);
	NumericVector FOI_R(max_age);
    
    double stoc1 = stoc[t/365];
    
    double birth = B[t/365];
    
    double N = 0.0;
    for( int j=0; j<(comp-2)*max_age; j++)
    {
        N += x[j];
    }
    
    // Calculate force of infection
    for(int a=0 ; a<max_age ; a++ )
    {
        FOI_A[a] = 0.0;
		FOI_R[a] = 0.0;
		
        for( int a2=0; a2<max_age; a2++ )
        {
            
            double Z = (1+epsbeta*cos(2*pi*t/365.0));
			
            FOI_A[a] += beta_A*Z*Cm(a,a2)*(x[index2(C1A,a2,max_age)] + x[index2(C2A,a2,max_age)] + x[index2(I1A,a2,max_age)] + x[index2(I2A,a2,max_age)])/N;
			
			FOI_R[a] += beta_R*Z*Cm(a,a2)*(x[index2(C1R,a2,max_age)] + x[index2(C2R,a2,max_age)] + x[index2(I1R,a2,max_age)] + x[index2(I2R,a2,max_age)])/N;
        }
    }
    
    
    
    dxdt[ index2(S , 0, max_age )] = birth*N + phi*x[index2(RA,0,max_age)] + phi*x[index2(RR,0,max_age)] + phi*x[index2(R,0,max_age)] - stoc1*FOI_A[0]*x[index2(S,0,max_age)] - stoc1*FOI_R[0]*x[index2(S,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(S,0,max_age)];
    
    dxdt[ index2(C1A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1A,0,max_age)] + stoc1*FOI_A[0]*x[index2(S,0,max_age)];
	
	dxdt[ index2(C1R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1R,0,max_age)] + stoc1*FOI_R[0]*x[index2(S,0,max_age)];
    
    dxdt[ index2(I1A , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1A,0,max_age)];
	
	dxdt[ index2(I1R , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1R,0,max_age)];
    
    dxdt[ index2(RA , 0, max_age )] = rho*x[index2(I1A,0,max_age)]+alpha_A*x[index2(C1A,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RA,0,max_age)] - stoc1*FOI_R[0]*x[index2(RA,0,max_age)];
    
	dxdt[ index2(RR , 0, max_age )] = rho*x[index2(I1R,0,max_age)]+alpha_R*x[index2(C1R,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RR,0,max_age)] - stoc1*FOI_A[0]*x[index2(RR,0,max_age)];
	
	dxdt[ index2(C2A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2A,0,max_age)] + stoc1*FOI_A[0]*x[index2(RR,0,max_age)];
	
	dxdt[ index2(C2R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2R,0,max_age)] + stoc1*FOI_R[0]*x[index2(RA,0,max_age)];
	
	dxdt[ index2(I2A , 0, max_age )] = ir_A[0]*x[index2(C2A,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2A,0,max_age)];
	
	dxdt[ index2(I2R , 0, max_age )] = ir_R[0]*x[index2(C2R,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2R,0,max_age)];
	
	dxdt[ index2(R, 0, max_age )] = alpha_A*x[index2(C2A,0,max_age)] + alpha_R*x[index2(C2R,0,max_age)] + rho*x[index2(I2A,0,max_age)] + rho*x[index2(I2R,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(R,0,max_age)]; 
	
	 
   
   // Incidence
    
    dxdt[ index2(IncA , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)] + ir_A[0]*x[index2(C2A,0,max_age)];
	
	dxdt[ index2(IncR , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)] + ir_R[0]*x[index2(C2R,0,max_age)];
	
    
	
	//Second age group
	
	    
    dxdt[ index2(S , 1, max_age )] = age*x[index2(S,0,max_age)]+phi*x[index2(RA,1,max_age)] + phi*x[index2(RR,1,max_age)] + phi*x[index2(R,1,max_age)] - stoc1*FOI_A[1]*x[index2(S,1,max_age)] - stoc1*FOI_R[1]*x[index2(S,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(S,1,max_age)];
    
    dxdt[ index2(C1A , 1, max_age )] = age*x[index2(C1A,0,max_age)]-(ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1A,1,max_age)] + stoc1*FOI_A[1]*x[index2(S,1,max_age)];
	
	dxdt[ index2(C1R , 1, max_age )] = age*x[index2(C1R,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1R,1,max_age)] + stoc1*FOI_R[1]*x[index2(S,1,max_age)];
    
    dxdt[ index2(I1A , 1, max_age )] = age*x[index2(I1A,0,max_age)]+ir_A[1]*x[index2(C1A,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1A,1,max_age)];
	
	dxdt[ index2(I1R , 1, max_age )] = age*x[index2(I1R,0,max_age)]+ir_R[1]*x[index2(C1R,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1R,1,max_age)];
    
    dxdt[ index2(RA , 1, max_age )] = age*x[index2(RA,0,max_age)]+rho*x[index2(I1A,1,max_age)]+alpha_A*x[index2(C1A,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RA,1,max_age)] - stoc1*FOI_R[1]*x[index2(RA,1,max_age)];
    
	dxdt[ index2(RR , 1, max_age )] = age*x[index2(RR,0,max_age)]+rho*x[index2(I1R,1,max_age)]+alpha_R*x[index2(C1R,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RR,1,max_age)] - stoc1*FOI_A[1]*x[index2(RR,1,max_age)];
	
	dxdt[ index2(C2A , 1, max_age )] = age*x[index2(C2A,0,max_age)]-(ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2A,1,max_age)] + stoc1*FOI_A[1]*x[index2(RR,1,max_age)];
	
	dxdt[ index2(C2R , 1, max_age )] = age*x[index2(C2R,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2R,1,max_age)] + stoc1*FOI_R[1]*x[index2(RA,1,max_age)];
	
	dxdt[ index2(I2A , 1, max_age )] = age*x[index2(I2A,0,max_age)]+ir_A[1]*x[index2(C2A,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2A,1,max_age)];
	
	dxdt[ index2(I2R , 1, max_age )] = age*x[index2(I2R,0,max_age)]+ir_R[1]*x[index2(C2R,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2R,1,max_age)];
	
	dxdt[ index2(R, 1, max_age )] = age*x[index2(R,0,max_age)]+alpha_A*x[index2(C2A,1,max_age)] + alpha_R*x[index2(C2R,1,max_age)] + rho*x[index2(I2A,1,max_age)] + rho*x[index2(I2R,1,max_age)]
	 - (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(R,1,max_age)]; 
	
	 // Incidence
	   
	 dxdt[ index2(IncA , 1, max_age )] = ir_A[1]*x[index2(C1A,1,max_age)] + ir_A[1]*x[index2(C2A,1,max_age)];
	
	dxdt[ index2(IncR , 1, max_age )] = ir_R[1]*x[index2(C1R,1,max_age)] + ir_R[1]*x[index2(C2R,1,max_age)];
	 
    
    // Annual age groups
    for( size_t a=2 ; a<(max_age-1) ; ++a )
    {
        
            dxdt[ index2(S , a, max_age )] = age*(x[index2(S,a-1,max_age)]-x[index2(S,a,max_age)]) + phi*x[index2(RA,a,max_age)] + phi*x[index2(RR,a,max_age)] + phi*x[index2(R,a,max_age)] - stoc1*FOI_A[a]*x[index2(S,a,max_age)]
			 - stoc1*FOI_R[a]*x[index2(S,a,max_age)] - (mu[index2(floor(t/365),a,max_age)])*x[index2(S,a,max_age)];
    
    dxdt[ index2(C1A , a, max_age )] = age*(x[index2(C1A,a-1,max_age)]-x[index2(C1A,a,max_age)]) -(ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)])*x[index2(C1A,a,max_age)] + stoc1*FOI_A[a]*x[index2(S,a,max_age)];
	
	dxdt[ index2(C1R , a, max_age )] = age*(x[index2(C1R,a-1,max_age)]-x[index2(C1R,a,max_age)]) -(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)])*x[index2(C1R,a,max_age)] + stoc1*FOI_R[a]*x[index2(S,a,max_age)];
    
    dxdt[ index2(I1A , a, max_age )] = age*(x[index2(I1A,a-1,max_age)]-x[index2(I1A,a,max_age)]) +ir_A[a]*x[index2(C1A,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I1A,a,max_age)];
	
	dxdt[ index2(I1R , a, max_age )] = age*(x[index2(I1R,a-1,max_age)]-x[index2(I1R,a,max_age)]) +ir_R[a]*x[index2(C1R,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I1R,a,max_age)];
    
    dxdt[ index2(RA , a, max_age )] = age*(x[index2(RA,a-1,max_age)]-x[index2(RA,a,max_age)]) +rho*x[index2(I1A,a,max_age)]+alpha_A*x[index2(C1A,a,max_age)]-(phi+mu[index2(floor(t/365),a,max_age)])*x[index2(RA,a,max_age)]
	 - stoc1*FOI_R[a]*x[index2(RA,a,max_age)];
    
	dxdt[ index2(RR , a, max_age )] = age*(x[index2(RR,a-1,max_age)]-x[index2(RR,a,max_age)]) +rho*x[index2(I1R,a,max_age)]+alpha_R*x[index2(C1R,a,max_age)]-(phi+mu[index2(floor(t/365),a,max_age)])*x[index2(RR,a,max_age)]
	 - stoc1*FOI_A[a]*x[index2(RR,a,max_age)];
	
	dxdt[ index2(C2A , a, max_age )] = age*(x[index2(C2A,a-1,max_age)]-x[index2(C2A,a,max_age)]) -(ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)])*x[index2(C2A,a,max_age)] + stoc1*FOI_A[a]*x[index2(RR,a,max_age)];
	
	dxdt[ index2(C2R , a, max_age )] = age*(x[index2(C2R,a-1,max_age)]-x[index2(C2R,a,max_age)]) -(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)])*x[index2(C2R,a,max_age)] + stoc1*FOI_R[a]*x[index2(RA,a,max_age)];
	
	dxdt[ index2(I2A , a, max_age )] = age*(x[index2(I2A,a-1,max_age)]-x[index2(I2A,a,max_age)]) +ir_A[a]*x[index2(C2A,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I2A,a,max_age)];
	
	dxdt[ index2(I2R , a, max_age )] = age*(x[index2(I2R,a-1,max_age)]-x[index2(I2R,a,max_age)]) +ir_R[a]*x[index2(C2R,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I2R,a,max_age)];
	
	dxdt[ index2(R, a, max_age )] = age*(x[index2(R,a-1,max_age)]-x[index2(R,a,max_age)]) +alpha_A*x[index2(C2A,a,max_age)] + alpha_R*x[index2(C2R,a,max_age)] + rho*x[index2(I2A,a,max_age)] + rho*x[index2(I2R,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)])*x[index2(R,a,max_age)]; 
	 
	 
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , a, max_age )] = ir_A[a]*x[index2(C1A,a,max_age)] + ir_A[a]*x[index2(C2A,a,max_age)];
	
	dxdt[ index2(IncR , a, max_age )] = ir_R[a]*x[index2(C1R,a,max_age)] + ir_R[a]*x[index2(C2R,a,max_age)];
	    
        
    }
	
    
    // Final Age Group
    dxdt[ index2(S , max_age-1, max_age )] = age*(x[index2(S,(max_age-2),max_age)])+phi*x[index2(RA,max_age-1,max_age)] + phi*x[index2(RR,max_age-1,max_age)] + phi*x[index2(R,max_age-1,max_age)] - stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)]
			 - stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)] - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(S,max_age-1,max_age)];
    
    dxdt[ index2(C1A , max_age-1, max_age )] = age*(x[index2(C1A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)];
	
	dxdt[ index2(C1R , max_age-1, max_age )] = age*(x[index2(C1R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)];
    
    dxdt[ index2(I1A , max_age-1, max_age )] = age*(x[index2(I1A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1A,max_age-1,max_age)];
	
	dxdt[ index2(I1R , max_age-1, max_age )] = age*(x[index2(I1R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1R,max_age-1,max_age)];
    
    dxdt[ index2(RA , max_age-1, max_age )] = age*(x[index2(RA,(max_age-2),max_age)])+rho*x[index2(I1A,max_age-1,max_age)]+alpha_A*x[index2(C1A,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RA,max_age-1,max_age)]
	 - stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)];
    
	dxdt[ index2(RR , max_age-1, max_age )] = age*(x[index2(RR,(max_age-2),max_age)])+rho*x[index2(I1R,max_age-1,max_age)]+alpha_R*x[index2(C1R,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RR,max_age-1,max_age)]
	 - stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)];
	
	dxdt[ index2(C2A , max_age-1, max_age )] = age*(x[index2(C2A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)];
	
	dxdt[ index2(C2R , max_age-1, max_age )] = age*(x[index2(C2R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)];
	
	dxdt[ index2(I2A , max_age-1, max_age )] = age*(x[index2(I2A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2A,max_age-1,max_age)];
	
	dxdt[ index2(I2R , max_age-1, max_age )] = age*(x[index2(I2R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2R,max_age-1,max_age)];
	
	dxdt[ index2(R, max_age-1, max_age )] = age*(x[index2(R,(max_age-2),max_age)])+alpha_A*x[index2(C2A,max_age-1,max_age)] + alpha_R*x[index2(C2R,max_age-1,max_age)] + rho*x[index2(I2A,max_age-1,max_age)] + rho*x[index2(I2R,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(R,max_age-1,max_age)]; 
	 
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , max_age-1, max_age )] = ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)] + ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)];

	dxdt[ index2(IncR , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)];
        
    
    return List::create(dxdt, stoc1,B, N);
}



// [[Rcpp::export]]



List PentaLongV(double t, NumericVector x, NumericVector params, NumericMatrix Cm, NumericVector mu, double max_age, NumericVector B, NumericVector ir_A, NumericVector ir_R, NumericVector stoc, NumericVector epi_A, NumericVector epi_R)
{
    
    NumericVector dxdt(x.length());
    
    double phi = params[0]; //loss of immunity
    double rho = params[1]; // recovery of disease
    double alpha_A = params[2]; // duration of carriage menA
	double alpha_R = params[3]; // duration of carriage menCWYX
    double epsbeta = params[4]; // seasonal amplitude
    double beta_A = params[5]; // contact parameter
	double beta_R = params[6];
	double ksi = params[7];
	double delta = params[8];
	double wl = params[9];
    
    double pi = M_PI;
    
    double age = 1.0/365.0; // Daily Ageing
    
    int comp=30;
    
    enum  state_variables {S,C1A,C1R,I1A,I1R,RA,RR,C2A,C2R,I2A,I2R,R,SVlong,C1AVlong,C1RVlong,I1AVlong,I1RVlong,RAVlong,RRVlong,C2AVlong,C2RVlong,I2AVlong,I2RVlong,RVlong,IncA,IncAV,IncA_all,IncR,IncRV,IncR_all};
    
    NumericVector FOI_A(max_age);
	NumericVector FOI_R(max_age);
    
    double stoc1 = stoc[t/365];
    
    double birth = B[t/365];
    
    double N = 0.0;
    for( int j=0; j<(comp-6)*max_age; j++)
    {
        N += x[j];
    }
    
    // Calculate force of infection
    for(int a=0 ; a<max_age ; a++ )
    {
        FOI_A[a] = 0.0;
		FOI_R[a] = 0.0;
		
        for( int a2=0; a2<max_age; a2++ )
        {
            
            double Z = (1+epsbeta*cos(2*pi*t/365.0));
			
            FOI_A[a] += beta_A*Z*Cm(a,a2)*(x[index2(C1A,a2,max_age)] + x[index2(C2A,a2,max_age)] + x[index2(I1A,a2,max_age)] + x[index2(I2A,a2,max_age)] 
			+ x[index2(C1AVlong,a2,max_age)] + x[index2(C2AVlong,a2,max_age)] + x[index2(I1AVlong,a2,max_age)] + x[index2(I2AVlong,a2,max_age)])/N;
			
			FOI_R[a] += beta_R*Z*Cm(a,a2)*(x[index2(C1R,a2,max_age)] + x[index2(C2R,a2,max_age)] + x[index2(I1R,a2,max_age)] + x[index2(I2R,a2,max_age)] 
			+ x[index2(C1RVlong,a2,max_age)] + x[index2(C2RVlong,a2,max_age)] + x[index2(I1RVlong,a2,max_age)] + x[index2(I2RVlong,a2,max_age)])/N;
        }
    }
    
    
    
    dxdt[ index2(S , 0, max_age )] = birth*N + phi*x[index2(RA,0,max_age)] + phi*x[index2(RR,0,max_age)] + phi*x[index2(R,0,max_age)] - stoc1*FOI_A[0]*x[index2(S,0,max_age)] - stoc1*FOI_R[0]*x[index2(S,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(S,0,max_age)] + wl*x[index2(SVlong,0,max_age)];
    
    dxdt[ index2(C1A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1A,0,max_age)] + stoc1*FOI_A[0]*x[index2(S,0,max_age)]
	+ wl*x[index2(C1AVlong,0,max_age)];
	
	dxdt[ index2(C1R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1R,0,max_age)] + stoc1*FOI_R[0]*x[index2(S,0,max_age)]
	+ wl*x[index2(C1RVlong,0,max_age)];
    
    dxdt[ index2(I1A , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1A,0,max_age)];
	
	dxdt[ index2(I1R , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1R,0,max_age)];
    
    dxdt[ index2(RA , 0, max_age )] = rho*x[index2(I1A,0,max_age)]+alpha_A*x[index2(C1A,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RA,0,max_age)] - stoc1*FOI_R[0]*x[index2(RA,0,max_age)]
	+ wl*x[index2(RAVlong,0,max_age)];
    
	dxdt[ index2(RR , 0, max_age )] = rho*x[index2(I1R,0,max_age)]+alpha_R*x[index2(C1R,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RR,0,max_age)] - stoc1*FOI_A[0]*x[index2(RR,0,max_age)]
	+ wl*x[index2(RRVlong,0,max_age)];
	
	dxdt[ index2(C2A , 0, max_age )] = -(ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2A,0,max_age)] + stoc1*FOI_A[0]*x[index2(RR,0,max_age)]
	+ wl*x[index2(C2AVlong,0,max_age)];
	
	dxdt[ index2(C2R , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2R,0,max_age)] + stoc1*FOI_R[0]*x[index2(RA,0,max_age)]
	+ wl*x[index2(C2RVlong,0,max_age)];
	
	dxdt[ index2(I2A , 0, max_age )] = ir_A[0]*x[index2(C2A,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2A,0,max_age)];
	
	dxdt[ index2(I2R , 0, max_age )] = ir_R[0]*x[index2(C2R,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2R,0,max_age)];
	
	dxdt[ index2(R, 0, max_age )] = alpha_A*x[index2(C2A,0,max_age)] + alpha_R*x[index2(C2R,0,max_age)] + rho*x[index2(I2A,0,max_age)] + rho*x[index2(I2R,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(R,0,max_age)] + wl*x[index2(RVlong,0,max_age)]; 
	
	// MenAfriVac
	 
	dxdt[ index2(SVlong , 0, max_age )] = phi*x[index2(RAVlong,0,max_age)] + phi*x[index2(RRVlong,0,max_age)] + phi*x[index2(RVlong,0,max_age)] - stoc1*(1-delta)*FOI_A[0]*x[index2(SVlong,0,max_age)] - stoc1*FOI_R[0]*x[index2(SVlong,0,max_age)]
	 - (age+mu[index2(floor(t/365),0,max_age)])*x[index2(SVlong,0,max_age)] - wl*x[index2(SVlong,0,max_age)];
    
    dxdt[ index2(C1AVlong , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1AVlong,0,max_age)] 
	+ stoc1*(1-delta)*FOI_A[0]*x[index2(SVlong,0,max_age)] - wl*x[index2(C1AVlong,0,max_age)];
	
	dxdt[ index2(C1RVlong , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C1RVlong,0,max_age)] + stoc1*FOI_R[0]*x[index2(SVlong,0,max_age)]
	- wl*x[index2(C1RVlong,0,max_age)];
    
    dxdt[ index2(I1AVlong , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1AVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1AVlong,0,max_age)];
	
	dxdt[ index2(I1RVlong , 0, max_age )] = ir_R[0]*x[index2(C1RVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I1RVlong,0,max_age)];
    
    dxdt[ index2(RAVlong , 0, max_age )] = rho*x[index2(I1AVlong,0,max_age)]+alpha_A*x[index2(C1AVlong,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RAVlong,0,max_age)] 
	- stoc1*FOI_R[0]*x[index2(RAVlong,0,max_age)] - wl*x[index2(RAVlong,0,max_age)];
    
	dxdt[ index2(RRVlong , 0, max_age )] = rho*x[index2(I1RVlong,0,max_age)]+alpha_R*x[index2(C1RVlong,0,max_age)]-(phi+age+mu[index2(floor(t/365),0,max_age)])*x[index2(RRVlong,0,max_age)]
	 - stoc1*(1-delta)*FOI_A[0]*x[index2(RRVlong,0,max_age)] - wl*x[index2(RRVlong,0,max_age)];
	
	dxdt[ index2(C2AVlong , 0, max_age )] = -((1-ksi)*ir_A[0]+alpha_A+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2AVlong,0,max_age)] + stoc1*(1-delta)*FOI_A[0]*x[index2(RRVlong,0,max_age)]
	- wl*x[index2(C2AVlong,0,max_age)];
	
	dxdt[ index2(C2RVlong , 0, max_age )] = -(ir_R[0]+alpha_R+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(C2RVlong,0,max_age)] + stoc1*FOI_R[0]*x[index2(RAVlong,0,max_age)]
	- wl*x[index2(C2RVlong,0,max_age)];
	
	dxdt[ index2(I2AVlong , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C2AVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2AVlong,0,max_age)];
	
	dxdt[ index2(I2RVlong , 0, max_age )] = ir_R[0]*x[index2(C2RVlong,0,max_age)]-(rho+age+mu[index2(floor(t/365),0,max_age)])*x[index2(I2RVlong,0,max_age)];
	
	dxdt[ index2(RVlong, 0, max_age )] = alpha_A*x[index2(C2AVlong,0,max_age)] + alpha_R*x[index2(C2RVlong,0,max_age)] + rho*x[index2(I2AVlong,0,max_age)] + rho*x[index2(I2RVlong,0,max_age)]
	 - (phi+mu[index2(floor(t/365),0,max_age)]+age)*x[index2(RVlong,0,max_age)] - wl*x[index2(RVlong,0,max_age)]; 
   

	 
   
    
   
   // Incidence
    
    dxdt[ index2(IncA , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)] + ir_A[0]*x[index2(C2A,0,max_age)];
	dxdt[ index2(IncAV , 0, max_age )] = (1-ksi)*ir_A[0]*x[index2(C1AVlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2AVlong,0,max_age)];

	dxdt[ index2(IncA_all , 0, max_age )] = ir_A[0]*x[index2(C1A,0,max_age)] + ir_A[0]*x[index2(C2A,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C1AVlong,0,max_age)] + (1-ksi)*ir_A[0]*x[index2(C2AVlong,0,max_age)];
	
	dxdt[ index2(IncR , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)] + ir_R[0]*x[index2(C2R,0,max_age)];
	dxdt[ index2(IncRV , 0, max_age )] = ir_R[0]*x[index2(C1RVlong,0,max_age)] + ir_R[0]*x[index2(C2RVlong,0,max_age)];

	dxdt[ index2(IncR_all , 0, max_age )] = ir_R[0]*x[index2(C1R,0,max_age)] + ir_R[0]*x[index2(C2R,0,max_age)] + ir_R[0]*x[index2(C1RVlong,0,max_age)] + ir_R[0]*x[index2(C2RVlong,0,max_age)];
	
    
	
	//Second age group
	
	    
    dxdt[ index2(S , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(S,0,max_age)]+phi*x[index2(RA,1,max_age)] + phi*x[index2(RR,1,max_age)] + phi*x[index2(R,1,max_age)] - stoc1*FOI_A[1]*x[index2(S,1,max_age)] - stoc1*FOI_R[1]*x[index2(S,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(S,1,max_age)] + wl*x[index2(SVlong,1,max_age)];
    
    dxdt[ index2(C1A , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C1A,0,max_age)]-(ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1A,1,max_age)] + stoc1*FOI_A[1]*x[index2(S,1,max_age)]
	+ wl*x[index2(C1AVlong,1,max_age)];
	
	dxdt[ index2(C1R , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C1R,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1R,1,max_age)] + stoc1*FOI_R[1]*x[index2(S,1,max_age)]
	+ wl*x[index2(C1RVlong,1,max_age)];
    
    dxdt[ index2(I1A , 1, max_age )] = age*x[index2(I1A,0,max_age)]+ir_A[1]*x[index2(C1A,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1A,1,max_age)];
	
	dxdt[ index2(I1R , 1, max_age )] = age*x[index2(I1R,0,max_age)]+ir_R[1]*x[index2(C1R,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1R,1,max_age)];
    
    dxdt[ index2(RA , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(RA,0,max_age)]+rho*x[index2(I1A,1,max_age)]+alpha_A*x[index2(C1A,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RA,1,max_age)] - stoc1*FOI_R[1]*x[index2(RA,1,max_age)]
	+ wl*x[index2(RAVlong,1,max_age)];
	
	dxdt[ index2(RR , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(RR,0,max_age)]+rho*x[index2(I1R,1,max_age)]+alpha_R*x[index2(C1R,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RR,1,max_age)] - stoc1*FOI_A[1]*x[index2(RR,1,max_age)]
	+ wl*x[index2(RRVlong,1,max_age)];
	
	dxdt[ index2(C2A , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C2A,0,max_age)]-(ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2A,1,max_age)] + stoc1*FOI_A[1]*x[index2(RR,1,max_age)]
	+ wl*x[index2(C2AVlong,1,max_age)];
	
	dxdt[ index2(C2R , 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(C2R,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2R,1,max_age)] + stoc1*FOI_R[1]*x[index2(RA,1,max_age)]
	+ wl*x[index2(C2RVlong,1,max_age)];
	
	dxdt[ index2(I2A , 1, max_age )] = age*x[index2(I2A,0,max_age)]+ir_A[1]*x[index2(C2A,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2A,1,max_age)];
	
	dxdt[ index2(I2R , 1, max_age )] = age*x[index2(I2R,0,max_age)]+ir_R[1]*x[index2(C2R,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2R,1,max_age)];
	
	dxdt[ index2(R, 1, max_age )] = (1-epi_A[t/365]-epi_R[t/365])*age*x[index2(R,0,max_age)]+alpha_A*x[index2(C2A,1,max_age)] + alpha_R*x[index2(C2R,1,max_age)] + rho*x[index2(I2A,1,max_age)] + rho*x[index2(I2R,1,max_age)]
	 - (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(R,1,max_age)] + wl*x[index2(RVlong,1,max_age)]; 
	
	// MenAfriVac
	 
	dxdt[ index2(SVlong , 1, max_age )] = age*x[index2(SVlong,0,max_age)]+phi*x[index2(RAVlong,1,max_age)] + phi*x[index2(RRVlong,1,max_age)] + phi*x[index2(RVlong,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(SVlong,1,max_age)] - stoc1*FOI_R[1]*x[index2(SVlong,1,max_age)]
	 - (age+mu[index2(floor(t/365),1,max_age)])*x[index2(SVlong,1,max_age)] - wl*x[index2(SVlong,1,max_age)];
    
    dxdt[ index2(C1AVlong , 1, max_age )] = age*x[index2(C1AVlong,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1AVlong,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(SVlong,1,max_age)]
	-wl*x[index2(C1AVlong,1,max_age)];
	
	dxdt[ index2(C1RVlong , 1, max_age )] = age*x[index2(C1RVlong,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C1RVlong,1,max_age)] + stoc1*FOI_R[1]*x[index2(SVlong,1,max_age)]
	-wl*x[index2(C1RVlong,1,max_age)];
    
    dxdt[ index2(I1AVlong , 1, max_age )] = age*x[index2(I1AVlong,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C1AVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1AVlong,1,max_age)];
	
	dxdt[ index2(I1RVlong , 1, max_age )] = age*x[index2(I1RVlong,0,max_age)]+ir_R[1]*x[index2(C1RVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I1RVlong,1,max_age)];
    
    dxdt[ index2(RAVlong , 1, max_age )] = age*x[index2(RAVlong,0,max_age)]+rho*x[index2(I1AVlong,1,max_age)]+alpha_A*x[index2(C1AVlong,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RAVlong,1,max_age)] - stoc1*FOI_R[1]*x[index2(RAVlong,1,max_age)]
	-wl*x[index2(RAVlong,1,max_age)];
    
	dxdt[ index2(RRVlong , 1, max_age )] = age*x[index2(RRVlong,0,max_age)]+rho*x[index2(I1RVlong,1,max_age)]+alpha_R*x[index2(C1RVlong,1,max_age)]-(phi+age+mu[index2(floor(t/365),1,max_age)])*x[index2(RRVlong,1,max_age)] - stoc1*(1-delta)*FOI_A[1]*x[index2(RRVlong,1,max_age)]
	-wl*x[index2(RRVlong,1,max_age)];
	
	dxdt[ index2(C2AVlong , 1, max_age )] = age*x[index2(C2AVlong,0,max_age)]-((1-ksi)*ir_A[1]+alpha_A+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2AVlong,1,max_age)] + stoc1*(1-delta)*FOI_A[1]*x[index2(RRVlong,1,max_age)]
	-wl*x[index2(C2AVlong,1,max_age)];
	
	dxdt[ index2(C2RVlong , 1, max_age )] = age*x[index2(C2RVlong,0,max_age)]-(ir_R[1]+alpha_R+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(C2RVlong,1,max_age)] + stoc1*FOI_R[1]*x[index2(RAVlong,1,max_age)]
	-wl*x[index2(C2RVlong,1,max_age)];
	
	dxdt[ index2(I2AVlong , 1, max_age )] = age*x[index2(I2AVlong,0,max_age)]+(1-ksi)*ir_A[1]*x[index2(C2AVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2AVlong,1,max_age)];
	
	dxdt[ index2(I2RVlong , 1, max_age )] = age*x[index2(I2RVlong,0,max_age)]+ir_R[1]*x[index2(C2RVlong,1,max_age)]-(rho+age+mu[index2(floor(t/365),1,max_age)])*x[index2(I2RVlong,1,max_age)];
	
	dxdt[ index2(RVlong, 1, max_age )] = age*x[index2(RVlong,0,max_age)]+alpha_A*x[index2(C2AVlong,1,max_age)] + alpha_R*x[index2(C2RVlong,1,max_age)] + rho*x[index2(I2AVlong,1,max_age)] + rho*x[index2(I2RVlong,1,max_age)]
	 	- (phi+mu[index2(floor(t/365),1,max_age)]+age)*x[index2(RVlong,1,max_age)] - wl*x[index2(RVlong,1,max_age)]; 
    
	 
	 
    
	 // Incidence
	   
	 dxdt[ index2(IncA , 1, max_age )] = ir_A[1]*x[index2(C1A,1,max_age)] + ir_A[1]*x[index2(C2A,1,max_age)];
	
	dxdt[ index2(IncAV , 1, max_age )] = (1-ksi)*ir_A[1]*x[index2(C1AVlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2AVlong,1,max_age)];
		dxdt[ index2(IncA_all , 1, max_age )] = ir_A[1]*x[index2(C1A,1,max_age)] + ir_A[1]*x[index2(C2A,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C1AVlong,1,max_age)] + (1-ksi)*ir_A[1]*x[index2(C2AVlong,1,max_age)];
	
	dxdt[ index2(IncR , 1, max_age )] = ir_R[1]*x[index2(C1R,1,max_age)] + ir_R[1]*x[index2(C2R,1,max_age)];
	dxdt[ index2(IncRV , 1, max_age )] = ir_R[1]*x[index2(C1RVlong,1,max_age)] + ir_R[1]*x[index2(C2RVlong,1,max_age)];
	
	dxdt[ index2(IncR_all , 1, max_age )] = ir_R[1]*x[index2(C1R,1,max_age)] + ir_R[1]*x[index2(C2R,1,max_age)] + ir_R[1]*x[index2(C1RVlong,1,max_age)] + ir_R[1]*x[index2(C2RVlong,1,max_age)];
	
	 
    
    // Annual age groups
    for( size_t a=2 ; a<(max_age-1) ; ++a )
    {
        
            dxdt[ index2(S , a, max_age )] = age*(x[index2(S,a-1,max_age)]-x[index2(S,a,max_age)]) + phi*x[index2(RA,a,max_age)] + phi*x[index2(RR,a,max_age)] + phi*x[index2(R,a,max_age)] - stoc1*FOI_A[a]*x[index2(S,a,max_age)]
			 - stoc1*FOI_R[a]*x[index2(S,a,max_age)] - (mu[index2(floor(t/365),a,max_age)])*x[index2(S,a,max_age)] + wl*x[index2(SVlong,a,max_age)];
    
    dxdt[ index2(C1A , a, max_age )] = age*(x[index2(C1A,a-1,max_age)]-x[index2(C1A,a,max_age)]) -(ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)])*x[index2(C1A,a,max_age)] + stoc1*FOI_A[a]*x[index2(S,a,max_age)]
	+ wl*x[index2(C1AVlong,a,max_age)];
	
	dxdt[ index2(C1R , a, max_age )] = age*(x[index2(C1R,a-1,max_age)]-x[index2(C1R,a,max_age)]) -(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)])*x[index2(C1R,a,max_age)] + stoc1*FOI_R[a]*x[index2(S,a,max_age)]
	+ wl*x[index2(C1RVlong,a,max_age)];
    
    dxdt[ index2(I1A , a, max_age )] = age*(x[index2(I1A,a-1,max_age)]-x[index2(I1A,a,max_age)]) +ir_A[a]*x[index2(C1A,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I1A,a,max_age)];
	
	dxdt[ index2(I1R , a, max_age )] = age*(x[index2(I1R,a-1,max_age)]-x[index2(I1R,a,max_age)]) +ir_R[a]*x[index2(C1R,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I1R,a,max_age)];
    
    dxdt[ index2(RA , a, max_age )] = age*(x[index2(RA,a-1,max_age)]-x[index2(RA,a,max_age)]) +rho*x[index2(I1A,a,max_age)]+alpha_A*x[index2(C1A,a,max_age)]-(phi+mu[index2(floor(t/365),a,max_age)])*x[index2(RA,a,max_age)]
	 - stoc1*FOI_R[a]*x[index2(RA,a,max_age)] + wl*x[index2(RAVlong,a,max_age)];
    
	dxdt[ index2(RR , a, max_age )] = age*(x[index2(RR,a-1,max_age)]-x[index2(RR,a,max_age)]) +rho*x[index2(I1R,a,max_age)]+alpha_R*x[index2(C1R,a,max_age)]-(phi+mu[index2(floor(t/365),a,max_age)])*x[index2(RR,a,max_age)]
	 - stoc1*FOI_A[a]*x[index2(RR,a,max_age)] + wl*x[index2(RRVlong,a,max_age)];
	
	dxdt[ index2(C2A , a, max_age )] = age*(x[index2(C2A,a-1,max_age)]-x[index2(C2A,a,max_age)]) -(ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)])*x[index2(C2A,a,max_age)] + stoc1*FOI_A[a]*x[index2(RR,a,max_age)]
	+ wl*x[index2(C2AVlong,a,max_age)];
	
	dxdt[ index2(C2R , a, max_age )] = age*(x[index2(C2R,a-1,max_age)]-x[index2(C2R,a,max_age)]) -(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)])*x[index2(C2R,a,max_age)] + stoc1*FOI_R[a]*x[index2(RA,a,max_age)]
	+ wl*x[index2(C2RVlong,a,max_age)];
	
	dxdt[ index2(I2A , a, max_age )] = age*(x[index2(I2A,a-1,max_age)]-x[index2(I2A,a,max_age)]) +ir_A[a]*x[index2(C2A,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I2A,a,max_age)];
	
	dxdt[ index2(I2R , a, max_age )] = age*(x[index2(I2R,a-1,max_age)]-x[index2(I2R,a,max_age)]) +ir_R[a]*x[index2(C2R,a,max_age)]-(rho+mu[index2(floor(t/365),a,max_age)])*x[index2(I2R,a,max_age)];
	
	dxdt[ index2(R, a, max_age )] = age*(x[index2(R,a-1,max_age)]-x[index2(R,a,max_age)]) +alpha_A*x[index2(C2A,a,max_age)] + alpha_R*x[index2(C2R,a,max_age)] + rho*x[index2(I2A,a,max_age)] + rho*x[index2(I2R,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)])*x[index2(R,a,max_age)] + wl*x[index2(RVlong,a,max_age)]; 
	 
	 
	 	// MenAfriVac
	 
	dxdt[ index2(SVlong , a, max_age )] = age*x[index2(SVlong,a-1,max_age)]+phi*x[index2(RAVlong,a,max_age)] + phi*x[index2(RRVlong,a,max_age)] + phi*x[index2(RVlong,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(SVlong,a,max_age)] - stoc1*FOI_R[a]*x[index2(SVlong,a,max_age)]
	 - (age+mu[index2(floor(t/365),a,max_age)])*x[index2(SVlong,a,max_age)] - wl*x[index2(SVlong,a,max_age)];
    
    dxdt[ index2(C1AVlong , a, max_age )] = age*x[index2(C1AVlong,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1AVlong,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(SVlong,a,max_age)]
	- wl*x[index2(C1AVlong,a,max_age)];
	
	dxdt[ index2(C1RVlong , a, max_age )] = age*x[index2(C1RVlong,a-1,max_age)]-(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C1RVlong,a,max_age)] + stoc1*FOI_R[a]*x[index2(SVlong,a,max_age)]
	- wl*x[index2(C1RVlong,a,max_age)];
    
    dxdt[ index2(I1AVlong , a, max_age )] = age*x[index2(I1AVlong,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C1AVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1AVlong,a,max_age)];
	
	dxdt[ index2(I1RVlong , a, max_age )] = age*x[index2(I1RVlong,a-1,max_age)]+ir_R[a]*x[index2(C1RVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I1RVlong,a,max_age)];
    
    dxdt[ index2(RAVlong , a, max_age )] = age*x[index2(RAVlong,a-1,max_age)]+rho*x[index2(I1AVlong,a,max_age)]+alpha_A*x[index2(C1AVlong,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RAVlong,a,max_age)] - stoc1*FOI_R[a]*x[index2(RAVlong,a,max_age)]
	- wl*x[index2(RAVlong,a,max_age)];
    
	dxdt[ index2(RRVlong , a, max_age )] = age*x[index2(RRVlong,a-1,max_age)]+rho*x[index2(I1RVlong,a,max_age)]+alpha_R*x[index2(C1RVlong,a,max_age)]-(phi+age+mu[index2(floor(t/365),a,max_age)])*x[index2(RRVlong,a,max_age)] - stoc1*(1-delta)*FOI_A[a]*x[index2(RRVlong,a,max_age)]
	- wl*x[index2(RRVlong,a,max_age)];
	
	dxdt[ index2(C2AVlong , a, max_age )] = age*x[index2(C2AVlong,a-1,max_age)]-((1-ksi)*ir_A[a]+alpha_A+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2AVlong,a,max_age)] + stoc1*(1-delta)*FOI_A[a]*x[index2(RRVlong,a,max_age)]
	- wl*x[index2(C2AVlong,a,max_age)];
	
	dxdt[ index2(C2RVlong , a, max_age )] = age*x[index2(C2RVlong,a-1,max_age)]-(ir_R[a]+alpha_R+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(C2RVlong,a,max_age)] + stoc1*FOI_R[a]*x[index2(RAVlong,a,max_age)]
	- wl*x[index2(C2RVlong,a,max_age)];
	
	dxdt[ index2(I2AVlong , a, max_age )] = age*x[index2(I2AVlong,a-1,max_age)]+(1-ksi)*ir_A[a]*x[index2(C2AVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2AVlong,a,max_age)];
	
	dxdt[ index2(I2RVlong , a, max_age )] = age*x[index2(I2RVlong,a-1,max_age)]+ir_R[a]*x[index2(C2RVlong,a,max_age)]-(rho+age+mu[index2(floor(t/365),a,max_age)])*x[index2(I2RVlong,a,max_age)];
	
	dxdt[ index2(RVlong, a, max_age )] = age*x[index2(RVlong,a-1,max_age)]+alpha_A*x[index2(C2AVlong,a,max_age)] + alpha_R*x[index2(C2RVlong,a,max_age)] + rho*x[index2(I2AVlong,a,max_age)] + rho*x[index2(I2RVlong,a,max_age)]
	 - (phi+mu[index2(floor(t/365),a,max_age)]+age)*x[index2(RVlong,a,max_age)] - wl*x[index2(RVlong,a,max_age)]; 
   
   
    
	 
	 
    
	 
    
  	 // Incidence
	   
	 dxdt[ index2(IncA , a, max_age )] = ir_A[a]*x[index2(C1A,a,max_age)] + ir_A[a]*x[index2(C2A,a,max_age)];
	dxdt[ index2(IncAV , a, max_age )] = (1-ksi)*ir_A[a]*x[index2(C1AVlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2AVlong,a,max_age)];
	
	dxdt[ index2(IncA_all , a, max_age )] = ir_A[a]*x[index2(C1A,a,max_age)] + ir_A[a]*x[index2(C2A,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C1AVlong,a,max_age)] + (1-ksi)*ir_A[a]*x[index2(C2AVlong,a,max_age)];
	
	dxdt[ index2(IncR , a, max_age )] = ir_R[a]*x[index2(C1R,a,max_age)] + ir_R[a]*x[index2(C2R,a,max_age)];
	dxdt[ index2(IncRV , a, max_age )] = ir_R[a]*x[index2(C1RVlong,a,max_age)] + ir_R[a]*x[index2(C2RVlong,a,max_age)];
	
	dxdt[ index2(IncR_all , a, max_age )] = ir_R[a]*x[index2(C1R,a,max_age)] + ir_R[a]*x[index2(C2R,a,max_age)] + ir_R[a]*x[index2(C1RVlong,a,max_age)] + ir_R[a]*x[index2(C2RVlong,a,max_age)];
        
        
        
        
    }
	
    
    // Final Age Group
    dxdt[ index2(S , max_age-1, max_age )] = age*(x[index2(S,(max_age-2),max_age)])+phi*x[index2(RA,max_age-1,max_age)] + phi*x[index2(RR,max_age-1,max_age)] + phi*x[index2(R,max_age-1,max_age)] - stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)]
			 - stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)] - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(S,max_age-1,max_age)]
			 + wl*x[index2(SVlong,max_age-1,max_age)];
    
    dxdt[ index2(C1A , max_age-1, max_age )] = age*(x[index2(C1A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(S,max_age-1,max_age)]
	+ wl*x[index2(C1AVlong,max_age-1,max_age)];
	
	dxdt[ index2(C1R , max_age-1, max_age )] = age*(x[index2(C1R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(S,max_age-1,max_age)]
	+ wl*x[index2(C1RVlong,max_age-1,max_age)];
    
    dxdt[ index2(I1A , max_age-1, max_age )] = age*(x[index2(I1A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1A,max_age-1,max_age)];
	
	dxdt[ index2(I1R , max_age-1, max_age )] = age*(x[index2(I1R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1R,max_age-1,max_age)];
    
    dxdt[ index2(RA , max_age-1, max_age )] = age*(x[index2(RA,(max_age-2),max_age)])+rho*x[index2(I1A,max_age-1,max_age)]+alpha_A*x[index2(C1A,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RA,max_age-1,max_age)]
	 - stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)] + wl*x[index2(RAVlong,max_age-1,max_age)];
    
	dxdt[ index2(RR , max_age-1, max_age )] = age*(x[index2(RR,(max_age-2),max_age)])+rho*x[index2(I1R,max_age-1,max_age)]+alpha_R*x[index2(C1R,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RR,max_age-1,max_age)]
	 - stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)] + wl*x[index2(RAVlong,max_age-1,max_age)];
	
	dxdt[ index2(C2A , max_age-1, max_age )] = age*(x[index2(C2A,(max_age-2),max_age)])-(ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2A,max_age-1,max_age)] + stoc1*FOI_A[max_age-1]*x[index2(RR,max_age-1,max_age)]
	+ wl*x[index2(C2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(C2R , max_age-1, max_age )] = age*(x[index2(C2R,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2R,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RA,max_age-1,max_age)]
	+ wl*x[index2(C2RVlong,max_age-1,max_age)];
	
	dxdt[ index2(I2A , max_age-1, max_age )] = age*(x[index2(I2A,(max_age-2),max_age)])+ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2A,max_age-1,max_age)];
	
	dxdt[ index2(I2R , max_age-1, max_age )] = age*(x[index2(I2R,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2R,max_age-1,max_age)];
	
	dxdt[ index2(R, max_age-1, max_age )] = age*(x[index2(R,(max_age-2),max_age)])+alpha_A*x[index2(C2A,max_age-1,max_age)] + alpha_R*x[index2(C2R,max_age-1,max_age)] + rho*x[index2(I2A,max_age-1,max_age)] + rho*x[index2(I2R,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(R,max_age-1,max_age)] + wl*x[index2(RVlong,max_age-1,max_age)]; 
	 
	 	 	// MenAfriVac
	 
	dxdt[ index2(SVlong , max_age-1, max_age )] = age*(x[index2(SVlong,(max_age-2),max_age)])+phi*x[index2(RAVlong,max_age-1,max_age)] + phi*x[index2(RRVlong,max_age-1,max_age)] + phi*x[index2(RVlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SVlong,max_age-1,max_age)] - stoc1*FOI_R[max_age-1]*x[index2(SVlong,max_age-1,max_age)]
	 - (mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(SVlong,max_age-1,max_age)] - wl*x[index2(SVlong,max_age-1,max_age)];
    
    dxdt[ index2(C1AVlong , max_age-1, max_age )] = age*(x[index2(C1AVlong,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1AVlong,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(SVlong,max_age-1,max_age)]
	- wl*x[index2(C1AVlong,max_age-1,max_age)];
	
	dxdt[ index2(C1RVlong , max_age-1, max_age )] = age*(x[index2(C1RVlong,(max_age-2),max_age)])+-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C1RVlong,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(SVlong,max_age-1,max_age)]
	- wl*x[index2(C1RVlong,max_age-1,max_age)];
    
    dxdt[ index2(I1AVlong , max_age-1, max_age )] = age*(x[index2(I1AVlong,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C1AVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1AVlong,max_age-1,max_age)];
	
	dxdt[ index2(I1RVlong , max_age-1, max_age )] = age*(x[index2(I1RVlong,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C1RVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I1RVlong,max_age-1,max_age)];
    
    dxdt[ index2(RAVlong , max_age-1, max_age )] = age*(x[index2(RAVlong,(max_age-2),max_age)])+rho*x[index2(I1AVlong,max_age-1,max_age)]+alpha_A*x[index2(C1AVlong,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RAVlong,max_age-1,max_age)] - stoc1*FOI_R[max_age-1]*x[index2(RAVlong,max_age-1,max_age)]
	- wl*x[index2(RAVlong,max_age-1,max_age)];
    
	dxdt[ index2(RRVlong , max_age-1, max_age )] = age*(x[index2(RRVlong,(max_age-2),max_age)])+rho*x[index2(I1RVlong,max_age-1,max_age)]+alpha_R*x[index2(C1RVlong,max_age-1,max_age)]-(phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RRVlong,max_age-1,max_age)] - stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRVlong,max_age-1,max_age)]
	- wl*x[index2(RRVlong,max_age-1,max_age)];
	
	dxdt[ index2(C2AVlong , max_age-1, max_age )] = age*(x[index2(C2AVlong,(max_age-2),max_age)])-((1-ksi)*ir_A[max_age-1]+alpha_A+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2AVlong,max_age-1,max_age)] + stoc1*(1-delta)*FOI_A[max_age-1]*x[index2(RRVlong,max_age-1,max_age)]
	- wl*x[index2(C2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(C2RVlong , max_age-1, max_age )] = age*(x[index2(C2RVlong,(max_age-2),max_age)])-(ir_R[max_age-1]+alpha_R+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(C2RVlong,max_age-1,max_age)] + stoc1*FOI_R[max_age-1]*x[index2(RAVlong,max_age-1,max_age)]
	- wl*x[index2(C2RVlong,max_age-1,max_age)];
	
	dxdt[ index2(I2AVlong , max_age-1, max_age )] = age*(x[index2(I2AVlong,(max_age-2),max_age)])+(1-ksi)*ir_A[max_age-1]*x[index2(C2AVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(I2RVlong , max_age-1, max_age )] = age*(x[index2(I2RVlong,(max_age-2),max_age)])+ir_R[max_age-1]*x[index2(C2RVlong,max_age-1,max_age)]-(rho+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(I2RVlong,max_age-1,max_age)];
	
	dxdt[ index2(RVlong, max_age-1, max_age )] = age*(x[index2(RVlong,(max_age-2),max_age)])+alpha_A*x[index2(C2AVlong,max_age-1,max_age)] + alpha_R*x[index2(C2RVlong,max_age-1,max_age)] + rho*x[index2(I2AVlong,max_age-1,max_age)] + rho*x[index2(I2RVlong,max_age-1,max_age)]
	 - (phi+mu[index2(floor(t/365),max_age-1,max_age)])*x[index2(RVlong,max_age-1,max_age)] - wl*x[index2(RVlong,max_age-1,max_age)]; 
   
    
	 
  	 // Incidence
	   
	 dxdt[ index2(IncA , max_age-1, max_age )] = ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)] + ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)];
	dxdt[ index2(IncAV , max_age-1, max_age )] = (1-ksi)*ir_A[max_age-1]*x[index2(C1AVlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(IncA_all , max_age-1, max_age )] = ir_A[max_age-1]*x[index2(C1A,max_age-1,max_age)] + ir_A[max_age-1]*x[index2(C2A,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C1AVlong,max_age-1,max_age)] + (1-ksi)*ir_A[max_age-1]*x[index2(C2AVlong,max_age-1,max_age)];
	
	dxdt[ index2(IncR , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)];
	dxdt[ index2(IncRV , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1RVlong,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2RVlong,max_age-1,max_age)];
	
	dxdt[ index2(IncR_all , max_age-1, max_age )] = ir_R[max_age-1]*x[index2(C1R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2R,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C1RVlong,max_age-1,max_age)] + ir_R[max_age-1]*x[index2(C2RVlong,max_age-1,max_age)];
        
    
    return List::create(dxdt, stoc1,birth, N);
}
