#ifndef LIBRARY_myfits_H
#define LIBRARY_myfits_H
#include <TF1.h>
#include "TMath.h"

using namespace std;


/*CALICE-AHCAL thickeness
 Radiation Length,  X0   = 1.05;
 Interaction Length, LI   = 0.11;
*/


//******************************************************************
// Fit Function for longitudinal Profile of EM-Showers: 3 parameters
//******************************************************************
Double_t lngProfileFunction_Electrons(Double_t *x, Double_t *par) {
  // Gamma function (https://en.wikipedia.org/wiki/Incomplete_gamma_function) with parameters:
  // [0] - Energy
  // [1] - alpha
  // [2] - beta

  double z_X0     = x[0];
  double X0       = 1.05;
  double Energy   = par[0];
  double Alpha    = par[1];
  double Beta     = par[2];

  double shortComponent = (1/TMath::Gamma(Alpha)) * pow((z_X0/Beta),(Alpha - 1)) * (exp(-z_X0/Beta)/Beta); // In terms of X_0
  return Energy * X0 * shortComponent;
}



//**********************************************************************
// Fit Function for longitudinal Profile of Hadron Showers: 6 parameters
//**********************************************************************
Double_t lngProfileFunction(Double_t *x, Double_t *par) {

  //-- parameters --
  // [0] - Energy
  // [1] - Alpha short
  // [2] - Beta short
  // [3] - Alpha long
  // [4] - Beta long
  // [5] - fraction

  double z_LI = x[0];
  double X0   = 1.05;
  double LI   = 0.11;
  double z_X0 = (X0/LI)*z_LI;

  double Energy   = par[0];
  double Alpha_S  = par[1];
  double Beta_S   = par[2];
  double Alpha_L  = par[3];
  double Beta_L   = par[4];
  double Frac     = par[5];


  double shortComponent = (Frac/TMath::Gamma(Alpha_S)) * pow((z_X0/Beta_S),(Alpha_S - 1)) * (exp(-z_X0/Beta_S)/(Beta_S*LI/X0));
  double longComponent = ((1-Frac)/TMath::Gamma(Alpha_L)) * pow((z_LI/Beta_L),(Alpha_L - 1)) * (exp(-z_LI/Beta_L)/Beta_L);

  return Energy * LI * (shortComponent + longComponent);
}


//******************************************************************
// Fit Function for Longitudinal Long Component: 4 parameters
//******************************************************************
Double_t lngProfileFunction_SHORT(Double_t *x, Double_t *par) {

  //-- parameters --
  // [0] - scaling
  // [1] - alpha short
  // [2] - beta  short
  // [3] - fraction

  double z_LI = x[0];
  double X0   = 1.05;
  double LI   = 0.11;
  double z_X0 = (X0/LI)*z_LI;

  double Energy   = par[0];
  double Alpha_S  = par[1];
  double Beta_S   = par[2];
  double Frac     = par[3];

  double shortComponent = (Frac/TMath::Gamma(Alpha_S)) * pow((z_X0/Beta_S),(Alpha_S - 1)) * (exp(-z_X0/Beta_S)/(Beta_S*LI/X0));
  return Energy * LI * shortComponent;
}


//******************************************************************
//  Fit Function for Longitudinal Short Component: 4 parameters
//******************************************************************
Double_t lngProfileFunction_LONG(Double_t *x, Double_t *par) {

  //-- parameters --
  // [0] - scaling
  // [1] - alpha long
  // [2] - beta  long
  // [3] - fraction

  double z_LI     = x[0];
  double LI       = 0.11;
  double Energy   = par[0];
  double Alpha_L  = par[1];
  double Beta_L   = par[2];
  double Frac     = par[3];

  double longComponent = ((1-Frac)/TMath::Gamma(Alpha_L)) * pow((z_LI/Beta_L),(Alpha_L - 1)) * (exp(-z_LI/Beta_L)/Beta_L);

  return Energy * LI * longComponent;
}


//******************************************************************
// Fit Function for Radial Shower Shape: 4 parameters
//******************************************************************
Double_t radProfileFunction(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - beta_core
  // [2] - fraction
  // [3] - beta_halo

  double r          = x[0];
  double energy     = par[0];
  double beta_core  = par[1];
  double f          = par[2];
  double beta_halo  = par[3];

  double coreComponent = (f /pow(beta_core, 2) ) *  exp(-r/beta_core);
  double haloComponent = (1 - f)/pow(beta_halo, 2)  *  exp(-r/beta_halo);

  return energy/(2*M_PI) * (coreComponent + haloComponent);
}


//******************************************************************
//  Fit Function for Radial Core: 3 parameters
//******************************************************************
Double_t radProfileFunction_CORE(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - fraction
  // [2] - beta_core

  double r          = x[0];
  double energy     = par[0];
  double f          = par[1];
  double beta_core  = par[2];

  double coreComponent = (f /pow(beta_core, 2) ) *  exp(-r/beta_core) ;
  return energy/(2*M_PI) * coreComponent;
}

//******************************************************************
//  Fit Function for Radial Halo: 3 parameters
//******************************************************************
Double_t radProfileFunction_HALO(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - fraction
  // [2] - beta_halo

  double r          = x[0];
  double energy     = par[0];
  double f          = par[1];
  double beta_halo  = par[2];

  double haloComponent = (1 - f)/pow(beta_halo, 2)  *  exp(-r/beta_halo) ; // r/beta
  return energy/(2*M_PI)  * haloComponent;
}



//**************************************************************************
// Improvised Fit Function for Radial Shapes of Hadron Showers: 7 parameters
//**************************************************************************
Double_t EFK_FITRAD(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - K constant: E_inf/E_vis; (Leakage)
  // [2] - fraction
  // [3] - beta_core_EM
  // [4] - beta_core_HAD
  // [5] - beta_halo_HAD
  // [6] - f, f_h (core and halo fraction of Hadronic Shower)

  double r               = x[0];
  double E               = par[0];
  double K               = par[1];
  double f               = par[2];
  double beta_core_EM    = par[3];
  double beta_core_HAD   = par[4];
  double beta_halo_HAD   = par[5];
  double f_h             = par[6];

  double coreComponent =  f * K *  (exp(-r/beta_core_EM)/pow(beta_core_EM, 2));
  double haloComponent =  (1-f*K) * ( f_h * (exp(-r/beta_core_HAD)/pow(beta_core_HAD,2)) + (1-f_h) * (exp(-r/beta_halo_HAD)/pow(beta_halo_HAD,2)) );

  return E/(2*M_PI*K) * (coreComponent + haloComponent);

}


//******************************************************************
//Integral of Effective Radial function (x_lower bin, x_upperbin)
//******************************************************************

Double_t EFK_FITRADINT(Double_t *x, Double_t *par) {

  TF1 radfunc = TF1("radfunc",EFK_FITRAD, 0, 300, 7);
  radfunc.SetParameters(par);

  return radfunc.Integral(x[0], x[1])/(x[1] - x[0]);
}

//******************************************************************
//Derivative of Effective Radial function (x_lower bin, x_upperbin)
//******************************************************************

Double_t EFK_FITRADDERV(Double_t *x, Double_t *par) {

  TF1 radfunc = TF1("radfunc",EFK_FITRAD, 0, 300, 7);
  radfunc.SetParameters(par);

  return radfunc.Derivative(x[0], par, 0.01);
}



//******************************************************************
// Improvised Fit Function for CORE of Hadron Showers: 4 parameters
//******************************************************************
Double_t EFK_FITRAD_CORE(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - K constant: E_inf/E_vis
  // [2] - fraction
  // [3] - beta_core_EM

  double r               = x[0];
  double E               = par[0];
  double K               = par[1];
  double f               = par[2];
  double beta_core_EM    = par[3];

  double coreComponent =  f * K *  (exp(-r/beta_core_EM)/pow(beta_core_EM, 2));

  return E/(2*M_PI*K) * coreComponent;
}


//*****************************************************************************
// Improvised Fit Function for Layer-Wise EM-Core of Hadron Showers: 3 parameters
//*****************************************************************************
Double_t EFK_FITRAD_CORE_L2L(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - fraction
  // [2] - beta_core_EM

  double r               = x[0];
  double E               = par[0];
  double f               = par[1];
  double beta_core_EM    = par[2];

  double coreComponent =  f * (exp(-r/beta_core_EM)/pow(beta_core_EM, 2));

  return E/(2*M_PI) * coreComponent;
}


//***************************************************************************************
// Improvised Fit Function for Layer-Wise Hadronic Core of Hadron Showers: 5 parameters
//******************************************************************
Double_t EFK_FITRAD_HADCORE(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - K constant: E_inf/E_vis
  // [2] - fraction
  // [3] - beta_core_HAD
  // [4] - fraction, f_h

  double r               = x[0];
  double E               = par[0];
  double K               = par[1];
  double f               = par[2];
  double beta_core_HAD   = par[3];
  double f_h             = par[4];

  double hadCoreComponent =  (1-f*K) * ( f_h * (exp(-r/beta_core_HAD)/pow(beta_core_HAD,2)) );

  return E/(2*M_PI*K) * hadCoreComponent;
}


//*******************************************************************
// Improvised Fit Function for Layer-Wise Hadronic Halo: 5 parameters
//*******************************************************************
Double_t EFK_FITRAD_HADHALO(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - K constant: E_inf/E_vis
  // [2] - fraction
  // [3] - beta_halo_HAD
  // [4] - fraction, f_h

  double r               = x[0];
  double E               = par[0];
  double K               = par[1];
  double f               = par[2];
  double beta_halo_HAD   = par[3];
  double f_h             = par[4];

  double hadHaloComponent =  (1-f*K) * ((1-f_h) * (exp(-r/beta_halo_HAD)/pow(beta_halo_HAD,2)) );

  return E/(2*M_PI*K) * hadHaloComponent;
}


//******************************************************
// Improvised Fit Function for Layer-Wise: 6 parameters
//******************************************************
Double_t EFK_FITRAD_L2L(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - fraction
  // [2] - beta_core_EM
  // [3] - beta_core_HAD
  // [4] - beta_halo_HAD
  // [5] - fraction, f_h

  double r               = x[0];
  double E               = par[0];
  double f               = par[1];
  double beta_core_EM    = par[2];
  double beta_core_HAD   = par[3];
  double beta_halo_HAD   = par[4];
  double f_h             = par[5];

  double coreComponent =  f * (exp(-r/beta_core_EM)/pow(beta_core_EM, 2));
  double haloComponent =  (1-f) * ( f_h * (exp(-r/beta_core_HAD)/pow(beta_core_HAD,2)) + (1-f_h) * (exp(-r/beta_halo_HAD)/pow(beta_halo_HAD,2)) );

  return E/(2*M_PI) * (coreComponent + haloComponent);

}

//********************************************************************
// Improvised Fit Function for Layer-Wise Hadronic Core: 4 parameters
//********************************************************************
Double_t EFK_FITRAD_HADCORE_L2L(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - fraction
  // [2] - beta_core_HAD
  // [3] - fraction, f_h

  double r               = x[0];
  double E               = par[0];
  double f               = par[1];
  double beta_core_HAD   = par[2];
  double f_h             = par[3];

  double hadCoreComponent =  (1-f) * ( f_h * (exp(-r/beta_core_HAD)/pow(beta_core_HAD,2)) );

  return E/(2*M_PI) * hadCoreComponent;
}



//******************************************************************
// Improvised Fit Function for Layer-Wise Hadronic Halo: 4 parameter
//******************************************************************
Double_t EFK_FITRAD_HADHALO_L2L(Double_t *x, Double_t *par) {

  //--- parameters ---
  // [0] - Energy
  // [1] - fraction
  // [2] - beta_halo_HAD
  // [3] - fraction, f_h

  double r               = x[0];
  double E               = par[0];
  double f               = par[1];
  double beta_halo_HAD   = par[2];
  double f_h             = par[3];

  double hadHaloComponent =  (1-f) * ((1-f_h) * (exp(-r/beta_halo_HAD)/pow(beta_halo_HAD,2)) );

  return E/(2*M_PI) * hadHaloComponent;
}


//******************************************************************
// CHI2-Method for Simulatanous fits of Longitudinal and
// Radial Showers
// https://root.cern.ch/doc/master/classROOT_1_1Fit_1_1Fitter.html
// fcn returns f, the function to be minimized.
// par – number of currently variable parameters.
//******************************************************************

//GLOBAL VARIABLES BEAWARE!!!!!!!
TH1D *h_longitudinal_avg, *h_radial_energydensity_sys;
TGraphErrors *ge;

const int N_LONGPARS = 6;//Number of Longitudinal Pars
const int N_RADPARS  = 7;//Number of Radial Pars
const double MAX_Z = 4.0,  MIN_R = 0.0 ,  MAX_R = 300.0;//Fit Ranges
const int NRINGS = 39;//Number of Rings
double x[NRINGS], y[NRINGS], ex[NRINGS], ey[NRINGS];

double CHI2NDF_L =0, CHI2NDF_R =0, CHI2NDF_TOTAL =0;

double fcn(const double* par)
{
  // definition of parameters
  // par[0]  - Longitudinal_Energy
  // par[1]  - alpha_short
  // par[2]  - beta_short
  // par[3]  - alpha_long
  // par[4]  - beta_long
  // par[5]  - L_fraction

  // par[6]  - Radial_Energy
  // par[7]  - K constant: E_inf/E_vis
  // par[8]  - R_fraction
  // par[9]  - beta_core_EM
  // par[10] - beta_core_HAD
  // par[11] - beta_halo_HAD
  // par[12] - Halo_fraction

  double  f = 0.;

  const int PAR_LFREE = 6, PAR_RFREE = 5;
  const double NFITTED_BINS_LONG = h_longitudinal_avg->GetXaxis()->FindBin(4.0);
  const double NFITTED_BINS_RAD  = h_radial_energydensity_sys->GetXaxis()->FindBin(300.0) -1;

  int ipar_l[N_LONGPARS] = {0,1,2,3,4,5};
  int ipar_r[N_RADPARS]  = {0,6,5,7,8,9,10};//(Parameter 5 and 8 in common)

  double par_l[N_LONGPARS];
  for (int i = 0; i < N_LONGPARS; ++i)
      par_l[i] = par[ipar_l[i] ];
  double par_r[N_RADPARS];
  for (int i = 0; i < N_RADPARS; ++i)
      par_r[i] = par[ipar_r[i] ];


  //==============CHI2_LONGITUDINAL=================//
  double chisq_long = 0;

  for (int i_long=1; i_long<=h_longitudinal_avg->GetNbinsX(); i_long++)
    {
      double x_long       = h_longitudinal_avg->GetBinCenter(i_long);
      if(x_long > MAX_Z ) continue;

      double measure_long = h_longitudinal_avg->GetBinContent(i_long);
      double error_long   = h_longitudinal_avg->GetBinError(i_long);
      double func_long    = lngProfileFunction( &x_long, par_l);

      if (error_long<=0) continue;
      double delta_long   = (func_long - measure_long)/error_long;

      chisq_long += delta_long * delta_long;
    }

  //==============CHI2_RADIAL=================//
  double chisq_radial = 0;

  for (int i_rad=1;i_rad<=h_radial_energydensity_sys->GetNbinsX(); i_rad++)
    {
      double x_radial       = h_radial_energydensity_sys->GetBinCenter(i_rad);
      if( x_radial < MIN_R || x_radial > MAX_R ) continue;
      double measure_radial = h_radial_energydensity_sys->GetBinContent(i_rad);
      double error_radial   = h_radial_energydensity_sys->GetBinError(i_rad);
      double binwidthx      = h_radial_energydensity_sys->GetBinWidth(i_rad)/ sqrt(12.);

      double x_bin[2] = {h_radial_energydensity_sys->GetBinLowEdge(i_rad), h_radial_energydensity_sys->GetBinLowEdge(i_rad +1)};
      double func_radial    = EFK_FITRADINT(x_bin, par_r);

      if ( error_radial<=0) continue;
      double delta_radial   =  (func_radial - measure_radial)/(error_radial +  (binwidthx * EFK_FITRADDERV(x_bin, par_r)));

      chisq_radial += delta_radial * delta_radial;
    }

  //=============COMBINED CHI2==================//
  f = chisq_long + chisq_radial;

  /* NDF:  No. of fitted bins - No. of free parameters */
  CHI2NDF_L = chisq_long/(NFITTED_BINS_LONG - PAR_LFREE); //LONGITUDINAL: Total 6 Free parameters,
  CHI2NDF_R = chisq_radial/(NFITTED_BINS_RAD  - PAR_RFREE); //RADIAL: Total 5 Free parameters + 2 FIXED (from EM Showers),
  CHI2NDF_TOTAL = CHI2NDF_L + CHI2NDF_R;

  /* cout << "____________________________________________________________________________________________" <<endl; */
  /* cout << " CHI2_L: " << chisq_long << " |NDF_L: " << NFITTED_BINS_LONG-6 <<" |CHI2NDF_L: " << CHI2NDF_L */
  /*   << " |CHI2_R: " << chisq_radial << " |NDF_R: " << NFITTED_BINS_RAD-5 << " |CHI2NDF_R: " << CHI2NDF_R */
  /*   << " |CHI2_TOTAL: " << f  << " |CHI2NDF_TOTAL: "  << CHI2NDF_TOTAL<< endl; */
  /* cout << "____________________________________________________________________________________________" <<endl; */

  return f;
}


//******************************************************************
// Three-Dimensional Hadronic Shower Parametrization: 12 parameters
//******************************************************************
Double_t _2Dfitting(Double_t *x, Double_t *par) {

  // definition of parameters
  // par[0]  - EM-Electromagnetic
  // par[1]  - alpha_short
  // par[2]  - beta_short
  // par[3]  - beta_EM_core
  // par[4]  - EM-Hadronic
  // par[5]  - alpha_long
  // par[6]  - beta_long
  // par[7]  - N
  // par[8]  - Alpha f_h
  // par[9]  - Beta f_h
  // par[10] - beta_core
  // par[11] - beta_halo

  double z_LI = x[0];
  double r    = abs(x[1]);
  double X0   = 1.05;
  double LI   = 0.11;
  double z_X0 = (X0/LI)*z_LI;

  double E_EM       = par[0];
  double Alpha_S    = par[1];
  double Beta_S     = par[2];
  double Beta_C_EM  = par[3];
  double E_HAD      = par[4];
  double Alpha_L    = par[5];
  double Beta_L     = par[6];
  double N          = par[7];
  double Alpha_T    = par[8];
  double Beta_T     = par[9];
  double Beta_C_HAD = par[10];
  double Beta_H_HAD = par[11];

  double shortComponent = (E_EM/TMath::Gamma(Alpha_S)) * (pow((z_X0/Beta_S),(Alpha_S - 1))) * (exp(-z_X0/Beta_S)/(Beta_S*LI/X0)) * (exp(-r/Beta_C_EM)/pow(Beta_C_EM, 2));
  double f_H = (N/TMath::Gamma(Alpha_T)) * (pow((z_LI/Beta_T),(Alpha_T - 1))) * (exp(-z_LI/Beta_T)/Beta_T);
  double longComponent =  (E_HAD/TMath::Gamma(Alpha_L)) * (pow((z_LI/Beta_L),(Alpha_L - 1))) * (exp(-z_LI/Beta_L)/Beta_L) * (   f_H * (exp(-r/Beta_C_HAD)/pow(Beta_C_HAD,2)) + (1-f_H) * (exp(-r/Beta_H_HAD)/pow(Beta_H_HAD,2)) );

  return   1/(2*M_PI) * (shortComponent + longComponent);
}


#endif
