#ifndef LIBRARY_GLOBAL_H
#define LIBRARY_GLOBAL_H

#include "../include/myFits.h"
#include <TGraphErrors.h>
#include <vector>

using namespace std;


//*******************************************************************************
// iterative fit of graph
//*******************************************************************************
TF1* fitProfile(TGraphErrors* gr,  int npar, int nIter,  double chi2ndfLowLimit)
{
  double chi2ndf = 100;
  double chi2ndfMin = chi2ndf;
  int iter = 0;

  TFitResultPtr fitres;
  char funName[20];
  sprintf(funName,"radFun");


  //---------- Random number generator -----------------------------------
  TRandom2 rg;

  TF1 *ffit = new TF1(funName, EFK_FITRAD, 0, 300, npar);
  ffit->SetParameters(100, 0.5, 0.2, 10, 30, 80, 0.4);
  //---- iterations to get best fit
  while ( chi2ndf > chi2ndfLowLimit && iter < nIter )
    {
      iter++;
      for ( int ip = 0; ip < npar; ip++ )
      	{
      	  /* ffit->SetParameter(ip, par[ip]); */
      	  /* ffit->SetParErrors(ip, parErr[ip]); */
      	  if (ip ==3 )
      	    ffit->FixParameter(ip,10.3);
      	  else
      	    ffit->ReleaseParameter(ip);
      	}

      fitres = gr->Fit(funName,"QRM0S");
      chi2ndf = ffit->GetChisquare()/ffit->GetNDF();

      double iterFrac = ffit->GetParameter(2);
      if ( iterFrac < 0 || iterFrac > 1 ) {
	       chi2ndf = 100;
      }

      else
    	{ // check if slopes > 0
    	  if ( ffit->GetParameter(4) < 0 || ffit->GetParameter(5) < 0 ) {
    	  chi2ndf = 100;
    	  }
	    }

      if ( chi2ndf < chi2ndfMin) {
	       chi2ndfMin = chi2ndf;
	       cout << chi2ndf << " : " << ffit->GetChisquare() << " : "   << ffit->GetParameter(0) << endl;
	       ffit->Draw("same");
      }
    }
  return ffit;
}

#endif
