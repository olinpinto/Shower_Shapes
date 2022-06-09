#ifndef LIBRARY_minimizer_H
#define LIBRARY_minimizer_H

#include "Fit/Fitter.h"
#include "Fit/ParameterSettings.h"
#include "myFits.h"
#include <vector>


using namespace std;
using ParamsSettings = std::vector<ROOT::Fit::ParameterSettings>;

//====Global variables===//
/* Vectors for CHI2 */
vector <double> Chi2overNDF_LONG, Chi2overNDF_RAD, Chi2overNDF_COMBINED;

ROOT::Fit::Fitter fitter;

TH2D *h_cor, *h_cov;

/* Function to fix the parameter value */
void fixParams(int Parindex, ParamsSettings& parameters)
{
  auto &param = parameters[Parindex];
  param.Fix();
}


void runFit(const std::vector<double>& params, const std::vector<double>& paramsErr)
{
  ParamsSettings startParams;


  // EM fraction is the same for long and radial. Calculate the starting value as average of the long and radial fractions from the preliminary fit
  double longFrac = 0.5 * (params[5] + params[8]);
  double longFracErr = 0.5 * (paramsErr[5] + paramsErr[8]);

  double Energy = 0.5 * (params[0] + params[6]);
  double EnergyErr = 0.5 * (paramsErr[0] + paramsErr[6]);


  startParams.push_back(ROOT::Fit::ParameterSettings("Energy", Energy, EnergyErr));  // 0
  startParams.push_back(ROOT::Fit::ParameterSettings("Alpha_short", params[1], paramsErr[1]));  // 1 // 0, 6.5
  startParams.push_back(ROOT::Fit::ParameterSettings("Beta_short", params[2], paramsErr[2]));  // 2
  startParams.push_back(ROOT::Fit::ParameterSettings("Alpha_long", params[3], paramsErr[3]));  // 3
  startParams.push_back(ROOT::Fit::ParameterSettings("Beta_long", params[4], paramsErr[4]));  // 4
  startParams.push_back(ROOT::Fit::ParameterSettings("EM_fraction", longFrac, longFracErr));   // 5
  startParams.push_back(ROOT::Fit::ParameterSettings("K_value", params[7], paramsErr[7]));  // 6
  startParams.push_back(ROOT::Fit::ParameterSettings("B_EM-core", params[9], paramsErr[9])); // 7
  startParams.push_back(ROOT::Fit::ParameterSettings("B_HAD-core", params[10], paramsErr[10])); // 8
  startParams.push_back(ROOT::Fit::ParameterSettings("B_HAD-halo", params[11], paramsErr[11], 50 , 100)); // 9 //, 60, 100
  startParams.push_back(ROOT::Fit::ParameterSettings("R_Halo_fraction", params[12], paramsErr[12])); // 10


  /* Parameter to be fixed */
  fixParams(6, startParams);
  fixParams(7, startParams);

  //for (int p =0; p<startParams.size(); ++p)
    //cout << startParams[p].IsFixed() << " : "  << startParams[p].Value() <<endl;

  /* set the start parameters from where the fit should start */
  fitter.Config().SetParamsSettings(startParams);

  /* Fitter Options */
  ROOT::Math::MinimizerOptions minOpt;
  minOpt.SetMinimizerType("Minuit2");
  //minOpt.SetPrintLevel(1);
  fitter.Config().SetMinimizerOptions(minOpt);

  // Set the data options
  ROOT::Fit::DataOptions options;
  options.fCoordErrors = false;// NO EFFECT!! use errors on the x coordinates when available (default is true)

  /* fit function, needs to know how many paramters and your fcn */
  fitter.FitFCN(startParams.size(), fcn);

}

/* Function storing fit results */
void FitResults()
{
  const auto& fitResult = fitter.Result();

  const auto chi2 = fitResult.Chi2();

  const auto NDF = fitResult.Ndf();

  const auto pvalue = fitResult.Prob();

  const auto MinFcn = fitResult.MinFcnValue();

  const auto nFreeParams = fitResult.NFreeParameters();

  const auto nTotalParams = fitResult.NTotalParameters();

  const auto& par = fitResult.Parameters();

  const auto& err = fitResult.Errors();

  fitter.GetMinimizer()->PrintResults();
  //fitResult.PrintCovMatrix(cout);

  h_cor = new TH2D("h_correlation", "", nTotalParams, 0, nTotalParams, nTotalParams, 0, nTotalParams);
  h_cov = new TH2D("h_covariance", "", nTotalParams, 0, nTotalParams, nTotalParams, 0, nTotalParams);

  const char *x_llabels[11] = {"E", "#alpha_{s}", "#beta_{s}","#alpha_{l}", "#beta_{l}","f_{EM}", "K","#beta_{c}^{EM}","#beta_{c}^{HAD}", "#beta_{h}^{HAD}", "f_{h}"};
  const char *y_llabels[11] = {"E", "alpha_{s}", "#beta_{s}","#alpha_{l}", "#beta_{l}","f_{EM}", "K","#beta_{c}^{EM}","#beta_{c}^{HAD}", "#beta_{h}^{HAD}", "f_{h}"};


  for (int i_par =0; i_par<nTotalParams; ++i_par)
    for (int j_par =0; j_par<nTotalParams; ++j_par)
      {
      	h_cor->Fill(i_par, j_par, fitResult.Correlation(i_par,j_par));
      	h_cov->Fill(i_par, j_par, fitResult.CovMatrix(i_par,j_par));
      }

  for (Int_t i = 0; i < nTotalParams; i++)
    {
      h_cor->GetXaxis()->SetBinLabel(i + 1, x_llabels[i]);
      h_cor->GetYaxis()->SetBinLabel(i + 1, x_llabels[i]);
      h_cov->GetXaxis()->SetBinLabel(i + 1, x_llabels[i]);
      h_cov->GetYaxis()->SetBinLabel(i + 1, x_llabels[i]);
    }

  h_cor->Write();
  h_cov->Write();

  /* Storing CHI2NDF values */
  Chi2overNDF_LONG.push_back(CHI2NDF_L);
  Chi2overNDF_RAD.push_back(CHI2NDF_R);
  Chi2overNDF_COMBINED.push_back(CHI2NDF_TOTAL);
}




/* Needed for 2D contour plots */
/* 1-Sigma = 68.27% */
/* 2-Sigma = 95.45% */
/* 3-Sigma = 99.73% */
TGraph get2DContour(double confLevel=0.683, unsigned nPoints=50)
{
  auto *minimizer = fitter.GetMinimizer();
  if (!minimizer) {
    cout << "Could not get minimizer from fitter\n";
  }

  cout << "Recalculating the errors for the desired confLevel = " << confLevel << "\n";
  const double oldErrDef = minimizer->ErrorDef();
  minimizer->SetErrorDef(oldErrDef * TMath::ChisquareQuantile(confLevel, 2));
  const bool hesse = fitter.CalculateHessErrors();
  minimizer->PrintResults();
  const bool minos = fitter.CalculateMinosErrors();
  minimizer->PrintResults();
  if (!hesse || !minos) {
    std::cout << "HESSE (" << hesse << ") or MINOS (" << minos << ") errors could not be determined\n";
    return TGraph();
  }

  std::vector<double> x_c(50, 0.0);
  std::vector<double> y_c(50, 0.0);

  if (!minimizer->Contour(0, 1, nPoints, x_c.data(), y_c.data())) {
    std::cout << "Could not get contour\n";
    return TGraph();
  }

  minimizer->SetErrorDef(oldErrDef); // reset

  return TGraph(nPoints, x_c.data(), y_c.data());
}


#endif
