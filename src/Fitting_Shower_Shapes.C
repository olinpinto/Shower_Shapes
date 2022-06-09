//myIncludes
#include "../include/myFits.h"
#include "../include/Minimizer.h"
#include "../include/Global.h"

void Fitting_Shower_Shapes() {

  //===============================================
  // Fitting Shower shapes using the CALICE-AHCAL
  //      Sept 2020
  //   @author O. Pinto
  //===============================================


  bool isElectrons = false;
  bool isPions = true;

  //=================================================================================================================//
  //                                           Fitting Shower Shapes for Electrons                                   //
  //=================================================================================================================//

  if (isElectrons)
    {
	  //===================fitting longitudinal profile=================//
	  //-- parameters --
	  // [0] - Energy
	  // [1] - alpha_short
	  // [2] - beta_ short

	  longprof = new TF1("longprof",lngProfileFunction_Electrons, 0, 38, 3);
	  longprof->SetParNames ("E", "alpha", "beta");
	  longprof->SetParameters(h_energySum->GetMean(), 2.0, 1.2);

	  //===Radial fit for EM-Showers====//
	  TF1 *radhist = new TF1("radhist",radProfileFunction, 0, 300,4);
	  TF1 *radhist_CORE = new TF1("radhist_CORE",radProfileFunction_CORE, 0, 300,3);
	  TF1 *radhist_HALO = new TF1("radhist_HALO",radProfileFunction_HALO, 0, 300,3);

	  radhist->SetParNames ("E", "b_{c}","f_{EM}", "b_{h}");
	  radhist->SetParameters(h_energySum->GetMean(), 5.,0.8, 20.);

	  double x[NRINGS], y[NRINGS], ex[NRINGS], ey[NRINGS];
	  for (int i_rad=1;i_rad<=h_radial_energydensity_sys->GetNbinsX(); i_rad++)
	    {
	      x[i_rad]  = h_radial_energydensity_sys->GetBinCenter(i_rad);
	      y[i_rad]  = h_radial_energydensity_sys->GetBinContent(i_rad);
	      ex[i_rad] = h_radial_energydensity_sys->GetBinWidth(i_rad)/ sqrt(12.);
	      ey[i_rad] = h_radial_energydensity_sys->GetBinError(i_rad);
	    }

	  TGraphErrors *ge = new TGraphErrors(NRINGS,x,y,ex,ey);
	  ge->Draw("AP");

    TFitResultPtr rFit = ge->Fit("radhist", "EMRS");
    TMatrixDSym rcov = rFit->GetCovarianceMatrix();
  }

  //=================================================================================================================//
  //                                           Fitting Shower Shapes for Pions                                       //
  //=================================================================================================================//

  if (isPions)
    {
  	  h_longitudinal_avg =(TH1D*)inputFile->Get("h_longitudinal_avg");
  	  h_radial_energydensity_sys =(TH1D*)inputFile->Get("h_radial_energydensity_sys"); // Radial histogram
  	  ge =(TGraphErrors*)inputFile->Get("Graph"); // Radial Graph

  	  //===================fitting longitudinal profile=================//
  	  TF1 *longprof = new TF1("longprof",lngProfileFunction, 0, 4.0,6);
  	  TF1 *longprof_short = new TF1("longprof_short",lngProfileFunction_SHORT, 0, 4.0,4);
  	  TF1 *longprof_long = new TF1("longprof_long",lngProfileFunction_LONG, 0, 4.0,4);

  	  longprof->SetParNames ("E", "a_{s}", "b_{s}","a_{l}", "b_{l}","f_{EM}");
  	  longprof_short->SetParNames ("E", "a_{s}", "b_{s}","f_{EM}");
  	  longprof_long->SetParNames ("E", "a_{l}", "b_{l}","f_{EM}");
  	  longprof->SetParameters(h_energySum->GetMean(), 3.5, 1.5, 1.0, 1.0, 0.25);

  	  TFitResultPtr lFit = h_longitudinal_avg->Fit("longprof","QEMRS");
  	  TMatrixDSym lcov = lFit->GetCovarianceMatrix();
  	  TMatrixDSym lcor = lFit->GetCorrelationMatrix();
  	  int n_lpar = longprof->GetNpar();

  	  //=====SETTING PARAMETERS FOR THE TWO COMPONENTS OF LONGITUDINAL FUNCTION====//
  	  longprof_short->SetParameters(longprof->GetParameter(0), longprof->GetParameter(1), longprof->GetParameter(2), longprof->GetParameter(5));
  	  longprof_long->SetParameters(longprof->GetParameter(0), longprof->GetParameter(3), longprof->GetParameter(4), longprof->GetParameter(5));

  	  h_longitudinal_avg->Draw("E");
  	  longprof->Draw("SAME");
  	  longprof_short->SetLineColor(kBlue+3);
  	  longprof_short->Draw("SAME");
  	  longprof_long->SetLineColor(kGreen+3);
  	  longprof_long->Draw("SAME");

  	  //====================FITTING RADIAL PROFILES========================//

  	  //____________________PRELIMINARY METHOD_________________________________//
  	  //--- parameters ---
  	  // [0] - Total Energy
  	  // [1] - beta_core
  	  // [2] - fraction
  	  // [3] - beta_halo

  	  TF1 *radhist = new TF1("radhist",radProfileFunction, MIN_R, MAX_R,4);// the range of fit and #pars
  	  TF1 *radhist_CORE = new TF1("radhist_CORE",radProfileFunction_CORE, MIN_R, MAX_R,3);
  	  TF1 *radhist_HALO = new TF1("radhist_HALO",radProfileFunction_HALO, MIN_R, MAX_R,3);

  	  radhist->SetParNames ("E", "b_{c}","f_{EM}", "b_{h}");
  	  radhist->SetParameters(h_energySum->GetMean(), 12.,0.6, 65.);

  	  //============TGraphErrors for RADIAL profile===============//
  	  ge->Draw("AP");
  	  TFitResultPtr rFit = ge->Fit("radhist", "QEMRS");

  	  //____________________EFFECTIVE RADIAL FITTING_________________________________//
  	  TF1 *radhist_EFK = new TF1("radhist_EFK",EFK_FITRAD, MIN_R, MAX_R, 7);
  	  TF1 *radhist_EFKCORE = new TF1("radhist_EFKCORE",EFK_FITRAD_CORE, MIN_R, MAX_R,4);
  	  TF1 *radhist_EFKHADCORE = new TF1("radhist_EFKHADCORE",EFK_FITRAD_HADCORE, MIN_R, MAX_R,5);
  	  TF1 *radhist_EFKHADHALO = new TF1("radhist_EFKHADHALO",EFK_FITRAD_HADHALO, MIN_R, MAX_R,5);

  	  radhist_EFK->SetParNames ("E", "K", "f_{EM}","b_{c}^{EM}", "b_{c}^{HAD}", "b_{h}^{HAD}", "f_{h}");
  	  radhist_EFKCORE->SetParNames ("E", "K", "f_{EM}", "b_{c}^{EM}");
  	  radhist_EFKHADCORE->SetParNames ("E", "K", "f_{EM}", "b_{c}^{HAD}", "f_{h}");
  	  radhist_EFKHADHALO->SetParNames ("E", "K", "f_{EM}", "b_{h}^{HAD}", "f_{h}");

  	  double E_EM  = longprof->GetParameter(5)*longprof->GetParameter(0);// E_EM = f*E
  	  double E_HAD = (1-longprof->GetParameter(5))*longprof->GetParameter(0);// E_Had = (1-f)*E
  	  double K = longprof->GetParameter(0)/h_energySum->GetMean();

      radhist_EFK->SetParameters(h_energySum->GetMean(), K, 0.1, 10.0, radhist->GetParameter(1), radhist->GetParameter(3), 0.3);
      radhist_EFK->FixParameter(1, K);
      radhist_EFK->FixParameter(3, 11.30);// Fixed for Data
      radhist_EFK->FixParameter(3, 10.55);// Fixed for MC's

  	  ge->Fit("radhist_EFK", "QEMRNO");

  	  //=====SETTING PARAMETERS FOR THE TWO COMPONENTS OF RADIAL FUNCTION====//
  	  radhist_EFKCORE->SetParameters(radhist_EFK->GetParameter(0), radhist_EFK->GetParameter(1), radhist_EFK->GetParameter(2), radhist_EFK->GetParameter(3));
  	  radhist_EFKHADCORE->SetParameters(radhist_EFK->GetParameter(0), radhist_EFK->GetParameter(1), radhist_EFK->GetParameter(2), radhist_EFK->GetParameter(4), radhist_EFK->GetParameter(6));
  	  radhist_EFKHADHALO->SetParameters(radhist_EFK->GetParameter(0), radhist_EFK->GetParameter(1), radhist_EFK->GetParameter(2), radhist_EFK->GetParameter(5), radhist_EFK->GetParameter(6));

  	  radhist_EFK->Draw("SAME");
  	  radhist_EFKCORE->Draw("SAME");
  	  radhist_EFKCORE->SetLineColor(kBlue+3);
  	  radhist_EFKHADCORE->Draw("SAME");
  	  radhist_EFKHADCORE->SetLineColor(kMagenta+3);
  	  radhist_EFKHADHALO->Draw("SAME");
  	  radhist_EFKHADHALO->SetLineColor(kGreen+3);
  	  legend_radialEFK->Draw("SAME");

  	  TMatrixDSym rcov = rFit->GetCovarianceMatrix();
  	  TMatrixDSym rcor = rFit->GetCorrelationMatrix();
  	  int n_rpar = radhist->GetNpar();

  	  for (int i_rpar =0; i_rpar<n_rpar; ++i_rpar)
  	    for (int j_rpar =0; j_rpar<n_rpar; ++j_rpar)
  		     h_rcorr->Fill(i_rpar, j_rpar, rcor[i_rpar][j_rpar]);

  	  const char *xr_labels[4] = {"E", "b_{c}","f_{EM}", "b_{h}"};
  	  const char *yr_labels[4] = {"E", "b_{c}","f_{EM}", "b_{h}"};
  	  for (Int_t i = 0; i < n_rpar; i++)
  	    {
  	      h_rcorr->GetXaxis()->SetBinLabel(i + 1, xr_labels[i]);
  	      h_rcorr->GetYaxis()->SetBinLabel(i + 1, xr_labels[i]);
  	    }

  	  //=====SETTING PARAMETERS FOR THE TWO COMPONENTS OF RADIAL FUNCTION====//
  	  radhist_CORE->SetParameters(radhist->GetParameter(0), radhist->GetParameter(2), radhist->GetParameter(1));
  	  radhist_HALO->SetParameters(radhist->GetParameter(0), radhist->GetParameter(2), radhist->GetParameter(3));

  	  radhist_CORE->SetLineColor(kBlue+3);
  	  radhist_HALO->SetLineColor(kGreen+3);
  	  radhist->Draw("SAME");
  	  legend_radial->Draw("SAME");
  	  radhist_CORE->Draw("SAME");
  	  radhist_HALO->Draw("SAME");

	  //-------------------------------------------------------------------------------------------------------------------------------------------//
    //==========================================SIMULTANEOUS FITTING LONGITUDINAL & RADIAL PROFILE ==============================================//
	  //-------------------------------------------------------------------------------------------------------------------------------------------//

      vector <double> params, paramsErr;

      for (int p=0; p<=5; ++p)
    		{
    		  params.push_back(longprof->GetParameter(p));
    		  paramsErr.push_back(longprof->GetParError(p));
    		}
      for (int p=0; p<=6; ++p)
    		{
    		  params.push_back(radhist_EFK->GetParameter(p));
    		  paramsErr.push_back(radhist_EFK->GetParError(p));
    		}

      runFit(params, paramsErr);
      FitResults();
      params = fitter.Result().Parameters();
      paramsErr = fitter.Result().Errors();

      TF1 *myCHI2_global_long;
      TF1 *myCHI2_short_part;
      TF1 *myCHI2_long_part;
      TF1 *myCHI2_global_rad;
      TF1 *myCHI2_rad_core;
      TF1 *myCHI2_rad_hadcore;
      TF1 *myCHI2_rad_hadhalo;

		  myCHI2_global_long = new TF1("myCHI2_global_long",lngProfileFunction, 0, MAX_Z,N_LONGPARS);
		  myCHI2_global_long->SetParameters(params[0], params[1],params[2], params[3], params[4], params[5]);

		  myCHI2_short_part = new TF1("myCHI2_short_part",lngProfileFunction_SHORT, 0, MAX_Z,4);
		  myCHI2_short_part->SetParameters(params[0], params[1],params[2], params[5]);

		  myCHI2_long_part = new TF1("myCHI2_long_part",lngProfileFunction_LONG, 0, MAX_Z,4);
		  myCHI2_long_part->SetParameters(params[0], params[3], params[4], params[5]);

		  myCHI2_global_rad = new TF1("myCHI2_global_rad",EFK_FITRAD, MIN_R, MAX_R,N_RADPARS);
		  myCHI2_global_rad->SetParameters(params[0], params[6], params[5], params[7], params[8], params[9], params[10]);

		  myCHI2_rad_core = new TF1("myCHI2_rad_core",EFK_FITRAD_CORE, MIN_R, MAX_R,4);
		  myCHI2_rad_core->SetParameters(params[0], params[6], params[5], params[7]);

		  myCHI2_rad_hadcore = new TF1("myCHI2_rad_hadcore",EFK_FITRAD_HADCORE, MIN_R, MAX_R,5);
		  myCHI2_rad_hadcore->SetParameters(params[0], params[6], params[5], params[8], params[10]);

		  myCHI2_rad_hadhalo = new TF1("myCHI2_rad_hadhalo",EFK_FITRAD_HADHALO, MIN_R, MAX_R,5);
		  myCHI2_rad_hadhalo->SetParameters(params[0], params[6], params[5], params[9], params[10]);

      //******************************************************************
      //  Plot profiles
      //  longitudinal and radial
      //******************************************************************
      //LONGITUDINAL//
      h_longitudinal_avg->Draw("E");
      myCHI2_global_long->Draw("SAME");
      myCHI2_short_part->SetLineColor(kBlue+3);
      myCHI2_short_part->Draw("SAME");
      myCHI2_long_part->SetLineColor(kGreen+3);
      myCHI2_long_part->Draw("SAME");

      //RADIAL//
      h_radial_energydensity_sys->Draw("E");
      myCHI2_global_rad->Draw("SAME");
      myCHI2_rad_core->SetLineColor(kBlue+3);
      myCHI2_rad_core->Draw("SAME");
      myCHI2_rad_hadcore->SetLineColor(kMagenta+3);
      myCHI2_rad_hadcore->Draw("SAME");
      myCHI2_rad_hadhalo->SetLineColor(kGreen+3);
      myCHI2_rad_hadhalo->Draw("SAME");


	  //=================================================================================================================//
	  //                                           LayerWise-Fitting Shower Shapes for Pions                             //
	  //=================================================================================================================//

  	  TF1 *radhistL2L_EFK = new TF1("radhistL2L_EFK",EFK_FITRAD_L2L, MIN_R, MAX_R, 6);
  	  TF1 *radhistL2L_EFKCORE = new TF1("radhistL2L_EFKCORE",EFK_FITRAD_CORE_L2L, MIN_R, MAX_R,3);
  	  TF1 *radhistL2L_EFKHADCORE = new TF1("radhistL2L_EFKHADCORE",EFK_FITRAD_HADCORE_L2L, MIN_R, MAX_R,4);
  	  TF1 *radhistL2L_EFKHADHALO = new TF1("radhistL2L_EFKHADHALO",EFK_FITRAD_HADHALO_L2L, MIN_R, MAX_R,4);

  	  radhistL2L_EFK->SetParNames ("E", "f_{EM}","b_{c}^{EM}", "b_{c}^{HAD}", "b_{h}^{HAD}", "f_{h}");
  	  radhistL2L_EFKCORE->SetParNames ("E", "f_{EM}", "b_{c}^{EM}");
  	  radhistL2L_EFKHADCORE->SetParNames ("E", "f_{EM}", "b_{c}^{HAD}", "f_{h}");
  	  radhistL2L_EFKHADHALO->SetParNames ("E", "f_{EM}", "b_{h}^{HAD}", "f_{h}");

  	  radhistL2L_EFK->SetParameters(h_longitudinal_avg->GetBinContent(i_hist+1), radhist_EFK->GetParameter(2)*0.2, radhist_EFK->GetParameter(3), radhist_EFK->GetParameter(4), radhist_EFK->GetParameter(5), radhist_EFK->GetParameter(6)*0.5);
  	  TFitResultPtr radfit_L2L = gr_layer2layer[i_hist]->Fit("radhistL2L_EFK", "QEMRS");
  	  radhistL2L_EFKCORE->SetParameters(radhistL2L_EFK->GetParameter(0), radhistL2L_EFK->GetParameter(1), radhistL2L_EFK->GetParameter(2));
  	  radhistL2L_EFKHADCORE->SetParameters(radhistL2L_EFK->GetParameter(0), radhistL2L_EFK->GetParameter(1), radhistL2L_EFK->GetParameter(3), radhistL2L_EFK->GetParameter(5));
  	  radhistL2L_EFKHADHALO->SetParameters(radhistL2L_EFK->GetParameter(0), radhistL2L_EFK->GetParameter(1), radhistL2L_EFK->GetParameter(4), radhistL2L_EFK->GetParameter(5));

      //=================================================================================================================//
      //                                           Shower Shapes 3D fitting                                              //
      //=================================================================================================================//

      cout << "===========Performing 3D fitting ===========" << endl;

      TF2 *fit3D = new TF2("fit3D",_2Dfitting, 0, 4.0, -300, 300, 12);

      double E_EM = gr_radial_L2L[0][7][i_file]->GetFunction("f1")->GetParameter(0);
      double A_S = gr_radial_L2L[0][7][i_file]->GetFunction("f1")->GetParameter(1);
      double B_S = gr_radial_L2L[0][7][i_file]->GetFunction("f1")->GetParameter(2);

      double E_HAD = gr_radial_L2L[0][8][i_file]->GetFunction("f1")->GetParameter(0);
      double A_L = gr_radial_L2L[0][8][i_file]->GetFunction("f1")->GetParameter(1);
      double B_L = gr_radial_L2L[0][8][i_file]->GetFunction("f1")->GetParameter(2);

      double N = gr_radial_L2L[0][5][i_file]->GetFunction("f1")->GetParameter(0);
      double A_T = gr_radial_L2L[0][5][i_file]->GetFunction("f1")->GetParameter(1);
      double B_T = gr_radial_L2L[0][5][i_file]->GetFunction("f1")->GetParameter(2);

      const std::vector<std::pair<double, std::string> > parameters = {
        {E_EM, "E_{EM}"},
        {A_S, "a_{s}"},
        {B_S, "b_{s}"},
        {11.0, "b_{c}^{EM}"},
        {E_HAD,  "E_{HAD}"},
        {A_L,  "a_{l}"},
        {B_L,  "b_{l}"},
        {N,  "N"},
        {A_T,  "a_{t}"},
        {B_T,  "b_{t}"},
        {Beta_Core,  "b_{c}^{HAD}"}, //value from EM showers
        {Beta_Halo,  "b_{h}^{HAD}"}// value from EM Showers
      };
      for (size_t i =0; i< parameters.size(); ++i) {
        fit3D->SetParameter(i, parameters[i].first);
        fit3D->SetParName(i, parameters[i].second.c_str());
      }
      fit3D->FixParameter(10, Beta_Core);
      fit3D->FixParameter(11, Beta_Halo);


      for (int i=0; i<3; ++i)
        {
          gr_3Dshapes->Draw("lego2Z");//("lego2Z");//with option Z it will draw the color scale
          fit3D->FixParameter(3, 11.30);// Data
          fit3D->FixParameter(3, 10.55);// MC's
          gr_3Dshapes->Fit("fit3D", "QEMRS");
          TFitResultPtr FitResult3D = gr_3Dshapes->Fit("fit3D", "EMRS");
          fit3D->Draw("lego2same");
          gr_3Dshapes->SetTitle("; Layer (#lambda_{I}) ; Radius (mm); Energy (MIP)");

          TMatrixDSym FitResult3Dcov = FitResult3D->GetCovarianceMatrix();
          gr_3Dshapes->Draw("COLZ");
          fit3D->Draw("same");

          TMatrixDSym rfcor = FitResult3D->GetCorrelationMatrix();
          int nf_par = fit3D->GetNpar();
          h_fcorr = new TH2D("", "", nf_par, 0, nf_par, nf_par, 0, nf_par);
          h_fcorr->SetMarkerSize(1.5);
          h_fcorr->GetZaxis()->SetTitle("Correlation");

          for (int i_rpar =0; i_rpar<nf_par; ++i_rpar)
            {
              h_fcorr->GetXaxis()->SetBinLabel(i_rpar + 1, parameters[i_rpar].second.c_str());
              h_fcorr->GetYaxis()->SetBinLabel(i_rpar + 1, parameters[i_rpar].second.c_str());
              for (int j_rpar =0; j_rpar<nf_par; ++j_rpar)
                  h_fcorr->Fill(i_rpar, j_rpar, rfcor[i_rpar][j_rpar]);
            }
        }

        // ==============Plotting the residual for 3D modelling================//
      TString nameres; nameres.Form("Relativedeviation%d",i_file);
      c_3DRes.emplace_back(new TCanvas(nameres, nameres, 2000,1500));
      c_3DRes[i_file]->SetRightMargin(0.24);
      c_3DRes[i_file]->SetTheta(45.);
      c_3DRes[i_file]->SetPhi(60.);
      c_3DRes[i_file]->SetLogz(0);

        for (int i=0; i<3; ++i)
          {
            const auto* xVals = gr_3Dshapes->GetX();
            const auto* yVals = gr_3Dshapes->GetY();
            const auto* zVals = gr_3Dshapes->GetZ();
            const auto* EzVals = gr_3Dshapes->GetEZ();

            TGraph2D *gr2_residual = new TGraph2D();
            gr2_residual->SetTitle("; Layer (#lambda_{I}) ; Radius (mm); Relative deviation");

            for (int i = 0; i < gr_3Dshapes[0][0]->GetN(); ++i) {
              const Double_t fitVal = fit3D->Eval(xVals[i], yVals[i]);
              const Double_t residual = (zVals[i] - fitVal)/fitVal;
              gr2_residual->SetPoint(i,xVals[i],yVals[i],residual);
            }

            gr2_residual->Draw("lego2Z");

            gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
            TPaletteAxis *palette2 = (TPaletteAxis*)gr2_residual->GetHistogram()->GetListOfFunctions()->FindObject("palette");
            if (palette2) {
              palette2->SetY1NDC(0.45);
              gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn
            }
          }
        }
