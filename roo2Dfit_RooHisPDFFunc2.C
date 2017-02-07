#include "roo2Dfit_RooHisPDF.h"


void roo2Dfit_RooHisPDFFunc2(TString nModel = "RttbCon"){

  //setTDRStyle();
  gSystem->Load("libRooFit") ;
  using namespace RooFit;
  using namespace RooStats;

  HistoFit InFile[18];
  float ratio[15];

  // Data
  InFile[data] = LoadSample("DataSingleLep");

  // MC: Signal
  InFile[ttbb]   = LoadSample("ttbar_LepJetsPowhegPythiattbb");
  InFile[ttb]    = LoadSample("ttbar_LepJetsPowhegPythiattbj");
  InFile[ttcc]   = LoadSample("ttbar_LepJetsPowhegPythiattcc");
  InFile[ttLF]   = LoadSample("ttbar_LepJetsPowhegPythiattLF"); // Includes ttc
  InFile[ttccLF] = LoadSample("ttbar_LepJetsPowhegPythiattccLF");

  // MC: Backgrounds
  InFile[Bkgtt]    = LoadSample("ttbar_PowhegPythiaBkgtt"); // ttbarBkg + tt
  InFile[BkgOther] = LoadSample("BkgOther"); 

  TString dirfigname_pdf;
  dirfigname_pdf = dirnameIn + "figures_" + fl + "/ttbbFITNuisance" + nModel + "/pdf/";
  // make a dir if it does not exist!!
  gSystem->mkdir(dirfigname_pdf,       kTRUE);
  TString Rfile =  dirnameIn + "figures_" + fl + "/ttbbFITNuisance" + nModel + "/FitResultImp.log";
  FILE* fResult = fopen(Rfile, "w");


  // ----------------------------------------------------------------------------------
  // Initial Parameters
  // ----------------------------------------------------------------------------------
  // -- Cross Sections from MC 
  RooRealVar *Vis_Xsecttjj   = new RooRealVar("Vis_Xsecttjj",  "ttjj cross section Vis Ph-Sp",   26.40,    10.0,    50.0);// mu
  RooRealVar *Vis_Xsecttbb   = new RooRealVar("Vis_Xsecttbb",  "ttbb cross section Vis Ph-Sp",    0.375,    0.001,   1.5);// mu
  RooRealVar *Full_Xsecttjj  = new RooRealVar("Full_Xsecttjj", "ttjj cross section Full Ph-Sp", 290.2,    180.0,   300.0);
  RooRealVar *Full_Xsecttbb  = new RooRealVar("Full_Xsecttbb", "ttbb cross section Full Ph-Sp",   3.9,      1.0,     6.0);
  // Ratios Xsecttbb/Xsecttjj    
  RooFormulaVar *Vis_Xsecttbb_Xsecttjj  = new RooFormulaVar("Vis_Xsecttbb_Xsecttjj",  "Xsecttbb/Xsecttjj Vis-PhSp",  
							    "Vis_Xsecttbb/Vis_Xsecttjj",   RooArgList(*Vis_Xsecttbb,*Vis_Xsecttjj));
  RooFormulaVar *Full_Xsecttbb_Xsecttjj = new RooFormulaVar("Full_Xsecttbb_Xsecttjj", "Xsecttbb/Xsecttjj Full-PhSp", 
							    "Full_Xsecttbb/Full_Xsecttjj", RooArgList(*Full_Xsecttbb,*Full_Xsecttjj));
  RooRealVar *Vis_C_Xsecttbb_Xsecttjj   = new RooRealVar("Vis_C_Xsecttbb_Xsecttjj",  "Xsecttbb/Xsecttjj Vis-PhSp",  0.375/26.40,  0.001, 0.5); // 0.0142
  RooRealVar *Full_C_Xsecttbb_Xsecttjj  = new RooRealVar("Full_C_Xsecttbb_Xsecttjj", "Xsecttbb/Xsecttjj Full-PhSp", 3.90/290.2,  0.0,  1.0); // 0.0134

  
  // -- Efficiency
  vEffttjj[muJets] = 0.05;  vEffttjj[eJets] = 0.05;  vEffttjj[LepJets] = 0.0561;
  vEffttbb[muJets] = 0.22;  vEffttbb[eJets] = 0.22;  vEffttbb[LepJets] = 0.2216;
  // -- Acceptance
  vAccttjj[muJets] = 0.25*0.34;  vAccttjj[eJets] = 0.25*0.34;  vAccttjj[LepJets] = 0.2413*0.34;
  vAccttbb[muJets] = 0.29*0.34;  vAccttbb[eJets] = 0.29*0.34;  vAccttbb[LepJets] = 0.2914*0.34;
  // -- Luminosity 
  vLumi = 36814.; // Full 2016 Dataset

  for (unsigned int ch=2; ch<3; ch++){
    // Efficiencies ttjj
    RooRealVar *Effttjj_nom = new RooRealVar("Effttjj_nom", "ttjj Nom. Efficiency", vEffttjj[ch]);
    // Efficiencies ttbb
    RooRealVar *Effttbb_nom = new RooRealVar("Effttbb_nom", "ttbb Nom. Efficiency", vEffttbb[ch]);
    
    // Aceptancies ttjj
    RooRealVar *Accttjj_nom = new RooRealVar("Accttjj_nom",  "ttjj Nom. Acceptancy", vAccttjj[ch]);
    // Aceptancies ttbb
    RooRealVar *Accttbb_nom = new RooRealVar("Accttbb_nom",  "ttbb Nom. Acceptancy", vAccttbb[ch]);
    
    // Luminosity
    RooRealVar *Lumi_nom  = new RooRealVar("Lumi_nom",  "Luminosity", vLumi);
    
    // Variables CSV2 and CSV3
    // + 0.1 to avoid empty bins
    RooRealVar *CSV2 = new RooRealVar("CSV2", "CSV for Jet 3", 0.0,1.0); 
    RooRealVar *CSV3 = new RooRealVar("CSV3", "CSV for Jet 4", 0.0,1.0);  

    float n_ttjj = InFile[ttbb].events[ch] + 
      InFile[ttcc].events[ch] + 
      InFile[ttb].events[ch]  + 
      InFile[ttLF].events[ch];    

    ratio[ttbb]    = InFile[ttbb].events[ch]/n_ttjj; 
    ratio[ttb]     = InFile[ttb].events[ch]/n_ttjj; 
    ratio[ttccLF]  = InFile[ttccLF].events[ch]/n_ttjj; 

    cout << "----- Ratios from MC -----" << endl;
    cout << "ttbb/ttjj   = " << ratio[ttbb]   << endl;
    cout << "ttbj/ttjj   = " << ratio[ttb]    << endl;
    cout << "ttccLF/ttjj = " << ratio[ttccLF] << endl;
    cout << "--------------------------" << endl;
    
   
    // Ratio ttbj/ttbb
    RooRealVar *CRatio_ttbjttbb = new RooRealVar("CRatio_ttbjttbb", "MC ratio ttbj/ttbb",  InFile[ttb].events[ch]/InFile[ttbb].events[ch],InFile[ttb].events[ch]/InFile[ttbb].events[ch],InFile[ttb].events[ch]/InFile[ttbb].events[ch]);
    // Fitted ratio
    RooRealVar *FRatio_ttbbttjj = new RooRealVar("FRatio_ttbbttjj", "FITTED ratio ttbb/ttjj",   ratio[ttbb], 0.0, 0.2); 

    // Signal: N events
    RooRealVar *n_ttjj_var   = new RooRealVar("n_ttjj_var",   "number of ttjj events",   n_ttjj);    
    RooRealVar *n_ttbb_var   = new RooRealVar("n_ttbb_var",   "number of ttbb events",   InFile[ttbb].events[ch]);    
    RooRealVar *n_ttbj_var   = new RooRealVar("n_ttbj_var",   "number of ttbj events",   InFile[ttb].events[ch]);    
    RooRealVar *n_ttccLF_var = new RooRealVar("n_ttccLF_var", "number of ttccLF events", InFile[ttccLF].events[ch]);    
    // Background
    RooRealVar *n_Bkgtt_var    = new RooRealVar("n_Bkgtt_var",    "number of tt background events",    InFile[Bkgtt].events[ch]);
    RooRealVar *n_BkgOther_var = new RooRealVar("n_BkgOther_var", "number of Other background events", InFile[BkgOther].events[ch]);

    // 2D Arguments
    RooArgList *arg_CSV = new RooArgList(*CSV2, *CSV3);

    // // ----------------------------------------------------------------------------------
    // // ----------------------    Systematic Uncertainties   -----------------------------
    // // ----------------------------------------------------------------------------------

    RooDataHist *DataHisSys[18][3][14][3]; // [sample][case][syst][var]
    RooHistPdf  *HisPdfSys[18][3][14][3];  // [sample][case][syst][var]
    RooHistFunc  *HisFuncSys[18][3][14][3];  // [sample][case][syst][var]
    RooArgList  *arg_PdfVar[18][3][3],*arg_FuncVar[18][3][3];     // [sample][case][var]

    RooRealVar *AlphaSys[14][3]; // [syst][case]   
    RooArgList *arg_AlphaSys[3]; // [case]

    PiecewiseInterpolation *PInterSyst[18][3]; // [sample][case]
    

    for(int isam=ttbb; isam<=BkgOther; isam++){
      for(int isys=0; isys<14; isys++){
	if(isys!=btagcfI) continue;
    	for(int ivar=0; ivar<3; ivar++){
	  
    	  //////////
    	  // CSV2 //
    	  //////////
    	  // histograms
    	  DataHisSys[isam][0][isys][ivar] = new RooDataHist(SamNam[isam] + "_his2" + SystNam[isys] + VarNam[ivar],   SamNam[isam] + " Histogram",
    							    *CSV2, InFile[isam].hsyst1D[isys][ivar][0][ch]);
    	  // PDF
    	  HisPdfSys[isam][0][isys][ivar]   = new RooHistPdf(SamNam[isam] + "_hisPdf2" + SystNam[isys] + VarNam[ivar],   SamNam[isam] + " PDF",
    							    RooArgSet(*CSV2), *DataHisSys[isam][0][isys][ivar]);
	  // Func
    	  HisFuncSys[isam][0][isys][ivar]   = new RooHistFunc(SamNam[isam] + "_hisFunc2" + SystNam[isys] + VarNam[ivar],   SamNam[isam] + " PDF",
    							    RooArgSet(*CSV2), *DataHisSys[isam][0][isys][ivar]);

	  std::cout<< "INFO--> Histogram and PDF for " << SamNam[isam] + "_" + SystNam[isys] + VarNam[ivar] << " have been created for 2" << std::endl;

    	} // for(ivar)
      } // for(isys)



      
      // ----------------------------------------------------------------------------------
      // ----------------------------------------------------------------------------------

      RooRealVar *alpha = new RooRealVar("alpha", "Alpha Syst ",-500.0, 500.0);
      
      //////////
      // CSV2 //
      //////////      
      PInterSyst[isam][0] = new PiecewiseInterpolation(SamNam[isam] + "_PInterSyst2", "Interpolation for "  + SamNam[isam],
      						       *HisPdfSys [isam][0][btagcfI][Nom],// Nom [sample][Scen][syst][var] 
      						       *HisPdfSys [isam][0][btagcfI][Down],  // Down[sample][Scen][var]  
      						       *HisPdfSys [isam][0][btagcfI][Up],    // Up  [sample][Scen][var]  
      						       *alpha);           // Nuisance[Scen]

      // PInterSyst[isam][0] = new PiecewiseInterpolation(SamNam[isam] + "_PInterSyst2", "Interpolation for "  + SamNam[isam],
      // 						       *HisFuncSys [isam][0][btagcfI][Nom],// Nom [sample][Scen][syst][var] 
      // 						       *HisFuncSys [isam][0][btagcfI][Down],  // Down[sample][Scen][var]  
      // 						       *HisFuncSys [isam][0][btagcfI][Up],    // Up  [sample][Scen][var]  
      // 						       *alpha);           // Nuisance[Scen]
      
      PInterSyst[isam][0]->setPositiveDefinite(kTRUE);
      
      
    }// for(isam)
    
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // RooFormulaVar *CRttbjttbbFR_ttbbttjj = new RooFormulaVar("CRttbjttbbFR_ttbbttjj", "ttbj term",
    // 							     "CRatio_ttbjttbb*Full_C_Xsecttbb_Xsecttjj",
    // 							     RooArgList(*CRatio_ttbjttbb,*Full_C_Xsecttbb_Xsecttjj));
    
    // RooFormulaVar *Nttjj = new RooFormulaVar("Nttjj","Nttjj","Full_Xsecttjj*Lumi_nom",RooArgList(*Full_Xsecttjj,*Lumi_nom));
    
    // RooFormulaVar *Full_Nttjj = new RooFormulaVar("Full_Nttjj","Full_Nttjj","Nttjj*Effttjj_nom*Accttjj_nom",RooArgList(*Nttjj,*Effttjj_nom,*Accttjj_nom));
    
    // RooFormulaVar *Full_k = new RooFormulaVar("Full_k","Norm k","Full_Xsecttjj*0.00345",RooArgList(*Full_Xsecttjj)); // 0.00345 from 1/290.2 
    
    // RooFormulaVar *Full_Bkgtt = new RooFormulaVar("Full_Bkgtt","Norm k x n_Bkgtt_var","Full_k*n_Bkgtt_var",RooArgList(*Full_k,*n_Bkgtt_var));
      

    // // RooAddPdf *m_ttjj = new RooAddPdf("m_ttjj","ttjj Model", 
    // // 				      RooArgList(*HisPdfSys[ttbb][0][btagcfI][Nom],*HisPdfSys[ttb][0][btagcfI][Nom],*HisPdfSys[ttccLF][0][btagcfI][Nom]), 
    // // 				      RooArgList(*Full_C_Xsecttbb_Xsecttjj,*CRttbjttbbFR_ttbbttjj));
    // // RooAddPdf *m_total = new RooAddPdf("m_total","Total Model",
    // // 				       RooArgList(*m_ttjj,*HisPdfSys[Bkgtt][0][btagcfI][Nom],*HisPdfSys[BkgOther][0][btagcfI][Nom]), 
    // // 				       RooArgList(*Full_Nttjj,*Full_Bkgtt,*n_BkgOther_var));

    // RooAddPdf *m_ttjj = new RooAddPdf("m_ttjj","ttjj Model", 
    // 				      RooArgList(*PInterSyst[ttbb][0],*PInterSyst[ttb][0],*PInterSyst[ttccLF][0]), 
    // 				      RooArgList(*Full_C_Xsecttbb_Xsecttjj,*CRttbjttbbFR_ttbbttjj));
    // RooAddPdf *m_total = new RooAddPdf("m_total","Total Model",
    // 				       RooArgList(*m_ttjj,*PInterSyst[Bkgtt][0],*PInterSyst[BkgOther][0]), 
    // 				       RooArgList(*Full_Nttjj,*Full_Bkgtt,*n_BkgOther_var));

    // //RooGaussian *GaussCon = new RooGaussian("GaussCon","Gauss Con",*alpha,0,1);
    // // RooGaussian GaussCon("GaussCon","Gauss Con",,0,1);

    // // RooProdPdf *m_totalCon = new RooProdPdf("m_totalCon","Total Model Constrained",
    // // 					    RooArgList(*m_total,*GaussCon));

    // // ----------------------------------------------------------------------------------
    // // ----------------------------------------------------------------------------------


    /////////////////////////////////////
    // WorkSpace 2: Fit with Nuisance  //
    /////////////////////////////////////

    RooWorkspace *WS_NuMo;
    WS_NuMo = new RooWorkspace("Fit CSV2  with Nuisance Parameters");

    for(int isam=ttbb; isam<=BkgOther; isam++) WS_NuMo->import(*PInterSyst[isam][0]);

    // Parameters
    WS_NuMo->import(*CRatio_ttbjttbb);

    WS_NuMo->import(*n_ttjj_var);
    WS_NuMo->import(*n_ttbb_var);
    WS_NuMo->import(*n_ttbj_var);
    WS_NuMo->import(*n_ttccLF_var);
    WS_NuMo->import(*n_Bkgtt_var);
    WS_NuMo->import(*n_BkgOther_var);

    // Visible and Full Ph-Sp
    WS_NuMo->import(*Lumi_nom);
    WS_NuMo->import(*Accttjj_nom);
    WS_NuMo->import(*Accttbb_nom);
    WS_NuMo->import(*Effttjj_nom);
    WS_NuMo->import(*Effttbb_nom);
    WS_NuMo->import(*Full_Xsecttjj);
    WS_NuMo->import(*Full_Xsecttbb_Xsecttjj);
    WS_NuMo->import(*Full_C_Xsecttbb_Xsecttjj);

    // // -- Efficiency with Nuisance Parameters
    // WS_NuMo->factory("prod::Effttbb(Effttbb_nom,EttbbbetaScale,EttbbbetaPS,EttbbbetaPDF)");
    // WS_NuMo->factory("prod::Effttjj(Effttjj_nom,EttjjbetaScale,EttjjbetaPS,EttjjbetaPDF)");
    // // -- Efficiencies and Acceptance rations 
    WS_NuMo->factory("expr::Effttbbttjj('Effttbb_nom/Effttjj_nom',Effttbb_nom,Effttjj_nom)");
    WS_NuMo->factory("expr::Accttbbttjj('Accttbb_nom/Accttjj_nom',Accttbb_nom,Accttjj_nom)");

    WS_NuMo->factory("expr::SFbb('1.0/n_ttbb_var',n_ttbb_var)");
    WS_NuMo->factory("expr::SFbj('1.0/n_ttbj_var',n_ttbj_var)");
    WS_NuMo->factory("expr::SFccLF('1.0/n_ttccLF_var',n_ttccLF_var)");

    // // Full Ph-Sp
    WS_NuMo->factory("prod::Full_C_XsecRatiottbbttjj(Effttbbttjj,Accttbbttjj,Full_C_Xsecttbb_Xsecttjj)");
    WS_NuMo->factory("expr::Full_Nttjj('Full_Xsecttjj*Lumi_nom', Full_Xsecttjj,Lumi_nom)");
    WS_NuMo->factory("prod::Full_k(Full_Xsecttjj,0.00345)"); // 0.00345 from 1/290.2

    // // -- Scenario2

    // // WS_NuMo->factory("ASUM::Full_C_ttjj_PInterSyst2(prod(Full_C_XsecRatiottbbttjj,SFbb)*ttbb_PInterSyst2,prod(SFbj,Full_C_XsecRatiottbbttjj,CRatio_ttbjttbb)*ttb_PInterSyst2, ttccLF_PInterSyst2)");
    // // WS_NuMo->factory("ASUM::Full_C_TotModel2_sys(prod(Full_Nttjj,Effttjj_nom,Accttjj_nom)*Full_C_ttjj_PInterSyst2, prod(Full_k)*Bkgtt_PInterSyst2,BkgOther_PInterSyst2)");


    WS_NuMo->factory("ASUM::Full_C_ttjj_NoNu(Full_C_XsecRatiottbbttjj*ttbb_hisPdf2btagcfINom,prod(Full_C_XsecRatiottbbttjj,CRatio_ttbjttbb)*ttb_hisPdf2btagcfINom, ttccLF_hisPdf2btagcfINom)");
    WS_NuMo->factory("ASUM::Full_C_TotModel2_NoNu(prod(Full_Nttjj,Effttjj_nom,Accttjj_nom)*Full_C_ttjj_NoNu, prod(Full_k,n_Bkgtt_var)*Bkgtt_hisPdf2btagcfINom, n_BkgOther_var*BkgOther_hisPdf2btagcfINom)");


    WS_NuMo->factory("ASUM::Full_C_ttjj_PInterSyst2(Full_C_XsecRatiottbbttjj*ttbb_PInterSyst2,prod(Full_C_XsecRatiottbbttjj,CRatio_ttbjttbb)*ttb_PInterSyst2, ttccLF_PInterSyst2)");
    WS_NuMo->factory("ASUM::Full_C_TotModel2_sys(prod(Full_Nttjj,Effttjj_nom,Accttjj_nom)*Full_C_ttjj_PInterSyst2, prod(Full_k,n_Bkgtt_var)*Bkgtt_PInterSyst2, n_BkgOther_var*BkgOther_PInterSyst2)");

    // Constrain the models
    TString sTotModel2_sys  = "PROD::TotModel2_SysCons(Full_C_TotModel2_sys"; 
    
    // for(int isys=0; isys<14; isys++) sTotModel2_sys  += ",Gaussian(0,AlphaSys2"  + SystNam[isys] + ",1)"; 
    sTotModel2_sys  += ",Gaussian(0,alpha,100)"; 
    // To close the model definition
    sTotModel2_sys  += ")"; 
    
    std::cout << sTotModel2_sys  << std::endl;

    WS_NuMo->factory(sTotModel2_sys);

    // without NP

    // WS_NuMo->factory("ASUM::Full_C_ttjj_PDF(Full_C_XsecRatiottbbttjj*ttbb_PInterSyst2,prod(Full_C_XsecRatiottbbttjj,CRatio_ttbjttbb)*ttb_PInterSyst2, ttccLF_PInterSyst2)");
    // WS_NuMo->factory("ASUM::Full_C_TotModel2_sys(prod(Full_Nttjj,Effttjj_nom,Accttjj_nom)*Full_C_ttjj_PInterSyst2, prod(Full_k,n_Bkgtt_var)*Bkgtt_PInterSyst2, n_BkgOther_var*BkgOther_PInterSyst2)");



    WS_NuMo->Print();

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    // histograms: Data 
    RooDataHist *DataHis;
    DataHis = new RooDataHist("data_his2",   "data Histogram",   
    			      *CSV2, InFile[data].hist1D[0][ch]);

    // //m_total->fitTo(*DataHis);

    // // WS_NuMo->var("alpha")->setConstant(kTRUE);

    // RooAbsReal *nll_NuMo = m_total->createNLL(*DataHis);
    RooAbsReal *nll_NuMo = WS_NuMo->pdf("TotModel2_SysCons")->createNLL(*DataHis);
    // RooAbsReal *nll_NuMo = WS_NuMo->pdf("Full_C_TotModel2_NoNu")->createNLL(*DataHis);
    // // Create MINUIT interface object
    RooMinimizer m_NuMo(*nll_NuMo);
    // // Call MIGRAD to minimize the likelihood
    m_NuMo.migrad();
    // // Run HESSE to calculate errors from d2L/dp2
    m_NuMo.hesse();
    // // More precise estimation of the parameter uncertainties
    // // MINOS uses a likelihood scan
    // //m_NuMo.minos();

    // // RooProfileLL pll("pll", "", *nll_NuMo, RooArgSet(*WS_NuMo->var("Full_Xsecttjj"),*WS_NuMo->var("Full_C_Xsecttbb_Xsecttjj")));


    // //RooAbsReal *pll_NuMo = nll_NuMo->createProfile(*WS_NuMo->var("Full_Xsecttjj"));
    //  RooPlot *plot_xsec_NuMo = WS_NuMo->var("Full_Xsecttjj")->frame(); 
    // // pll_NuMo->plotOn(plot_xsec_NuMo,LineColor(kRed)) ;
    // nll_NuMo->plotOn(plot_xsec_NuMo,LineColor(kBlue),ShiftToZero()) ;
    // // // // // Save the current fit result
    // // // // RooFitResult* r_NuMo = m_NuMo.save() ;

    //  RooPlot *plot_ratio_NuMo = WS_NuMo->var("Full_C_Xsecttbb_Xsecttjj")->frame(); 
    // // pll_NuMo->plotOn(plot_xsec_NuMo,LineColor(kRed)) ;
    // nll_NuMo->plotOn(plot_ratio_NuMo,LineColor(kBlue),ShiftToZero()) ;
    // // // // // Save the current fit result
    // // // // RooFitResult* r_NuMo = m_NuMo.save() ;

    // TCanvas *c1 = new TCanvas("Ca","Canvas");
    // c1->Divide(2,1);
    // c1->SetGrid();

    // c1->cd(1);

    // // plot_xsec_NuMo->SetMaximum(10);
    // // TLine *line_1 = new TLine(180,5,210,5);
    // plot_xsec_NuMo->SetMaximum(5.);
    // plot_xsec_NuMo->SetMinimum(0);
    // plot_xsec_NuMo->GetXaxis()->SetLimits(180,220);
    // // line_1->SetLineColor(kBlue);
    // plot_xsec_NuMo->Draw();
    // //plot_xsec_NuMo->GetXaxis()->SetRangeUser(195,200);
    // // line_1->Draw("SAME");
    // // // float nuival = WS_NuMo->var("alpha")->getVal();

    // c1->cd(2);

    // plot_ratio_NuMo->SetMaximum(5.);
    // plot_ratio_NuMo->SetMinimum(0);
    // plot_ratio_NuMo->GetXaxis()->SetLimits(0.005,0.03);
    // // line_1->SetLineColor(kBlue);
    // plot_ratio_NuMo->Draw();

    // c1->SaveAs("Canvas_Test.pdf");


    // WS_NuMo->var("alpha")->setVal(0);
    // WS_NuMo->var("alpha")->setConstant(kTRUE);


    // // // Call MIGRAD to minimize the likelihood
    // m_NuMo.migrad();
    // // // // // Run HESSE *to calculate errors from d2L/dp2
    // m_NuMo.hesse();

    // RooStats::ModelConfig *mc = (RooStats::ModelConfig*) WS_NuMo->obj("TotModel2_SysCons");
    // RooAbsData* data =  WS_NuMo->data("data_his2");

    // ProfileLikelihoodCalculator pl(*DataHis,*mc);
    // pl.SetConfidenceLevel(0.683); // 68% interval
    // LikelihoodInterval* interval = pl.GetInterval();
  } // for(ch)
  
  fclose(fResult);
  
}


HistoFit LoadSample(TString FileName){

  
  TFile* fInput  = new TFile(dirnameIn + fl + "_" + FileName + ".root");
  if(!fInput->GetFile()){
    std::cerr << dirnameIn + fl +  FileName << " not Found!!!" << std::endl;
    std::exit(0);
  }  

  TString Cut = "2btag";   
  TString nch[3] = {"mujets","ejets","ljets"};
  TString nca[3] = {"hKinAdd1CSV","hKinAdd2CSV","h2DKinAddCSV"};
  
  HistoFit Output;
  
  for(int irca=0; irca<3;irca++){
    for(int irch=0; irch<2;irch++){
      if (irca<2){
	Output.hist1D[irca][irch] = (TH1D*)fInput->Get(Cut + "/" + nch[irch] + "/" + nca[irca] + "_" + nch[irch] + "_" + Cut)->Clone();
	// Avoid negative entries (from aMC@NLO MC)
	for(int ibin = 1; ibin<= Output.hist1D[irca][irch]->GetNbinsX();ibin++){
	  double binCont = Output.hist1D[irca][irch]->GetBinContent(ibin);
	  if (binCont < 0.0) Output.hist1D[irca][irch]->SetBinContent(ibin,0.0);
	}

	// Number of Events
	if (irca == 0) Output.events[irch] = Output.hist1D[0][irch]->Integral(0,1000);
	
	// Normalization
	//Output.hist1D[irca][irch]->Scale(1.0/Output.hist1D[irca][irch]->Integral());

      } // if(irca)
      else {
	Output.hist2D[irch] = (TH2D*)fInput->Get(Cut + "/" + nch[irch] + "/" + nca[irca] + "_" + nch[irch] + "_" + Cut)->Clone();
	// Avoid negative entries (from aMC@NLO MC)
	for(int ibinX = 1; ibinX<= Output.hist2D[irch]->GetNbinsX();ibinX++){
	  for(int ibinY = 1; ibinY<= Output.hist2D[irch]->GetNbinsX();ibinY++){
	    double binCont = Output.hist2D[irch]->GetBinContent(ibinX,ibinY);
	    if (binCont < 0.0) Output.hist2D[irch]->SetBinContent(ibinX,ibinY,0.0);
	  }
	}
	// Normalization
	//Output.hist2D[irch]->Scale(1.0/Output.hist2D[irch]->Integral());

      } // else

    } // for(irch)
 
    if (irca<2){
      Output.hist1D[irca][2] = (TH1D*)Output.hist1D[irca][0]->Clone();
      Output.hist1D[irca][2]->Add(Output.hist1D[irca][1]);
      // Normalization
      //Output.hist1D[irca][2]->Scale(1.0/Output.hist1D[irca][2]->Integral());

    }
    else{
      Output.hist2D[2] = (TH2D*)Output.hist2D[0]->Clone();
      Output.hist2D[2]->Add(Output.hist2D[1]);    
      // Normalization
      //Output.hist2D[2]->Scale(1.0/Output.hist2D[2]->Integral());
    }
  } // for(irca)

  Output.events[2] = Output.events[0] + Output.events[1];

  cout << "YIELDS = " << Output.events[0] << "  " << Output.events[1] << "  " << Output.events[2] << endl;
  
  cout << "All the histograms from " << dirnameIn + fl +  FileName << " have been loaded successfully!!!" << endl;

  // ----------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------
  
  if(!FileName.Contains("Data")){
        
    for(unsigned int isys = 0; isys<14; isys++){

      if(isys != btagcfI) continue;
      
      TFile* fInputSys[3];
      fInputSys[Nom]  = new TFile(dirnameIn + fl + "_" + FileName + ".root");
      
      fInputSys[Up]    = new TFile(dirnameIn + fl + "_" + FileName + "_SYS_" + SystNam[isys] + "_Up.root");
      fInputSys[Down]  = new TFile(dirnameIn + fl + "_" + FileName + "_SYS_" + SystNam[isys] + "_Down.root");
      
      if(!fInputSys[Nom]->GetFile() || !fInputSys[Up]->GetFile() || !fInputSys[Down]->GetFile()){
	std::cerr << SystNam[isys] << " systematic variation for " << dirnameIn + fl +  FileName << " not Found!!!" << std::endl;
	std::exit(0);
      }  

      for(unsigned int ivar = 0; ivar<3; ivar++){

	std::cout << "Loading histograms for " << fInputSys[ivar]->GetName() << std::endl; 
	
	for(int irca=0; irca<3;irca++){
	  for(int irch=0; irch<2;irch++){

	    if (irca<2){
	      Output.hsyst1D[isys][ivar][irca][irch] = (TH1D*)fInputSys[ivar]->Get(Cut + "/" + nch[irch] + "/" +  nca[irca] + "_" + nch[irch] + "_" + Cut)->Clone();
	      // Avoid negative entries (from aMC@NLO MC)
	      for(int ibin = 1; ibin<= Output.hsyst1D[isys][ivar][irca][irch]->GetNbinsX();ibin++){
		double binCont = Output.hsyst1D[isys][ivar][irca][irch]->GetBinContent(ibin);
		if (binCont < 0.0) Output.hsyst1D[isys][ivar][irca][irch]->SetBinContent(ibin,0.0);
	      }
	      // Normalize wrt the number of event
	      Output.hsyst1D[isys][ivar][irca][irch]->Scale(1.0/Output.hsyst1D[isys][ivar][irca][irch]->Integral());
	    } // if(irca)

	    else {
	      Output.hsyst2D[isys][ivar][irch] = (TH2D*)fInputSys[ivar]->Get(Cut + "/" + nch[irch] + "/" + nca[irca] + "_" + nch[irch] + "_" + Cut)->Clone();
	      // Avoid negative entries (from aMC@NLO MC)
	      for(int ibinX = 1; ibinX<= Output.hsyst2D[isys][ivar][irch]->GetNbinsX();ibinX++){
		for(int ibinY = 1; ibinY<= Output.hsyst2D[isys][ivar][irch]->GetNbinsX();ibinY++){
		  double binCont = Output.hsyst2D[isys][ivar][irch]->GetBinContent(ibinX,ibinY);
		  if (binCont < 0.0) Output.hsyst2D[isys][ivar][irch]->SetBinContent(ibinX,ibinY,0.0);
		}
	      }
	      // Normalize wrt the number of events
	      Output.hsyst2D[isys][ivar][irch]->Scale(1.0/Output.hsyst2D[isys][ivar][irch]->Integral());
	    } // else
	  } // for(irch)
	  
	  if (irca<2){
	    Output.hsyst1D[isys][ivar][irca][2] = (TH1D*)Output.hsyst1D[isys][ivar][irca][0]->Clone();
	    Output.hsyst1D[isys][ivar][irca][2]->Add(Output.hsyst1D[isys][ivar][irca][1]);
	    // Normalize wrt the number of events
	    Output.hsyst1D[isys][ivar][irca][2]->Scale(1.0/Output.hsyst1D[isys][ivar][irca][2]->Integral());
	  }
	  else{
	    Output.hsyst2D[isys][ivar][2] = (TH2D*)Output.hsyst2D[isys][ivar][0]->Clone();
	    Output.hsyst2D[isys][ivar][2]->Add(Output.hsyst2D[isys][ivar][1]);    
	    // Normalize wrt the number of events
	    Output.hsyst2D[isys][ivar][2]->Scale(1.0/Output.hsyst2D[isys][ivar][2]->Integral());    
	  }
	} // for(irca)
      } // for(ivar)
    } // for(isys)
  } // if(!data)
  
  //fInput->Close();
  
  return Output;
}
