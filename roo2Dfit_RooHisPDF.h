#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include <vector>
#include <string>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <set>
#include "tdrstyle.C"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"


typedef struct HistoFit{
  TH1D     *hist1D[2][3];
  TH2D     *hist2D[3];
  double    events[3];
  TH1     *hsyst1D[14][3][2][3];
  TH2D     *hsyst2D[14][3][3];

} HistoFit;


TString dirnameIn= "TopResults/";
TString fl  = "hSF-NuiPar10M-v0_Tree_LepJets_PileUp17_v8-0-4_Spring16-80X_36814pb-1";  
//TString fl  = "hSF-NuiPar-v1_Tree_LepJets_KFCSVOrder01NoSkim_v8-0-1_Spring16-80X_15920pb-1";  

float vEffttjj[3];
float vEffttbb[3];
float vAccttjj[3];
float vAccttbb[3];
float vLumi;

enum Channel{muJets,eJets,LepJets};
TString name_ch[3] = {"muJets","eJets","LepJets"};

enum FileSample{data, 
		ttbb, ttb, ttccLF, 
		Bkgtt, BkgOther,
		ttcc, ttLF, ttjj,
		BkgFull};
TString SamNam[18] = {"data", 
		      "ttbb", "ttb", "ttccLF", 
		      "Bkgtt", "BkgOther",
		      "ttcc", "ttLF", "ttjj",
		      "BkgFull"};

int colors[15] = { 1 ,
		   TColor::GetColor("#660000"), TColor::GetColor("#ffcc00"), TColor::GetColor("#cc6600"),
		   TColor::GetColor("#ff6565"), TColor::GetColor("#8000ff"),
		   TColor::GetColor("#cc6600"), TColor::GetColor("#ff0000"), TColor::GetColor("#3366ff"),
		   TColor::GetColor("#33cc33")};


enum FileSystematic{PileUp,
		    JER, JES,
		    LES, LepSF,
		    btagjes, btaglf, btaghf, btaghfsI, btaghfsII, btaglfsI, btaglfsII, btagcfI, btagcfII};
TString SystNam[14]= {"PileUp",
		      "JER", "JES",
		      "LES", "LepSF",
		      "btagjes", "btaglf", "btaghf", "btaghfsI", "btaghfsII", "btaglfsI", "btaglfsII", "btagcfI", "btagcfII"};

enum TheorySystematic{Lumi,Scale,PS,PDF};
TString SystThNam[3]= {"Scale","PS","PDF"};

enum FileVariation{Nom, Down, Up};
TString VarNam[3]= {"Nom", "Down", "Up"};

HistoFit LoadSample(TString FileName);
RooPlot *PlotPDFModel(RooRealVar *var, RooHistPdf *Model, RooDataHist *DataHis, bool PlotData, TLegend *leg);
RooPlot *PlotPDF_NuMo(RooRealVar *var, RooWorkspace *WS,  RooDataHist *DataHis, bool PlotData, TLegend *leg);
TH1 *HistoPDF_NuMo(TString var, RooWorkspace *WS, RooDataHist *DataHis, bool PlotData, TLegend *leg);

void NuiParameterPlot_2 (TCanvas *canvas_sys_NuMo, RooWorkspace *WS, TString cat = "ttbb", TString par = "btagcfI");
void NuiParameterPlot_3 (TCanvas *canvas_sys_NuMo, RooWorkspace *WS, TString cat = "ttbb", TString par = "btagcfI");
//TCanvas *NuiParameterPlot_23(RooWorkspace *WS, TString cat = "ttbb", TString par = "btagcfI");
void NuiParameterPlot_23(TCanvas *canvas_sys_NuMo, RooWorkspace *WS, TString cat = "ttbb", TString par = "btagcfI");
 
