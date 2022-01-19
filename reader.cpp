#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

#include <TROOT.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include "TStyle.h"

#include "DVCSEvent.h"

#include <iostream>
#include <fstream>

using namespace HepMC3;

//main
int main(int argc, char **argv) {


    double Q2Bins[13] = {1.0, 1.77828, 3.16228, 5.62341, 10.0, 17.7828, 31.6228, 56.2341, 100.0, 177.828, 316.228, 562.341, 1000.0};

    double xBjBins[20] = {0.0001, 0.000158489, 0.000251189, 0.000398107, 0.000630957, 0.001, 0.00158489, 0.00251189, 0.00398107, 0.00630957, 0.01, 0.0158489, 0.0251189, 0.0398107, 0.0630957, 0.1, 0.158489, 0.251189, 0.398107, 0.630957};

    double tBins[12] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3};

    int i,j,k,l,m,n;

    double photonEnergyCut = 0.0;
    double electronEnergyCut = 0.0;

    double phi_low = 0.03;
    double phi_high = 6.2531;

    //int Q2BinNumber = 13;
    //int xBjBinNumber = 20;

    int Q2BinNumber = sizeof(Q2Bins)/sizeof(Q2Bins[0]);
    int xBjBinNumber = sizeof(xBjBins)/sizeof(xBjBins[0]);
    int tBinNumber = sizeof(tBins)/sizeof(tBins[0]);

    int photonEBinNumber = 39;
    int electronEBinNumber = 39;
    int PhiBinNumber = 39;

    //histograms

    std::vector<TH1*> h_ALU(4000, nullptr);

    std::vector<TH1D*> h_phi_parallel(4000, nullptr);
    std::vector<TH1D*> h_phi_antiparallel(4000, nullptr);

    // 1D histograms for ALL (DVCS + BH + INT)  
    TH1F *xB_all = new TH1F("xB_all", "", 19, xBjBins);
    TH1F *Q2_all = new TH1F("Q2_all", "", 12, Q2Bins);

    TH1F *t_parallel = new TH1F("t_parallel", "", 11, 0.2, 1.3);
    TH1F *t_antiparallel = new TH1F("t_antiparallel", "", 11, 0.2, 1.3);
    TH1F *t_all = new TH1F("t_all", "", 11, 0.2, 1.3);

    // 2D histograms for ALL (DVCS + BH + INT) 
    TH2F *photonEanglePhi_parallel = new TH2F("Ephoton_phi_parallel","",photonEBinNumber,0,10,PhiBinNumber,phi_low,phi_high);
    TH2F *photonEanglePhi_antiparallel = new TH2F("Ephoton_phi_antiparallel","",photonEBinNumber,0,10,PhiBinNumber,phi_low,phi_high);

    TH2F *electronEanglePhi_parallel = new TH2F("Ee_phi_parallel","",electronEBinNumber,0,10,PhiBinNumber,phi_low,phi_high);
    TH2F *electronEanglePhi_antiparallel = new TH2F("Ee_phi_antiparallel","",electronEBinNumber,0,10,PhiBinNumber,phi_low,phi_high);

    TH2F *photonE_Q2_parallel = new TH2F("Ephoton_Q2_parallel","",photonEBinNumber,0,10,12,Q2Bins);
    TH2F *photonE_Q2_antiparallel = new TH2F("Ephoton_Q2_antiparallel","",photonEBinNumber,0,10,12,Q2Bins);

    TH2F *electronE_Q2_parallel = new TH2F("Ee_Q2_parallel","",electronEBinNumber,0,10,12,Q2Bins);
    TH2F *electronE_Q2_antiparallel = new TH2F("Ee_Q2_antiparallel","",electronEBinNumber,0,10,12,Q2Bins);

    TH2F *photonE_xB_parallel = new TH2F("photonE_xB_parallel","",photonEBinNumber,0,10,19,xBjBins);
    TH2F *photonE_xB_antiparallel = new TH2F("photonE_xB_antiparallel","",photonEBinNumber,0,10,19,xBjBins);

    TH2F *electronE_xB_parallel = new TH2F("electronE_xB_parallel","",electronEBinNumber,0,10,19,xBjBins);
    TH2F *electronE_xB_antiparallel = new TH2F("electronE_xB_antiparallel","",electronEBinNumber,0,10,19,xBjBins);

    TH2F *photonE_t_parallel = new TH2F("Ephoton_t_parallel","",photonEBinNumber,0,10,11,tBins);
    TH2F *photonE_t_antiparallel = new TH2F("Ephoton_t_antiparallel","",photonEBinNumber,0,10,11,tBins);

    TH2F *electronE_t_parallel = new TH2F("Ee_t_parallel","",electronEBinNumber,0,10,11,tBins);
    TH2F *electronE_t_antiparallel = new TH2F("Ee_t_antiparallel","",electronEBinNumber,0,10,11,tBins);

    TH2D *xBQ2_parallel = new TH2D("xBQ2 parallel helicity","Parallel Helicity",19,xBjBins,12,Q2Bins);
    TH2D *xBQ2_antiparallel = new TH2D("xBQ2 antiparallel helicity","Antiparallel Helicity",19,xBjBins,12,Q2Bins);
    TH2D *xBQ2_all = new TH2D("xBQ2 all","",19,xBjBins,12,Q2Bins);

    TH2D *xB_t_parallel = new TH2D("xB_t parallel helicity","",19,xBjBins,11,tBins);
    TH2D *xB_t_antiparallel = new TH2D("xB_t antiparallel helicity","",19,xBjBins,11,tBins);

    TH2D *xB_phi_parallel = new TH2D("xB_phi parallel helicity","",19,xBjBins,PhiBinNumber,phi_low,phi_high);
    TH2D *xB_phi_antiparallel = new TH2D("xB_phi antiparallel helicity","",19,xBjBins,PhiBinNumber,phi_low,phi_high);

    TH2D *Q2_t_parallel = new TH2D("Q2_t parallel helicity","",12,Q2Bins,11,tBins);
    TH2D *Q2_t_antiparallel = new TH2D("Q2_t antiparallel helicity","",12,Q2Bins,11,tBins);

    TH2D *Q2_phi_parallel = new TH2D("Q2_phi parallel helicity","",12,Q2Bins,PhiBinNumber,phi_low,phi_high);
    TH2D *Q2_phi_antiparallel = new TH2D("Q2_phi antiparallel helicity","",12,Q2Bins,PhiBinNumber,phi_low,phi_high);

    TH2D *t_phi_parallel = new TH2D("t_phi parallel helicity","",11,tBins,PhiBinNumber,phi_low,phi_high);
    TH2D *t_phi_antiparallel = new TH2D("t_phi antiparallel helicity","",11,tBins,PhiBinNumber,phi_low,phi_high);

    // 3D histograms for ALL (DVCS + BH + INT) 
    TH3D *xB_Q2_t_all = new TH3D("xB_Q2_t all","",19,xBjBins,12,Q2Bins,11,tBins);

    // 1D histograms for BH 
    TH1F *xB_BH = new TH1F("xB_BH", "", 19, xBjBins);
    TH1F *Q2_BH = new TH1F("Q2_BH", "", 12, Q2Bins);
    TH1F *t_BH = new TH1F("t_BH", "", 11, 0.2, 1.3);

    // 2D histograms for BH 
    TH2D *xBQ2_BH = new TH2D("xBQ2 BH","",19,xBjBins,12,Q2Bins);

    // 3D histograms for BH 
    TH3D *xB_Q2_t_BH = new TH3D("xB_Q2_t BH","",19,xBjBins,12,Q2Bins,11,tBins);

    // 1D histograms for BH subtracted ALL  
    TH1F *xB_BH_sub = new TH1F("xB_BH_sub", "", 19, xBjBins);
    TH1F *Q2_BH_sub = new TH1F("Q2_BH_sub", "", 12, Q2Bins);
    TH1F *t_BH_sub = new TH1F("t_BH_sub", "", 11, 0.2, 1.3);
    
    // 2D histograms for BH subtracted ALL  
    TH2D *xBQ2_BH_sub = new TH2D("xBQ2 BH_sub","",19,xBjBins,12,Q2Bins);

    // 3D histograms for BH subtracted ALL  
    TH3D *xB_Q2_t_BH_sub = new TH3D("xB_Q2_t BH_sub","",19,xBjBins,12,Q2Bins,11,tBins);

    // histograms for A_LU

    for(k=0;k<4000; k++) {
        //h_ALU[k] = new TH1(TString::Format("h_ALU_%d", k), "", 10, 0.0, 6.2831);

        h_phi_parallel[k] = new TH1D(TString::Format("h_phi_parallel_%d", k), "", 11, phi_low, phi_high);
        h_phi_antiparallel[k] = new TH1D(TString::Format("h_phi_antiparallel_%d", k), "", 11, phi_low, phi_high);

    }


//****************************************************** Fill the Histograms for ALL ******************************************************

    for(m=11; m<=38; m++) {
        ReaderAscii inputFile("10_100_ALL_parallel_3deg/even." + to_string(m) + ".txt");

        //loop over event
        size_t iEvent = 0;

        while(! inputFile.failed()) {

                //event
                GenEvent evt(Units::GEV,Units::MM);
                //read event from input file
                inputFile.read_event(evt);
                //if reading failed - exit loop
                if(inputFile.failed() ) break;

                //DVCS event
                DVCSEvent dvcsEvent(evt, 1);

                //fill

            if(dvcsEvent.getGammaOut_E() > photonEnergyCut && dvcsEvent.getEOut_E() > electronEnergyCut) { //Apply a cut on the photon and electron energies

	       xB_all->Fill(dvcsEvent.getXb());

	       Q2_all->Fill(dvcsEvent.getQ2());

	       t_parallel->Fill(-1. * dvcsEvent.getT());

  	       t_all->Fill(-1. * dvcsEvent.getT());

                photonEanglePhi_parallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getPhi());

                electronEanglePhi_parallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getPhi());

                photonE_xB_parallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getXb());

                electronE_xB_parallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getXb());

                photonE_Q2_parallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getQ2());

                electronE_Q2_parallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getQ2());

                photonE_t_parallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getT());

                electronE_t_parallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getT());

                xBQ2_parallel->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2());

	        xBQ2_all->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2());

                xB_t_parallel->Fill(dvcsEvent.getXb(),dvcsEvent.getT());

                xB_phi_parallel->Fill(dvcsEvent.getXb(),dvcsEvent.getPhi());

                Q2_t_parallel->Fill(dvcsEvent.getQ2(),dvcsEvent.getT());

                Q2_phi_parallel->Fill(dvcsEvent.getQ2(),dvcsEvent.getPhi());

                t_phi_parallel->Fill(dvcsEvent.getT(),dvcsEvent.getPhi());

	        xB_Q2_t_all->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2(),-1.0 * dvcsEvent.getT());

             }    

	        // Fill the histograms for BSA (parallel helicity) 
                if(dvcsEvent.getGammaOut_E() > photonEnergyCut && dvcsEvent.getEOut_E() > electronEnergyCut) { //Apply a cut on the photon and electron energies
                    for(i=0; i<(Q2BinNumber-1); i++) {
                        for(j=0; j<(xBjBinNumber-1); j++){
                            for(k=0; k<(tBinNumber-1); k++){
                                if((Q2Bins[i] < dvcsEvent.getQ2()) && (dvcsEvent.getQ2() < Q2Bins[i+1])) {
                                    if((xBjBins[j] < dvcsEvent.getXb()) && (dvcsEvent.getXb() < xBjBins[j+1])) {
                                        if((tBins[k] < (-1.0 * dvcsEvent.getT())) && ((-1.0 * dvcsEvent.getT()) < tBins[k+1])) {
                                            h_phi_parallel[i*(xBjBinNumber-1)*(tBinNumber-1)+j*(tBinNumber-1)+k]->Fill(dvcsEvent.getPhi());
                                        }
                                    }
                                }
                            }

                        }

                    }
                }

                //id
                iEvent++;
        }
        //close file
        inputFile.close();
    }


    for(m=11; m<=38; m++) {
        ReaderAscii inputFile("10_100_ALL_antiparallel_3deg/even." + to_string(m) + ".txt");

        //loop over event
        size_t iEvent = 0;

        while(! inputFile.failed()) {

            //event
            GenEvent evt(Units::GEV,Units::MM);

            //read event from input file
            inputFile.read_event(evt);
            //if reading failed - exit loop
            if(inputFile.failed() ) break;
            //DVCS event
            DVCSEvent dvcsEvent(evt, 1);

            //fill

        if(dvcsEvent.getGammaOut_E() > photonEnergyCut && dvcsEvent.getEOut_E() > electronEnergyCut) { //Apply a cut on the photon and electron energies    

	   xB_all->Fill(dvcsEvent.getXb());

	   Q2_all->Fill(dvcsEvent.getQ2());

	   t_antiparallel->Fill(-1. * dvcsEvent.getT());

	   t_all->Fill(-1. * dvcsEvent.getT());

            photonEanglePhi_antiparallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getPhi());

            electronEanglePhi_antiparallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getPhi());

            photonE_xB_antiparallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getXb());

            electronE_xB_antiparallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getXb());

            photonE_Q2_antiparallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getQ2());

            electronE_Q2_antiparallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getQ2());

            photonE_t_antiparallel->Fill(dvcsEvent.getGammaOut_E(),dvcsEvent.getT());

            electronE_t_antiparallel->Fill(dvcsEvent.getEOut_E(),dvcsEvent.getT());

            xBQ2_antiparallel->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2());

	    xBQ2_all->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2());

            xB_t_antiparallel->Fill(dvcsEvent.getXb(),dvcsEvent.getT());

            xB_phi_antiparallel->Fill(dvcsEvent.getXb(),dvcsEvent.getPhi());

            Q2_t_antiparallel->Fill(dvcsEvent.getQ2(),dvcsEvent.getT());

            Q2_phi_antiparallel->Fill(dvcsEvent.getQ2(),dvcsEvent.getPhi());

            t_phi_antiparallel->Fill(dvcsEvent.getT(),dvcsEvent.getPhi());

	    xB_Q2_t_all->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2(),-1.0 * dvcsEvent.getT());

        }    
	    // Fill the histograms for BSA (antiparallel helicity)
            if(dvcsEvent.getGammaOut_E() > photonEnergyCut && dvcsEvent.getEOut_E() > electronEnergyCut) { //Apply a cut on the photon and electron energies
                for(i=0; i<(Q2BinNumber-1); i++) {
                    for(j=0; j<(xBjBinNumber-1); j++){
                        for(k=0; k<(tBinNumber-1); k++){
                            if((Q2Bins[i] < dvcsEvent.getQ2()) && (dvcsEvent.getQ2() < Q2Bins[i+1])) {
                                if((xBjBins[j] < dvcsEvent.getXb()) && (dvcsEvent.getXb() < xBjBins[j+1])) {
                                    if((tBins[k] < (-1.0 * dvcsEvent.getT())) && ((-1.0 * dvcsEvent.getT()) < tBins[k+1])) {
                                        h_phi_antiparallel[i*(xBjBinNumber-1)*(tBinNumber-1)+j*(tBinNumber-1)+k]->Fill(dvcsEvent.getPhi());
                                    }
                                }
                            }
                        }

                    }

                }
            }

            //id
            iEvent++;
        }
        //close file
        inputFile.close();
    }

//****************************************************** Fill the Histograms for BH to be used in calculating the number of events (not for asymmetries) ******************************************************

    for(m=11; m<=33; m++) {
        ReaderAscii inputFile("10_100_BH_parallel_3deg/even." + to_string(m) + ".txt");

        //loop over event
        size_t iEvent = 0;

        while(! inputFile.failed()) {

                //event
                GenEvent evt(Units::GEV,Units::MM);
                //read event from input file
                inputFile.read_event(evt);
                //if reading failed - exit loop
                if(inputFile.failed() ) break;

                //DVCS event
                DVCSEvent dvcsEvent(evt, 1);

                //fill

	        xB_BH->Fill(dvcsEvent.getXb());

	        Q2_BH->Fill(dvcsEvent.getQ2());

	        t_BH->Fill(-1. * dvcsEvent.getT());

 	        xBQ2_BH->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2());

	        xB_Q2_t_BH->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2(),-1.0 * dvcsEvent.getT());

                //id
                iEvent++;
        }
        //close file
        inputFile.close();
    }


    for(m=11; m<=33; m++) {
        ReaderAscii inputFile("10_100_BH_antiparallel_3deg/even." + to_string(m) + ".txt");

        //loop over event
        size_t iEvent = 0;

        while(! inputFile.failed()) {

            //event
            GenEvent evt(Units::GEV,Units::MM);

            //read event from input file
            inputFile.read_event(evt);
            //if reading failed - exit loop
            if(inputFile.failed() ) break;
            //DVCS event
            DVCSEvent dvcsEvent(evt, 1);

            //fill

	   xB_BH->Fill(dvcsEvent.getXb());

	   Q2_BH->Fill(dvcsEvent.getQ2());

	   t_BH->Fill(-1. * dvcsEvent.getT());

	   xBQ2_BH->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2());

	   xB_Q2_t_BH->Fill(dvcsEvent.getXb(),dvcsEvent.getQ2(),-1.0 * dvcsEvent.getT());

            //id
            iEvent++;
        }
        //close file
        inputFile.close();
    }

//************************************************************************************************************

    // Subtract the BH events from the ALL
    xB_BH_sub->Add(xB_all,xB_BH,1.0,-1.0);

    Q2_BH_sub->Add(Q2_all,Q2_BH,1.0,-1.0);

    t_BH_sub->Add(t_all,t_BH,1.0,-1.0);

    xBQ2_BH_sub->Add(xBQ2_all,xBQ2_BH,1.0,-1.0);

    xB_Q2_t_BH_sub->Add(xB_Q2_t_all,xB_Q2_t_BH,1.0,-1.0);


    double limunosity = 1.0;


    // calculate the ratio of events phi=pi / phi=0 for each energy bin
    n = 40;
    double phiF, phiL, E[n], phiRatio[n];

    for(i = 1; i <= 40; i++) {
        phiF = photonEanglePhi_parallel->GetBinContent(i,1);
        phiL = photonEanglePhi_parallel->GetBinContent(i,20);

        phiRatio[i-1] = phiL/phiF;

        // cout << "RatioParallel: " << phiRatio[i-1] << " \n";

    }

    TGraph *graphRatioParallel = new TGraph(n,E,phiRatio);

    for(i = 1; i <= 40; i++) {
        phiF = photonEanglePhi_antiparallel->GetBinContent(i,1);
        phiL = photonEanglePhi_antiparallel->GetBinContent(i,20);

        phiRatio[i-1] = phiL/phiF;

        // cout << "RatioAntiparallel: " << phiRatio[i-1] << " \n";

    }

    TGraph *graphRatioAntiparallel = new TGraph(n,E,phiRatio);

    TCanvas* can_phi_ratio = new TCanvas("can_phi_ratio");
    can_phi_ratio->Divide(1, 2);

    can_phi_ratio->cd(1);
    can_phi_ratio->cd(1)->SetLogy();
    graphRatioParallel->Draw("APL");

    can_phi_ratio->cd(2);
    can_phi_ratio->cd(2)->SetLogy();
    graphRatioAntiparallel->Draw("APL");

    can_phi_ratio->Print("plots/2D/plots_phi_ratio.pdf", "pdf");


    // Write the number of events for 1D distributions

    n = tBinNumber-1;
    m = xBjBinNumber-1;
    l = Q2BinNumber-1;

    double t[n];
    double xB[m]; 
    double Q2[l];

    ofstream myfile;
    myfile.open ("1D_distributions.txt");

    myfile << " Number of events in each xBj-bin: \n\n";

    for(i = 1; i <= m; i++) {

        //xB[i-1] = xB_all->GetBinContent(i);
        xB[i-1] = xB_BH_sub->GetBinContent(i);

	myfile << "bin range: " << xBjBins[i-1] << "-" << xBjBins[i] << " Number of events: " << xB[i-1] << "\n"; 

    }

    myfile << "\n Number of events in each Q2-bin: \n\n";

    for(i = 1; i <= l; i++) {
        //Q2[i-1] = Q2_all->GetBinContent(i);
        Q2[i-1] = Q2_BH_sub->GetBinContent(i);

	myfile << "bin range: " << Q2Bins[i-1] << "-" << Q2Bins[i] << " Number of events: " << Q2[i-1] << "\n"; 

    }

    myfile << "\n Number of events in each t-bin: \n\n";
    
    for(i = 1; i <= n; i++) {
        //t[i-1] = t_all->GetBinContent(i);
        t[i-1] = t_BH_sub->GetBinContent(i);

	myfile << "bin range: " << tBins[i-1] << "-" << tBins[i] << " Number of events: " << t[i-1] << "\n"; 

    }

    myfile.close();


    // Write the number of events for 2D distributions

    m = xBjBinNumber-1;
    l = Q2BinNumber-1;

    double xB_Q2[m][l];

    // ofstream myfile;
    myfile.open ("xB_Q2_distribution.txt");
    
    for(i = 1; i <= m; i++) {

	myfile << "xB bin range: " << setw(11) << xBjBins[i-1] << "-" << setw(11) << xBjBins[i] << setw(12) << "    Number of events for each Q2 bin:  " << setw(12); 

	for(j = 1; j <= l; j++) {

	    //xB_Q2[i-1][j-1] = xBQ2_all->GetBinContent(i,j);
	    xB_Q2[i-1][j-1] = xBQ2_BH_sub->GetBinContent(i,j);

	    myfile << setw(6) << xB_Q2[i-1][j-1] << "\t";  
	}

        myfile << "\n"; 

    }

    myfile.close();



    // Write the number of events for 3D distributions

    m = xBjBinNumber-1;
    l = Q2BinNumber-1;
    n = tBinNumber-1;

    double xB_Q2_t[m][l][n];

    // ofstream myfile;
    myfile.open ("xB_Q2_t_distribution.txt");
    
    for(i = 1; i <= m; i++) {

	for(j = 1; j <= l; j++) {

	    myfile << "xB range: " << setw(11) << xBjBins[i-1] << "-" << setw(11) << xBjBins[i] << setw(12) << " Q2 range: " << setw(11) << Q2Bins[j-1] << "-" << setw(11) << Q2Bins[j] << setw(12) << "    Number of events for each t bin:  " << setw(12); 

	    for(k = 1; k <= n; k++) {

	         //xB_Q2_t[i-1][j-1][k-1] = xB_Q2_t_all->GetBinContent(i,j,k);
		xB_Q2_t[i-1][j-1][k-1] = xB_Q2_t_BH_sub->GetBinContent(i,j,k);

	        myfile << setw(6) << xB_Q2_t[i-1][j-1][k-1] << "\t";  
	    }

	myfile << "\n";

	}

        myfile << "\n\n"; 

    }

    myfile.close();


    // Plots for t-distribution
    TCanvas* can_t = new TCanvas("can_t");
    can_t->Divide(1, 3);

    gStyle->SetOptStat(0);

    can_t->cd(1);
    can_t->cd(1)->SetLogy();
    t_parallel->SetMinimum(1);
    t_parallel->Draw("hist");
    t_parallel->GetXaxis()->SetTitle("|t| (GeV^{2})");
    t_parallel->GetYaxis()->SetTitle("Number of events");
    t_parallel->GetYaxis()->CenterTitle(true);
    
    can_t->cd(2);
    can_t->cd(2)->SetLogy();
    t_antiparallel->SetMinimum(1);
    t_antiparallel->Draw("hist");
    t_antiparallel->GetXaxis()->SetTitle("|t| (GeV^{2})");
    t_antiparallel->GetYaxis()->SetTitle("Number of events");
    t_antiparallel->GetYaxis()->CenterTitle(true);

    can_t->cd(3);
    can_t->cd(3)->SetLogy();
    t_all->SetMinimum(1);
    t_all->Draw("hist");
    t_all->GetXaxis()->SetTitle("|t| (GeV^{2})");
    t_all->GetYaxis()->SetTitle("Number of events");
    t_all->GetYaxis()->CenterTitle(true);

    can_t->Print("plots/1D/plots_1D_t.pdf", "pdf");

    // Plots for photonE-Phi correlation
    TCanvas* can_E_phi = new TCanvas("can_E_phi");
    can_E_phi->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_E_phi->cd(1);
    photonEanglePhi_parallel->Draw("colz");
    can_E_phi->cd(1)->SetLogz();

    can_E_phi->cd(2);
    photonEanglePhi_antiparallel->Draw("colz");
    can_E_phi->cd(2)->SetLogz();

    can_E_phi->Print("plots/2D/plots_2D_Ephoton_phi.pdf", "pdf");

    // Plots for electronE-Phi correlation
    TCanvas* can_Ee_phi = new TCanvas("can_Ee_phi");
    can_Ee_phi->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_Ee_phi->cd(1);
    electronEanglePhi_parallel->Draw("colz");
    can_Ee_phi->cd(1)->SetLogz();

    can_Ee_phi->cd(2);
    electronEanglePhi_antiparallel->Draw("colz");
    can_Ee_phi->cd(2)->SetLogz();

    can_Ee_phi->Print("plots/2D/plots_2D_Eelectron_phi.pdf", "pdf");

    // Plots for photonE-xB correlation
    TCanvas* can_photonE_xB = new TCanvas("can_photonE_xB");
    can_photonE_xB->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_photonE_xB->cd(1);
    can_photonE_xB->cd(1)->SetLogy();
    can_photonE_xB->cd(1)->SetLogz();
    photonE_xB_parallel->Draw("colz");

    can_photonE_xB->cd(2);
    can_photonE_xB->cd(2)->SetLogy();
    can_photonE_xB->cd(2)->SetLogz();
    photonE_xB_antiparallel->Draw("colz");

    can_photonE_xB->Print("plots/2D/plots_2D_photonE_xB.pdf", "pdf");

    // Plots for electronE-xB correlation
    TCanvas* can_electronE_xB = new TCanvas("can_electronE_xB");
    can_electronE_xB->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_electronE_xB->cd(1);
    can_electronE_xB->cd(1)->SetLogy();
    can_electronE_xB->cd(1)->SetLogz();
    electronE_xB_parallel->Draw("colz");

    can_electronE_xB->cd(2);
    can_electronE_xB->cd(2)->SetLogy();
    can_electronE_xB->cd(2)->SetLogz();
    electronE_xB_antiparallel->Draw("colz");

    can_electronE_xB->Print("plots/2D/plots_2D_electronE_xB.pdf", "pdf");

    // Plots for photonE-Q2 correlation
    TCanvas* can_photonE_Q2 = new TCanvas("can_photonE_Q2");
    can_photonE_Q2->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_photonE_Q2->cd(1);
    can_photonE_Q2->cd(1)->SetLogy();
    can_photonE_Q2->cd(1)->SetLogz();
    photonE_Q2_parallel->Draw("colz");

    can_photonE_Q2->cd(2);
    can_photonE_Q2->cd(2)->SetLogy();
    can_photonE_Q2->cd(2)->SetLogz();
    photonE_Q2_antiparallel->Draw("colz");

    can_photonE_Q2->Print("plots/2D/plots_2D_photonE_Q2.pdf", "pdf");

    // Plots for electronE-Q2 correlation
    TCanvas* can_electronE_Q2 = new TCanvas("can_electronE_Q2");
    can_electronE_Q2->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_electronE_Q2->cd(1);
    can_electronE_Q2->cd(1)->SetLogy();
    can_electronE_Q2->cd(1)->SetLogz();
    electronE_Q2_parallel->Draw("colz");

    can_electronE_Q2->cd(2);
    can_electronE_Q2->cd(2)->SetLogy();
    can_electronE_Q2->cd(2)->SetLogz();
    electronE_Q2_antiparallel->Draw("colz");

    can_electronE_Q2->Print("plots/2D/plots_2D_electronE_Q2.pdf", "pdf");

    // Plots for photonE-t correlation
    TCanvas* can_photonE_t = new TCanvas("can_photonE_t");
    can_photonE_t->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_photonE_t->cd(1);
    can_photonE_t->cd(1)->SetLogy();
    can_photonE_t->cd(1)->SetLogz();
    photonE_t_parallel->Draw("colz");

    can_photonE_t->cd(2);
    can_photonE_t->cd(2)->SetLogy();
    can_photonE_t->cd(2)->SetLogz();
    photonE_t_antiparallel->Draw("colz");

    can_photonE_t->Print("plots/2D/plots_2D_photonE_t.pdf", "pdf");

    // Plots for electronE-t correlation
    TCanvas* can_electronE_t = new TCanvas("can_electronE_t");
    can_electronE_t->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_electronE_t->cd(1);
    can_electronE_t->cd(1)->SetLogy();
    can_electronE_t->cd(1)->SetLogz();
    electronE_t_parallel->Draw("colz");

    can_electronE_t->cd(2);
    can_electronE_t->cd(2)->SetLogy();
    can_electronE_t->cd(2)->SetLogz();
    electronE_t_antiparallel->Draw("colz");

    can_electronE_t->Print("plots/2D/plots_2D_electronE_t.pdf", "pdf");

    // Plots for xB-Q2 correlation
    TCanvas* can_xB_Q2 = new TCanvas("can_xB_Q2");
    can_xB_Q2->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_xB_Q2->cd(1);
    can_xB_Q2->cd(1)->SetLogx();
    can_xB_Q2->cd(1)->SetLogy();
    can_xB_Q2->cd(1)->SetLogz();
    xBQ2_parallel->Draw("colz");
    xBQ2_parallel->GetXaxis()->SetTitle("xB");
    xBQ2_parallel->GetYaxis()->SetTitle("Q^{2}");
    xBQ2_parallel->GetYaxis()->CenterTitle(true);


    can_xB_Q2->cd(2);
    can_xB_Q2->cd(2)->SetLogx();
    can_xB_Q2->cd(2)->SetLogy();
    can_xB_Q2->cd(2)->SetLogz();
    xBQ2_antiparallel->Draw("colz");
    xBQ2_antiparallel->GetXaxis()->SetTitle("xB");
    xBQ2_antiparallel->GetYaxis()->SetTitle("Q^{2}");
    xBQ2_antiparallel->GetYaxis()->CenterTitle(true);

    can_xB_Q2->Print("plots/2D/plots_2D_xB_Q2.pdf", "pdf");

    // Plots for xB-t correlation
    TCanvas* can_xB_t = new TCanvas("can_xB_t");
    can_xB_t->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_xB_t->cd(1);
    can_xB_t->cd(1)->SetLogx();
    can_xB_t->cd(1)->SetLogy();
    can_xB_t->cd(1)->SetLogz();
    xB_t_parallel->Draw("colz");

    can_xB_t->cd(2);
    can_xB_t->cd(2)->SetLogx();
    can_xB_t->cd(2)->SetLogy();
    can_xB_t->cd(2)->SetLogz();
    xB_t_antiparallel->Draw("colz");

    can_xB_t->Print("plots/2D/plots_2D_xB_t.pdf", "pdf");

    // Plots for xB-phi correlation
    TCanvas* can_xB_phi = new TCanvas("can_xB_phi");
    can_xB_phi->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_xB_phi->cd(1);
    can_xB_phi->cd(1)->SetLogx();
    can_xB_phi->cd(1)->SetLogy();
    can_xB_phi->cd(1)->SetLogz();
    xB_phi_parallel->Draw("colz");

    can_xB_phi->cd(2);
    can_xB_phi->cd(2)->SetLogx();
    can_xB_phi->cd(2)->SetLogy();
    can_xB_phi->cd(2)->SetLogz();
    xB_phi_antiparallel->Draw("colz");

    can_xB_phi->Print("plots/2D/plots_2D_xB_phi.pdf", "pdf");

    // Plots for Q2-t correlation
    TCanvas* can_Q2_t = new TCanvas("can_Q2_t");
    can_Q2_t->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_Q2_t->cd(1);
    can_Q2_t->cd(1)->SetLogx();
    can_Q2_t->cd(1)->SetLogy();
    can_Q2_t->cd(1)->SetLogz();
    Q2_t_parallel->Draw("colz");

    can_Q2_t->cd(2);
    can_Q2_t->cd(2)->SetLogx();
    can_Q2_t->cd(2)->SetLogy();
    can_Q2_t->cd(2)->SetLogz();
    Q2_t_antiparallel->Draw("colz");

    can_Q2_t->Print("plots/2D/plots_2D_Q2_t.pdf", "pdf");

    // Plots for Q2-phi correlation
    TCanvas* can_Q2_phi = new TCanvas("can_Q2_phi");
    can_Q2_phi->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_Q2_phi->cd(1);
    can_Q2_phi->cd(1)->SetLogx();
    can_Q2_phi->cd(1)->SetLogy();
    can_Q2_phi->cd(1)->SetLogz();
    Q2_phi_parallel->Draw("colz");

    can_Q2_phi->cd(2);
    can_Q2_phi->cd(2)->SetLogx();
    can_Q2_phi->cd(2)->SetLogy();
    can_Q2_phi->cd(2)->SetLogz();
    Q2_phi_antiparallel->Draw("colz");

    can_Q2_phi->Print("plots/2D/plots_2D_Q2_phi.pdf", "pdf");

    // Plots for t-phi correlation
    TCanvas* can_t_phi = new TCanvas("can_t_phi");
    can_t_phi->Divide(1, 2);

    gStyle->SetOptStat(0);

    can_t_phi->cd(1);
    can_t_phi->cd(1)->SetLogx();
    can_t_phi->cd(1)->SetLogy();
    can_t_phi->cd(1)->SetLogz();
    t_phi_parallel->Draw("colz");

    can_t_phi->cd(2);
    can_t_phi->cd(2)->SetLogx();
    can_t_phi->cd(2)->SetLogy();
    can_t_phi->cd(2)->SetLogz();
    t_phi_antiparallel->Draw("colz");

    can_t_phi->Print("plots/2D/plots_2D_t_phi.pdf", "pdf");

    // Plots for ALU asymmetry
    for(k=1;k<=(Q2BinNumber-1); k++){

        for(l=1;l<=(xBjBinNumber-1); l++) {

            TCanvas* can_all = new TCanvas("can_all");
            can_all->Divide(4, 3);

            for(m=1;m<=(tBinNumber-1); m++) {

                can_all->cd(m);

                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1] = h_phi_parallel[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->GetAsymmetry(h_phi_antiparallel[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]);


                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->SetMinimum(-1.0);
                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->Draw("E");
                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->GetXaxis()->SetTitle("phi");
                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->GetYaxis()->SetTitle("A_LU");
                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->GetYaxis()->CenterTitle(true);
                h_ALU[(k-1) * (xBjBinNumber-1) * (tBinNumber-1) + (l-1) * (tBinNumber-1) + m-1]->Scale(limunosity);


            }

            // Print the Asymmetries for all |t| bins at a given Q2 and xBj bins
            can_all->Print(TString::Format("plots/ALU/plots_all_%d_%d.pdf", k,l), "pdf");

        }


    }


    //return
    return 0;
}
















