#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include <TH1.h>
#include <TROOT.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "DVCSEvent.h"

#include <iostream>

using namespace HepMC3;

//main
int main(int argc, char **argv) {


    double Q2Bins[] = {1.0, 1.77828, 3.16228, 5.62341, 10, 17.7828, 31.6228, 56.2341, 100, 177.828, 316.228, 562.341, 1000.0};

    double xBjBins[] = {0.0001, 0.000158489, 0.000251189, 0.000398107, 0.000630957, 0.001, 0.00158489, 0.00251189, 0.00398107, 0.00630957,
            0.01, 0.0158489, 0.0251189, 0.0398107, 0.0630957, 0.1, 0.158489, 0.251189, 0.398107, 0.630957};

    int i,j,k,l,m;

    //int Q2BinNumber = 13;
    //int xBjBinNumber = 20;

    int Q2BinNumber = sizeof(Q2Bins)/sizeof(Q2Bins[0]);
    int xBjBinNumber = sizeof(xBjBins)/sizeof(xBjBins[0]);


    //histograms

    std::vector<TH1D*> h_t(300, nullptr);

    for(k=0;k<300; k++) {
        h_t[k] = new TH1D(TString::Format("h_t_%d", k), "", 10, 0.2, 1.3);
        h_t[k]-> Sumw2();
    }


    for(m=1; m<=3; m++) {
        ReaderAscii inputFile("10_100_DVCS_parallel/events_DVCS_parallel_" + to_string(m) + ".txt");

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
                for(i=0; i<(Q2BinNumber-1); i++) {
                    for(j=0; j<(xBjBinNumber-1); j++){

                        if((Q2Bins[i] < dvcsEvent.getQ2()) && (dvcsEvent.getQ2() < Q2Bins[i+1])) {
                            if((xBjBins[j] < dvcsEvent.getXb()) && (dvcsEvent.getXb() < xBjBins[j+1])) {
                                h_t[i*(xBjBinNumber-1)+j]->Fill(-1. * dvcsEvent.getT());

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

    for(m=1; m<=3; m++) {
        ReaderAscii inputFile("10_100_DVCS_antiparallel/events_DVCS_antiparallel_" + to_string(m) + ".txt");

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
            for(i=0; i<(Q2BinNumber-1); i++) {
                for(j=0; j<(xBjBinNumber-1); j++){

                    if((Q2Bins[i] < dvcsEvent.getQ2()) && (dvcsEvent.getQ2() < Q2Bins[i+1])) {
                        if((xBjBins[j] < dvcsEvent.getXb()) && (dvcsEvent.getXb() < xBjBins[j+1])) {
                            h_t[i*(xBjBinNumber-1)+j]->Fill(-1. * dvcsEvent.getT());

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

    double limunosity =  0.0;


    for(k=1;k<=(Q2BinNumber-1); k++){

        TCanvas* can_all = new TCanvas("can_all");
        can_all->Divide(5, 4);

        for(l=1;l<=(xBjBinNumber-1); l++) {

            can_all->cd(l);
            can_all->cd(l)->SetLogy();
            h_t[(k-1) * (xBjBinNumber-1) + l-1]->SetMinimum(1);
            h_t[(k-1) * (xBjBinNumber-1) + l-1]->Draw("E");
            h_t[(k-1) * (xBjBinNumber-1) + l-1]->GetXaxis()->SetTitle("|t| (GeV^{2})");
            h_t[(k-1) * (xBjBinNumber-1) + l-1]->GetYaxis()->SetTitle("Number of events");
            h_t[(k-1) * (xBjBinNumber-1) + l-1]->GetYaxis()->CenterTitle(true);

        }

        can_all->Print(TString::Format("plots/plots_all_%d.pdf", k), "pdf");

    }

    //return
    return 0;
}

