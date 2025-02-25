//plot two particles(neutron/proton) separation energy for an isotone chain
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "../bindingEnergy/ameMass2020.h"
#include "./shared.h"

void plotIsotones(short n0=50){
	cout<<"calculating alpha-particle Q-value for isotopes: N="<<GREEN<<n0<<RESET<<"\n";
	auto *gra=new TGraph();
	gra->SetMarkerSize(2);
	gra->SetLineWidth(2);
	gra->SetMarkerStyle(20);
	gra->SetMarkerColor(4);
	gra->SetLineColor(4);
	gra->SetLineStyle(9);

	TLine *ls[nlines];
	double bindingEnergyPerNucleon0, bindingEnergyPerNucleon1;
	string elementName;

	getNucleusBindingE2020(4,2,elementName,bindingEnergyPerNucleon0);
	double bindingEnergyAlpha=bindingEnergyPerNucleon0*4.;//total binding energy of alpha

	for(short a=n0;a<300;a++){//mass number
		const short z=a-n0;
		if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon0)
			&&getNucleusBindingE2020(a-4,z-2,elementName,bindingEnergyPerNucleon1)){
			double qvalue=bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-4)-bindingEnergyAlpha;//alpha is bound
			//qvalue=-qvalue;
			gra->SetPoint(gra->GetN(),z,qvalue);//alpha is bound
		}
	}
	auto *c=new TCanvas("cgap","cgap",1200,800);

	gra->Draw("apl");
	gra->SetTitle(Form("#alpha-particle Q-value of (N=%d)",n0));
	gra->GetYaxis()->SetTitle("Q_{#alpha}(MeV)");
	gra->GetXaxis()->SetTitle("Proton number Z");
	float xLow,xHigh,yLow,yHigh;
	setLines(gra,ls,xLow,xHigh,yLow,yHigh);
	for(short i=0;i<nlines;i++){
		float x=ls[i]->GetX1();
		if(x>xLow && x<xHigh) ls[i]->Draw();
	}
}
