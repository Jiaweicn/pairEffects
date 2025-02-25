//plot one particle(neutron/proton) separation energy for an isotone chain
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "../bindingEnergy/ameMass2020.h"
#include "./shared.h"

void plotIsotones(short n0=50){
	cout<<"calculating separation energy for isotopes: N="<<GREEN<<n0<<RESET<<"\n";
	auto *grn=new TGraph();
	auto *grp=new TGraph();

	auto *lgd=new TLegend(0.45,0.75,0.55,0.9);
	auto *mgr=new TMultiGraph();
	TLine *ls[nlines];
	initGrs(mgr,grn,grp,lgd);
	double bindingEnergyPerNucleon0, bindingEnergyPerNucleon1;
	string elementName;

	for(short a=n0;a<300;a++){//mass number
		const short z=a-n0;
		if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon0)
			&&getNucleusBindingE2020(a-1,z,elementName,bindingEnergyPerNucleon1))
			grn->SetPoint(grn->GetN(),z,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-1));
		if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon0)
			&&getNucleusBindingE2020(a-1,z-1,elementName,bindingEnergyPerNucleon1))
			grp->SetPoint(grp->GetN(),z,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-1));
	}
	auto *c=new TCanvas("cgap","cgap",1200,800);

	mgr->Draw("apl");
	mgr->SetTitle(Form("One-particle Separation Energy of (N=%d)",n0));
	mgr->GetYaxis()->SetTitle("S_{n}/S_{p}(MeV)");
	mgr->GetXaxis()->SetTitle("Proton number Z");
	float xLow,xHigh,yLow,yHigh;
	setLines(grn,grp,ls,xLow,xHigh,yLow,yHigh);
	for(short i=0;i<nlines;i++){
		float x=ls[i]->GetX1();
		if(x>xLow && x<xHigh) ls[i]->Draw();
	}
	mgr->Draw("p");
	lgd->Draw();
}
