//plot two particles(neutron/proton) separation energy for an isotope chain
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "../bindingEnergy/ameMass2020.h"
#include "./shared.h"

void plotIsotopes(short z0=50){
	cout<<"calculating two-particles separation energy for isotopes: Z="<<GREEN<<z0<<RESET<<"\n";
	auto *gr2n=new TGraph();
	auto *gr2p=new TGraph();
	auto *grnp=new TGraph();

	auto *lgd=new TLegend(0.45,0.72,0.55,0.9);
	auto *mgr=new TMultiGraph();
	TLine *ls[nlines];
	initGrs(mgr,gr2n,gr2p,grnp,lgd);
	double bindingEnergyPerNucleon0, bindingEnergyPerNucleon1;
	string elementName,elementName0;
	bool nameGot=false;

	for(short a=z0;a<300;a++){//mass number
		const short n=a-z0;
		if(getNucleusBindingE2020(a,z0,elementName0,bindingEnergyPerNucleon0)//nn
			&&getNucleusBindingE2020(a-2,z0,elementName0,bindingEnergyPerNucleon1))
			gr2n->SetPoint(gr2n->GetN(),n,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-2));

		if(getNucleusBindingE2020(a,z0,elementName,bindingEnergyPerNucleon0)//pp
			&&getNucleusBindingE2020(a-2,z0-2,elementName,bindingEnergyPerNucleon1))
			gr2p->SetPoint(gr2p->GetN(),n,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-2));
		if(getNucleusBindingE2020(a,z0,elementName,bindingEnergyPerNucleon0)//np
			&&getNucleusBindingE2020(a-2,z0-1,elementName,bindingEnergyPerNucleon1))
			grnp->SetPoint(grnp->GetN(),n,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-2));
	}
	auto *c=new TCanvas("cgap","cgap",1200,800);

	mgr->Draw("apl");
	mgr->SetTitle(Form("Two-particles Separation Energy of %s(Z=%d)",elementName0.c_str(),z0));
	mgr->GetYaxis()->SetTitle("S_{2n}/S_{2p}/S_{np}(MeV)");
	mgr->GetXaxis()->SetTitle("Neutron number N");
	float xLow,xHigh,yLow,yHigh;
	setLines(gr2n,gr2p,ls,xLow,xHigh,yLow,yHigh);
	for(short i=0;i<nlines;i++){
		float x=ls[i]->GetX1();
		if(x>xLow && x<xHigh) ls[i]->Draw();
	}
	mgr->Draw("p");
	lgd->Draw();
}
