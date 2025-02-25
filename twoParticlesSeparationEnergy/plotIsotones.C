//plot two particles(neutron/proton) separation energy for an isotone chain
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "../bindingEnergy/ameMass2020.h"
#include "./shared.h"

void plotIsotones(short n0=50){
	cout<<"calculating two-particles separation energy for isotopes: N="<<GREEN<<n0<<RESET<<"\n";
	auto *gr2n=new TGraph();
	auto *gr2p=new TGraph();
	auto *grnp=new TGraph();

	auto *lgd=new TLegend(0.45,0.72,0.55,0.9);
	auto *mgr=new TMultiGraph();
	TLine *ls[nlines];
	initGrs(mgr,gr2n,gr2p,grnp,lgd);
	double bindingEnergyPerNucleon0, bindingEnergyPerNucleon1;
	string elementName;

	for(short a=n0;a<300;a++){//mass number
		const short z=a-n0;
		if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon0)//nn
			&&getNucleusBindingE2020(a-2,z,elementName,bindingEnergyPerNucleon1))
			gr2n->SetPoint(gr2n->GetN(),z,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-2));
		if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon0)//pp
			&&getNucleusBindingE2020(a-2,z-2,elementName,bindingEnergyPerNucleon1))
			gr2p->SetPoint(gr2p->GetN(),z,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-2));
		if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon0)//np
			&&getNucleusBindingE2020(a-2,z-1,elementName,bindingEnergyPerNucleon1))
			grnp->SetPoint(grnp->GetN(),z,bindingEnergyPerNucleon0*a-bindingEnergyPerNucleon1*(a-2));
	}
	auto *c=new TCanvas("cgap","cgap",1200,800);

	mgr->Draw("apl");
	mgr->SetTitle(Form("Two-particles Separation Energy of (N=%d)",n0));
	mgr->GetYaxis()->SetTitle("S_{2n}/S_{2p}/S_{np}(MeV)");
	mgr->GetXaxis()->SetTitle("Proton number Z");
	float xLow,xHigh,yLow,yHigh;
	setLines(gr2n,gr2p,ls,xLow,xHigh,yLow,yHigh);
	for(short i=0;i<nlines;i++){
		float x=ls[i]->GetX1();
		if(x>xLow && x<xHigh) ls[i]->Draw();
	}
	mgr->Draw("p");
	lgd->Draw();
}
