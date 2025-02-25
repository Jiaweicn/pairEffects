//plot binding energy per nucleon
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "ameMass2020.h"
#include "./shared.h"

void plotBE(short z0=50,short n0=0){//3 kind of nucleon-pair gaps of a series of isotopes/isotones
	short a0=0;
	auto *gro=new TGraph();
	auto *gre=new TGraph();
	auto *mgr=new TMultiGraph();
	if(z0>0 && n0==0){//isotopes
		a0=z0;
		cout<<"calculating binding energy per nucleon for isotopes: Z="<<GREEN<<z0<<RESET<<"\n";
	}
	else if(n0>0 && z0==0){//isotones
		cout<<"calculating binding energy per nucleon for isotones: N="<<GREEN<<n0<<RESET<<"\n";
		a0=n0;
	}
	else{
		cout<<"wrong input!\n";
		exit(1);
	}
	gro->SetMarkerSize(2);	gre->SetMarkerSize(2);
	gro->SetLineWidth(2);	gre->SetLineWidth(2);
	gro->SetMarkerStyle(20);gre->SetMarkerStyle(22);
	gro->SetMarkerColor(1);	gre->SetMarkerColor(2);
	gro->SetLineColor(1);	gre->SetLineColor(2);
	gro->SetLineStyle(9);	gre->SetLineStyle(9);

	double bindingEnergyPerNucleon;
	string elementName;

	for(short a=a0;a<300;a++){//mass number
		if(z0>0 && n0==0){//isotopes
			const short n=a-z0;
			if(getNucleusBindingE2020(a,z0,elementName,bindingEnergyPerNucleon)){
				if(n%2==1) gro->SetPoint(gro->GetN(),n,bindingEnergyPerNucleon);
				else if(n%2==0) gre->SetPoint(gre->GetN(),n,bindingEnergyPerNucleon);
			}
		}
		if(n0>0 && z0==0){//isotones
			const short z=a-n0;
			if(getNucleusBindingE2020(a,z,elementName,bindingEnergyPerNucleon)){
				if(z%2==1) gro->SetPoint(gro->GetN(),z,bindingEnergyPerNucleon);
				else if(z%2==0) gre->SetPoint(gre->GetN(),z,bindingEnergyPerNucleon);
			}
		}
	}
	auto *c=new TCanvas("cgap","cgap",1200,800);

	mgr->Add(gro);
	mgr->Add(gre);
	mgr->Draw("apl");
	if(z0>0 && n0==0){
		mgr->SetTitle(Form("Binding Energy Per Nucleon of %s(Z=%d)",elementName.c_str(),z0));
		mgr->GetXaxis()->SetTitle("Neutron number N");
	}
	else if(n0>0 && z0==0){
		mgr->SetTitle(Form("Binding Energy Per Nucleon of (N=%d)",n0));
		mgr->GetXaxis()->SetTitle("Proton number Z");
	}
	float xLow,xHigh,yLow,yHigh;
	TLine *ls[nlines];
	setLines(gro,gre,ls,xLow,xHigh,yLow,yHigh);
	for(short i=0;i<nlines;i++){
		float x=ls[i]->GetX1();
		if(x>xLow && x<xHigh) ls[i]->Draw();
	}
	mgr->GetYaxis()->SetTitle("BE/A(MeV)");
}
