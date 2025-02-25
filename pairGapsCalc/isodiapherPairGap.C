//calculation nucleon-pair gap of an isodiapher chain
//using mass to calc gap
//using data from mass table of artemis
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "../constant.h"
#include "ameMass2020.h"
const short nxgrs0=3,nygrs0=4;
#include "shared.h"

void initMGR(const short nMinusZ,TGraph *grs[nxgrs0][nygrs0],TMultiGraph *mgrs[nxgrs0]){
	const short colors[]={4,2,1,6};
	const short markers[]={89,90,91,95};
	const short lines[]={9,5,7,2};
	for(short i=0;i<nxgrs0;i++){//different cantegories of nuclei
		mgrs[i]=new TMultiGraph();
		for(short j=0;j<nygrs0;j++){//different cantegories of nuclei
			grs[i][j]=new TGraph();
			grs[i][j]->SetMarkerSize(2);
			grs[i][j]->SetLineWidth(2);
			grs[i][j]->SetMarkerStyle(markers[j]);
			grs[i][j]->SetLineStyle(lines[j]);
			grs[i][j]->SetMarkerColor(colors[i]);
			grs[i][j]->SetLineColor(colors[i]);
			if(i==0){//nn
				grs[i][j]->SetTitle(Form("neutron-neutron pair gap(N-Z=%d)",nMinusZ));
				grs[i][j]->GetXaxis()->SetTitle("Neutron number N");
				grs[i][j]->GetYaxis()->SetTitle("#Delta_{nn}/MeV");
			}
			else if(i==1){//pp
				grs[i][j]->SetTitle(Form("proton-proton pair gap(N-Z=%d)",nMinusZ));
				grs[i][j]->GetXaxis()->SetTitle("Proton number Z");
				grs[i][j]->GetYaxis()->SetTitle("#Delta_{pp}/MeV");
			}
			else if(i==2){//np
				grs[i][j]->SetTitle(Form("neutron-proton interaction energy(N-Z=%d)",nMinusZ));
				grs[i][j]->GetXaxis()->SetTitle("Mass number A");
				//grs[i][j]->GetXaxis()->SetTitle("Proton number Z");
				grs[i][j]->GetYaxis()->SetTitle("#delta_{np}/MeV");
			}
			//mgrs[i]->Add(grs[i][j]);
		}
	}
}
void plot2(const short &nMinusZ,TGraph *grs[nxgrs0][nygrs0]){
	string pairLgd[nxgrs0]={"#Delta_{nn}","#Delta_{pp}","#delta_{np}"};
	string typeLgd[nygrs0]={"eZ-eN","eZ-oN","oZ-eN","oZ-oN"};
	auto *cgap=new TCanvas("cgap","cgap",900,600);
	auto *lgd=new TLegend(0.65,0.65,0.9,0.9);
	lgd->SetFillStyle(0);
	auto *mgr=new TMultiGraph();
	for(short i=0;i<nxgrs0-1;i++){//different pairs
		for(short j=0;j<nygrs0;j++){//different cantegories of nuclei
			if(grs[i][j]->GetN()>0){
				mgr->Add(grs[i][j]);
				lgd->AddEntry(grs[i][j],(pairLgd[i]+": "+typeLgd[j]).c_str());
			}
		}
	}
	mgr->Draw("apl");
	mgr->SetTitle(Form("nucleon-pair gaps(N-Z=%d)",nMinusZ));
	mgr->GetXaxis()->SetTitle("p/n number");
	mgr->GetYaxis()->SetTitle("#Delta/MeV");
	TLine *ls[nlines];
	float xLow,xHigh,yLow,yHigh;
	setLines(grs,nxgrs0-1,nygrs0,ls,xLow,xHigh,yLow,yHigh);
	for(short i=0;i<nlines;i++){
		float x=ls[i]->GetX1();
		if(x>xLow && x<xHigh) ls[i]->Draw();
	}
	mgr->Draw("p");
	lgd->Draw();
}
void plot3(const short &nMinusZ,TGraph *grs[nxgrs0][nygrs0],TMultiGraph *mgrs[nxgrs0]){
	string typeLgd[nygrs0]={"eZ-eN","eZ-oN","oZ-eN","oZ-oN"};
	TCanvas *cgaps[nxgrs0];
	TLegend *lgds[nxgrs0];
	TLine *ls[nxgrs0][nlines];
	float xLow,xHigh,yLow,yHigh;
	for(short i=0;i<nxgrs0;i++){//different pairs
		cgaps[i]=new TCanvas("","",900,600);
		lgds[i]=new TLegend(0.75,0.75,0.9,0.9);
		lgds[i]->SetFillStyle(0);
		for(short j=0;j<nygrs0;j++){//different cantegories of nuclei
			if(grs[i][j]->GetN()>0){
				mgrs[i]->Add(grs[i][j]);
				lgds[i]->AddEntry(grs[i][j],typeLgd[j].c_str());
			}
		}
		mgrs[i]->Draw("apl");
		mgrs[i]->SetTitle(grs[i][0]->GetTitle());
		mgrs[i]->GetXaxis()->SetTitle(grs[i][0]->GetXaxis()->GetTitle());
		mgrs[i]->GetYaxis()->SetTitle(grs[i][0]->GetYaxis()->GetTitle());
		setLines(grs[i],ls[i],xLow,xHigh,yLow,yHigh);
		for(short j=0;j<nlines;j++){
			float x=ls[i][j]->GetX1();
			if(x>xLow && x<xHigh) ls[i][j]->Draw();
		}
		mgrs[i]->Draw("p");
		lgds[i]->Draw();
	}
}
void isodiapherPairGap(const short nMinusZ=20){//pair gaps of an isodiapher chain
	TGraph *grs[nxgrs0][nygrs0];
	TMultiGraph *mgrs[nxgrs0];
	initMGR(nMinusZ,grs,mgrs);

	for(short z=1;z<120;z++){//charge of nuclei
		const short n=z+nMinusZ;//neutron number
		const short a=z+n;//mass number
		double nnGap,ppGap,npDelta;
		short type;
		if(!nucleiType(z,n,type)) continue;
		//if(calcNPgap1993(a,z,npDelta))//np pair gap, MeV
		if(calcNPgap2020(a,z,npDelta))//np pair gap, MeV
			grs[2][type]->SetPoint(grs[2][type]->GetN(),a,npDelta);
		else continue;
		//if(calcNgap1993(a,z,nnGap,npDelta))//nn pair gap, MeV
		if(calcNgap2020(a,z,nnGap,npDelta))//nn pair gap, MeV
			grs[0][type]->SetPoint(grs[0][type]->GetN(),n,nnGap);
		//if(calcPgap1993(a,z,ppGap,npDelta))//pp pair gap, MeV
		if(calcPgap2020(a,z,ppGap,npDelta))//pp pair gap, MeV
			grs[1][type]->SetPoint(grs[1][type]->GetN(),z,ppGap);
	}
	plot2(nMinusZ,grs);//plot nn & pp gaps in one canvas
	plot3(nMinusZ,grs,mgrs);//plot nn, pp & np gaps in 3 canvases
}
