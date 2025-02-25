const short magics[]={2,8,20,28,50,82,126};
const short nlines=sizeof(magics)/sizeof(magics[0]);
void initGrs(TMultiGraph *mgr,TGraph *gr2n,TGraph *gr2p,TGraph *grnp,TLegend *lgd){
	gr2n->SetMarkerSize(2);
	gr2n->SetLineWidth(2);
	gr2n->SetMarkerStyle(22);
	gr2n->SetMarkerColor(4);
	gr2n->SetLineColor(4);
	gr2n->SetLineStyle(9);
	gr2p->SetMarkerSize(2);
	gr2p->SetLineWidth(2);
	gr2p->SetMarkerStyle(23);
	gr2p->SetMarkerColor(2);
	gr2p->SetLineColor(2);
	gr2p->SetLineStyle(9);
	grnp->SetMarkerSize(2);
	grnp->SetLineWidth(2);
	grnp->SetMarkerStyle(20);
	grnp->SetMarkerColor(1);
	grnp->SetLineColor(1);
	grnp->SetLineStyle(9);

	lgd->AddEntry(gr2n,"S_{2n}");
	lgd->AddEntry(gr2p,"S_{2p}");
	lgd->AddEntry(grnp,"S_{np}");
	mgr->Add(gr2n);
	mgr->Add(gr2p);
	mgr->Add(grnp);

}

void setLines(TGraph *grn,TGraph *grp,TLine *ls[nlines],float &xLow,float &xHigh,float &yLow,float &yHigh){
	for(short i=0;i<grn->GetN();i++){
		float x=grn->GetPointX(i);
		float y=grn->GetPointY(i);
		if(i==0){xLow=x;xHigh=x;yLow=y;yHigh=y;}
		else{
			if(xLow>x) xLow=x;
			if(xHigh<x) xHigh=x;
			if(yLow>y) yLow=y;
			if(yHigh<y) yHigh=y;
		}
	}
	for(short i=0;i<grp->GetN();i++){
		float x=grp->GetPointX(i);
		float y=grp->GetPointY(i);
		if(xLow>x) xLow=x;
		if(xHigh<x) xHigh=x;
		if(yLow>y) yLow=y;
		if(yHigh<y) yHigh=y;
	}
	yLow=yLow-(yHigh-yLow)*0.03;
	yHigh=yHigh+(yHigh-yLow)*0.03;
	for(short i=0;i<nlines;i++){
		ls[i]=new TLine(magics[i],yLow,magics[i],yHigh);
		ls[i]->SetLineWidth(3);
		ls[i]->SetLineStyle(2);
		ls[i]->SetLineColor(1);
	}
}
