const short magics[]={2,8,20,28,50,82,126};
const short nlines=sizeof(magics)/sizeof(magics[0]);
void initGrs(TMultiGraph *mgr,TGraph *grn,TGraph *grp,TLegend *lgd){
	grn->SetMarkerSize(2);
	grn->SetLineWidth(2);
	grn->SetMarkerStyle(20);
	grn->SetMarkerColor(4);
	grn->SetLineColor(4);
	grn->SetLineStyle(9);
	grp->SetMarkerSize(2);
	grp->SetLineWidth(2);
	grp->SetMarkerStyle(22);
	grp->SetMarkerColor(2);
	grp->SetLineColor(2);
	grp->SetLineStyle(9);

	lgd->AddEntry(grn,"S_{n}");
	lgd->AddEntry(grp,"S_{p}");
	mgr->Add(grn);
	mgr->Add(grp);

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
