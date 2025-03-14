const short magics[]={2,8,20,28,50,82,126};
const short nlines=sizeof(magics)/sizeof(magics[0]);

void setLines(TGraph *gra,TLine *ls[nlines],float &xLow,float &xHigh,float &yLow,float &yHigh){
	for(short i=0;i<gra->GetN();i++){
		float x=gra->GetPointX(i);
		float y=gra->GetPointY(i);
		if(i==0){xLow=x;xHigh=x;yLow=y;yHigh=y;}
		else{
			if(xLow>x) xLow=x;
			if(xHigh<x) xHigh=x;
			if(yLow>y) yLow=y;
			if(yHigh<y) yHigh=y;
		}
	}
	xLow=xLow-(xHigh-xLow)*0.03;
	xHigh=xHigh+(xHigh-xLow)*0.03;
	yLow=yLow-(yHigh-yLow)*0.03;
	yHigh=yHigh+(yHigh-yLow)*0.03;
	for(short i=0;i<nlines;i++){
		ls[i]=new TLine(magics[i],yLow,magics[i],yHigh);
		ls[i]->SetLineWidth(3);
		ls[i]->SetLineStyle(2);
		ls[i]->SetLineColor(1);
	}
}
