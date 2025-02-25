const short magics[]={2,8,20,28,50,82,126};
const short nlines=sizeof(magics)/sizeof(magics[0]);

void setLines(TGraph *gro,TGraph *gre,TLine *ls[nlines],float &xLow,float &xHigh,float &yLow,float &yHigh){
	for(short i=0;i<gro->GetN();i++){
		float x=gro->GetPointX(i);
		float y=gro->GetPointY(i);
		if(i==0){xLow=x;xHigh=x;yLow=y;yHigh=y;}
		else{
			if(xLow>x) xLow=x;
			if(xHigh<x) xHigh=x;
			if(yLow>y) yLow=y;
			if(yHigh<y) yHigh=y;
		}
	}
	for(short i=0;i<gre->GetN();i++){
		float x=gre->GetPointX(i);
		float y=gre->GetPointY(i);
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
