//calculating nucleon-pair gaps
/*bool nucleiRuleOut(const short &z,const short &n){
	const short magics[]={8,14,28,50,82,126};//z, n
	for(short i=0;i<sizeof(magics)/sizeof(magics[0]);i++)
		if(z==magics[i]) return false;
		else if(n==magics[i]) return false;
	return true;

}*/
bool calcNPgap1993(const short &a,const short &z,double &npDelta){//calc np pair gap
	const short n=a-z;
	if(z<0 || n<0) {cout<<"wrong z="<<z<<" & n="<<n<<" input!\n";return false;}

	double mnp0s[9]={0.};
	if(getNucleusMass1993(a+1,z,mnp0s[0]) &&
		getNucleusMass1993(a-1,z,mnp0s[1]) &&
		getNucleusMass1993(a-1,z-1,mnp0s[2]) &&
		getNucleusMass1993(a+1,z+1,mnp0s[3]) &&
		getNucleusMass1993(a,z,mnp0s[4]) &&
		getNucleusMass1993(a+2,z+1,mnp0s[5]) &&
		getNucleusMass1993(a,z-1,mnp0s[6]) &&
		getNucleusMass1993(a,z+1,mnp0s[7]) &&
		getNucleusMass1993(a-2,z-1,mnp0s[8])) NULL;
	else return false;//else the execution is failed

	double tmp=(2.*(mnp0s[0]+mnp0s[1]+mnp0s[2]+mnp0s[3]) -4.*mnp0s[4]- (mnp0s[5]+mnp0s[6]+mnp0s[7]+mnp0s[8])) /4.;//MeV
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		npDelta=tmp;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		npDelta=-tmp;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		npDelta=-tmp;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		npDelta=tmp;//MeV
	else return false;
	return true;
}
bool calcNPgap2020(const short &a,const short &z,double &npDelta){//calc np pair gap
	const short n=a-z;
	if(z<0 || n<0) {cout<<"wrong z="<<z<<" & n="<<n<<" input!\n";return false;}

	double mnp0s[9]={0.};
	if(getNucleusMass2020(a+1,z,mnp0s[0]) &&
		getNucleusMass2020(a-1,z,mnp0s[1]) &&
		getNucleusMass2020(a-1,z-1,mnp0s[2]) &&
		getNucleusMass2020(a+1,z+1,mnp0s[3]) &&
		getNucleusMass2020(a,z,mnp0s[4]) &&
		getNucleusMass2020(a+2,z+1,mnp0s[5]) &&
		getNucleusMass2020(a,z-1,mnp0s[6]) &&
		getNucleusMass2020(a,z+1,mnp0s[7]) &&
		getNucleusMass2020(a-2,z-1,mnp0s[8])) NULL;
	else return false;//else the execution is failed

	double tmp=(2.*(mnp0s[0]+mnp0s[1]+mnp0s[2]+mnp0s[3]) -4.*mnp0s[4]- (mnp0s[5]+mnp0s[6]+mnp0s[7]+mnp0s[8])) /4.;//MeV
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		npDelta=tmp;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		npDelta=-tmp;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		npDelta=-tmp;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		npDelta=tmp;//MeV
	else return false;
	return true;
}
bool calcNgap1993(const short &a,const short &z,double &nnGap,const double &npDelta){//calc nn pair gap
	const short n=a-z;
	if(z<0 || n<0) {cout<<"wrong z="<<z<<" & n="<<n<<" input!\n";return false;}

	double mz0s[5]={0.};//z is constant for these nuclei
	if(getNucleusMass1993(a-2,z,mz0s[0]) && //mass has to be obtained
		getNucleusMass1993(a-1,z,mz0s[1]) && 
		getNucleusMass1993(a,z,mz0s[2]) && 
		getNucleusMass1993(a+1,z,mz0s[3]) && 
		getNucleusMass1993(a+2,z,mz0s[4])) NULL;
	else return false;//else the execution is failed

	double tmp= -(mz0s[0]-4.*mz0s[1]+6.*mz0s[2]-4.*mz0s[3]+mz0s[4]) /8.;//MeV
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		nnGap=tmp;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		nnGap=-tmp;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		nnGap=tmp;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		nnGap=-tmp+ npDelta;//MeV
	else return false;
	return true;
}
bool calcNgap2020(const short &a,const short &z,double &nnGap,const double &npDelta){//calc nn pair gap
	const short n=a-z;
	if(z<0 || n<0) {cout<<"wrong z="<<z<<" & n="<<n<<" input!\n";return false;}

	double mz0s[5]={0.};//z is constant for these nuclei
	if(getNucleusMass2020(a-2,z,mz0s[0]) && //mass has to be obtained
		getNucleusMass2020(a-1,z,mz0s[1]) && 
		getNucleusMass2020(a,z,mz0s[2]) && 
		getNucleusMass2020(a+1,z,mz0s[3]) && 
		getNucleusMass2020(a+2,z,mz0s[4])) NULL;
	else return false;//else the execution is failed

	double tmp= -(mz0s[0]-4.*mz0s[1]+6.*mz0s[2]-4.*mz0s[3]+mz0s[4]) /8.;//MeV
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		nnGap=tmp;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		nnGap=-tmp;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		nnGap=tmp;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		nnGap=-tmp+ npDelta;//MeV
	else return false;
	return true;
}
bool calcPgap1993(const short &a,const short &z,double &ppGap,const double &npDelta){//calc pp pair gap
	const short n=a-z;
	if(z<0 || n<0) {cout<<"wrong z="<<z<<" & n="<<n<<" input!\n";return false;}

	double mn0s[5]={0.};//n is constant for these nuclei
	if(getNucleusMass1993(a-2,z-2,mn0s[0]) &&
		getNucleusMass1993(a-1,z-1,mn0s[1]) &&
		getNucleusMass1993(a,z,mn0s[2]) &&
		getNucleusMass1993(a+1,z+1,mn0s[3]) &&
		getNucleusMass1993(a+2,z+2,mn0s[4])) NULL;
	else return false;//else the execution is failed

	double tmp= -(mn0s[0]-4.*mn0s[1]+6.*mn0s[2]-4.*mn0s[3]+mn0s[4]) /8.;//MeV
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		ppGap=tmp;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		ppGap=tmp +npDelta;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		ppGap=-tmp;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		ppGap=-tmp +npDelta;//MeV
	else return false;
	return true;
}
bool calcPgap2020(const short &a,const short &z,double &ppGap,const double &npDelta){//calc pp pair gap
	const short n=a-z;
	if(z<0 || n<0) {cout<<"wrong z="<<z<<" & n="<<n<<" input!\n";return false;}

	double mn0s[5]={0.};//n is constant for these nuclei
	if(getNucleusMass2020(a-2,z-2,mn0s[0]) &&
		getNucleusMass2020(a-1,z-1,mn0s[1]) &&
		getNucleusMass2020(a,z,mn0s[2]) &&
		getNucleusMass2020(a+1,z+1,mn0s[3]) &&
		getNucleusMass2020(a+2,z+2,mn0s[4])) NULL;
	else return false;//else the execution is failed

	double tmp= -(mn0s[0]-4.*mn0s[1]+6.*mn0s[2]-4.*mn0s[3]+mn0s[4]) /8.;//MeV
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		ppGap=tmp;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		ppGap=tmp +npDelta;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		ppGap=-tmp;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		ppGap=-tmp +npDelta;//MeV
	else return false;
	return true;
}

const short magics[]={2,8,20,28,50,82,126};
const short nlines=sizeof(magics)/sizeof(magics[0]);
void setLines(TGraph *grs[nygrs0],TLine *ls[nlines],float &xLow,float &xHigh,const float &yLow,const float &yHigh){
	for(short i=0;i<nygrs0;i++){
		for(short j=0;j<grs[i]->GetN();j++){
			float x=grs[i]->GetPointX(j);
			float y=grs[i]->GetPointY(j);
			if(i==0 && j==0){xLow=x;xHigh=x;}
			else{
				if(xLow>x) xLow=x;
				if(xHigh<x) xHigh=x;
			}
		}
	}
	xLow=xLow-(xHigh-xLow)*0.03;
	xHigh=xHigh+(xHigh-xLow)*0.03;
	for(short i=0;i<nlines;i++){
		ls[i]=new TLine(magics[i],yLow,magics[i],yHigh);
		ls[i]->SetLineWidth(1);
		ls[i]->SetLineStyle(9);
		ls[i]->SetLineColor(1);
	}
}
void setLines(TGraph *grs[nygrs0],TLine *ls[nlines],float &xLow,float &xHigh,float &yLow,float &yHigh){
	for(short i=0;i<nygrs0;i++){
		for(short j=0;j<grs[i]->GetN();j++){
			float x=grs[i]->GetPointX(j);
			float y=grs[i]->GetPointY(j);
			if(i==0 && j==0){xLow=x;xHigh=x;yLow=y;yHigh=y;}
			else{
				if(xLow>x) xLow=x;
				if(xHigh<x) xHigh=x;
				if(yLow>y) yLow=y;
				if(yHigh<y) yHigh=y;
			}
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
void setLines(TGraph *grs[nxgrs0][nygrs0],short nxgrs,short nygrs,TLine *ls[nlines],float &xLow,float &xHigh,float &yLow,float &yHigh){
	for(short i=0;i<nxgrs;i++){
		for(short j=0;j<nygrs;j++){
			for(short k=0;k<grs[i][j]->GetN();k++){
				float x=grs[i][j]->GetPointX(k);
				float y=grs[i][j]->GetPointY(k);
				if(i==0 && j==0 && k==0){xLow=x;xHigh=x;yLow=y;yHigh=y;}
				else{
					if(xLow>x) xLow=x;
					if(xHigh<x) xHigh=x;
					if(yLow>y) yLow=y;
					if(yHigh<y) yHigh=y;
				}
			}
		}
	}
	xLow=xLow-(xHigh-xLow)*0.03;
	xHigh=xHigh+(xHigh-xLow)*0.03;
	yLow=yLow-(yHigh-yLow)*0.02;
	yHigh=yHigh+(yHigh-yLow)*0.03;
	for(short i=0;i<nlines;i++){
		ls[i]=new TLine(magics[i],yLow,magics[i],yHigh);
		ls[i]->SetLineWidth(3);
		ls[i]->SetLineStyle(2);
		ls[i]->SetLineColor(1);
	}
}
