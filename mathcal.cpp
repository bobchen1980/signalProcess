#include "mathcal.h"
# define PI 3.1415926535897932384626433832795
#include <COMPLEX>

using namespace std;


/*---------------------------------------------------------------------------------------------------*/

/*FFT*/
/*----基2快速傅立叶变换---*/

void FFT(float xin[],int N)
{

	int f,m,LH,nm,i,k,j,L;		//定义变量
	float p , ps ;
	int le,B,ip;
	complex<float>  w,t;
	//---------------------
	LH=N/2;						//最高位值的权值
	f=N;						//f为中间变量
	for(m=1;(f=f/2)!=1;m++){;}	//求得点数为2的几次幂，即级数m
	nm=N-2;						//求倒位序时，点0和点N-1不必参与
	j=N/2;
	/****************************/
	complex<float>   *din = new complex<float> [N];
	for(i=0;i<N;i++)
	{
		din[i].real(xin[i]);
		din[i].imag(0);
	}							//将进来的实数转化成复数


	/*倒位序变址运算*/
	for(i=1;i<=nm;i++)
	{
		if(i<j) {t=din[j];din[j]=din[i];din[i]=t;}
		
		k=LH;
		while(j>=k) {j=j-k;k=k/2;}
		j=j+k;
	}
	/***************************/
	
	/*----------FFT计算--------*/

	for(L=1;L<=m;L++)								//第三圈（外圈）循环：按级
	{  
		le=pow(2,L);					//le为2的L次方
		B=le/2;							//蝶形两节点间距，由小到大
	
		for(j=0;j<=B-1;j++)							//第二圈（中圈）循环：按相同旋转因子
		{
		   p=pow(2,m-L)*j;
		   ps=2*PI/N*p;
		   w.real(cos(ps));
		   w.imag(-sin(ps));					//加权系数计算
		   for(i=j;i<=N-1;i=i+le)					//第三圈（里圈）循环：按蝶形
			 { 
			  ip=i+B;
			  t = din[ip]*w;
			  din[ip] = din[i] - t;			 
			  din[i]  = din[i] + t;			  
			 }
		}
	}

		xin[0] = sqrt(din[0].real()*din[0].real()+din[0].imag()*din[0].imag())/N;
		if (din[0].real()<0) {
			xin[0]=-xin[0];
		}
	for (int h=1; h<N/2; h++)
	{
		xin[h] = 2*sqrt(din[h].real()*din[h].real()+din[h].imag()*din[h].imag())/N;
	}
	//释放内存
	delete []din;
} 
/*-------------------------------------FFT-----------------------------------------------*/

/************************************************************************/
/* 函数：Inte_cal                                                       */
/* 梯形求积法                                                           */
/* X(N) = X(N) + (Ts/2)*[x(N-1)+x(N)]									*/
/************************************************************************/
void Inte_cal(float xin[], int N, float f)
{
	float xout[4096];
	xout[0] = 0;
	for (int k=0; k<N-1; k++)
	{
		xout[k+1] = xout[k] + (1/(2*f))*(xin[k+1]+xin[k]);
	}
	/*memcpy(xin, xout, N);*/
	for ( k=0; k<N; k++)
	{
		xin[k] = xout[k];
	}

	return;
}
/*-----------------------------------Inte_cal--------------------------------------------*/

/*audioline*/
float audioline(float xin[],int N, int fs)
{
	float p0 = (float)2e-5;
	for (int c=0;c<N;c++)
	{
		xin[c] = pow((xin[c]/p0),2);
	}
	Inte_cal(xin,N,4096);
	float leq;
	leq = (float)10*log10(xin[N-1]/(N/fs));

	return leq;
	
}
/*-----------------------------------audioline------------------------------------------*/

/*A_weight*/
float A_weight(float xin[],int N,int fs)
{	
//-----------------------------
	char* lpszFileName;
	CFile aweight;
	CFileException fileException;
	//-----------------------------
	lpszFileName = "C:\\Program Files\\Evndaq\\AWeightdata.dat";/**/
	aweight.Open(lpszFileName,CFile::modeRead | CFile::typeBinary, &fileException);
//------------------------------------------
	FFT(xin,N);
	float p00 = (float)2e-5;
	int fp;
	float aw;
	float av = 0;
	float tem;
	for (int d = 0; d < N/2; d++)
	{	
		if (xin[d] == 0)
		{
			xin[d] = (float)2e-5;
		}
		tem = xin[d];
		xin[d]=20*log10(tem/p00);
		fp = d*fs/N;
		aweight.Seek(2*fp*sizeof(float),CFile::begin);
		aweight.Read((void*)&aw,sizeof(float));
		xin[d] =xin[d] + aw;
		av += pow(10,0.1*xin[d]);

	}
	av = 10*log10(av)-14;//修正14db
	aweight.Close();
	return av;

}
/*-----------------------------------A_weight--------------------------------------------*/

/*octave_anal*/
void octave_anal(float xin[],float yout[9], int N, int fs)
{
	int fscale[10]
		= {22, 45, 90, 180, 355, 710, 1400, 2800, 5600, 11200};//频率界线点

	/*float octave[9]={0}; //各中心频率声压值*/
	for (int i=0;i<9;i++)
	{
		yout[i]=0;
	}
//-------------------------------------//	
	FFT(xin,N);
	xin[0]=0;
	for (int q=0; q<N/2;q++)
	{
		for (int fo=0;fo<9;fo++)
		{
			if (((q*fs/N)>fscale[fo]) && ((q*fs/N)<=fscale[fo+1]))
			{
				/*octave[fo] = octave[fo] + (xin[q])*(xin[q]);*/
				yout[fo] = yout[fo] + (xin[q])*(xin[q]);
			}
		}
 
	}
	for (int fo1=0;fo1<9;fo1++)
	{
		if (yout[fo1]==0)
		{
			yout[fo1] = (float)4e-10;
		}
		/*octave[fo1] = 20*log10(sqrt(octave[fo1])/(2e-5));*/
		yout[fo1] = 20*log10(sqrt(yout[fo1])/(2e-5));
	}
}

/*----------------------------------octave_anal------------------------------------------*/


/*one_third_octave*/
void one_third_octave(float xin[], int N, int fs)
{
	int trifs[27]
		= { 45, 56, 71, 90, 112, 140, 180, 224, 280, 355, 450, 560, 710, 900, 1120, 1400,
			1800, 2240, 2800, 3550, 4500, 5600, 7100, 9000, 11200, 14000, 18000};

	float trioctave[26]={0}; //各中心频率声压值
	
	FFT(xin,N);
	for (int q; q<N/2;q++)
	{
		for (int fo=0;fo<9;fo++)
		{
			if (((q*fs/N)>trifs[fo]) && ((q*fs/N)<=trifs[fo+1]))
			{
				trioctave[fo] = trioctave[fo] + (xin[q])*(xin[q]);
			}
		}
 
	}
	for (int fo1=0;fo1<9;fo1++)
	{
		if (trioctave[fo1]==0)
		{
			trioctave[fo1] =(float) 4e-10;
		}
		trioctave[fo1] = 20*log10(sqrt(trioctave[fo1])/(2e-5));
	}


}

/*-------------------------------one_third_octave----------------------------------------*/

/*statisticcal*/
void statisticcal(float *xin,long counts, double ave, float adev, float sedv,float var,float dmax,
				  float dmin,float absmax)
{
	ave = 0;
	adev = 0;
	sedv = 0;
	var = 0;
	dmax = xin[0];
	dmin =xin[0];
	absmax = xin[0];


//------平均值--------//
	for (long i=0; i<counts;i++)
	{
		ave = ave + xin[i];
		if (dmax<xin[i])
		{
			dmax = xin[i];
		}
		if (dmin>xin[i])
		{
			dmin = xin[i];
		}
		if (absmax<fabs(xin[i]))
		{
			absmax = fabs(xin[i]);
		}
	}
	ave = ave/counts;
//------平均差-------//
	for (i=0; i<counts; i++)
	{
		adev = adev + fabs(xin[i]-ave);
	}
	adev = adev/counts;
//-------方差-------//
	for (i=0; i<counts; i++)
	{
		var = var + (xin[i]-ave)*(xin[i]-ave);
	}
	var = var/(counts-1);
//------标准差------//
	sedv = sqrt(var);


	return;
}
/*---------------------------------statisticcal------------------------------------------*/

/*array_cpy*/
void array_cpy(float* array_target,float* array_sourece,int data_num)
{
	for(int cnt=0; cnt<data_num; cnt++)
	{
		array_target[cnt] = array_sourece[cnt];
	}
}
/*---------------------------------array_cpy---------------------------------------------*/

