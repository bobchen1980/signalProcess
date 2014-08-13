#include "mathcal.h"
# define PI 3.1415926535897932384626433832795
#include <COMPLEX>

using namespace std;


/*---------------------------------------------------------------------------------------------------*/

/*FFT*/
/*----��2���ٸ���Ҷ�任---*/

void FFT(float xin[],int N)
{

	int f,m,LH,nm,i,k,j,L;		//�������
	float p , ps ;
	int le,B,ip;
	complex<float>  w,t;
	//---------------------
	LH=N/2;						//���λֵ��Ȩֵ
	f=N;						//fΪ�м����
	for(m=1;(f=f/2)!=1;m++){;}	//��õ���Ϊ2�ļ����ݣ�������m
	nm=N-2;						//��λ��ʱ����0�͵�N-1���ز���
	j=N/2;
	/****************************/
	complex<float>   *din = new complex<float> [N];
	for(i=0;i<N;i++)
	{
		din[i].real(xin[i]);
		din[i].imag(0);
	}							//��������ʵ��ת���ɸ���


	/*��λ���ַ����*/
	for(i=1;i<=nm;i++)
	{
		if(i<j) {t=din[j];din[j]=din[i];din[i]=t;}
		
		k=LH;
		while(j>=k) {j=j-k;k=k/2;}
		j=j+k;
	}
	/***************************/
	
	/*----------FFT����--------*/

	for(L=1;L<=m;L++)								//����Ȧ����Ȧ��ѭ��������
	{  
		le=pow(2,L);					//leΪ2��L�η�
		B=le/2;							//�������ڵ��࣬��С����
	
		for(j=0;j<=B-1;j++)							//�ڶ�Ȧ����Ȧ��ѭ��������ͬ��ת����
		{
		   p=pow(2,m-L)*j;
		   ps=2*PI/N*p;
		   w.real(cos(ps));
		   w.imag(-sin(ps));					//��Ȩϵ������
		   for(i=j;i<=N-1;i=i+le)					//����Ȧ����Ȧ��ѭ����������
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
	//�ͷ��ڴ�
	delete []din;
} 
/*-------------------------------------FFT-----------------------------------------------*/

/************************************************************************/
/* ������Inte_cal                                                       */
/* ���������                                                           */
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
	av = 10*log10(av)-14;//����14db
	aweight.Close();
	return av;

}
/*-----------------------------------A_weight--------------------------------------------*/

/*octave_anal*/
void octave_anal(float xin[],float yout[9], int N, int fs)
{
	int fscale[10]
		= {22, 45, 90, 180, 355, 710, 1400, 2800, 5600, 11200};//Ƶ�ʽ��ߵ�

	/*float octave[9]={0}; //������Ƶ����ѹֵ*/
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

	float trioctave[26]={0}; //������Ƶ����ѹֵ
	
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


//------ƽ��ֵ--------//
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
//------ƽ����-------//
	for (i=0; i<counts; i++)
	{
		adev = adev + fabs(xin[i]-ave);
	}
	adev = adev/counts;
//-------����-------//
	for (i=0; i<counts; i++)
	{
		var = var + (xin[i]-ave)*(xin[i]-ave);
	}
	var = var/(counts-1);
//------��׼��------//
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

