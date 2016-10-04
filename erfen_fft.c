                                         */
/*****************************************************************************/



/********************************************************************************************/
/*    此函数主要目的:计算10^-20到10^20之间的浮点值的cos值(正负数均可)，在时间、空间足够的情 */
/*    况下可精确到小数点后任意一位，并可在一般配置的电脑上1分钟内精确到100万位              */
/*    的值。                                                                                */
/*    主要数据结构：用int[]数组存储大数，其中int[0]表示正负号，正号用1表                    */
/*    示，符号用-1表示，0表示0；int[1]表示存储数据的radix的幂次	；                          */
/*    int[i]里每一项存储数字的位数由radix决定，每个int型结点里都存有Log10_radix             */
/*    位数，其中radix由radix test这段代码确定，以防止每个int型结点溢出。                    */
/********************************************************************************************/
/*
---- calculation of PI(= 3.14159...) using FFT ----
    by T.Ooura, ver. LG1.1.2-MP1.5a Sep. 2001.

This is a test program to estimate the performance of
the FFT routines: fft*g.c.

Example compilation:
    GNU      : gcc -O6 -ffast-math pi_fft.c fftsg.c -lm -o pi_fftsg
    SUN      : cc -fast -xO5 pi_fft.c fft8g.c -lm -o pi_fft8g
    Microsoft: cl /O2 /G6 pi_fft.c fft4g.c /Fepi_fft4g.exe
    ...
    etc.
*/

/* Please check the following macros before compiling */
#ifndef DBL_ERROR_MARGIN
#define DBL_ERROR_MARGIN 0.3  /* must be < 0.5 */
#endif


#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mem.h>
#include <assert.h>

void mp_load_0(int n, int radix, int out[]);/**/
void mp_load_1(int n, int radix, int out[]);/**/
void mp_copy(int n, int radix, int in[], int out[]);/*将一个数据从in[]里读入复制给out[]，即令out[]=in[]，n表示out[]的长度*/
void mp_round(int n, int radix, int m, int inout[]);/*将数据倒序，例如将123变成321，n表示位数，inout[]为数据的读入和倒序后保存的地址*/
int mp_cmp(int n, int radix, int in1[], int in2[]);/*比较数据in1[]和数据in2[]的大小*/
void mp_add(int n, int radix, int in1[], int in2[], int out[]);/*将数据in1[]和数据in2[]相加存入out[]中，即加法函数*/
void mp_sub(int n, int radix, int in1[], int in2[], int out[]);/*数据in1[]-数据in2[]，存入out[]中，即减法函数*/
void mp_imul(int n, int radix, int in1[], int in2, int out[]);/*数乘：将数据in1[]*一个int型整数in2得到out[]*/
int mp_idiv(int n, int radix, int in1[], int in2, int out[]);/*数除：将数据in1[]除以一个int型整数in2，得到数据out[]*/
void mp_idiv_2(int n, int radix, int in[], int out[]);/*为了提高程序效率节省时间，数除时除数是2时而编的一个函数，即数据in[]/2得到ut[]*/
double mp_mul_radix_test(int n, int radix, int nfft,
        double tmpfft[], int ip[], double w[]);     /*计算radix的值*/
void mp_mul(int n, int radix, int in1[], int in2[], int out[],
        int tmp[], int nfft, double tmp1fft[], double tmp2fft[],
        double tmp3fft[], int ip[], double w[]);     /*乘法运算，in1[]*in2[]，输出到out[]里*/
void mp_squ(int n, int radix, int in[], int out[], int tmp[],
        int nfft, double tmp1fft[], double tmp2fft[],
        int ip[], double w[]);     /*平方运算，out[]=(in[])^2*/
void mp_mulh(int n, int radix, int in1[], int in2[], int out[],
        int nfft, double in1fft[], double outfft[],
        int ip[], double w[]);     /*同样是乘法运算，但速度比 mp_mul函数快，精确度稍微降低，容易溢出，属于快乘，即in1[]*in2[]=out[]*/
void mp_squh(int n, int radix, int in[], int out[],
        int nfft, double inoutfft[], int ip[], double w[]);/*同样是平方运算，但速度比 mp_squ函数快，精确度稍微降低，容易溢出，属于快乘方*/
int mp_inv(int n, int radix, int in[], int out[],
        int tmp1[], int tmp2[], int nfft,
        double tmp1fft[], double tmp2fft[], int ip[], double w[]);/*求倒数，out[]=1/in[]，计算除法时将除数变成倒数，即将除法变成乘法，提升了效率即为了使in1[]/in2[]变成in1[]*(1/in2[])*/
int mp_sqrt(int n, int radix, int in[], int out[],
        int tmp1[], int tmp2[], int nfft,
        double tmp1fft[], double tmp2fft[], int ip[], double w[]);/*开方运算,out[]=(in[])^(0.5)*/
void mp_sprintf(int n, int log10_radix, int in[], char out[]);/*调入int型数组in[]，将最后浮点值的cos的值从字符串out[]中输出*/
void mp_sscanf(int n, int log10_radix, char in[], int out[]);/*调入字符串in[]，从字符串将数据读入到int型数组out[]中*/
void mp_fprintf(int n, int log10_radix, int in[], FILE *fout);/*调入int型数组in[]，将最后浮点值的cos的值写入fout文件中*/

	char nn[2000],kk[2000];  /*全局变量，nn存储输入的浮点值，kk将nn的里输入部分除了小数点均变为0*/
	int jj,log10_radix,flag;      /*全局变量 jj为记录基底为10负多少次方若为整数，jj=0，
	                                 log10_radix为记录每个int型结点里存的位数,flag判断是否为整数*/

/*
	binary的作用是将浮点值里除以基底取出整数，将其化成倒写的二进制，
	如nn[]为1.0010，则先除以0.0001再将其化成倒写的二进制0101100011100100，
	储存在mm中，返回main函数*/
int * binary()
{

	int * mm=(int*)malloc(2000*sizeof(int));/*mm为存储倒写的二进制数据*/
	int ii=0;
	memset(kk,'0',strlen(nn));
	kk[1]='.';
	kk[strlen(nn)+1]='\0';
	while(strcmp(kk,nn)!=0)
	{
		int j=0,k=0,i;
		for(i=2;i<(int)strlen(nn);i++)
		{
			if(nn[i]==nn[i]/2*2)
				k=0;//k是余数
			else k=1;
			nn[i]=(nn[i]-48)/2+j*5+48;
			j=k;
		}
		mm[ii++]=k;

	}
	mm[ii]=2;
	return mm;
}

void check()//使负数变为正数，因为cos（-a）=cos（a）
{
	int i;
	for(i=0;i<=(int)strlen(nn);i++)
	{
		nn[i-1]=nn[i];
	}
}

void check1(int nfft)//检验输入合法
{
	if(nfft<=0)
	{
		printf("输入有误！！\n");
		exit(0);
	}
}

void check2(char *nn)//检验输入合法
{
	int i,j=0;
	for(i=0;i<strlen(nn);i++)
	{
		if(((nn[i]=='+')||(nn[i]=='-'))&&(i==0)) continue;
		if(nn[i]=='.')
		{
			j++;
			continue;
		}
		if(nn[i]>'9'||nn[i]<'0')
			check1(-1);//等于退出程序
	}
	if(j>1)//禁止两个或多个小数点
		check1(-1);
}

////////////////////////////////////////////////
// 二分法计算cos值
 /* m=[(a+b)/2] }
	P(a,b)=P(a,m)Q(m,b)+P(m,b), Q(a,b)=Q(a,m)Q(m,b)
	R(a,a+1)=(2a+1)(2a+2)    Q(a,a+1)=(2a+1)(2a+2)x*x
	P(a,a+1)=(-1)^(a+1)*/
/*
   用二分法不断递归这种求e的方法来求cos10^(-jj)的值。
  cos(x)=1-1/(2!*x(-2))+1/(4!*x(-4))……调用参数中P代表分子，Q代表分母，
  根据公式P(N0,N1)=P(N0,NMid)*Q(NMid,N1)+P(NMid,N1),
  Q(N0,N1)=Q(N0,NMid)*Q(NMid,N1)如此不断通分递归得到基底的cos值
*/
void BinarySplittingCOS(long N0,long N1,int radix, int *P, int *Q,int *ip,double *w)
{
	int *PP, *QQ;
	long NMid; /*NMid=(N0+N1)/2*/
	if (N1-N0 == 1)
	{
		if(N0==(int)(N0/2)*2)
		{
			P[2]=0;
			P[0]=0;
		}
		else
			if(N0==(int)((N0-1)/4)*4+1)
			{
				P[2]=1;
				P[0]=-1;
			}
			else
			{
				P[0]=1;
				P[2]=1;
			}
		P[1] = 0;
		P[3] = 0;
		if(flag!=1)
		{
			Q[0] = 1;
			if(N1 < radix)
			{
				Q[1] = jj/log10_radix;
				Q[2] = N1;
				memset(Q+3,0,jj/log10_radix);
			}
			else
			{
				Q[1] = jj/log10_radix+1;
				Q[2] = N1/radix;
				Q[3] = N1 - Q[2]*radix;
				memset(Q+4,0,jj/log10_radix);
			}
		}
		else
		{
			Q[0] = 1;
			if(N1 < radix)
			{
				Q[1] = 0;
				Q[2] = N1;
				Q[3] = 0;
			}
			else
			{
				Q[1] = 1;
				Q[2] = N1/radix;
				Q[3] = N1 - Q[2]*radix;
			}
		}
		return;
	}
	NMid = (N0+N1) >> 1;
	BinarySplittingCOS(N0, NMid, radix, P, Q, ip, w);
	//开始计算继续运算所需空间
	long nMaxSize;
	nMaxSize = N1-NMid;
	if(nMaxSize >= 10)
		nMaxSize = nMaxSize* (log10(nMaxSize)+(jj/log10_radix+2)+3);//key
	else
		nMaxSize += (jj/log10_radix+2)*(8+jj/log10_radix*2);//10 or 8?8 is not enough.

	PP =(int *)malloc(nMaxSize*sizeof(int));// new int[nMaxSize];
	assert(PP);
	QQ = (int *)malloc(nMaxSize*sizeof(int));
	assert(QQ);
	memset(PP, 0, nMaxSize*sizeof(int));
	memset(QQ, 0, nMaxSize*sizeof(int));
	BinarySplittingCOS(NMid, N1, radix, PP, QQ, ip, w);
	int nfft;
	nfft = 1;
	long n = P[1]>=QQ[1]?P[1]:QQ[1];//比较p与qq谁长
	n += 1;
	while(nfft < n)
	nfft <<= 1;
	nfft <<= 1;
	n = nfft+2;
	double *tmp1fft, *tmp2fft;
	tmp1fft = (double *)malloc((n+2)*sizeof(double));//new double[n+2];
	assert(tmp1fft);
	tmp2fft = (double *)malloc((n+2)*sizeof(double));
	assert(tmp2fft);
	//计算：公式P(a,b)=P(a,m)Q(m,b)+P(m,b), Q(a,b)=Q(a,m)Q(m,b)
	mp_mulh(n, radix, P, QQ, P, nfft, (double*)tmp1fft, (double*)tmp2fft, ip, (double*)w);
	mp_add(P[1]+3, radix, P, PP, P);
	mp_mulh(n, radix, Q, QQ, Q, nfft, (double*)tmp1fft, (double*)tmp2fft, ip, (double*)w);
	if(PP)
		free(PP);
	if(QQ)
		free(QQ);
	if(tmp1fft)
		free(tmp1fft);
	if(tmp2fft)
		free(tmp2fft);
}

void radix_test(int *radix,double *err,int *n,int *nfft,double *d1,int *ip,double * w)
{//通过误差计算来获得radix的最大值，保证不溢出
	log10_radix = 1;
    *radix = 10;
    *err = mp_mul_radix_test(*n, *radix, *nfft, d1, ip, w);
    *err += DBL_EPSILON * ((*n) * (*radix) * (*radix) / 4);
    while (100 * (*err) < DBL_ERROR_MARGIN && (*radix) <= INT_MAX / 20)
	{
        (*err) *= 100;
        log10_radix++;
        (*radix) *= 10;
    }
}

void adjust(int *log2_nfft,int *nfft)
{//调整nfft
	for ((*log2_nfft) = 1; (1 << (*log2_nfft)) < (*nfft); (*log2_nfft)++);
    *nfft = 1 << (*log2_nfft);
}

void result(char *zhizhen)
{
	int ii;
	if(flag==1)//把整数整理为小数
	{
		for(ii=strlen(nn)+2;ii>=2;ii--)
		{
			nn[ii]=nn[ii-2];//后移2位
		}
		nn[0]='0';
		nn[1]='.';
		return;
	}
	if(flag==2)//把整数与小数的共有体整理为小数
	{
		for(ii=strlen(nn)+2;ii>=zhizhen-nn+2;ii--)//如1.2334
		{
			nn[ii]=nn[ii-1];//后移1位
		}
		for(ii=zhizhen-nn+1;ii>=2;ii--)
		{
			nn[ii]=nn[ii-2];//后移2位
		}
		nn[0]='0';
		nn[1]='.';
	}
	for(ii=strlen(nn);ii<jj+2+(zhizhen-nn);ii++)
	{
		nn[ii]='0';
	}
	nn[ii]='\0';
}
/*
    main函数：

   计算cos的值采用了泰勒级数展开的基本思想，但为了提高算法的效率，先计算
   10^-(log10_radix*k)（(log10_radix*k)即为里面的变量jj，为基底，jj>=0）的cos值，
  如计算cos1.00123的100万位的值，则先求出cos0.00000001的值，再根据
  cos(a+b)=cos(a)*cos(b)-sin(a)*sin(b)，sin(a+b)=sin(a)*cos(b)+cos(a)*sin(b)
  以及二倍角公式进行迭代求出所得结果
*/
int main()
{
	char *zhizhen;
    int nfft, log2_nfft, radix, n,ii;
	/*nfft表示计算到的位数，若要精确到相应的位数，需要更大的位数，先求出log2_nfft，再取2^(log2_nfft)*/
    double err, d_time, n_op;/*err是计算误差，d_time是运行时间， n_op是操作次数*/
    int *a, *b, *c, *g,*h,*d,*f,*e, *i1, *i2, *ip,MK=202682;  /*MK为级数展开时n！的n的最大值，其余为临时变量*/
    double *d1, *d2, *d3, *w,sum=0;
    time_t t_1, t_2;
    FILE  *f_out;/*f_out为输出文件*/
	//input
    printf("COS calculation to estimate the FFT benchmarks\n");
    printf("Accuracy=?\n");
    scanf("%d", &nfft);
	check1(nfft);
    printf("Float=?\n");
    scanf("%s", nn);
	check2(nn);
	for(ii=0;ii<strlen(nn);ii++)
	{
		if(nn[ii]!='0'&&nn[ii]!='.'&&nn[ii]!='-')
			break;
	}
	if(ii==strlen(nn))
	{
		printf("Don't you know cos0=1?\n");
		return 0;//cos0不需要保存，都知道是1
	}
	if(nn[0]=='-')
		check();
	memset(kk,'0',strlen(nn));
	kk[1]='.';
	kk[strlen(nn)+1]='\0';
	ii=strlen(nn);
	if(ii>22)
	{
		printf("A float cannot be so accurate!\n");
		return 0;
	}
	if(nn[1]!='.'||nn[0]!='0')
	{
		zhizhen=strstr(nn,".");
		if (zhizhen==NULL)
		{
			flag=1;
			MK=205023;
			//nfft=1000000;
			goto Init;
		}
		else
		{
			flag=2;
			jj=ii-(zhizhen-nn)-1;
		}
	}
	else
		jj=ii-2;//jj记录基底是10的负多少次方
Init:
	for(ii=1;;ii++)
	{
		sum+=log10((double)ii);
		if ((int)sum>nfft)
			break;
	}
	MK=ii;
	ii=nfft;


    printf("initializing...\n");
	//adjust nfft to 2^K
	adjust(&log2_nfft,&nfft);
	//开空间
    n = nfft + 2;
    ip = (int *) malloc((3 + (int) sqrt(0.5 * nfft)) * sizeof(int));
    w = (double *) malloc(nfft / 2 * sizeof(double));
    a = (int *) malloc((n + 2) * sizeof(int));
    b = (int *) malloc((n + 2) * sizeof(int));
    c = (int *) malloc((n + 2) * sizeof(int));
	d = (int *) malloc((n + 2) * sizeof(int));
	f = (int *) malloc((n + 2) * sizeof(int));
	g = (int *) malloc((n + 2) * sizeof(int));
	h = (int *) malloc((n + 2) * sizeof(int));
    e = (int *) malloc((n + 2) * sizeof(int));
    i1 = (int *) malloc((n + 2) * sizeof(int));
    i2 = (int *) malloc((n + 2) * sizeof(int));
    d1 = (double *) malloc((nfft + 2) * sizeof(double));
    d2 = (double *) malloc((nfft + 2) * sizeof(double));
    d3 = (double *) malloc((nfft + 2) * sizeof(double));
    if (d3 == NULL)
	{
        printf("Allocation Failure!\n");
        exit(1);
    }
    ip[0] = 0;
    /* ---- radix test ---- */
	radix_test(&radix,&err,&n,&nfft,d1,ip,w);
	if(flag!=1)
	{//调整jj=log10_radix*N使得算法可基于结点
		if(jj!=jj/log10_radix*log10_radix)
			jj=jj/log10_radix*log10_radix+log10_radix;
		if((ii/jj+1)<MK)
			MK=(ii/jj+1);
	}
    printf("MK=%d\nnfft= %d\nradix= %d\nerror_margin= %g\n",MK, nfft, radix, err);
    printf("calculating %d digits of PI...\n", log10_radix * (n - 2));
    /* ---- time check ---- */
    time(&t_1);

	memset(a, 0, (n + 2)*sizeof(int));
	memset(b, 0, (n + 2)*sizeof(int));
	BinarySplittingCOS(0,MK,radix,a,b,ip,w);
	mp_inv(n, radix, b, c, i1, i2, nfft, d1, d2, ip, w);
	mp_mulh(n, radix, c, a, c, nfft, d1, d2, ip, w);
	mp_sscanf(n, log10_radix, "1", e);
	mp_add(n, radix, e, c, c);
	//get cos
	mp_squh(n, radix, c, a, nfft, d1, ip, w);//a=c*c
	mp_sub(n, radix, e, a, b);//b=1-a
	mp_sqrt(n, radix, b, a, i1, i2, nfft, d1, d2, ip, w);//a=sqrt(b)
	//get sin
	//begin accumulating
	mp_sscanf(n, log10_radix, "1", f);
	mp_sscanf(n, log10_radix, "0", d);
	int *mm;
	result(zhizhen);
	mm=binary();
	nfft>>=1;
	n=nfft+2;
    ii=0;
	while(mm[ii]!=2)
	{
		if(mm[ii]==1)
		{
			mp_mulh(n, radix, c, f, g,nfft, d1, d2, ip, w);//g=f*c
			mp_mulh(n, radix, a, d, h,  nfft, d1, d2,  ip, w);//h=a*d
			mp_sub(n, radix, g, h, b);//b=g-h
			mp_mulh(n, radix, a, f, g,  nfft, d1, d2, ip, w);//g=f*a
			mp_mulh(n, radix, c, d, h, nfft, d1, d2,  ip, w);//h=c*d
			mp_add(n, radix, g, h, d);//d=g+h
			mp_copy(n,radix, b, f);//f=b
		}
		mp_mulh(n, radix, c, a, b, nfft, d1, d2, ip, w);//b=a*c
		mp_add(n, radix, b, b, a);//a=b*2
		mp_squh(n, radix, c, b, nfft, d1, ip, w);//b=c*c
		mp_add(n, radix, b, b, b);//b=b*2
		mp_sub(n, radix, b, e, c);//c=b-1
		ii++;

	}
    /* ---- time check ---- */
    time(&t_2);
    /* ---- output ---- */
    f_out = fopen("f:\\yzl.dat", "w");
    printf("writing cos.dat...\n");
    mp_fprintf(n, log10_radix, f, f_out);
    fclose(f_out);
    free(d3);
    free(d2);
    free(d1);
    free(i2);
    free(i1);
    free(e);
    free(c);
    free(b);
	free(d);
	free(f);
	free(g);
	free(h);
    free(a);
    free(w);
    free(ip);
    /* ---- benchmark ---- */
    n_op = 50.0 * nfft * log2_nfft * log2_nfft;
    printf("floating point operation: %g op.\n", n_op);
    d_time = difftime(t_2, t_1);
    printf("execution time: %g sec. (real time)\n", d_time);
    return 0;
}


/* -------- multiple precision routines -------- */


#include <math.h>
#include <float.h>
#include <stdio.h>

/* ---- floating point format ----
    data := data[0] * pow(radix, data[1]) *
            (data[2] + data[3]/radix + data[4]/radix/radix + ...),
    data[0]       : sign (1;data>0, -1;data<0, 0;data==0)
    data[1]       : exponent (0;data==0)
    data[2...n+1] : digits
   ---- function prototypes ----
    void mp_load_0(int n, int radix, int out[]);
    void mp_load_1(int n, int radix, int out[]);
    void mp_copy(int n, int radix, int in[], int out[]);
    void mp_round(int n, int radix, int m, int inout[]);
    int mp_cmp(int n, int radix, int in1[], int in2[]);
    void mp_add(int n, int radix, int in1[], int in2[], int out[]);
    void mp_sub(int n, int radix, int in1[], int in2[], int out[]);
    void mp_imul(int n, int radix, int in1[], int in2, int out[]);
    int mp_idiv(int n, int radix, int in1[], int in2, int out[]);
    void mp_idiv_2(int n, int radix, int in[], int out[]);
    double mp_mul_radix_test(int n, int radix, int nfft,
            double tmpfft[], int ip[], double w[]);
    void mp_mul(int n, int radix, int in1[], int in2[], int out[],
            int tmp[], int nfft, double tmp1fft[], double tmp2fft[],
            double tmp3fft[], int ip[], double w[]);
    void mp_squ(int n, int radix, int in[], int out[], int tmp[],
            int nfft, double tmp1fft[], double tmp2fft[],
            int ip[], double w[]);
    void mp_mulh(int n, int radix, int in1[], int in2[], int out[],
            int nfft, double in1fft[], double outfft[],
            int ip[], double w[]);
    void mp_squh(int n, int radix, int in[], int out[],
            int nfft, double inoutfft[], int ip[], double w[]);
    int mp_inv(int n, int radix, int in[], int out[],
            int tmp1[], int tmp2[], int nfft,
            double tmp1fft[], double tmp2fft[], int ip[], double w[]);
    int mp_sqrt(int n, int radix, int in[], int out[],
            int tmp1[], int tmp2[], int nfft,
            double tmp1fft[], double tmp2fft[], int ip[], double w[]);
    void mp_sprintf(int n, int log10_radix, int in[], char out[]);
    void mp_sscanf(int n, int log10_radix, char in[], int out[]);
    void mp_fprintf(int n, int log10_radix, int in[], FILE *fout);
   ----
*/


/* -------- mp_load routines -------- */

/**/
void mp_load_0(int n, int radix, int out[])
{
    int j;

    for (j = 0; j <= n + 1; j++) {
        out[j] = 0;
    }
}


void mp_load_1(int n, int radix, int out[])
{
    int j;

    out[0] = 1;
    out[1] = 0;
    out[2] = 1;
    for (j = 3; j <= n + 1; j++) {
        out[j] = 0;
    }
}


void mp_copy(int n, int radix, int in[], int out[])
{
    int j;

    for (j = 0; j <= n + 1; j++) {
        out[j] = in[j];
    }
}


void mp_round(int n, int radix, int m, int inout[])
{
    int j, x;

    if (m < n) {
        for (j = n + 1; j > m + 2; j--) {
            inout[j] = 0;
        }
        x = 2 * inout[m + 2];
        inout[m + 2] = 0;
        if (x >= radix) {
            for (j = m + 1; j >= 2; j--) {
                x = inout[j] + 1;
                if (x < radix) {
                    inout[j] = x;
                    break;
                }
                inout[j] = 0;
            }
            if (x >= radix) {
                inout[2] = 1;
                inout[1]++;
            }
        }
    }
}


/* -------- mp_add routines -------- */


int mp_cmp(int n, int radix, int in1[], int in2[])
{
    int mp_unsgn_cmp(int n, int in1[], int in2[]);

    if (in1[0] > in2[0]) {
        return 1;
    } else if (in1[0] < in2[0]) {
        return -1;
    }
    return in1[0] * mp_unsgn_cmp(n, &in1[1], &in2[1]);
}


void mp_add(int n, int radix, int in1[], int in2[], int out[])
{
    int mp_unsgn_cmp(int n, int in1[], int in2[]);
    int mp_unexp_add(int n, int radix, int expdif,
            int in1[], int in2[], int out[]);
    int mp_unexp_sub(int n, int radix, int expdif,
            int in1[], int in2[], int out[]);
    int outsgn, outexp, expdif;

    expdif = in1[1] - in2[1];
    outexp = in1[1];
    if (expdif < 0) {
        outexp = in2[1];
    }
    outsgn = in1[0] * in2[0];
    if (outsgn >= 0) {
        if (outsgn > 0) {
            outsgn = in1[0];
        } else {
            outsgn = in1[0] + in2[0];
            outexp = in1[1] + in2[1];
            expdif = 0;
        }
        if (expdif >= 0) {
            outexp += mp_unexp_add(n, radix, expdif,
                    &in1[2], &in2[2], &out[2]);
        } else {
            outexp += mp_unexp_add(n, radix, -expdif,
                    &in2[2], &in1[2], &out[2]);
        }
    } else {
        outsgn = mp_unsgn_cmp(n, &in1[1], &in2[1]);
        if (outsgn >= 0) {
            expdif = mp_unexp_sub(n, radix, expdif,
                    &in1[2], &in2[2], &out[2]);
        } else {
            expdif = mp_unexp_sub(n, radix, -expdif,
                    &in2[2], &in1[2], &out[2]);
        }
        outexp -= expdif;
        outsgn *= in1[0];
        if (expdif == n) {
            outsgn = 0;
        }
    }
    if (outsgn == 0) {
        outexp = 0;
    }
    out[0] = outsgn;
    out[1] = outexp;
}


void mp_sub(int n, int radix, int in1[], int in2[], int out[])
{
    int mp_unsgn_cmp(int n, int in1[], int in2[]);
    int mp_unexp_add(int n, int radix, int expdif,
            int in1[], int in2[], int out[]);
    int mp_unexp_sub(int n, int radix, int expdif,
            int in1[], int in2[], int out[]);
    int outsgn, outexp, expdif;

    expdif = in1[1] - in2[1];
    outexp = in1[1];
    if (expdif < 0) {
        outexp = in2[1];
    }
    outsgn = in1[0] * in2[0];
    if (outsgn <= 0) {
        if (outsgn < 0) {
            outsgn = in1[0];
        } else {
            outsgn = in1[0] - in2[0];
            outexp = in1[1] + in2[1];
            expdif = 0;
        }
        if (expdif >= 0) {
            outexp += mp_unexp_add(n, radix, expdif,
                    &in1[2], &in2[2], &out[2]);
        } else {
            outexp += mp_unexp_add(n, radix, -expdif,
                    &in2[2], &in1[2], &out[2]);
        }
    } else {
        outsgn = mp_unsgn_cmp(n, &in1[1], &in2[1]);
        if (outsgn >= 0) {
            expdif = mp_unexp_sub(n, radix, expdif,
                    &in1[2], &in2[2], &out[2]);
        } else {
            expdif = mp_unexp_sub(n, radix, -expdif,
                    &in2[2], &in1[2], &out[2]);
        }
        outexp -= expdif;
        outsgn *= in1[0];
        if (expdif == n) {
            outsgn = 0;
        }
    }
    if (outsgn == 0) {
        outexp = 0;
    }
    out[0] = outsgn;
    out[1] = outexp;
}


/* -------- mp_add child routines -------- */


int mp_unsgn_cmp(int n, int in1[], int in2[])
{
    int j, cmp;

    cmp = 0;
    for (j = 0; j <= n && cmp == 0; j++) {
        cmp = in1[j] - in2[j];
    }
    if (cmp > 0) {
        cmp = 1;
    } else if (cmp < 0) {
        cmp = -1;
    }
    return cmp;
}


int mp_unexp_add(int n, int radix, int expdif,
        int in1[], int in2[], int out[])
{
    int j, x, carry;

    carry = 0;
    if (expdif == 0 && in1[0] + in2[0] >= radix) {
        x = in1[n - 1] + in2[n - 1];
        carry = x >= radix ? -1 : 0;
        for (j = n - 1; j > 0; j--) {
            x = in1[j - 1] + in2[j - 1] - carry;
            carry = x >= radix ? -1 : 0;
            out[j] = x - (radix & carry);
        }
        out[0] = -carry;
    } else {
        if (expdif > n) {
            expdif = n;
        }
        for (j = n - 1; j >= expdif; j--) {
            x = in1[j] + in2[j - expdif] - carry;
            carry = x >= radix ? -1 : 0;
            out[j] = x - (radix & carry);
        }
        for (j = expdif - 1; j >= 0; j--) {
            x = in1[j] - carry;
            carry = x >= radix ? -1 : 0;
            out[j] = x - (radix & carry);
        }
        if (carry != 0) {
            for (j = n - 1; j > 0; j--) {
                out[j] = out[j - 1];
            }
            out[0] = -carry;
        }
    }
    return -carry;
}


int mp_unexp_sub(int n, int radix, int expdif,
        int in1[], int in2[], int out[])
{
    int j, x, borrow, ncancel;

    if (expdif > n) {
        expdif = n;
    }
    borrow = 0;
    for (j = n - 1; j >= expdif; j--) {
        x = in1[j] - in2[j - expdif] + borrow;
        borrow = x < 0 ? -1 : 0;
        out[j] = x + (radix & borrow);
    }
    for (j = expdif - 1; j >= 0; j--) {
        x = in1[j] + borrow;
        borrow = x < 0 ? -1 : 0;
        out[j] = x + (radix & borrow);
    }
    ncancel = 0;
    for (j = 0; j < n && out[j] == 0; j++) {
        ncancel = j + 1;
    }
    if (ncancel > 0 && ncancel < n) {
        for (j = 0; j < n - ncancel; j++) {
            out[j] = out[j + ncancel];
        }
        for (j = n - ncancel; j < n; j++) {
            out[j] = 0;
        }
    }
    return ncancel;
}


/* -------- mp_imul routines -------- */

/*数乘：将数据in1[]*一个int型整数in2得到out[]*/
void mp_imul(int n, int radix, int in1[], int in2, int out[])
{
    void mp_unsgn_imul(int n, double dradix, int in1[], double din2,
            int out[]);

    if (in2 > 0) {
        out[0] = in1[0];
    } else if (in2 < 0) {
        out[0] = -in1[0];
        in2 = -in2;
    } else {
        out[0] = 0;
    }
    mp_unsgn_imul(n, radix, &in1[1], in2, &out[1]);
    if (out[0] == 0) {
        out[1] = 0;
    }
}

/*数除：将数据in1[]除以一个int型整数in2，得到数据out[]*/
int mp_idiv(int n, int radix, int in1[], int in2, int out[])
{
    void mp_load_0(int n, int radix, int out[]);
    void mp_unsgn_idiv(int n, double dradix, int in1[], double din2,
            int out[]);

    if (in2 == 0) {
        return -1;
    }
    if (in2 > 0) {
        out[0] = in1[0];
    } else {
        out[0] = -in1[0];
        in2 = -in2;
    }
    if (in1[0] == 0) {
        mp_load_0(n, radix, out);
        return 0;
    }
    mp_unsgn_idiv(n, radix, &in1[1], in2, &out[1]);
    return 0;
}

/*为了提高程序效率节省时间，数除时除数是2时而编的一个函数，即数据in[]/2得到out[]*/
void mp_idiv_2(int n, int radix, int in[], int out[])
{
    int j, ix, carry, shift;

    out[0] = in[0];
    shift = 0;
    if (in[2] == 1) {
        shift = 1;
    }
    out[1] = in[1] - shift;
    carry = -shift;
    for (j = 2; j <= n + 1 - shift; j++) {
        ix = in[j + shift] + (radix & carry);
        carry = -(ix & 1);
        out[j] = ix >> 1;
    }
    if (shift > 0) {
        out[n + 1] = (radix & carry) >> 1;
    }
}


/* -------- mp_imul child routines -------- */


void mp_unsgn_imul(int n, double dradix, int in1[], double din2,
        int out[])
{
    int j, carry, shift;
    double x, d1_radix;

    d1_radix = 1.0 / dradix;
    carry = 0;
    for (j = n; j >= 1; j--) {
        x = din2 * in1[j] + carry + 0.5;
        carry = (int) (d1_radix * x);
        out[j] = (int) (x - dradix * carry);
    }
    shift = 0;
    x = carry + 0.5;
    while (x > 1) {
        x *= d1_radix;
        shift++;
    }
    out[0] = in1[0] + shift;
    if (shift > 0) {
        while (shift > n) {
            carry = (int) (d1_radix * carry + 0.5);
            shift--;
        }
        for (j = n; j >= shift + 1; j--) {
            out[j] = out[j - shift];
        }
        for (j = shift; j >= 1; j--) {
            x = carry + 0.5;
            carry = (int) (d1_radix * x);
            out[j] = (int) (x - dradix * carry);
        }
    }
}


void mp_unsgn_idiv(int n, double dradix, int in1[], double din2,
        int out[])
{
    int j, ix, carry, shift;
    double x, d1_in2;

    d1_in2 = 1.0 / din2;
    shift = 0;
    x = 0;
    do {
        shift++;
        x *= dradix;
        if (shift <= n) {
            x += in1[shift];
        }
    } while (x < din2 - 0.5);
    x += 0.5;
    ix = (int) (d1_in2 * x);
    carry = (int) (x - din2 * ix);
    out[1] = ix;
    shift--;
    out[0] = in1[0] - shift;
    if (shift >= n) {
        shift = n - 1;
    }
    for (j = 2; j <= n - shift; j++) {
        x = in1[j + shift] + dradix * carry + 0.5;
        ix = (int) (d1_in2 * x);
        carry = (int) (x - din2 * ix);
        out[j] = ix;
    }
    for (j = n - shift + 1; j <= n; j++) {
        x = dradix * carry + 0.5;
        ix = (int) (d1_in2 * x);
        carry = (int) (x - din2 * ix);
        out[j] = ix;
    }
}


/* -------- mp_mul routines -------- */

/*计算*/
double mp_mul_radix_test(int n, int radix, int nfft,
        double tmpfft[], int ip[], double w[])
{
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_csqu(int nfft, double dinout[]);
    double mp_mul_d2i_test(int radix, int nfft, double din[]);
    int j, ndata, radix_2;

    ndata = (nfft >> 1) + 1;
    if (ndata > n) {
        ndata = n;
    }
    tmpfft[nfft + 1] = radix - 1;
    for (j = nfft; j > ndata; j--) {
        tmpfft[j] = 0;
    }
    radix_2 = (radix + 1) / 2;
    for (j = ndata; j > 2; j--) {
        tmpfft[j] = radix_2;
    }
    tmpfft[2] = radix;
    tmpfft[1] = radix - 1;
    tmpfft[0] = 0;
    rdft(nfft, 1, &tmpfft[1], ip, w);
    mp_mul_csqu(nfft, tmpfft);
    rdft(nfft, -1, &tmpfft[1], ip, w);
    return 2 * mp_mul_d2i_test(radix, nfft, tmpfft);
}


void mp_mul(int n, int radix, int in1[], int in2[], int out[],
        int tmp[], int nfft, double tmp1fft[], double tmp2fft[],
        double tmp3fft[], int ip[], double w[])
{
    void mp_copy(int n, int radix, int in[], int out[]);
    void mp_add(int n, int radix, int in1[], int in2[], int out[]);
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_i2d(int n, int radix, int nfft, int shift,
            int in[], double dout[]);
    void mp_mul_cmul(int nfft, double din[], double dinout[]);
    void mp_mul_cmuladd(int nfft, double din1[], double din2[],
            double dinout[]);
    void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[]);
    int n_h, shift;

    shift = (nfft >> 1) + 1;
    while (n > shift) {
        if (in1[shift + 2] + in2[shift + 2] != 0) {
            break;
        }
        shift++;
    }
    n_h = n / 2 + 1;
    if (n_h < n - shift) {
        n_h = n - shift;
    }
    /* ---- tmp3fft = (upper) in1 * (lower) in2 ---- */
    mp_mul_i2d(n, radix, nfft, 0, in1, tmp1fft);
    rdft(nfft, 1, &tmp1fft[1], ip, w);
    mp_mul_i2d(n, radix, nfft, shift, in2, tmp3fft);
    rdft(nfft, 1, &tmp3fft[1], ip, w);
    mp_mul_cmul(nfft, tmp1fft, tmp3fft);
    /* ---- tmp = (upper) in1 * (upper) in2 ---- */
    mp_mul_i2d(n, radix, nfft, 0, in2, tmp2fft);
    rdft(nfft, 1, &tmp2fft[1], ip, w);
    mp_mul_cmul(nfft, tmp2fft, tmp1fft);
    rdft(nfft, -1, &tmp1fft[1], ip, w);
    mp_mul_d2i(n, radix, nfft, tmp1fft, tmp);
    /* ---- tmp3fft += (upper) in2 * (lower) in1 ---- */
    mp_mul_i2d(n, radix, nfft, shift, in1, tmp1fft);
    rdft(nfft, 1, &tmp1fft[1], ip, w);
    mp_mul_cmuladd(nfft, tmp1fft, tmp2fft, tmp3fft);
    /* ---- out = tmp + tmp3fft ---- */
    rdft(nfft, -1, &tmp3fft[1], ip, w);
    mp_mul_d2i(n_h, radix, nfft, tmp3fft, out);
    if (out[0] != 0) {
        mp_add(n, radix, out, tmp, out);
    } else {
        mp_copy(n, radix, tmp, out);
    }
}


void mp_squ(int n, int radix, int in[], int out[], int tmp[],
        int nfft, double tmp1fft[], double tmp2fft[],
        int ip[], double w[])
{
    void mp_add(int n, int radix, int in1[], int in2[], int out[]);
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_i2d(int n, int radix, int nfft, int shift,
            int in[], double dout[]);
    void mp_mul_cmul(int nfft, double din[], double dinout[]);
    void mp_mul_csqu(int nfft, double dinout[]);
    void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[]);
    int n_h, shift;

    shift = (nfft >> 1) + 1;
    while (n > shift) {
        if (in[shift + 2] != 0) {
            break;
        }
        shift++;
    }
    n_h = n / 2 + 1;
    if (n_h < n - shift) {
        n_h = n - shift;
    }
    /* ---- tmp = (upper) in * (lower) in ---- */
    mp_mul_i2d(n, radix, nfft, 0, in, tmp1fft);
    rdft(nfft, 1, &tmp1fft[1], ip, w);
    mp_mul_i2d(n, radix, nfft, shift, in, tmp2fft);
    rdft(nfft, 1, &tmp2fft[1], ip, w);
    mp_mul_cmul(nfft, tmp1fft, tmp2fft);
    rdft(nfft, -1, &tmp2fft[1], ip, w);
    mp_mul_d2i(n_h, radix, nfft, tmp2fft, tmp);
    /* ---- out = 2 * tmp + ((upper) in)^2 ---- */
    mp_mul_csqu(nfft, tmp1fft);
    rdft(nfft, -1, &tmp1fft[1], ip, w);
    mp_mul_d2i(n, radix, nfft, tmp1fft, out);
    if (tmp[0] != 0) {
        mp_add(n_h, radix, tmp, tmp, tmp);
        mp_add(n, radix, out, tmp, out);
    }
}


void mp_mulh(int n, int radix, int in1[], int in2[], int out[],
        int nfft, double in1fft[], double outfft[], int ip[], double w[])
{
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_i2d(int n, int radix, int nfft, int shift,
            int in[], double dout[]);
    void mp_mul_cmul(int nfft, double din[], double dinout[]);
    void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[]);

    mp_mul_i2d(n, radix, nfft, 0, in1, in1fft);
    rdft(nfft, 1, &in1fft[1], ip, w);
    mp_mul_i2d(n, radix, nfft, 0, in2, outfft);
    rdft(nfft, 1, &outfft[1], ip, w);
    mp_mul_cmul(nfft, in1fft, outfft);
    rdft(nfft, -1, &outfft[1], ip, w);
    mp_mul_d2i(n, radix, nfft, outfft, out);
}


void mp_mulh_use_in1fft(int n, int radix, double in1fft[],
        int shift, int in2[], int out[], int nfft, double outfft[],
        int ip[], double w[])
{
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_i2d(int n, int radix, int nfft, int shift,
            int in[], double dout[]);
    void mp_mul_cmul(int nfft, double din[], double dinout[]);
    void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[]);
    int n_h;

    while (n > shift) {
        if (in2[shift + 2] != 0) {
            break;
        }
        shift++;
    }
    n_h = n / 2 + 1;
    if (n_h < n - shift) {
        n_h = n - shift;
    }
    mp_mul_i2d(n, radix, nfft, shift, in2, outfft);
    rdft(nfft, 1, &outfft[1], ip, w);
    mp_mul_cmul(nfft, in1fft, outfft);
    rdft(nfft, -1, &outfft[1], ip, w);
    mp_mul_d2i(n_h, radix, nfft, outfft, out);
}


void mp_squh(int n, int radix, int in[], int out[],
        int nfft, double inoutfft[], int ip[], double w[])
{
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_i2d(int n, int radix, int nfft, int shift,
            int in[], double dout[]);
    void mp_mul_csqu(int nfft, double dinout[]);
    void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[]);

    mp_mul_i2d(n, radix, nfft, 0, in, inoutfft);
    rdft(nfft, 1, &inoutfft[1], ip, w);
    mp_mul_csqu(nfft, inoutfft);
    rdft(nfft, -1, &inoutfft[1], ip, w);
    mp_mul_d2i(n, radix, nfft, inoutfft, out);
}


void mp_squh_use_in1fft(int n, int radix, double inoutfft[], int out[],
        int nfft, int ip[], double w[])
{
    void rdft(int n, int isgn, double *a, int *ip, double *w);
    void mp_mul_csqu(int nfft, double dinout[]);
    void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[]);

    mp_mul_csqu(nfft, inoutfft);
    rdft(nfft, -1, &inoutfft[1], ip, w);
    mp_mul_d2i(n, radix, nfft, inoutfft, out);
}


/* -------- mp_mul child routines -------- */


void mp_mul_i2d(int n, int radix, int nfft, int shift,
        int in[], double dout[])
{
    int j, x, carry, ndata, radix_2, topdgt;

    ndata = 0;
    topdgt = 0;
    if (n > shift) {
        topdgt = in[shift + 2];
        ndata = (nfft >> 1) + 1;
        if (ndata > n - shift) {
            ndata = n - shift;
        }
    }
    dout[nfft + 1] = in[0] * topdgt;
    for (j = nfft; j > ndata; j--) {
        dout[j] = 0;
    }
    /* ---- abs(dout[j]) <= radix/2 (to keep FFT precision) ---- */
    if (ndata > 1) {
        radix_2 = radix / 2;
        carry = 0;
        for (j = ndata + 1; j > 3; j--) {
            x = in[j + shift] - carry;
            carry = x >= radix_2 ? -1 : 0;
            dout[j - 1] = x - (radix & carry);
        }
        dout[2] = in[shift + 3] - carry;
    }
    dout[1] = topdgt;
    dout[0] = in[1] - shift;
}


void mp_mul_cmul(int nfft, double din[], double dinout[])
{
    int j;
    double xr, xi, yr, yi;

    dinout[0] += din[0];
    dinout[1] *= din[1];
    dinout[2] *= din[2];
    for (j = 3; j < nfft; j += 2) {
        xr = din[j];
        xi = din[j + 1];
        yr = dinout[j];
        yi = dinout[j + 1];
        dinout[j] = xr * yr - xi * yi;
        dinout[j + 1] = xr * yi + xi * yr;
    }
    dinout[nfft + 1] *= din[nfft + 1];
}


void mp_mul_cmuladd(int nfft, double din1[], double din2[],
        double dinout[])
{
    int j;
    double xr, xi, yr, yi;

    dinout[1] += din1[1] * din2[1];
    dinout[2] += din1[2] * din2[2];
    for (j = 3; j < nfft; j += 2) {
        xr = din1[j];
        xi = din1[j + 1];
        yr = din2[j];
        yi = din2[j + 1];
        dinout[j] += xr * yr - xi * yi;
        dinout[j + 1] += xr * yi + xi * yr;
    }
    dinout[nfft + 1] += din1[nfft + 1] * din2[nfft + 1];
}


void mp_mul_csqu(int nfft, double dinout[])
{
    int j;
    double xr, xi;

    dinout[0] *= 2;
    dinout[1] *= dinout[1];
    dinout[2] *= dinout[2];
    for (j = 3; j < nfft; j += 2) {
        xr = dinout[j];
        xi = dinout[j + 1];
        dinout[j] = xr * xr - xi * xi;
        dinout[j + 1] = 2 * xr * xi;
    }
    dinout[nfft + 1] *= dinout[nfft + 1];
}


void mp_mul_d2i(int n, int radix, int nfft, double din[], int out[])
{
    int j, carry, carry1, carry2, shift, ndata;
    double x, scale, d1_radix, d1_radix2, pow_radix, topdgt;

    scale = 2.0 / nfft;
    d1_radix = 1.0 / radix;
    d1_radix2 = d1_radix * d1_radix;
    topdgt = din[nfft + 1];
    x = topdgt < 0 ? -topdgt : topdgt;
    shift = x + 0.5 >= radix ? 1 : 0;
    /* ---- correction of cyclic convolution of din[1] ---- */
    x *= nfft * 0.5;
    din[nfft + 1] = din[1] - x;
    din[1] = x;
    /* ---- output of digits ---- */
    ndata = n;
    if (n > nfft + 1 + shift) {
        ndata = nfft + 1 + shift;
        for (j = n + 1; j > ndata + 1; j--) {
            out[j] = 0;
        }
    }
    x = 0;
    pow_radix = 1;
    for (j = ndata + 1 - shift; j <= nfft + 1; j++) {
        x += pow_radix * din[j];
        pow_radix *= d1_radix;
        if (pow_radix < DBL_EPSILON) {
            break;
        }
    }
    x = d1_radix2 * (scale * x + 0.5);
    carry2 = ((int) x) - 1;
    carry = (int) (radix * (x - carry2) + 0.5);
    for (j = ndata; j > 1; j--) {
        x = d1_radix2 * (scale * din[j - shift] + carry + 0.5);
        carry = carry2;
        carry2 = ((int) x) - 1;
        x = radix * (x - carry2);
        carry1 = (int) x;
        out[j + 1] = (int) (radix * (x - carry1));
        carry += carry1;
    }
    x = carry + ((double) radix) * carry2 + 0.5;
    if (shift == 0) {
        x += scale * din[1];
    }
    carry = (int) (d1_radix * x);
    out[2] = (int) (x - ((double) radix) * carry);
    if (carry > 0) {
        for (j = n + 1; j > 2; j--) {
            out[j] = out[j - 1];
        }
        out[2] = carry;
        shift++;
    }
    /* ---- output of exp, sgn ---- */
    x = din[0] + shift + 0.5;
    shift = ((int) x) - 1;
    out[1] = shift + ((int) (x - shift));
    out[0] = topdgt > 0.5 ? 1 : -1;
    if (out[2] == 0) {
        out[0] = 0;
        out[1] = 0;
    }
}


double mp_mul_d2i_test(int radix, int nfft, double din[])
{
    int j, carry, carry1, carry2;
    double x, scale, d1_radix, d1_radix2, err;

    scale = 2.0 / nfft;
    d1_radix = 1.0 / radix;
    d1_radix2 = d1_radix * d1_radix;
    /* ---- correction of cyclic convolution of din[1] ---- */
    x = din[nfft + 1] * nfft * 0.5;
    if (x < 0) {
        x = -x;
    }
    din[nfft + 1] = din[1] - x;
    /* ---- check of digits ---- */
    err = 0;
    carry = 0;
    carry2 = 0;
    for (j = nfft + 1; j > 1; j--) {
        x = d1_radix2 * (scale * din[j] + carry + 0.5);
        carry = carry2;
        carry2 = ((int) x) - 1;
        x = radix * (x - carry2);
        carry1 = (int) x;
        x = radix * (x - carry1);
        carry += carry1;
        x = x - 0.5 - ((int) x);
        if (x > err) {
            err = x;
        } else if (-x > err) {
            err = -x;
        }
    }
    return err;
}


/* -------- mp_inv routines -------- */


int mp_inv(int n, int radix, int in[], int out[],
        int tmp1[], int tmp2[], int nfft,
        double tmp1fft[], double tmp2fft[], int ip[], double w[])
{
    int mp_get_nfft_init(int radix, int nfft_max);
    void mp_inv_init(int n, int radix, int in[], int out[]);
    int mp_inv_newton(int n, int radix, int in[], int inout[],
            int tmp1[], int tmp2[], int nfft, double tmp1fft[],
            double tmp2fft[], int ip[], double w[]);
    int n_nwt, nfft_nwt, thr, prc;

    if (in[0] == 0) {
        return -1;
    }
    nfft_nwt = mp_get_nfft_init(radix, nfft);
    n_nwt = nfft_nwt + 2;
    if (n_nwt > n) {
        n_nwt = n;
    }
    mp_inv_init(n_nwt, radix, in, out);
    thr = 8;
    do {
        n_nwt = nfft_nwt + 2;
        if (n_nwt > n) {
            n_nwt = n;
        }
        prc = mp_inv_newton(n_nwt, radix, in, out,
                tmp1, tmp2, nfft_nwt, tmp1fft, tmp2fft, ip, w);
        if (thr * nfft_nwt >= nfft) {
            thr = 0;
            if (2 * prc <= n_nwt - 2) {
                nfft_nwt >>= 1;
            }
        } else {
            if (3 * prc < n_nwt - 2) {
                nfft_nwt >>= 1;
            }
        }
        nfft_nwt <<= 1;
    } while (nfft_nwt <= nfft);
    return 0;
}


int mp_sqrt(int n, int radix, int in[], int out[],
        int tmp1[], int tmp2[], int nfft,
        double tmp1fft[], double tmp2fft[], int ip[], double w[])
{
    void mp_load_0(int n, int radix, int out[]);
    int mp_get_nfft_init(int radix, int nfft_max);
    void mp_sqrt_init(int n, int radix, int in[], int out[], int out_rev[]);
    int mp_sqrt_newton(int n, int radix, int in[], int inout[],
            int inout_rev[], int tmp[], int nfft, double tmp1fft[],
            double tmp2fft[], int ip[], double w[], int *n_tmp1fft);
    int n_nwt, nfft_nwt, thr, prc, n_tmp1fft;

    if (in[0] < 0) {
        return -1;
    } else if (in[0] == 0) {
        mp_load_0(n, radix, out);
        return 0;
    }
    nfft_nwt = mp_get_nfft_init(radix, nfft);
    n_nwt = nfft_nwt + 2;
    if (n_nwt > n) {
        n_nwt = n;
    }
    mp_sqrt_init(n_nwt, radix, in, out, tmp1);
    n_tmp1fft = 0;
    thr = 8;
    do {
        n_nwt = nfft_nwt + 2;
        if (n_nwt > n) {
            n_nwt = n;
        }
        prc = mp_sqrt_newton(n_nwt, radix, in, out,
                tmp1, tmp2, nfft_nwt, tmp1fft, tmp2fft,
                ip, w, &n_tmp1fft);
        if (thr * nfft_nwt >= nfft) {
            thr = 0;
            if (2 * prc <= n_nwt - 2) {
                nfft_nwt >>= 1;
            }
        } else {
            if (3 * prc < n_nwt - 2) {
                nfft_nwt >>= 1;
            }
        }
        nfft_nwt <<= 1;
    } while (nfft_nwt <= nfft);
    return 0;
}


/* -------- mp_inv child routines -------- */


int mp_get_nfft_init(int radix, int nfft_max)
{
    int nfft_init;
    double r;

    r = radix;
    nfft_init = 1;
    do {
        r *= r;
        nfft_init <<= 1;
    } while (DBL_EPSILON * r < 1 && nfft_init < nfft_max);
    return nfft_init;
}


void mp_inv_init(int n, int radix, int in[], int out[])
{
    void mp_unexp_d2mp(int n, int radix, double din, int out[]);
    double mp_unexp_mp2d(int n, int radix, int in[]);
    int outexp;
    double din;

    out[0] = in[0];
    outexp = -in[1];
    din = 1.0 / mp_unexp_mp2d(n, radix, &in[2]);
    while (din < 1) {
        din *= radix;
        outexp--;
    }
    out[1] = outexp;
    mp_unexp_d2mp(n, radix, din, &out[2]);
}


void mp_sqrt_init(int n, int radix, int in[], int out[], int out_rev[])
{
    void mp_unexp_d2mp(int n, int radix, double din, int out[]);
    double mp_unexp_mp2d(int n, int radix, int in[]);
    int outexp;
    double din;

    out[0] = 1;
    out_rev[0] = 1;
    outexp = in[1];
    din = mp_unexp_mp2d(n, radix, &in[2]);
    if (outexp % 2 != 0) {
        din *= radix;
        outexp--;
    }
    outexp /= 2;
    din = sqrt(din);
    if (din < 1) {
        din *= radix;
        outexp--;
    }
    out[1] = outexp;
    mp_unexp_d2mp(n, radix, din, &out[2]);
    outexp = -outexp;
    din = 1.0 / din;
    while (din < 1) {
        din *= radix;
        outexp--;
    }
    out_rev[1] = outexp;
    mp_unexp_d2mp(n, radix, din, &out_rev[2]);
}


void mp_unexp_d2mp(int n, int radix, double din, int out[])
{
    int j, x;

    for (j = 0; j < n; j++) {
        x = (int) din;
        if (x >= radix) {
            x = radix - 1;
            din = radix;
        }
        din = radix * (din - x);
        out[j] = x;
    }
}


double mp_unexp_mp2d(int n, int radix, int in[])
{
    int j;
    double d1_radix, dout;

    d1_radix = 1.0 / radix;
    dout = 0;
    for (j = n - 1; j >= 0; j--) {
        dout = d1_radix * dout + in[j];
    }
    return dout;
}


int mp_inv_newton(int n, int radix, int in[], int inout[],
        int tmp1[], int tmp2[], int nfft, double tmp1fft[],
        double tmp2fft[], int ip[], double w[])
{
    void mp_load_1(int n, int radix, int out[]);
    void mp_round(int n, int radix, int m, int inout[]);
    void mp_add(int n, int radix, int in1[], int in2[], int out[]);
    void mp_sub(int n, int radix, int in1[], int in2[], int out[]);
    void mp_mulh(int n, int radix, int in1[], int in2[], int out[],
            int nfft, double in1fft[], double outfft[],
            int ip[], double w[]);
    void mp_mulh_use_in1fft(int n, int radix, double in1fft[],
            int shift, int in2[], int out[], int nfft, double outfft[],
            int ip[], double w[]);
    int n_h, shift, prc;

    shift = (nfft >> 1) + 1;
    n_h = n / 2 + 1;
    if (n_h < n - shift) {
        n_h = n - shift;
    }
    /* ---- tmp1 = inout * (upper) in (half to normal precision) ---- */
    mp_round(n, radix, shift, inout);
    mp_mulh(n, radix, inout, in, tmp1,
            nfft, tmp1fft, tmp2fft, ip, w);
    /* ---- tmp2 = 1 - tmp1 ---- */
    mp_load_1(n, radix, tmp2);
    mp_sub(n, radix, tmp2, tmp1, tmp2);
    /* ---- tmp2 -= inout * (lower) in (half precision) ---- */
    mp_mulh_use_in1fft(n, radix, tmp1fft, shift, in, tmp1,
            nfft, tmp2fft, ip, w);
    mp_sub(n_h, radix, tmp2, tmp1, tmp2);
    /* ---- get precision ---- */
    prc = -tmp2[1];
    if (tmp2[0] == 0) {
        prc = nfft + 1;
    }
    /* ---- tmp2 *= inout (half precision) ---- */
    mp_mulh_use_in1fft(n_h, radix, tmp1fft, 0, tmp2, tmp2,
            nfft, tmp2fft, ip, w);
    /* ---- inout += tmp2 ---- */
    if (tmp2[0] != 0) {
        mp_add(n, radix, inout, tmp2, inout);
    }
    return prc;
}


int mp_sqrt_newton(int n, int radix, int in[], int inout[],
        int inout_rev[], int tmp[], int nfft, double tmp1fft[],
        double tmp2fft[], int ip[], double w[], int *n_tmp1fft)
{
    void mp_round(int n, int radix, int m, int inout[]);
    void mp_add(int n, int radix, int in1[], int in2[], int out[]);
    void mp_sub(int n, int radix, int in1[], int in2[], int out[]);
    void mp_idiv_2(int n, int radix, int in[], int out[]);
    void mp_mulh(int n, int radix, int in1[], int in2[], int out[],
            int nfft, double in1fft[], double outfft[],
            int ip[], double w[]);
    void mp_squh(int n, int radix, int in[], int out[],
            int nfft, double inoutfft[], int ip[], double w[]);
    void mp_squh_use_in1fft(int n, int radix, double inoutfft[], int out[],
            int nfft, int ip[], double w[]);
    int n_h, nfft_h, shift, prc;

    nfft_h = nfft >> 1;
    shift = nfft_h + 1;
    if (nfft_h < 2) {
        nfft_h = 2;
    }
    n_h = n / 2 + 1;
    if (n_h < n - shift) {
        n_h = n - shift;
    }
    /* ---- tmp = inout_rev^2 (1/4 to half precision) ---- */
    mp_round(n_h, radix, (nfft_h >> 1) + 1, inout_rev);
    if (*n_tmp1fft != nfft_h) {
        mp_squh(n_h, radix, inout_rev, tmp,
                nfft_h, tmp1fft, ip, w);
    } else {
        mp_squh_use_in1fft(n_h, radix, tmp1fft, tmp,
                nfft_h, ip, w);
    }
    /* ---- tmp = inout_rev - inout * tmp (half precision) ---- */
    mp_round(n, radix, shift, inout);
    mp_mulh(n_h, radix, inout, tmp, tmp,
            nfft, tmp1fft, tmp2fft, ip, w);
    mp_sub(n_h, radix, inout_rev, tmp, tmp);
    /* ---- inout_rev += tmp ---- */
    mp_add(n_h, radix, inout_rev, tmp, inout_rev);
    /* ---- tmp = in - inout^2 (half to normal precision) ---- */
    mp_squh_use_in1fft(n, radix, tmp1fft, tmp,
            nfft, ip, w);
    mp_sub(n, radix, in, tmp, tmp);
    /* ---- get precision ---- */
    prc = in[1] - tmp[1];
    if (in[2] > tmp[2]) {
        prc++;
    }
    if (tmp[0] == 0) {
        prc = nfft + 1;
    }
    /* ---- tmp = tmp * inout_rev / 2 (half precision) ---- */
    mp_round(n_h, radix, shift, inout_rev);
    mp_mulh(n_h, radix, inout_rev, tmp, tmp,
            nfft, tmp1fft, tmp2fft, ip, w);
    *n_tmp1fft = nfft;
    mp_idiv_2(n_h, radix, tmp, tmp);
    /* ---- inout += tmp ---- */
    if (tmp[0] != 0) {
        mp_add(n, radix, inout, tmp, inout);
    }
    return prc;
}


/* -------- mp_io routines -------- */


void mp_sprintf(int n, int log10_radix, int in[], char out[])
{
    int j, k, x, y, outexp, shift;

    if (in[0] < 0) {
        *out++ = '-';
    }
    x = in[2];
    shift = log10_radix;
    for (k = log10_radix; k > 0; k--) {
        y = x % 10;
        x /= 10;
        out[k] = '0' + y;
        if (y != 0) {
            shift = k;
        }
    }
    out[0] = out[shift];
    out[1] = '.';
    for (k = 1; k <= log10_radix - shift; k++) {
        out[k + 1] = out[k + shift];
    }
    outexp = log10_radix - shift;
    out += outexp + 2;
    for (j = 3; j <= n + 1; j++) {
        x = in[j];
        for (k = log10_radix - 1; k >= 0; k--) {
            y = x % 10;
            x /= 10;
            out[k] = '0' + y;
        }
        out += log10_radix;
    }
    *out++ = 'e';
    outexp += log10_radix * in[1];
    sprintf(out, "%d", outexp);
}


void mp_sscanf(int n, int log10_radix, char in[], int out[])
{
    char *s;
    int j, x, outexp, outexp_mod;

    while (*in == ' ') {
        in++;
    }
    out[0] = 1;
    if (*in == '-') {
        out[0] = -1;
        in++;
    } else if (*in == '+') {
        in++;
    }
    while (*in == ' ' || *in == '0') {
        in++;
    }
    outexp = 0;
    for (s = in; *s != '\0'; s++) {
        if (*s == 'e' || *s == 'E' || *s == 'd' || *s == 'D') {
            if (sscanf(++s, "%d", &outexp) != 1) {
                outexp = 0;
            }
            break;
        }
    }
    if (*in == '.') {
        do {
            outexp--;
            while (*++in == ' ');
        } while (*in == '0' && *in != '\0');
    } else if (*in != '\0') {
        s = in;
        while (*++s == ' ');
        while (*s >= '0' && *s <= '9' && *s != '\0') {
            outexp++;
            while (*++s == ' ');
        }
    }
    x = outexp / log10_radix;
    outexp_mod = outexp - log10_radix * x;
    if (outexp_mod < 0) {
        x--;
        outexp_mod += log10_radix;
    }
    out[1] = x;
    x = 0;
    j = 2;
    for (s = in; *s != '\0'; s++) {
        if (*s == '.' || *s == ' ') {
            continue;
        }
        if (*s < '0' || *s > '9') {
            break;
        }
        x = 10 * x + (*s - '0');
        if (--outexp_mod < 0) {
            if (j > n + 1) {
                break;
            }
            out[j++] = x;
            x = 0;
            outexp_mod = log10_radix - 1;
        }
    }
    while (outexp_mod-- >= 0) {
        x *= 10;
    }
    while (j <= n + 1) {
        out[j++] = x;
        x = 0;
    }
    if (out[2] == 0) {
        out[0] = 0;
        out[1] = 0;
    }
}


void mp_fprintf(int n, int log10_radix, int in[], FILE *fout)
{
    int j, k, x, y, outexp, shift;
    char out[256];

    if (in[0] < 0) {
        putc('-', fout);
    }
    x = in[2];
    shift = log10_radix;
    for (k = log10_radix; k > 0; k--) {
        y = x % 10;
        x /= 10;
        out[k] = '0' + y;
        if (y != 0) {
            shift = k;
        }
    }
    putc(out[shift], fout);
    putc('.', fout);
    for (k = 1; k <= log10_radix - shift; k++) {
        putc(out[k + shift], fout);
    }
    outexp = log10_radix - shift;
    for (j = 3; j <= n + 1; j++) {
        x = in[j];
        for (k = log10_radix - 1; k >= 0; k--) {
            y = x % 10;
            x /= 10;
            out[k] = '0' + y;
        }
        for (k = 0; k < log10_radix; k++) {
            putc(out[k], fout);
        }
    }
    putc('e', fout);
    outexp += log10_radix * in[1];
    sprintf(out, "%d", outexp);
    for (k = 0; out[k] != '\0'; k++) {
        putc(out[k], fout);
    }
}
