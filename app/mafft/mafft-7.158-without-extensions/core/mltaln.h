#define USE_XCED 0

#if USE_XCED
#include "config.h"
#include "share.h"
#else
#endif

#include "mafft.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include "mtxutl.h"
#include <float.h>
#include <stdarg.h>
#ifdef enablemultithread
#include <pthread.h>
#endif

#define VERSION "7.158b"
#define SHOWVERSION reporterr( "%s (%s) Version " VERSION " alg=%c, model=%s, amax=%3.1f\n%d thread(s)\n", progName( argv[0] ), (dorp=='d')?"nuc":((nblosum==-2)?"text":"aa"), alg, modelname, specificityconsideration, nthread )

#define FFT_THRESHOLD  80
#define FFT_WINSIZE_P   20
#define FFT_WINSIZE_D   100
#define DISPSEQF  60
#define DISPSITEI    0
#define MAXITERATION 500
#define M  500000       /* njob no saidaiti */
#define N 5000000       /* nlen no saidaiti */
#define MAXSEG 100000
#define B     256
#define C     60       /*  1 gyou no mojisuu */
#define D      6
#define rnd() ( ( 1.0 / ( RAND_MAX + 1.0 ) ) * rand() )
#define MAX(X,Y)    ( ((X)>(Y))?(X):(Y) )
#define MIN(X,Y)    ( ((X)<(Y))?(X):(Y) )
#define G(X)        ( ((X)>(0))?(X):(0) )
#define BEFF   1.0   /* 0.6 ni suruto zureru */
#define WIN   3
#define SGAPP -1000
#define GETA2 0.001
#define GETA3 0.001
#define NOTSPECIFIED  100009
#define SUEFF   0.1  /* upg/(spg+upg)  -> sueff.sed */ 
#define DIVLOCAL 0
#define INTMTXSCALE 1000000.0
#define JTT 201
#define TM  202

extern char modelname[500];
extern int njob, nlenmax;
extern int amino_n[0x80];
extern char amino_grp[0x80];
extern int amino_dis[0x80][0x80];
extern double **n_disLN;
extern double amino_dis_consweight_multi[0x80][0x80];
extern int **n_dis;
extern int **n_disFFT;
extern double **n_dis_consweight_multi;
extern char amino[0x80];
extern double polarity[0x80];
extern double volume[0x80];
extern int ribosumdis[37][37];

extern int ppid;
extern double thrinter;
extern double fastathreshold;
extern int pslocal, ppslocal;
extern int constraint;
extern int divpairscore;
extern int fmodel; // 1-> fmodel 0->default -1->raw
extern int nblosum; // 45, 50, 62, 80
extern int kobetsubunkatsu;
extern int bunkatsu;
extern int dorp;
extern int niter;
extern int contin;
extern int calledByXced;
extern int devide;
extern int scmtd;
extern int weight;
extern int utree;
extern int tbutree;
extern int refine;
extern int check;
extern double cut;
extern int cooling;
extern int trywarp;
extern int penalty, ppenalty, penaltyLN;
extern int penalty_dist, ppenalty_dist;
extern int RNApenalty, RNAppenalty;
extern int RNApenalty_ex, RNAppenalty_ex;
extern int penalty_ex, ppenalty_ex, penalty_exLN;
extern int penalty_EX, ppenalty_EX;
extern int penalty_OP, ppenalty_OP;
extern int penalty_shift;
extern double penalty_shift_factor;
extern int offset, poffset, offsetLN, offsetFFT;
extern int RNAthr, RNApthr;
extern int scoremtx;
extern int TMorJTT;
extern char use_fft;
extern char force_fft;
extern int nevermemsave;
extern int fftscore;
extern int fftWinSize;
extern int fftThreshold;
extern int fftRepeatStop;
extern int fftNoAnchStop;
extern int divWinSize;
extern int divThreshold;
extern int disp;
extern int outgap;
extern char alg;
extern int cnst;
extern int mix;
extern int tbitr;
extern int tbweight;
extern int tbrweight;
extern int disopt;
extern int pamN;
extern int checkC;
extern float geta2;
extern int treemethod;
extern int kimuraR;
extern char *swopt;
extern int fftkeika;
extern int score_check;
extern char *inputfile;
extern char *addfile;
extern int addprofile;
extern float consweight_multi;
extern float consweight_rna;
extern char RNAscoremtx;

extern char *signalSM;
extern FILE *prep_g;
extern FILE *trap_g;
extern char **seq_g;
extern char **res_g;
extern int rnakozo;
extern char rnaprediction;

/* sengen no ichi ha koko dake de ha nai */
extern void constants();
extern char **Calignm1();
extern char **Dalignm1();
extern char **align0();
extern double Cscore_m_1( char **, int, int, double ** );
extern double score_m_1(  char **, int, int, double ** );
extern double score_calc0( char **, int, double **, int );
extern char seqcheck( char ** );
extern float substitution( char *, char * );
extern float substitution_score( char *, char * );
extern float substitution_nid( char *, char * );
extern float substitution_hosei( char *, char * );
extern double ipower( double, int );
extern float translate_and_Calign();
extern float A__align();
extern float A__align11();
extern float A__align_gapmap();
extern float partA__align();
extern float L__align11( double **scoringmtx, float scoreoffset, char **seq1, char **seq2, int alloclen, int *off1pt, int *off2pt );
extern float G__align11();
extern float Falign();
extern float Falign_localhom();
extern float Conalign();
extern float Aalign();
extern float imp_match_out_sc( int, int );
extern float part_imp_match_out_sc( int, int );
extern void ErrorExit();
extern void cpmx_calc();
extern void intergroup_score( char **, char **, double *, double *, int, int, int, double * );
extern int conjuctionfortbfast();
extern int fastconjuction();
extern char seqcheck( char ** );

typedef struct _LocalHom
{
	int nokori;
	struct _LocalHom *next;
	struct _LocalHom *last;
	int start1;
	int end1;
	int start2;
	int end2;
	double opt;
	int overlapaa;
	int extended;
	double importance;
	float fimportance;
	double wimportance;
	char korh;
} LocalHom;

typedef struct _NodeInCub
{
	int step;
	int LorR;
} NodeInCub;

typedef struct _Node
{
	struct _Node *children[3];
	int tmpChildren[3];
	double length[3];
	double *weightptr[3];
	int top[3];
	int *members[3];
} Node;

typedef struct _Segment
{
	int start;
	int end;
	int center;
	double score;
	int skipForeward;
	int skipBackward;
	struct _Segment *pair;
	int  number;
} Segment;

typedef struct _Segments
{
	Segment group1;
	Segment group2;
	int number1;
	int number2;
} Segments;

typedef struct _Bchain
{
	struct _Bchain *next;
	struct _Bchain *prev;
	int pos;
} Bchain;

typedef struct _Achain
{
	int next;
	int prev;
//	int curr;
} Achain;


typedef struct _Fukusosuu
{
    double R;
    double I;
} Fukusosuu;

typedef struct _Gappattern
{
	int len;
	float freq;
} Gappat;

typedef struct _RNApair
{
	int uppos;
	float upscore;
	int downpos;
	float downscore;
	int bestpos;
	float bestscore;
} RNApair;

typedef struct _Treedep
{
	int child0;
	int child1;
	int done;
	float distfromtip;
} Treedep;

typedef struct _Addtree
{
	int nearest;
	float dist1;
	char *neighbors;
	float dist2;
} Addtree;

#include "fft.h"
#include "dp.h"
#include "functions.h"

#ifdef enablemultithread
#define TLS __thread
#else
#define TLS 
#endif

extern TLS int commonAlloc1;
extern TLS int commonAlloc2;
extern TLS int **commonIP;
extern TLS int **commonJP;
extern int nthread;
extern int randomseed;
extern int parallelizationstrategy;
#define BESTFIRST 0
#define BAATARI0 1
#define BAATARI1 2
#define BAATARI2 3

extern int scoreout;
extern int outnumber;

extern int legacygapcost;

extern TLS char *newgapstr;

extern int nalphabets;
extern int nscoredalphabets;
extern double specificityconsideration;
extern int ndistclass, maxdistclass;
extern char **gmsg;
extern int gmsglen;
