//**************//
// Important   //
//  Constants //
//***********//

// Number of parallel MOPAC
// (one node operation only)
#define MAXPARA 8
// Number of cores per Q-Chem instance
#define DFTPARA 1
// Number of iterations to wait for queue to finish
#define MAX_TIME_WAIT 10000

// Safe Mode Options
#define SKIP_ISOMERS 677
#define SKIPDFT 1
#define SKIPMOPAC 1

// Screening Thresholds 
#define DFT_EMAX 25
#define LST_EMAX_TS 250  //250
#define GSM_EMAX_TS 20  //90
//#define DFT_EMAX_TS 25
#define PM6_EMAX 220     //60
#define PM6_EMAX_TS 150

// Correction for PM6 overestimating H2 stability
#define H2PENALTY 15

// Thresholds for PM6 structures
#define MMVSMOPAC 4 //at least 2

#define SOLVENTGRAD 1 //xyz force constant constraint
#define EVSORIGTOL 4.5
#define MINTSE 10

//MOPAC/DFT Settings
#define NOMOOPT 1
#define QSTART 0 //note this overrides options below
#define UNRESTRICTED 0
#define DFTB3LYP 1
#define B631GSS 0
#define LANL2DZ 1
#define BASISMIX 0
#define BASISGEN 0
#define DFT_OPT_CYCLES 250
#define CATION 0
#define TRIPLET 0

