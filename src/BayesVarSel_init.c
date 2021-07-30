#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void flsBF(void *, void *, void *, void *, void *, void *);
extern void flsConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void flsSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void flsUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gBF(void *, void *, void *, void *, void *);
extern void gConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsflsConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsflsSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsflsUser(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsgConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsgSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsgUser(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsLiangConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsLiangSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsLiangUser(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsRobustConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsRobustSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsRobustUser(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsRobust2Const(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsRobust2SB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsRobust2User(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsintrinsicConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsintrinsicSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsintrinsicUser(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsZSConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsZSSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsZSUser(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void intrinsicBF(void *, void *, void *, void *, void *);
extern void LiangBF(void *, void *, void *, void *, void *);
extern void LiangConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void LiangSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void LiangUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RobustBF(void *, void *, void *, void *, void *);
extern void RobustConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RobustSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RobustUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Robust2Const(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Robust2SB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Robust2User(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void intrinsicConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void intrinsicSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void intrinsicUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geointrinsicConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geointrinsicSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geointrinsicUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geointrinsic2Const(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geointrinsic2SB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geointrinsic2User(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ZSBF(void *, void *, void *, void *, void *);
extern void ZSConst(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ZSSB(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ZSUser(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsFSBSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsFSB(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsFConstConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsFConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GibbsFSBConst(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
    {"flsBF",            (DL_FUNC) &flsBF,             6},
    {"flsConst",         (DL_FUNC) &flsConst,          9},
    {"flsSB",            (DL_FUNC) &flsSB,             9},
    {"flsUser",          (DL_FUNC) &flsUser,           9},
    {"gBF",              (DL_FUNC) &gBF,               5},
    {"gConst",           (DL_FUNC) &gConst,            9},
    {"GibbsflsConst",    (DL_FUNC) &GibbsflsConst,    10},
    {"GibbsflsSB",       (DL_FUNC) &GibbsflsSB,       10},
    {"GibbsflsUser",     (DL_FUNC) &GibbsflsUser,     10},
    {"GibbsgConst",      (DL_FUNC) &GibbsgConst,      10},
    {"GibbsgSB",         (DL_FUNC) &GibbsgSB,         10},
    {"GibbsgUser",       (DL_FUNC) &GibbsgUser,       10},
    {"GibbsLiangConst",  (DL_FUNC) &GibbsLiangConst,  10},
    {"GibbsLiangSB",     (DL_FUNC) &GibbsLiangSB,     10},
    {"GibbsLiangUser",   (DL_FUNC) &GibbsLiangUser,   10},
    {"GibbsRobustConst", (DL_FUNC) &GibbsRobustConst, 10},
    {"GibbsRobustSB",    (DL_FUNC) &GibbsRobustSB,    10},
    {"GibbsRobustUser",  (DL_FUNC) &GibbsRobustUser,  10},
    {"GibbsRobust2Const",(DL_FUNC) &GibbsRobust2Const, 10},
    {"GibbsRobust2SB",   (DL_FUNC) &GibbsRobust2SB,    10},
    {"GibbsRobust2User", (DL_FUNC) &GibbsRobust2User,  10},
    {"GibbsintrinsicConst",(DL_FUNC) &GibbsintrinsicConst, 10},
    {"GibbsintrinsicSB",   (DL_FUNC) &GibbsintrinsicSB,    10},
    {"GibbsintrinsicUser", (DL_FUNC) &GibbsintrinsicUser,  10},				
    {"GibbsZSConst",     (DL_FUNC) &GibbsZSConst,     10},
    {"GibbsZSSB",        (DL_FUNC) &GibbsZSSB,        10},
    {"GibbsZSUser",      (DL_FUNC) &GibbsZSUser,      10},
    {"gSB",              (DL_FUNC) &gSB,               9},
    {"gUser",            (DL_FUNC) &gUser,             9},
    {"intrinsicBF",      (DL_FUNC) &intrinsicBF,       5},
    {"LiangBF",          (DL_FUNC) &LiangBF,           5},
    {"LiangConst",       (DL_FUNC) &LiangConst,        9},
    {"LiangSB",          (DL_FUNC) &LiangSB,           9},
    {"LiangUser",        (DL_FUNC) &LiangUser,         9},
    {"RobustBF",         (DL_FUNC) &RobustBF,          5},
    {"RobustConst",      (DL_FUNC) &RobustConst,       9},
    {"RobustSB",         (DL_FUNC) &RobustSB,          9},
    {"RobustUser",       (DL_FUNC) &RobustUser,        9},
    {"Robust2Const",     (DL_FUNC) &Robust2Const,      9},
    {"Robust2SB",        (DL_FUNC) &Robust2SB,         9},
    {"Robust2User",      (DL_FUNC) &Robust2User,       9},
    {"intrinsicConst",   (DL_FUNC) &intrinsicConst,    9},
    {"intrinsicSB",      (DL_FUNC) &intrinsicSB,       9},
    {"intrinsicUser",    (DL_FUNC) &intrinsicUser,     9},	
	{"geointrinsicConst",   (DL_FUNC) &geointrinsicConst,    9},
    {"geointrinsicSB",      (DL_FUNC) &geointrinsicSB,       9},
    {"geointrinsicUser",    (DL_FUNC) &geointrinsicUser,     9},	
    {"geointrinsic2Const",   (DL_FUNC) &geointrinsic2Const,    9},
    {"geointrinsic2SB",      (DL_FUNC) &geointrinsic2SB,       9},
    {"geointrinsic2User",    (DL_FUNC) &geointrinsic2User,     9},	
	{"ZSBF",             (DL_FUNC) &ZSBF,              5},
    {"ZSConst",          (DL_FUNC) &ZSConst,           9},
    {"ZSSB",             (DL_FUNC) &ZSSB,              9},
    {"ZSUser",           (DL_FUNC) &ZSUser,            9},
    {"GibbsFSBSB", (DL_FUNC) &GibbsFSBSB, 10},
    {"GibbsFSB",   (DL_FUNC) &GibbsFSB, 10},		
    {"GibbsFConstConst", (DL_FUNC) &GibbsFConstConst, 10},
    {"GibbsFSBConst", (DL_FUNC) &GibbsFSBConst, 10},			
    {"GibbsFConst",(DL_FUNC) &GibbsFConst, 10},		
		{NULL, NULL, 0}
};

void R_init_BayesVarSel(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
