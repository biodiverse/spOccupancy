#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {
  SEXP tIntPGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, SEXP consts_r, 
                 SEXP pDetLong_r, SEXP JLong_r, SEXP nObsLong_r, 
                 SEXP nOccRELong_r, SEXP nDetRELong_r, SEXP betaStarting_r, 
                 SEXP alphaStarting_r, SEXP sigmaSqPsiStarting_r, 
                 SEXP sigmaSqPStarting_r, SEXP betaStarStarting_r, 
                 SEXP alphaStarStarting_r, SEXP zStarting_r, 
                 SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, 
                 SEXP zYearIndx_r, SEXP zDatIndx_r, 
                 SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
                 SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
                 SEXP alphaNREIndx_r, SEXP alphaColIndx_r, 
                 SEXP muBeta_r, SEXP SigmaBeta_r, 
                 SEXP muAlpha_r, SEXP sigmaAlpha_r, 
                 SEXP sigmaSqPsiA_r, SEXP sigmaSqPsiB_r, 
                 SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
                 SEXP ar1_r, SEXP ar1Vals_r, SEXP tuning_r,
                 SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r,
                 SEXP nThreads_r, SEXP verbose_r, 
                 SEXP nReport_r, SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r, 
                 SEXP chainInfo_r, SEXP waicNObsIndx_r, SEXP waicCellIndx_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, r, q, rr, ll,  t, info, nProtect=0, l;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";

    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xp = REAL(Xp_r);
    int *XpRE = INTEGER(XpRE_r); 
    int *XRE = INTEGER(XRE_r); 
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pOcc = INTEGER(consts_r)[2];
    int pOccRE = INTEGER(consts_r)[3];
    int nOccRE = INTEGER(consts_r)[4];
    int pDet = INTEGER(consts_r)[5];
    int pDetRE = INTEGER(consts_r)[6];
    int nDetRE = INTEGER(consts_r)[7];
    int nYearsMax = INTEGER(consts_r)[8];
    int nData = INTEGER(consts_r)[9];
    int nWAIC = INTEGER(consts_r)[10];
    int *pDetLong = INTEGER(pDetLong_r); 
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc;
    /**********************************
     * Priors
     * *******************************/
    double *muBeta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppOcc, sizeof(double));   
    F77_NAME(dcopy)(&ppOcc, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *sigmaAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(sigmaAlpha_r), &inc, sigmaAlpha, &inc);
    double *sigmaSqPsiA = REAL(sigmaSqPsiA_r); 
    double *sigmaSqPsiB = REAL(sigmaSqPsiB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    int *JLong = INTEGER(JLong_r);
    // Total number of sites across all data sets, including
    // sites that are sampled by multiple data sources
    int JSum = 0; 
    for (i = 0; i < nData; i++) {
      JSum += JLong[i];
    }
    int *nObsLong = INTEGER(nObsLong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nOccRELong = INTEGER(nOccRELong_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *dataIndx = INTEGER(dataIndx_r); 
    int *alphaIndx = INTEGER(alphaIndx_r); 
    int *zYearIndx = INTEGER(zYearIndx_r); 
    int *zDatIndx = INTEGER(zDatIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *alphaNREIndx = INTEGER(alphaNREIndx_r);
    int *alphaColIndx = INTEGER(alphaColIndx_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    /**********************************
     * AR1 Parameters 
     * *******************************/
    int ar1 = INTEGER(ar1_r)[0];
    double rhoA = REAL(ar1Vals_r)[0];
    double rhoB = REAL(ar1Vals_r)[1];
    double sigmaSqTA = REAL(ar1Vals_r)[2];
    double sigmaSqTB = REAL(ar1Vals_r)[3];
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    double acceptRate = REAL(acceptRate_r)[0];
    double *tuning = REAL(tuning_r); 
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int thinIndx = 0; 
    int sPost = 0; 
    int status = 0;
    // For looping through data sets
    int stNObs = 0; 
    int stAlpha = 0; 
    // For getting likelihood values for WAIC
    int *waicCellIndx = INTEGER(waicCellIndx_r);
    int *waicNObsIndx = INTEGER(waicNObsIndx_r);
  
    // Some constants
    int JnYears = J * nYearsMax;

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    /**********************************************************************
     * Print Information 
     * *******************************************************************/
    if (verbose) {
      if (currChain == 1) {
        Rprintf("----------------------------------------\n");
        Rprintf("\tModel description\n");
        Rprintf("----------------------------------------\n");
        Rprintf("Integrated multi-season Occupancy Model with Polya-Gamma latent variable\nfit with %i sites and %i primary time periods.\n\n", J, nYearsMax);
        Rprintf("Integrating %i occupancy data sets.\n\n", nData); 
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        if (ar1) {
          Rprintf("Using an AR(1) temporal autocorrelation matrix in the occurrence sub-model.\n\n");
        }
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      }
     Rprintf("----------------------------------------\n");
     Rprintf("\tChain %i\n", currChain);
     Rprintf("----------------------------------------\n");
     Rprintf("Sampling ... \n");
     #ifdef Win32
       R_FlushConsole();
     #endif
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Occurrence fixed effects
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqPsi = (double *) R_alloc(pOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&pOccRE, REAL(sigmaSqPsiStarting_r), &inc, sigmaSqPsi, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nOccRE, sizeof(double)); 
    F77_NAME(dcopy)(&nOccRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Detection fixed effects
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent Occurrence
    double *z = (double *) R_alloc(JnYears, sizeof(double));   
    F77_NAME(dcopy)(&JnYears, REAL(zStarting_r), &inc, z, &inc);
    // PG auxiliary variables 
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaOcc = (double *) R_alloc(JnYears, sizeof(double)); zeros(omegaOcc, JnYears);
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(kappaDet, nObs);
    double *kappaOcc = (double *) R_alloc(JnYears, sizeof(double)); zeros(kappaOcc, JnYears);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pOcc * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDet * nPost);
    SEXP zSamples_r; 
    PROTECT(zSamples_r = Rf_allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    zeros(REAL(zSamples_r), JnYears * nPost);
    SEXP etaSamples_r; 
    if (ar1) {
      PROTECT(etaSamples_r = Rf_allocMatrix(REALSXP, nYearsMax, nPost)); nProtect++; 
      zeros(REAL(etaSamples_r), nYearsMax * nPost);
    }
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = Rf_allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    zeros(REAL(psiSamples_r), JnYears * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = Rf_allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = Rf_allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetRE * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqPsiSamples_r; 
    SEXP betaStarSamples_r; 
    if (pOccRE > 0) {
      PROTECT(sigmaSqPsiSamples_r = Rf_allocMatrix(REALSXP, pOccRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPsiSamples_r), pOccRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nOccRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nOccRE * nPost);
    }
    SEXP pSamples_r; 
    PROTECT(pSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
    zeros(REAL(pSamples_r), nObs * nPost);
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, nWAIC, nPost)); nProtect++;
    zeros(REAL(likeSamples_r), nWAIC * nPost);
    
    /**********************************************************************
     * Some constants and temporary variables to be used later
     * *******************************************************************/
    int JnYearspOccRE = J * nYearsMax * pOccRE; 
    int nObspDet = nObs * pDet;
    int nObspDetRE = nObs * pDetRE;
    int nnYears = nYearsMax * nYearsMax;
    double tmp_0 = 0.0;
    double tmp_02 = 0.0;
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_ppOcc2 = (double *) R_alloc(ppOcc, sizeof(double));
    zeros(tmp_ppOcc2, ppOcc);
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_one = (double *) R_alloc(1, sizeof(double)); 
    int *tmp_JnYearsInt = (int *) R_alloc(JnYears, sizeof(int));
    for (j = 0; j < JnYears; j++) {
      tmp_JnYearsInt[j] = zero; 
    }
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J1, J);
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double));
    zeros(tmp_nObs, nObs);
    double *tmp_JnYears = (double *) R_alloc(JnYears, sizeof(double));
    zeros(tmp_JnYears, JnYears);
    double *tmp_JnYearspOcc = (double *) R_alloc(JnYears * pOcc, sizeof(double));
    zeros(tmp_JnYearspOcc, JnYears * pOcc);

    // For latent occupancy
    double psiNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(JnYears, sizeof(double)); 
    zeros(psi, JnYears); 
    double *piProd = (double *) R_alloc(JnYears, sizeof(double)); 
    ones(piProd, JnYears); 
    double *ySum = (double *) R_alloc(JnYears, sizeof(double)); zeros(ySum, JnYears);
    // Stuff for WAIC
    double *yWAIC = (double *) R_alloc(nWAIC, sizeof(double)); ones(yWAIC, nWAIC);
    double *piProdWAIC = (double *) R_alloc(nWAIC, sizeof(double)); 
    ones(piProdWAIC, nWAIC); 
    double *yWAICSum = (double *) R_alloc(nWAIC, sizeof(double)); zeros(yWAICSum, nWAIC);

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info FCONE); 
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    // Have "separate" multivariate normal priors for the different sets of coefficients
    // that vary across the data sets. 
    // Get size of vector
    int currSize = 0; 
    // Index of starting prior values. 
    int *alphaSigmaIndx = (int *) R_alloc(nData, sizeof(int)); 
    int *alphaMuIndx = (int *) R_alloc(nData, sizeof(int)); 
    int tmp0 = 0; 
    int tmp02 = 0; 
    for (q = 0; q < nData; q++) {
      currSize += pDetLong[q] * pDetLong[q];  
      alphaSigmaIndx[q] = tmp0; 
      tmp0 += pDetLong[q] * pDetLong[q]; 
      alphaMuIndx[q] = tmp02; 
      tmp02 += pDetLong[q]; 
    } // q
    double *SigmaAlphaInv = (double *) R_alloc(currSize, sizeof(double)); zeros(SigmaAlphaInv, currSize); 
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    // Fill SigmaAlpha
    for (q = 0, j = 0; q < nData; q++) {
      for (i = 0; i < pDetLong[q]; i++, j++) {
        SigmaAlphaInv[alphaSigmaIndx[q] + i * pDetLong[q] + i] = sigmaAlpha[j]; 
        // Rprintf("Index: %i\n", alphaSigmaIndx[q] + i * pDetLong[q] + i); 
      } // i
      F77_NAME(dpotrf)(lower, &pDetLong[q], &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaInv failed\n");}
      F77_NAME(dpotri)(lower, &pDetLong[q], &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
      if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaInv failed\n");}
      F77_NAME(dsymv)(lower, &pDetLong[q], &one, &SigmaAlphaInv[alphaSigmaIndx[q]], &pDetLong[q], &muAlpha[alphaMuIndx[q]], &inc, &zero, 
                      &SigmaAlphaInvMuAlpha[alphaMuIndx[q]], &inc FCONE);
    } // q
   

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site/year-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JnYears, sizeof(double)); 
    zeros(betaStarSites, JnYears); 
    int *betaStarLongIndx = (int *) R_alloc(JnYearspOccRE, sizeof(int));
    // Initial sums
    for (t = 0; t < nYearsMax; t++) {
      for (j = 0; j < J; j++) {
        for (l = 0; l < pOccRE; l++) {
          betaStarLongIndx[l * JnYears + t * J + j] = which(XRE[l * JnYears + t * J + j], 
                                                            betaLevelIndx, nOccRE);
          betaStarSites[t * J + j] += betaStar[betaStarLongIndx[l * JnYears + t * J + j]];
        } // l
      } // j
    } // t
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    int *alphaStarLongIndx = (int *) R_alloc(nObspDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < alphaNREIndx[i]; l++) {
        alphaStarLongIndx[l * nObs + i] = which(XpRE[l * nObs + i], alphaLevelIndx, nDetRE);
        alphaStarObs[i] += alphaStar[alphaStarLongIndx[l * nObs + i]];
      }
    }
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pOccRE, sizeof(int)); 
    for (l = 0; l < pOccRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nOccRE); 
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
    }

    /**********************************************************************
     * Set up AR1 and MH stuff
     *********************************************************************/
    int nTheta = 2, sigmaSqTIndx = 0, rhoIndx = 1;
    double *accept = (double *) R_alloc(nTheta, sizeof(double)); 
    zeros(accept, nTheta); 
    double *accept2 = (double *) R_alloc(nTheta, sizeof(double)); 
    zeros(accept2, nTheta); 
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    double logMHRatio, logPostCurr = 0.0, logPostCand = 0.0, detCand = 0.0, detCurr = 0.0;
    double rhoCand = 0.0; 
    SEXP thetaSamples_r; 
    if (ar1) {
      PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nTheta, nPost)); nProtect++; 
      zeros(REAL(thetaSamples_r), nTheta * nPost);
    }
    // Initiate values
    theta[sigmaSqTIndx] = REAL(ar1Vals_r)[5];
    theta[rhoIndx] = REAL(ar1Vals_r)[4];
    double rho = theta[rhoIndx];
    double *SigmaEta = (double *) R_alloc(nnYears, sizeof(double));
    double *SigmaEtaCand = (double *) R_alloc(nnYears, sizeof(double));
    double *tmp_nYearsMax = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nYearsMax2 = (double *) R_alloc(nYearsMax, sizeof(double));
    double *tmp_nnYears = (double *) R_alloc(nnYears, sizeof(double));
    if (ar1) {
      AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
      // expCov(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
      clearUT(SigmaEta, nYearsMax);
      F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){Rf_error("c++ error: Cholesky failed in initial time covariance matrix\n");}
      F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
      if(info != 0){Rf_error("c++ error: Cholesky inverse failed in initial time covariance matrix\n");}
    }
    double *eta = (double *) R_alloc(nYearsMax, sizeof(double)); zeros(eta, nYearsMax);
    // For sigmaSqT sampler
    double aSigmaSqTPost = 0.5 * nYearsMax + sigmaSqTA;
    double bSigmaSqTPost = 0.0;
    double *etaTRInv = (double *) R_alloc(nYearsMax, sizeof(double));

    GetRNGstate();

    /**********************************************************************
     * Begin Sampler 
     * *******************************************************************/
    for (s = 0, ll = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, ll++) {
        /********************************************************************
         *Update Occupancy Auxiliary Variables 
         *******************************************************************/
        zeros(tmp_JnYears, JnYears);
        for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            // Only calculate omegaOcc when there are observations at that 
            // site/year combo. Otherwise, you're just wasting time. 
            if (zDatIndx[t * J + j] == 1) { 
              tmp_JnYears[t * J + j] = F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, 
                 	                      beta, &inc);
              omegaOcc[t * J + j] = rpg(1.0, tmp_JnYears[t * J + j] + betaStarSites[t * J + j] + 
                                        eta[t]);
              // Update kappa values along the way. 
              kappaOcc[t * J + j] = z[t * J + j] - 1.0 / 2.0; 
            }
          } // j 
        } // t
        zeros(kappaDet, nObs);
        /********************************************************************
         *Update Detection Auxiliary Variables 
         *******************************************************************/
        // Note that all of the variables are sampled, but only those at 
        // locations with z[j] == 1 actually effect the results. 
        for (i = 0; i < nObs; i++) {
          stAlpha = which(dataIndx[i], alphaIndx, pDet); 
          // If the site is occupied and data were collected at that site. 
          if (z[zLongIndx[i]] == 1.0 && zDatIndx[zLongIndx[i]] == 1) {
            omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, 
                                            &alpha[stAlpha], &inc) + alphaStarObs[i]);
          }
        } // i
        
        /********************************************************************
         *Update Occupancy Regression Coefficients
         *******************************************************************/
        zeros(tmp_JnYears, JnYears);
        for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            if (zDatIndx[t * J + j] == 1) {
              tmp_JnYears[t * J + j] = kappaOcc[t * J + j] - omegaOcc[t * J + j] * (betaStarSites[t * J + j] + eta[t]); 
            }
          } // j
        } // t
        /********************************
         * Compute b.beta
         *******************************/
        // This is fine, because the elements in tmp_JnYears corresponding
        // to unobserve site/time locations is set to 0 and not changed. 
        F77_NAME(dgemv)(ytran, &JnYears, &pOcc, &one, X, &JnYears, 
        		tmp_JnYears, &inc, &zero, tmp_pOcc, &inc FCONE); 	 
        for (j = 0; j < pOcc; j++) {
          tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
        // This is fine, because omegaOcc == 0 for the site/year combos
        // that don't have any observations at them, which will cause this
        // whole product to go to 0.  
        for (t = 0; t < nYearsMax; t++) {
          for (j = 0; j < J; j++) {
            if (zDatIndx[t * J + j] == 1) {
              for(i = 0; i < pOcc; i++){
                tmp_JnYearspOcc[i * JnYears + t * J + j] = X[i * JnYears + t * J + j] * omegaOcc[t * J + j];
              } // i
            }
          } // j
        } // t

        F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &JnYears, &one, X, &JnYears, tmp_JnYearspOcc, &JnYears, &zero, tmp_ppOcc, &pOcc FCONE FCONE);
        for (j = 0; j < ppOcc; j++) {
          tmp_ppOcc[j] += SigmaBetaInv[j]; 
        } // j


        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf A.beta failed\n");}
        F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotri A.beta failed\n");}
        F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info FCONE); 
        if(info != 0){Rf_error("c++ error: dpotrf A.beta2 failed\n");}
        mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
        
        /********************************************************************
         *Update Detection Regression Coefficients
         *******************************************************************/
        for (q = 0; q < nData; q++) {
          zeros(tmp_nObs, nObs);
          // Starting locations
          stNObs = which(q, dataIndx, nObs); 
          stAlpha = which(q, alphaIndx, pDet); 
          /********************************
           * Compute b.alpha
           *******************************/
          // First multiply kappaDet * the current occupied values, such that values go 
          // to 0 if z == 0 and values go to kappaDet if z == 1
          for (i = 0; i < nObsLong[q]; i++) {
            // 1.0 is currently hardcoded in for occupancy data
            kappaDet[stNObs + i] = (y[stNObs + i] - 1.0/2.0) * z[zLongIndx[stNObs + i]];
            tmp_nObs[stNObs + i] = kappaDet[stNObs + i] - omegaDet[stNObs + i] * alphaStarObs[stNObs + i];
            tmp_nObs[stNObs + i] *= z[zLongIndx[stNObs + i]];
          } // i
          // Xp * kappaDet + 0 * tmp_pDet. Output is stored in tmp_pDet
          F77_NAME(dgemv)(ytran, &nObsLong[q], &pDetLong[q], &one, &Xp[stNObs], &nObs, 
                          &tmp_nObs[stNObs], &inc, &zero, &tmp_pDet[stAlpha], &inc FCONE); 	  
          for (j = 0; j < pDetLong[q]; j++) {
            tmp_pDet[stAlpha + j] += SigmaAlphaInvMuAlpha[stAlpha + j]; 
          } // j
          
          /********************************
           * Compute A.alpha
           * *****************************/
          for (j = 0; j < nObsLong[q]; j++) {
            for (i = 0; i < pDetLong[q]; i++) {
              tmp_nObspDet[stNObs + i*nObs + j] = Xp[stNObs + i * nObs + j] * omegaDet[stNObs + j] * z[zLongIndx[stNObs + j]];
            } // i
          } // j
          
          // This finishes off A.alpha
          // 1 * Xp * tmp_nObspDet + 0 * tmp_ppDet = tmp_ppDet
          F77_NAME(dgemm)(ytran, ntran, &pDetLong[q], &pDetLong[q], &nObsLong[q], &one, &Xp[stNObs], 
                          &nObs, &tmp_nObspDet[stNObs], &nObs, &zero, &tmp_ppDet[alphaSigmaIndx[q]], 
                          &pDetLong[q] FCONE FCONE);
          
          for (j = 0; j < pDetLong[q] * pDetLong[q]; j++) {
            tmp_ppDet[alphaSigmaIndx[q] + j] += SigmaAlphaInv[alphaSigmaIndx[q] + j]; 
            // Rprintf("tmp_ppDet: %f\n", tmp_ppDet[alphaSigmaIndx[q] + j]); 
          } // j
          
          // This gives the Cholesky of A.alpha
          // Computes cholesky of tmp_ppDet. Output stored in tmp_ppOcc
          F77_NAME(dpotrf)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf A.alpha failed\n");}
          // Computes the inverse tmp_ppOcc. Stored in tmp_ppOcc. This is A.beta.inv. 
          F77_NAME(dpotri)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri A.alpha failed\n");}
          // A.alpha.inv %*% b.alpha
          // 1 * tmp_ppDet * tmp_pDet + 0 * tmp_pDet2 
          // (which is currently nothing) = tmp_pDet2
          F77_NAME(dsymv)(lower, &pDetLong[q], &one, &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &tmp_pDet[stAlpha], &inc, &zero, &tmp_pDet2[stAlpha], &inc FCONE);
          // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
          F77_NAME(dpotrf)(lower, &pDetLong[q], &tmp_ppDet[alphaSigmaIndx[q]], &pDetLong[q], &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(&alpha[stAlpha], &tmp_pDet2[stAlpha], &tmp_ppDet[alphaSigmaIndx[q]], pDetLong[q]);
        } // q

        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pOccRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nOccRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqPsi[l] = rigamma(sigmaSqPsiA[l] + nOccRELong[l] / 2.0, sigmaSqPsiB[l] + tmp_0); 
        }

        /********************************************************************
         *Update Detection random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nDetRELong[l], &alphaStar[alphaStarStart[l]], &inc, &alphaStar[alphaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] / 2.0, sigmaSqPB[l] + tmp_0); 
        }

        /********************************************************************
         *Update Occupancy random effects
         *******************************************************************/
        if (pOccRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nOccRE; l++) {
            /********************************
             * Compute b.beta.star
             *******************************/
            zeros(tmp_one, inc);
            tmp_0 = 0.0;	      
            // Only allow information to come from when XRE == betaLevelIndx[l]. 
            // aka information only comes from the sites with any given level 
            // of a random effect. 
            for (t = 0; t < nYearsMax; t++) {
              for (j = 0; j < J; j++) {
                if (XRE[betaStarIndx[l] * JnYears + t * J + j] == betaLevelIndx[l]) {
                  if (zDatIndx[t * J + j] == 1) {
                    tmp_02 = 0.0;
                    for (rr = 0; rr < pOccRE; rr++) {
                      tmp_02 += betaStar[betaStarLongIndx[rr * JnYears + t * J + j]];
                      } 
                    tmp_one[0] += kappaOcc[t * J + j] - (F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, beta, &inc) + 
                    	    tmp_02 - betaStar[l] + eta[t]) * omegaOcc[t * J + j];
                    tmp_0 += omegaOcc[t * J + j];
                    }
                }
              }
            }
            /********************************
             * Compute A.beta.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqPsi[betaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            betaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
        
          // Update the RE sums for the current species
          zeros(betaStarSites, JnYears);
          for (t = 0; t < nYearsMax; t++) {
            for (j = 0; j < J; j++) {
              for (l = 0; l < pOccRE; l++) {
                betaStarSites[t * J + j] += betaStar[betaStarLongIndx[l * JnYears + t * J + j]];
              }
            }
          }
        }
        
        /********************************************************************
         *Update Detection random effects
         *******************************************************************/
        if (pDetRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nDetRE; l++) {
            /********************************
             * Compute b.alpha.star
             *******************************/
            // Only allow information to come from when z[r] == 1 and XpRE == alphaLevelIndx[l]
            zeros(tmp_one, inc);
            tmp_0 = 0.0;
            for (i = 0; i < nObs; i++) {
              stAlpha = which(dataIndx[i], alphaIndx, pDet); 
              if ((z[zLongIndx[i]] == 1.0) && (XpRE[alphaColIndx[l] * nObs + i] == alphaLevelIndx[l])) {
                tmp_02 = 0.0;
                for (rr = 0; rr < alphaNREIndx[i]; rr++) {
                  tmp_02 += alphaStar[alphaStarLongIndx[rr * nObs + i]];
                }
                tmp_one[0] += kappaDet[i] - (F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc) + tmp_02 - alphaStar[l]) * omegaDet[i];
                tmp_0 += omegaDet[i];
              }
            }
            /********************************
             * Compute A.alpha.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            alphaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
          zeros(alphaStarObs, nObs); 
          // Update the RE sums
          for (i = 0; i < nObs; i++) {
            for (l = 0; l < alphaNREIndx[i]; l++) {
              alphaStarObs[i] += alphaStar[alphaStarLongIndx[l * nObs + i]]; 
            }
          }
        }

        if (ar1) {
          /********************************************************************
           *Update sigmaSqT
           *******************************************************************/
          // Form correlation matrix. 
          AR1(nYearsMax, theta[rhoIndx], 1.0, SigmaEta);
          // expCov(nYearsMax, theta[rhoIndx], 1.0, SigmaEta);
          clearUT(SigmaEta, nYearsMax);
          F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
          F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
          fillUTri(SigmaEta, nYearsMax);
          // Compute t(eta) %*% SigmaEta^-1 %*% eta
          for (t = 0; t < nYearsMax; t++) {
            etaTRInv[t] = F77_NAME(ddot)(&nYearsMax, &SigmaEta[t], &nYearsMax, 
                                   eta, &inc); 
          }
          bSigmaSqTPost = F77_NAME(ddot)(&nYearsMax, etaTRInv, &inc, eta, &inc);	
          bSigmaSqTPost /= 2.0;
          bSigmaSqTPost += sigmaSqTB;
          theta[sigmaSqTIndx] = rigamma(aSigmaSqTPost, bSigmaSqTPost);
          
          /********************************************************************
           *Update rho
           *******************************************************************/
          rho = theta[rhoIndx];
          rhoCand = logitInv(rnorm(logit(rho, rhoA, rhoB), exp(tuning[rhoIndx])), rhoA, rhoB); 
          theta[rhoIndx] = rhoCand; 

          // Construct proposal covariance matrix. 
          AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEtaCand);
          // expCov(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEtaCand);
          clearUT(SigmaEtaCand, nYearsMax);

          /********************************
           * Proposal
           *******************************/
          // Invert SigmaEtaCand and log det cov. 
          detCand = 0.0;
          F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky failed in proposal covariance matrix\n");}
          // Get log of the determinant of the covariance matrix. 
          for (k = 0; k < nYearsMax; k++) {
            detCand += 2.0 * log(SigmaEtaCand[k*nYearsMax+k]);
          } // k
          F77_NAME(dpotri)(lower, &nYearsMax, SigmaEtaCand, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky inverse failed in proposal covariance matrix\n");}
          logPostCand = 0.0; 
          // Jacobian and Uniform prior. 
          logPostCand += log(rhoCand - rhoA) + log(rhoB - rhoCand); 
          F77_NAME(dsymv)(lower, &nYearsMax, &one,  SigmaEtaCand, &nYearsMax, eta, &inc, &zero, tmp_nYearsMax, &inc FCONE);
          logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&nYearsMax, eta, &inc, tmp_nYearsMax, &inc);
          /********************************
           * Current
           *******************************/
          theta[rhoIndx] = rho; 
          AR1(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
          // expCov(nYearsMax, theta[rhoIndx], theta[sigmaSqTIndx], SigmaEta);
          clearUT(SigmaEta, nYearsMax);
          detCurr = 0.0;
          F77_NAME(dpotrf)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky failed in covariance matrix\n");}
          for (k = 0; k < nYearsMax; k++) {
            detCurr += 2.0 * log(SigmaEta[k*nYearsMax+k]);
          } // k
          F77_NAME(dpotri)(lower, &nYearsMax, SigmaEta, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: Cholesky inverse failed in covariance matrix\n");}
          logPostCurr = 0.0; 
          logPostCurr += log(rho - rhoA) + log(rhoB - rho); 
          // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
          F77_NAME(dsymv)(lower, &nYearsMax, &one, SigmaEta, &nYearsMax, eta, &inc, &zero, 
                   tmp_nYearsMax, &inc FCONE);
          logPostCurr += -0.5*detCurr-0.5*F77_NAME(ddot)(&nYearsMax, eta, &inc, tmp_nYearsMax, &inc);
          // MH Accept/Reject
          logMHRatio = logPostCand - logPostCurr; 
          if (runif(0.0, 1.0) <= exp(logMHRatio)) {
            theta[rhoIndx] = rhoCand;
            accept[rhoIndx]++;
            F77_NAME(dcopy)(&nnYears, SigmaEtaCand, &inc, SigmaEta, &inc); 
          }
          /********************************************************************
           *Update eta 
           *******************************************************************/
          /********************************
           * Compute b.w
           *******************************/
          zeros(tmp_nYearsMax, nYearsMax);
          for (j = 0; j < J; j++) {
            for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_nYearsMax[t] += kappaOcc[t * J + j] - omegaOcc[t * J + j] * (F77_NAME(ddot)(&pOcc, &X[t * J + j], &JnYears, beta, &inc) + betaStarSites[t * J + j]);
                }
              }
            }
          /********************************
           * Compute A.w
           *******************************/
          // Copy inverse covariance matrix into tmp_JJ
          F77_NAME(dcopy)(&nnYears, SigmaEta, &inc, tmp_nnYears, &inc); 
          for (j = 0; j < J; j++) {
            for (t = 0; t < nYearsMax; t++) {
              if (zDatIndx[t * J + j] == 1) {
                tmp_nnYears[t * nYearsMax + t] += omegaOcc[t * J + j];
                }
              } // t
            } // j
        
          // Cholesky of A.eta
          F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf on A.eta failed\n");}
          // Inverse of A.eta
          F77_NAME(dpotri)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotri on A.eta failed\n");}
          // A.eta.inv %*% b.eta. Stored in tmp_
          F77_NAME(dsymv)(lower, &nYearsMax, &one, tmp_nnYears, &nYearsMax, 
                   tmp_nYearsMax, &inc, &zero, tmp_nYearsMax2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &nYearsMax, tmp_nnYears, &nYearsMax, &info FCONE); 
          if(info != 0){Rf_error("c++ error: dpotrf on A.eta failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(eta, tmp_nYearsMax2, tmp_nnYears, nYearsMax);
        }

        /********************************************************************
         *Update Latent Occupancy
         *******************************************************************/
        // Compute detection probability 
        for (i = 0; i < nObs; i++) {
          stAlpha = which(dataIndx[i], alphaIndx, pDet); 
          detProb[i] = logitInv(F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, &alpha[stAlpha], &inc) + alphaStarObs[i], zero, one);
          if (tmp_JnYearsInt[zLongIndx[i]] == 0) {
            psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &JnYears, 
                                                        beta, &inc) + betaStarSites[zLongIndx[i]] + 
          			                         eta[zYearIndx[zLongIndx[i]]], zero, one); 
          }
          piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
          ySum[zLongIndx[i]] += y[i]; 	
          tmp_JnYearsInt[zLongIndx[i]]++;
          // Update stuff for WAIC
          piProdWAIC[waicNObsIndx[i]] *= pow(detProb[i], y[i]);
          piProdWAIC[waicNObsIndx[i]] *= pow(1.0 - detProb[i], 1 - y[i]);
          yWAICSum[waicNObsIndx[i]] += y[i];
        } // i
        for (t = 0; t < nYearsMax; t++) {
          // Compute occupancy probability 
          for (j = 0; j < J; j++) {
            psiNum = psi[t * J + j] * piProd[t * J + j]; 
            // If the site j was sampled in year t
            if (zDatIndx[t * J + j] == 1) {
              if (ySum[t * J + j] == zero) {
                z[t * J + j] = rbinom(one, psiNum / (psiNum + (1.0 - psi[t * J + j])));
              } else {
                z[t * J + j] = one; 
              }
            } else {
              psi[t * J + j] = logitInv(F77_NAME(ddot)(&pOcc, &X[t * J + j], 
                                                       &JnYears, beta, &inc) + 
                                                       betaStarSites[t * J + j] + eta[t], 
                                        zero, one); 
              z[t * J + j] = rbinom(one, psi[t * J + j]);
            }
            piProd[t * J + j] = one;
            ySum[t * J + j] = zero; 
            tmp_JnYearsInt[t * J + j] = 0; 
          } // j
        } // t
        // Calculating WAIC
        for (j = 0; j < nWAIC; j++) {
          if (yWAICSum[j] == zero) {
            yWAIC[j] = (1.0 - psi[waicCellIndx[j]]) + psi[waicCellIndx[j]] * piProdWAIC[j];
          } else {
            yWAIC[j] = psi[waicCellIndx[j]] * piProdWAIC[j];
          }
          piProdWAIC[j] = one;
          yWAICSum[j] = zero;
        } // j
            

        /********************************************************************
         *Save samples
         *******************************************************************/
        if (ll >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
            F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&JnYears, psi, &inc, &REAL(psiSamples_r)[sPost*JnYears], &inc); 
            F77_NAME(dcopy)(&JnYears, z, &inc, &REAL(zSamples_r)[sPost*JnYears], &inc); 
            if (ar1) {
              F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[sPost*nTheta], &inc); 
              F77_NAME(dcopy)(&nYearsMax, eta, &inc, &REAL(etaSamples_r)[sPost*nYearsMax], &inc);
            }
            if (pOccRE > 0) {
              F77_NAME(dcopy)(&pOccRE, sigmaSqPsi, &inc, 
                              &REAL(sigmaSqPsiSamples_r)[sPost*pOccRE], &inc);
              F77_NAME(dcopy)(&nOccRE, betaStar, &inc, 
                              &REAL(betaStarSamples_r)[sPost*nOccRE], &inc);
            }
            if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, 
                              &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetRE, alphaStar, &inc, 
                              &REAL(alphaStarSamples_r)[sPost*nDetRE], &inc);
            }
            for (i = 0; i < nObs; i++) {
              REAL(pSamples_r)[sPost * nObs + i] = detProb[i]; 
            } // i
            F77_NAME(dcopy)(&nWAIC, yWAIC, &inc, &REAL(likeSamples_r)[sPost*nWAIC], &inc);
            sPost++; 
            thinIndx = 0; 
          }
        }
        R_CheckUserInterrupt();
      } // r (end batch)
      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      if (ar1) {
        for (k = 0; k < nTheta; k++) {
          accept2[k] = accept[k] / batchLength;
          if (accept[k] / batchLength > acceptRate) {
            tuning[k] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          } else{
              tuning[k] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            }
          accept[k] = 0;
        }
      }
      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          if (ar1) {
            Rprintf("\tParameter\tAcceptance\tTuning\n");	  
            Rprintf("\trho\t\t%3.1f\t\t%1.5f\n", 100.0*accept2[rhoIndx], exp(tuning[rhoIndx]));
          }
          Rprintf("-------------------------------------------------\n");
          #ifdef Win32
          R_FlushConsole();
          #endif
          status = 0;
        }
        status++;        
      }
    } // all batches
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }
    PutRNGstate();
 
    SEXP result_r, resultName_r;
    int nResultListObjs = 6;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pOccRE > 0) {
      nResultListObjs += 2;
    }
    if (ar1) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, likeSamples_r);
    SET_VECTOR_ELT(result_r, 5, pSamples_r);
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 6, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 7, alphaStarSamples_r);
    }
    if (pOccRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 8;
      } else {
        tmp_0 = 6;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqPsiSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    int ar1Ind = 0;
    if (ar1) {
      if (pOccRE > 0) {
        ar1Ind = tmp_0 + 2;
      } else if (pDetRE > 0) {
        ar1Ind = 8; 
      } else {
        ar1Ind = 6;
      }
      SET_VECTOR_ELT(result_r, ar1Ind, etaSamples_r);
      SET_VECTOR_ELT(result_r, ar1Ind + 1, thetaSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("like.samples"));
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("p.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("alpha.star.samples")); 
    }
    if (pOccRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.psi.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples")); 
    }
    if (ar1) {
      SET_VECTOR_ELT(resultName_r, ar1Ind, Rf_mkChar("eta.samples")); 
      SET_VECTOR_ELT(resultName_r, ar1Ind + 1, Rf_mkChar("theta.samples")); 

    }
    
    // Set the names of the output list.  
    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

