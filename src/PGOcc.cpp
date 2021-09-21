#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
  SEXP PGOcc(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP pocc_r, SEXP pdet_r, 
             SEXP J_r, SEXP K_r, SEXP betaStarting_r, SEXP alphaStarting_r, 
	     SEXP zStarting_r, SEXP zLongIndx_r, SEXP muBeta_r, 
	     SEXP muAlpha_r, SEXP SigmaBeta_r, SEXP SigmaAlpha_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	     SEXP nBurn_r, SEXP nThin_r, SEXP nPost_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *upper = "U";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
    char const *lside = "L";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xp = REAL(Xp_r);
    double *muBeta = REAL(muBeta_r); 
    double *muAlpha = REAL(muAlpha_r); 
    double *SigmaBetaInv = REAL(SigmaBeta_r); 
    double *SigmaAlphaInv = REAL(SigmaAlpha_r); 
    int pOcc = INTEGER(pocc_r)[0];
    int pDet = INTEGER(pdet_r)[0];
    int J = INTEGER(J_r)[0];
    int *K = INTEGER(K_r); 
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int nObs = 0;
    for (j = 0; j < J; j++) {
      nObs += K[j]; 
    } // j
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nThin = INTEGER(nThin_r)[0]; 
    int nBurn = INTEGER(nBurn_r)[0]; 
    int nPost = INTEGER(nPost_r)[0]; 
    int status = 0; 
    double *z = REAL(zStarting_r); 
    int thinIndx = 0;
    int sPost = 0;  

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
      nThreads = 1;
    }
#endif
    
    /**********************************************************************
     * Print Information 
     * *******************************************************************/
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Occupancy model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
      Rprintf("Number of MCMC samples: %i \n", nSamples);
      Rprintf("Burn-in: %i \n", nBurn); 
      Rprintf("Thinning Rate: %i \n", nThin); 
      Rprintf("Total Posterior Samples: %i \n\n", nPost); 
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      Rprintf("Sampling ... \n");
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Occupancy covariates
    double *beta = (double *) R_alloc(pOcc, sizeof(double));   
    F77_NAME(dcopy)(&pOcc, REAL(betaStarting_r), &inc, beta, &inc);
    // Detection covariates
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Auxiliary variables
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double));
    double *omegaOcc = (double *) R_alloc(J, sizeof(double));
    double *kappaDet = (double *) R_alloc(nObs, sizeof(double)); 
    double *kappaOcc = (double *) R_alloc(J, sizeof(double)); 

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pOcc, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP zSamples_r; 
    PROTECT(zSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    SEXP psiSamples_r; 
    PROTECT(psiSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(INTSXP, nObs, nPost)); nProtect++; 
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int ppDet = pDet * pDet;
    int ppOcc = pOcc * pOcc; 
    int JpOcc = J * pOcc; 
    int nObspDet = nObs * pDet;
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppOcc = (double *) R_alloc(ppOcc, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc = (double *) R_alloc(pOcc, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pOcc2 = (double *) R_alloc(pOcc, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = 0; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpOcc = (double *) R_alloc(JpOcc, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
   
    // For latent occupancy
    double psiNum; 
    double psiNew; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); 
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); 
    int *yRep = (int *) R_alloc(nObs, sizeof(int)); 

    // For normal priors
    // Occupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pOcc, SigmaBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pOcc, SigmaBetaInv, &pOcc, &info); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pOcc, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pOcc, &pOcc, &one, SigmaBetaInv, &pOcc, muBeta, &inc, &zero, SigmaBetaInvMuBeta, &inc); 	  
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, SigmaAlphaInvMuAlpha, &inc); 	  


    GetRNGstate(); 

    for (s = 0; s < nSamples; s++) {
      /********************************************************************
       *Update Occupancy Auxiliary Variables 
       *******************************************************************/
      for (j = 0; j < J; j++) {
        omegaOcc[j] = rpg(1.0, F77_NAME(ddot)(&pOcc, &X[j], &J, beta, &inc));
      } // j

      /********************************************************************
       *Update Detection Auxiliary Variables 
       *******************************************************************/
      for (i = 0; i < nObs; i++) {
        omegaDet[i] = rpg(1.0, F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc));
      } // i

      /********************************************************************
       *Update Occupancy Regression Coefficients
       *******************************************************************/
      for (j = 0; j < J; j++) {
	kappaOcc[j] = z[j] - 1.0 / 2.0; 
      } // j

      /********************************
       * Compute b.beta
       *******************************/
      F77_NAME(dgemv)(ytran, &J, &pOcc, &one, X, &J, kappaOcc, &inc, &zero, tmp_pOcc, &inc); 	 
      for (j = 0; j < pOcc; j++) {
        tmp_pOcc[j] += SigmaBetaInvMuBeta[j]; 
      } // j 

      /********************************
       * Compute A.beta
       * *****************************/
      // tmp_JpOcc is X %*% omegaOcc. 
      for(j = 0; j < J; j++){
        for(i = 0; i < pOcc; i++){
          tmp_JpOcc[i*J+j] = X[i*J+j]*omegaOcc[j];
        }
      }

      // This finishes off A.beta
      // 1 * X * tmp_JpOcc + 0 * tmp_ppOcc = tmp_ppOcc
      F77_NAME(dgemm)(ytran, ntran, &pOcc, &pOcc, &J, &one, X, &J, tmp_JpOcc, &J, &zero, tmp_ppOcc, &pOcc);
      for (j = 0; j < ppOcc; j++) {
        tmp_ppOcc[j] += SigmaBetaInv[j]; 
      } // j

      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      F77_NAME(dpotri)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotri here failed\n");}
      F77_NAME(dsymv)(lower, &pOcc, &one, tmp_ppOcc, &pOcc, tmp_pOcc, &inc, &zero, tmp_pOcc2, &inc);
      F77_NAME(dpotrf)(lower, &pOcc, tmp_ppOcc, &pOcc, &info); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      mvrnorm(beta, tmp_pOcc2, tmp_ppOcc, pOcc);
      
      /********************************************************************
       *Update Detection Regression Coefficients
       *******************************************************************/
      // /********************************
      //  * Compute b.alpha
      //  *******************************/
      // First multiply kappDet * the current occupied values, such that values go 
      // to 0 if z == 0 and values go to kappaDet if z == 1
      for (i = 0; i < nObs; i++) {
        kappaDet[i] = (y[i] - 1.0/2.0) * z[zLongIndx[i]];
      } // i
      
      F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, kappaDet, &inc, &zero, tmp_pDet, &inc); 	  
      for (j = 0; j < pDet; j++) {
        tmp_pDet[j] += SigmaAlphaInvMuAlpha[j]; 
      } // j

      /********************************
       * Compute A.alpha
       * *****************************/
      for (j = 0; j < nObs; j++) {
        for (i = 0; i < pDet; i++) {
          tmp_nObspDet[i*nObs + j] = Xp[i * nObs + j] * omegaDet[j] * z[zLongIndx[j]];
        } // i
      } // j

      F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet);

      for (j = 0; j < ppDet; j++) {
        tmp_ppDet[j] += SigmaAlphaInv[j]; 
      } // j

      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
      F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
      F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc);
      // Computes cholesky of tmp_ppDet again stored back in tmp_ppDet. This chol(A.alpha.inv)
      F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info); 
      if(info != 0){error("c++ error: dpotrf here failed\n");}
      mvrnorm(alpha, tmp_pDet2, tmp_ppDet, pDet);

      /********************************************************************
       *Update Latent Occupancy
       *******************************************************************/
      // Compute detection probability 
      for (i = 0; i < nObs; i++) {
        detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc), zero, one);
        if (tmp_J[zLongIndx[i]] == 0) {
          psi[zLongIndx[i]] = logitInv(F77_NAME(ddot)(&pOcc, &X[zLongIndx[i]], &J, beta, &inc), zero, one); 
        }
        piProd[zLongIndx[i]] *= (1.0 - detProb[i]);
        ySum[zLongIndx[i]] += y[i]; 	
        tmp_J[zLongIndx[i]]++;
      } // i
      // Compute occupancy probability 
      for (j = 0; j < J; j++) {
        psiNum = psi[j] * piProd[j]; 
        if (ySum[j] == zero) {
          z[j] = rbinom(one, psiNum / (psiNum + (1.0 - psi[j])));           
        } else {
          z[j] = one; 
        }
        piProd[j] = one;
        ySum[j] = zero; 
        tmp_J[j] = 0; 
      } // j

      /********************************************************************
       *Save samples
       *******************************************************************/
      if (s >= nBurn) {
        thinIndx++; 
	if (thinIndx == nThin) {
          F77_NAME(dcopy)(&pOcc, beta, &inc, &REAL(betaSamples_r)[sPost*pOcc], &inc);
          F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
          F77_NAME(dcopy)(&J, psi, &inc, &REAL(psiSamples_r)[sPost*J], &inc); 
          F77_NAME(dcopy)(&J, z, &inc, &REAL(zSamples_r)[sPost*J], &inc); 
	  // Replicate data set for GoF
          for (i = 0; i < nObs; i++) {
            yRep[i] = rbinom(one, detProb[i] * z[zLongIndx[i]]);
            INTEGER(yRepSamples_r)[sPost * nObs + i] = yRep[i]; 
          } // i
          sPost++; 
	  thinIndx = 0; 
	}
      }

      /********************************************************************
       * Report
       *******************************************************************/
      if (status == nReport){
        if(verbose){
          Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
      	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
      	  R_FlushConsole();
          #endif
      	}
        status = 0;
      }
      
      status++;

      R_CheckUserInterrupt();
    }
    if (verbose) {
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    }
    PutRNGstate();

    SEXP result_r, resultName_r;
    int nResultListObjs = 5;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, zSamples_r); 
    SET_VECTOR_ELT(result_r, 3, psiSamples_r);
    SET_VECTOR_ELT(result_r, 4, yRepSamples_r);

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("z.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("psi.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("y.rep.samples")); 
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

