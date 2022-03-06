#include <math.h>
#include <string.h>

double ccfpix_c(double* M_L, double* M_H, double* X, double* THAR, const double DELTA, const int N, const int M){
  // This function calculates the cross-correlation function

  double CCFPIX;

  int COND;
  double FRACTION, PIX_INIT, PIX_END;

  double M_LLOC[N], M_HLOC[N];

  int i, j;

  // SHIFT MASK BY DELTA
  for(i=0; i<N; i++){
    M_LLOC[i] = M_L[i] + DELTA;
    M_HLOC[i] = M_H[i] + DELTA;
  }

  // I marks where we are in terms of masks
  i = 0;
  CCFPIX = 0.0;

  COND = 0;
  for(j=1; j<M-1; j++){
     PIX_INIT = 0.5*(X[j-1] + X[j]);
     PIX_END  = 0.5*(X[j] + X[j+1]);
     while ((M_HLOC[i] < PIX_INIT) && (COND == 0)){
       if (i == N-1){
         COND = 1;
       }
       if (COND == 0){
         i = i + 1;
       }
     }

     if ((PIX_END < M_HLOC[i]) && (PIX_INIT > M_LLOC[i])){
       CCFPIX = CCFPIX + THAR[j];
     }
     else if ((PIX_END < M_HLOC[i]) && (PIX_INIT < M_LLOC[i]) && (PIX_END > M_LLOC[i])){
       FRACTION = (PIX_END - M_LLOC[i]) / (PIX_END - PIX_INIT);
       CCFPIX = CCFPIX + THAR[j] * FRACTION;
     }
     else if ((PIX_END > M_HLOC[i]) && (PIX_INIT > M_LLOC[i]) && (PIX_INIT < M_HLOC[i])){
       FRACTION = (M_HLOC[i] - PIX_INIT) / (PIX_END - PIX_INIT);
       CCFPIX = CCFPIX + THAR[j] * FRACTION;
     }
     else if ((PIX_END > M_HLOC[i]) && (PIX_INIT < M_LLOC[i])){
       FRACTION = (M_HLOC[i] - M_LLOC[i]) / (PIX_END - PIX_INIT);
       CCFPIX = CCFPIX + THAR[j] * FRACTION;
     }
  }

  return CCFPIX;

}

void ccfcos_c(double* M_L, double* M_H, double* WAV, double* SPEC, double* WEIGHT, double* SN, const double V_R, const int N, const int M, double* ccfcos, double* snw){
  // This function calculates the cross-correlation function
  // of a spectrum (SPE, at wavelengths WAV) with a mask. M_L
  // and M_H contain the beginning and end of regions where the mask
  // equals 1, ordered.
  // This version uses a sine function instead of a box --> better lobes
  // The function used is sin((lambda-lambda_init)*pi/(2*half_width))

  // Speed of light, km/s
  double C=2.99792458e5;

  int COND;
  double GAMMA;
  double PIX_INIT, PIX_END;
  double M_LLOC[N], M_HLOC[N], HW[N];
  double WD, NORM, ARG;
  double CCFCOS,SNW;

  int i, j;

  // DOPPLER FACTOR
  GAMMA = sqrt(1. + V_R / C) / sqrt(1. - V_R / C);

  // DOPPLER SHIFT MASK
  for(i=0; i<N; i++){
    M_LLOC[i] = M_L[i] * GAMMA;
    M_HLOC[i] = M_H[i] * GAMMA;
    HW[i]     = 0.5*( M_HLOC[i] - M_LLOC[i] );
  }


  // I marks where we are in terms of masks
  i = 0;
  CCFCOS = 0.0;
  SNW = 0.0;
  COND = 0;

  for(j=1; j<M-1; j++){
    PIX_INIT = 0.5*(WAV[j-1] + WAV[j]);
    PIX_END  = 0.5*(WAV[j] + WAV[j+1]);
    while ((M_HLOC[i] < PIX_INIT) && (COND == 0)){
      if (i == N-1){
        COND = 1;
      }
      if (COND == 0){
        i = i + 1;
      }
    }

    NORM = HW[i];
    ARG  = M_PI / (2.0*HW[i]);

    if ((PIX_END < M_HLOC[i]) && (PIX_INIT > M_LLOC[i])){
      WD = NORM*( cos((PIX_INIT-M_LLOC[i])*ARG) - cos((PIX_END-M_LLOC[i])*ARG) );
      WD = WD / (PIX_END-PIX_INIT);
      CCFCOS = CCFCOS + WD*SPEC[j] * WEIGHT[i] * SN[j];
      SNW = SNW + WD*SN[j]*WEIGHT[i];
    }
    else if ((PIX_END < M_HLOC[i]) && (PIX_INIT < M_LLOC[i]) && (PIX_END > M_LLOC[i])){
      WD = NORM*( 1.0 - cos((PIX_END-M_LLOC[i])*ARG) );
      WD = WD / (PIX_END-PIX_INIT);
      CCFCOS = CCFCOS + SPEC[j] * WEIGHT[i] * WD * SN[j];
      SNW = SNW + WD*SN[j]*WEIGHT[i];
    }
    else if ((PIX_END > M_HLOC[i]) && (PIX_INIT > M_LLOC[i]) && (PIX_INIT < M_HLOC[i])){
      WD = NORM*( cos((PIX_INIT-M_LLOC[i])*ARG) + 1.0);
      WD = WD / (PIX_END-PIX_INIT);
      CCFCOS = CCFCOS + SPEC[j] * WEIGHT[i] * WD * SN[j];
      SNW = SNW + WD*SN[j]*WEIGHT[i];
    }
    else if ((PIX_END > M_HLOC[i]) && (PIX_INIT < M_LLOC[i])){
      WD = 2 * NORM;
      WD = WD / (PIX_END-PIX_INIT);
      CCFCOS = CCFCOS + SPEC[j] * WEIGHT[i] * WD * SN[j];
      SNW = SNW + WD*SN[j]*WEIGHT[i];
    }
  }

  memcpy(ccfcos,&CCFCOS,sizeof(CCFCOS));
  memcpy(snw,&SNW,sizeof(SNW));
}