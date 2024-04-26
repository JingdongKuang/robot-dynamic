//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Chebyshef_40hz_Order3.h
//
// MATLAB Coder version            : 5.6
// C/C++ source code generated on  : 14-Mar-2024 10:31:47
//

#ifndef CHEBYSHEF_40HZ_ORDER3_H
#define CHEBYSHEF_40HZ_ORDER3_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

struct dsp_BiquadFilter_0 {
	int S0_isInitialized;
	double W0_FILT_STATES[4];
	int W1_PreviousNumChannels;
	double P0_ICRTP;
	double P1_RTP1COEFF[6];
	double P2_RTP2COEFF[4];
	double P3_RTP3COEFF[3];
	boolean_T P4_RTP_COEFF3_BOOL[3];
};

namespace coder {
	namespace dspcodegen {
		class BiquadFilter {
		public:
			BiquadFilter();
			~BiquadFilter();
			boolean_T matlabCodegenIsDeleted;
			dsp_BiquadFilter_0 cSFunObject;

		private:
			boolean_T isSetupComplete;
		};

	} // namespace dspcodegen
} // namespace coder

struct Chebyshef_40hz_Order3PersistentData {
	coder::dspcodegen::BiquadFilter Hd;
};

struct Chebyshef_40hz_Order3StackData {
	Chebyshef_40hz_Order3PersistentData* pd;
};
// Type Definitions
class ChebyshefFilter {
public:
  ChebyshefFilter();
  ~ChebyshefFilter();
  double dofilter(double x);
  Chebyshef_40hz_Order3StackData *getStackData();

private:
  Chebyshef_40hz_Order3PersistentData pd_;
  Chebyshef_40hz_Order3StackData SD_;
};

#endif
//
// File trailer for Chebyshef_40hz_Order3.h
//
// [EOF]
//
