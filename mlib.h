
#ifndef MLIB_H
#define MLIB_H


float true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength);

float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues);

void cs_computations(struct phase_cs_in  p_in, struct phase_cs_out *p_out );

void sym_comp(struct phase_cs_in pa, struct phase_cs_in pb,struct phase_cs_in pc,struct sym_out *sym);

void sym_mag(struct sym_out sym, struct sym_out *sym_back, struct sym_out *sym_rms);

void signal_spectra(
	
	float rtInput, 
	struct spectra *h,
	unsigned int qBufferLength,		//updated buffer length
	float *twBufferReal,			//twiddle factor Real coeffs
	float *twBufferImag,			//twiddle factor Imag coeffs    
	unsigned int pCounter);

float signal_thd(struct spectra h);

void pvp_filter(struct pvp_data in,struct pvp_data *in_back,struct pvp_data *out, struct pvp_data *out_back,float ts);

float peak_detect_rms(float rtInput, float *pData,unsigned int pDataCounter, unsigned int dataLength);

float thermal_status(struct thermal_parameters therm, float mem);


#endif
