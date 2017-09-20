float true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength);
float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues);


struct spectra {

	float qBuffer[100];
	float foutReal[13];
	float foutImag[13];
	float foutMag[13];

};



void signal_spectra(
	
	float rtInput, 
	struct spectra *h,
	unsigned int qBufferLength,	//updated buffer length
	float *twBufferReal,			//twiddle factor Real coeffs
	float *twBufferImag,			//twiddle factor Imag coeffs    
	unsigned int pCounter);
