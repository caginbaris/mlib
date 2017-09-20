float true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength);
float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues);


void signal_spectra(
	
	float rtInput, 
	float *qBuffer,	//updated buffer
	unsigned int qBufferLength,	//updated buffer length
	float *twBufferReal,	//twiddle factor Real coeffs
	float *twBufferImag,	//twiddle factor Imag coeffs    
	float *foutReal,	//back real for hist
	float *foutImag,	//back imag for hist
	float *foutMag	//magnitude output

	);
