
#define pi       3.1415926535897932f
#define _2pi  	6.2831853071795864f
#define sqrt2  1.4142135623730950f
#define _i2     0.5f


#define sym_i3 0.333333333333333f
#define sym_r  0.5f
#define sym_i  0.86602540378443864f
#define sym_rms_scale 1.0f



float true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength);
float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues);



struct phase_cs_in{

		float Vc;
		float Vs;
		float Ic;
		float Is;
};

struct phase_cs_out{

		float rms_V;
		float rms_V2;
		float rms_I;
		float rms_I2;
		float P;
		float Q;
		float X;
		float R;	
		float phase_V;
		float phase_I;
		float phase_disp;

};

void cs_computations(struct phase_cs_in  p_in, struct phase_cs_out *p_out );


struct sym_out{

	
	float V1;
	float V2;
	float V0;

	float I1;
	float I2;
	float I0;

};

void sym_comp(struct phase_cs_in pa, struct phase_cs_in pb,struct phase_cs_in pc,struct sym_out *sym);



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
