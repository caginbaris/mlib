

#include <math.h>
#include "measurement_functions.h"



// function-1
// True RMS half cycle
// delayLineArray contains half period circular data of input
// delayLineCounter global counter for true rms calculation
// length of delayLineArray - mult. inverse can be u

float true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength)
{


	unsigned int i;
	float rms = 0, rms_sum = 0, rms_data=0;

	rms_sum = 0;

	*(delayLineArray+delayLineCounter) = rtInput;

	for (i = 0; i < arrayLength; i++)
	{
			
		rms_data=*delayLineArray ++;
		rms_sum = rms_sum + rms_data*rms_data;
		

	}

	rms = sqrtf(rms_sum / arrayLength);

	return rms;

}



// function-2
// cs element generation
//cos->c, sin->s
//fs:2.5/5 e3... and scale value should be considered

float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues){

	unsigned int i;
     float *z1_ptr,*z2_ptr,*coeff_ptr;
     float output;

	z1_ptr=zValues; 		//background data
	z2_ptr=z1_ptr; 		//data update
	coeff_ptr=coeff+coeffLength-1;  //last element


	output=(*z1_ptr++) *(*coeff_ptr--);

	for(i=2;i<coeffLength;i++){
	
	*z2_ptr++ =*z1_ptr;
	output+=(*z1_ptr++) *(*coeff_ptr--);

	}


	output+=rtInput *(*coeff_ptr);
	*z2_ptr=rtInput;

	return(output);

}


// function-3
//spectral analysis
// caution: input structure has to be initialized and used only ones
// external pCounter is needed
//twfactors truncated to 13th
 


void signal_spectra(
	
	float rtInput, 
	struct spectra *h,
	unsigned int qBufferLength,	//updated buffer length
	float *twBufferReal,			//twiddle factor Real coeffs
	float *twBufferImag,			//twiddle factor Imag coeffs    
	unsigned int pCounter)

{

	float x_error;
	float temp_real,temp_imag;
	float out_scale;
	unsigned int i;

	out_scale=1.41421356f/(float)qBufferLength;


	x_error=h->qBuffer[pCounter]-rtInput;
	h->qBuffer[pCounter]=rtInput;

	for(i=0;i<13;i++){

	temp_real =twBufferReal[i]* (h->foutReal[i]+x_error)-twBufferImag[i]*h->foutImag[i];
	temp_imag=twBufferImag[i]* (h->foutReal[i]+x_error)+twBufferReal[i]*h->foutImag[i];

	h->foutMag[i]=out_scale*sqrtf(temp_real*temp_real+temp_imag*temp_imag);

	h->foutReal[i]=temp_real;
	h->foutImag[i]=temp_imag;
	
	}


}





