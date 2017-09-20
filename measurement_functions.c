

#include <math.h>



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
//

void signal_spectra(	float rtInput, 
						float *qBuffer, 				//updated buffer
						unsigned int qBufferLength, 	//updated buffer length
						float *twBufferReal,              //twiddle factor Real coeffs
						float *twBufferImag,			//twiddle factor Imag coeffs	
						float *foutReal,				//back real for hist
 						float *foutImag,				//back imag for hist
					     float *foutMag					//magnitude output
	)

{

	float x_error;
	float temp_real,temp_imag;
	float out_scale;
	static unsigned int f_count=0;
	

	out_scale=1.41421356f/(float)qBufferLength;

	//x_error=*qBuffer-rtInput;
	//*qBuffer++=rtInput;

	x_error=qBuffer[f_count]-rtInput;
	qBuffer[f_count]=rtInput;



	if(++f_count==qBufferLength){f_count=0;}

	temp_real =twBufferReal[0]* (foutReal[0]+x_error)-twBufferImag[0]*foutImag[0];
	temp_imag=twBufferImag[0]* (foutReal[0]+x_error)+twBufferReal[0]*foutImag[0];

	foutMag[0]=out_scale*sqrtf(temp_real*temp_real+temp_imag*temp_imag);

	foutReal[0]=temp_real;
	foutImag[0]=temp_imag;


}





