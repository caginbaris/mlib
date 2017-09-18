// function-1
// True RMS half cycle
// delayLineArray contains half period circular data of input
// delayLineCounter global counter for true rms calculation
// length of delayLineArray - mult. inverse can be u

#include <math.h>

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

	rms = sqrt(rms_sum / arrayLength);

	return rms;

}



// function-2
// fundamental RMS full cycle

float fundamental_rms(float rtInput,float *coeff, unsigned int coeffLength, float *zValues){

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




