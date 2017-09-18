
// True RMS for distorted data
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
