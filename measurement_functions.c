

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

//function-3
void cs_computations(struct phase_cs_in p_in, struct phase_cs_out *p_out ){

	
   	p_out->rms_V2 =(p_in.Vc)*(p_in.Vc) + (p_in.Vs)*(p_in.Vs);
    	p_out->rms_I2 =(p_in.Ic)*(p_in.Ic) + (p_in.Is)*(p_in.Is);

   	p_out->rms_V =sqrtf((p_in.Vc)*(p_in.Vc) + (p_in.Vs)*(p_in.Vs));
    	p_out->rms_I =sqrtf((p_in.Ic)*(p_in.Ic) + (p_in.Is)*(p_in.Is));

	p_out->P=(p_in.Vc)*(p_in.Ic) +(p_in.Vs)*(p_in.Is);
	p_out->Q=(p_in.Vs)*(p_in.Ic) -(p_in.Vc)*(p_in.Is);


	if(p_out->rms_I2>1.0f){

	p_out->X=2.0f*p_out->Q/p_out->rms_I2;
	p_out->R=2.0f*p_out->X/p_out->rms_I2;
	
	}


	p_out->phase_V= -atan2f(p_in.Vc,p_in.Vs)+pi;
	p_out->phase_I= -atan2f(p_in.Ic,p_in.Is)+pi;

	p_out->phase_disp=p_out->phase_I-p_out->phase_V;

	if(p_out->phase_disp<0)	{p_out->phase_disp=p_out->phase_disp+_2pi;}
	if(p_out->phase_disp>pi)	{p_out->phase_disp=p_out->phase_disp-_2pi;}

     
	


}


//function-4


void sym_comp(struct phase_cs_in pa, struct phase_cs_in pb,struct phase_cs_in pc,struct sym_out*sym){

	
	sym->V0  =(pa.Vc + pb.Vc + 	pc.Vc)*sym_3;
	sym->V1  =(pa.Vc + pb.Vc*sym_r +   pc.Vc*sym_r - pb.Vs*sym_i + pc.Vs*sym_i)*sym_3;
	sym->V2 =(pa.Vc + pb.Vc*sym_r +   pc.Vc*sym_r + pb.Vs*sym_i - pc.Vs*sym_i)*sym_3;


	sym->I0  =(pa.Ic + pb.Ic + 	pc.Ic)*sym_3;
	sym->I1  =(pa.Ic + pb.Ic*sym_r +   pc.Ic*sym_r - pb.Is*sym_i + pc.Is*sym_i)*sym_3;
	sym->I2 =(pa.Ic + pb.Ic*sym_r +   pc.Ic*sym_r + pb.Is*sym_i - pc.Is*sym_i)*sym_3;
	
	
	}


//function-5

void sym_mag(struct sym_out sym, struct sym_out *sym_back, struct sym_out *sym_rms  ){

	float temp;

	temp = (sym.V0-sym_back->V0)*sym_rms_scale;
	sym_rms->V0=sqrtf(temp*temp+sym.V0*sym.V0);
	sym_back->V0=sym.V0;

	temp = (sym.V1-sym_back->V1)*sym_rms_scale;
	sym_rms->V1=sqrtf(temp*temp+sym.V1*sym.V1);
	sym_back->V1=sym.V1;

	temp = (sym.V2-sym_back->V2)*sym_rms_scale;
	sym_rms->V2=sqrtf(temp*temp+sym.V2*sym.V2);
	sym_back->V2=sym.V2;


	temp = (sym.I0-sym_back->I0)*sym_rms_scale;
	sym_rms->I0=sqrtf(temp*temp+sym.I0*sym.I0);
	sym_back->I0=sym.I0;

	temp = (sym.I1-sym_back->I1)*sym_rms_scale;
	sym_rms->I1=sqrtf(temp*temp+sym.I1*sym.I1);
	sym_back->I1=sym.I1;

	temp = (sym.I2-sym_back->I2)*sym_rms_scale;
	sym_rms->I2=sqrtf(temp*temp+sym.I2*sym.I2);
	sym_back->I2=sym.I2;


}


// function-6
//spectral analysis
// caution: input structure has to be initialized and used only once
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

	out_scale=sqrt2/(float)qBufferLength;


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








