SQLite format 3   @    6                                                           6 -�   � zA�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  C_indexFilenameBrowseCREATE INDEX Filename ON Browse(Filename)4KindexTagBrowseCREATE INDEX Tag ON Browse(Tag)7OindexNameBrowseCREATE INDEX Name ON Browse(Name)��atableBrowseBrowseCREATE TABLE Browse (Kind INTEGER,Name TEXT,Tag TEXT,Filename TEXT,Lineno INTEGER,Text TEXT,Extra INTEGER)   �    ��������� � D�                     =;-ou@�B;=rms_summeasurement_functions.cfloat rms = 0, rms_sum =/�A;#rmsmeasurement_functions.cfloat rms =1�@;+imeasurement_functions.cunsigned int i;A�6;1foutMagspectrameasurement_functions.h
float foutMag[13];C�7;3foutImagspectrameasurement_functions.h	float foutImag[13];C�8;3foutRealspectrameasurement_functions.hfloat foutReal[13];B�9;3qBufferspectrameasurement_functions.hfloat qBuffer[100];(�:;spectrameasurement_functions.h   		';�Acs_generationme��?		;�ktrue_rmsmeasurement_functions.cfloat true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength) {6�>;)rtInputmeasurement_functions.cfloat rtInput,E�=);9delayLineArraymeasurement_functions.cfloat *delayLineArray,O�<-;IdelayLineCountermeasurement_functions.cunsigned int delayLineCounter,E�;#;?arrayL� �}   �:   �.   �   �H   �   �   �i   �Y   �    ��������������� � ����[! �H�(6CQ]X 
 � � ��������!d.9EPo{�}�d � � ������m<0�papbpcsymatan2fatan2f
sqrtf
sqrtf+cs_computations	p_in
p_outoutputcoeff_ptrz2_ptrz1_ptr'cs_generation
#coeffLengthI0�I1�I2�Ic�Is�P�Q�R�V0�V1�V2�Vc�Vs�X�	_2pi�_i2�a�b�c�-delayLineCounter�#arrayLength�foutImag�foutMag�foutReal�iisqrt2�rmstrue_rms rtInput�)delayLineArray�phase_I�phase_V�#phase_cs_in�%phase_cs_out�!phase_disp�pi�%pn_rms_scale�pvp_data�\ rms_sumqBuffer�
rms_I�rms_I2�
rms_V�rms_V2� i
sqrtfrms_dataspectra�
sqrt2�zValues
sym_i�sym_i3��
coeffsym_out�
sym_r�'sym_rms_scale�    pc   �    �������������������|vpjd^XRLF@:4.("
����������������������ztnhb\VPJD>82,& �������|jXF4" � � � �   sym_outasym_out`sym_out_%phase_cs_out]%phase_cs_out\%phase_cs_out[%phase_cs_outZ%phase_cs_outY%phase_cs_outX%phase_cs_outW%phase_cs_outV%phase_cs_outU%phase_cs_outT%phase_cs_outS#phase_cs_inQ#phase_cs_inP#phase_cs_inO#phase_cs_inN^RMLKJIHGFEDCBA@?>=<;:9876543210/.-,+*)('&%$#"! 
	 ��   <   �    ��� } `+He�����0Mj�����5Ro�����           ;measure;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_f   
;measurement_functions.hN   	;measurement_functions.c5   ;measurement_functions.c   
% %m�@�B����������%;5twBufferImag    >;5temp_imeasurement_functio   ';
_i2measurement_functions.h   *;
sym_i3measurement_functions.h   );
sy   5;apvp_datameasurement_functions.hUfloat a;   ,;pvp_datameasurement_functions.hSE�|#;?arrayLengthmeasurement_functions.cunsigned int arrayLength)O�}-;IdelayLineCountermeasurement_functions.cunsigned int delayLineCounter,E�~);9delayLineArraymeasurement_functions.cfloat *delayLineArray,6�;)rtInputmeasurement_functions.cfloat rtInput,�� 		;�ktrue_rmsmeasurement_functions.cfloat true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength) {1�;+imeasurement_functions.cunsigned int i;/�;#rmsmeasurement_functions.cfloat rms =@�;=rms_summeasurement_functions.cfloat rms = 0, rms_sum =O�;Yrms_datameasurement_functions.cfloat rms = 0, rms_sum = 0, rms_data=0&�;sqrtfmeasurement_functions.c"
    �  � �+He�����0Mj�����5Ro��������ionmeasurA�8%C�k;;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c ;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c	;measurement_functions.c
;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c    G������6�������7�h#;Isphase_cs_inmeasurement_functions.hfloat Is;7�i#;7�=;)x_errormeasurement_functions.6�f%;Pphase_cs_outmeasurement_functions.h"float P;@�g%;'rms_I2phase_cs_out(�w;spectrameasurement_f7�;+zValuesmeasurement_functions.c/float *zValues)E�#;?coeffLengthmeasurement_functions.c/unsigned int coeffLength,3�;'coeffmeasurement_functions.c/float *coeff,6�	;)rtInputmeasurement_functions.c/float rtInput,��
		';�Acs_generationmeasurement_functions.c/float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues){1�;+imeasurement_functions.c1unsigned int i;5�;)z1_ptrmeasurement_functions.c2float *z1_ptr,=�;9z2_ptrmeasurement_functions.c2float *z1_ptr,*z2_ptr,K�;Ocoeff_ptrmeasurement_functions.c2float *z1_ptr,*z2_ptr,*coeff_ptr;4�;'outputmeasurement_functions.c3float output;B�;Ep_outmeasurement_functions.cLstruct phase_cs_out *p_out )
    �  � �+He�����0Mj�����5Ro��������_scaleC�l;3foutIm;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c ;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c;measurement_functions.c ;measurement_functions.c!;measurement_functions.c";measurement_functions.c#;measurement_functions.c$;measurement_functions.c%;measurement_functions.c&;measurement_functions.c';measurement_functions.c(;measurement_functions.c);measurement_functions.c*;measurement_functions.c+;measurement_functions.c,;measurement_functions.c-;measurement_functions.c.;measurement_functions.c/;measurement_functions.c0;measurement_functions.c1;measurement_functions.c2;measurement_functions.c3;measurement_functions.c4
   H He�����0Mj�����5Ro����5Ro����imeasurement_functions;measurement_functions.c>;measurement_functions.c?;measurement_functions.c@;measurement_functions.cA;measurement_functions.cB;measurement_functions.cC;measurement_functions.hD;measurement_functions.hE;measurement_functions.c6;measurement_functions.c7;measurement_functions.c8;measurement_functions.c9;measurement_functions.c:;measurement_functions.c;;measurement_functions.c<;measurement_functions.c=;measurement_functions.c>;measurement_functions.c?;measurement_functions.c@;measurement_functions.cA;measurement_functions.cB;measurement_functions.cC;measurement_functions.hD;measurement_functions.hE;measurement_functions.hF;measurement_functions.hG;measurement_functions.hH;measurement_functions.hI;measurement_functions.hJ;measurement_functions.hK;measurement_functions.hL;measurement_functions.hM
 6�uem��}��]U�����������+:O�bv���=��'����.�D�������Po{�}�d � � ������m<0�papbpcsymatan2fatan2f
sqrtf
sqrtf+cs_computations	p_in
p_outoutputcoeff_ptrz2_ptrz1_ptr'cs_generation
#coeffLengthI0�I1�I2�Ic�Is�P�Qambln_peakCp_peakBiA
pData>%pDataCounter=!dataLength<in:in_back9out8out_back7i4out_scale3h-pCounter) V2cV0bI1aI2`I0_PYQXXWRVVcQVsPIcOIsN_i2Hisqrt2G	_2piE#arrayLength�atan2fatan2f� bV1d
coeff#coeffLengthcoeff_ptr+cs_computations'cs_generation
)delayLineArray�-delayLineCounter�foutRealhfoutImaggfoutMagfii   ickoutput	p_in
p_outpapb
M < � ���d�T�u ����H��<!.t` ��S^i������*t}����(��55AD�"��S����d � � ������m<0�pa4pb3+peak_detect_rms!pvp_filterpvp_datanspectraj+peak_detect_rms@rtInput?!pvp_filter;ts6
sqrtf5temp_imag2temp_real1x_error0)signal_spectra/rtInput.'qBufferLength,%twBufferReal+%twBufferImag*
sqrtf(
sqrtf'
sqrtf&
sqrtf%
sqrtf$
sqrtf#	temp"sym_mag!sym sym_backsym_rmstemp_itemp_r� pha!phase_dispS#phase_cs_inR%pn_rms_scaleM'sym_rms_scaleL
sym_iK
sym_rJsym_i3I
sqrt2FpiDrms� rms_rms_I2Zphase_VUphase_ITrms_datarms_sumrtInput�rtInput	rms_V2\
rms_I[
sqrtf
sqrtf
sqrtfsymsym_comp   sym_i�qBufferisym_oute%phase_cs_out^
rms_V]true_rms z1_ptrz2_ptrzValues
     ` 0Mj�����5Ro���������eH+ � � � � } `y����������;measurement_functions.hn;measurement_functions.hm;measurement_functions.hl;measurement_functions.hk;measurement_functions.hj;measurement_functions.hi;measurement_functions.hh;measurement_functions.hg;measurement_functions.hf;measurement_functions.he;measurement_functions.hd;measurement_functions.hc;measurement_functions.hb;measurement_functions.ha;measurement_functions.h`;measurement_functions.h_;measurement_functions.hO;measurement_functions.hP;measurement_functions.hQ;measurement_functions.hR;measurement_functions.hS;measurement_functions.hT;measurement_functions.hU;measurement_functions.hV;measurement_functions.hW;measurement_functions.hX;measurement_functions.hY;measurement_functions.hZ;measurement_functions.h[;measurement_functions.h\;measurement_functions.h];measurement_functions.h^    B ����vL"���F � | B������eH+ � � �7�.;)rtInputmeasurement_functions.c �float rtInput,5�-;1hmeasurement_functions.c �struct spectra *h,J�,';CqBufferLengthmeasurement_functions.c �unsigned int qBufferLength,B�+%;5twBufferRealmeasurement_functions.c �float *twBufferReal,B�*%;5twBufferImagmeasurement_functions.c �float *twBufferImag,@�);9pCountermeasurement_functions.c �unsigned int pCounter)'�(;sqrtfmeasurement_functions.c �'�';sqrtfmeasurement_functions.c �'�&;sqrtfmeasurement_functions.c �'�%;sqrtfmeasurement_functions.c �'�$;sqrtfmeasurement_functions.c �'�#;sqrtfmeasurement_functions.c �1�";#tempmeasurement_functions.c �float temp;~�!		;�7sym_magmeasurement_functions.c �void sym_mag(struct sym_out sym, struct sym_out *sym_back, struct sym_out *sym_rms ){8� ;3symmeasurement_functions.c �struct sym_out sym,C�;?sym_backmeasurement_functions.c �struct sym_out *sym_back,
   @� ����������������������
"(.4:@FLRX^djpv|���������������������ztnhb\VPJD>82,& �������|jXF4" � � � � trsym_outasym_out`sym_out_%phase_cs_out]%phase_cs_out\%phase_cs_out[%phase_cs_outZ%phase_cs_outY%phase_cs_outX%phase_cs_outW%phase_cs_outV%phase_cs_outU%phase_cs_outT%phase_cs_outS#phase_cs_inQ#phase_cs_inP#phase_cs_inO#phase_cs_inN^RMLKJIHGFEDCBA@?>=<���� 	
 !"#$%&'()*+,-./0123456789:;    D �|9���\&���Y � � D         �����;B�i;3qBufferspectrameasurement_functions.hBfloat qBuffer[100];C�h;3foutRealspectrameasurement_functions.hCfloat foutReal[13];C�g;3foutImagspectrameasurement_functions.hDfloat foutImag[13];A�f;1foutMagspectrameasurement_functions.hEfloat foutMag[13];(�e;sym_outmeasurement_functions.h/3�d;V1sym_outmeasurement_functions.h2float V1;3�c;V2sym_outmeasurement_functions.h3float V2;3�b;V0sym_outmeasurement_functions.h4float V0;3�a;I1sym_outmeasurement_functions.h6float I1;3�`;I2sym_outmeasurement_functions.h7float I2;3�_;I0sym_outmeasurement_functions.h8float I0;-�^%;phase_cs_outmeasurement_functions.h>�]%;%rms_Vphase_cs_outmeasurement_functions.hfloat rms_V;@�\%;'rms_V2phase_cs_outmeasurement_functions.hfloat rms_V2;>�[%;%rms_Iphase_cs_outmeasurement_functions.h float rms_I;@�Z%;'rms_I2phase_cs_outmeasurement_functions.h!float rms_I2;
   2� PV\bhntz�����������0������*<N`r�������#	�����C6�����C60BTbp~���������JD>82,
�������|jXF4" � � � �s_outmsym_outsym_out
sym_out	sym_outsym_outsym_out%phase_cs_out%phase_cs_o�����   measurement_functions.h�   measurement_functions.h�   measurement_functions.h�   measurnpvp_datampvp_datalpvp_datakjspectraispectrahspectragspectrafesym_outdsym_outc=>?@ABCDEFGHIJKLMR^#phase_cs_inN#phase_cs_inO#phase_cs_inP#phase_cs_inQ%phase_cs_outS%phase_cs_outT%phase_cs_outU%phase_cs_outV%phase_cs_outW%phase_cs_outX%phase_cs_outY%phase_cs_outZ%phase_cs_out[%phase_cs_out\%phase_cs_out]sym_out_sym_out`sym_outasym_outb    S ���S#��u;�|7 � � � S�������Ypeak_detect_rmsme3�p6�Y%;Pphase_cs_outmeasurement_functions.h"float P;6�X%;Qphase_cs_outmeasurement_functions.h#float Q;6�W%;Xphase_cs_outmeasurement_functions.h$float X;6�V%;Rphase_cs_outmeasurement_functions.h%float R;B�U%;)phase_Vphase_cs_outmeasurement_functions.h&float phase_V;B�T%;)phase_Iphase_cs_outmeasurement_functions.h'float phase_I;H�S!%;/phase_dispphase_cs_outmeasurement_functions.h(float phase_disp;,�R#;phase_cs_inmeasurement_functions.h7�Q#;Vcphase_cs_inmeasurement_functions.hfloat Vc;7�P#;Vsphase_cs_inmeasurement_functions.hfloat Vs;7�O#;Icphase_cs_inmeasurement_functions.hfloat Ic;7�N#;Isphase_cs_inmeasurement_functions.hfloat Is;-�M%;
pn_rms_scalemeasurement_functions.h.�L';
sym_rms_scalemeasurement_functions.h&�K;
sym_imeasurement_functions.h
&�J;
sym_rmeasurement_functions.h	'�I;
sym_i3measurement_functions.h                                         sqrtfme�7�B�;3qBufferspectrame7�$;3symmeasurement_functions.csstruct sym_out*sym)'�#;atan2fmeasurement_functions.cb'�";atan2fmeasurement_functions.ca&�!;sqrtfmeasurement_functions.cS&� ;sqrtfmeasurement_functions.cR{�		+;�#cs_computationsmeasurement_functions.cLvoid cs_computations(struct phase_cs_in p_in, struct phase_cs_out *p_out ){=�;=p_inmeasurement_functions.cLstruct phase_cs_in p_in,B�;Ep_outmeasurement_fuB�;Ep_outmeasurement_functions.cLstruct phase_cs_out *p_out )=�;=p_inmeasurement_functions.cLstruct phase_cs_in p_in,   ~	+;�#cs_computationsmeasurement_functions.cLvoid cs_computations(struct phase_cs_in p_in, struct phase_cs_out *p_out ){&� ;sqrtfmeasurement_functions.cR&�!;sqrtfmeasurement_functions.cS'�";atan2fmeasurement_functions.ca'�#;atan2fmeasurement_functions.cb7�$;3symmeasurement_functions.csstruct sym_out*sym)    � 7��w9��c& � ��y? � � Qo�����������"4FXH�`!%;/phase_dispphase_cs_outmeasurement_functions.h(float phase_disp;,�_#;7�:;3inmeasurement_functions.c �struct pvp_data in,B�9;?in_backmeasurement_functions.c �struct pvp_data *in_back,:�8;7outmeasurement_functions.c �struct pvp_data *out,D�7;Aout_backmeasurement_functions.c �struct pvp_data *out_back,-�6;tsmeasurement_functions.c �float ts)'�5;sqrtfmeasurement_functions.c �2�4;+imeasurement_functions.c �unsigned int i;;�3;-out_scalemeasurement_functions.c �float out_scale;E�2;Atemp_imagmeasurement_functions.c �float temp_real,temp_imag;;�1;-temp_realmeasurement_functions.c �float temp_real,7�0;)x_errormeasurement_functions.c �float x_error;�E�/		);�7signal_spectrameasurement_functions.c �void signal_spectra( float rtInput, struct spectra *h, unsigned int qBufferLength, float *twBufferReal, float *twBufferImag, unsigned int pCounter) {    b �B���b&�� � � b               Vphase_@�)7�H;)rB�*;Ep_outmeB�;?sym_rmsmeasurement_functions.c �struct sym_out *sym_rms );�;5temp_imeasurement_functions.cufloat temp_r,temp_i;4�;'temp_rmeasurement_functions.cufloat temp_r,��		;�Wsym_compmeasurement_functions.csvoid sym_comp(struct phase_cs_in pa, struct phase_cs_in pb,struct phase_cs_in pc,struct sym_out*sym){9�;9pameasurement_functions.csstruct phase_cs_in pa,9�;9pbmeasurement_functions.csstruct phase_cs_in pb,9�;9pcmeasurement_functions.csstruct phase_cs_in pc,7�;3symmeasurement_functions.csstruct sym_out*sym)'�;atan2fmeasurement_functions.cb'�;atan2fmeasurement_functions.ca&�;sqrtfmeasurement_functions.cS&�;sqrtfmeasurement_functions.cR{�		+;�#cs_computationsmeasurement_functions.cLvoid cs_computations(struct phase_cs_in p_in, struct phase_cs_out *p_out ){=�;=p_inmeasurement_functions.cLstruct phase_cs_in p_in,   
 ��k6
���f� � n < B B         6�;/�;#rmsmeasurement_functions.cfloat rms =1�;+imeasurement_functions.cunsigned int i;�� 		;�ktrue_rmsmeasurement_functions.cfloat true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength) {C�u;3foutRealspectrameasurement_functions.hCfloat foutReal[13];B�v;3qBufferspectrameasurement_functions.hBfloat qBuffer[100];(�w;spectrameasureE�|#;?arrayLengthmeasurement_functions.cunsigned int arrayLength)O�}-;IdelayLineCountermeasurement_functions.cunsigned int delayLineCounter,E�~);9delayLineArraymeasurement_functions.cfloat *delayLineArray,6�;)rtInputmeasurement_functions.cfloat rtInput,�� 		;�k)�n;pvp_datameasurement_functions.hS2�m;apvp_datameasurement_functions.hUfloat a;2�l;bpvp_datameasurement_functions.hVfloat b;2�k;cpvp_datameasurement_functions.hWfloat c;(�j;spectrameasurement_functions.h@    L Z��W��M � � � s L 3 4�s� �		+;�Ypeak_detect_rmsmeasur$�H;
_i2measurement_functions.h'�G;
isqrt2measurement_functions.h&�F;
sqrt2measurement_functions.h%�E;
_2pimeasurement_functions.h#�D;
pimeasurement_functions.h6�C;)n_peakmeasurement_functions.c �float n_peak=06�B;)p_peakmeasurement_functions.c �float p_peak=03�A;-imeasurement_functions.c �unsigned int i=0��@		+;�Ypeak_detect_rmsmeasurement_functions.c �float peak_detect_rms(float rtInput, float *pData,unsigned int pDataCounter, unsigned int dataLength){7�?;)rtInputmeasurement_functions.c �float rtInput,4�>;'pDatameasurement_functions.c �float *pData,H�=%;ApDataCountermeasurement_functions.c �unsigned int pDataCounter,D�<!;=dataLengthmeasurement_functions.c �unsigned int dataLength)�"�;		!;�ypvp_filtermeasurement_functions.c �void pvp_filter(struct pvp_data in,struct pvp_data *in_back,struct pvp_data *out, struct pvp_data *out_back,float ts){