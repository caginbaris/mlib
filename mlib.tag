SQLite format 3   @    b                                                             b -�   � zA�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  C_indexFilenameBrowseCREATE INDEX Filename ON Browse(Filename)4KindexTagBrowseCREATE INDEX Tag ON Browse(Tag)7OindexNameBrowseCREATE INDEX Name ON Browse(Name)��atableBrowseBrowseCREATE TABLE Browse (Kind INTEGER,Name TEXT,Tag TEXT,Filename TEXT,Lineno INTEGER,Text TEXT,Extra INTEGER)   �    ������������                     =;-ou@�B;=rms_summeasurement_functions.cfloat rms = 0, rms_sum =/�A;#rmsmeasurement_functions.cfloat rms =1�@;+imeasurement_functions.cunsigned int i;A�6;1foutMagspectrameasurement_functions.h
float foutMag[13];C�7;3foutImagspectrameasurement_functions.h	float foutImag[13];C�8;3foutRealspectrameasurement_functions.hfloat foutReal[13];B�9;3qBufferspectrameasurement_functions.hfloat qBuffer[100];(�:;spectrameasurement_functions.h   		';�Acs_generationme��?		;�ktrue_rmsmeasurement_functions.cfloat true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength) {6�>;)rtInputmeasurement_functions.cfloat rtInput,E�=);9delayLineArraymeasurement_functions.cfloat *delayLineArray,O�<-;IdelayLineCountermeasurement_functions.cunsigned int delayLineCounter,   �x   �f   �S   �C   �5   �(   �   �O   �k   �[   �   �~   �    ���������~vnf^WMD=,%��������zsg[ �NA0����������~qdWLA6* � � � � � � � �
p_out�zValues�z2_ptr�z1_ptr�tsktrue_rms1thermal_parametersotaul'sym_rms_scaley
sym_rwsym_out`sym_i3v
sym_ix
sqrtf�
sqrt2sspectraertInput�rtInput~rms_sum�rms_data�rms_V2W
rms_VXrms_I2U
rms_IVrms�rmsnqBufferdpvp_datai%pn_rms_scalezpiq!phase_dispN%phase_cs_outY#phase_cs_inMphase_VPphase_IOoutput�isqrt2ti�i�freezejfoutRealcfoutMagafoutImagb-delayLineCounter|)delayLineArray}'cs_generation�coeff_ptr�#coeffLength�
coeff�cfbg#arrayLength{ah_i2u	_2pirXRVsKVcLV2^V1_V0]RQQSPT-MLIB_CONSTANTS_Hp;MEASUREMENT_FUNCTIONS_HHIsI	Inom   
sqrtf�   out_scale�   �    ��BHNTZ`flrx~�������������$  � � � � � � � � � � � � � ����(:L^p�������� -:GTan{��������������������������HMY`eiopqrstuvwxyz�����~}|{#phase_cs_inI#phase_cs_inJ#phase_cs_inK#phase_cs_inL%phase_cs_outN%phase_cs_outO%phase_cs_outP%phase_cs_outQ%phase_cs_outR%phase_cs_outS%phase_cs_outT%phase_cs_outU%phase_cs_outV%phase_cs_outW%phase_cs_outXpvp_datafpvp_datagpvp_datahspectraaspectrabspectracspectradsym_outZsym_out[sym_out\sym_out]sym_out^sym_out_1thermal_parametersj1thermal_parametersk1thermal_parametersl1thermal_parameter   spectraa   �   `    ���|`+He�����0Mj�����5Ro�����           ;measure;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measurement_functions.h�;measuremen   1mlib_definitions.h�   ;measurement_functions.hn   
;measurement_functions.hQ   	;measurement_functions.c�   ;measurement_functions.c�    Y��Ap�����������%;5twBufferImag   @;=p_inmeasurement_functions.cMstruct phase_cs_in p_in,   );sqrtfmeasurement_   H;Atemp_imagmeasurement   9;)p_peakmeasurement_functions.cfloat p_peak=0   9;)n_peakmeasurement_functions.cfloat n_peak=0   2;!memmeasurement_functions.c$float mem)   J;Mthermmeasurement_functions.c$struct thermal_parameters therm,   4;#tempmeasurement_functions.c'float temp;   @!;/t_constantmeasurement_functions.c(float t_constant;8�H;;
MEASUREMENT_FUNCTIONS_Hmeasurement_functions.h7�I#;Isphase_cs_inmeasurement_functions.hfloat Is;7�J#;Icphase_cs_inmeasurement_functions.hfloat Ic;7�K#;Vsphase_cs_inmeasurement_functions.hfloat Vs;7�L#;Vcphase_cs_inmeasurement_functions.hfloat Vc;,�M#;phase_cs_inmeasurement_functions.h	H�N!%;/phase_dispphase_cs_outmeasurement_functions.hfloat phase_disp;B�O%;)phase_Iphase_cs_outmeasurement_functions.hfloat phase_I;
     `  ` } � � � �+He�����0Mj�����5Ro����ionmeasurA�8%C�k;;measurement_functions.c{;measurement_functions.c|;measurement_functions.c};measurement_functions.c~;measurement_functions.c;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�   3 3x��#\��#d��������7�:;)rtInputmeasurement_functions.cfloat rtIn'�;sqrtfmeasurement_functions.c �'� ;sqrtfmeasurement_functions.c �'�!;sH�<%;ApDataCountermeasurement_functions.cunsigned int pDataCounter,4�=;'pDatameasurem=�G!;/B�P%;)phase_Vphase_cs_outmeasurement_functions.hfloat phase_V;6�Q%;Rphase_cs_outmeasurement_functions.hfloat R;6�R%;Xphase_cs_outmeasurement_functions.hfloat X;6�S%;Qphase_cs_outmeasurement_functions.hfloat Q;6�T%;Pphase_cs_outmeasurement_functions.hfloat P;@�U%;'rms_I2phase_cs_outmeasurement_functions.hfloat rms_I2;>�V%;%rms_Iphase_cs_outmeasurement_functions.hfloat rms_I;@�W%;'rms_V2phase_cs_outmeasurement_functions.hfloat rms_V2;>�X%;%rms_Vphase_cs_outmeasurement_functions.hfloat rms_V;-�Y%;phase_cs_outmeasurement_functions.h3�Z;I0sym_outmeasurement_functions.h-float I0;3�[;I2sym_outmeasurement_functions.h,float I2;
    �  � � �+He�����0Mj�����5Ro����������C�l;;measur;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�
    �  � � �+He�����0Mj�����5Ro����������surement;measur;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.c�;measurement_functions.hH;measurement_functions.hI;measurement_functions.hJ;measurement_functions.hK;measurement_functions.hL;measurement_functions.hM;measurement_functions.hN;measurement_functions.hO;measurement_functions.hP
   K � v`~P�X��� ������h
oH'8/@7?GOW^eox ����������	-CQ+^l �x�������������outR�I0ZI1\I2[IcJIc�	InommIs	Inom�freeze�a�b�c�foutReal�foutImag�foutMag�V1�V2�V0�I1�I2�I0�P�Q�I0ZI1\I2[IcJIc�	InommIsIIs�;MEASUREMENT_FUNCTIONS_HH-MLIB_CONSTANTS_H�1MLIB_DEFINITIONS_H�PTQSRQR�V0]V1_V2^VcLVc�VsKVs�XRX�	_2pi�_i2�ah#arrayLength{atan2f�atan2f�bgcf
coeff�#coeffLength�coeff_ptr�+cs_computations�'cs_generation�!dataLength�)delayLineArray}-delayLineCounter|foutImagbfoutMagafoutRealcfreezejh�h�	hsum�i�i�i�i�i�in�in_back�isqrt2�mem�n_peak�out�out_back�
   ; x�����������
$1BS8eu�����+����a�lJU+8ER_l��������������������%phase_cs_outY!phase_dispN!phase_disp�pi�%pn_rms_scale�pvp_datai!pvp_filter�qBufferrms�pvp_data�spectra�qBuffer�%phase_cs_out�
rms_V�rms_V2�
rms_I�rms_I2�output�pCounter�
pData�%pDataCounter�	p_in�
p_out�p_peak�pa�pb�pc�+peak_detect_rms�phase_IOphase_I�phase_VPphase_V�#phase_cs_inM#phase_cs_in�%phase_cs_outY!phase_dispN!phase_disp�pi�%pn_rms_scale�pvp_datai!pvp_filter�qBufferd'qBufferLength�rmsnrms�
rms_IVrms_I2U
rms_VXrms_V2Wrms_data�rms_sum�rtInput~rtInput�rtInput�rtInput�)signal_spectra�!signal_thd�spectrae
sqrt2�
sqrtf�
sqrtf�
sqrtf�
sqrtf�
sqrtf�
sqrtf�
sqrtf�
sqrtf�
   `� ����������� &,28>DJPV\bhntz����������������������
"(.4:@FLRX^djpv|���������������������6<BHNTZ`flrx��������#phHMY`eio{|}~���������������������HMY`eio{|}~������������������������������������������������������������������������������������
   /V \bhntz����V������� 2DVhz�������
.@Rdv��������Xp����� {�1thermal%phase_cs_outO%phase_cs_outP%phase_cs_outQ%phase_cs_outR%phase_cs_outS%phase_cs_outT%phase_cs_outU%phase_cs_outV%phase_cs_outW%phase_cs_outX%phase_cs_out�%phase_ ����������#phase_cs_inI#phase_cs_inJ#phase_cs_inK#phase_cs_inL#phase_cs_in�#phase_cs_in�#phase_cs_in�#phase_cs_in�%phase_cs_outN%phase_cs_outO%phase_cs_outP%phase_cs_outQ%phase_cs_outR%phase_cs_outS%phase_cs_outT%phase_cs_outU%phase_cs_outV%phase_cs_outW%phase_cs_outX%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�%phase_cs_out�pvp_datafpvp_datagpvp_datahpvp_data�pvp_data�pvp_data�    Z �v)����v. � � � Z�������� -:GTan{����� � � � � �4�5;/hmeasurement_functions.c �struct spectra h)'�4;sqrtfmeasurement_functions.c �2�3;+imeasurement_functions.c �unsigned int i;;�2;-out_scalemeasurement_functions.c �float out_scale;E�1;Atemp_imagmeasurement_functions.c �float temp_real,temp_imag;;�0;-temp_realmeasurement_functions.c �float temp_real,7�/;)x_errormeasurement_functions.c �float x_error;�E�.		);�7signal_spectrameasurement_functions.c �void signal_spectra( float rtInput, struct spectra *h, unsigned int qBufferLength, float *twBufferReal, float *twBufferImag, unsigned int pCounter) {7�-;)rtInputmeasurement_functions.c �float rtInput,5�,;1hmeasurement_functions.c �struct spectra *h,J�+';CqBufferLengthmeasurement_functions.c �unsigned int qBufferLength,B�*%;5twBufferRealmeasurement_functions.c �float *twBufferReal,B�)%;5twBufferImagmeasurement_functions.c �float *twBufferImag,
    �  � �+He�����0Mj�����5Ro����|�����tameasur;measurement_functions.hU;measurement_functions.hV;measurement_functions.hW;measurement_functions.hX;measu;measurement_functions.hR;measurement_functions.hS;measurement_functions.hT;measurement_functions.hU;measurement_functions.hV;measurement_functions.hW;measurement_functions.hX;measurement_functions.hY;measurement_functions.hZ;measurement_functions.h[;measurement_functions.h\;measurement_functions.h];measurement_functions.h^;measurement_functions.h_;measurement_functions.h`;measurement_functions.ha;measurement_functions.hb;measurement_functions.hc;measurement_functions.hd;measurement_functions.he;measurement_functions.hf;measurement_functions.hg;measurement_functions.hh;measurement_functions.hi;measurement_functions.hj;measurement_functions.hk;measurement_functions.hl;measurement_functions.hm    W d0��i@��O � � W ? ? ,   ''�I;
sym_i3m@�";9pCountermeas5�;)z1_ptrmeasurement_functions.c2float *z1_ptr,1�
;+imeasurement_functions.c1unsigned int i;��			';�Acs_generationmeasurement_functions.c/float cs_generation(float rtInput,float *coeff, unsigned int coeffLength, float *zValues){6�;)rtInputmeasurement_functions.c/float rtInput,3�;'coeffmeasurement_functions.c/float *coeff,E�#;?coeffLengthmeasurement_functions.c/unsigned int coeffLength,7�;+zValuesmeasurement_functions.c/float *zValues)&�;sqrtfmeasurement_functions.c"O�;Yrms_datameasurement_functions.cfloat rms = 0, rms_sum = 0, rms_data=0@�;=rms_summeasurement_functions.cfloat rms = 0, rms_sum =/�;#rmsmeasurement_functions.cfloat rms =1� ;+imeasurement_functions.cunsigned int i;��		;�ktrue_rmsmeasurement_functions.cfloat true_rms(float rtInput, float *delayLineArray, unsigned int delayLineCounter, unsigned int arrayLength) {    ? �zE��j-�� � v ?�����     7�K#;4�C;'pDatameasurement_functions.cfloat *pData,H�B%;ApDataCountermeasurement_functions.cunsigned int pDataCounter,D�A!;=dataLengthmeasurement_functions.cunsigned int dataLength)�"�@		!;�ypvp_filtermeasurement_functions.c �void pvp_filter(struct pvp_data in,struct pvp_data *in_back,struct pvp_data *out, struct pvp_data *out_back,float ts){7�?;3inmeasurement_functions.c �struct pvp_data in,B�>;?in_backmeasurement_functions.c �struct pvp_data *in_back,:�=;7outmeasurement_functions.c �struct pvp_data *out,D�<;Aout_backmeasurement_functions.c �struct pvp_data *out_back,-�;;tsmeasurement_functions.c �float ts)'�:;sqrtfmeasurement_functions.c �7�9;1thdmeasurement_functions.c �float hsum=0,thd=02�8;%hsummeasurement_functions.c �float hsum=02�7;+imeasurement_functions.c �unsigned int i;N�6	!;Ssignal_thdmeasurement_functions.c �float signal_thd(struct spectra h){    ? n7��n3�~T*  � � � ?��.�tatusmeasurem@�(;9pCountermeasurement_functions.c �unsigned int pCounter)'�';sqrtfmeasurement_functions.c �'�&;sqrtfmeasurement_functions.c �'�%;sqrtfmeasurement_functions.c �'�$;sqrtfmeasurement_functions.c �'�#;sqrtfmeasurement_functions.c �'�";sqrtfmeasurement_functions.c �1�!;#tempmeasurement_functions.c �float temp;~� 		;�7sym_magmeasurement_functions.c �void sym_mag(struct sym_out sym, struct sym_out *sym_back, struct sym_out *sym_rms ){8�;3symmeasurement_functions.c �struct sym_out sym,C�;?sym_backmeasurement_functions.c �struct sym_out *sym_back,B�;?sym_rmsmeasurement_functions.c �struct sym_out *sym_rms );�;5temp_imeasurement_functions.cufloat temp_r,temp_i;4�;'temp_rmeasurement_functions.cufloat temp_r,��		;�Wsym_compmeasurement_functions.csvoid sym_comp(struct phase_cs_in pa, struct phase_cs_in pb,struct phase_cs_in pc,struct sym_out*sym){    a  a � �9d��4y��Co�                                                         3�\;I1sym_outmeasurement_functions.h+float I1;3�];V0sym_outmeasurement_functions.h)float V0;3�^;V2sym_outmeasurement_functions.h(float V2;3�_;V1sym_outmeasurement_functions.h'float V1;(�`;sym_outmeasurement_functions.h$A�a;1foutMagspectrameasurement_functions.h:float foutMag[13];C�b;3foutImagspectrameasurement_functions.h9float foutImag[13];C�c;3foutRealspectrameasurement_functions.h8float foutReal[13];B�d;3qBufferspectrameasurement_functions.h7float qBuffer[100];(�e;spectrameasurement_functions.h52�f;cpvp_datameasurement_functions.hOfloat c;2�g;bpvp_datameasurement_functions.hNfloat b;2�h;apvp_datameasurement_functions.hMfloat a;)�i;pvp_datameasurement_functions.hKM�j1;5freezethermal_parametersmeasurement_functions.h]unsigned int freeze;>�k1;tsthermal_parametersmeasurement_functions.h\float ts;s  X �x5�+ � � X X X X X X X X X X X X                                          6�~;)rtInputmeasurement_functions.cfloat rtInput,E�});9delayLineArraymeasurement_functions.cfloat *delayLineArray,O�|-;IdelayLineCountermeasurement_functions.cunsigned int delayLineCounter,E�{#;?arrayLengthmeasurement_functions.cunsigned int arrayLength)  �%-
pn_rms_scalemlib_constants.h  c'-
sym_rms_scalemlib_constants.h  9-
sym_imlib_constants.h  -
sym_rmlib_constants.h   �-
sym_i3mlib_constants.h   �-
_i2mlib_constants.h   �-
isqrt2mlib_constants.h   �-
sqrt2mlib_constants.h
   m-
_2pimlib_constants.h	   L-
pimlib_constants.h   ---
MLIB_CONSTANTS_Hmlib_constants.h3�o1;thermal_parametersmeasurement_functions.hW@�n1;!rmsthermal_parametersmeasurement_functions.hYfloat rms;B�m1;#Inomthermal_parametersmeasurement_functions.hZfloat Inom;@�l1;!tauthermal_parametersmeasurement_functions.h[float tau;    � �r;��8���X � �s K K                                   (		';�Acs_generationmeasurement_functions.c/float cs_generation(float rtInput,float9�;9pameasurement_functions.csstruct phase_cs_in pa,9�;9pbmeasurement_functions.csstruct phase_cs_in pb,9�;9pcmeasurement_functions.csstruct phase_cs_in pc,7�;3symmeasurement_functions.csstruct sym_out*sym)'�;atan2fmeasurement_functions.cb'�;atan2fmeasurement_functions.ca&�;sqrtfmeasurement_functions.cS&�;sqrtfmeasurement_functions.cR{�		+;�#cs_computationsmeasurement_functions.cLvoid cs_computations(struct phase_cs_in p_in, struct phase_cs_out *p_out ){=�;=p_inmeasurement_functions.cLstruct phase_cs_in p_in,B�;Ep_outmeasurement_functions.cLstruct phase_cs_out *p_out )4�;'outputmeasurement_functions.c3float output;K�;Ocoeff_ptrmeasurement_functions.c2float *z1_ptr,*z2_ptr,*coeff_ptr;=�;9z2_ptrmeasurement_functions.c2float *z1_ptr,*z2_ptr,    C ���oF��vA��Q � � � C                ;�f%1'rms_I2phase_cs_outmlib_definitions.hfloat rms_I2;1�e%1Pphase_cs_outmlib_definitions.hfloat P;1�d%1Qphase_cs_outmlib_definitions.hfloat Q;1�c%1Xphase_cs_outmlib_definitions.hfloat X;1�b%1Rphase_cs_outmlib_definitions.hfloat R;=�a%1)phase_Vphase_cs_outmlib_definitions.hfloat phase_V;=�`%1)phase_Iphase_cs_outmlib_definitions.hfloat phase_I;C�_!%1/phase_dispphase_cs_outmlib_definitions.hfloat phase_disp;'�^#1phase_cs_inmlib_definitions.h2�]#1Vcphase_cs_inmlib_definitions.h	float Vc;2�\#1Vsphase_cs_inmlib_definitions.h
float Vs;2�[#1Icphase_cs_inmlib_definitions.hfloat Ic;2�Z#1Isphase_cs_inmlib_definitions.hfloat Is;.�Y11
MLIB_DEFINITIONS_Hmlib_definitions.h&�X%-
pn_rms_scalemlib_constants.h'�W'-
sym_rms_scalemlib_constants.h�V-
sym_imlib_constants.h�U-
sym_rmlib_constants.h �T-
sym_i3mlib_constants.h
   & IT_hq����<�����+"1=LU`x���4������_outmeasurement_functions.h!float rms_I2;6�Y%;Pphase_cs_outmeasurement_fu9�;9pameasurement_functions.ctstruct phase_cs_in pa,9�;9pbmeasurement_functions.ctstruct phase_cs_in pb,9�;9pcmeasurement_functions.ctstruct phase_cs_in pc,7�;3symmeasurement_functions.ctstruct sym_out*sym)'�;atan2fmeasurement_functions.cc'�;atan2fmeasurement_functions.cb&�;sqrtfmeasurement_functions.cT&�1thermal_parameters tau�ts�sym_out�
sqrtf�
sqrtf�sym�sym�sym_back�sym_comp�
sym_i�sym_i3�sym_mag�sym_out`
sym_r�sym_rms�'sym_rms_scale�!t_constant�taul	temp�	temp�temp_i�temp_imag�temp_r�temp_real�thd�
therm�1thermal_parameterso)thermal_status�true_rmstskts�%twBufferImag�%twBufferReal�x_error�z1_ptr�z2_ptr�zValues�
   Y Yv������&<Rh������(@Xp�����(@Xp�����ment_fun-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definiti;measurement_functions.ho-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�-mlib_constants.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�    L �+���Q�^ � � � � l L :                @�V4�=;'p�S-
_i2mlib_constants.h �R-
isqrt2mlib_constants.h�Q-
sqrt2mlib_constants.h
�P-
_2pimlib_constants.h	�O-
pimlib_constants.h*�N--
MLIB_CONSTANTS_Hmlib_constants.h=�M!;/t_constantmeasurement_functions.c'float t_constant;1�L;#tempmeasurement_functions.c&float temp;r�K		);�thermal_statusmeasurement_functions.c#float thermal_status(struct thermal_parameters therm, float mem) {G�J;Mthermmeasurement_functions.c#struct thermal_parameters therm,/�I;!memmeasurement_functions.c#float mem)6�H;)n_peakmeasurement_functions.cfloat n_peak=06�G;)p_peakmeasurement_functions.cfloat p_peak=03�F;-imeasurement_functions.cunsigned int i=0��E		+;�Ypeak_detect_rmsmeasurement_functions.cfloat peak_detect_rms(float rtInput, float *pData,unsigned int pDataCounter, unsigned int dataLength){7�D;)rtInputmeasurement_functions.cfloat rtInput,    L ��J���[*���S � � | L                                -�x1bpvp_datamlib_definitions.hDfloat b;-�w1cpvp_datamlib_definitions.hEfloat c;#�v1spectramlib_definitions.h-=�u13qBufferspectramlib_definitions.h/float qBuffer[100];>�t13foutRealspectramlib_definitions.h0float foutReal[13];>�s13foutImagspectramlib_definitions.h1float foutImag[13];<�r11foutMagspectramlib_definitions.h2float foutMag[13];#�q1sym_outmlib_definitions.h.�p1V1sym_outmlib_definitions.h"float V1;.�o1V2sym_outmlib_definitions.h#float V2;.�n1V0sym_outmlib_definitions.h$float V0;.�m1I1sym_outmlib_definitions.h&float I1;.�l1I2sym_outmlib_definitions.h'float I2;.�k1I0sym_outmlib_definitions.h(float I0;(�j%1phase_cs_outmlib_definitions.h9�i%1%rms_Vphase_cs_outmlib_definitions.hfloat rms_V;;�h%1'rms_V2phase_cs_outmlib_definitions.hfloat rms_V2;9�g%1%rms_Iphase_cs_outmlib_definitions.hfloat rms_I;   5 ��^"��f5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             .� 11thermal_parametersmlib_definitions.hI;�11!rmsthermal_parametersmlib_definitions.hKfloat rms;=�~11#Inomthermal_parametersmlib_definitions.hLfloat Inom;;�}11!tauthermal_parametersmlib_definitions.hMfloat tau;9�|11tsthermal_parametersmlib_definitions.hNfloat ts;H�{115freezethermal_parametersmlib_definitions.hOunsigned int freeze;$�z1pvp_datamlib_definitions.hA-�y1apvp_datamlib_definitions.hCfloat a;
   � 8Ph������(@Xp����� �                                                                                                                                                                                                                                                                                                                                                                                                                                                            1mlib_definitions.h 1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�1mlib_definitions.h�
    IVcp}����������&3@Xp�����1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       1thermal_parameters�1thermal_parameters�spectrabspectracspectradspectra�spectra�spectra�spectra�sym_outZsym_out[sym_out\sym_out]sym_out^sym_out_sym_out�sym_out�sym_out�sym_out�sym_out�sym_out�1thermal_parametersj1thermal_parametersk1thermal_parametersl1thermal_parametersm1thermal_parametersn1thermal_parameters�1thermal_parameters�1thermal_parameters�