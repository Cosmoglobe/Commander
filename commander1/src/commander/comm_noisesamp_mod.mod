  �  C   k820309    l          18.0        �\                                                                                                          
       comm_noisesamp_mod.f90 COMM_NOISESAMP_MOD              BAND MY_CHISQ NU                      @                              
                            @                              
                            @                              
                         @                                '                     #X    #Y    #Z    #W    #GSET 	   #EMPTY 
                � $                                                             � $                                                            � $                                                            � $                                                            � $                             	               
                � $                             
                                    @               @                '8                   #CMB_AMP    #FG_AMP    #TEMP_AMP              �                                                            
            &                   &                                                     �                                          `                 
            &                   &                   &                                                     �                                          �                 
            &                   &                                                             @                               '                                                                                                                                                                                                                                                                                            #         @                                                      #INTEGER    #STRING              
                                                                                                          1 #         @                                                   
   #PARFILE    #PARNAME    #PAR_INT    #PAR_CHAR    #PAR_STRING    #PAR_SP    #PAR_DP    #PAR_LGT    #PAR_PRESENT    #DESC                                                                   1                                                               1                                                                                                                     1                                                               1                                                	                                                      
                                                                                                                                                                                  1 #         @                                  !                    #MAP_ID "   #COEFF #   #OUTPUT_STATS $   #CHISQ_FULLSKY %   #CHISQ_HIGHLAT &   #CHISQ_MAP '   #CHISQ_RMS (   #CHISQ_BAND )   #NU_BAND *   #CHAIN +   #ITER ,             
                                 "                     
                                  #     8             #GENVEC              
                                $                                                     %     
                                                 &     
                      �                           '                   
               & p                   & p                                                                                    (     
                      �                           )                   
               & p                                                         �                           *                   
               & p                                                    
                                +                     
                                ,           %         @                               -                 
   
       #HANDLE .   #X_IN /   #LNL 0   #PRIOR 2   #STATUS 3   #N_EVAL 4   #LNL_IN 5   #OPTIMIZE 6   #USE_PRECOMPUTED_GRID 7   #TOLERANCE_ 8                                              .                     #PLANCK_RNG              
     �                           /                   
              & p                                          %         @                                0                    
       #X 1             
                                1     
                                                2                   
     p          p            p                                                                    3                                                      4                      
    �                           5                   
              & p                                                    
                                6                     
                                7                     
                                8     
               @ @                              9                   
                &                   &                                           #         @                                   :                    #PARAMFILE ;             
@ @                             ;                    1 #         @                                   <                    #HANDLE =   #S >   #NOISEAMP ?             
D @                               =                     #PLANCK_RNG              
  @                               >     8             #GENVEC              
D     �                           ?                   
               & p                                          %         @  @                            @                    
       #X A                                                        
  @                              A     
         �   2      fn#fn (   �   !   b   uapp(COMM_NOISESAMP_MOD     �   @   J   COMM_N_MULT_MOD    3  @   J   ARS_MOD    s  @   J   COMM_MP_MOD "   �  �       PLANCK_RNG+RNGMOD $   4  H   a   PLANCK_RNG%X+RNGMOD $   |  H   a   PLANCK_RNG%Y+RNGMOD $   �  H   a   PLANCK_RNG%Z+RNGMOD $     H   a   PLANCK_RNG%W+RNGMOD '   T  H   a   PLANCK_RNG%GSET+RNGMOD (   �  H   a   PLANCK_RNG%EMPTY+RNGMOD '   �  w       GENVEC+COMM_GENVEC_MOD /   [  �   a   GENVEC%CMB_AMP+COMM_GENVEC_MOD .     �   a   GENVEC%FG_AMP+COMM_GENVEC_MOD 0   �  �   a   GENVEC%TEMP_AMP+COMM_GENVEC_MOD '   w  P       #UNLPOLY+ISO_C_BINDING "   �  p       I4B+HEALPIX_TYPES !   7  p       DP+HEALPIX_TYPES $   �  @       NUMBAND+COMM_BP_MOD &   �  a       INT2STRING+COMM_UTILS .   H  @   a   INT2STRING%INTEGER+COMM_UTILS -   �  L   a   INT2STRING%STRING+COMM_UTILS )   �  �       GET_PARAMETER+COMM_UTILS 1   �	  L   a   GET_PARAMETER%PARFILE+COMM_UTILS 1   �	  L   a   GET_PARAMETER%PARNAME+COMM_UTILS 1   9
  @   a   GET_PARAMETER%PAR_INT+COMM_UTILS 2   y
  L   a   GET_PARAMETER%PAR_CHAR+COMM_UTILS 4   �
  L   a   GET_PARAMETER%PAR_STRING+COMM_UTILS 0     @   a   GET_PARAMETER%PAR_SP+COMM_UTILS 0   Q  @   a   GET_PARAMETER%PAR_DP+COMM_UTILS 1   �  @   a   GET_PARAMETER%PAR_LGT+COMM_UTILS 5   �  @   a   GET_PARAMETER%PAR_PRESENT+COMM_UTILS .     L   a   GET_PARAMETER%DESC+COMM_UTILS *   ]  �       COMPUTE_CHISQ+COMM_MP_MOD 1   D  @   a   COMPUTE_CHISQ%MAP_ID+COMM_MP_MOD 0   �  T   a   COMPUTE_CHISQ%COEFF+COMM_MP_MOD 7   �  @   a   COMPUTE_CHISQ%OUTPUT_STATS+COMM_MP_MOD 8     @   a   COMPUTE_CHISQ%CHISQ_FULLSKY+COMM_MP_MOD 8   X  @   a   COMPUTE_CHISQ%CHISQ_HIGHLAT+COMM_MP_MOD 4   �  �   a   COMPUTE_CHISQ%CHISQ_MAP+COMM_MP_MOD 4   D  @   a   COMPUTE_CHISQ%CHISQ_RMS+COMM_MP_MOD 5   �  �   a   COMPUTE_CHISQ%CHISQ_BAND+COMM_MP_MOD 2     �   a   COMPUTE_CHISQ%NU_BAND+COMM_MP_MOD 0   �  @   a   COMPUTE_CHISQ%CHAIN+COMM_MP_MOD /   �  @   a   COMPUTE_CHISQ%ITER+COMM_MP_MOD +   $  �       SAMPLE_INVSAMP+INVSAMP_MOD 2   �  X   a   SAMPLE_INVSAMP%HANDLE+INVSAMP_MOD 0   R  �   a   SAMPLE_INVSAMP%X_IN+INVSAMP_MOD /   �  W      SAMPLE_INVSAMP%LNL+INVSAMP_MOD 1   9  @   a   SAMPLE_INVSAMP%LNL%X+INVSAMP_MOD 1   y  �   a   SAMPLE_INVSAMP%PRIOR+INVSAMP_MOD 2     @   a   SAMPLE_INVSAMP%STATUS+INVSAMP_MOD 2   M  @   a   SAMPLE_INVSAMP%N_EVAL+INVSAMP_MOD 2   �  �   a   SAMPLE_INVSAMP%LNL_IN+INVSAMP_MOD 4     @   a   SAMPLE_INVSAMP%OPTIMIZE+INVSAMP_MOD @   ]  @   a   SAMPLE_INVSAMP%USE_PRECOMPUTED_GRID+INVSAMP_MOD 6   �  @   a   SAMPLE_INVSAMP%TOLERANCE_+INVSAMP_MOD    �  �       N_PRIOR )   �  W       INITIALIZE_NOISESAMP_MOD 3   �  L   a   INITIALIZE_NOISESAMP_MOD%PARAMFILE !   $  i       SAMPLE_NOISEAMPS (   �  X   a   SAMPLE_NOISEAMPS%HANDLE #   �  T   a   SAMPLE_NOISEAMPS%S *   9  �   a   SAMPLE_NOISEAMPS%NOISEAMP    �  �       LNL_NOISEAMP    K  @   a   LNL_NOISEAMP%X 