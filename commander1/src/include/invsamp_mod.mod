  �  B   k820309    l          18.0        3�\                                                                                                          
       InvSamp_mod.f90 INVSAMP_MOD              N_SPLINE DELTA_LNL MIN_NUM_ACTIVE_POINT TOLERANCE                                                     
                            @                              
                                                           
                         @                                '                     #X    #Y    #Z    #W    #GSET 	   #EMPTY 
                � $                                                             � $                                                            � $                                                            � $                                                            � $                             	               
                � $                             
                                                                                                                                                                                                                                                                                                                                                                                                                                          %         @                                                   
       #HANDLE              
                                                      #PLANCK_RNG    #         @                                                       #S    #X    #Y    #BOUNDARY    #REGULAR    #LINEAR              
                                      �               #SPLINE_TYPE              
                                                   
              &                                                     
                                                    
              &                                                     
                                                   
    p          p            p                                    
                                                     
                                           #         @                                                      #X    #Y    #YP1    #YPN    #Y2              
                                                   
 
             &                                                     
                                                    
              &                                                     
                                      
                
                                      
                                                                   
               &                                           %         @                                                     
       #S     #X !             
                                        �              #SPLINE_TYPE              
                                 !     
      %         @                               "                    
       #XA #   #YA $   #Y2A %   #X &             
                                #                   
              &                                                     
                                 $                   
              &                                                     
                                 %                   
              &                                                     
                                 &     
                                                  '                                       �              1000%         @                                (                 
   
       #HANDLE )   #X_IN *   #LNL +   #PRIOR -   #STATUS .   #N_EVAL /   #LNL_IN 0   #OPTIMIZE 1   #USE_PRECOMPUTED_GRID 2   #TOLERANCE_ 3             D @                               )                     #PLANCK_RNG              
 @   �                           *                   
              & p                                          %         @                               +                    
       #X ,                                        
                                ,     
                 @                              -                   
     p          p            p                                    F @                              .                      F @                              /                      
@   �                           0                   
              & p                                                    
 @                              1                     
 @                              2                     
 @                              3     
      #         @                                  4                    #X_NEW 5   #Y_NEW 6   #X 7   #Y 8   #N 9   #STAT :             
                                 5     
                
                                 6     
                
D                                7                   
               &                                                     
D                                8                   
               &                                                     
D                                9                      
D                                :                              @               @                '�                    #X ;   #Y <   #Y2 =   #BOUNDARY >   #REGULAR ?   #LINEAR @             �                              ;                              
            &                                                      �                              <            H                 
            &                                                      �                              =            �                 
            &                                                        �                              >            �                 
  p          p            p                                       �                              ?     �                          �                              @     �                �   $      fn#fn !   �   B   b   uapp(INVSAMP_MOD      @   J   HEALPIX_TYPES    F  @   J   RNGMOD    �  @   J   SPLINE_1D_MOD "   �  �       PLANCK_RNG+RNGMOD $   G  H   a   PLANCK_RNG%X+RNGMOD $   �  H   a   PLANCK_RNG%Y+RNGMOD $   �  H   a   PLANCK_RNG%Z+RNGMOD $     H   a   PLANCK_RNG%W+RNGMOD '   g  H   a   PLANCK_RNG%GSET+RNGMOD (   �  H   a   PLANCK_RNG%EMPTY+RNGMOD "   �  p       I4B+HEALPIX_TYPES !   g  p       DP+HEALPIX_TYPES "   �  p       LGT+HEALPIX_TYPES !   G  p       SP+HEALPIX_TYPES     �  \       RAND_UNI+RNGMOD '     X   a   RAND_UNI%HANDLE+RNGMOD ,   k  �       SPLINE_SIMPLE+SPLINE_1D_MOD .   �  Y   a   SPLINE_SIMPLE%S+SPLINE_1D_MOD .   H  �   a   SPLINE_SIMPLE%X+SPLINE_1D_MOD .   �  �   a   SPLINE_SIMPLE%Y+SPLINE_1D_MOD 5   `  �   a   SPLINE_SIMPLE%BOUNDARY+SPLINE_1D_MOD 4   �  @   a   SPLINE_SIMPLE%REGULAR+SPLINE_1D_MOD 3   4	  @   a   SPLINE_SIMPLE%LINEAR+SPLINE_1D_MOD +   t	  p       SPLINE_PLAIN+SPLINE_1D_MOD -   �	  �   a   SPLINE_PLAIN%X+SPLINE_1D_MOD -   p
  �   a   SPLINE_PLAIN%Y+SPLINE_1D_MOD /   �
  @   a   SPLINE_PLAIN%YP1+SPLINE_1D_MOD /   <  @   a   SPLINE_PLAIN%YPN+SPLINE_1D_MOD .   |  �   a   SPLINE_PLAIN%Y2+SPLINE_1D_MOD ,     ^       SPLINT_SIMPLE+SPLINE_1D_MOD .   f  Y   a   SPLINT_SIMPLE%S+SPLINE_1D_MOD .   �  @   a   SPLINT_SIMPLE%X+SPLINE_1D_MOD +   �  p       SPLINT_PLAIN+SPLINE_1D_MOD .   o  �   a   SPLINT_PLAIN%XA+SPLINE_1D_MOD .   �  �   a   SPLINT_PLAIN%YA+SPLINE_1D_MOD /   �  �   a   SPLINT_PLAIN%Y2A+SPLINE_1D_MOD -     @   a   SPLINT_PLAIN%X+SPLINE_1D_MOD &   S  t       INVSAMP_MAX_NUM_EVALS    �  �       SAMPLE_INVSAMP &   �  X   a   SAMPLE_INVSAMP%HANDLE $   �  �   a   SAMPLE_INVSAMP%X_IN #   �  r      SAMPLE_INVSAMP%LNL %   �  @   a   SAMPLE_INVSAMP%LNL%X %   7  �   a   SAMPLE_INVSAMP%PRIOR &   �  @   a   SAMPLE_INVSAMP%STATUS &     @   a   SAMPLE_INVSAMP%N_EVAL &   K  �   a   SAMPLE_INVSAMP%LNL_IN (   �  @   a   SAMPLE_INVSAMP%OPTIMIZE 4     @   a   SAMPLE_INVSAMP%USE_PRECOMPUTED_GRID *   [  @   a   SAMPLE_INVSAMP%TOLERANCE_ *   �  }       UPDATE_INVSAMP_SAMPLE_SET 0     @   a   UPDATE_INVSAMP_SAMPLE_SET%X_NEW 0   X  @   a   UPDATE_INVSAMP_SAMPLE_SET%Y_NEW ,   �  �   a   UPDATE_INVSAMP_SAMPLE_SET%X ,   $  �   a   UPDATE_INVSAMP_SAMPLE_SET%Y ,   �  @   a   UPDATE_INVSAMP_SAMPLE_SET%N /   �  @   a   UPDATE_INVSAMP_SAMPLE_SET%STAT *   0  �       SPLINE_TYPE+SPLINE_1D_MOD ,   �  �   a   SPLINE_TYPE%X+SPLINE_1D_MOD ,   Q  �   a   SPLINE_TYPE%Y+SPLINE_1D_MOD -   �  �   a   SPLINE_TYPE%Y2+SPLINE_1D_MOD 3   y  �   a   SPLINE_TYPE%BOUNDARY+SPLINE_1D_MOD 2     H   a   SPLINE_TYPE%REGULAR+SPLINE_1D_MOD 1   ]  H   a   SPLINE_TYPE%LINEAR+SPLINE_1D_MOD 