  >4  q   k820309    l          18.0        ­¹Ù[                                                                                                          
       math_tools.f90 MATH_TOOLS                                                     
                                                              u #INVERT_MATRIX_DPC    #INVERT_MATRIX_DP    #INVERT_MATRIX_SP                                                           u #INVERT_MATRIX_WITH_MASK_DPC    #INVERT_MATRIX_WITH_MASK_DP                                                           u #CONVERT_FRACT2SIGMA_SP    #CONVERT_FRACT2SIGMA_DP                                                 	                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                   
                  ÿÿÿÿÿÿï                                                          HUGE                                                  
                
          )       -DTû!	@        3.141592653589793238462643383279502884197                                                 
                
          *       mBP×Ò?        0.2820947917738781434740397257803862929220#         @      X                                                 #MATRIX              
D@                                                                & p                  & p                                          #         @      X                                                 #MATRIX    #CHOLESKY    #STATUS              
D@                                                 
               & p                  & p                                                    
 @                                                   F @                                          #         @      X                                                 #MATRIX              
D@                                                 	               & p                  & p                                          #         @      X                                                 #MATRIX              
D@                                                                & p                  & p                                          #         @      X                                                 #MATRIX              
D@                                                 
               & p                  & p                                          #         @      X                                                 #SIGMA    #FRACT              D                                     	                 
                                      	      #         @      X                                                 #SIGMA    #FRACT              D                                     
                 
                                      
      #         @                                                       #A    #ACROSS               
 @                                                 
 
             &                   &                                                     D                                                    
               &                   &                                           #         @                                   !                    #A "   #ACROSS #             
 @                              "                   
              &                   &                                                     D @                              #                   
               &                   &                                           #         @                                   $                    #MATRIX %             
D@                              %                   	                & p                  & p                                          #         @                                  &                    #MATRIX '   #EIGENVALS (   #EIGENVECTORS )             
                                 '                   
 $             & p                  & p                                                    D@                              (                   
 %              & p                                                    D @                              )                   
 &              & p                  & p                                          #         @                                   *                    #A +   #EIGENVALS ,             
                                 +                   
 +             & p                  & p                                                    D@                              ,                   
 ,              & p                                          #         @                                   -                    #MATRIX .   #EIGENVALS /   #EIGENVECTORS 0   #STATUS 1             @                              .                   
 1              &                   &                                                     D @                              /                   
 2              &                                                     D @                              0                   
 3              &                   &                                                     D @                              1            #         @                                   2                    #A 3   #POW 4             
D                               3                   
 6              & p                  & p                                                    
                                 4     
      #         @                                   5                    #A 6   #X 7   #B 8                                             6                    ;              &                   &                                                     D                                7                    <              &                                                  0  @                              8                    =              &                                           #         @                                   9                    #MATRIX :   #THRESHOLD ;             
D@                              :                   
 A              & p                  & p                                                    
                                 ;     
      #         @                                   <                    #A =   #X >   #B ?             
                                 =                   
 E             & p                  & p                                                    D@                              >                   
 G              & p                                                    
                                 ?                   
 F             & p                                          #         @                                   @                    #NLMAX A   #M B   #THETA C   #PLM D             
                                 A                     
  @                              B                     
  @                              C     
               D                                D                    
 L    p           & p         5  p        r A         5  p        r A   p         p                          #         @                                   E                    #FRAC F   #SIGMA G             D                                F     	                 
  @                              G     	      %         @                               H                    
       #X I             
  @                              I     
      #         @                                   J                    #A K   #L L   #IERR M             
                                K                   
 M             &                   &                                                     D @                              L                   
 N              &                   &                                                     F @                              M            #         @                                   N                    #A O   #IERR P             
D@                              O                   
 P              &                   &                                                     F @                              P            #         @                                   Q                    #L R   #B S   #X T             
@@                              R                   
 R             &                   &                                                     
                                 S                   
 S             &                                                     D                                T                   
 T              &                                           #         @                                   U                    #A V   #B W   #C X             
@@                              V                   
 V             &                   &                                                     
@ @                              W                   
 W             &                   &                                                     D @                              X                   
 X              &                   &                                           #         @                                   Y                    #A Z   #B [   #C \             
@@                              Z                   
 Y             &                   &                                                     
@ @                              [                   
 Z             &                   &                                                     D @                              \                   
 [              &                   &                                           %         @                                ]                    
       #ARRAY ^             @                              ^                   
 \              &                                           %         @                                 _                    
       #ARRAY `             D@                              `                   
 ]              &                                           #         @                                   a                    #A b   #B c   #C d   #R e   #U f             
                                 b                   
 ^             &                                                  0  
 @                              c                   
 _             &                                                     
                                 d                   
 `             &                                                     
                                 e                   
 a             &                                                     D                                f                   
 b              &                                           #         @                                   g                    #A1 h   #B1 i   #C1 j   #A2 k   #B2 l   #C2 m             
D                                h     
                 
D                                i     
                 
D                                j     
                 
D                                k     
                 
D                                l     
                 
D                                m     
              "      fn#fn    Â   @   J   HEALPIX_TYPES "            gen@INVERT_MATRIX ,            gen@INVERT_MATRIX_WITH_MASK (     x       gen@CONVERT_FRACT2SIGMA "   ~  p       DPC+HEALPIX_TYPES "   î  p       I4B+HEALPIX_TYPES !   ^  p       DP+HEALPIX_TYPES "   Î  p       LGT+HEALPIX_TYPES !   >  p       SP+HEALPIX_TYPES "   ®  p       I8B+HEALPIX_TYPES %     p       MAX_DP+HEALPIX_TYPES #     =       HUGE+HEALPIX_TYPES !   Ë         PI+HEALPIX_TYPES (   d         SQ4PI_INV+HEALPIX_TYPES "   þ  T       INVERT_MATRIX_DPC )   R  ¬   a   INVERT_MATRIX_DPC%MATRIX !   þ  n       INVERT_MATRIX_DP (   l  ¬   a   INVERT_MATRIX_DP%MATRIX *   	  @   a   INVERT_MATRIX_DP%CHOLESKY (   X	  @   a   INVERT_MATRIX_DP%STATUS !   	  T       INVERT_MATRIX_SP (   ì	  ¬   a   INVERT_MATRIX_SP%MATRIX ,   
  T       INVERT_MATRIX_WITH_MASK_DPC 3   ì
  ¬   a   INVERT_MATRIX_WITH_MASK_DPC%MATRIX +     T       INVERT_MATRIX_WITH_MASK_DP 2   ì  ¬   a   INVERT_MATRIX_WITH_MASK_DP%MATRIX '     ^       CONVERT_FRACT2SIGMA_SP -   ö  @   a   CONVERT_FRACT2SIGMA_SP%SIGMA -   6  @   a   CONVERT_FRACT2SIGMA_SP%FRACT '   v  ^       CONVERT_FRACT2SIGMA_DP -   Ô  @   a   CONVERT_FRACT2SIGMA_DP%SIGMA -     @   a   CONVERT_FRACT2SIGMA_DP%FRACT '   T  [       COMPUTE_PSEUDO_INVERSE )   ¯  ¤   a   COMPUTE_PSEUDO_INVERSE%A .   S  ¤   a   COMPUTE_PSEUDO_INVERSE%ACROSS (   ÷  [       COMPUTE_PSEUDO_INVERSE2 *   R  ¤   a   COMPUTE_PSEUDO_INVERSE2%A /   ö  ¤   a   COMPUTE_PSEUDO_INVERSE2%ACROSS +     T       INVERT_MATRIX_WITH_MASK_SP 2   î  ¬   a   INVERT_MATRIX_WITH_MASK_SP%MATRIX (     u       GET_EIGEN_DECOMPOSITION /     ¬   a   GET_EIGEN_DECOMPOSITION%MATRIX 2   »     a   GET_EIGEN_DECOMPOSITION%EIGENVALS 5   K  ¬   a   GET_EIGEN_DECOMPOSITION%EIGENVECTORS     ÷  ^       GET_EIGENVALUES "   U  ¬   a   GET_EIGENVALUES%A *        a   GET_EIGENVALUES%EIGENVALS             EIGEN_DECOMP $     ¤   a   EIGEN_DECOMP%MATRIX '   ¶     a   EIGEN_DECOMP%EIGENVALS *   B  ¤   a   EIGEN_DECOMP%EIGENVECTORS $   æ  @   a   EIGEN_DECOMP%STATUS '   &  X       COMPUTE_HERMITIAN_ROOT )   ~  ¬   a   COMPUTE_HERMITIAN_ROOT%A +   *  @   a   COMPUTE_HERMITIAN_ROOT%POW    j  ]       SOLVE_SYSTEM    Ç  ¤   a   SOLVE_SYSTEM%A    k     a   SOLVE_SYSTEM%X    ÷     a   SOLVE_SYSTEM%B '     c       INVERT_SINGULAR_MATRIX .   æ  ¬   a   INVERT_SINGULAR_MATRIX%MATRIX 1     @   a   INVERT_SINGULAR_MATRIX%THRESHOLD "   Ò  ]       SOLVE_SYSTEM_REAL $   /  ¬   a   SOLVE_SYSTEM_REAL%A $   Û     a   SOLVE_SYSTEM_REAL%X $   k     a   SOLVE_SYSTEM_REAL%B $   û  n       COMP_NORMALISED_PLM *   i   @   a   COMP_NORMALISED_PLM%NLMAX &   ©   @   a   COMP_NORMALISED_PLM%M *   é   @   a   COMP_NORMALISED_PLM%THETA (   )!  ä   a   COMP_NORMALISED_PLM%PLM #   "  ]       CONVERT_SIGMA2FRAC (   j"  @   a   CONVERT_SIGMA2FRAC%FRAC )   ª"  @   a   CONVERT_SIGMA2FRAC%SIGMA    ê"  W       CORR_ERF    A#  @   a   CORR_ERF%X #   #  `       CHOLESKY_DECOMPOSE %   á#  ¤   a   CHOLESKY_DECOMPOSE%A %   $  ¤   a   CHOLESKY_DECOMPOSE%L (   )%  @   a   CHOLESKY_DECOMPOSE%IERR *   i%  Y       CHOLESKY_DECOMPOSE_SINGLE ,   Â%  ¤   a   CHOLESKY_DECOMPOSE_SINGLE%A /   f&  @   a   CHOLESKY_DECOMPOSE_SINGLE%IERR    ¦&  ]       CHOLESKY_SOLVE !   '  ¤   a   CHOLESKY_SOLVE%L !   §'     a   CHOLESKY_SOLVE%B !   3(     a   CHOLESKY_SOLVE%X    ¿(  ]       MATMUL_SYMM    )  ¤   a   MATMUL_SYMM%A    À)  ¤   a   MATMUL_SYMM%B    d*  ¤   a   MATMUL_SYMM%C    +  ]       MATMUL_GEN    e+  ¤   a   MATMUL_GEN%A    	,  ¤   a   MATMUL_GEN%B    ­,  ¤   a   MATMUL_GEN%C    Q-  [       MEAN    ¬-     a   MEAN%ARRAY    8.  [       VARIANCE    .     a   VARIANCE%ARRAY    /  k       TRIDAG    /     a   TRIDAG%A    0     a   TRIDAG%B    ¢0     a   TRIDAG%C    .1     a   TRIDAG%R    º1     a   TRIDAG%U    F2  x       MOV3    ¾2  @   a   MOV3%A1    þ2  @   a   MOV3%B1    >3  @   a   MOV3%C1    ~3  @   a   MOV3%A2    ¾3  @   a   MOV3%B2    þ3  @   a   MOV3%C2 