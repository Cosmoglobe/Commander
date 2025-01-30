 Function FLA_Apod_Rec,  Mode,                    $
                         SciMode    = SciMode,    $
                         Linearized = Linearized, $
                         FakeIt     = FakeIt
;
; Author: K. Turpie (based on code by Alice R. Trenholme)
;         General Sciences Corporation
;
;         Reference look-up tables used are defined by the FLA_Init_FIR routine.
;         They have the following values and meanings:
;
;         speeds        = 0  slow
;                       = 1  fast
;
;         lengths       = 0  short
;                       = 1  long
;
;         channels      = 0  RH
;                       = 1  RL
;                       = 2  LH
;                       = 3  LL
;
;         NGroups       = the correct number of adds/group
;
;         SciMode       = microprocessor mode (4 is the usual)
;
;         Linearized    = 0  no
;                       = 1  yes
;
;         FakeIt        = 0  no
;                       = 1  yes
;
;  Changes:
;
;    Converted to FLA.
;    Gene Eplee, 14 September 1993, SER 11413
;
;
 On_Error, 2 ; return to caller if an error occurs
;
 Common FIRAS_Modes,   Modes,      NGroups,     Peaks,   $
                       Channels,   Lengths,     Speeds,  $
                       Nyquists
;
 Common FIRAS_Gain,    GainMode,   Gain,                 $
                       FreqArray,  Bin0,        BinN
;
 Common FIRAS_Apod,    Apod
;
;Begin,
;
    If (N_Elements( FakeIt )    eq 0) then FakeIt     = 0
    If (N_Elements( SciMode )   eq 0) then SciMode    = 4
    If (N_Elements( Linearize ) eq 0) then Linearized = 1
;
    ModeIndx = Where( StrUpcase( Mode ) eq Modes )
    Ptr      = ModeIndx(0)
    If (Ptr lt 0) then Message, 'ERROR: ' + Mode + ' is an invalid mode.'
;
    If ((Channels(Ptr) eq 0) or (Channels(Ptr) eq 2)) then Ch  = 1 else Ch  = 0
    If ((SciMode       eq 1) or (SciMode       eq 3)) then Sci = 0 else Sci = 1
;
    If (FakeIt eq 1) then IRec = 384                       $
                     else IRec = 32 * (NGroups(Ptr)-1)   + $
                                 16 *  Speeds(Ptr)       + $
                                  8 *  Lengths(Ptr)      + $
                                  4 *  Ch                + $
                                  2 *  Sci               + $
                                       Linearized 
    Return, IRec
 End
