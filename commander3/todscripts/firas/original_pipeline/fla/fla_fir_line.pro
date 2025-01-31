;-------------------------------------------------------------------------------
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    FLA_FIR_LINE is an IDL function that generates a FIRAS synthetic line
;    profile for a given channel, scan mode, and line frequency.
;
; CALLING SEQUENCE:
; Function FLA_FIR_LINE, Freq, Mode, Solution,        $
;                              Phase = Phase,         $
;                              Amplitude = Amplitude, $
;                              Start = Start,         $
;                              Stop  = Stop,          $
;                              Uncal = Uncal,         $
;                              gauss = gauss,         $
;                              sinc = sinc
;
; ARGUMENTS (I=input, []=optional):
;    freq         I  frequency of line profile
;    mode         I  instrument channel and scan mode
;    solution     I  calibration model solution filename extension
;    [phase]      I  phase shift applied to line profile (Default=0)
;    [amplitude]  I  amplitude of line profile (Default=1)
;    [start]      I  Initial frequency bin of profile (Default=1)
;    [stop]       I  Final frequency bin of profile (Default=167 for high freq)
;    [uncal]      I  Uncalibrated prorfile (Default = calibrated)
;    gauss        I  Gaussian self-apodization function width
;                         The standard value is 1.0/156.0
;    sinc         I  Sinc self-apodizaiton function width
;                         The standard value is 0.0
;
;  WARNING:
;
;    Subroutines of this routine use the logical pointer CSDR$FIRAS_REF to point
;    to the reference datasets FEX_NYQUIST and FEX_APOD and the logical pointer
;    CSDR$FIRAS_IN to point to the gain function file FEX_GAIN_CCSS.
;
;  EXAMPLE:
;    To generate a line profile for the 63 icm [CII] for the RHSS scan mode
;    using the F16_93HYBRID calibration model solution and a gaussian
;    self-apodization with a width of 1/156, the invocation is:
;
;    UIDL> CII = fla_fir_line (63.395, 'rhss', 'f16_93hybrid', $
;                              gauss=1.0/156.0, sinc=0.0)
;
;#
;
;  Author:  Kevin Turpie
;           General Sciences Corp.
;  Converted to FLA by:  Gene Eplee
;                        General Sciences Corp.
;                        14 September 1993
;                        SER 11413
;
;  Changes:
;
;    Added Self-apodization.
;    Joel Gales, ARC
;
;    Added calibration model solution file name extension.
;    Gene Eplee, GSC, 14 September 1993.
;
;    Converted ifg sample array to IDL indexing.
;    Gene Eplee, GSC, 21 September 1993.
;
;-
;-------------------------------------------------------------------------------
;
 Function FLA_FIR_LINE, Freq, Mode, Solution,        $
                              Phase = Phase,         $
                              Amplitude = Amplitude, $
                              Start = Start,         $
                              Stop  = Stop,          $
                              Uncal = Uncal,         $
			      gauss = gauss,         $
			      sinc = sinc
;
;
 Common FIRAS_Modes,   Modes,      NGroups,     Peaks,   $
                       Channels,   Lengths,     Speeds,  $
                       Nyquists

 Common FIRAS_Gain,    GainMode,   Gain,                 $
                       FreqArray,  Bin0,        BinN

 Common FIRAS_Apod,    Apod

    FLA_Init_FIR, Mode, Solution

    F     = FreqArray           ; FreqArray, Bin0, and BinN come
    Start = Bin0                ; from the common block FIRAS_GAIN
    Stop  = BinN                ;

    If (N_Elements( Phase ) eq 0) then Phase = 0.0D0
    If (N_Elements( Amplitude ) eq 0) then Amplitude = 1.0
    If (Amplitude eq 0) then Amplitude = 1.0
    If (N_Elements( Gauss ) eq 0) then begin
        Message, "Keyword GAUSS must be set to some value."
    EndIf
    If (N_Elements( Sinc ) eq 0) then begin
        Message, "Keyword SINC must be set to some value."
    EndIf


    ModeIndx = Where( StrUpcase( Mode ) eq Modes )
    If (ModeIndx(0) lt 0) then Message, Mode + " is an invalid mode."

;                Mode    Freq    Length   Speed
;              +---------------------------------
;            0 | HSS     high    short    slow
;            1 | HSF     high    short    fast
;            2 | HLS     high    long     slow
;            3 | HLF     high    long     fast
;            4 | LSS     low     short    slow
;            5 | LSF     low     short    fast
;            6 | LLS     low     long     slow
;            7 | LLF     low     long     fast
;              |
;

    nyquist = Nyquists( ModeIndx(0) mod 8 )
	  ;  from common block FIRAS_MODES in file fir_parm.inc


    Peak    = Peaks( ModeIndx(0) )


;   Create line signal and apodization function :

;   LineIFG models the IFG for a single line at the input frequency Freq,
;   centered at the point of zero path difference (given by Peak), and
;   sampled at 1/(2*Nyquist) intervals.

;   Compute angular frequency (computed) and define independent variable :

    w = Double( Freq * !DPi ) / Double( Nyquist )

    x = DIndGen( 512 )

    u = 0.5 * w * sinc * (x - peak) + 1e-12
    v = gauss * w * (x - peak) / (2 * !dpi)
    apod_self = 0*x + 1
	IF (sinc NE 0) THEN apod_self = sin(u)/u

	IF (gauss NE 0) THEN apod_self = exp(-v*v)

    signal = Amplitude * Cos( w * (x - peak) + Phase ) * apod_self

;   Apodize line signal :                  ; Apod comes from the common block
    Signal  = Signal * Apod                ; FIRAS_APOD in fir_parm.inc

;   Rotate synthetic line signal :
    Signal  = Shift( Signal, -Peak )

;   Creat synthetic line profile in the frequency domain :
    Trnsfrm = FFT( Signal, -1 )

    Calibrate = (not Keyword_Set( Uncal ))

    If (Calibrate) then begin
       Center = Fix( Freq * 256.0 / Nyquist + 0.5 ) - Start   ; Gain comes from
       Gv     = Abs( Gain( Center ) )                         ; the common block
       Normlz = 512. * Abs( Gain ) / (Gv * Nyquist)           ; FIRAS_GAIN in
    EndIf                                                     ; fir_parm.inc

;   Truncate unwanted channels and calibrate the line with the gain function :

    If (Calibrate) then Line = Trnsfrm(Start:Stop) * Normlz $
                   else Line = Trnsfrm(Start:Stop)

    Return, Line
 End
