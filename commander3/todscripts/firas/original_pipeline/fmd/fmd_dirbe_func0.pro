;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_DIRBE_FUNC0 creates the individual channel/scan mode and the combined
;  channel/scan mode versions of the DIRBE gradients.
;
;DESCRIPTION:
;     An IDL procedure to write the 13 channel/scan mode versions of the 
;     DIRBE gradients using the input files CSDR$FIRAS_REF:FEX_GRAD.DAT and 
;     CSDR$FIRAS_IN:ccss.ISS (where ccss is the channel/scan mode) into the 
;     output file CSDR$FIRAS_OUT:FMD_DIRBE_FUNC0_ccss.ISS.
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS:
;     chanscan -- string specifying the channel/scan mode; one of
;                 RHS, RHF, RLS, RSF, RLF, LHS, LHF, LLS, LSF, LLF,
;                 HIGH, LOWF, HRES
;
;WARNINGS:
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_REF   directory containing reference dataset FEX_GRAD.DAT
;       CSDR$FIRAS_IN    directory containing input dataset ccss.ISS
;	CSDR$FIRAS_OUT	 directory where output file will be written
;
;EXAMPLE:
;     $ UIDL
;     setlog,'csdr$firas_ref','e'
;     setlog,'csdr$firas_in','f'
;     setlog,'csdr$firas_out','g'
;     fmd_dirbe_func0, 'chanscan'
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;
;MODIFICATION HISTORY:
;     Written by D. Fixsen, HSTX
;     Modified by S. Brodd, HSTX, 5/6/97, for Pass 4, SER
;     Modified by K.Jensen, HSTX, 5/29/97, include FSL_IDX and SKY_LBL arrays.
;-
;______________________________________________________________________________
;
pro fmd_dirbe_func0, chanscan
;
; Restore reference and input datasets.
;
restore,'csdr$firas_ref:fex_grad.dat'
restore,'csdr$firas_in:'+chanscan+'.iss'
;
; Set unused restored variables to zero to avoid warning messages.
;
tm=0 & nifgs=0 & scan=0 & time=0 & st_sub=0 & f=0 & sl_weight_rat=0
n_rhs=0 & b_rhs=0 & l_rhs=0 & d_rhs=0 & n_rhf=0 & b_rhf=0 & l_rhf=0 & d_rhf=0
n_rls=0 & b_rls=0 & l_rls=0 & d_rls=0 & n_rsf=0 & b_rsf=0 & l_rsf=0 & d_rsf=0
n_lhs=0 & b_lhs=0 & l_lhs=0 & d_lhs=0 & n_lhf=0 & b_lhf=0 & l_lhf=0 & d_lhf=0
n_lls=0 & b_lls=0 & l_lls=0 & d_lls=0 & n_lsf=0 & b_lsf=0 & l_lsf=0 & d_lsf=0
n_rlf=0 & b_rlf=0 & l_rlf=0 & d_rlf=0 & n_llf=0 & b_llf=0 & l_llf=0 & d_llf=0
galcut=0 & cal_nifgs=0 & cal_tm=0 & xcal=0 & ical=0 & refh=0 & skyh=0 & dihd=0
sky_glitch=0 & cal_glitch=0 & cal_wgts=0 & sky_dihd=0 & sky_s0=0 & cal_s0=0
cal_lbl=0 & solution=0 & chan_label=0 & sky_idx=0 & cal_idx=0
;
; Calculate intermediate arrays B, L, D, DU, DW, Q, PW, and Z.
;
b=glat*!pi/180
l=glon*!pi/180

d=[[cos(b)*cos(l)],[cos(b)*sin(l)],[sin(b)]]
k=long(total(px gt -1))-1

du=fltarr(k+1) 
dw=du

for i=long(0),k do begin
    du(i)=total(d(i,*)*u(*,px(i)))
    dw(i)=total(d(i,*)*w(*,px(i)))
end

q=[[du*0+1],[du],[dw],[du^2],[dw^2],[du*dw]]

pw=dblarr(6144) 

for i=long(0),k do begin
    pw(px(i))=pw(px(i))+sky_wgts(i)
end

z=sky_wgts/pw(px) 
;
; Calculate G8, G9, and G10.
;
g8=fltarr(k+1) 
g9=g8 
g10=g8

for i=long(0),k do begin
    g8(i)=q(i,*)#a8(*,px(i))
    g9(i)=q(i,*)#a9(*,px(i))
    g10(i)=q(i,*)#a10(*,px(i))
end
;
; Calculate F8, F9 and F10;
;
f8=fltarr(6144)
f9=f8 
f10=f8

for i=long(0),k do begin
    f8(px(i))=f8(px(i))+z(i)*g8(i)
    f9(px(i))=f9(px(i))+z(i)*g9(i)
    f10(px(i))=f10(px(i))+z(i)*g10(i)
end
;
; Subtract Fs to get final Gs.
;
g8=g8-f8(px) 
g9=g9-f9(px) 
g10=g10-f10(px)
;
; Save output file.
;
save,file='csdr$firas_out:fmd_dirbe_func0_'+chanscan+'.iss', $
     g8,g9,g10,f8,f9,f10,fsl_idx,sky_lbl

end
