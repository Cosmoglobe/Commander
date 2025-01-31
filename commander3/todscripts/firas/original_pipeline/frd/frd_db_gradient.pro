;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FRD_DB_GRADIENT creates the FIRAS version of the DIRBE gradients.
;
;DESCRIPTION:
;     An IDL procedure to convert the DIRBE gradients from the input file
;     CSDR$FIRAS_REF:FEX_DB_BANDS.DAT into the FIRAS gradients output file 
;     CSDR$FIRAS_OUT:FEX_GRAD.DAT.
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS:
;     None
;
;WARNINGS:
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_REF   directory containing dataset FEX_DB_BANDS.DAT
;	CSDR$FIRAS_OUT	 directory where FEX_GRAD.DAT file will be written
;
;EXAMPLE:
;     $ UIDL
;     setlog,'csdr$firas_ref','e'
;     setlog,'csdr$firas_out','f'
;     frd_db_gradient
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
;     Modified by S. Brodd, HSTX, 5/6/97, for Pass 4, SPR 12348
;-
;______________________________________________________________________________
;
pro frd_db_gradient
;
; Create intermediate arrays P, U, V, and W.
;
p=intarr(9,6144)

v=coorconv(indgen(6144),infmt='p',outfmt='u',inco='f',outco='g')

for i=0,6143 do begin
    k=reverse(sort(v#reform(v(i,*))))
    p(*,i)=k(0:8)
end

r=sqrt(total(double(v(*,0:1))^2,2))
u=double(transpose(v)) 
w=u

for i=0,6143 do begin
    u(*,i)=[reform(-u(0:1,i)*(u(2,i)/r(i))),r(i)]
end

for i=0,6143 do begin
    w(*,i)=[w(1,i)*u(2,i)-w(2,i)*u(1,i), $
    w(2,i)*u(0,i)-w(0,i)*u(2,i),         $
    w(0,i)*u(1,i)-w(1,i)*u(0,i)]
end
;
; Create intermediate array A.
;
one=fltarr(9)+1. 
a=dblarr(6,9,6144)

for i=0,6143 do begin 
    ud=v(p(*,i),*)#u(*,i) 
    wd=v(p(*,i),*)#w(*,i)
    d=transpose([[one],[ud],[wd],[ud^2],[wd^2],[ud*wd]])
    a(*,*,i)=invert(d#transpose(d))#d 
end
;
; Restore DIRBE gradients file.
;
restore,'csdr$firas_ref:fex_db_bands.dat'
;
; Calculate A8, A9, and A10.
;
a8=fltarr(6,6144) 
a9=a8 
a10=a8

for i=0,6143 do begin
    a8(*,i)=a(*,*,i)#b8(p(*,i))
    a9(*,i)=a(*,*,i)#b9(p(*,i))
    a10(*,i)=a(*,*,i)#b10(p(*,i))
end
; 
; Save output file.
;
save,file='csdr$firas_out:fex_grad.dat',a8,a9,a10,u,w

end
