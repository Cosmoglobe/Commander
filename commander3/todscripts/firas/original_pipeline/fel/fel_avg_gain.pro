pro fel_avg_gain,chanscan
;
;  Average the Gain JCJ's
;
restore,'csdr$firas_in:jcj_spec.'+chanscan
nsize = size(jcj(*,*,*))
ag = complexarr(nsize(1)-1,nsize(2))
uj = nsize(1)-2
uk = nsize(2)-1
for j=0,uj do for k=0,uk do $
    ag(j,k) = total(jcj(j+1,k,*)*conj(jcj(0,k,*)))/total(abs(jcj(0,k,*))^2)
save,file='csdr$firas_out:jcj_spec.'+chanscan,ao,ag,jcj,f,idh,idl,pixel,gl

return
end
