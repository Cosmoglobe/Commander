pro fel_avg_off,chanscan
;
;  Average the Offset JCJ's
;
restore,'csdr$firas_in:jcj_spec.'+chanscan
nsize = size(jcj(*,*,idh))
ao = complexarr(nsize(1)-1,nsize(2))
uj = nsize(1)-2
uk = nsize(2)-1
norm = float(nsize(3))
for j=0,uj do for k=0,uk do ao(j,k)=total(jcj(j+1,k,idh))/norm
save,file='csdr$firas_out:jcj_spec.'+chanscan,ao,jcj,f,idh,idl,pixel,gl

return
end
