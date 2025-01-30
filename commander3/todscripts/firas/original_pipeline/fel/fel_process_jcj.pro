pro fel_process_jcj,chanscan
;
;  Run Read_JCJ and save the results.
;
fel_read_jcj,chanscan,pixel,f,jcj
gl = coorconv(pixel,infmt='P',outfmt='L',inco='R6',outco='G')
idl = where(abs(gl(*,1)) lt 10.0)
idh = where(abs(gl(*,1)) gt 10.0)
save,filename='csdr$firas_out:jcj_spec.'+chanscan,pixel,f,jcj,gl,idh,idl

return
end
