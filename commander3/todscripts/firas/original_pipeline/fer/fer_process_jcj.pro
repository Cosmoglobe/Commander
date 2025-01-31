pro fer_process_jcj,smode

;
;  Run Read_JCJ and save the results.
;
fer_read_jcj,smode,pixel,model,f,jcj
gl = coorconv(pixel,infmt='P',outfmt='L',inco='R6',outco='G')
idl = where(abs(gl(*,1)) lt 10.0)
idh = where(abs(gl(*,1)) gt 10.0)
save,filename='csdr$firas_out:jcj_spec.'+smode,pixel,model,f,jcj,gl,idh,idl

return
end
