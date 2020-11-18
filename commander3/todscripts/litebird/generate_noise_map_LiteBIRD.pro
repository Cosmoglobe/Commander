pro make_beam, fwhm, lmax, filename
  b = gaussbeam(fwhm, lmax)
  cl2fits, b, filename
end

pro make_rms, sigma0, nside, filename, seed
  noise = fltarr(12*nside^2,3)
  noise[*,*] = sigma0*randomn(seed, 12*nside^2,3)
  write_tqu, filename, noise, /ring
  print, stddev(noise)
end

pro make_all_inputs

  ; Les inn tabell
  nside       = 256L
  lmax        = 3*nside
  freqs       = ['040', '050', '060', '068_1', '068_2', '078_1', '078_2', '089_1', '089_2', '100_1', '100_2', '119_1', '119_2', '140_1', '140_2', '166', '195_1', '195_2', '235', '280', '337', '402']
  sensitivity = []*nside/512L
  fwhms       = []

  seed = 736563
  for i = 0, n_elements(fwhms)-1 do begin
;     make_beam, fwhms [i], lmax, "LB_v28_beam_"+freqs[i] +".fits"
     make_rms, sensitivity[i], nside, "LB_v28_noise_256_"+freqs[i] +".fits", seed
     print, 'sequence number: ', freqs[i]
  end

end
