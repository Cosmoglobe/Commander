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
  sensitivity = [8.629, 4.771, 3.749, 2.316, 2.316, 1.907, 1.907, 1.637, 1.637, 1.126, 1.126, 0.782, 0.782, 0.822, 0.822, 0.846, 0.943, 0.943, 2.206, 2.617, 3.637, 7.262]*nside/512L
  fwhms       = [69.3, 56.8, 49.0, 41.6, 44.5, 36.9, 40.0, 33.0, 36.7, 30.2, 37.8, 26.3, 33.6, 23.7, 30.8, 28.9, 28.0, 28.6, 24.7, 22.5, 20.9, 17.9]

  seed = 736563
  for i = 0, n_elements(fwhms)-1 do begin
;     make_beam, fwhms [i], lmax, "LB_v28_beam_"+freqs[i] +".fits"
     make_rms, sensitivity[i], nside, "LB_v28_noise_256_"+freqs[i] +".fits", seed
     print, 'sequence number: ', freqs[i]
  end

end
