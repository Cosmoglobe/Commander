files = file_search('/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/*.fits')

amp_drain_current = [$
'DFK111B8DNI',  'DFK112B8DNI',  'DFK121B8DNI',  'DFK122B8DNI', $
'DRK111B8DNI',  'DRK112B8DNI',  'DRK121B8DNI',  'DRK122B8DNI', $
'DFKA111B8DNI', 'DFKA112B8DNI', 'DFKA121B8DNI', 'DFKA122B8DNI', $
'DRKA111B8DNI', 'DRKA112B8DNI', 'DRKA121B8DNI', 'DRKA122B8DNI', $
'DFQ111B8DNI',  'DFQ112B8DNI',  'DFQ121B8DNI',  'DFQ122B8DNI', $
'DRQ111B8DNI',  'DRQ112B8DNI',  'DRQ121B8DNI',  'DRQ122B8DNI', $
'DFQ211B8DNI',  'DFQ212B8DNI',  'DFQ221B8DNI',  'DFQ222B8DNI', $
'DRQ211B8DNI',  'DRQ212B8DNI',  'DRQ221B8DNI',  'DRQ222B8DNI', $
'DFV111B8DNI',  'DFV112B8DNI',  'DFV121B8DNI',  'DFV122B8DNI', $
'DRV111B8DNI',  'DRV112B8DNI',  'DRV121B8DNI',  'DRV122B8DNI', $
'DFV211B8DNI',  'DFV212B8DNI',  'DFV221B8DNI',  'DFV222B8DNI', $
'DRV211B8DNI',  'DRV212B8DNI',  'DRV221B8DNI',  'DRV222B8DNI', $
'DFW111B8DNI',  'DFW112B8DNI',  'DFW121B8DNI',  'DFW122B8DNI', $
'DRW111B8DNI',  'DRW112B8DNI',  'DRW121B8DNI',  'DRW122B8DNI', $
'DFW211B8DNI',  'DFW212B8DNI',  'DFW221B8DNI',  'DFW222B8DNI', $
'DRW211B8DNI',  'DRW212B8DNI',  'DRW221B8DNI',  'DRW222B8DNI', $
'DFW311B8DNI',  'DFW312B8DNI',  'DFW321B8DNI',  'DFW322B8DNI', $
'DRW311B8DNI',  'DRW312B8DNI',  'DRW321B8DNI',  'DRW322B8DNI', $
'DFW411B8DNI',  'DFW412B8DNI',  'DFW421B8DNI',  'DFW422B8DNI', $
'DRW411B8DNI',  'DRW412B8DNI',  'DRW421B8DNI',  'DRW422B8DNI']

radiometer_rf_bias = [$
'DRK113RFBI0',   'DRK114RFBI1',   'DRK123RFBI2',   'DRK124RFBI3', $
'DRKA113RFBI36', 'DRKA114RFBI37', 'DRKA123RFBI38', 'DRKA124RFBI39', $
'DRQ113RFBI20',  'DRQ114RFBI21',  'DRQ123RFBI22',  'DRQ124RFBI23', $
'DRQ213RFBI28',  'DRQ214RFBI29',  'DRQ223RFBI30',  'DRQ224RFBI31', $
'DRV113RFBI32',  'DRV114RFBI33',  'DRV123RFBI34',  'DRV124RFBI35', $
'DRV213RFBI12',  'DRV214RFBI13',  'DRV223RFBI14',  'DRV224RFBI15', $
'DRW113RFBI4',   'DRW114RFBI5',   'DRW123RFBI6',   'DRW124RFBI7', $
'DRW213RFBI24',  'DRW214RFBI25',  'DRW223RFBI26',  'DRW224RFBI27', $
'DRW313RFBI16',  'DRW314RFBI17',  'DRW323RFBI18',  'DRW324RFBI19', $
'DRW413RFBI8',   'DRW414RFBI9',   'DRW423RFBI10',  'DRW424RFBI11']

trs_temperatures = [$
'DTATOPPRIT', 'DTAMIDPRIT', 'DTATOPSECT', 'DTAMIDSECT', 'DTABOTSECT', $
'DTBTOPPRIT', 'DTBMIDPRIT', 'DTBTOPSECT', 'DTBMIDSECT', $
'DTAPXMIDRADT', 'DTBPXMIDRADT', 'DTAMXTOPRADT', 'DTBMXBOTRADT']

fpa_temps = [$
'DFK1AFEEDT', 'DFKA1BFEEDT', 'DFQ1AFEEDT', 'DFQ2BFEEDT', $
'DFW3AFEEDT', 'DFW3BFEEDT', 'DFV11FPATEET', 'DFV22FPATEET', $
'DFW11FPATEET', 'DFW22FPATEET', 'DFW32FPATEET', 'DFW3AOMTT', $
'DFW3BOMTT', 'DFK1BOMTT', 'DFKA1AOMTT', 'DFQ1BOMTT', 'DFQ2AOMTT']

rxb_temps = [$
'DRV111RXBAMPT', 'DRV222RXBAMPT', 'DRW111RXBAMPT', 'DRW221RXBAMPT',$
'DRW321RXBAMPT', 'DRK12RXBRIBT', 'DRKA12RXBRIBT', 'DRQ1RXBRIBT', $
'DRQ2RXBRIBT', 'DRW3RXBRIBT', 'DRPYPSHPRTKT', 'DRMYPSHPRTKT']

aeu_temps = [$
'DAW323_4AMPT', 'DAW2_14_23AMP_ADT', 'DAV113_4ADT', 'WDAW113_4ADT',$
'DAV223_4AMPT', 'DAQ113_4ADT', 'DAIHK1BDT', 'DAIHK2BDT', 'DACONVBDT']

pdu_temps = [$
'DPPINTT1', 'DPPINTT2', 'DPPINTT3',$
'DPV111_2FPAT', 'DPW221_2FPAT', 'DPW321_2FPAT',$
'DPV221_2RXBT', 'DPW111_2RXBT', 'DPW321_2RXBT']

aeu_voltages = [$
'DAP15VBD1', 'DAM15VBD1', 'DAP12VBD1', 'DAM12VBD1', 'DAP5VBD1',$
'DAP15VBD2', 'DAM15VBD2', 'DAP12VBD2', 'DAM12VBD2', 'DAP5VBD2',$
'DABD1V', $
'DABD2V']
;'DARREF1BD1', 'DARREF2BD1', $
;'DARREF1BD2', 'DARREF2BD2', $
;'DASPARE1']

pdu_voltages = [$
'DPFP7_2V', 'DPFM7_2V', 'DPRP7_2V', 'DPRM7_2V', $
'DFPLEDP10V', 'DPPHSWCONVP9V', 'DPPHSWCONVM9V',$
'DPFLDAP6_2V', 'DPFLDAM6_2V', 'DPFLDBP6_2V', 'DPFLDBM6_2V', $
'DPHKP15V', 'DPHKP5V', 'DPHKM15V']

fname_ampdrain = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/ampdrain.txt'
fname_rfbias   = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/rfbias.txt'
fname_trstemp  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/trstemp.txt'
fname_fpatemp  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/fpatemp.txt'
fname_rxbtemp  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/rxbtemp.txt'
fname_aeutemp  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/aeutemp.txt'
fname_pdutemp  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/pdutemp.txt'
fname_aeuvolt  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/aeuvolt.txt'
fname_pduvolt  = '/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/housekeeping/pduvolt.txt'


openw, 1, fname_ampdrain, width=1100
openw, 2, fname_rfbias, width=800
openw, 3, fname_trstemp, width=200
openw, 4, fname_fpatemp, width=300
openw, 5, fname_rxbtemp, width=200
openw, 6, fname_aeutemp, width=200
openw, 7, fname_pdutemp, width=200
openw, 8, fname_aeuvolt, width=200
openw, 9, fname_pduvolt, width=200
for i=0, n_elements(files)-1 do begin
    file = files(i)
    print, file, i+1, n_elements(files)
    fits_read_tod, file, arch
    T_arr = fltarr(n_elements(amp_drain_current))
    for n=0,n_elements(amp_drain_current)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, amp_drain_current(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 1, T_arr

    T_arr = fltarr(n_elements(radiometer_rf_bias))
    for n=0,n_elements(radiometer_rf_bias)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, radiometer_rf_bias(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 2, T_arr

    T_arr = fltarr(n_elements(trs_temperatures))
    for n=0,n_elements(trs_temperatures)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, trs_temperatures(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 3, T_arr

    T_arr = fltarr(n_elements(fpa_temps))
    for n=0,n_elements(fpa_temps)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, fpa_temps(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 4, T_arr

    T_arr = fltarr(n_elements(rxb_temps))
    for n=0,n_elements(rxb_temps)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, rxb_temps(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 5, T_arr

    T_arr = fltarr(n_elements(aeu_temps))
    for n=0,n_elements(aeu_temps)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, aeu_temps(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 6, T_arr

    T_arr = fltarr(n_elements(pdu_temps))
    for n=0,n_elements(pdu_temps)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, pdu_temps(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 7, T_arr

    T_arr = fltarr(n_elements(aeu_voltages))
    for n=0,n_elements(aeu_voltages)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, aeu_voltages(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 8, T_arr

    T_arr = fltarr(n_elements(pdu_voltages))
    for n=0,n_elements(pdu_voltages)-1 do begin
        T = AIHK_Arch2Mnemonic(arch, pdu_voltages(n), 0)
        T_arr(n) = T(n_elements(T)/2)
    endfor
    printf, 9, T_arr
endfor
close, 1
close, 2
close, 3
close, 4
close, 5
close, 6
close, 7
close, 8
close, 9


end
