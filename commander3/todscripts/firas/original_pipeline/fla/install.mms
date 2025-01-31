!
!   Purpose: Install facility fla_lines_analysis.
!
!   Author: Gene Eplee, GSC, 7/1/93, SER 11413
!
!   Modifications: S. Alexander, HSTX, 4/94, SER 11703, Split line profiles.
!                  K. Jensen, HSTX, 10/26/94, SER 11935, Add code for
!                             line maps, dust, and temperature.
!                  K. Jensen, HSTX, 2/27/95, SPR 12117, Modifications for
!                             errors uncovered during validation.
!                  K. Jensen, HSTX, 3/3/95, SPR 12124, Procedures not loaded.

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg .tlb .exe .pro

fla : csdr$cld:fla.cld,-
      csdr$library:csdrmsg.olb(fla_msg),-
      csdr$help:csdrhelp.hlb(fla),-
      csdr$system:fla_gain.exe,-
      csdr$idl:fla.pro,-
      csdr$idl:fla_apod_rec.pro,-
      csdr$idl:fla_csq_high.pro,-
      csdr$idl:fla_csq_hres.pro,-
      csdr$idl:fla_csq_lowf.pro,-
      csdr$idl:fla_csq_lres.pro,-
      csdr$idl:fla_dipole.pro,-
      csdr$idl:fla_dst_st.pro,-
      csdr$idl:fla_fex_cvs_st.pro,-
      csdr$idl:fla_fip_csq.pro,-
      csdr$idl:fla_fip_lmh.pro,-
      csdr$idl:fla_fip_lml.pro,-
      csdr$idl:fla_fip_tcb.pro,-
      csdr$idl:fla_fir_line.pro,-
      csdr$idl:fla_fit.pro,-
      csdr$idl:fla_fitdust.pro,-
      csdr$idl:fla_fit_hi.pro,-
      csdr$idl:fla_fit_lhr.pro,-
      csdr$idl:fla_fit_llr.pro,-
      csdr$idl:fla_fms_sky_st.pro,-
      csdr$idl:fla_gal_tmpl.pro,-
      csdr$idl:fla_get_tcbr.pro,-
      csdr$idl:fla_init_fir.pro,-
      csdr$idl:fla_match.pro,-
      csdr$idl:fla_tcmbr.pro
 ! Facility fla has been installed.

csdr$cld:fla.cld : fla.cld
 copy fla.cld csdr$cld:fla.cld;0

csdr$system:fla_gain.exe : fla_gain.exe
 copy fla_gain.exe csdr$system:fla_gain.exe;0

csdr$idl:fla.pro : fla.pro
 copy fla.pro csdr$idl:fla.pro;0

csdr$idl:fla_apod_rec.pro : fla_apod_rec.pro
 copy fla_apod_rec.pro csdr$idl:fla_apod_rec.pro;0

csdr$idl:fla_csq_high.pro : fla_csq_high.pro
 copy fla_csq_high.pro csdr$idl:fla_csq_high.pro;0

csdr$idl:fla_csq_hres.pro : fla_csq_hres.pro
 copy fla_csq_hres.pro csdr$idl:fla_csq_hres.pro;0

csdr$idl:fla_csq_lowf.pro : fla_csq_lowf.pro
 copy fla_csq_lowf.pro csdr$idl:fla_csq_lowf.pro;0

csdr$idl:fla_csq_lres.pro : fla_csq_lres.pro
 copy fla_csq_lres.pro csdr$idl:fla_csq_lres.pro;0

csdr$idl:fla_dipole.pro : fla_dipole.pro
 copy fla_dipole.pro csdr$idl:fla_dipole.pro;0

csdr$idl:fla_dst_st.pro : fla_dst_st.pro
 copy fla_dst_st.pro csdr$idl:fla_dst_st.pro;0

csdr$idl:fla_fex_cvs_st.pro : fla_fex_cvs_st.pro
 copy fla_fex_cvs_st.pro csdr$idl:fla_fex_cvs_st.pro;0

csdr$idl:fla_fip_csq.pro : fla_fip_csq.pro
 copy fla_fip_csq.pro csdr$idl:fla_fip_csq.pro;0

csdr$idl:fla_fip_lmh.pro : fla_fip_lmh.pro
 copy fla_fip_lmh.pro csdr$idl:fla_fip_lmh.pro;0

csdr$idl:fla_fip_lml.pro : fla_fip_lml.pro
 copy fla_fip_lml.pro csdr$idl:fla_fip_lml.pro;0

csdr$idl:fla_fip_tcb.pro : fla_fip_tcb.pro
 copy fla_fip_tcb.pro csdr$idl:fla_fip_tcb.pro;0

csdr$idl:fla_fir_line.pro : fla_fir_line.pro
 copy fla_fir_line.pro csdr$idl:fla_fir_line.pro;0

csdr$idl:fla_fit.pro : fla_fit.pro
 copy fla_fit.pro csdr$idl:fla_fit.pro;0

csdr$idl:fla_fitdust.pro : fla_fitdust.pro
 copy fla_fitdust.pro csdr$idl:fla_fitdust.pro;0

csdr$idl:fla_fit_hi.pro : fla_fit_hi.pro
 copy fla_fit_hi.pro csdr$idl:fla_fit_hi.pro;0

csdr$idl:fla_fit_lhr.pro : fla_fit_lhr.pro
 copy fla_fit_lhr.pro csdr$idl:fla_fit_lhr.pro;0

csdr$idl:fla_fit_llr.pro : fla_fit_llr.pro
 copy fla_fit_llr.pro csdr$idl:fla_fit_llr.pro;0

csdr$idl:fla_fms_sky_st.pro : fla_fms_sky_st.pro
 copy fla_fms_sky_st.pro csdr$idl:fla_fms_sky_st.pro;0

csdr$idl:fla_gal_tmpl.pro : fla_gal_tmpl.pro
 copy fla_gal_tmpl.pro csdr$idl:fla_gal_tmpl.pro;0

csdr$idl:fla_get_tcbr.pro : fla_get_tcbr.pro
 copy fla_get_tcbr.pro csdr$idl:fla_get_tcbr.pro;0

csdr$idl:fla_init_fir.pro : fla_init_fir.pro
 copy fla_init_fir.pro csdr$idl:fla_init_fir.pro;0

csdr$idl:fla_match.pro : fla_match.pro
 copy fla_match.pro csdr$idl:fla_match.pro;0

csdr$idl:fla_tcmbr.pro : fla_tcmbr.pro
 copy fla_tcmbr.pro csdr$idl:fla_tcmbr.pro;0
