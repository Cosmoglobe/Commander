'''Definition of data types for structured arrays
'''

import numpy as np

dt_adt = '<i8' #64-byte quadword. number of 100ns intervals since 1858-11-17
dt_word = '<i2'
dt_uword = '<u2'
dt_byte = '<i1'
dt_long = '<i4'
dt_float = '<f4'
dt_double = '<f8'
dt_complex = '<c8'
dt_doublecomplex = '<c16'

#names_cthead = ['gmt', 'time', 'space_time', 'mjr_frm_no', 'orbit', 'hskp1_tlm_fmt', 'hskp2_tlm_fmt',
#                'ingest_spares', 'dataset_id', 'instr_spares']

#fmt_cthead = ['<S14', '<i8', ('<i1', (6,)), '<i4', '<i4', '<i1', '<i1', ('<i1', (18,)), '<i2', ('<i1', (6,))]

#dt_ct_head = np.dtype({'names' : names_cthead, 'formats' : fmt_cthead})

dt_ct_head = np.dtype([('gmt', '<S14'),
                       ('time', dt_adt),
                       ('space_time', dt_byte, (6,)),
                       ('mjr_frm_no', dt_long),
                       ('orbit', dt_long),
                       ('hskp1_tlm_fmt', dt_byte),
                       ('hskp2_tlm_fmt', dt_byte),
                       ('ingest_spares', dt_byte, (18,)),
                       ('dataset_id', dt_word),
                       ('instr_spares', dt_byte, (6,))])

dt_sci_head = np.dtype([('chan_id', '<i2'),
                        ('gain', '<i2'),
                        ('mtm_speed', '<i2'),
                        ('mtm_length', '<i2'),
                        ('data_qual', '<i1', (60,)),
                        ('data_ready', '<i2', (8,)),
                        ('sc_head0', '<i2'),
                        ('sc_head1a', '<i1'),
                        ('sc_head1b', '<i1'),
                        ('sc_head2', '<i2'),
                        ('sc_head3', '<i2'),
                        ('sc_head4', '<i2'),
                        ('sc_head5', '<i2'),
                        ('sc_head6', '<i2'),
                        ('sc_head7', '<i2'),
                        ('sc_head8', '<i2'),
                        ('sc_head9', '<i2'),
                        ('sc_head10', '<i2'),
                        ('sc_head11', '<i2'),
                        ('sc_head12', '<i2'),
                        ('sc_head13', '<i2'),
                        ('sc_head14', '<i2'),
                        ('sc_head15', '<i2'),
                        ('sc_head16', '<i2'),
                        ('sc_head17', '<i2'),
                        ('sc_head18', '<i2'),
                        ('sc_head19', '<i2'),
                        ('sc_head20', '<i2'),
                        ('sc_head21', '<i2'),
                        ('sc_head22', '<i2'),
                        ('sc_head23', '<i2'),
                        ('sc_head24', '<i2'),
                        ('sc_head25', '<S2'),])

dt_ifg_data = np.dtype([('ifg', '<i2', (512,)),
                        ('gltch', '<u2', (32,))])

dt_dq_data = np.dtype([('fake', '<i2'),
                       ('xcal_pos', '<i2'),
                       ('iref_temp', '<f4'),
                       ('eng_time', dt_adt),
                       ('eng_rec', '<i4'),
                       ('data_quality', '<i1', (110,)),
                       ('ifg_no', '<i4'),
                       ('dq_spares', '<i1', (24,))])

dt_collect_time = np.dtype([('midpoint_time', '<i8'),
                            ('badtime_flag', '<i1'),
                            ('fpp_spare', '<i1')])

dt_attitude = np.dtype([('pixel_no', '<i4'),
                        ('equatorial', '<f4', (3,)),
                        ('ra', '<i2'),
                        ('dec', '<i2'),
                        ('terr_pixel_no', '<i4'),
                        ('terr_latitude', '<i2'),
                        ('terr_longitude', '<i2'),
                        ('earth_limb', '<i2'),
                        ('earth_limb_azimuth', '<i2'),
                        ('sun_angle', '<i2'),
                        ('moon_angle', '<i2'),
                        ('moon_az_angle', '<i2'),
                        ('moon_phase', '<i2'),
                        ('sun_moon_dist', '<f4'),
                        ('cobe_moon_dist', '<f4'),
                        ('altitude', '<i2'),
                        ('projected_barycentric_velocity', '<i2'),
                        ('mcilwain_l_param', '<i2'),
                        ('galactic_longitude', '<i2'),
                        ('galactic_latitude', '<i2'),
                        ('ecliptic_longitude', '<i2'),
                        ('ecliptic_latitude', '<i2'),
                        ('orbital_phase', '<i2'),
                        ('projected_geocentric_velocity', '<i2'),
                        ('scan_angle', '<i2'),
                        ('sc_rotation_angle', '<i2'),
                        ('solution', '<i1'),
                        ('pixel_definition', '<S1'),
                        ('skymap_index', '<i2'),
                        ('exc_galactic_lat', '<i2'),
                        ('terr_rad_byte', '<i1'),
                        ('outside_galaxy_cut', '<i1'),
                        ('att_spares', '<i1', (2,))])


dt_fdq_sdf = np.dtype([('ct_head', dt_ct_head),
                       ('sci_head', dt_sci_head),
                       ('ifg_data', dt_ifg_data),
                       ('dq_data', dt_dq_data),
                       ('collect_time', dt_collect_time),
                       ('attitude', dt_attitude)])


dt_sci_time = np.dtype([('bin_time', dt_adt)])
#dt_sci_time = dt_adt

dt_en_head = np.dtype([('limits', '<S23'),
                       ('convert', '<S23'),
                       ('a_spare', '<i4'),
                       ('key_id', dt_word),
                       ('number_of_records', dt_long),
                       ('enghead_spares', dt_byte, (8,)),
                       ('sci_time', dt_sci_time, (4,)),
                       ('first_eng_time', dt_adt),
                       ('last_eng_time', dt_adt),
                       ('dq_summary_flag', dt_byte, (4,)),
                       ('ifg_no', dt_long, (4,)),
                       ('att_summary_flag', dt_byte, (4,)),
                       ('enghead_spares2', dt_byte, (34,))])

#Must give names, formats, offsets separately if I want offsets so I can have 
#multiple names referring to the same area in memory (i.e. unions)
names_en_stat = ['group1',
                 'stat_word_1',
                 'int_ref_temp_a',
                 'ref_hrn_temp_a',
                 'stat_word_4',
                 'stat_word_5',
                 'sky_hrn_temp_a',
                 'ext_cal_temp_a',
                 'stat_word_8',
                 'stat_word_9',
                 'int_ref_temp_b',
                 'ref_hrn_temp_b',
                 'stat_word_12',
                 'stat_word_13',
                 'sky_hrn_temp_b',
                 'ext_cal_temp_b',
                 'stat_word_16',
                 'group2',
                 'grt_addr',
                 'micro_stat_bus',
                 'bol_cmd_bias',
                 'dwell_stat',
                 'lvdt_stat',
                 'hot_spot_cmd',
                 'engstat_spares',
                 'power_a_status',
                 'power_b_status',
                 'engstat_spares2']

fmts_en_stat = [(dt_uword, (16,)),
                dt_uword, dt_uword, dt_uword, dt_uword, dt_uword, dt_uword, dt_uword, dt_uword,
                dt_uword, dt_uword, dt_uword, dt_uword, dt_uword, dt_uword, dt_uword, dt_uword,
                (dt_byte, (14,)),
                (dt_byte, (2,)), (dt_byte, (4,)), (dt_byte, (4,)), (dt_byte, (2,)), (dt_byte, (2,)),
                (dt_byte, (2,)),
                (dt_byte, (10,)),
                (dt_byte, (2,)),
                (dt_byte, (2,)),
                (dt_byte, (6,))]

off_en_stat = [0,
               0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
               32,
               32, 34, 38, 42, 44,
               46, 48, 58, 60, 62]

dt_en_stat = np.dtype({'names' : names_en_stat,
                       'formats' : fmts_en_stat,
                       'offsets' : off_en_stat})

#Engineering status type without the unions
dt_en_stat2 = np.dtype([('group1', dt_uword, (16,)),
                        ('group2', dt_byte, (14,)),
                        ('hot_spot_cmd', dt_byte, (2,)),
                        ('engstat_spares', dt_byte, (10,)),
                        ('power_a_status', dt_byte, (2,)),
                        ('power_b_status', dt_byte, (2,)),
                        ('engstat_spares2', dt_byte, (6,))])

names_en_analog = ['grt',
                   'a_lo_grt', 'a_hi_grt', 'b_lo_grt', 'b_hi_grt',
                   'a_lo_xcal_tip', 'a_lo_skyhorn', 'a_lo_refhorn', 'a_lo_ical', 'a_lo_dihedral',
                   'a_lo_bol_assem', 'a_lo_mirror', 'a_lo_cal_resistors', 'a_lo_xcal_cone', 'a_lo_collimator',
                   'a_hi_xcal_tip', 'a_hi_skyhorn', 'a_hi_refhorn', 'a_hi_ical', 'a_hi_dihedral',
                   'a_hi_bol_assem', 'a_hi_mirror', 'a_hi_cal_resistors', 'a_hi_xcal_cone', 'a_hi_collimator',
                   'b_lo_xcal_tip', 'b_lo_skyhorn', 'b_lo_refhorn', 'b_lo_ical', 'b_lo_dihedral',
                   'b_lo_bol_assem', 'b_lo_mirror', 'b_lo_cal_resistors', 'b_lo_xcal_cone', 'b_lo_collimator',
                   'b_hi_xcal_tip', 'b_hi_skyhorn', 'b_hi_refhorn', 'b_hi_ical', 'b_hi_dihedral',
                   'b_hi_bol_assem', 'b_hi_mirror', 'b_hi_cal_resistors', 'b_hi_xcal_cone', 'b_hi_collimator',
                   'group1',
                   'temp_ctl', 'ipdu_temp', 'cna_temp', 'dbx_tmp', 'stat_mon_temp',
                   'pamp_chan', 'pamp_op', 'hot_spot', 'mtm_cal_mtr', 'mtm_pos', 'bol_volt',
                   'ipdu_bolt', 'ipdu_amp']

fmts_en_analog = [(dt_float, (64,)),
                  (dt_float, (16,)), (dt_float, (16,)), (dt_float, (16,)), (dt_float, (16,)),
                  dt_float, dt_float, dt_float, dt_float, dt_float,
                  (dt_float, (4,)), dt_float, (dt_float, (4,)), dt_float, dt_float,
                  dt_float, dt_float, dt_float, dt_float, dt_float,
                  (dt_float, (4,)), dt_float, (dt_float, (4,)), dt_float, dt_float,
                  dt_float, dt_float, dt_float, dt_float, dt_float,
                  (dt_float, (4,)), dt_float, (dt_float, (4,)), dt_float, dt_float,
                  dt_float, dt_float, dt_float, dt_float, dt_float,
                  (dt_float, (4,)), dt_float, (dt_float, (4,)), dt_float, dt_float,
                  (dt_float, (62,)),
                  (dt_float, (8,)), (dt_float, (2,)), (dt_float, (4, )), (dt_float, (2,)), (dt_float, (2)),
                  dt_float, dt_float, (dt_float, (2,)), (dt_float, (2,)), (dt_float, (2,)), (dt_float, (4,)),
                  (dt_float, (20,)), (dt_float, (12,))]

off_en_analog = [0,
                 0, 64, 128, 192,
                 0, 4, 8, 12, 16,
                 20, 36, 40, 56, 60,
                 64, 68, 72, 76, 80,
                 84, 100, 104, 120, 124,
                 128, 132, 136, 140, 144,
                 148, 164, 168, 184, 188,
                 192, 196, 200, 204, 208,
                 212, 228, 232, 248, 252,
                 256,
                 256, 288, 296, 312, 320,
                 328, 332, 336, 344, 352, 360,
                 376, 456]

dt_en_analog = np.dtype({'names' : names_en_analog,
                         'formats' : fmts_en_analog,
                         'offsets' : off_en_analog})

dt_en_analog2 = np.dtype([('grt', dt_float, (64,)),
                          ('group1', dt_float, (62,))])

dt_en_xcal = np.dtype([('pos', dt_word, (2,)),
                       ('xcal_spares', dt_byte, (50,))])

dt_chan = np.dtype([('up_sci_mode', dt_byte),
                    ('fakeit', dt_byte),
                    ('up_adds_per_group', dt_byte),
                    ('up_swps_per_ifg', dt_byte),
                    ('xmit_mtm_speed', dt_byte),
                    ('xmit_mtm_len', dt_byte),
                    ('sci_gain', dt_word),
                    ('dither', dt_byte),
                    ('setup_spares', dt_byte, (15,))])

dt_en_tempdiff = np.dtype([('xcal', dt_word),
                           ('ical', dt_word),
                           ('skyhorn', dt_word),
                           ('refhorn', dt_word),
                           ('dihedral', dt_word),
                           ('collimator_mirror', dt_word),
                           ('bol_assem', dt_word, (4,))])

dt_en_tail = np.dtype([('engtail_spares', dt_byte, (8,)),
                       ('hskp_flag', dt_byte),
                       ('lmac_analog_temp', dt_float),
                       ('lmac_digital_temp', dt_float),
                       ('tlm_qual_maj_frm', dt_byte, (2,)),
                       ('eng_spares', dt_byte, (9,))])

dt_fdq_eng = np.dtype([('ct_head', dt_ct_head),
                       ('en_head', dt_en_head),
                       ('en_stat', dt_en_stat),
                       ('en_analog', dt_en_analog),
                       ('en_xcal', dt_en_xcal),
                       ('chan', dt_chan, (4,)),
                       ('en_tempdiff', dt_en_tempdiff, (2,)),
                       ('en_tail', dt_en_tail)])

dt_fdq_eng2 = np.dtype([('ct_head', dt_ct_head),
                        ('en_head', dt_en_head),
                        ('en_stat', dt_en_stat2),
                        ('en_analog', dt_en_analog2),
                        ('en_xcal', dt_en_xcal),
                        ('chan', dt_chan, (4,)),
                        ('en_tempdiff', dt_en_tempdiff, (2,)),
                        ('en_tail', dt_en_tail)])

dt_short_science = np.dtype([('time', dt_adt),
                             ('channel_id', dt_byte),
                             ('pixel_no', dt_long),
                             ('mtm_scan_speed', dt_byte),
                             ('mtm_scan_length', dt_byte),
                             ('sci_mode', dt_byte),
                             ('adds_per_group', dt_byte),
                             ('transition_flag', dt_byte),
                             ('pixel_definition', '<S1'),
                             ('skymap_index', dt_word),
                             ('data_quality', dt_byte),
                             ('attitude_quality', dt_byte),
                             ('xcal_temp', dt_float),
                             ('skyhorn_temp', dt_float),
                             ('refhorn_temp', dt_float),
                             ('ical_temp', dt_float),
                             ('dihedral_temp', dt_float),
                             ('bolometer_temp', dt_float),
                             ('bol_cmd_bias', dt_byte),
                             ('exc_galactic_lat', dt_word),
                             ('sun_angle', dt_word),
                             ('moon_angle', dt_word),
                             ('earth_limb', dt_word),
                             ('galactic_latitude', dt_word),
                             ('galactic_longitude', dt_word),
                             ('glitch_count', dt_word), 
                             ('original_pixel', dt_word)]) 

#No saved files, but reproducing con_check.rdl (fic_scc.rdl)
dt_fic_scc = np.dtype([('gmt', '<S14'),
                       ('time', dt_adt),
                       ('coadd_gmt', '<S14'),
                       ('coadd_time', dt_adt),
                       ('chan_id', dt_byte),
                       ('mtm_speed', dt_byte),
                       ('mtm_length', dt_byte),
                       ('fakeit', dt_byte),
                       ('sci_mode', dt_byte),
                       ('adds_per_group', dt_byte),
                       ('xcal_pos', dt_byte),
                       ('sweeps', dt_byte),
                       ('gain', dt_float),
                       ('glitch_rate', dt_float),
                       ('pixel_no', dt_long),
                       ('earth_limb', dt_word),
                       ('earth_limb_azimuth', dt_word),
                       ('sun_angle', dt_word),
                       ('moon_angle', dt_word),
                       ('moon_az_angle', dt_word),
                       ('moon_phase', dt_word),
                       ('sun_moon_dist', dt_float),
                       ('cobe_moon_dist', dt_float),
                       ('b', dt_float),
                       ('c', dt_float),
                       ('glitch_pos', dt_word),
                       ('glitch_amp', dt_float),
                       ('glitch_iter', dt_word),
                       ('glitch_signal', dt_float),
                       ('noise', dt_float),
                       ('con_check', dt_long),
                       ('spares', dt_byte, (16,))])

dt_fil_scc = np.dtype([('gmt', '<S14'),
                       ('time', dt_adt),
                       ('coadd_gmt', '<S14'),
                       ('coadd_time', dt_adt),
                       ('chan_id', dt_byte),
                       ('mtm_speed', dt_byte),
                       ('mtm_length', dt_byte),
                       ('fakeit', dt_byte),
                       ('sci_mode', dt_byte),
                       ('adds_per_group', dt_byte),
                       ('xcal_pos', dt_byte),
                       ('sweeps', dt_byte),
                       ('gain', dt_float),
                       ('glitch_rate', dt_float),
                       ('pixel_no', dt_long),
                       ('earth_limb', dt_word),
                       ('earth_limb_azimuth', dt_word),
                       ('sun_angle', dt_word),
                       ('moon_angle', dt_word),
                       ('moon_az_angle', dt_word),
                       ('moon_phase', dt_word),
                       ('sun_moon_dist', dt_float),
                       ('cobe_moon_dist', dt_float),
                       ('b', dt_float),
                       ('c', dt_float),
                       ('pre_dg_median', dt_float),
                       ('pre_dg_noise', dt_float),
                       ('glitch_pos', dt_word),
                       ('glitch_amp', dt_float),
                       ('glitch_iter', dt_word),
                       ('glitch_signal', dt_float),
                       ('post_dg_median', dt_float),
                       ('noise', dt_float),
                       ('con_check', dt_long),
                       ('spares', dt_byte, (4,))])

dt_con_check = dt_fil_scc

#Fortran code has primary and secondary templates as separate structures, but I am 
#combining them. Primary template will just ignore the 3 extra variables

#12 bytes
dt_ifgs = np.dtype([('time', dt_adt),
                    ('pixel_no', dt_long)])

#1205 bytes
dt_template = np.dtype([('num_ifgs', dt_word),
                        ('neighbors', dt_byte),
                        ('neighbor_num_ifgs', dt_word),
                        ('ifgs', dt_ifgs, (100,))])

#21 bytes total
#LOOK AT BOOLEAN FOR SUBTRACTED
dt_pstemplates = np.dtype([('subtracted', dt_byte),
                           ('amplitude', dt_float),
                           ('snr', dt_float),
                           ('variance', dt_float),
                           ('b_average', dt_float),
                           ('b_variance', dt_float)])

#12 bytes
dt_transient = np.dtype([('c_average', dt_float),
                         ('c_variance', dt_float),
                         ('bl_trans_coeff', dt_float)])

#6 bytes
dt_deglitch = np.dtype([('glitch_iter', dt_word),
                        ('glitch_signal', dt_float)])

names_coad_spec_data = ['chan_id', 'mtm_speed', 'mtm_length', 'fakeit', 'sci_mode',
                        'adds_per_group', 'xcal_pos', 'sweeps', 'gain_sum', 'glitch_rate',
                        'dq_summary_flag', 'att_summary_flag', 'peak_pos',
                        'nyquist_hertz', 'nyquist_icm',
                        'prim_template', 'sec_template', 'transient', 'deglitch',
                        'noise', 'bin_info', 'bl_coeffs',
                        'bol_cmd_bias', 'bol_volt', 'bol_volt_min', 'bol_volt_max',
                        'temp',
                        'xcal', 'ical', 'skyhorn', 'refhorn', 'dihedral', 'collimator_mirror', 'bol_assem',
                        'temp_sigma',
                        'xcal_sigma', 'ical_sigma', 'skyhorn_sigma', 'refhorn_sigma',
                        'dihedral_sigma', 'collimator_mirror_sigma', 'bol_assem_sigma',
                        'temp_min',
                        'xcal_min', 'ical_min', 'skyhorn_min', 'refhorn_min',
                        'dihedral_min', 'collimator_mirror_min', 'bol_assem_min',
                        'temp_max',
                        'xcal_max', 'ical_max', 'skyhorn_max', 'refhorn_max',
                        'dihedral_max', 'collimator_mirror_max', 'bol_assem_max',
                        'spares']

fmts_coad_spec_data = [dt_byte, dt_byte, dt_byte, dt_byte, dt_byte,
                       dt_byte, dt_byte, dt_byte, dt_long, dt_float,
                       dt_byte, dt_byte, dt_word,
                       dt_float, dt_float,
                       dt_pstemplates, dt_pstemplates, dt_transient, dt_deglitch,
                       dt_float, (dt_word, (256,)), (dt_float, (5,)),
                       dt_word, dt_float, (dt_float, (4,)), (dt_float, (4,)),
                       (dt_float, (10,)),
                       dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                       (dt_float, (10,)),
                       dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                       (dt_float, (10,)),
                       dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                       (dt_float, (10,)),
                       dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                       (dt_byte, (14,))]

off_coad_spec_data = [0, 1, 2, 3, 4,
                      5, 6, 7, 8, 12,
                      16, 17, 18,
                      20, 24,
                      28, 49, 70, 82, #dt_template ...
                      88, 92, 604, 
                      624, 626, 630, 646,
                      662,
                      662, 666, 670, 674, 678, 682, 686,
                      702,
                      702, 706, 710, 714, 718, 722, 726,
                      742,
                      742, 746, 750, 754, 758, 762, 766,
                      782,
                      782, 786, 790, 794, 798, 802, 806,
                      822]

names_coad_spec_datal = ['chan_id', 'mtm_speed', 'mtm_length', 'fakeit', 'sci_mode',
                         'adds_per_group', 'xcal_pos', 'sweeps', 'gain_sum', 'glitch_rate',
                         'dq_summary_flag', 'att_summary_flag', 'peak_pos',
                         'nyquist_hertz', 'nyquist_icm', 'orphans',
                         'template', 'prim_template', 'sec_template', 'transient', 'deglitch',
                         'noise', 'bin_info', 'bl_coeffs',
                         'bol_cmd_bias', 'bol_volt', 'bol_volt_min', 'bol_volt_max',
                         'temp',
                         'xcal', 'ical', 'skyhorn', 'refhorn', 'dihedral', 'collimator_mirror', 'bol_assem',
                         'temp_sigma',
                         'xcal_sigma', 'ical_sigma', 'skyhorn_sigma', 'refhorn_sigma',
                         'dihedral_sigma', 'collimator_mirror_sigma', 'bol_assem_sigma',
                         'temp_min',
                         'xcal_min', 'ical_min', 'skyhorn_min', 'refhorn_min',
                         'dihedral_min', 'collimator_mirror_min', 'bol_assem_min',
                         'temp_max',
                         'xcal_max', 'ical_max', 'skyhorn_max', 'refhorn_max',
                         'dihedral_max', 'collimator_mirror_max', 'bol_assem_max',
                         'spares']

fmts_coad_spec_datal = [dt_byte, dt_byte, dt_byte, dt_byte, dt_byte,
                        dt_byte, dt_byte, dt_byte, dt_long, dt_float,
                        dt_byte, dt_byte, dt_word,
                        dt_float, dt_float, dt_byte,
                        dt_template, dt_pstemplates, dt_pstemplates, dt_transient, dt_deglitch,
                        dt_float, (dt_word, (256,)), (dt_float, (5,)),
                        dt_word, dt_float, (dt_float, (4,)), (dt_float, (4,)),
                        (dt_float, (10,)),
                        dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                        (dt_float, (10,)),
                        dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                        (dt_float, (10,)),
                        dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                        (dt_float, (10,)),
                        dt_float, dt_float, dt_float, dt_float, dt_float, dt_float, (dt_float, (4,)),
                        (dt_byte, (44,))]

off_coad_spec_datal = [0, 1, 2, 3, 4,
                       5, 6, 7, 8, 12,
                       16, 17, 18,
                       20, 24, 28,
                       29, 1234, 1255, 1276, 1288,
                       1294, 1298, 1810,
                       1830, 1832, 1836, 1852,
                       1868,
                       1868, 1872, 1876, 1880, 1884, 1888, 1892,
                       1908,
                       1908, 1912, 1916, 1920, 1924, 1928, 1932,
                       1948,
                       1948, 1952, 1956, 1960, 1964, 1968, 1972,
                       1988,
                       1988, 1992, 1996, 2000, 2004, 2008, 2012,
                       2028]

#836 bytes
dt_coad_spec_data = np.dtype({'names' : names_coad_spec_data,
                              'formats' : fmts_coad_spec_data,
                              'offsets' : off_coad_spec_data})

#2042 bytes
dt_coad_spec_datal = np.dtype({'names' : names_coad_spec_datal,
                               'formats' : fmts_coad_spec_datal,
                               'offsets' : off_coad_spec_datal})

dt_coad_spec_datal2 = np.dtype([('chan_id', dt_byte),
                                ('mtm_speed', dt_byte),
                                ('mtm_length', dt_byte),
                                ('fakeit', dt_byte),
                                ('sci_mode', dt_byte),
                                ('adds_per_group', dt_byte),
                                ('xcal_pos', dt_byte),
                                ('sweeps', dt_byte),
                                ('gain_sum', dt_long),
                                ('glitch_rate', dt_float),
                                ('dq_summary_flag', dt_byte),
                                ('att_summary_flag', dt_byte),
                                ('peak_pos', dt_word),
                                ('nyquist_hertz', dt_float),
                                ('nyquist_icm', dt_float),
                                ('orphans', dt_byte),
                                ('template', dt_template),
                                ('prim_template', dt_pstemplates),
                                ('sec_template', dt_pstemplates),
                                ('transient', dt_transient),
                                ('deglitch', dt_deglitch),
                                ('noise', dt_float),
                                ('bin_info', dt_word, (256, )),
                                ('bl_coeffs', dt_float, (5, )),
                                ('bol_cmd_bias', dt_word),
                                ('bol_volt', dt_float),
                                ('bol_volt_min', dt_float, (4, )),
                                ('bol_volt_max', dt_float, (4, )),
                                ('temp', dt_float, (10, )),
                                ('temp_sigma', dt_float, (10, )),
                                ('temp_min', dt_float, (10, )),
                                ('temp_max', dt_float, (10, )),
                                ('spares', dt_byte, (44, ))])

dt_coad_spec_head = np.dtype([('first_gmt', '<S14'),
                              ('first_time', dt_adt),
                              ('first_space_time', dt_byte, (6,)),
                              ('first_mjr_frm_no', dt_long),
                              ('last_gmt', '<S14'),
                              ('last_time', dt_adt),
                              ('last_space_time', dt_byte, (6,)),
                              ('last_mjr_frm_no', dt_long),
                              ('coadd_no', dt_long),
                              ('num_ifgs', dt_word),
                              ('label', '<S60'),
                              ('comb_num_ifgs', dt_float),
                              ('real_comb_chi_square', dt_float),
                              ('imag_comb_chi_square', dt_float),
                              ('spares', dt_byte, (4,))])

dt_coad_spec_headl = np.dtype([('first_gmt', '<S14'),
                               ('first_time', dt_adt),
                               ('first_space_time', dt_byte, (6,)),
                               ('first_mjr_frm_no', dt_long),
                               ('last_gmt', '<S14'),
                               ('last_time', dt_adt),
                               ('last_space_time', dt_byte, (6,)),
                               ('last_mjr_frm_no', dt_long),
                               ('coadd_no', dt_long),
                               ('num_ifgs', dt_word),
                               ('adj_num_ifgs', dt_float),
                               ('times', dt_word, (100,)),
                               ('label', '<S60'),
                               ('comb_num_ifgs', dt_float),
                               ('adj_comb_num_ifgs', dt_float),
                               ('real_comb_chi_square', dt_float),
                               ('imag_comb_chi_square', dt_float),
                               ('spares', dt_byte, (24,))])

dt_coad_data = np.dtype([('ifg', dt_float, (512,)),
                         ('real_var', dt_float, (257,)),
                         ('imag_var', dt_float, (257,)),
                         ('real_imag_var', dt_float, (257,)),
                         ('spares', dt_byte, (318,))])

dt_coad_datal = np.dtype([('ifg', dt_float, (512,)),
                          ('fft_length', dt_long),
                          ('real_var', dt_float, (361,)),
                          ('imag_var', dt_float, (361,)),
                          ('real_imag_var', dt_float, (361,)),
                          ('spares', dt_byte, (1186,))])

names_en_sigma = ['sig_grt',
                   'sig_a_lo_grt', 'sig_a_hi_grt', 'sig_b_lo_grt', 'sig_b_hi_grt',
                   'sig_a_lo_xcal_tip', 'sig_a_lo_skyhorn', 'sig_a_lo_refhorn', 'sig_a_lo_ical',
                   'sig_a_lo_dihedral', 'sig_a_lo_bol_assem', 'sig_a_lo_mirror', 'sig_a_lo_cal_resistors',
                   'sig_a_lo_xcal_cone', 'sig_a_lo_collimator',
                   'sig_a_hi_xcal_tip', 'sig_a_hi_skyhorn', 'sig_a_hi_refhorn', 'sig_a_hi_ical',
                   'sig_a_hi_dihedral', 'sig_a_hi_bol_assem', 'sig_a_hi_mirror', 'sig_a_hi_cal_resistors',
                   'sig_a_hi_xcal_cone', 'sig_a_hi_collimator',
                   'sig_b_lo_xcal_tip', 'sig_b_lo_skyhorn', 'sig_b_lo_refhorn', 'sig_b_lo_ical',
                   'sig_b_lo_dihedral', 'sig_b_lo_bol_assem', 'sig_b_lo_mirror', 'sig_b_lo_cal_resistors',
                   'sig_b_lo_xcal_cone', 'sig_b_lo_collimator',
                   'sig_b_hi_xcal_tip', 'sig_b_hi_skyhorn', 'sig_b_hi_refhorn', 'sig_b_hi_ical',
                   'sig_b_hi_dihedral', 'sig_b_hi_bol_assem', 'sig_b_hi_mirror', 'sig_b_hi_cal_resistors',
                   'sig_b_hi_xcal_cone', 'sig_b_hi_collimator',
                   'group2',
                   'sig_temp_ctl', 'sig_ipdu_temp', 'sig_cna_temp', 'sig_dbx_tmp', 'sig_stat_mon_temp',
                   'sig_pamp_chan', 'sig_pamp_op', 'sig_hot_spot', 'sig_mtm_cal_mtr', 'sig_mtm_pos',
                   'sig_bol_volt', 'sig_ipdu_bolt', 'sig_ipdu_amp']

#Since en_sigma is just the sigmas of en_analog, we can use the formats and offsets of en_analog
dt_en_sigma = np.dtype({'names' : names_en_sigma,
                        'formats' : fmts_en_analog,
                        'offsets' : off_en_analog})

dt_en_sigma2 = np.dtype([('sig_grt', dt_float, (64,)),
                         ('group2', dt_float, (62,))])

dt_fic_sky = np.dtype([('ct_head', dt_ct_head),
                       ('coad_spec_head', dt_coad_spec_head),
                       ('coad_spec_data', dt_coad_spec_data),
                       ('coad_data', dt_coad_data),
                       ('en_stat', dt_en_stat),
                       ('en_analog', dt_en_analog),
                       ('en_sigma', dt_en_sigma),
                       ('en_tempdiff', dt_en_tempdiff, (2,)),
                       ('attitude', dt_attitude)])

dt_fil_sky = np.dtype([('ct_head', dt_ct_head),
                       ('coad_spec_head', dt_coad_spec_headl),
                       ('coad_spec_data', dt_coad_spec_datal),
                       ('coad_data', dt_coad_datal),
                       ('en_stat', dt_en_stat),
                       ('en_analog', dt_en_analog),
                       ('en_sigma', dt_en_sigma),
                       ('en_tempdiff', dt_en_tempdiff, (2,)),
                       ('attitude', dt_attitude)])

#fil_sky data type that can be written to a HDF5 file
dt_fil_sky2 = np.dtype([('ct_head', dt_ct_head),
                        ('coad_spec_head', dt_coad_spec_headl),
                        ('coad_spec_data', dt_coad_spec_datal2),
                        ('coad_data', dt_coad_datal),
                        ('en_stat', dt_en_stat2),
                        ('en_analog', dt_en_analog2),
                        ('en_sigma', dt_en_sigma2),
                        ('en_tempdiff', dt_en_tempdiff, (2,)),
                        ('attitude', dt_attitude)])

dt_coa_rec = dt_fil_sky

dt_gltchpro = np.dtype([('gltchpro', dt_float, (512,))])
dt_gtrf = np.dtype([('trans', dt_float, (128,))])
dt_basis = np.dtype([('leg_poly', dt_double, (5, 512))])

dt_ident = np.dtype([('chan_id', dt_byte),
                     ('mtm_speed', dt_byte),
                     ('mtm_length', dt_byte),
                     ('sci_mode', dt_byte),
                     ('adds_per_group', dt_byte)])

#This does not use unions that the Fortran code uses, but we are not strapped for memory so we
#don't need to do this
dt_fil_cov = np.dtype([('ct_head', dt_ct_head),
               ('ident', dt_ident),
               ('num_ifgs', dt_long),
               ('adj_num_ifgs', dt_float),
               ('num_coadds', dt_long),
               ('avg_rvoltspec', dt_double, (361,)),
               ('avg_ivoltspec', dt_double, (361,)),
               ('volt_rsqsum', dt_double, (361,)),
               ('volt_isqsum', dt_double, (361,)),
               ('sum_quots', dt_double),
               ('sum_wts', dt_double),
               ('sum_wts_sq', dt_double),
               ('wtd_disp_tot', dt_double, (334,)),
               ('spares', dt_byte, (7,)),
               ('bin_total', dt_long),
               ('wtd_bin_total', dt_float),
               ('disp', dt_doublecomplex, (361,)),
               ('norms', dt_double, (2,)),
               ('temp', dt_double, (8,)),
               ('covar', dt_double, (1050,)),
               ('bin_spares', dt_byte, (3,))])

#TODO: rewrite dt_fil_cov so that everything is in one record and covariance matrices are more
#easily accessed
dt_fil_cov2 = np.dtype([('ct_head', dt_ct_head),
                        ('ident', dt_ident),
                        ('num_ifgs', dt_long),
                        ('adj_num_ifgs', dt_float),
                        ('num_coadds', dt_long),
                        ('avg_rvoltspec', dt_double, (361,)),
                        ('avg_ivoltspec', dt_double, (361,)),
                        ('volt_rsqsum', dt_double, (361,)),
                        ('volt_isqsum', dt_double, (361,)),
                        ('sum_quots', dt_double),
                        ('sum_wts', dt_double),
                        ('sum_wts_sq', dt_double),
                        ('wtd_disp_tot', dt_double, (334,)),
                        ('bin_total', dt_long, (256,)),
                        ('wtd_bin_total', dt_float, (256,)),
                        ('disp', dt_doublecomplex, (256, 361)),
                        ('norms', dt_double, (256, 2)),
                        ('temp', dt_double, (256, 8)),
                        ('realrealcovar', dt_double, (369, 369)),
                        ('realimagcovar', dt_double, (369, 361)),
                        ('imagimagcovar', dt_double, (361, 361))])

dt_combinations = np.dtype([('order_0_input', '<S4', (8,)),
                            ('order_1_input', '<S4', (8,)),
                            ('order_2_input', '<S4', (4,)),
                            ('order_3_input', '<S4', (2,)),
                            ('order_0_nifgs', dt_word, (8,)),
                            ('order_1_nifgs', dt_word, (8,)),
                            ('order_2_nifgs', dt_word, (4,)),
                            ('order_3_nifgs', dt_word, (2,)),
                            ('current_order', dt_byte),
                            ('combination_quality', dt_byte),
                            ('spares', dt_byte)])

dt_spec_datal = np.dtype([('model_ttag', '<S14'),
                          ('model_label', '<S40'),
                          ('combine_script', '<S80'),
                          ('responsivity', dt_float),
                          ('resp_sigma', dt_float),
                          ('time_constant', dt_float),
                          ('tc_sigma', dt_float),
                          ('fft_length', dt_long),
                          ('lofreq_bin', dt_long),
                          ('hifreq_bin', dt_long),
                          ('spec', dt_complex, (361,)),
                          ('real_var', dt_float, (361,)),
                          ('imag_var', dt_float, (361,)),
                          ('real_imag_var', dt_float, (361,)),
                          ('phase_corr', dt_float),
                          ('pc_sigma', dt_float),
                          ('qrad', dt_float),
                          ('qrad_sigma', dt_float),
                          ('ir_power', dt_float),
                          ('ir_power_sigma', dt_float),
                          ('calibrated', dt_byte),
                          ('coadd_vars', dt_byte),
                          ('fil_vars', dt_byte),
                          ('fsl_vars', dt_byte),
                          ('destriped', dt_byte),
                          ('spares', dt_byte, (24,)),
                          ('combinations', dt_combinations)])

#ct_head = 64
#coad_spec_headl = 374
#coad_spec_datal = 2042
dt_fsl_sky = np.dtype([('ct_head', dt_ct_head),
                       ('coad_spec_head', dt_coad_spec_headl),
                       ('coad_spec_data', dt_coad_spec_datal),
                       ('spec_data', dt_spec_datal),
                       ('en_stat', dt_en_stat),
                       ('en_analog', dt_en_analog),
                       ('en_sigma', dt_en_sigma),
                       ('en_tempdiff', dt_en_tempdiff, (2,)),
                       ('attitude', dt_attitude)])

dt_fsl_sky2 = np.dtype([('ct_head', dt_ct_head),
                        ('coad_spec_head', dt_coad_spec_headl),
                        ('coad_spec_data', dt_coad_spec_datal2),
                        ('spec_data', dt_spec_datal),
                        ('en_stat', dt_en_stat2),
                        ('en_analog', dt_en_analog2),
                        ('en_sigma', dt_en_sigma2),
                        ('en_tempdiff', dt_en_tempdiff, (2,)),
                        ('attitude', dt_attitude)])

dt_fex_dtrf = np.dtype([('trans', dt_float, (128,))])

dt_fishspec = np.dtype([('nifgs', dt_long),
                        ('adj_nifgs', dt_float),
                        ('time', dt_double),
                        ('temp', dt_float, (7,)),
                        ('tsigma', dt_float, (7,)),
                        ('volt', dt_double),
                        ('bias', dt_double),
                        ('gain', dt_long),
                        ('vspec', dt_doublecomplex, (360,)),
                        ('vsigma', dt_double, (360,))])

dt_fex_var = np.dtype([('gmt', '<S14'),
                       ('time', dt_adt),
                       ('channel', dt_byte),
                       ('scan_length', dt_byte),
                       ('scan_speed', dt_byte),
                       ('model_label', '<S40'),
                       ('galat_exc', dt_float),
                       ('nsky_ifgs', dt_float),
                       ('variances', dt_double, (361, )),
                       ('spares', dt_byte, (47,))])

#or should I just write names, size, offsets with variances/cvar pointing to same area
dt_fex_cvs = np.dtype([('gmt', '<S14'),
                       ('time', dt_adt),
                       ('channel', dt_byte),
                       ('scan_length', dt_byte),
                       ('scan_speed', dt_byte),
                       ('model_label', '<S40'),
                       ('galat_exc', dt_float),
                       ('ncal_ifgs', dt_float),
                       ('cvar', dt_double, (361, )),
                       ('spares', dt_byte, (47,))])

dt_hskp_head = np.dtype([('stat_monitor_cmd', dt_uword, (8, )),
                         ('ipdu_stat', dt_byte, (8, )),
                         ('dwell_stat', dt_byte, (2, )),
                         ('lvdt_stat', dt_byte, (2, )),
                         ('u_proc_stat', dt_byte, (4, )),
                         ('bol_cmd_bias', dt_byte, (4, ))])

dt_hskp_end = np.dtype([('ntch_flt_a', dt_byte, (5, )),
                        ('ntch_flt_b', dt_byte, (5, )),
                        ('tlm_qual_maj_frm', dt_byte),
                        ('power_a_status', dt_byte),
                        ('power_b_status', dt_byte),
                        ('lmac_analog_temp', dt_byte),
                        ('lmac_digital_temp', dt_byte),
                        ('spares', dt_byte, (17, ))])

names_temp_counts = ['ex_cal', 'sky_horn', 'ref_horn', 'iref_source', 'dihedral',
                     'bol_assem', 'mirror', 'cal_resist', 'ex_cal_segment', 'colimator',
                     'grt']

fmts_temp_counts = [dt_uword, dt_uword, dt_uword, dt_uword, dt_uword,
                    (dt_uword, (4, )), dt_uword, (dt_uword, (4, )), dt_uword, dt_uword,
                    (dt_uword, (16, ))]

off_temp_counts = [0, 2, 4, 6, 8,
                   10, 18, 20, 28, 30,
                   0]

dt_temp_counts = np.dtype({'names' : names_temp_counts,
                           'formats' : fmts_temp_counts,
                           'offsets' : off_temp_counts})

dt_temp_counts2 = np.dtype([('ex_cal', dt_uword),
                            ('sky_horn', dt_uword),
                            ('ref_horn', dt_uword),
                            ('iref_source', dt_uword),
                            ('dihedral', dt_uword),
                            ('bol_assem', dt_uword, (4, )),
                            ('mirror', dt_uword),
                            ('cal_resist', dt_uword, (4, )),
                            ('ex_cal_segment', dt_uword),
                            ('colimator', dt_uword)])

dt_temperature = np.dtype([('side_amp', dt_temp_counts, (4, )),
                           ('ipdu', dt_byte, (2, )),
                           ('drive_box_a', dt_byte),
                           ('chan_pre_amp', dt_byte),
                           ('chan_temp', dt_byte, (4, )),
                           ('stat_mon', dt_byte, (2, )),
                           ('hot_spot_current', dt_byte, (2, )),
                           ('optical_preamp', dt_byte),
                           ('drive_box_b', dt_byte)])

dt_temperature2 = np.dtype([('side_amp', dt_temp_counts2, (4, )),
                            ('ipdu', dt_byte, (2, )),
                            ('drive_box_a', dt_byte),
                            ('chan_pre_amp', dt_byte),
                            ('chan_temp', dt_byte, (4, )),
                            ('stat_mon', dt_byte, (2, )),
                            ('hot_spot_current', dt_byte, (2, )),
                            ('optical_preamp', dt_byte),
                            ('drive_box_b', dt_byte)])

names_v_and_i = ['dig_conv_n_15', 'dig_conv_p_15', 'dig_conv_p_5', 'ana_conv_p_15', 'ana_conv_n_15',
                 'bias_pre_reg', 'int_ps_p_28', 'int_ps_p_15', 'int_ps_n_15', 'int_ps_p_5',
                 'ipdu_volt',
                 'cur_bias_pre_reg', 'cur_ana_conv', 'cur_dig_conv', 'con_current', 'con_int_current',
                 'ipdu_amp',
                 'mtm_cal_motor', 'bol_bias_volt']

fmts_v_and_i = [(dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )),
                (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )),
                (dt_byte, (20, )),
                (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (2, )), (dt_byte, (4, )), (dt_byte, (2, )),
                (dt_byte, (12, )),
                (dt_byte, (2, )), (dt_byte, (4, ))]

off_v_and_i = [0, 2, 4, 6, 8,
               10, 12, 14, 16, 18,
               0,
               20, 22, 24, 26, 30,
               20,
               32, 34]

dt_v_and_i = np.dtype({'names' : names_v_and_i,
                       'formats' : fmts_v_and_i,
                       'offsets' : off_v_and_i})

dt_v_and_i2 = np.dtype([('ipdu_volt', dt_byte, (20, )),
                        ('ipdu_amp', dt_byte, (12, )),
                        ('mtm_cal_motor', dt_byte, (2, )),
                        ('bol_bias_volt', dt_byte, (4, ))])

dt_hskp_struct = np.dtype([('hskp_head', dt_hskp_head),
                           ('temps', dt_temperature),
                           ('v_and_i', dt_v_and_i)])

dt_hskp_struct2 = np.dtype([('hskp_head', dt_hskp_head),
                            ('temps', dt_temperature2),
                            ('v_and_i', dt_v_and_i2)])

dt_hskp_tail = np.dtype([('gmt_mjf', '<S14'),
                         ('spares', dt_word)])

dt_nfs_hkp = np.dtype([('ct_head', dt_ct_head),
                       ('frame', dt_hskp_struct, (2, )),
                       ('hskp_tail', dt_hskp_tail),
                       ('mj_frm', dt_hskp_end, (2, ))])

dt_nfs_hkp2 = np.dtype([('ct_head', dt_ct_head),
                        ('frame', dt_hskp_struct2, (2, )),
                        ('hskp_tail', dt_hskp_tail),
                        ('mj_frm', dt_hskp_end, (2, ))])

dt_fex_av_calrs = np.dtype([('ct_head', dt_ct_head),
                            ('data_stop', dt_byte, (14, )),
                            ('data_stop_time', dt_adt),
                            ('prev_data_start', dt_byte, (14, )),
                            ('ave_period', dt_float),
                            ('num_good_record', dt_long),
                            ('num_bad_record', dt_long),
                            ('calres_ave_a_lo', dt_float, (4, )),
                            ('calres_ave_a_hi', dt_float, (4, )),
                            ('calres_ave_b_lo', dt_float, (4, )),
                            ('calres_ave_b_hi', dt_float, (4, )),
                            ('calres_dev_a_lo', dt_float, (4, )),
                            ('calres_dev_a_hi', dt_float, (4, )),
                            ('calres_dev_b_lo', dt_float, (4, )),
                            ('calres_dev_b_hi', dt_float, (4, )),
                            ('calres_bad_a_lo', dt_float, (4, )),
                            ('calres_bad_a_hi', dt_float, (4, )),
                            ('calres_bad_b_lo', dt_float, (4, )),
                            ('calres_bad_b_hi', dt_float, (4, )),
                            ('spares', dt_byte, (48, ))])


