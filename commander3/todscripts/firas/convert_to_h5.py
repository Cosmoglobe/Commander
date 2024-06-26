# Inputs include the engineering data and scientific data:
# FDQ_ENG=/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng.h5
# FDQ_SDF=/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf.h5
# Current files are direct numpy array dumps

import h5py

data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng.h5')
fname_out = '/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5'

d = data['fdq_eng'][()]

groups = ['ct_head', 'en_head', 'en_stat', 'en_analog', 'en_xcal', 'chan', 'en_tempdiff',
        'en_tail']

subgroups = [
        ['gmt', 'time', 'space_time', 'mjr_frm_no', 'orbit', 'hskp1_tlm_fmt',
        'hskp2_tlm_fmt', 'ingest_spares', 'dataset_id', 'instr_spares'],
        ['limits', 'convert', 'a_spare', 'key_id', 'number_of_records',
        'enghead_spares', 'sci_time/bin_time', 'first_eng_time',
        'last_eng_time', 'dq_summary_flag', 'ifg_no', 'att_summary_flag',
        'enghead_spares2'],
        ['group1', 'group2', 'hot_spot_cmd', 'engstat_spares', 'power_a_status',
            'power_b_status', 'engstat_spares2'],
        ['grt', 'group1'],
        ['pos', 'xcal_spares'],
        ['up_sci_mode', 'fakeit', 'up_adds_per_group', 'up_swps_per_ifg',
        'xmit_mtm_speed', 'xmit_mtm_len', 'sci_gain', 'dither', 'setup_spares'],
        ['xcal', 'ical', 'skyhorn', 'refhorn', 'dihedral', 'collimator_mirror',
            'bol_assem'],
        ['engtail_spares', 'hskp_flag', 'lmac_analog_temp', 'lmac_digital_temp',
            'tlm_qual_maj_frm', 'eng_spares']
        ]

names_group1 = [
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
                 'stat_word_16']
names_group2 = [
                 'grt_addr',
                 'micro_stat_bus',
                 'bol_cmd_bias',
                 'dwell_stat',
                 'lvdt_stat',
                 'hot_spot_cmd'
                 ]

names_en_analog_grt = [
                   'a_lo_xcal_tip', 'a_lo_skyhorn', 'a_lo_refhorn', 'a_lo_ical', 'a_lo_dihedral',
                   'a_lo_bol_assem', 'a_lo_mirror', 'a_lo_cal_resistors', 'a_lo_xcal_cone', 'a_lo_collimator',
                   'a_hi_xcal_tip', 'a_hi_skyhorn', 'a_hi_refhorn', 'a_hi_ical', 'a_hi_dihedral',
                   'a_hi_bol_assem', 'a_hi_mirror', 'a_hi_cal_resistors', 'a_hi_xcal_cone', 'a_hi_collimator',
                   'b_lo_xcal_tip', 'b_lo_skyhorn', 'b_lo_refhorn', 'b_lo_ical', 'b_lo_dihedral',
                   'b_lo_bol_assem', 'b_lo_mirror', 'b_lo_cal_resistors', 'b_lo_xcal_cone', 'b_lo_collimator',
                   'b_hi_xcal_tip', 'b_hi_skyhorn', 'b_hi_refhorn', 'b_hi_ical', 'b_hi_dihedral',
                   'b_hi_bol_assem', 'b_hi_mirror', 'b_hi_cal_resistors',
                   'b_hi_xcal_cone', 'b_hi_collimator']
names_en_analog_group1 = [
                   'temp_ctl', 'ipdu_temp', 'cna_temp', 'dbx_tmp', 'stat_mon_temp',
                   'pamp_chan', 'pamp_op', 'hot_spot', 'mtm_cal_mtr', 'mtm_pos', 'bol_volt',
                   'ipdu_bolt', 'ipdu_amp']


lens_grt =  [
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1]
lens_group1 =    [
                  8,2,4,2,2,
                  1,1,2,2,2,4,
                  20, 12]



f = h5py.File(fname_out, 'w')
for i in range(len(groups)):
    for j in range(len(subgroups[i])):
        if subgroups[i][j] == 'sci_time/bin_time':
            f.create_dataset('en_head/sci_time/bin_time',
                    data=d['en_head']['sci_time']['bin_time'])
        elif (i == 2) & (subgroups[i][j] == 'group1'):
            for k in range(len(names_group1)):
                f.create_dataset(f'{groups[i]}/{names_group1[k]}',
                        data=d['en_stat']['group1'][:,k])
        elif (i == 2) & (subgroups[i][j] == 'group2'):
            f.create_dataset(f'{groups[i]}/grt_addr',
                    data=d['en_stat']['group2'][:,0:2])
            f.create_dataset(f'{groups[i]}/micro_stat_bus',
                    data=d['en_stat']['group2'][:,2:6])
            f.create_dataset(f'{groups[i]}/bol_cmd_bias',
                    data=d['en_stat']['group2'][:,6:10])
            f.create_dataset(f'{groups[i]}/dwell_stat',
                    data=d['en_stat']['group2'][:,10:12])
            f.create_dataset(f'{groups[i]}/lvdt_stat',
                    data=d['en_stat']['group2'][:,12:14])
        elif (i == 3) & (subgroups[i][j] == 'group1'):
            ind0 = 0
            for k in range(len(lens_group1)):
                f.create_dataset(f'{groups[i]}/{names_en_analog_group1[k]}',
                        data=d['en_analog']['group1'][:,ind0:ind0+lens_group1[k]])
                ind0 += lens_group1[k]
        elif (i == 3) & (subgroups[i][j] == 'grt'):
            ind0 = 0
            for k in range(len(lens_grt)):
                f.create_dataset(f'{groups[i]}/{names_en_analog_grt[k]}',
                        data=d['en_analog']['grt'][:,ind0:ind0+lens_grt[k]])
                ind0 += lens_grt[k]

        else:
            f.create_dataset(f'{groups[i]}/{subgroups[i][j]}',
                    data=d[groups[i]][subgroups[i][j]])
f.close()


data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf.h5')
fname_out = '/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5'
#     lh, ll, rh, rl, refer to combinations of left/right and high/low
keys = ['fdq_sdf_lh', 'fdq_sdf_ll', 'fdq_sdf_rh', 'fdq_sdf_rl']
groups = ['ct_head', 'sci_head', 'ifg_data', 'dq_data', 'collect_time', 'attitude']
subgroups = [
        ['gmt', 'time', 'space_time', 'mjr_frm_no', 'orbit', 'hskp1_tlm_fmt',
        'hskp2_tlm_fmt', 'ingest_spares', 'dataset_id', 'instr_spares'],
        ['chan_id', 'gain', 'mtm_speed', 'mtm_length', 'data_qual',
            'data_ready', 'sc_head0', 'sc_head1a', 'sc_head1b'] + [f'sc_head{i}' for i in range(2, 26)],
        ['ifg', 'gltch'],
        ['fake', 'xcal_pos', 'iref_temp', 'eng_time', 'eng_rec', 'data_quality',
            'ifg_no', 'dq_spares'],
        ['midpoint_time', 'badtime_flag', 'fpp_spare'],
        ['pixel_no', 'equatorial', 'ra', 'dec',
         'terr_pixel_no', 'terr_latitude', 'terr_longitude', 'earth_limb',
         'earth_limb_azimuth', 'sun_angle', 'moon_angle', 'moon_az_angle',
         'moon_phase', 'sun_moon_dist', 'cobe_moon_dist', 'altitude',
         'projected_barycentric_velocity', 'mcilwain_l_param',
         'galactic_longitude', 'galactic_latitude',
         'ecliptic_longitude', 'ecliptic_latitude',
         'orbital_phase', 'projected_geocentric_velocity', 'scan_angle',
         'sc_rotation_angle', 'solution', 'pixel_definition', 'skymap_index',
         'exc_galactic_lat', 'terr_rad_byte', 'outside_galaxy_cut',
         'att_spares']
        ]

f = h5py.File(fname_out, 'w')
for k in keys:
    d = data[k][()]
    for i in range(len(groups)):
        for j in range(len(subgroups[i])):
            f.create_dataset(f'{k}/{groups[i]}/{subgroups[i][j]}',
                    data=d[groups[i]][subgroups[i][j]])
f.close()
