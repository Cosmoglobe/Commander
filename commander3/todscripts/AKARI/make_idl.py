month = "2006_07"

text_content = """
i = {i}
file_list = File_Search('/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/{month}/', '*LW*.fits.gz')
tsd = Obj_New('tsd')
r = tsd->Read_TSD(file_list[i])
pos_bs2detpix, tsd, elon_tmp, elat_tmp, /ECLIPTIC
aftime_tmp = tsd->Get_Value_TSD('FIS_OBS', 'AFTIME')
blankData_tmp = tsd->get_Flag_tsd('FIS_OBS', 'BLANK')
aocuFlag_tmp = tsd->get_status_tsd('AOCU', 'ADS_RATE_MODE_B0') OR $
                tsd->get_status_tsd('AOCU', 'ADS_RATE_MODE_B1') OR $
                tsd->get_status_tsd('GADS', 'BLANK')
flags_tmp = blankData_tmp OR aocuFlag_tmp
ccf = CCF_Read_CCF('OFSTL_001_002', tsd, CCF_FILE='/mn/stornext/d23/cmbco/cg/AKARI/ASTRO-F/reduction/CCF/data/OFSTL_001_002.ccf', HEADER=hdr, ERROR=err)
pos_bs2detpix, tsd, lon, lat
Path = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_pkl/lon/txt/'
OPENW, 1, Path + FILE_BASENAME(file_list[{i}], '.fits.gz') + '_lon_gads.txt'
printf, 1, lon, format='(4(D, TR1))'
CLOSE, 1
Path = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_pkl/lat/txt/'
OPENW, 3, Path + FILE_BASENAME(file_list[{i}], '.fits.gz') + '_lat_gads.txt'
printf, 3, lat, format='(4(D, TR1))'
CLOSE, 3
"""
file_name = 'lw_'+str(month)+'_gads.pro'
with open("/mn/stornext/d23/cmbco/cg/AKARI/code/idl/"+str(file_name), 'w') as file:
    for i in range(800):
        file.write(text_content.format(i=i, month=month))
        file.write('\n\n')  # 各文書の間に空行を追加
    file.write("END")


text_content = """
i = {i}
file_list = File_Search('/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/{month}/', '*SW*.fits.gz')
tsd = Obj_New('tsd')
r = tsd->Read_TSD(file_list[i])
pos_bs2detpix, tsd, elon_tmp, elat_tmp, /ECLIPTIC
aftime_tmp = tsd->Get_Value_TSD('FIS_OBS', 'AFTIME')
blankData_tmp = tsd->get_Flag_tsd('FIS_OBS', 'BLANK')
aocuFlag_tmp = tsd->get_status_tsd('AOCU', 'ADS_RATE_MODE_B0') OR $
                tsd->get_status_tsd('AOCU', 'ADS_RATE_MODE_B1') OR $
                tsd->get_status_tsd('PR', 'BLANK')
flags_tmp = blankData_tmp OR aocuFlag_tmp
ccf = CCF_Read_CCF('OFSTS_001_002', tsd, CCF_FILE='/mn/stornext/d23/cmbco/cg/AKARI/ASTRO-F/reduction/CCF/data/OFSTS_001_002.ccf', HEADER=hdr, ERROR=err)
pos_bs2detpix, tsd, lon, lat
Path = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_pkl/lon/txt/'
OPENW, 1, Path + FILE_BASENAME(file_list[{i}], '.fits.gz') + '_lon_gads.txt'
printf, 1, lon, format='(4(D, TR1))'
CLOSE, 1
Path = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_pkl/lat/txt/'
OPENW, 3, Path + FILE_BASENAME(file_list[{i}], '.fits.gz') + '_lat_gads.txt'
printf, 3, lat, format='(4(D, TR1))'
CLOSE, 3
"""

file_name = 'sw_'+str(month)+'_gads.pro'
with open("/mn/stornext/d23/cmbco/cg/AKARI/code/idl/"+str(file_name), 'w') as file:
    for i in range(800):
        file.write(text_content.format(i=i, month=month))
        file.write('\n\n')  # 各文書の間に空行を追加
    file.write("END")
print('DONE')

print(str(month))