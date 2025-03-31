import os
import pickle
import numpy as np
import time
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import healpy as hp
from astropy.coordinates import SkyCoord, Galactic
from multiprocessing import Pool
import multiprocessing as mp
from functools import partial
import math
from astropy.wcs import WCS
import astropy.units as u
import glob
import itertools
import fitsio
from pathlib import Path
from datetime import datetime, timezone
start = time.time()

# # # IDLで作成したtxtファイルをpklファイルに変換する
# # # lat
# # def process_file(file_path, month):
# #     print(file_path)
# #     with open(file_path) as f:
# #         l = [s.rstrip() for s in f.readlines()]
# #     data_str = []
# #     data_2dlist = []
# #     lat = []
# #     for a in range(len(l)):
# #         if l[a] == ' ':
# #             pass
# #         else:
# #             data_str.append(l[a])
# #     for m in range(len(l)):
# #         data_2dlist.append([float(num) for num in data_str[m].split()])
# #     lat_list = list(itertools.chain.from_iterable(data_2dlist))
# #     for y in range(100):  # LW,SWで変更点あり
# #         lat.append(lat_list[y::100])  # LW,SWで変更点あり
# #     result_file_name = file_path.replace('.txt', '.pkl')
# #     with open(result_file_name, 'wb') as p:
# #         pickle.dump(lat, p)

# # def main():
# #     month = '2007_01'
# #     directory_path = f'/mn/stornext/d23/cmbco/cg/AKARI/test/lat/all'
# #     txt_files = glob.glob(os.path.join(directory_path, 'FIS_SW_*.txt'))  # LW,SWで変更点あり
# #     txt_files = sorted(txt_files)
# #     with Pool(processes=25) as pool:
# #         pool.starmap(process_file, [(file_path, month) for file_path in txt_files])
# #     end = time.time()
# #     time_diff = end - start
# #     print(time_diff)

# # if __name__ == '__main__':
# #     main()
# # # end #

# # # lon
# # def process_file(file_path, month):
# #     print(file_path)
# #     with open(file_path) as f:
# #         l = [s.rstrip() for s in f.readlines()]
# #     data_str = []
# #     data_2dlist = []
# #     lat = []
# #     for a in range(len(l)):
# #         if l[a] == ' ':
# #             pass
# #         else:
# #             data_str.append(l[a])
# #     for m in range(len(l)):
# #         data_2dlist.append([float(num) for num in data_str[m].split()])
# #     lat_list = list(itertools.chain.from_iterable(data_2dlist))
# #     for y in range(100):  # LW,SWで変更点あり
# #         lat.append(lat_list[y::100])  # LW,SWで変更点あり
# #     result_file_name = file_path.replace('.txt', '.pkl')
# #     with open(result_file_name, 'wb') as p:
# #         pickle.dump(lat, p)

# # def main():
# #     month = '2006_05'
# #     directory_path = f'/mn/stornext/d23/cmbco/cg/AKARI/test/lon/all'
# #     txt_files = glob.glob(os.path.join(directory_path, 'FIS_SW_*.txt'))  # LW,SWで変更点あり
# #     txt_files = sorted(txt_files)
# #     with Pool(processes=25) as pool:
# #         pool.starmap(process_file, [(file_path, month) for file_path in txt_files])
# #     end = time.time()
# #     time_diff = end - start
# #     print(time_diff)

# # if __name__ == '__main__':
# #     main()
# # # end #





# # # 「あかり」FISデータから直接、各種データを読み込みpklファイルとして保存する
# # # monthとLW,SWは別のターミナルで実行させる
# def process_file(file_name):  # 処理を行う関数を定義
#     file_path = os.path.join(path_add, file_name)
#     data = fitsio.FITS(file_path)
#     # d = data[1].read()
#     # status = d['STATUS']
#     # flag = d['FLAG']
#     # pixel = d['PIX_FLAG']
#     # det = d['DET']
#     # flux = d['FLUX']

#     print(file_name)

#     # t = np.array(status).reshape(len(status), len(status[0]))  # status
#     # status = t.T
#     # shtop, calalon, calason, calbon, sinalon, sinason = [], [], [], [], [], []
#     # for i in range(100):  # LWの場合は75, SWの場合は100に変更する
#     #     # shtop.append(status[1])
#     #     calalon.append(status[16])
#     #     calason.append(status[17])
#     #     # calbon.append(status[18])
#     #     # sinalon.append(status[19])
#     #     # sinason.append(status[20])
#     # # save_to_file(file_path, month, 'shtop', 'status', shtop)
#     # # save_to_file(file_path, month, 'calbon', 'status', calbon)
#     # # save_to_file(file_path, month, 'sinalon', 'status', sinalon)
#     # # save_to_file(file_path, month, 'sinason', 'status', sinason)
#     # for i in range(len(calalon)):
#     #     new_row = calalon[i].copy()
#     #     for j in range(len(calalon[i]) - 1):
#     #         if calalon[i][j] == True and calalon[i][j + 1] == False:
#     #             if j + 1 < len(new_row):
#     #                 new_row[j + 1] = True
#     #             if j + 2 < len(new_row):
#     #                 new_row[j + 2] = True
#     #             if j + 3 < len(new_row):
#     #                 new_row[j + 3] = True
#     #     calalon[i] = new_row
#     # save_to_file(file_path, month, 'calalon_re', 'status', calalon)
#     # for i in range(len(calason)):
#     #     new_row = calason[i].copy()
#     #     for j in range(len(calason[i]) - 1):
#     #         if calason[i][j] == True and calason[i][j + 1] == False:
#     #             if j + 1 < len(new_row):
#     #                 new_row[j + 1] = True
#     #             if j + 2 < len(new_row):
#     #                 new_row[j + 2] = True
#     #             if j + 3 < len(new_row):
#     #                 new_row[j + 3] = True
#     #             if j + 4 < len(new_row):
#     #                 new_row[j + 4] = True
#     #             if j + 5 < len(new_row):
#     #                 new_row[j + 5] = True
#     #             if j + 6 < len(new_row):
#     #                 new_row[j + 6] = True
#     #             if j + 7 < len(new_row):
#     #                 new_row[j + 7] = True
#     #             if j + 8 < len(new_row):
#     #                 new_row[j + 8] = True
#     #             if j + 9 < len(new_row):
#     #                 new_row[j + 9] = True
#     #     calason[i] = new_row
#     # save_to_file(file_path, month, 'calason_re', 'status', calason)



#     # t = np.array(flag).reshape(len(flag), len(flag[0]))  # flag
#     # flag = t.T
#     # # bad_frame = [flag[0] for _ in range(100)]  # LWの場合は75, SWの場合は100に変更する
#     # # save_to_file(file_path, month, 'bad_frame', 'flag', bad_frame)
#     # in_saa = [flag[0] for _ in range(75)]  # LWの場合は75, SWの場合は100に変更する
#     # save_to_file(file_path, month, 'in_saa', 'flag', in_saa)
#     # near_moon = [flag[0] for _ in range(75)]  # LWの場合は75, SWの場合は100に変更する
#     # save_to_file(file_path, month, 'near_moon', 'flag', near_moon)



#     # t = np.array(pixel).reshape(len(pixel), (len(pixel[0]) * len(pixel[0][0])))  # pixel_flag
#     # pixel = t.T
#     # mtgl_tail,flutter,no_rp_corr,blank,bad,dead,saturate,reset,rstanom,no_diff,gpgl_type1,gpgl_type2,gpgl_type3,gpgl_type4,gpgl_tail,mtgl_type1,mtgl_type2,mtgl_type3,mtgl_type4,no_peri_corr = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
#     # for i in range(75):  # LWの場合は75, SWの場合は100に変更する
#     #     # bad.append(pixel[0])
#     #     # dead.append(pixel[3])
#     #     # saturate.append(pixel[4])
#     #     reset.append(pixel[5])
#     #     # rstanom.append(pixel[6])
#     #     # no_diff.append(pixel[7])
#     #     # gpgl_type1.append(pixel[16])
#     #     # gpgl_type2.append(pixel[17])
#     #     # gpgl_type3.append(pixel[18])
#     #     # gpgl_type4.append(pixel[19])
#     #     # gpgl_tail.append(pixel[20])
#     #     # mtgl_type1.append(pixel[21])
#     #     # mtgl_type2.append(pixel[22])
#     #     # mtgl_type3.append(pixel[23])
#     #     # mtgl_type4.append(pixel[24])
#     #     # no_peri_corr.append(pixel[26])
#     #     # mtgl_tail.append(pixel[25])
#     #     # flutter.append(pixel[30])
#     #     # no_rp_corr.append(pixel[8])
#     #     # blank.append(pixel[31])
#     # # save_to_file(file_path, month, 'bad', 'pix_flag', bad)
#     # # save_to_file(file_path, month, 'dead', 'pix_flag', dead)
#     # # save_to_file(file_path, month, 'saturate', 'pix_flag', saturate)
#     # # save_to_file(file_path, month, 'rstanom', 'pix_flag', rstanom)
#     # # save_to_file(file_path, month, 'no_diff', 'pix_flag', no_diff)
#     # # save_to_file(file_path, month, 'gpgl_type1', 'pix_flag', gpgl_type1)
#     # # save_to_file(file_path, month, 'gpgl_type2', 'pix_flag', gpgl_type2)
#     # # save_to_file(file_path, month, 'gpgl_type3', 'pix_flag', gpgl_type3)
#     # # save_to_file(file_path, month, 'gpgl_type4', 'pix_flag', gpgl_type4)
#     # # save_to_file(file_path, month, 'gpgl_tail', 'pix_flag', gpgl_tail)
#     # # save_to_file(file_path, month, 'mtgl_type1', 'pix_flag', mtgl_type1)
#     # # save_to_file(file_path, month, 'mtgl_type2', 'pix_flag', mtgl_type2)
#     # # save_to_file(file_path, month, 'mtgl_type3', 'pix_flag', mtgl_type3)
#     # # save_to_file(file_path, month, 'mtgl_type4', 'pix_flag', mtgl_type4)
#     # # save_to_file(file_path, month, 'no_peri_corr', 'pix_flag', no_peri_corr)
#     # # save_to_file(file_path, month, 'mtgl_tail', 'pix_flag', mtgl_tail)
#     # # save_to_file(file_path, month, 'flutter', 'pix_flag', flutter)
#     # # save_to_file(file_path, month, 'no_rp_corr', 'pix_flag', no_rp_corr)
#     # # save_to_file(file_path, month, 'blank', 'pix_flag', blank)

#     # # for i in range(len(reset)):  # SWの場合のみ必要、LWの場合はコメントアウトしていい
#     # #     if i >= 60:
#     # #     new_row = reset[i].copy()
#     # #     for j in range(len(reset[i]) - 1):
#     # #         if reset[i][j] == True and reset[i][j + 1] == False:
#     # #             new_row[j + 1] = True
#     # #     reset[i] = new_row
#     # # save_to_file(file_path, month, 'reset_re', 'pix_flag', reset)

#     # for i in range(len(reset)):
#     #     new_row = reset[i].copy()
#     #     for j in range(len(reset[i]) - 1):
#     #         if reset[i][j] == True and reset[i][j + 1] == False:
#     #             if j + 1 < len(new_row):
#     #                 new_row[j + 1] = True
#     #             if j + 2 < len(new_row):
#     #                 new_row[j + 2] = True
#     #             if j + 3 < len(new_row):
#     #                 new_row[j + 3] = True
#     #     reset[i] = new_row
#     # save_to_file(file_path, month, 'reset_re_LW', 'pix_flag', reset)


#     # t = np.array(det).reshape(len(det), len(det[0]))  # det
#     # det = t.T
#     # save_to_file(file_path, month, 'det', 'det', det)


#     # t = np.array(flux).reshape(len(flux), len(flux[0]))  # flux
#     # flux = t.T
#     # for i in range(len(flux)):
#     #     for j in range(len(flux[0])):
#     #         if -100 < flux[i][j] < 500e3:
#     #             pass
#     #         else:
#     #             flux[i][j] = 0
#     # save_to_file(file_path, month, 'flux_re', 'flux', flux)


#     # ga = data[6].read()
#     # aa_lun = ga['AA_LUN']  # GADS aa_lun
#     # for i in range(len(aa_lun)):
#     #     if aa_lun[i] < 21.5:  # ImageData.hでは35deg
#     #         aa_lun[i] = 1
#     #     else:
#     #         aa_lun[i] = 0
#     # lune = []
#     # for i in range(75):  # LWの場合は75, SWの場合は100に変更する
#     #     lune.append(aa_lun)  # 0: survey/maneuver 1:pointing
#     # save_to_file(file_path, month, 'aa_lun2', 'gads', lune)


#     pr = data[7].read()
#     # radec = pr['RA']  # PR RA
#     # for i in range(len(radec)):
#     #     if radec[i] == 0:
#     #         radec[i] = 1
#     #     else:
#     #         radec[i] = 0
#     # lune = []
#     # for i in range(75):  # LWの場合は75, SWの場合は100に変更する
#     #     lune.append(radec)  # 0: survey/maneuver 1:pointing
#     # save_to_file(file_path, month, 'radec', 'gads', lune)

#     aa_lun = ga['AA_LUN']  # GADS aa_lun
#     for i in range(len(aa_lun)):
#         if aa_lun[i] < 21.5:  # ImageData.hでは35deg
#             aa_lun[i] = 1
#         else:
#             aa_lun[i] = 0
#     lune = []
#     for i in range(75):  # LWの場合は75, SWの場合は100に変更する
#         lune.append(aa_lun)  # 0: survey/maneuver 1:pointing
#     save_to_file(file_path, month, 'aa_lun2', 'gads', lune)


# def save_to_file(file_path, month, suffix, direc, data):  # ファイルの保存を行う関数
#     result_file_name = file_path.replace('_gb.fits.gz', f'_{suffix}.pkl').replace(
#         f'akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/{month}', f'akari_TSD_pkl/{direc}/all')
#     with open(result_file_name, 'wb') as p:  # suffixだけではなくて、ディレクトリとファイル名を別に選択させる
#         pickle.dump(data, p)

# if __name__ == '__main__':
#     directory = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/'
#     months = ['2006_04','2006_05','2006_06','2006_07','2006_08','2006_09','2006_10','2006_11','2006_12','2007_01','2007_02','2007_03','2007_04','2007_05','2007_06','2007_07','2007_08']
#     for month in months:
#         path_add = str(directory) + str(month)
#         files = [f for f in os.listdir(path_add) if f.startswith('FIS_LW_')]  # LW,SWを変更 -> 上にも変更点がある

#         with Pool(processes=50) as pool:
#             pool.map(process_file, files)

#     end = time.time()
#     time_diff = end - start
#     print(time_diff)

#     print("done "+str(month))










def moonag(mjd, DEG2RAD):
    """
    月の計算に必要な各角度・パラメータを算出する。
    (mjd: Modified Julian Date)
    戻り値: ta, a, b, c, d, e, g, j, l, m, n, v, w (すべてラジアン)
    """
    ta = (mjd - 15019.5) / 36525.0
    tb = ta * ta

    a = DEG2RAD * (ta * 4.77e5   + 296.1044608 + ta * 198.849108  + tb * 0.009192)
    b = DEG2RAD * (ta * 483120.0  + 11.250889   + ta * 82.02515    - tb * 0.003211)
    c = DEG2RAD * (ta * 480960.0  + 270.434164  + ta * 307.883142  - tb * 0.001133)
    d = DEG2RAD * (ta * 444960    + 350.737486  + ta * 307.114217  - tb * 0.001436)
    e = DEG2RAD * (ta * 35640     + 98.998753   + ta * 359.372886)
    g = DEG2RAD * (ta * 35999.04975 + 358.475833 - tb * 1.5e-4)
    j = DEG2RAD * (ta * 2880      + 225.444651  + ta * 154.906654)
    l = DEG2RAD * (ta * 36000.76892 + 279.696678 + tb * 3.03e-4)
    m = DEG2RAD * (ta * 19080     + 319.529425  + ta * 59.8585    + tb * 1.81e-4)
    n = DEG2RAD * (259.183275     - ta * 1800   - ta * 134.142008 + tb * 0.002078)
    v = DEG2RAD * (ta * 58320     + 212.603219 + ta * 197.803875  + tb * 0.001286)
    w = DEG2RAD * (ta * 58320     + 342.767053 + ta * 199.211911 * 3.1e-4 * tb)
    return ta, a, b, c, d, e, g, j, l, m, n, v, w


def moonth(ta, a, b, c, d, e, g, j, l, m, n, v, w):
    """
    月の位置ベクトルに関する内部変数 (mx, my, mz) を計算する。
    すべての引数、戻り値はラジアン（または無次元値）である。
    """
    # MOON THETA (mx)
    mx  = -0.00101 * math.sin(a + b - 4 * d)
    mx -= 0.00102 * math.sin(a - b - 4 * d - n)
    mx -= ta * 0.00103 * math.sin(a - b - n)
    mx -= 0.00107 * math.sin(a - g - b - 2 * d - n)
    mx -= 0.00121 * math.sin(2 * a - b - 4 * d - n)
    mx += 0.0013  * math.sin(3 * a + b + n)
    mx -= 0.00131 * math.sin(a + b - n)
    mx += 0.00136 * math.sin(a + b - d + n)
    mx -= 0.00145 * math.sin(g + b)
    mx -= 0.00149 * math.sin(a + g - b - 2 * d)
    mx += 0.00157 * math.sin(g - b + d - n)
    mx -= 0.00159 * math.sin(g - b)
    mx += 0.00184 * math.sin(a - g + b - 2 * d + n)
    mx -= 0.00194 * math.sin(b - 2 * d - n)
    mx -= 0.00196 * math.sin(g - b + 2 * d - n)
    mx += 0.002    * math.sin(b - d)
    mx -= 0.00205 * math.sin(a + g - b)
    mx += 0.00235 * math.sin(a - g - b)
    mx += 0.00246 * math.sin(a - 3 * b - n)
    mx -= 0.00262 * math.sin(2 * a + b - 2 * d)
    mx -= 0.00283 * math.sin(a + g + b - 2 * d)
    mx -= 0.00339 * math.sin(g - b - 2 * d - n)
    mx += 0.00345 * math.sin(a - b + n)
    mx -= 0.00347 * math.sin(g - b + 2 * d)
    mx -= 0.00383 * math.sin(b + d + n)
    mx -= 0.00411 * math.sin(a + g + b + n)
    mx -= 0.00442 * math.sin(2 * a - b - 2 * d - n)
    mx += 0.00449 * math.sin(a - b + 2 * d)
    mx -= 0.00456 * math.sin(3 * b - 2 * d + n)
    mx += 0.00466 * math.sin(a + b + 2 * d + n)
    mx += 0.0049  * math.sin(2 * a - b)
    mx += 0.00561 * math.sin(2 * a + b)
    mx += 0.00564 * math.sin(a - g + b + n)
    mx -= 0.00638 * math.sin(a + g - b - n)
    mx -= 0.00713 * math.sin(a + g - b - 2 * d - n)
    mx -= 0.00929 * math.sin(g + b - 2 * d)
    mx -= 0.00947 * math.sin(2 * a - b - n)
    mx += 0.00965 * math.sin(a - g - b - n)
    mx += 0.0097  * math.sin(b + 2 * d)
    mx += 0.01064 * math.sin(b - d + n)
    mx -= ta * 0.0125 * math.sin(b + n)
    mx -= 0.01434 * math.sin(g + b - 2 * d + n)
    mx -= 0.01652 * math.sin(a + g + b - 2 * d + n)
    mx -= 0.01868 * math.sin(2 * a + b - 2 * d + n)
    mx += 0.027   * math.sin(2 * a + b + n)
    mx -= 0.02994 * math.sin(a - b - 2 * d)
    mx -= 0.03759 * math.sin(g + b + n)
    mx -= 0.03982 * math.sin(g - b - n)
    mx += 0.04732 * math.sin(b + 2 * d + n)
    mx -= 0.04771 * math.sin(b - n)
    mx -= 0.06505 * math.sin(a + b - 2 * d)
    mx += 0.13622 * math.sin(a + b)
    mx -= 0.14511 * math.sin(a - b - 2 * d - n)
    mx -= 0.18354 * math.sin(b - 2 * d)
    mx -= 0.20017 * math.sin(b - 2 * d + n)
    mx -= 0.38899 * math.sin(a + b - 2 * d + n)
    mx += 0.40248 * math.sin(a - b)
    mx += 0.65973 * math.sin(a + b + n)
    mx += 1.96763 * math.sin(a - b - n)
    mx += 4.95372 * math.sin(b)
    mx += 23.89684 * math.sin(b + n)
    
    # MOON RHO (my)
    my  =  0.05491 * math.cos(2 * a + g)
    my +=  0.0629  * math.cos(a + d)
    my -=  0.06444 * math.cos(4 * d)
    my -=  0.06652 * math.cos(2 * a - g)
    my -=  0.07369 * math.cos(g - 4 * d)
    my +=  0.08119 * math.cos(a - 3 * d)
    my -=  0.09261 * math.cos(a + 4 * d)
    my +=  0.10177 * math.cos(a - 2 * b + 2 * d)
    my +=  0.10225 * math.cos(a + g + 2 * d)
    my -=  0.10243 * math.cos(a + 2 * g - 2 * d)
    my -=  0.12291 * math.cos(2 * b)
    my -=  0.12291 * math.cos(2 * a - 2 * b)
    my -=  0.12428 * math.cos(a + g - 4 * d)
    my -=  0.14986 * math.cos(3 * a)
    my -=  0.1607  * math.cos(a - g + 2 * d)
    my -=  0.16949 * math.cos(a - d)
    my +=  0.17697 * math.cos(a + 2 * b - 2 * d)
    my -=  0.18815 * math.cos(2 * a - 4 * d)
    my -=  0.19482 * math.cos(2 * g - 2 * d)
    my +=  0.22383 * math.cos(2 * b - 2 * d)
    my +=  0.22594 * math.cos(3 * a - 2 * d)
    my +=  0.24454 * math.cos(2 * a + g - 2 * d)
    my -=  0.31717 * math.cos(g + d)
    my -=  0.36333 * math.cos(a - 4 * d)
    my +=  0.47999 * math.cos(a - g - 2 * d)
    my +=  0.63844 * math.cos(g + 2 * d)
    my +=  0.8617  * math.cos(g)
    my +=  1.50534 * math.cos(a - 2 * b)
    my -=  1.67417 * math.cos(a + 2 * d)
    my +=  1.99463 * math.cos(a + g)
    my +=  2.07579 * math.cos(d)
    my -=  2.455   * math.cos(a - g)
    my -=  2.74067 * math.cos(a + g - 2 * d)
    my -=  3.83002 * math.cos(g - 2 * d)
    my -=  5.37817 * math.cos(2 * a)
    my +=  6.60763 * math.cos(2 * a - 2 * d)
    my -= 53.97626 * math.cos(2 * d)
    my -= 68.62152 * math.cos(a - 2 * d)
    my -= 395.13669 * math.cos(a)
    my += 3649.33705
    
    # MOON PHI (mz)
    mz  = -0.001    * math.sin(a - g - 2 * b - 2 * n)
    mz -= 0.001    * math.sin(a + g - 4 * d)
    mz += 0.001    * math.sin(2 * a - g)
    mz += 0.00102  * math.sin(a - g + 2 * d)
    mz -= 0.00106  * math.sin(2 * a - 2 * b - n)
    mz -= 0.00106  * math.sin(2 * a + n)
    mz -= 0.00109  * math.sin(a + 2 * b - 2 * d)
    mz -= 0.0011   * math.sin(2 * b - d + 2 * n)
    mz += 0.00112  * math.sin(4 * d)
    mz -= 0.00122  * math.sin(2 * a - n)
    mz -= 0.00122  * math.sin(2 * a + 2 * b + n)
    mz += 0.00149  * math.sin(g + 2 * b - 2 * d + 2 * n)
    mz -= 0.00157  * math.sin(2 * a - 4 * d)
    mz += 0.00171  * math.sin(a + g + 2 * b - 2 * d + 2 * n)
    mz -= 0.00175  * math.sin(2 * a + 2 * b - 2 * d)
    mz -= 0.0019   * math.sin(2 * g - 2 * d)
    mz += 0.00193  * math.sin(a + 16 * e - 18 * w)
    mz += 0.00194  * math.sin(2 * a + 2 * b - 2 * d + 2 * n)
    mz += 0.00201  * math.sin(g - 2 * d - n)
    mz += 0.00201  * math.sin(g + 2 * b - 2 * d + 2 * n)
    mz -= 0.00207  * math.sin(a + 2 * g - 2 * d)
    mz -= 0.0021   * math.sin(2 * g)
    mz -= 0.00213  * math.sin(2 * d - n)
    mz -= 0.00213  * math.sin(2 * b + 2 * d + n)
    mz -= 0.00215  * math.sin(3 * a - 2 * d)
    mz -= 0.00247  * math.sin(a - 4 * d)
    mz -= 0.00253  * math.sin(a - 2 * b + 2 * d)
    mz += 0.00279  * ta * math.sin(2 * b + 2 * n)
    mz -= 0.0028   * math.sin(2 * a + 2 * b + 2 * n)
    mz += 0.00312  * math.sin(3 * a)
    mz -= 0.00317  * math.sin(a + 2 * b)
    mz -= 0.0035   * math.sin(a + 16 * e - 18 * w)
    mz += 0.0039   * math.sin(g + 2 * b + 2 * n)
    mz += 0.00413  * math.sin(g - 2 * b - 2 * n)
    mz -= 0.0049   * math.sin(2 * n)
    mz -= 0.00491  * math.sin(2 * b + 2 * d + 2 * n)
    mz += 0.00504  * math.sin(g + d)
    mz += 0.00516  * math.sin(a - d)
    mz -= 0.00621  * math.sin(g + 2 * d)
    mz += 0.00648  * math.sin(a - 2 * b - 2 * d - n)
    mz += 0.00648  * math.sin(a - 2 * d + n)
    mz += 0.007    * math.sin(a - g - 2 * d)
    mz += 0.01122  * math.sin(a + 2 * d)
    mz += 0.0141   * math.sin(a - 2 * d - n)
    mz += 0.0141   * math.sin(a + 2 * b - 2 * d + n)
    mz += 0.01424  * math.sin(a - 2 * b)
    mz += 0.01506  * math.sin(a - 2 * b - 2 * d - 2 * n)
    mz -= 0.01567  * math.sin(2 * b - 2 * d)
    mz += 0.02077  * math.sin(2 * b - 2 * d + n)
    mz -= 0.02527  * math.sin(a + g)
    mz -= 0.02952  * math.sin(a - n)
    mz -= 0.02952  * math.sin(a + 2 * b + n)
    mz -= 0.03487  * math.sin(d)
    mz += 0.03684  * math.sin(a - g)
    mz -= 0.03983  * math.sin(2 * d + n)
    mz += 0.03983  * math.sin(2 * b - 2 * d + n)
    mz += 0.04037  * math.sin(a + 2 * b - 2 * d + 2 * n)
    mz += 0.04221  * math.sin(2 * a)
    mz -= 0.04273  * math.sin(g - 2 * d)
    mz -= 0.05566  * math.sin(2 * a - 2 * d)
    mz -= 0.05697  * math.sin(a + g - 2 * d)
    mz -= 0.06846  * math.sin(a + 2 * b + 2 * n)
    mz -= 0.08724  * math.sin(a - 2 * b - n)
    mz -= 0.08724  * math.sin(a + n)
    mz -= 0.11463  * math.sin(2 * b)
    mz -= 0.18647  * math.sin(g)
    mz -= 0.20417  * math.sin(a - 2 * b - 2 * n)
    mz += 0.59616  * math.sin(2 * d)
    mz += 1.07142  * math.sin(n)
    mz -= 1.07447  * math.sin(2 * b + n)
    mz -= 1.28658  * math.sin(a - d)
    mz -= 2.4797   * math.sin(2 * b + 2 * n)
    mz += 6.32962  * math.sin(a)

    return mx, my, mz


def at_moon(mjd, EARTH_RADIUS, MOON_RADIUS, DEG2RAD, TWO_PI):
    """
    月の位置を計算する関数
      入力:
         mjd: Modified Julian Date
      戻り値:
         pos: 地球中心から月へのベクトル [km]（J2000座標系）
         size: 月の見かけの角直径 [radian]
         phase: 月の位相 [radian] (0:new, pi:full)
         dist: 月までの距離 [km]
    """
    ta, a, b, c, d, e, g, j, l, m, n, v, w = moonag(mjd, DEG2RAD)
    mx, my, mz = moonth(ta, a, b, c, d, e, g, j, l, m, n, v, w)
    
    # 以下の計算は元コードの通り
    # ※注意: 元コード中、myは実際には月のρ²に関する値として用いられていると思われるため、そのまま再現
    r_xy = math.sqrt(my - mx * mx)
    sin_delta = mz / r_xy
    cos_delta = math.sqrt(1.0 - sin_delta * sin_delta)
    sin_c = math.sin(c)
    cos_c = math.cos(c)
    
    dist = EARTH_RADIUS * math.sqrt(my)
    # x_tod は大気層外での月の位置ベクトル（J2000座標系とする）\n    x_tod[0] = EARTH_RADIUS * r_xy * (cos_delta * cos(c) - sin_delta * sin(c))\n    x_tod[1] = EARTH_RADIUS * r_xy * (sin_delta * cos(c) + cos_delta * sin(c))\n    x_tod[2] = EARTH_RADIUS * mx
    x_tod = [
        EARTH_RADIUS * r_xy * (cos_delta * cos_c - sin_delta * math.sin(c)),
        EARTH_RADIUS * r_xy * (sin_delta * cos_c + cos_delta * math.sin(c)),
        EARTH_RADIUS * mx,
    ]
    
    size = math.atan(MOON_RADIUS / dist)
    phase = d % TWO_PI

    # atPrecessionの処理はここでは省略（J2000座標系のままとする）
    pos = x_tod

    return pos, size, phase, dist


def precession_rm(mjd, MJD_J2000, DEG2RAD):
    """
    与えられたMJDの歳差補正回転行列を計算する
    """
    t = (mjd - MJD_J2000) / 36525.0
    
    zeta  = (2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3) * DEG2RAD / 3600.0
    z     = (2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3) * DEG2RAD / 3600.0
    theta = (2004.3109 * t - 0.42665 * t**2 - 0.041833 * t**3) * DEG2RAD / 3600.0
    
    rm = [
        [math.cos(zeta) * math.cos(theta) * math.cos(z) - math.sin(zeta) * math.sin(z),
         -math.sin(zeta) * math.cos(theta) * math.cos(z) - math.cos(zeta) * math.sin(z),
         -math.sin(theta) * math.cos(z)],
        [math.cos(zeta) * math.cos(theta) * math.sin(z) + math.sin(zeta) * math.cos(z),
         -math.sin(zeta) * math.cos(theta) * math.sin(z) + math.cos(zeta) * math.cos(z),
         -math.sin(theta) * math.sin(z)],
        [math.cos(zeta) * math.sin(theta),
         -math.sin(zeta) * math.sin(theta),
         math.cos(theta)]
    ]
    return rm


def mat_mult(a, b):
    """
    3x3行列同士の積を計算する
    """
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat_vec_mult(mat, vec):
    """
    3x3行列と3次元ベクトルの積を計算する
    """
    return [sum(mat[i][j] * vec[j] for j in range(3)) for i in range(3)]


def at_precess_rm(mjd0, mjd, MJD_J2000, DEG2RAD):
    """
    MJD0 から MJD への歳差補正回転行列を計算する
    """
    rm_a_to_2000 = precession_rm(mjd0, MJD_J2000, DEG2RAD)
    rm_b_to_2000 = precession_rm(mjd, MJD_J2000, DEG2RAD)
    rm_2000_to_b = [[rm_b_to_2000[j][i] for j in range(3)] for i in range(3)]  # 転置で逆行列
    return mat_mult(rm_a_to_2000, rm_2000_to_b)


def at_precession(mjd0, x0, mjd, MJD_J2000, DEG2RAD):
    """
    MJD0 の時点での赤道座標 x0 を、MJD における赤道座標に変換する
    """
    rm = at_precess_rm(mjd0, mjd, MJD_J2000, DEG2RAD)
    return mat_vec_mult(rm, x0)


def at_kepler(g, eccent):
# ----------------------------
# atKepler: ケプラー方程式を解く関数
# g: 平均近点角 [rad]
# eccent: 離心率
# 戻り値: 解となる離心近点角 [rad]
# ----------------------------
    eps = 1e-15
    imax = 50
    E = g
    # 平均近点角が0なら解は0とする
    if g == 0.0:
        return 0.0
    for i in range(imax):
        deltaE = (g - E + eccent * math.sin(E)) / (1 - eccent * math.cos(E))
        E += deltaE
        # Eが0にならない場合の相対誤差
        error = abs(deltaE / E) if E != 0 else abs(deltaE)
        if error < eps:
            return E
    raise Exception("Kepler equation did not converge.")


def at_orb_plane(xys, omeltl, omebig, aincln):
# ----------------------------
# atOrbPlane: 軌道面上の座標を地心座標系に変換する
# xys: 軌道面座標 (3次元ベクトル、xys[2]は通常0)
# omeltl: 近地点引数（argument of perigee） [rad]
# omebig: 昇交点経度（right ascension of ascending node） [rad]
# aincln: 軌道傾斜角 [rad]
# 戻り値: 地心座標系での3次元ベクトル
# ----------------------------
    coso = math.cos(omebig)
    sino = math.sin(omebig)
    cosi = math.cos(aincln)
    sini = math.sin(aincln)
    
    # C++コードにおける変換式（各項は軌道要素に依存）:
    # xyz[0] = xys[0]*(coso*cos(omeltl) - sino*cosi*sin(omeltl)) - xys[1]*(coso*sin(omeltl) + sino*cosi*cos(omeltl));
    # xyz[1] = xys[0]*(sino*cos(omeltl) + coso*cosi*sin(omeltl)) + xys[1]*(coso*cosi*cos(omeltl) - sino*sin(omeltl));
    # xyz[2] = xys[0]*(sini*sin(omeltl)) + xys[1]*(sini*cos(omeltl)) + xys[2]*cosi;
    xyz0 = xys[0] * (coso * math.cos(omeltl) - sino * cosi * math.sin(omeltl)) - \
           xys[1] * (coso * math.sin(omeltl) + sino * cosi * math.cos(omeltl))
    xyz1 = xys[0] * (sino * math.cos(omeltl) + coso * cosi * math.sin(omeltl)) + \
           xys[1] * (coso * cosi * math.cos(omeltl) - sino * math.sin(omeltl))
    xyz2 = xys[0] * (sini * math.sin(omeltl)) + \
           xys[1] * (sini * math.cos(omeltl)) + \
           xys[2] * cosi
    return [xyz0, xyz1, xyz2]


def at_sat_pos(mjd, atElement):
# ----------------------------
# atSatPos: 指定時刻 mjd における衛星位置を計算する
# 入力:
#   mjd: Modified Julian Date
#   atElement: グローバル変数atElementに含まれる軌道要素（辞書）\n
# 戻り値:
#   xyz: 地球中心から見た衛星の位置ベクトル（J2000座標系、単位: km）
# ----------------------------
    # 時間差（基準MJDとの差）
    timed = mjd - atElement["mjdz"]
    timem = timed * 1440.0
    semaxs = atElement["semiax"] + atElement["adot"] * timed
    omebig = atElement["ragome"] + atElement["ragdot"] * timed
    omeltl = atElement["smaome"] + atElement["smodot"] * timed
    amean0 = atElement["omean0"] + (atElement["znbadt"] * 0.5 * timed + atElement["znbar"]) * timem

    e = atElement["eccent"]
    # ケプラー方程式を解いて離心近点角 u を求める
    u = at_kepler(amean0, e)

    # 軌道面上での位置（xys）を計算
    xys = [0.0, 0.0, 0.0]
    xys[0] = semaxs * (math.cos(u) - e)
    xys[1] = semaxs * math.sqrt(1 - e * e) * math.sin(u)
    xys[2] = 0.0  # 軌道面上では通常 z成分は0

    # 軌道面から地心座標系への変換
    xyz = at_orb_plane(xys, omeltl, omebig, atElement["aincln"])
    return xyz


def at_mjulian_d(attime):
# ----------------------------
# 日付からユリウス日 (JD) を計算する（グレゴリオ暦版）  
# JD = 367*year - int((7*(year+int((month+9)/12)))/4) + int(275*month/9) + day + 1721013.5 + (hour+minute/60+second/3600)/24  
# 修正ユリウス日 (MJD) = JD - 2400000.5  
# ----------------------------
    """
    attime: 辞書 {"yr":..., "mo":..., "dy":..., "hr":..., "mn":..., "sc":...}
    戻り値: mjd (float)
    """
    year = attime["yr"]
    month = attime["mo"]
    day = attime["dy"]
    hour = attime["hr"]
    minute = attime["mn"]
    second = attime["sc"]
    
    # 月が1月または2月の場合、年と月を補正
    if month <= 2:
        year -= 1
        month += 12
    A = int(year/100)
    B = 2 - A + int(A/4)
    JD = int(365.25*(year+4716)) + int(30.6001*(month+1)) + day + B - 1524.5
    frac_day = (hour + minute/60.0 + second/3600.0) / 24.0
    JD += frac_day
    MJD = JD - 2400000.5
    return MJD


def at_mj_date(mjd):
# ----------------------------
# atMJDate: mjdを日時（文字列など）に変換する（ここでは簡易に文字列化）
# ----------------------------
    # JD = MJD + 2400000.5
    JD = mjd + 2400000.5
    # 単純な計算ではなく、datetime.fromtimestamp()などを使う方法もありますが、ここではダミー実装
    # ※実際の変換が必要な場合はastropy.time.Timeなどを利用してください。
    return f"MJD={mjd}"


def at_set_element(filename, mjd0, kchk, atElement, DEG2RAD, EARTH_RADIUS, TWO_PI):
# ----------------------------
# at_set_element: 軌道要素ファイルからatElementを設定する関数  
# filename: 軌道要素ファイルのパス  
# mjd0: 基準とするMJD（float）  
# kchk: チェックフラグ（整数、0なら2次微分を理論値で計算）  
# 戻り値: 0 (正常終了) または例外発生（エラー）
# ----------------------------
    # 定数（static変数相当）
    coef1  = 0.0234375
    coef2  = 0.3515625
    coef3  = 1.7916667
    coef4  = 1.25
    coef5  = 4.375
    coef6  = 6.642857
    coef7  = 1.7142857
    coef8  = 5.0625
    coef9  = 6.75
    coef10 = 1.928571
    coef11 = 1.666667
    coef12 = 0.2083333
    coef13 = 0.8571428
    coef14 = 1.5
    gravc  = 0.0743666
    j2     = 0.00108264
    days   = 86400.0
    dycon  = 398599.0
    aj2    = 0.001082628
    aj4    = -2.12e-6
    eradi  = EARTH_RADIUS  # 6378.16 km

    # attimeの初期化（7個の整数を格納するための辞書）
    attime = {"yr":0, "mo":0, "dy":0, "hr":0, "mn":0, "sc":0}

    ibuff = [0]*7   # 整数バッファ
    dbuff = [0.0]*11  # 実数バッファ

    mjdz = 0.0
    jblk = 0
    code2 = 0

    # ファイルを開く
    try:
        with open(filename, "r") as f:
            # print(f.read())
            # ファイルの各ブロックを読み、mjd1がmjd0以下の最新のブロックを採用する
            while True:
                # 7個の整数を読み込む
                line = f.readline()
                if not line:
                    break
                parts = line.split()
                if len(parts) < 7:
                    break
                try:
                    # ibuff[0]～ibuff[6] を整数として読み込み
                    ibuff = [int(x) for x in parts[:7]]
                except Exception:
                    break
                # attimeはibuffの2番目～7番目（インデックス1～6）を利用
                attime["yr"] = ibuff[1]
                attime["mo"] = ibuff[2]
                attime["dy"] = ibuff[3]
                attime["hr"] = ibuff[4]
                attime["mn"] = ibuff[5]
                attime["sc"] = ibuff[6]
                # MJDを計算
                mjd1 = at_mjulian_d(attime)
                if jblk == 0 or mjd1 <= mjd0:
                    # 次の行から11個の実数を読み込む
                    line_dbuff = f.readline()
                    if not line_dbuff:
                        break
                    parts_dbuff = line_dbuff.split()
                    if len(parts_dbuff) < 11:
                        break
                    try:
                        dbuff = [float(x) for x in parts_dbuff[7:18]]
                    except Exception:
                        break
                    mjdz = mjd1
                    jblk += 1
                if mjd1 > mjd0:
                    break
    except Exception as e:
        raise Exception("OPEN_ERROR: " + str(e))
    
    if jblk == 0:
        raise Exception("FILE_FORMAT_ERROR: no valid block found")
    
    # atElementの各項目を設定（atElementはグローバル辞書）
    atElement["mjdz"] = mjdz
    atElement["itz"] = at_mj_date(mjdz)  # 簡易な日時文字列
    atElement["semiax"] = dbuff[0]
    atElement["eccent"] = dbuff[1]
    atElement["aincln"] = dbuff[2] * DEG2RAD
    atElement["ragome"] = dbuff[3] * DEG2RAD
    atElement["smaome"] = dbuff[4] * DEG2RAD
    atElement["omean0"] = dbuff[5] * DEG2RAD
    atElement["adot"]   = dbuff[6]
    atElement["ragdot"] = dbuff[7] * DEG2RAD
    atElement["smodot"] = dbuff[8] * DEG2RAD
    atElement["eccdot"] = 0.0  # 変化率はサポートしていない
    atElement["aindot"] = 0.0

    # CONVERT KOZAI'S MEAN TO BROWSLER' MEAN
    d1 = math.sin(atElement["aincln"])
    d2 = math.sqrt(1 - atElement["eccent"]**2)
    d4 = atElement["semiax"] / EARTH_RADIUS
    # d3 = 1 - j2 * 1.5 * (1 - d1*d1*1.5) / (d2*(d2*d2)) / (d4*d4)
    # ※ C++版の式中「d1*d1*1.5」は解釈に注意（恐らく 1.5*(d1*d1)）と解釈する）
    d3 = 1 - j2 * 1.5 * (1 - d1*d1) / (d2**3) / (d4**2)
    atElement["semiax"] /= d3
    atElement["adot"]   /= d3

    # THEORETICAL VALUE OF TIME DERIVATIVES OF ELEMENTS
    d1_val = atElement["semiax"]
    d2_val = d1_val
    pp1 = math.sqrt(dycon / (d2_val * (d1_val**2)))
    ppp = TWO_PI / pp1
    e2 = atElement["eccent"]**2
    semrec = (1 - e2) * atElement["semiax"] / eradi
    eee = math.sqrt(1 - e2)
    cosi = math.cos(atElement["aincln"])
    sini = math.sin(atElement["aincln"])
    cosi2 = cosi * cosi
    sini2 = sini * sini
    con1 = aj2 / (semrec**2)
    aj22 = con1 * con1
    con1 *= 1.5
    con2 = aj4 / (semrec**4)
    
    zzbar = pp1 * ( con1 * eee * (1 - sini2 * 1.5) + 1 \
             + coef1 * aj22 * ( (((eee * 25.0 + 144.0) * eee + 105.0) * cosi2 \
             - (eee * 90.0 + 96.0) * eee - 30.0) * cosi2 \
             + (eee * 25.0 + 16.0) * eee - 15.0 ) * eee \
             - coef2 * con2 * e2 * eee * ((cosi2 * 35.0 - 30.0) * cosi2 + 3.0) )
    atElement["znbar"] = zzbar * 60.0

    if atElement["smodot"] == 0.0:
        atElement["smodot"] = ( con1 * zzbar * (2 - sini2 * 2.5) \
             * ( con1 * (e2 * 0.5 + 2 - eee * 2 - (coef3 - e2 / 48.0 - eee * 3) * sini2) + 1) \
             - coef4 * aj22 * e2 * pp1 * cosi2 * cosi2 \
             - coef5 * con2 * pp1 * ((sini2 * 5.25 - coef6) * sini2 + coef7 \
             + e2 * ((coef8 * sini2 - coef9) * sini2 + coef10)) ) * days

    if atElement["ragdot"] == 0.0:
        atElement["ragdot"] = -cosi * ( zzbar * con1 * ( con1 * (e2 / 6.0 + 1.5 - eee * 2 - (coef11 - coef12 * e2 - eee * 3) * sini2) + 1) \
             + con2 * 4.375 * pp1 * (e2 * 1.5 + 1) * (coef13 - coef14 * sini2) ) * days

    if atElement["adot"] != 0.0 and kchk == 0:
        semax0 = atElement["semiax"] / EARTH_RADIUS
        atElement["znbadt"] = gravc * -1.5 / math.sqrt(semax0) / (semax0**2) * atElement["adot"] / EARTH_RADIUS

    atElement["perige"] = atElement["semiax"] * (1 - atElement["eccent"]) - eradi
    atElement["apoge"]  = atElement["semiax"] * (atElement["eccent"] + 1) - eradi

    return atElement


def at_vect_to_pol(x, TWO_PI):
    """
    直交座標系ベクトル x を極座標系に変換する。
    入力:
        x: 3次元ベクトル [x, y, z]
    戻り値:
        y: 辞書 {'r': 半径, 'lon': 経度, 'lat': 緯度}
    """
    norm01 = x[0]**2 + x[1]**2
    r = math.sqrt(norm01 + x[2]**2)
    
    if r == 0.0:
        return {'r': 0.0, 'lon': 0.0, 'lat': 0.0}  # NULL_VECTOR の場合
    
    norm01 = math.sqrt(norm01)
    lat = math.atan2(x[2], norm01)  # 緯度
    lon = math.atan2(x[1], x[0])    # 経度
    if lon < 0.0:
        lon += TWO_PI
    
    return {'r': r, 'lon': lon, 'lat': lat}


def eps(mjd):
    """
    修正ユリウス日 mjd における黄道傾斜角 ε をラジアンで返す。
    ここでは J2000.0 (mjd=51544.5) からの簡易的な線形近似式を用いる例です。
    (ε = 23.439291° - 0.0130042° * T, T: 世紀数)
    """
    T = (mjd - 51544.5) / 36525.0
    epsilon_deg = 23.439291 - 0.0130042 * T
    return math.radians(epsilon_deg)


def equatorial2ecliptic(mjd, equatorial):
    """
    赤道座標 (equatorial) を黄道座標 (ecliptic) に変換する。
    
    Parameters:
      mjd         : 修正ユリウス日 (Modified Julian Date)
      equatorial  : [RA, Dec] (ラジアン) のリストまたはタプル
    
    Returns:
      ecliptic    : [λ, β] (黄道経度, 黄道緯度、ラジアン) のリスト
    """
    epsilon = eps(mjd)
    RA, Dec = equatorial[0], equatorial[1]

    # 以下、C++版の式に基づく変換
    work0 = math.atan2(math.sin(RA) * math.cos(epsilon) + math.tan(Dec) * math.sin(epsilon),
                       math.cos(RA))
    work1 = math.asin(math.sin(Dec) * math.cos(epsilon) - math.cos(Dec) * math.sin(epsilon) * math.sin(RA))
    
    if work0 < 0.0:
        work0 += 2.0 * math.pi

    return [work0, work1]


def angle(coord1, coord2):
    """
    2点の赤道座標 (longitude, latitude) [ラジアン] から角距離を計算する。
    
    Parameters:
      coord1: [lon1, lat1]
      coord2: [lon2, lat2]
      
    Returns:
      anglerad: 2点間の角距離 [ラジアン]
    """
    # 入力座標の各要素を取得
    lon1, lat1 = coord1
    lon2, lat2 = coord2

    # 球面上の2点間の内積から角度を計算
    val = math.sin(lat1) * math.sin(lat2) + math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2)
    # acosで角度を求める
    anglerad = math.acos(val)
    
    # 数値誤差の影響でvalが -1 または 1 に非常に近い場合の補正
    if val < -0.999999 or val > 0.999999:
        diff_lon = lon1 - lon2
        diff_lat = lat1 - lat2
        val2 = math.cos(lat1) * math.cos(lat2) * (diff_lon ** 2) + (diff_lat ** 2)
        anglerad = math.sqrt(val2)
    
    return anglerad


# 月の回避角 (moon_avoid) が 21.5° 未満なら即座に False に。
# それ以外で、21.5°以上34°未満の場合、さらに以下のいずれかの条件で False に
# 三つのパターンのうち、少なくとも一つで「距離」計算の結果が設定された上下限の間に入る。
# もしくは、moon_avoid が 23.0°～26.5°の範囲にあり、かつ算出された位相角が 27°～63°の間にある。
def calculate_obsflag(moon_avoid, m_ecliptic, ecliptic, obsflag, sample, image_param, DEG2RAD, RAD2DEG):
    test = []
    # 角度の定義
    moon20 = 21.5 * DEG2RAD
    moon23 = 23.0 * DEG2RAD
    moon27 = 26.5 * DEG2RAD
    moon34 = 34.0 * DEG2RAD
    phase27 = 27.0 * DEG2RAD
    phase63 = 63.0 * DEG2RAD
    
    # moon_avoid の範囲チェック
    if moon_avoid < moon20:
        # obsflag[sample] = False
        test.append(1)
    elif moon_avoid < moon34:
        l_diff = m_ecliptic[0] - ecliptic[0]
        while l_diff < -math.pi:
            l_diff += 2.0 * math.pi
        while l_diff > math.pi:
            l_diff -= 2.0 * math.pi
        
        b_diff = m_ecliptic[1] - ecliptic[1]
        
        # パラメータの定義
        r1 = 24.5 * DEG2RAD
        theta1 = 59.5 * DEG2RAD
        dist1_l = 23.5
        dist1_u = 25.7
        r2 = 29.4 * DEG2RAD
        theta2 = 179.4 * DEG2RAD
        dist2_l = 27.9
        dist2_u = 29.2
        r3 = 24.5 * DEG2RAD
        theta3 = -59.0 * DEG2RAD
        dist3_l = 23.5
        dist3_u = 25.5
        
        # 走査方向に基づく角度計算
        if image_param == 1:  # 0: to_sep, 1: to_nep
            phase_angle_rad = math.atan2(math.tan(b_diff), math.sin(l_diff))
        else:
            phase_angle_rad = math.atan2(-math.tan(b_diff), -math.sin(l_diff))
        
        phase_angle = phase_angle_rad * RAD2DEG
        
        # 距離の計算
        dist1 = math.acos((math.cos(moon_avoid) * math.cos(r1)) + 
                            (math.sin(moon_avoid) * math.sin(r1) * math.cos(phase_angle_rad - theta1))) * RAD2DEG
        dist2 = math.acos((math.cos(moon_avoid) * math.cos(r2)) + 
                            (math.sin(moon_avoid) * math.sin(r2) * math.cos(phase_angle_rad - theta2))) * RAD2DEG
        dist3 = math.acos((math.cos(moon_avoid) * math.cos(r3)) + 
                            (math.sin(moon_avoid) * math.sin(r3) * math.cos(phase_angle_rad - theta3))) * RAD2DEG
        
        # 条件に基づくフラグ設定
        if ((dist1 > dist1_l and dist1 < dist1_u) or
            (dist2 > dist2_l and dist2 < dist2_u) or
            (dist3 > dist3_l and dist3 < dist3_u)):
            # obsflag[sample] = False
            test.append(1)
        elif moon_avoid < moon27 and moon_avoid > moon23:
            if phase_angle_rad > phase27 and phase_angle_rad < phase63:
                # obsflag[sample] = False
                test.append(1)
            else:
                test.append(0)
        else:
            test.append(0)
    else:
        test.append(0)
    return test


def save_to_file(file_path, month, suffix, direc, data):  # ファイルの保存を行う関数
    result_file_name = file_path.replace('_gb.fits.gz', f'_{suffix}.pkl').replace(
        f'akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/{month}', f'akari_TSD_pkl/{direc}/all')
    with open(result_file_name, 'wb') as p:  # suffixだけではなくて、ディレクトリとファイル名を別に選択させる
        pickle.dump(data, p)





def make_moon(files):
    # month = "2006_04"
    # dir_path = f"/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/{month}/FIS_SW*"
    # file_path = glob.glob(dir_path)
    filename = "/mn/stornext/d23/cmbco/cg/AKARI/code/FISimage/diffuse/trunk/src/atFunctions/2.8/data/orbit.data"
    file_path = os.path.join(path_add, files)
    data = fitsio.FITS(file_path)
    d = data[1].read()
    dd = data[7].read()
    aftime = d['AFTIME']
    ra = dd['RA']
    dec = dd['DEC']
    NEP = SkyCoord(270, 66.560722222, unit=(u.deg, u.deg), frame='icrs')  # 北天
    SEP = SkyCoord(90, -66.560722222, unit=(u.deg, u.deg), frame='icrs')  # 南天

    EARTH_RADIUS = 6378.137      # 地球半径 [km]
    MOON_RADIUS = 1737.4         # 月半径 [km]
    TWO_PI = 2 * math.pi
    DEG2RAD = math.pi / 180.0
    RAD2DEG = 180.0 / math.pi

    atElement = {}

    image_param = []
    for i in range(len(ra)-1):
        obs_coord_a = SkyCoord(ra[i+1], dec[i+1], unit=(u.deg, u.deg), frame='icrs')
        obs_coord_b = SkyCoord(ra[i], dec[i], unit=(u.deg, u.deg), frame='icrs')

        angle_to_nep_b = (obs_coord_b.separation(NEP))  # NEPと観測データとの角度を計算
        angle_to_nep_a = (obs_coord_a.separation(NEP))

        if angle_to_nep_b > angle_to_nep_a:
            image_param.append(1)  # 角度がNEPに近い場合は1（NEP方向観測）、SEPに近い場合は0（SEP方向観測）
        else:
            image_param.append(0)
    a = image_param[-1]
    image_param.append(a)

    # t0 = Time('2000-01-01T00:00:00', scale='utc')  # 基準時刻（J2000.0）
    t0 = 51544.5
    MJD_J2000 = t0

    obsflag = [0] * len(ra) # 初期化

    mjd = []
    for i in range(len(aftime)):
        result_mjd = t0 + aftime[i] / 86400.0  # オフセット秒を日数に変換してMJDに加算
        mjd.append(result_mjd)

    for j in range(len(ra)):
        pos, size, phase, dist = at_moon(mjd[j], EARTH_RADIUS, MOON_RADIUS, DEG2RAD, TWO_PI)

        mjd0 = t0  # 元の MJD
        mjd1 = mjd[j]  # 変換先の MJD
        x0 = pos  # サンプルベクトル
        x1 = at_precession(mjd0, x0, mjd1, MJD_J2000, DEG2RAD)

        kchk = 0
        atElement = at_set_element(filename, mjd0, kchk, atElement, DEG2RAD, EARTH_RADIUS, TWO_PI)

        sat_xyz = at_sat_pos(mjd[j], atElement)

        moonpos = []
        for i in range(len(sat_xyz)):
            moonpos.append((-1 * sat_xyz[i]) + pos[i])

        moon = at_vect_to_pol(moonpos, TWO_PI)
        moon_eq = [moon['lon'], moon['lat']]

        m_ecliptic = equatorial2ecliptic(MJD_J2000, moon_eq)  # 実際の月の位置

        radec = [ra[j], dec[j]]  # ボアサイトの方向
        ecliptic = equatorial2ecliptic(MJD_J2000, radec)
        moon_avoid = angle(ecliptic, m_ecliptic)

        flag_test = calculate_obsflag(moon_avoid, m_ecliptic, ecliptic, obsflag, j, image_param[j], DEG2RAD, RAD2DEG)
        obsflag[j] = flag_test[0]

    lune = []
    for i in range(100):  # If LW: 75, elif SW: 100
        lune.append(obsflag)  # 0: survey/maneuver 1:pointing
    save_to_file(file_path, month, 'near_moon', 'flag', lune)

    end = time.time()
    time_diff = end - start
    print(time_diff, files)

if __name__ == '__main__':
    directory = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/'
    months = ['2006_04', '2006_05','2006_06','2006_07','2006_08','2006_09','2006_10','2006_11','2006_12','2007_01','2007_02','2007_03','2007_04','2007_05','2007_06','2007_07','2007_08']
    for month in months:
        path_add = str(directory) + str(month)
        files = [f for f in os.listdir(path_add) if f.startswith('FIS_SW_')]  # LW or SW

        with Pool(processes=100) as pool:
            pool.map(make_moon, files)



# # # test # # 
# month = "2006_04"
# dir_path = f"/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD/www.ir.isas.jaxa.jp/~yamamura/DR2/{month}/FIS_LW*"
# file_path = glob.glob(dir_path)
# filename = "/mn/stornext/d23/cmbco/cg/AKARI/code/FISimage/diffuse/trunk/src/atFunctions/2.8/data/orbit.data"

# nside = 256
# n_pix = hp.nside2npix(nside)  # Healpix ピクセルを取得

# NEP = SkyCoord(270, 66.560722222, unit=(u.deg, u.deg), frame='icrs')  # 北天
# SEP = SkyCoord(90, -66.560722222, unit=(u.deg, u.deg), frame='icrs')  # 南天

# EARTH_RADIUS = 6378.137      # 地球半径 [km]
# MOON_RADIUS = 1737.4         # 月半径 [km]
# TWO_PI = 2 * math.pi
# DEG2RAD = math.pi / 180.0
# RAD2DEG = 180.0 / math.pi

# atElement = {}

# lon = []
# lat = []
# for pix in range(n_pix):
#     theta, phi = hp.pix2ang(nside, pix)  # θ (天頂角), φ (経度)
#     lon.append(theta - np.pi/2)
#     lat.append(phi)

# moon_avoid = 25 * DEG2RAD
# image_param = 0

# # 角度の定義
# moon20 = 21.5 * DEG2RAD
# moon23 = 23.0 * DEG2RAD
# moon27 = 26.5 * DEG2RAD
# moon34 = 34.0 * DEG2RAD
# phase27 = 27.0 * DEG2RAD
# phase63 = 63.0 * DEG2RAD
# m_ecliptic = [0,0]

# test = []
# obsflag = [0] * len(range(n_pix))
# for j in range(len(lon)):
#     ecliptic = [lat[j], lon[j]]
#     m_ecliptic = [0,0]
#     if moon_avoid < 21.5 * DEG2RAD:
#         distance = math.sqrt(lat[j]**2 + lon[j]**2)
#         distance1 = math.sqrt((2*np.pi-lat[j])**2 + lon[j]**2)
#         if distance < (21.5 * DEG2RAD):
#             pass
#         elif distance1 < (21.5 * DEG2RAD):
#             pass
#         else:
#             continue
#     elif (21.5 * DEG2RAD) < moon_avoid < (34 * DEG2RAD):
#         distance = math.sqrt(lat[j]**2 + lon[j]**2)
#         distance1 = math.sqrt((2*np.pi-lat[j])**2 + lon[j]**2)
#         if 21.5 * DEG2RAD < distance < 34 * DEG2RAD:
#             pass
#         elif 21.5 * DEG2RAD < distance1 < 34 * DEG2RAD:
#             pass
#         else:
#             continue

#     sample = j
#     flag_test = calculate_obsflag(moon_avoid, m_ecliptic, ecliptic, obsflag, j, image_param)
#     obsflag[j] = flag_test[0]

# moon_avoid = 2 * DEG2RAD
# test = [0] * len(range(n_pix))
# for j in range(len(lon)):
#     ecliptic = [lat[j], lon[j]]
#     m_ecliptic = [0,0]
#     if moon_avoid < 21.5 * DEG2RAD:
#         distance = math.sqrt(lat[j]**2 + lon[j]**2)
#         distance1 = math.sqrt((2*np.pi-lat[j])**2 + lon[j]**2)
#         if distance < (21.5 * DEG2RAD):
#             pass
#         elif distance1 < (21.5 * DEG2RAD):
#             pass
#         else:
#             continue
#     elif (21.5 * DEG2RAD) < moon_avoid < (34 * DEG2RAD):
#         distance = math.sqrt(lat[j]**2 + lon[j]**2)
#         distance1 = math.sqrt((2*np.pi-lat[j])**2 + lon[j]**2)
#         if 21.5 * DEG2RAD < distance < 34 * DEG2RAD:
#             pass
#         elif 21.5 * DEG2RAD < distance1 < 34 * DEG2RAD:
#             pass
#         else:
#             continue

#     sample = j
#     flag_test = calculate_obsflag(moon_avoid, m_ecliptic, ecliptic, test, j, image_param)
#     test[j] = flag_test[0]

# tt = np.array(test) + np.array(obsflag)
# hp.write_map("test3.fits", tt, overwrite=True)
# hp.mollview(tt,flip="geo")
# hp.graticule()
# plt.savefig('test.png')
