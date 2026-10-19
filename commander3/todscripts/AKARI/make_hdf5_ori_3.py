import numpy as np
import h5py as h5
import os
import glob

# hdf5ファイルを繋ぎ合わせる
def copy_datasets(src_file, dest_file, src_group_base="000001", dest_group_base="000001"):
    # コピー先のファイルが存在しない場合は新規作成
    if not os.path.exists(dest_file):
        with h5.File(dest_file, "w") as _:
            pass  # ファイルを新規作成のみ行う

    # ソースファイルが存在するか確認
    if not os.path.exists(src_file):
        print(f"指定したソースファイルが存在しません: {src_file}")
        return  # ファイルがなければスキップ

    with h5.File(src_file, "r") as hf1, h5.File(dest_file, "r+") as hf2:
        # ソースグループ内の全てのデータセットを走査
        def visit_and_copy(name, obj):
            if isinstance(obj, h5.Dataset):
                # データセットの名前を取得
                dest_dataset_path = name.replace(src_group_base, dest_group_base)
                
                # コピー先のデータセットが存在する場合は削除
                if dest_dataset_path in hf2:
                    del hf2[dest_dataset_path]
                
                # データセットをコピー
                h5.h5o.copy(hf1.id, name.encode(), hf2.id, dest_dataset_path.encode())
                print(f"コピー完了: {dest_dataset_path}")
            elif isinstance(obj, h5.Group):
                # グループが存在しない場合は作成
                dest_group_path = name.replace(src_group_base, dest_group_base)
                if dest_group_path not in hf2:
                    hf2.create_group(dest_group_path)

        # 指定したグループを訪問してコピー
        hf1.visititems(visit_and_copy)


# filelsit.txtの作成
def save_data_to_txt(n, data_path, file_path):
    # 出力するデータのリストを初期化
    data = []
    
    # データを生成
    for i in range(584, n + 1):
    # for i in range(1, n + 1):
        data_line = f'{i} {data_path} 1 0.0 0.0'
        # data_line = f'{i} "{src_dir}/merged_AKARI_{band}_n{nside_out}_v{version:02}_flux_reflag.h5" 1 0.0 0.0'
        # data_line = f'{i} "/mn/stornext/d23/cmbco/cg/AKARI/ryosukem/data/akari/merged_output_flux.h5" 1 0.0 0.0'
        data.append(data_line)

    # 最初の行に n を追加
    data.insert(0, 11395)
    # data.insert(0, str(n))
    
    # ファイルに書き込む
    with open(file_path, 'w') as f:
        for line in data:
            f.write(str(line) + '\n')

def make_hdf5(band, src_dir):
    nside_out = 2048
    version = 3

    # copy
    src_files = sorted(glob.glob(os.path.join(src_dir, f"data_{band}/200*_flux_reflag.h5")))
    dest_file_path = f"{src_dir}/merged_AKARI_{band}_n{nside_out}_v{version:02}_flux_reflag.h5" # output file
    # dest_file_path = "/mn/stornext/d23/cmbco/cg/AKARI/ryosukem/data/akari/merged_output_flux.h5"

    print(src_dir)

    # 最初のファイルを000001として保存
    first_file = src_files[0]
    copy_datasets(first_file, dest_file_path, src_group_base="000001", dest_group_base="000001")

    # 残りのファイルを連番で追加
    for i, src_file in enumerate(src_files[1:], start=2):
        src_group_base = "000001"  # 各ファイルのソースグループ名（固定の場合）
        dest_group_base = f"{i:06d}"  # コピー先の連番グループ名
        copy_datasets(src_file, dest_file_path, src_group_base, dest_group_base)

    print('Completed to copy all files')

    n = 11978
    data_path = f"{src_dir}/merged_AKARI_{band}_n{nside_out}_v{version:02}_flux_reflag.h5"
    file_path = f"{src_dir}/filelist_merged_AKARI_{band}_n{nside_out}_v{version:02}_flux_reflag.txt"
    # file_path = '/mn/stornext/d23/cmbco/cg/AKARI/ryosukem/data/akari/filelist_N160_startt20060413074458_endt20060816230000.txt'  # 出力ファイル名
    save_data_to_txt(n, data_path, file_path)

if __name__ == '__main__':
    nside_out = 2048
    version = 3
    
    band = 'N160'
    src_dir = f"/mn/stornext/d23/cmbco/cg/AKARI/tamakim3/data"
    # src_dir = "/mn/stornext/d23/cmbco/cg/AKARI/ryosukem/data/akari"
    make_hdf5(band, src_dir)
