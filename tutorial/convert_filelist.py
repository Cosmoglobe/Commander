import os

filelists = [os.getcwd() + "/auxcmd3/filelist_30_cm3.txt",
             os.getcwd() + "/auxcmd3/filelist_44_cm3.txt",
             os.getcwd() + "/auxcmd3/filelist_70_cm3.txt"]

data_path = os.getcwd() + "/L2Data/"

for filelist in filelists:
    with open(filelist, "r") as my_filelist:
        filedata = my_filelist.read()

    # Replace the target string
    filedata = filedata.replace("input/l2data/", data_path)

    # Write the file out again
    with open(filelist, 'w') as my_filelist:
        my_filelist.write(filedata)

    print(filedata)    
