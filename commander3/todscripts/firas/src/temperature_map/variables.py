import numpy as np

data = np.load(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/sky.npz",
    allow_pickle=True,
)

# print(data.files)

print(f"stat_word_16: {data["stat_word_16_lf"][137900]}")
print(f"stat_word_13: {data["stat_word_13_lf"][158400]}")
print(f"stat_word_5: {data["stat_word_5_lf"][1055]}")
print(f"stat_word_9: {data["stat_word_9_lf"][6253]}")
print(f"lvdt_stat_a: {data["lvdt_stat_a_lf"][93609]}")
print(f"lvdt_stat_b: {data["lvdt_stat_b_lf"][93609]}")
