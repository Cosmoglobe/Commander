import numpy as np

channel = "rh"
mode = "ss"

data = np.load(f"{g.PREPROCESSED_DATA_PATH}cal_{channel}.npz")
mtm_length = data["mtm_length"][:]
mtm_speed = data["mtm_speed"][:]

ss_filter = (mtm_length == 0) & (mtm_speed == 0)

ifg = data["ifg"][ss_filter]
