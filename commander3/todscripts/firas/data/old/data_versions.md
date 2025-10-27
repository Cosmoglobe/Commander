| Version | Description                                                                                                                                                                                                                                                                                     |
| ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| v1      | Test for relaxing constraints so that we don't have to have all channels available, and only need to match to engineering data. Constraint on having mtm_speed and mtm_length the same on all 4 channels in science and enginnering data. Only saving id, time, ifgs, mtm_length and mtm_speed. |
| v2      | Test for matching the ifgs with each other over their own `gmt` variable. Only has IFGs and gmt.                                                                                                                                                                                                |
| v3      | Matching the IFGs with each other over their own `gmt` variable. Only has IFGs, gmt and xcal_pos. Constrained over xcal_pos. Fixing the delay between left/right channels.                                                                                                                      |
| v4      | Making GMTs into datetimes before merging to fix 99 problem.                                                                                                                                                                                                                                    |
| v5      | Getting rid of anything less than seconds to match the IFGs over GMT.                                                                                                                                                                                                                           |
| v6      | Interpolating over ICAL temperatures with the constraint of GMT time being within one minute each way. Contains: GMT, IFGs, xcal_pos, (interpolated) ICAL temperatures.                                                                                                                         |
| v7      | Added MTM length and speed to the data.                                                                                                                                                                                                                                                         |
| v8      | Changed tolerance to 16 seconds and added XCAL back in.                                                                                                                                                                                                                                         |
| v9      | Added fake-it mode flag to the data.                                                                                                                                                                                                                                                            |
| v10     | Clean by upmode = 4.                                                                                                                                                                                                                                                                            |
| v11     | Changed the way to pick high/low GRT temps.                                                                                                                                                                                                                                                     |
| v12     | Added bolometer voltages and commaded voltages and binary time.                                                                                                                                                                                                                                 |
| v13     | Added sci_gain.                                                                                                                                                                                                                                                                                 |
| v14     | Added sweeps and switched out sci_gain for gain.                                                                                                                                                                                                                                                |
| v15     | Added galactic latitude and longitude.                                                                                                                                                                                                                                                          |
| v16     | Removed galactic latitude and longitude. Changing script to make sure it's possible to have records without all of the channels. Added rejects for sun angle, earth limb angle and moon angle.                                                                                                  |

Divided the data into calibration and sky data:

| Version | Description                                          |
| ------- | ---------------------------------------------------- |
| v1      | New dataset with sky and calibration data separated. |
| v2      | Added a lot of new flags.                            |


`data_prep_v4.py`

| Version | Description                                                            |
| ------- | ---------------------------------------------------------------------- |
| v4.1    | Copied over from other data_prep script but without interpolating now. |
| v4.2    | Changed possible values for sweeps (took out 1).                       |