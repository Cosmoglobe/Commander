.SUFFIXES
.SUFFIXES  .MEM .MEC .MEX .RNT .RNX .BRN .RNO

.RNO.BRN
 $(RUNOFF)/INTERMEDIATE/NOOUTPUT $(RFLAGS) $(MMS$SOURCE)

DBASE_MANUAL : FIRAS.MEM, FIRAS.BRN
 ! Database Manual Created

FIRAS.MEM, FIRAS.BRN : FIRAS.RNO, -
                       FEX_CTH.RNO, -
                       FEX_DTF.RNO, -
                       FEX_ETF.RNO, -
                       SETUP.RNO
