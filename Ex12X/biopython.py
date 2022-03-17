from Bio.SeqUtils.ProtParam import ProteinAnalysis as protAnal
import pandas as pd
from datetime import datetime

def proteinInfos(serie, save="n", outputName="proteinInfo"):
    """ Get more information about a protein Sequence:
    -> Isoelectric point
    -> Molecular weight
    -> Aromaticity
    -> Hydrophobicity (gravy)"""
#calculate specific information
    for i, entry in enumerate(serie.fastaSequence):
        if len(entry) < 2:
            pass
        else:
            print(f"Process: {i} / {len(serie.fastaSequence)}", end='\r')
            #replace uracil with thymidin due to error with gravy
            entry = entry.replace("U", "T")
            serie.loc[i, ["pI"]] = protAnal(entry).isoelectric_point()
            serie.loc[i, ["chargeAtpH8"]] = protAnal(entry).charge_at_pH(8)
            serie.loc[i, ["molecularWeight"]] = protAnal(entry).molecular_weight()
            serie.loc[i, ["aromaticity"]] = protAnal(entry).aromaticity()
            serie.loc[i, ["hydrophobicity"]] = protAnal(entry).gravy()

    if save == "y":
        ti = datetime.now().strftime("%Y%m%d-%H%M%S")
        serie.to_csv(ti + "_" + "proteinSeqInfo" + ".csv")

    return serie
