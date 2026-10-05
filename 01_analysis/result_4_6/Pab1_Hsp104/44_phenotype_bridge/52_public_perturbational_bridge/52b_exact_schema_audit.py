#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path("/")
OUT = ROOT / "44_phenotype_bridge/52_public_perturbational_bridge"

SHATTUCK = OUT / "Shattuck2019_Source_Data.xlsx"
S6 = OUT / "Andersson2021_S6_Guk1_HSP104.xlsx"
S10 = OUT / "Andersson2021_S10_Guk1_timelapse.xlsx"


def show_block(path, sheet, r0, r1, c0, c1):
    x = pd.read_excel(
        path,
        sheet_name=sheet,
        header=None,
        engine="openpyxl",
    )

    print()
    print("=" * 120)
    print(path.name, "::", sheet)
    print(f"ROWS {r0}:{r1-1}  COLS {c0}:{c1-1}")
    print("=" * 120)

    z = x.iloc[r0:r1, c0:c1].copy()

                                                                        
    z.columns = [
        f"C{j}"
        for j in range(c0, c1)
    ]

    z.index = [
        f"R{i}"
        for i in range(r0, r1)
    ]

    print(
        z.to_string(
            max_rows=None,
            max_cols=None,
        )
    )


                                                              
                   
 
                                                    
                                                              

show_block(
    SHATTUCK,
    "Figure 6",
    0, 39,
    0, 19,
)


                                                              
              
 
                                                      
                                                              

xl = pd.ExcelFile(
    S6,
    engine="openpyxl",
)

print()
print("=" * 120)
print("S6 SHEETS")
print("=" * 120)
print(xl.sheet_names)

combined = pd.read_excel(
    S6,
    sheet_name="Combined",
    header=None,
    engine="openpyxl",
)

print()
print("=" * 120)
print("ANDERSSON S6 :: COMBINED")
print("shape:", combined.shape)
print("=" * 120)

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.width", 1000)

print(
    combined.to_string(
        header=True,
        index=True,
    )
)


                                                              
                  
                                                           
                                                         
                                
                                                              

for sheet in [
    "190423",
    "190424",
    "190425",
]:

    x = pd.read_excel(
        S6,
        sheet_name=sheet,
        header=None,
        engine="openpyxl",
    )

    print()
    print("=" * 120)
    print("ANDERSSON S6 ::", sheet, ":: FIRST 10 COLS")
    print("=" * 120)

    print(
        x.iloc[:, :10].to_string(
            header=True,
            index=True,
        )
    )


                                                              
      
                                               
                                                              

xl10 = pd.ExcelFile(
    S10,
    engine="openpyxl",
)

print()
print("=" * 120)
print("S10 SHEETS")
print("=" * 120)

for s in xl10.sheet_names:
    x = pd.read_excel(
        S10,
        sheet_name=s,
        header=None,
        engine="openpyxl",
    )

    print(
        s,
        "shape=",
        x.shape,
    )


                                                              
                                                 
                                                              

fig6 = pd.read_excel(
    SHATTUCK,
    sheet_name="Figure 6",
    header=None,
    engine="openpyxl",
)

fig6.to_csv(
    OUT / "52b_Shattuck_Figure6_raw.tsv",
    sep="\t",
    index=True,
    header=True,
)

combined.to_csv(
    OUT / "52b_Andersson_S6_Combined_raw.tsv",
    sep="\t",
    index=True,
    header=True,
)

print()
print("=" * 120)
print("WRITTEN")
print("=" * 120)

print(OUT / "52b_Shattuck_Figure6_raw.tsv")
print(OUT / "52b_Andersson_S6_Combined_raw.tsv")
