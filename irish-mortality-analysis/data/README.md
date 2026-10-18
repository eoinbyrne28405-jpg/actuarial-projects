# Data

`irish_life_table_17.csv`: Irish Life Tables No. 17 (2015–2017), published by the Central Statistics Office (CSO), Ireland.

| Column | Meaning |
|---|---|
| `age` | Exact age x, 0–105 |
| `qx_male`, `qx_female` | Probability that a person aged exactly x dies before x+1 |
| `lx_male`, `lx_female` | Survivors to age x out of 100,000 births |

Source: CSO Irish Life Tables No. 17, 2015–2017, single-year table as published in the Human Life-Table Database: https://www.lifetable.de/File/GetDocument/data%5CIRL%5CIRL000020152017CU1.pdf

The values were transcribed from that document. `R/01_life_table.R` checks that every row satisfies l(x+1) = l(x)(1 − q(x)) to within rounding, and the resulting life expectancies at birth (79.6 male, 83.4 female) match the CSO's published figures. CSO material is reusable with attribution under its copyright and re-use policy.
