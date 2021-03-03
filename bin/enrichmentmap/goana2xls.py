#!/usr/bin/env python

import sys
import pandas as pd


goana = sys.argv[1]
out = sys.argv[2]

goana = pd.read_csv(goana, sep="\t")
goana = goana[goana.Ont=="BP"]

d = {"GO.ID":[], 
     "Description":[],
     "p.Val":[],
     "FDR":[]}
d["GO.ID"] = goana.index.tolist()
d["Description"] = goana["Term"]
d["p.Val"] = goana["P.DE"]
d["FDR"] = "NA"

df = pd.DataFrame(d)

df.to_csv(out, index=False, encoding='ISO-8859-1', sep="\t", na_rep="NA")

