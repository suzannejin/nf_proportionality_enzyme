#!/usr/bin/env python

import sys
import pandas as pd


graflex = sys.argv[1]
out = sys.argv[2]

graflex = pd.read_csv(graflex, sep="\t")
graflex = graflex[graflex.Ontology=="BP"]

d = {"GO.ID":[], 
     "Description":[],
     "p.Val":[],
     "FDR":[]}
d["GO.ID"] = graflex["Concept"]
d["Description"] = graflex["Term"]
d["p.Val"] = graflex["FDR.over"]
d["FDR"] = "NA"

df = pd.DataFrame(d)

df.to_csv(out, index=False, encoding='ISO-8859-1', sep="\t", na_rep="NA")

