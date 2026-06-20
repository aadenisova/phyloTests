#!/usr/bin/env python3
from pathlib import Path

import pandas as pd
from ete3 import Tree

DATA = Path(__file__).resolve().parent / "data"

ann = pd.read_csv(DATA / "annotations.tsv", sep="\t")
birds = ann.loc[ann["Lineage"] == "Birds", ["# accession", "ScientificName"]].copy()
birds["Species"] = birds["ScientificName"].str.replace(" ", "_", regex=False)

traits = pd.read_csv(DATA / "ALLBIRDTRAITS_Database_March_2026.csv")
traits = traits[traits["ResearchEffort"] >= 25]

birds = birds[birds["Species"].isin(traits["Species"])]

acc2species = dict(zip(birds["# accession"], birds["Species"]))
species_keep = set(acc2species.values())

tree = Tree(str(DATA / "roadies_v1.1.16b.nwk"), format=1)
tree.prune(list(acc2species), preserve_branch_length=True)
for leaf in tree:
    leaf.name = acc2species[leaf.name]

traits_out = traits[traits["Species"].isin(species_keep)].copy()

tree.write(format=1, outfile=str(DATA / "roadies_birds_allbirdtraits.nwk"))
traits_out.to_csv(DATA / "ALLBIRDTRAITS_intersect.csv", index=False)

print(f"intersect: {len(species_keep)} species, {len(tree)} leaves, {len(traits_out)} trait rows")
