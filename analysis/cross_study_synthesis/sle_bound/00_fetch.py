import datetime
import hashlib
import os
import sys

import cellxgene_census
import fsspec
import h5py
import numpy as np
import pandas as pd

CENSUS = "2025-11-08"
DATASET = "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"
H5AD = "https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad"
OUT = "cache/perez2022"
LINEAGES = {"classical monocyte": "cM", "natural killer cell": "NK"}

genes = pd.read_csv("analysis/inputs/core_isgs.tsv", sep="\t")["gene"].tolist()
os.makedirs(OUT, exist_ok=True)
t0 = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

with cellxgene_census.open_soma(census_version=CENSUS) as c:
    ad = cellxgene_census.get_anndata(
        c, "Homo sapiens", X_name="raw",
        obs_value_filter=f"dataset_id == '{DATASET}' and cell_type in {list(LINEAGES)}",
        var_value_filter=f"feature_name in {genes}",
        obs_column_names=["observation_joinid", "dataset_id", "donor_id", "disease", "self_reported_ethnicity",
                          "cell_type", "assay", "is_primary_data", "raw_sum", "nnz"],
    )
X = ad.X.toarray()
cells = ad.obs.reset_index(drop=True)
for j, g in enumerate(ad.var["feature_name"]):
    cells[g] = X[:, j]
cells["census_lineage"] = cells["cell_type"].map(LINEAGES)

f = fsspec.open(H5AD, block_size=2**20).open()
o = h5py.File(f, "r")["obs"]
dec = lambda a: [x.decode() if isinstance(x, bytes) else x for x in a]
meta = pd.DataFrame({"observation_joinid": dec(o["observation_joinid"][:])})
for k in ["author_cell_type", "ct_cov", "Processing_Cohort", "disease_state", "sample_uuid", "ind_cov", "library_uuid"]:
    meta[k] = pd.Categorical.from_codes(o[k]["codes"][:], dec(o[k]["categories"][:]))
n_author = meta["author_cell_type"].value_counts()
cells = cells.merge(meta, on="observation_joinid", how="left", validate="one_to_one")

missing = [g for g in genes if g not in cells.columns]
if missing:
    sys.exit(f"panel genes missing: {missing}")
cols = ["observation_joinid", "donor_id", "ind_cov", "sample_uuid", "library_uuid", "disease", "disease_state",
        "self_reported_ethnicity", "Processing_Cohort", "assay", "is_primary_data", "cell_type", "census_lineage",
        "author_cell_type", "ct_cov", "raw_sum", "nnz"] + genes
path = f"{OUT}/perez2022_cM_NK_panel_counts.csv.gz"
cells[cols].to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
pd.DataFrame({"author_cell_type": n_author.index.astype(str), "n_cells_h5ad": n_author.values}).to_csv(
    f"{OUT}/perez2022_author_cell_type_counts.csv", index=False)

t1 = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
md5 = hashlib.md5(open(path, "rb").read()).hexdigest()
rows = [
    ("source", "Perez et al. 2022 Science doi:10.1126/science.abf1970; GEO GSE174188; CELLxGENE collection 436154da-bcf1-4130-9c8b-120ff9a888f2"),
    ("census_version", CENSUS),
    ("census_uri", "s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/"),
    ("census_api", f"cellxgene_census {cellxgene_census.__version__}; X_name=raw; obs filter dataset_id and cell_type in {list(LINEAGES)}; var filter feature_name in panel"),
    ("dataset_id", DATASET),
    ("dataset_version_id", "c55dc602-d168-4d15-acc1-5de4f2f5d551"),
    ("h5ad_url", H5AD),
    ("h5ad_bytes_not_downloaded", "12218105530"),
    ("h5ad_access", "HTTP range reads of obs/{observation_joinid,author_cell_type,ct_cov,Processing_Cohort,disease_state,sample_uuid,ind_cov,library_uuid} only; joined on observation_joinid"),
    ("library_size", "Census obs raw_sum = per-cell sum of raw counts over all genes"),
    ("retrieval_start_utc", t0),
    ("retrieval_end_utc", t1),
    ("output_file", path),
    ("output_bytes", str(os.path.getsize(path))),
    ("output_md5", md5),
    ("output_n_cells", str(len(cells))),
]
pd.DataFrame(rows, columns=["key", "value"]).to_csv("analysis/cross_study_synthesis/sle_bound/data_source.tsv", sep="\t", index=False)
print(cells.groupby(["census_lineage", "author_cell_type"], observed=True).size())
