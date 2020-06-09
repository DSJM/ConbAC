import csv
import gzip
import os
import pandas as pd
import scipy.io


def importMM(matrix_dir):
  mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
  mat = pd.DataFrame.sparse.from_spmatrix(mat)
  features_path = os.path.join(matrix_dir, "genes.tsv")
  features = [row[0] for row in csv.reader(open(features_path), delimiter="\t")]
  barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
  barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
  mat.index = features
  mat.columns = barcodes
  
  return mat

