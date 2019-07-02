import pandas as pd

def readEfficiently(tree, **kwargs):
  dfs = list(tree.iterate(**kwargs))
  return pd.concat(dfs)

def convert2hdf5(dfs, filename, treenames, userFunctions=[], data_columns=True, complib='blosc:lz4', complevel=1):
  from pandas import HDFStore
  store = HDFStore(filename, mode='w')
  if data_columns == True:
    data_columns = [True] * len(treenames)
  elif len(data_columns) != len(treenames or len(userFunctions)):
    raise ValueError
  for df in dfs:
    for _treename, _function, _data_columns in zip(treenames, userFunctions, data_columns):
      df_modified = _function(df)
      if df_modified is not None:
        store.put(_treename, df_modified, append=True, format='table', complib=complib, complevel=complevel, data_columns=_data_columns)
  store.close()
