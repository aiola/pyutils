def readEfficiently(tree, **kwargs):
  dfs = list(tree.iterate(**kwargs))
  return pd.concat(dfs)

def convert2hdf5(treename, filename, dfs):
  from pandas import HDFStore
  store = HDFStore(filename, mode='w')
  for df in dfs:
    store.put(treename, df, append=True, format='table', complib='blosc:lz4', complevel=1, data_columns=True)
  store.close()
