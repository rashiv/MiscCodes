#!/usr/bin/env python2.7

################################### Averaging duplicates in column of 2D numpy array ########################
def avg_dupl(arr):
  arr = np.array(arr)
  a = arr[np.argsort(arr[:,0])] # sort by column 0
  a = np.unique(a[:,0])
  b = pd.DataFrame(arr).groupby(0).mean().values
  b = np.squeeze(b)
# reduce 1 dim since Pandas makes each entry a separate list
  c = [a,b]
  return c
# return format: [[column 0],[column 1]]
