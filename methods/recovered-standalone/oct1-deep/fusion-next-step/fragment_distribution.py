"""Describe kallisto fragment histogram; primary Pizzly threshold follows pinned helper."""
import json, pathlib, sys
import h5py
import numpy as np
p=pathlib.Path(sys.argv[1]).resolve()
with h5py.File(p,'r') as f:
    x=np.asarray(f['aux/fld'],dtype='float64')
if not np.isfinite(x).all() or x.sum()<=0:raise SystemExit('Invalid/empty fragment histogram')
c=np.cumsum(x)/x.sum()
q={str(z):int(np.argmax(c>z)) for z in [.5,.9,.95,.99]}
r={'source':str(p),'histogram_count':float(x.sum()),'nonzero_min':int(np.flatnonzero(x)[0]),'nonzero_max':int(np.flatnonzero(x)[-1]),'mean':float(np.dot(np.arange(len(x)),x)/x.sum()),'quantiles_strictly_greater_cumulative':q,'primary_pizzly_insert_size':q['0.95'],'method':'Uses np.argmax(cumsum(hist)/sum(hist) > 0.95), matching Pizzly v0.37.3 scripts/get_fragment_length.py; this is a quantification-derived fragment estimate, not direct molecule sizing.','h5py_version':h5py.__version__,'numpy_version':np.__version__,'histogram':x.astype(int).tolist()}
out=p.parent/'fragment-distribution.json';out.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:v for k,v in r.items() if k!='histogram'},indent=2))
