"""Independent VTK_HEXAHEDRON corner Jacobians using explicit edge pairs."""
import argparse
import json
import xml.etree.ElementTree as ET
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('vtu'); p.add_argument('--output')
p.add_argument('--bsm',help='Optional matching BSM for coplanar wall-face diagnosis')
args=p.parse_args()
piece=ET.parse(args.vtu).getroot().find('./UnstructuredGrid/Piece')
points=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(-1,3)
arrays={a.get('Name'):np.fromstring(a.text,sep=' ',dtype=np.int64) for a in piece.findall('./Cells/DataArray')}
hex_ids=np.flatnonzero(arrays['types']==12)
ends=arrays['offsets'][hex_ids]
indices=np.stack([arrays['connectivity'][e-8:e] for e in ends])
corners=points[indices]
xp=[(0,1),(0,1),(3,2),(3,2),(4,5),(4,5),(7,6),(7,6)]
yp=[(0,3),(1,2),(1,2),(0,3),(4,7),(5,6),(5,6),(4,7)]
zp=[(0,4),(1,5),(2,6),(3,7)]*2
values=[]
for x,y,z in zip(xp,yp,zp):
    u=corners[:,x[1]]-corners[:,x[0]]
    v=corners[:,y[1]]-corners[:,y[0]]
    w=corners[:,z[1]]-corners[:,z[0]]
    denominator=np.linalg.norm(u,axis=1)*np.linalg.norm(v,axis=1)*np.linalg.norm(w,axis=1)
    values.append(np.divide(np.einsum('ij,ij->i',u,np.cross(v,w)),denominator,
                            out=np.full(len(hex_ids),-1.),where=denominator>0))
values=np.array(values).T
mins=values.min(axis=1)
worst=np.argsort(mins)[:8]
result=dict(vtu=args.vtu,hex_cells=len(hex_ids),minimum=float(mins.min()),
    cells_below_point_one=int((mins<.1).sum()),
    cells_negative_beyond_roundoff=int((mins<-1e-10).sum()),
    cells_near_zero=int((abs(mins)<=1e-10).sum()),
    negative_tolerance=1e-10,
    worst=[dict(cell=int(hex_ids[i]),minimum=float(mins[i]),corner_jacobians=values[i].tolist(),
                point_ids=indices[i].tolist(),coordinates=corners[i].tolist()) for i in worst])
if args.bsm:
    boundary={}
    with open(args.bsm) as f:
        f.readline(); n=int(f.readline().split()[1])
        for _ in range(n): f.readline()
        for _ in range(int(f.readline().split()[1])):
            row=f.readline().split(); owner,neighbour=map(int,row[2:4])
            if neighbour>=0: continue
            ids=np.array(list(map(int,row[5:])))
            xyz=points[ids]; normal=np.cross(xyz-xyz[0],np.roll(xyz,-1,axis=0)-xyz[0]).sum(axis=0)
            length=np.linalg.norm(normal)
            if length>0: boundary.setdefault(owner,[]).append((normal/length,xyz[0],set(ids.tolist())))
    matches=[]
    for i in np.flatnonzero(mins<.1):
        cell=int(hex_ids[i]); fs=boundary.get(cell,[])
        if any(np.dot(a[0],b[0])>1-1e-10 and abs(np.dot(a[0],b[1]-a[1]))<1e-10
               and len(a[2]&b[2])>=2 for j,a in enumerate(fs) for b in fs[j+1:]):
            matches.append(cell)
    result['bad_hex_with_adjacent_coplanar_wall_faces']=len(matches)
    result['coplanar_wall_examples']=matches[:12]
text=json.dumps(result,indent=2)+'\n'
if args.output: open(args.output,'w').write(text)
print(text)
