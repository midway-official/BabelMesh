#!/usr/bin/env python3
"""Diagnostic only: oriented polygon sums, edge incidence, analytic 30deg box.

Area vectors use the polygon cross-product identity, not a triangulation-based
self-intersection test. The orientation prediction uses the reader's initial
face-incidence-weighted point mean; it is not a substitute for the real reader.
"""
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def mean(points):
    return tuple(sum(p[k] for p in points)/len(points) for k in range(3))


def audit(path):
    with path.open() as f:
        assert f.readline().strip() == 'BABELSIM_MESH 3'
        points = [tuple(map(float, f.readline().split())) for _ in range(int(f.readline().split()[1]))]
        faces = []
        cells = collections.defaultdict(list)
        for i in range(int(f.readline().split()[1])):
            row = f.readline().split()
            n, owner, neighbour, patch = map(int, row[1:5])
            ring = list(map(int, row[5:]))
            assert len(ring) == n
            p = [points[v] for v in ring]
            terms = [cross(sub(p[j], p[0]), sub(p[j+1], p[0])) for j in range(1, n-1)]
            area = tuple(0.5*sum(t[k] for t in terms) for k in range(3))
            faces.append((owner, neighbour, patch, ring, mean(p), area))
            cells[owner].append((i, 1))
            if neighbour >= 0:
                cells[neighbour].append((i, -1))
        patches = [f.readline().split()[1] for _ in range(int(f.readline().split()[1]))]
    centers = {c: mean([points[v] for i, _ in fs for v in faces[i][3]]) for c, fs in cells.items()}
    flipped = [dot(a, sub(center, centers[o])) < 0 for o, _, _, _, center, a in faces]
    bad_edges, bad_raw, bad_reoriented, bad_halfspace = [], [], [], []
    for c, fs in cells.items():
        edges = collections.defaultdict(lambda: [0, 0])
        raw, oriented = [0.0]*3, [0.0]*3
        total = 0.0
        inward = []
        for i, sign in fs:
            _, _, _, ring, center, area = faces[i]
            total += math.sqrt(dot(area, area))
            if sign*dot(area, sub(center, centers[c])) <= 0:
                inward.append(i)
            for k in range(3):
                raw[k] += sign*area[k]
                oriented[k] += sign*area[k]*(-1 if flipped[i] else 1)
            for a, b in zip(ring, ring[1:]+ring[:1]):
                e = edges[tuple(sorted((a, b)))]
                e[0] += 1
                e[1] += sign*(1 if a < b else -1)
        if any(n != 2 or s != 0 for n, s in edges.values()):
            bad_edges.append(c)
        if math.sqrt(dot(raw, raw)) > 1e-9*max(total, 1):
            bad_raw.append(c)
        if math.sqrt(dot(oriented, oriented)) > 1e-9*max(total, 1):
            bad_reoriented.append(c)
        if inward:
            bad_halfspace.append(c)
    # Solid is [0.5,1.5]^3 rotated about (1,1,1) by +30deg in xy.
    co, si = math.cos(math.pi/6), math.sin(math.pi/6)
    axes = ((co, si, 0), (-si, co, 0), (0, 0, 1))
    orthogonality, skewness = [], []
    orth_entities, skew_entities = [], []
    for face_index,(owner, neighbour, patch, ring, fc, area) in enumerate(faces):
        co=centers[owner]
        d=sub(centers[neighbour],co) if neighbour>=0 else sub(fc,co)
        da=math.sqrt(dot(area,area)); dd=math.sqrt(dot(d,d))
        if da>0 and dd>0:
            angle=math.degrees(math.acos(max(-1,min(1,dot(area,d)/(da*dd)))))
            orthogonality.append(angle); orth_entities.append((angle,face_index,owner,neighbour))
            if neighbour>=0:
                denom=dot(area,d)
                if abs(denom)>1e-30:
                    t=dot(area,sub(fc,co))/denom
                    intersection=tuple(co[k]+t*d[k] for k in range(3))
                    value=math.sqrt(dot(sub(fc,intersection),sub(fc,intersection)))/math.sqrt(da)
                    skewness.append(value); skew_entities.append((value,face_index,owner,neighbour))
                else:
                    skewness.append(float('inf')); skew_entities.append((float('inf'),face_index,owner,neighbour))
    orthogonality.sort(); skewness.sort()
    p99=lambda xs: xs[min(len(xs)-1,math.ceil(.99*len(xs))-1)] if xs else 0
    wrong_wall, off_wall, exterior_owner = [], [], []
    cell_vertices = {c: {v for i, _ in fs for v in faces[i][3]} for c, fs in cells.items()}
    for i, (owner, n, patch, ring, fc, area) in enumerate(faces):
        if n >= 0 or patches[patch] != 'cad_boundary':
            continue
        q = [tuple(dot(sub(points[v], (1, 1, 1)), ax) for ax in axes) for v in ring]
        choices = [(max(abs(p[k]-s*.5) for p in q), k, s) for k in range(3) for s in (-1, 1)]
        error, k, s = min(choices)
        if error > 1e-8 or any(max(map(abs, p)) > .5+1e-8 for p in q):
            off_wall.append(i)
        elif dot(area, axes[k])*(-s) <= 0:
            wrong_wall.append(i)
        else:
            outward = tuple(-s*x for x in axes[k])
            sides = [dot(sub(points[v], fc), outward) for v in cell_vertices[owner]]
            if min(sides) >= -1e-10 and max(sides) > 1e-10:
                exterior_owner.append(i)
    return dict(mesh=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                cells=len(cells), faces=len(faces), vertices=len(points),
                raw_edge_incidence_failures=len(bad_edges), raw_area_closure_failures=len(bad_raw),
                predicted_reader_flips=sum(flipped), predicted_reader_closure_failures=len(bad_reoriented),
                point_mean_halfspace_failures=len(bad_halfspace),
                analytic_wall_orientation_failures=len(wrong_wall), analytic_wall_membership_failures=len(off_wall),
                wall_faces_with_owner_entirely_on_exterior_side=len(exterior_owner),
                max_non_orthogonality_degrees=orthogonality[-1] if orthogonality else 0,
                p99_non_orthogonality_degrees= p99(orthogonality),
                max_internal_skewness=skewness[-1] if skewness else 0,
                p99_internal_skewness=p99(skewness),
                worst_non_orthogonality=max(orth_entities) if orth_entities else None,
                worst_skewness=max(skew_entities) if skew_entities else None,
                samples=dict(edge_cells=bad_edges[:16], closure_cells=bad_reoriented[:16],
                             wrong_wall_faces=wrong_wall[:16], off_wall_faces=off_wall[:16],
                             exterior_owner_wall_faces=exterior_owner[:16]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.dumps(audit(args.mesh), indent=2)+'\n'
    if args.output:
        args.output.write_text(result)
    print(result, end='')
