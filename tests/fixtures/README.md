# STEP fixture provenance

`gmsh-box-1000mm.step` is an independent CAD-import fixture created with Gmsh 4.8.4 using its OpenCASCADE geometry factory. The complete source was:

```geo
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1000, 1000, 1000};
Save "/tmp/babelmesh-gmsh-box.step";
```

It was exported with `gmsh -0 /tmp/babelmesh-gmsh-fixture.geo -v 2`; its STEP header declares `SI_UNIT(.MILLI.,.METRE.)`, and its expected normalized size is a 1 m cube with volume 1 m³. SHA-256:

```text
15f5afb33e0ff9b8326dbe0f1fdf3eee220ecfccc3abe42592409289524e7014  gmsh-box-1000mm.step
```

The fixture contains only the elementary box geometry needed to test STEP unit conversion. It is not a BabelMesh export and is retained here to make M1 acceptance repeatable.

## Tetrahedral topology fixtures

`topology-three-tet-star.json` and `topology-four-tet-star.json` contain small,
closed triangular/quadrilateral bipyramids with two apex points and an internal
apex-to-apex edge. Their coordinates and positively oriented connectivity were
constructed for this repository and then fixed as regression inputs. They are
independent of the production mesher and contain no third-party model.

The first fixture exercises a 3-to-2 swap that reduces cell count and skewness;
the second exercises a 4-to-4 swap that preserves cell count and passes the
previously failing face-quality gate. The C++ contract test also checks unchanged
points, cavity boundary triangles and patch IDs, preserved volume, all four
maximum/p99 FVM metrics, layer protection, and actual BabelSim reader validation.
