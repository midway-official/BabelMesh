# Third-party dependency inventory

This inventory records the dependency versions in the current development
environment and their upstream license notices. It is not a substitute for a
release-specific software bill of materials. Recheck package versions and
license files before redistributing binaries.

The repository's application source is MIT-licensed. The linked hybrid
application also uses Gmsh under GPL-2.0-or-later with its stated exception;
MIT for BabelMesh source alone does not describe the linked binary
distribution terms. See the dependency notices for the actual build. OCCT is a separate shared
library under LGPL 2.1 with its Open CASCADE exception; binary releases still
need the associated notices and source/relinking terms.

The tested Ubuntu Jammy package versions are frozen in
[`dependencies/ubuntu-jammy.lock.json`](dependencies/ubuntu-jammy.lock.json),
and `make deps` checks them against the installed system packages before running
the compile/link probe.

| Component | Role | Tested / planned version | License |
|---|---|---:|---|
| Open CASCADE Technology | CAD kernel, Boolean operations, STEP | 7.5.1 (Ubuntu Jammy package) | LGPL-2.1-only WITH OCCT-exception-1.0 |
| Eigen | Sparse/dense numerical support | 3.4.0 (Ubuntu Jammy package) | MPL-2.0 for used modules; review per-file notices if adding unsupported modules |
| nlohmann/json | C++ project JSON parsing | 3.10.5 (Ubuntu Jammy package) | MIT |
| Gmsh | Constrained curve/surface and cavity meshing | 4.8.4+ds2-2build1 (Ubuntu Jammy) | GPL-2.0-or-later with the Gmsh linking exception for listed libraries; see full package copyright |
| CGNS | Standard mixed and NGON/NFACE mesh exchange | 3.4.0-3 (Ubuntu Jammy local cache) | Upstream zlib-style notice; preserve attribution and mark modified source |
| Qt 6 Core/Gui/Widgets/OpenGLWidgets | Desktop CAD, viewport and background tasks | 6.2.4+dfsg-2ubuntu1.1 | LGPL-3.0/GPL/commercial options; this build uses system shared libraries |
| Qt 6 SVG runtime | Render the application's original SVG widget indicators | 6.2.4-1ubuntu1 | LGPL-3.0/GPL/commercial options; system shared runtime |
| Python jsonschema | Contract/schema tests only | 3.2.0 (Ubuntu Jammy package) | MIT upstream; Debian packaging files have separate GPL-2+ notices |
| GNU Make | Build orchestration only | 4.3 (current host) | GPL-3.0-or-later; build tool is not linked into BabelMesh |

The OCCT license is LGPL 2.1 with the exception published by OCCT. BabelMesh
must preserve the notice that it uses OCCT and include the applicable license
texts/instructions with a binary release. The exception covers object code
incorporating material from OCCT headers; it does not remove all LGPL duties.
The project currently links shared system libraries on Ubuntu.

The hybrid dependency snapshot and cache locations are recorded in
[`dependencies/hybrid-ubuntu-jammy.json`](dependencies/hybrid-ubuntu-jammy.json).
`tools/bootstrap_hybrid_dependencies.sh` extracts distribution development
packages for local development; those packages retain their own copyright
files. The Gmsh distribution used here excludes the Netgen contribution, so
this implementation does not enable its unavailable Netgen optimizer.

Upstream references:

- [OCCT license statement](https://github.com/Open-Cascade-SAS/OCCT),
  [OCCT LGPL exception](https://github.com/Open-Cascade-SAS/OCCT/blob/master/OCCT_LGPL_EXCEPTION.txt)
- [nlohmann/json license and integration](https://github.com/nlohmann/json)
- [Eigen licensing page](https://eigen.tuxfamily.org/index.php?title=Main_Page#License)
- [Gmsh copying conditions and full license](https://gmsh.info/doc/texinfo/#Copying-conditions),
  package copyright: `.cache/dev-packages/gmsh/usr/share/doc/libgmsh-dev/copyright`
- [Qt licensing](https://doc.qt.io/qt-6/licensing.html),
  package copyright: `/usr/share/doc/qt6-base-dev/copyright`
- CGNS package copyright: `.cache/dev-packages/cgns/usr/share/doc/libcgns-dev/copyright`
- Python package metadata: `/usr/share/doc/python3-jsonschema/copyright`
- Eigen package metadata: `/usr/share/doc/libeigen3-dev/copyright`
- OCCT package metadata: `/usr/share/doc/libocct-foundation-dev/copyright`

The build uses system-provided libraries and does not copy their source into
this repository. If a dependency is vendored later, include its full notices
and license texts in the release package and update this inventory.
