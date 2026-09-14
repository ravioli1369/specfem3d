#!/bin/bash
# Runs both SH FK-injection validation cases (homogeneous + layered) and the
# independent Thomson-Haskell check. Build the binaries first from the root:
#   make xmeshfem3D xgenerate_databases xspecfem3D
set -e
here=$(cd "$(dirname "$0")" && pwd)
bin=$here/../../../bin
for case in homogeneous layered; do
  echo "=== $case ==="
  cd "$here/$case"
  rm -rf OUTPUT_FILES DATABASES_MPI; mkdir -p OUTPUT_FILES DATABASES_MPI
  for x in xmeshfem3D xgenerate_databases xspecfem3D; do ln -sf "$bin/$x" .; done
  ./xmeshfem3D > OUTPUT_FILES/log_mesh.txt 2>&1
  ./xgenerate_databases > OUTPUT_FILES/log_gendb.txt 2>&1
  ./xspecfem3D > OUTPUT_FILES/log_solver.txt 2>&1
  cd "$here"
done
echo "=== validation ==="
python3 "$here/validate_sh_fk.py" "$here/homogeneous" "$here/layered"
