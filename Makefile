CC     = cc
CXX    = c++
CFLAGS = -O3

# SuperLU build settings for puffup_c's --solver sparse path.
# `make sparse` builds src/puffup_c_sparse linked against system SuperLU.
# Headers live at /usr/include/superlu (Debian/Ubuntu) or /opt/homebrew/include
# (macOS Homebrew). BLAS comes from libblas (Linux) or Accelerate.framework (macOS).
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
  SUPERLU_INC = -I/opt/homebrew/include
  SUPERLU_LIB = -L/opt/homebrew/lib -lsuperlu -framework Accelerate
  LAPACK_LIB  = -framework Accelerate
  CGAL_INC    = -I/opt/homebrew/include
  CGAL_LIB    = -L/opt/homebrew/lib -lgmp -lmpfr
else
  SUPERLU_INC = -I/usr/include/superlu
  SUPERLU_LIB = -lsuperlu -lblas
  LAPACK_LIB  = -llapack -lblas
  CGAL_INC    =
  CGAL_LIB    = -lgmp -lmpfr
endif

all: src/horou_c src/horoz_c src/proof_c src/prove_c src/puffup_c src/realize_c src/embed_check

src/horou_c: src/horou_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/horoz_c: src/horoz_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/proof_c: src/proof_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

# Euclidean existence prover (Matt Ellison's argument, arXiv:2312.05376).
# Requires LAPACK for SVD (dgesvd). Do NOT use -ffast-math — error bounds
# rely on IEEE 754 semantics.
src/prove_c: src/prove.c
	$(CC) $(CFLAGS) -o $@ $< $(LAPACK_LIB) -lm

src/puffup_c: src/puffup_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

# Sparse-LU variant: same source, built with -DHAVE_SUPERLU and linked to
# SuperLU + BLAS. With this build, sparse is the default solver; pass
# --solver dense to fall back to the dense path for A/B comparison.
src/puffup_c_sparse: src/puffup_c.c
	$(CC) $(CFLAGS) -DHAVE_SUPERLU $(SUPERLU_INC) -o $@ $< $(SUPERLU_LIB) -lm

sparse: src/puffup_c_sparse

src/realize_c: src/realize_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

# Certified embeddedness check via CGAL Polygon_mesh_processing.
# Uses Exact_predicates_inexact_constructions_kernel — verdicts are
# certified, not float-precision.
src/embed_check: src/embed_check.cpp
	$(CXX) -O2 -std=c++17 $(CGAL_INC) -o $@ $< $(CGAL_LIB)

clean:
	rm -f src/horou_c src/horoz_c src/proof_c src/prove_c src/puffup_c src/puffup_c_sparse src/realize_c src/embed_check
