CC     = cc
CFLAGS = -O3

# SuperLU build settings for puffup_c's --solver sparse path.
# `make sparse` builds src/puffup_c_sparse linked against system SuperLU.
# Headers live at /usr/include/superlu (Debian/Ubuntu) or /opt/homebrew/include
# (macOS Homebrew). BLAS comes from libblas (Linux) or Accelerate.framework (macOS).
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
  SUPERLU_INC = -I/opt/homebrew/include
  SUPERLU_LIB = -L/opt/homebrew/lib -lsuperlu -framework Accelerate
else
  SUPERLU_INC = -I/usr/include/superlu
  SUPERLU_LIB = -lsuperlu -lblas
endif

all: src/horou_c src/horoz_c src/proof_c src/puffup_c src/realize_c

src/horou_c: src/horou_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/horoz_c: src/horoz_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/proof_c: src/proof_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/puffup_c: src/puffup_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

# Sparse-LU variant: same source, built with -DHAVE_SUPERLU and linked to
# SuperLU + BLAS. Use ./puffup_c_sparse --solver sparse to engage.
src/puffup_c_sparse: src/puffup_c.c
	$(CC) $(CFLAGS) -DHAVE_SUPERLU $(SUPERLU_INC) -o $@ $< $(SUPERLU_LIB) -lm

sparse: src/puffup_c_sparse

src/realize_c: src/realize_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f src/horou_c src/horoz_c src/proof_c src/puffup_c src/puffup_c_sparse src/realize_c
