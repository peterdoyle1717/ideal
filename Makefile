CC     = cc
CFLAGS = -O3

all: src/horou_c src/horoz_c src/proof_c src/puffup_c src/realize_c

src/horou_c: src/horou_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/horoz_c: src/horoz_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/proof_c: src/proof_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/puffup_c: src/puffup_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

src/realize_c: src/realize_c.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f src/horou_c src/horoz_c src/proof_c src/puffup_c src/realize_c
