CC=			gcc
CFLAGS=		-std=gnu99 -Wall -O3
EXE=		mbmap-one mbmap-batch ex-lift-check ex-place-check ex-group-check

.PHONY:all clean

all:$(EXE)

../libminibwa.a:
	(cd ..; make libminibwa.a)

mbmap-one:ex-one.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

mbmap-batch:ex-batch.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

ex-lift-check:ex-lift-check.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

ex-place-check:ex-place-check.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

ex-group-check:ex-group-check.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

clean:
	rm -f $(EXE)
