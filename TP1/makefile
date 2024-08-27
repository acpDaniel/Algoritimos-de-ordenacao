SIZES = 5000 10000 15000 20000 25000 30000 35000 40000
ALGORITHMS = 9

all:
	g++ -o tp1v2.out src/tp1v2.c

run_alg3:
	@echo "Vetor com chaves que pode repetir, e ordenacao aleatoria"
	@for size in $(SIZES); do \
		for algo in $(ALGORITHMS); do \
			./tp1v2.out $$size 0 0 0 $$algo; \
		done; \
	done

run_alg2:
	@echo "Vetor com chaves que pode repetir, e ordem decrescente"
	@for size in $(SIZES); do \
		for algo in $(ALGORITHMS); do \
			./tp1v2.out $$size 0 1 0 $$algo; \
		done; \
	done

run_alg1:
	@echo "Vetor com chaves que pode repetir, e ordem crescente"
	@for size in $(SIZES); do \
		for algo in $(ALGORITHMS); do \
			./tp1v2.out $$size 0 1 1 $$algo; \
		done; \
	done

clean:
	rm -f tp1caso2.out
