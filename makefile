.PHONY: git

DO:
	clear
	make -B matrix

matrix:
	g++ -fopenmp -g main.cpp -Wall -o something
	./something

git: git_push

git_push:
	@echo "Enter your commit message:"
	@read commit_message; \
	git add .; \
	git commit -m "$$commit_message"; \
	git push;
