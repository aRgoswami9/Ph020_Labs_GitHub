PLOTS = explicit_spring_motion.jpg\
 explicit_spring_error.jpg\
 truncation_error.jpg\
 energy_evolution_explicit.jpg\
 energy_evolution_implicit_explicit.jpg\
 implicit_global_error.jpg\
 implicit_phasespace.jpg\
 explicit_phasespace.jpg\
 symplectic_phasespace.jpg\
 symplectic_energy_evolution.jpg

.phony: all plots

all: writeup.pdf

plots: $(PLOTS)

writeup.pdf: writeup.tex plots
	pdflatex writeup.tex

explicit_spring_motion.jpg: numDEQ.py 
	python numDEQ.py 1

explicit_spring_error.jpg: numDEQ.py 
	python numDEQ.py 2

truncation_error.jpg: numDEQ.py
	python numDEQ.py 3

energy_evolution_explicit.jpg: numDEQ.py
	python numDEQ.py 4

energy_evolution_implicit_explicit.jpg implicit_global_error.jpg: numDEQ.py
	python numDEQ.py 5

implicit_phasespace.jpg explicit_phasespace.jpg symplectic_phasespace.jpg: numDEQ.py
	python numDEQ.py phasespace

symplectic_energy_evolution.jpg: numDEQ.py
	python numDEQ 2.3
