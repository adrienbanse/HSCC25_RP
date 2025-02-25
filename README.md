# Repeatability Evaluation Package for "Memory-dependent abstractions of stochastic systems through the lens of transfer operators"

This repository is the Repeatability Evaluation Package (REP) for the paper 
> Adrien Banse, Giannis Delimpaltadakis, Luca Laurenti, Manuel Mazo Jr. and Raphaël M. Jungers
> "Memory-dependent abstractions of stochastic systems through the lens of transfer operators"
> in Hybrid Systems Computational & Control, ACM (2025).

This REP reproduces exactly Figure 4 and Figure 5 of the paper. After execution, both figures will be saved in `HSCC25_RP` as `output/figure_4.pdf` and `output/figure_5.pdf`.

## Running instructions

This REP is written in Julia and can be executed in two different ways: either (1) via the notebook `RP.ipynb`, or (2) directly by executing the replica `RP.jl` via Docker. The first option (notebook) allows for a better readability of the code and the documentation, while the second (replica) is fully reproducible (if Docker is used), simpler and faster.

### Option 1: run the notebook

1. Install [Julia](https://julialang.org/downloads/).
2. Navigate to the `HSCC25_RP` repository, and open the Julia REPL. 
3. (If needed) Install `IJulia` by running
   ```julia-repl
   julia> using Pkg
   julia> Pkg.add("IJulia")
   ```
4. Make sure that you have added the needed packages to your environment.
     - Either activate the `HSCC25_RP` environment. Open the REPL in `HSCC25_RP/` and run
       ```julia-repl
       julia> using Pkg
       julia> Pkg.activate(".")
       julia> Pkg.instantiate()
       ```
     - Or add yourself the needed packages by running
       ```julia-repl
       julia> using Pkg
       julia> Pkg.add("Distributions"); Pkg.add("LinearAlgebra"); Pkg.add("Random"); Pkg.add("MatrixEquations"); Pkg.add("Plots")
       ```
5. Open the notebook
   ```julia-repl
   julia> using IJulia; notebook()
   ```
   It is likely that it is asked to install `Jupyter` if it is the first time that you use `IJulia`, answer `y` to do so.
6. Open in Jupyter the notebook `RP.ipynb`, and execute it cell by cell.
  
### Option 2: run the replica (Docker)

1. Open Docker, and make sure that you allow for at least 16GB of memory use.
2. Navigate to the `HSCC25_RP` repository.
3. Run `make docker`

If, in addition, you want to make sure that `RP.jl` is indeed a replica of the notebook, you can re-execute the replica command. First, [install Jupytext](https://jupytext.readthedocs.io/en/latest/), and then run `jupytext --to jl RP.ipynb` in your terminal. This command will replace the existing `RP.jl` file with the output of `jupytext`.

## Documentation

All the code source needed for the execution of this REP is contained in the `RP.ipynb` notebook. The latter contains docstrings for the implemented functions, as well as further explanations/links with the paper written in Markdown.
