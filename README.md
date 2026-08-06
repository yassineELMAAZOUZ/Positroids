# Positroids.jl

Experimental Julia tool for positroids, decorated permutations, Grassmann
necklaces, total positivity. 
The combinatorial core is dependency-free; Plots.jl is used only by `draw_chords`.

This repository consolidates the Julia prototypes in `Positroids_code.jl` and
`decorated_permutations.jl` and the Mathematica implementation
(`positroids.m`) of Jacob L. Bourjaily.


### Install and load locally

Download the folder. Then do the following:

1- You first need to have Julia installed.

2- In your terminal: 
          > julia --project="path_to_positroids_folder"
   Then in the Julia REPL:
   
          Julia> using Positroids; interactive_session();



<span style="color: red;">Fair warning</span>: This is still work in progress!! If you encounter any bugs or errors feel free to let me know and I can try to fix it.

