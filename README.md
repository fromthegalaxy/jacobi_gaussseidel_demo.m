# jacobi_gaussseidel_demo.m
Implements Jacobi and Gauss-Seidel iterations, then for each n ∈ {5, 10, 20, 50, 100, 250, 500} solves the matrix with n on the diagonal and 1 elsewhere, logs cond(A), relative errors, and iteration counts, and prints a comparison table.


Solves Ax = b for a family of test matrices using **Jacobi** and
  **Gauss–Seidel** iterations implemented from scratch, compares their
  accuracy and convergence speed, and prints a summary table.

-------------------------------------------------------------------------
Problem setup
  • A ∈ ℝⁿˣⁿ has   A(i,i) = n   and   A(i,j) = 1   for i ≠ j
  • b = ones(n,1)
  • n ∈ {5, 10, 20, 50, 100, 250, 500}

Iterative solvers
  • Jacobi:  stop when  ∥x^{k+1} – x^{k}∥₂ < 1e-5  or k = 1000
  • Gauss-Seidel: same stopping rule

For each n the script
  1) generates A and b,
  2) gets the “exact” solution  x_e = A\b ,
  3) runs Jacobi  →  x_j ,  iterations i_j
  4) runs Gauss–Seidel  →  x_gs ,  iterations i_g
  5) records
        – cond(A)
        – relative errors  ∥x_e–x_j∥/∥x_e∥  and  ∥x_e–x_gs∥/∥x_e∥
        – iteration counts  i_j , i_g
  6) displays a results table.

Output
  A MATLAB table printed to the console with columns:
      n | cond(A) | relErr_Jacobi | iter_Jacobi | relErr_GS | iter_GS
