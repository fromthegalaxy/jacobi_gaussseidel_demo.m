function diff_norm = getNorm2(x_old, x_new)
% Calculate norm-2 of the difference
diff_norm = sqrt(sum((x_new - x_old).^2));

end

function G = iter_matrix(A, A_d_inv, n)
% calculate the iteration matrix
G = eye(n) - A_d_inv * A;

end

function G = iter_matrix2(A_M_inv, A_N, n)
% calculate the iteration matrix
G = -A_M_inv * A_N;
end

%-------------------------------------------------------
%-------------------------------------------------------
% a

function [x, num_iter] = Jacobi(A,b)
% Jacobi: Solves Ax=b using Jacobi iteration method
%   x[k+1] = -A_D^-1 (A_L + A_R) * x[k] + A_d^-1 * b
%
% Inputs:
%   A: n by n nonsingular matrix
%   b: a R^n vector
%
% Outputs:
%   x: the approximate solution to Ax=b
%   num_iter: the number of iterations happened
%
% Stops after 1000 iterations or when norm-2 of the difference
%       between two successive iterates is less than 10−5

max_iter = 1000;
tol = 10e-5;

num_iter = 0;
n = length(b);          % Dimension of the system
x_old = zeros(n, 1);    % initialize first guess = zero
A_d = zeros(n, n);      % initalize the diagonal matrix

for i = 1:n
    % gets the diagonal entries
    A_d(i, i) = A(i, i);
end

A_d_inv = inv(A_d); % inverse A_d
G = iter_matrix(A, A_d_inv, n); % iteration matrix G

for k = 1:max_iter
    % iterate on x using iteration matrix G
    x_new = zeros(n,1);
    x_new = G * x_old + A_d_inv * b;
    num_iter = k; % update number of iterations

    if getNorm2(x_old, x_new) < tol
        x_old = x_new; % update x
        %if tol is reached, exit the for loop
        break;
    else
        x_old = x_new; % update x
    end

end
x = x_old;

end

%-------------------------------------------------------
%-------------------------------------------------------
% b

function [x, num_iter] = Gauss_Seidel(A,b)
% GAUSS_SEIDEL  Solve Ax = b with the Gauss-Seidel iterative method.
%
%   [X, NUM_ITER] = GAUSS_SEIDEL(A, B) computes an approximate solution X
%   of the linear system A*X = B by Gauss-Seidel iteration, starting from
%   the zero vector.  The iteration stops when
%
%           ||x^{k+1} – x^{k}||_2 < 1 × 10⁻⁵   or   k = 1000 ,
%
%   whichever occurs first, and returns the final iterate X together with
%   the number of iterations NUM_ITER.
%
%   INPUTS
%   ------
%   A : n-by-n nonsingular matrix (double).
%   B : n-by-1 right-hand-side vector (double).
%
%   OUTPUTS
%   -------
%   X        : n-by-1 Gauss-Seidel approximation to the solution of A*X = B.
%   NUM_ITER : Number of iterations executed before termination.
%
%   ALGORITHM
%   ---------
%   Decompose A = (D + L) + U where
%       D = diag(A),  L = strict lower part,  U = strict upper part.
%   Iterate
%       (D + L) x^{k+1} = b – U x^{k}
%   which is implemented by successively overwriting the components of
%   x^{k+1} in place, so each new component immediately influences the next.
%
%   EXAMPLE
%   -------
%       A = [ 4 -1  0;
%             -1  4 -1;
%              0 -1  3];
%       b = [15; 10; 10];
%       [x, iters] = Gauss_Seidel(A, b);


max_iter = 1000;
tol = 10e-5;

num_iter = 0;
n = length(b);          % Dimension of the system
x_old = zeros(n, 1);    % initialize first guess = zero
A_M = zeros(n, n);      % initalize the M matrix (A_d + A_l)
A_N = zeros(n, n);      % initalize the N matrix (A_r)


for i = 1:n
    for j = 1:n
        if j <= i
            % gets the A_d and A_l entries
            A_M(i, j) = A(i, j);
        else
            % gets the A_r entries
            A_N(i,j) = A(i, j);
        end
    end
end

A_M_inv = inv(A_M); % inverse A_M
G = iter_matrix2(A_M_inv, A_N, n); % iteration matrix G

for k = 1:max_iter
    % iterate on x using iteration matrix G
    x_new = zeros(n,1);
    x_new = G * x_old + A_M_inv * b;
    num_iter = k; % update number of iterations

    if getNorm2(x_old, x_new) < tol
        x_old = x_new; % update x
        %if tol is reached, exit the for loop
        break;
    else
        x_old = x_new; % update x
    end

end

x = x_old;
end



[m,s] = Jacobi([ 4  -1   0; -1   4  -1; 0  -1   3], [15; 10; 10]);

disp(m);
disp(s);


%-------------------------------------------------------
%-------------------------------------------------------
% c

function A = gen_matrix(n)
% generate matrix of 1 in R_nxn, with n as the diagnoal entries
A = ones(n, n);

for i = 1:n
    % changing diagonal entries
    A(i, i) = n;
end

end


n = [5; 10; 20; 50; 100; 250; 500];

% Initialize tables
K_A = zeros(7, 1);          % table of condition numbers of A
K_r_Jacobi= zeros(7, 1);    % table of Kr with Jacobi
Jacobi_iter = zeros(7, 1);  % table of number of iter Jacobi performed
K_r_GS = zeros(7, 1);       % table of Kr with G-S
GS_iter = zeros(7, 1);      % table of number of iter G-S performed

for i = 1:length(n)

    A = gen_matrix(n(i));   % generate matrix
    b = ones(n(i), 1);      % generate b
    x = A\b;                % get the true output
    [x_Jacobi,n_Jacobi] = Jacobi(A,b);  % get the output with Jacobi
    [x_GS,n_GS] = Gauss_Seidel(A,b);    % get the ouput with G-S

    K_A(i) = cond(A);       % cond number of A
    K_r_Jacobi(i) = norm(x - x_Jacobi) / norm(x);   % Kr of Jacobi
    Jacobi_iter(i) = n_Jacobi;  
    K_r_GS(i) = norm(x - x_GS) / norm(x);   % Kr of G-S
    GS_iter(i) = n_GS;      
end


% Create the table
T = table(n, K_A, K_r_Jacobi, Jacobi_iter, K_r_GS, GS_iter);

% % Display the table
disp(T);


