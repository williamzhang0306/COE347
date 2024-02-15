%{
Part 1
Compute the eigenvalues of the matrix for N = 10 and confirm that they are
λi = −2(1 − cos(πi/(N + 1)) i = 1, . . . , N 
%}

N = 10;

tridiag = make_tridiagonal(N,1,-2,1)

% eigen values
disp('Eigen Values for N = 10')
eig(tridiag)

% given function\
disp('-2* (1 - cos(pi*i/(N+1) ) )')
x = 1:N;
disp( -2*(1 - cos(pi*x'./(N+1)))  )
           


%{
Part 2:
Next, plot max |λ|, i.e. the maximum absolute value of all eigenvalues of TN versus N for N =
1, . . . , 20. Comment on the behavior of max |λ| as N changes. Reconciliate your finding with the
expression for λi provided above in Eq. (3).
Turn in plots that support your answers and explanations.
%}



N_values = zeros(20,1);

for idx = 1:20
    tridiag = make_tridiagonal(idx, 1 ,-2, 1);
    eigen_vals = abs(eig(tridiag));
    N_values(idx) = max(eigen_vals);
end

figure
plot(1:20, N_values, 'o-')
xlabel('N')
ylabel('|λ|')
title('Matrix Size vs Largest Eigenvalue', 'FontSize', 12)
saveas(gcf, 'Eigenvalues_plot.png');
