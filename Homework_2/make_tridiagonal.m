function [tridiag] = make_tridiagonal(N,a,b,c)
%MAKE_TRIDIAGONAL Summary of this function goes here
%   N - size of matrix
%   a - upper diagonal values
%   b - main
%   c - lower
tridiag = zeros(N);

for i = 1:N;
    for j = 1:N;

        % upper diagonal;
        if i == j+1
            tridiag(i,j) = 1;
        
        % main diagonal
        elseif i == j
            tridiag(i,j) = -2;
           
        % lower diagonal
        elseif i == j-1
            tridiag(i,j) = 1;
        end
    end
end
end

