function [eigenvector, eigenvalue, iter] = powerMethod(A, tol)
    
    eigenvector = zeros(size(A, 1), 1);
    x = rand(size(A, 1), 1);
    eigenvector = 10 * x;
    iter = 0;
    maxIter = 10e+09;

    while iter < maxIter

        eigenvector = A * x;
        eigenvalue = norm(A * eigenvector) / norm(eigenvector);

        if (eigenvalue - norm(A * x) / norm(x)) < tol
            break;
        endif

        x = eigenvector;
        iter += 1;

    endwhile

    eigenvalue = norm(A * eigenvector) / norm(eigenvector);
    eigenvector = eigenvector ./ norm(eigenvector);

endfunction
