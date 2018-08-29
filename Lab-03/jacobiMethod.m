%====================== JACOBI ROTATION METHOD ======================%
%
%====================================================================%
%
%   Calculates the eigenvalues and the corresponding eigenvectors
%   of the input matrix A (square-symmetric).
%
%   Returns a row vector "eigenvalues", containing eigenvalues of A;
%   and a matrix "eigenvectors", with the corresponding eigenvectors
%   stored along the columns.
%
%====================================================================%

function [eigenvalues, eigenvectors] = jacobiMethod(A)

    % Dimension of the matrix A
    n = size(A, 1);

    tolerance = 10e-10; changes = 1;
    U = eye(n);

    while changes > 0

        changes = 0;

        % Zeroing the off-diagonal elements one-by-one
        for p = 1:n
            for q = p+1:n

                if abs(A(p, q)) < tolerance
                    continue;
                endif

                changes = 1;

                % Calculating the theta that makes A[p, q] = A[q, p] = 0
                tau = 2 * A(p, q) / (A(q, q) - A(p, p));
                theta = 0.5 * atan(tau);
                % Generating the rotation matrix required
                R = eye(n);
                R(p, p) = cos(theta);
                R(p, q) = sin(theta);
                R(q, p) = -sin(theta);
                R(q, q) = cos(theta);

                % Performing the rotation
                A = R' * A * R;

                % Adding the new rotation to the composite operator
                U = U * R;

            endfor
        endfor
    
    endwhile

    fprintf("\n")
    disp(U);
    fprintf("\n")
    disp(A);

    eigenvalues = sum(A);
    eigenvectors = U;

endfunction