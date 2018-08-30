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
                c = cos(theta);
                s = sin(theta);
                R(p, p) = c;
                R(p, q) = s;
                R(q, p) = -s;
                R(q, q) = c;

                % Performing the rotation
                newA = A;
                newA(p, p+1:end) = c * A(p, p+1:end) - s * A(q, p+1:end);
                newA(p+1:end, p) = newA(p, p+1:end)';
                newA(q, q+1:end) = c * A(q, q+1:end) + s * A(p, q+1:end);
                newA(q+1:end, q) = newA(q, q+1:end)';
                newA(p, p) = c^2 * A(p, p) + s^2 * A(q, q) - 2 * c * s * A(p, q);
                newA(q, q) = s^2 * A(p, p) + c^2 * A(q, q) + 2 * c * s * A(p, q);
                newA(p, q) = newA(q, p) = (c^2 - s^2) * A(p, q) + c * s * (A(p, p) - A(q, q));
                A = newA;
                
                % Adding the new rotation to the composite operator
                U = U * R;

            endfor
        endfor
    
    endwhile

    eigenvalues = sum(A);
    eigenvectors = U;

endfunction