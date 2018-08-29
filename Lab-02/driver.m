A = [5 8 4 5; 8 4 3 1; 4 3 6 2; 5 1 2 7];
tol = 10e-05;

fprintf("Working on A =\n");
disp(A);

[lEigenvector, lEigenvalue] = powerMethod(A, tol);
[sEigenvector, sEigenvalue] = powerMethod(inv(A), tol);
sEigenvalue = 1/sEigenvalue;

fprintf("\nThe largest eigenvalue of A is\n");
disp(lEigenvalue);
fprintf("And the corresponding eigenvector is:\n");
disp(lEigenvector);

fprintf("\nThe smallest eigenvalue of A is\n");
disp(sEigenvalue);
fprintf("And the corresponding eigenvector is:\n");
disp(sEigenvector);

[impEigenvector, impEigenvalue] = powerMethod(inv(A - eye(size(A, 1))*(-4.2)), 10e-05);
impEigenvalue = -4.2 + 1/impEigenvalue;
fprintf("\nGiven that another eigenvalue is = -4.2, its improved value is\n");
disp(impEigenvalue);
fprintf("And the corresponding eigenvector is:\n");
disp(impEigenvector);
