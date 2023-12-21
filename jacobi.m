% Jacobi Method

% Clears previous input.
clear

% Choose defaults or get user input.
value = input("Do you want to use default values? (No (0) or Yes (1)) ");

if value == 1
    A = [1 6 0 2 3; 2 8 4 8 1; 2 0 4 5 0; 1 0 2 4 5; 1 2 7 0 1];
    x = [1; 1; 1; 1; 1];
    b = [1; 4; 0; 3; 0];
    max = 100;
    tol = 1.e-5;
else
    A = input("Matrix A: ");
    x = input("Guess for solution: ");
    b = input("Vector b");
    max = input("Max number of iterations: ");
    tol = input("Error tolerance: ");
end

% Initialize the iteration.
it = 0;
r = b-A*x;
n = size(A, 1);
error = norm(r)/norm(b);
if error < tol
    output("The approximation is already good enough.");
end

while all(error > tol)
    x_1 = x;
    for i=1:n
        sum = 0;
        for j=1:n
            if j~=i
                sum = sum+A(i,j)*x_1(j);
            end
        end
        x(i) = (1/A(i,i))*(b(i)-sum);
    end
    it=it+1;
    y(it,:)=x;
    error=abs(x_1-x);
end %#ok<*SYNER>

if error > tol
    output("No convergence.");
end
