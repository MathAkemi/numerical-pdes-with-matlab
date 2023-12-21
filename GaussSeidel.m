% Gauss-Seidel Method

format long;
a0file = fopen('fa0error.out','w');
afile = fopen('faerror.out','w');

% Choose defaults or get user input.
value = input("Do you want to use default values? (No (0) or Yes (1)) ");

if value == 1
    A = [5 0 -2; 3 5 1; 0 -3 4];
    b = [7 2 -4]';
    tol = 1.e-3;
else
    A = input("Matrix A: ");
    b = input("Vector b");
    tol = input("Error tolerance: ");
end

it = 10;

% Initialize the relevant matricies.
D = diag(diag(A));
L = tril(A) - D;
U = triu(A) - D;

n = size(A,1);
x_k = zeros(n,1);
x_k_1 = zeros(n,1);
k=0;
a0error = norm(b-A*x_k);
fprintf(a0file, '%d %8.6f \n', k, a0error);
numit = 0;

% alpha = 1
while a0error > tol
    for i=1:n
        sum1=0;
        for j=1:i-1
            sum1 = sum1 + A(i,j)*x_k_1(j);
        end
        
        sum2=0;
        for j=i+1:n
            sum2 = sum2 + A(i,j)*x_k(j);
        end
        
        x_k_1(i) = (b(i) - sum1 - sum2)/A(i,i);
    end
    a0error = norm(b-A*x_k_1);
    k = k+1;
    fprintf(a0file, '%d %8.6f \n', k, log(a0error));
    x_k=x_k_1;
    numit = numit+1;
end
fprintf("Data for alpha = 1: ")
x_k
numit
a0error

% alpha != 1
k = 0;
alpha = input("Value of alpha: ");
x_k = zeros(n,1);
aerror = norm(b-A*x_k);
fprintf(afile, '%d %8.6f \n', k, aerror);
errorVals2 = zeros(10,1);
numit = 0;
while aerror > tol
    for i=1:n
        sum1=0;
        for j=1:i-1
            sum1 = sum1 + A(i,j)*x_k_1(j);
        end
        
        sum2=0;
        for j=i+1:n
            sum2 = sum2 + A(i,j)*x_k(j);
        end
        
        x_k_1(i) = (1-alpha)*x_k(i) + alpha*(b(i) - sum1 - sum2)/A(i,i);
    end
    aerror = norm(b-A*x_k_1);
    k = k+1;
    fprintf(afile, '%d %8.6f \n', k, log(aerror));
    x_k=x_k_1;
    numit = numit+1;
end

fprintf("Data for alpha: ")
x_k
numit
aerror

fclose(a0file);
fclose(afile);
load fa0error.out
load faerror.out

k_vec1 = fa0error(:,1);
k_vec2 = faerror(:,1);
a0error_vec = fa0error(:,2);
aerror_vec = faerror(:,2);

figure(1)
plot(k_vec1, a0error_vec, 'r-')
hold on
plot(k_vec2, aerror_vec, 'b-')
xlabel ('no of iterations')
ylabel ('Error')

title('Residual after each iteration using Gauss-Seidel method')
