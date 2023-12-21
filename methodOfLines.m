%% Method of Lines

% Solving U_t = beta*U_xx
global A b;

% Ask if the user wants default values
answer = input("Use defaults? (L=1, T=12000, n=39, m=600, c=1.e-5, mixed)? (No (0) or Yes (1)) ");
if answer == 1
    L = 1;
    T = 12000;
    n = 39;
    c = 1.e-5;
    cond1 = 0;
    cond2 = 1;
end

% Get user input for parameters.
if answer == 0
    L = input("Total length: ");
    T = input("Total time: ");
    n = input("Number of space intervals: ");
    c = input("Constant of diffusifity: ");
    cond1 = input("Type of boundary condition at x=0 (Dirichlet (0) or Neumann (1)): ");
    cond2 = input("Type of boundary condition at x=L (Dirichlet (0) or Neumann (1)): ");
end

% Length and spacing.
x1 = 0;
x2 = L;
dx = L/(n+1.);
t = 0;

% Gain parameter.
s = c/(dx^2);

% Construct the matrix A and vector b.
Adiag = (-2*s)*ones(n,1);
Asubs = s*ones(n,1);
Asuper = s*ones(n,1);
A = spdiags([Asubs,Adiag,Asuper],[-1 0 1], n, n);
b = zeros(n,1);

% Get the conditions depending on type.
if answer == 0
    if cond1 == 0
        g1 = input("Dirichlet condition at x = 0: ");
    end
    if cond1 == 1
        g1 = input("Neumann condition at x = 0: ");
    end
    if cond2 == 0
        g2 = input("Dirichlet condition at x = L: ");
    end
    if cond2 == 1
        g2 = input("Neumann condition at x = L: ");
    end
else
    g1 = 1;
    g2 = 2;
end

% Build the boundary conditions into the matrix depending on their type.
if cond1 == 0
    b(1)=s*g1;
end
if cond1 == 1
    b(1)=s*dx*g1;
    A(1,1)=-s;
end
if cond2 == 0
    b(n)=s*g2;
end
if cond2 == 1
    b(n)=s*dx*g2;
    A(n,n)=-s;
end

% Get and plot the initial condition (function).
x = linspace(x1+dx,x2-dx, n)';
if answer == 0
    U0 = input("Initial condition (function): ");
else
    U0 = 3.*x+sin(2.*pi*x)+1;
end

figure(3)
plot(x,U0,'-')
title("Initial condition for temperature distribution")
xlabel("x")
ylabel("U")

% Solve and plot the final function using ode45.
[t, U] = ode45("Uprime",[0,T], U0);

% Graph a specific point in time.
value = input("Do you want a graph at a specific point in time? (No (0) or Yes (1)) ");
if value == 1
    k = input("At what value of t from 1 to m? ") ;
    figure(3)
    vect = U(:,k)
    plot(x, vect, "-")
    title("Graph at specific point in time")
    xlabel("x")
    ylabel("U")
end

figure(4)
mesh(t,x,U')
title("Variation of temperature distribution with time (method of lines)")
xlabel("t")
ylabel("x")
zlabel("U")
