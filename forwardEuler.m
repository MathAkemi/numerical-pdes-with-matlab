% Clears the previous state of every variable to avoid errors.
clear

% Ask if the user wants default values
answer = input("Use defaults? (L=1, T=12000, n=39, m=600, c=1.e-5, mixed)? (No (0) or Yes (1)) ");
if answer == 1
    L = 1;
    T = 12000;
    n = 39;
    m = 600;
    c = 1.e-5;
    cond1 = 0;
    cond2 = 1;
end

% Get basic variables and the type of conditions.
if answer == 0
    L = input("Total length: ");
    T = input("Total time: ");
    n = input("Number of space intervals: ");
    m = input("Number of time intervals: ");
    c = input("Constant of diffusifity: ")
    cond1 = input("Type of boundary condition at x=0 (Dirichlet (0) or Neumann (1)): ");
    cond2 = input("Type of boundary condition at x=L (Dirichlet (0) or Neumann (1)): ");
end

% Get the increments in space and time, and s.
x0 = 0;
x1 = L;
t=0;
dx = L/(n+1);
dt = T/m;
s = c*dt/dx^2;

% Terminate the program if the method would be unstable (s > 0.5).
if s > 0.5
    fprintf("s > 0.5 - method is unstable.")
end

% Get the initial condition.
x=linspace(x0+dx,x1-dx,n)';
if answer == 0
    U_tk = input("Initial condition (function in terms of x): ")
else
    U_tk = 3.*x+sin(2.*pi*x)+1;
end

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

% Build a matrix to increment the finite difference solution.
Adiag = (1-2*s)*ones(n,1);
Asubs=s*ones(n,1);
Asuper=s*ones(n,1);
A = spdiags([Asubs,Adiag,Asuper],[-1 0 1],n,n);

b = zeros(n,1);

% Build the boundary conditions into the matrix depending on their type.
if cond1 == 0
    b(1)=s*g1;
end
if cond1 == 1
    b(1)=-s*dx*g1;
    A(1,1)=1-3*s;
end
if cond2 == 0
    b(n)=s*g2;
end
if cond2 == 1
    b(n)=s*dx*g2;
    A(n,n)=1-s;
end

z1 = 0;
z2 = 0;
if cond1 == 0 && cond2 == 0
    z2 = g1;
    g1 = (g2-z1)/L
end
if cond1 == 0 && cond2 == 1
    z2 = g1;
    z1 = g2;
end
if cond1 == 1 && cond2 == 0
    z1 = g1;
    z2 = g2 - g1*L;
end
if cond1 == 1 && cond2 == 0
    fprintf("Steady state has no time variance, conditions impossible")
end

% Build the display for the initial condition.
figure(1)
plot(x,U_tk,'-')
title("Initial condition for temperature distribution")
xlabel("x")
ylabel("U")

init = U_tk;
% Store the solution for each time in the solution matrix.
U=zeros(n,m);
tvec=zeros(m,1);

U(:,1)=U_tk;
tvec(1)=t;

steady = z1*x + z2;

% Move the solution forward for each step in time.
for k = 1:m
    t = t+dt;
    U_tk_1 = A*U_tk+b;
    
    U(:,k) = U_tk_1;
    U_tk = U_tk_1;
    tvec(k) = t;
end

% Graph a specific point in time.
value = input("Do you want a graph at a specific point in time? (No (0) or Yes (1)) ");
if value == 1
    k = input("At what value of t from 1 to m? ");
    figure(3)
    hold on
    vect = U(:,k);
    plot(x, vect, "r")
    plot(x, steady, "b")
    title("Graph at specific point in time")
    xlabel("x")
    ylabel("U")
end


% Plot the graph of the solution.
figure(2)
mesh(tvec,x,U)
title("Temperature distribution over time")
xlabel("t")
ylabel("x")
zlabel("U")