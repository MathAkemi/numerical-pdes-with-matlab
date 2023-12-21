% Clear previous values of the variables.
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

% Get user input for parameters.
if answer == 0
    L = input("Total length: ");
    T = input("Total time: ");
    n = input("Number of space intervals: ");
    m = input("Number of time intervals: ");
    c = input("Constant of diffusifity: ");
    cond1 = input("Type of boundary condition at x=0 (Dirichlet (0) or Neumann (1)): ");
    cond2 = input("Type of boundary condition at x=L (Dirichlet (0) or Neumann (1)): ");
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

% Set boundary values and compute steps.
x1 = 0;
x2 = L;
t = 0;
dx = L/(n+1.);
dt = T/m;

% Get the diffusivity constant and check stability.
s = c*dt/(dx^2);
if s > 0.5
    fprintf("Diffusivity > 0.5, solution is unstable.")
end

% Build the matrix and vector.
Adiag = (1+2*s)*ones(n+1,1);
Asubs = -s*ones(n+1,1);
Asuper = -s*ones(n+1,1);
A = spdiags([Asubs,Adiag,Asuper],[-1 0 1],n+1,n+1);
b = zeros(n+1,1);

z1 = 0;
z2 =0;

% Build the boundary conditions into the matrix depending on their type.
if cond1 == 0
    b(1)=s*g1;
    z2 = g1;
end
if cond1 == 1
    b(1)=-dx*g1;
    A(1,2)=-2*s;
end
if cond2 == 0
    b(n+1)=s*g2;
end
if cond2 == 1
    b(n+1)=4.*s*dx;
    A(n+1,n)=-2.*s;
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
if cond1 == 1 && cond2 == 1
    fprintf("Steady state has no time variance, conditions impossible")
end

% Graph the initial condition.
x = linspace(x1+dx,x2, n+1);
if answer == 0
    U_tk = input("Initial condition (function): ");
else
    U_tk = 2.*x+sin(2.*pi*x)+1;
end

figure(1)
plot(x,U_tk,'-')
title("Initial condition for temperature distribution")
xlabel("X")
ylabel("U")

% Compute the solution at each time step.
U = zeros(n+1,m);
U_tk=U_tk';
tvec=zeros(m,1);
U(:,1)=U_tk;
tvec(1)=t;

for k = 1:m
    t = t+dt;
    v = U_tk + b;
    U_tk_1 = A\v;

    % Move to the next step in time.
    U(:,k) = U_tk_1;
    U_tk = U_tk_1;
    tvec(k)=t;
end

steady = z1*x + z2;

% Graph a specific point in time.
value = input("Do you want a graph at a specific point in time? (No (0) or Yes (1)) ");
if value == 1
    k = input("At what value of t from 1 to m? ") ;
    figure(3)
    hold on
    vect = U(:,k);
    plot(x, vect, "r")
    plot(x, steady, "b")
    title("Plot against the steady-state solution at specific point in time")
    xlabel("x")
    ylabel("U")
end

% Graph the final figure.
figure(2)
mesh(tvec,x,U)
title("Variation of temperature distribution")
xlabel("T")
ylabel("X")
zlabel("U")