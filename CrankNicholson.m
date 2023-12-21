% Clear previous values of the variables.
clear

% Ask if the user wants default values
answer = input("Use defaults? (L=1, T=12000, n=39, m=600, c=1.e-5, mixed)? (No (0) or Yes (1)) ");
if answer == 1
    L = 1;
    T = 50000;
    n = 39;
    m = 1000;
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

bcase = "DN";

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
        if cond1 == 0
            bcase = "DD";
        else
            bcase = "ND";
        end
    end
    if cond2 == 1
        g2 = input("Neumann condition at x = L: ");
        if cond1 == 0
            bcase = "DN";
        else
            bcase = "NN"
        end
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
s = c*dt/(2*(dx.^2));
if s > 0.25
    fprintf("Diffusivity > 0.5, solution is unstable.")
end

% Build the matrix and vector.
if bcase == "DN"
    Adiag = (1+2*s)*ones(n+1,1);
    Asubs = -s*ones(n+1,1);
    Asuper = -s*ones(n+1,1);
    A = spdiags([Asubs,Adiag,Asuper],[-1 0 1],n+1,n+1);
    
    Bdiag = (1-2*s)*ones(n+1,1);
    Bsubs = s*ones(n+1,1);
    Bsuper = s*ones(n+1,1);
    B = spdiags([Bsubs,Bdiag,Bsuper],[-1 0 1],n+1,n+1);
    b = zeros(n+1,1);
end
if bcase == "DD"
    Adiag = (1+2*s)*ones(n,1);
    Asubs = -s*ones(n,1);
    Asuper = -s*ones(n,1);
    A = spdiags([Asubs,Adiag,Asuper],[-1 0 1],n,n);
    
    Bdiag = (1-2*s)*ones(n,1);
    Bsubs = s*ones(n,1);
    Bsuper = s*ones(n,1);
    B = spdiags([Bsubs,Bdiag,Bsuper],[-1 0 1],n,n);
    b = zeros(n,1);
end

z1 = g2;
z2 = g1;
% Build the boundary conditions into the matrix depending on their type.
if cond1 == 0
    b(1)=2.*s*g1;
    z2 = g1;
end
%if cond1 == 1
%    b(1)=s*dx*g1;
%    A(1,2)=-2.*s;
%end
if cond2 == 0
    b(n)=2.*s*g2;
    z1 = (g2 - g1)/L;
end
if cond2 == 1
    b(n+1)=4.*s*dx*g2;
    A(n+1,n)=-2.*s;
    B(n+1, n-1) = s;
    z1 = g2;
end

% Graph the initial condition.
if bcase == "DN"
    x = linspace(x1+dx,x2, n+1);
end
if bcase == "DD"
    x = linspace(x1+dx,x2-dx,n);
end

if answer == 0
    U_tk = input("Initial condition (function): ");
else
    U_tk = 3.*x+sin(2.*pi*x)+1;
end

figure(1)
plot(x,U_tk,'-')
title("Initial condition for temperature distribution")
xlabel("X")
ylabel("U")

% Compute the solution at each time step.
if bcase == "DN"
    U = zeros(n+1,m);
end
if bcase == "DD"
    U = zeros(n,m);
end

U_tk=U_tk';
tvec=zeros(m,1);
U(:,1)=U_tk;
tvec(1)=t;

for k = 1:m
    t = t+dt;
    v = B*U_tk + b;
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
title("Variation of temperature distribution (Backward Euler)")
xlabel("T")
ylabel("X")
zlabel("U")