% Ask if the user wants default values
answer = input("Do you want to use defaults (L = 2pi, n = 128, T = 9, c(x) = 1/5 + sin(x-1)^2, U_kt_1 = exp(-100*(x-1).^2)? (no = 0, yes = 1): ");
if answer == 0
    L = input("Length of medium: ");
    T = input("Maximum time: ");
    n = input ("Number of x intervals: ");
else
    L = 2.*pi;
    T = 9.;
    n = 512;
end

dx = L/(2*n);
dt = dx/4;
x = linspace(0,L-dx,2*n)';

wavenu = [0:n -n+1:-1]';
if answer == 0
    c = input("Function for the wave speed: ");
    U_tk_1 = input("Initial condition at t = 0: ");
    U_tk_2 = exp(-100*(x - (dt/5.) - 1).^2);
else
    c = 1/5 + sin(x-1).^2;
    U_tk_1 = exp(-100*(x - 1).^2);
    U_tk_2 = exp(-100*(x - (dt/5.) - 1).^2);
end

figure(1)
plot(x, U_tk_1, '-')
title("Initial condition for U")
xlabel("x")
ylabel("U")

% Store the solution at every 100th timestep to plot later.
m = round(T/dt);
U = zeros(2*n, round(m/100));
tvec = zeros(round(m/100), 1);

t = 0;
U(:, 1) = U_tk_1;
tvect(1) = t;

U_tk = zeros(2*n, 1);

% Move solution forward for each step in time.

num = 1;
for k = 1:m-1
    t = t + dt;

    uhat = fft(U_tk_1);
    what = 1i*wavenu.*(uhat);
    w = (ifft(what));
    U_tk = U_tk_2 - 2*dt*c.*w;
    U_tk_2 = U_tk_1;
    U_tk_1 = U_tk;

    if mod(k, 100) == 0
        U(:, num + 1) = abs(U_tk);
        tvec(num + 1) = t;
        num = num + 1;
    end
end

% Plot the final figures.
figure(2)
plot(x, c, '-')
title("Variable wave speed")
xlabel("x")
ylabel("Variable wave speed c(x)")

figure(3)
waterfall(x, tvec, U')
title("Variation of U with time")
xlabel("x")
ylabel("t")
zlabel("U")

figure(4)
plot(x, U_tk)
title("Solution at t=9")
xlabel("x")
ylabel("U")