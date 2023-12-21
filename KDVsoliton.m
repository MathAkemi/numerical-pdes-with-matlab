% Set the values for the KDV equation.
A = 100;
V = A/3.;
L = 2.*pi;
x1 = L;
T = 80.*pi/V;
n = 128;

dx = L/(2.*n);
dt = T/320000; %dx^3/(2*pi^2);
m = 40; %80;
t = 0.;
x = linspace(0,x1-dx,2*n)';

wavenu = [0:n -n+1:-1]';
U_tk_1 = A*sech(sqrt(A/12.)*(x-pi)).^2;
U_tk_2 = A*sech(sqrt(A/12.)*(x + V*dt - pi)).^2;

figure(1)
plot(x, U_tk_1, '-')
title("Initial condition for U")
xlabel("x")
ylabel("U")

% Store the solution at every 100th timestep to plot later.
U = zeros(2*n, round(m/100));
tvec = zeros(round(m/100), 1);

U(:, 1) = U_tk_1;
tvec(1) = t;
xi = linspace(0,x1-dx,2*n)';

U_tk = zeros(2*n, 1);

% Move solution forward for each step in time.

num = 1;
for k = 1:m-1
    t = t + dt;

    uhat = fft(U_tk_1);
    what = 1i*wavenu.*(uhat);
    w = real(ifft(what));
    what = -1i*sin(dt*wavenu.^3).*(uhat);
    www = real(ifft(what));

    U_tk = U_tk_2 - 2*dt*(U_tk_1.*w + www);
    U_tk_2 = U_tk_1;
    U_tk_1 = U_tk;

    U(:, num + 1) = U_tk;
    tvec(num + 1) = t;
    num = num + 1;
end



% Plot the final figures.
figure(2)
waterfall(x, tvec, U')
title("Variation of U with time")
xlabel("x")
ylabel("t")
zlabel("U")

figure(3)
plot(x, U_tk)
title("Solution at end")
xlabel("x")
ylabel("U")