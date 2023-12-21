% Set the values for the KDV equation.
A = 100;
V = A/3.;
L = 2.*pi;
T = 2.*pi/V;
n = 128;

dx = L/(2.*n);
dt = (dx^3)/(2*pi^2);
m = round(T/dt);
t = 0.;
x = linspace(0,L-dx,2*n)';

wavenu = [0:n -n+1:-1]';
U_tk_1 = A*sech(sqrt(A/12.)*(x-3*pi/2)).^2 + 2*A*sech(sqrt(2.*A/12.)*(x-pi/2)).^2;
U_tk_2 = A*sech(sqrt(A/12.)*(x-3*pi/2 + V*dt)).^2 + 2*A*sech(sqrt(2.*A/12.)*(x-pi/2 + V*dt)).^2;

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

    U_tk = U_tk_2 - 2*dt*U_tk_1.*w - 2.*www;
    U_tk_2 = U_tk_1;
    U_tk_1 = U_tk;

    if mod(k, 10000) == 0
        U(:, num + 1) = U_tk;
        tvec(num + 1) = t;
        num = num + 1;
    end
end

% Plot the final figures.
figure(2)
waterfall(x, tvec, U')
title("Variation of U with time")
xlabel("x")
ylabel("t")
zlabel("U")

val = length(U); 
value = input("Do you want a graph at a specific point in time? (No (0) or Yes (1)) ");
if value == 1
    k = round(input("At what value of t from 1 to val? "));
    figure(3)
    hold on
    vect = U(:,k)';
    plot(x, vect, "r")
    title("Plot at the given step in time")
    xlabel("x")
    ylabel("U")
end

figure(4)
plot(x, U_tk)
title("Solution at end")
xlabel("x")
ylabel("U")