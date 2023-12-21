L = 2.*pi;
T = 5.;
n = 512;

dx = L/(2.*n);
dt = dx/4.;
x = linspace(0,L-dx,2*n)';
c = 1 + 0*x;

wavenu = [0:n -n+1:-1]';
U_tk_1 = zeros(2*n,1);
for i = 1:2*n
    if x(i) >= 0.5
        if x(i) <= 1.
            U_tk_1(i) = 1;
        end
    end
end

U_tk_2 = zeros(2*n,1);
for i = 1:2*n
    if x(i) - dt >= 0.5
        if x(i) - dt <= 1.
            U_tk_2(i) = 1;
        end
    end
end

figure(1)
plot(x, U_tk_1, '-')
axis ([0 2*pi 0 2])
title("Initial condition for U")
xlabel("x")
ylabel("U")

% Store the solution at every 100th timestep to plot later.
m = round(T/dt);
U = zeros(2*n, round(m/400));
tvec = zeros(round(m/400), 1);

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

    if mod(k, 400) == 0
        U(:, num + 1) = abs(U_tk);
        tvec(num + 1) = t;
        num = num + 1;
    end
end

% Plot the final figures.
figure(2)
plot(x, c, '-')
title("Constant wave speed c = 1")
xlabel("x")
ylabel("Constant wave speed c = 1")

figure(3)
waterfall(x, tvec, U')
axis ([0 1 0 2 0 1.1])
title("Variation of U with time")
xlabel("x")
ylabel("t")
zlabel("U")

figure(4)
plot(x, U_tk)
axis ([0 2 0 1.1])
title("Solution at t=5")
xlabel("x")
ylabel("U")