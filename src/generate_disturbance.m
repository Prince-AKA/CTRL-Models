function d = generate_disturbance(t,dt,A1,f1,A2,f2,sigma_d,tau_d)

Nt = length(t);

% Deterministic components
d_det = ...
    A1*sin(2*pi*f1*t) + ...
    A2*sin(2*pi*f2*t);

% Correlated stochastic component
d_colored = zeros(1,Nt);

noise_gain = sqrt(2*sigma_d^2/tau_d);

for k = 2:Nt

    d_colored(k) = ...
        d_colored(k-1) ...
        + dt*(-d_colored(k-1)/tau_d) ...
        + sqrt(dt)*noise_gain*randn;

end

% Total disturbance
d = d_det + d_colored;

end