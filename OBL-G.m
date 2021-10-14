clc; clear all;

%% Part 1: setting upparameters and notations
% PARAMETERS TO BE TUNED
N       = 5;  
L       = 1;
%l       = 3;   %want l th energy and l+1 th energy comparison
% END OF TUNABLE ZONE
verbose = 0;
gk       = zeros(N+2,N+1); % each column is a g_i: gk = [g_0 g_1 ... g_{N}]
gk(2:end,:) = eye(N+1);
fk       = eye(N+1);       % each column is a f_i: fk = [f_0 f_1 ... f_{N}]
g        = @(i)(gk(:,i+1));% short cut for starting i-indexing of g at 0
f        = @(i)(fk(:,i+1));% short cut for starting i-indexing of f at 0
G   = sdpvar(N+2,N+2,'symmetric');  %Gram-schdmit matrix, positive definite
F   = sdpvar(N+1,1);         %F matrix, all are positive, f(x)

%% Part 3: dual problem and alternative optimization

% PARAMETERS TO BE TUNED
lambda = zeros(1,N);
beta = zeros(1,N);
hxk = zeros(N,N);

for iter = 1:1000
if mod(iter, 2) ==1
    lambda = sdpvar(1,N);  
    beta= sdpvar(1,N);  
    hxk= double(hxk);
else
    hxk = sdpvar(N,N, 'full');
    lambda = double(lambda);
    beta = double(beta);
end
% for notational convenience, we shift the indices
lam         = @(i)(lambda(i));
bet         = @(i)(beta(i+1));
tau         = sdpvar(1);
hx          = @(i,j)(hxk(i,j+1)); % short cut for starting j-indexing of ap at 0


S = (tau / (2*L) - 1 ) *g(N) * g(N)';
for k = 1 :N
    coef = 0;
    for t = 0: (k-1)
        coef = coef + hx(k,t) / L * g(t);
    end
    S = S + lam(k)/2 * (g(k) * coef' + coef * g(k)' + 1/L * (g(k-1) - g(k)) * (g(k-1) - g(k))');
end

for k = 0:(N-1)
    coef = 0;
    for j = (k+1):N
        for t = 0:(j-1)
            coef = coef + hx(j,t)/L * g(t);
        end
    end
    S = S + bet(k)/2 * (-coef * g(k)' - g(k)* coef');
end

obj_dual = tau;
cons = (tau - lam(1) + bet(0) == 0);


if mod(iter, 2) == 1
    for k = 1:(N-1)
        cons = cons + (lam(k) - lam(k+1) + bet(k) == 0);
    end
    cons = cons + (lambda >=0);
	cons = cons + (beta >=0);
end
cons = cons + (S>=0);
cons = cons + (tau >= 0);
solver_opt      = sdpsettings('solver','mosek','verbose',verbose);
solverDetails   = optimize(cons,obj_dual,solver_opt);
end

%% Part 4: summary

hxk_real = zeros(N,N);
for i = 1:N
    hxk_real(i,i) = 3 * (N-i+1)/(N-i+3);
end

for i = 1:N
    for j = 0:(i-2)
        hxk_real(i, j+1) = 2 * ((N-i) * (N-i+1) * (N-i+2))/ ((N-j) * (N-j+1) * (N-j+2));
    end
end

alpha = 1/2 * ( 1 + sqrt(N * (N+1)/2));

h0_real_add =zeros(N);
h0_real_add(1) = alpha * 4/(N+2) + (N-2) / (N+2);
for k = 2:N
    h0_real_add(k) = alpha * 4/(N+3-k) + (N-1-k) / (N+3-k) * h0_real_add(k-1);
end

hxk_real(1,1) = h0_real_add(1);

for i = 2:N
    hxk_real(i, 1) = h0_real_add(i) - h0_real_add(i-1);
end
% compare two answer
hxk_real- double(hxk)
