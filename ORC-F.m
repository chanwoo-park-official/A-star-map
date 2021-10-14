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
x0       = zeros(N+2,1);
x0(1,1)  = 1;
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
    hxk= double(hxk);
else
    hxk = sdpvar(N,N, 'full');
    lambda = double(lambda);
end
% for notational convenience, we shift the indices
lam         = @(i)(lambda(i));
tau         = sdpvar(1);
hx          = @(i,j)(hxk(i,j+1)); % short cut for starting j-indexing of ap at 0


S = (tau) * (x0 * x0') + lam(1) / (2 *L) * g(0) * g(0)' ;

for k = 1 :(N-1)
    S = S - (2*lam(k) - lam(k+1))/(2*L) * g(k) * g(k)';
end

S = S + (2 - 2 * lam(N))/(2*L) * g(N) * g(N)';

for k = 1 :N
    S = S + lam(k)/(2*L) * (g(k-1) - g(k)) * (g(k-1) - g(k))';
end

for k = 1 :(N-1)
    for t = 0 : (k-1)
        coef = lam(k) * hx(k,t) /(2 * L);
        coef2 = 0;
        for j = (t+1):k
            coef2 = coef2 + hx(j,t)/L;
        end
        coef = coef + (lam(k+1) - lam(k))/2 * coef2;
        S = S + coef * ((g(k)* g(t)') + (g(t) * g(k)'));
    end
end

for t = 0:(N-1)
    coef = lam(N) /2 * hx(N,t)/L;
    for j = (t+1): N
        coef = coef + (1-lam(N))/2 * hx(j,t)/L;
    end
    S = S + coef * ((g(N)* g(t)') + (g(t) * g(N)'));
end

for k = 1:(N-1) 
    S = S - (lam(k+1) - lam(k))/2 * (x0 * g(k)' + g(k) *x0');
end

S = S - (lam(1))/2 *  (x0 * g(0)' + g(0) *x0');

S = S - (1 - lam(N))/2 * (x0 * g(N)' + g(N) *x0');

obj_dual = tau;
cons = (S>=0);


if mod(iter, 2) == 1
    cons = cons + (lambda >=0);
end
cons = cons + (tau >= 0);
solver_opt      = sdpsettings('solver','mosek','verbose',verbose);
solverDetails   = optimize(cons,obj_dual,solver_opt);
end

%% Part 4: summary

sxk_real = zeros(N,N);
syk_real = zeros(N,N);
szk_real = zeros(N,N);
phi = zeros(N+3,1);
phi_zero = 0;
phi(1) = 2;
for k = 2:(N+3)
    phi(k) = (phi(k-1) + 1) + sqrt(phi(k-1) + 1);
end
for k = 1:N
    if k ==1
        syk_real(1,1) = 1/L;
        szk_real(1,1) =(phi(1) - phi_zero)/L;
        sxk_real(1,1) =  phi(1)/phi(2) * syk_real(1,1)+ (phi(2) - phi(1) )/ phi(2) *szk_real(1,1);
    else
        syk_real(k,: ) = sxk_real(k-1,:);
        syk_real(k, k) = 1/L;
        szk_real(k, :) = szk_real(k-1, :);
        szk_real(k, k) = (phi(k) - phi(k-1))/L;
        sxk_real(k, :) = phi(k)/phi(k+1) * syk_real(k,:)+ (phi(k+1) - phi(k)) / phi(k+1) *szk_real(k,:);
    end
end

hxk_real = zeros(N,N);
hxk_real(1, :) = sxk_real(1, :);
for k = 2:N
    hxk_real(k, :) = sxk_real(k, :) - sxk_real(k-1, :);
end
% compare two answer
hxk_real- double(hxk)
