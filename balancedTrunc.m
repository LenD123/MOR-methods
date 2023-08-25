function [A1hat, B1hat, C1hat, D1hat] = balancedTrunc(A,B,C,m)
% m - number of reduced states
% p - number of inputs
% q - number of outputs
% n - number of states
[n,p] = size(B);
q = size(C,1);

% Solve lyapunov equations
P = lyap(A,B*B');
Q = lyap(A', C'*C);

% Cholesky decomp of Q
R = chol(Q);

% Diagonalization
[V,~] = eigs(R*P*R');
Sigma = sqrtm(V'*R*P*R'*V);

% Transformation matrix
T = inv(sqrtm(Sigma)) * V'*R;

% Inverse Transformation matrix
Tinv = inv(T);

% balanced system
Ahat = T*A*Tinv;
Bhat = T*B;
Chat = C*Tinv;

%Phat = TPT', Qhat = Tinv' Q Tinv, Phat = Qhat = Sigma

% reduced system
 A1hat = Ahat(1:m, 1:m);
 B1hat = Bhat(1:m, :);
 C1hat = Chat(:, 1:m);
 D1hat = zeros(q, p); 

end

