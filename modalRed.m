function [A1_, B1_, C1_, D1_] = modalRedLitz(A, B, C, m)
% m - Anzahl modale reduzierte Zustände
% p - Anzahl Eingänge
% q - Anzahl Ausgänge
% n - Anzahl Zustände
[n,p] = size(B);
q = size(C,1);

    % Berechnung Strukturdominanzmaße
    % Trafo auf Diagonalform
    [V, A_] = eigs(A);
    disp(diag(A_));
    B_ = V \ B;
    C_ = C * V;

    % Normierung
    u_j0 = 1;
    G_stat = - C_ / A_ * B_;
    mu = max(abs(G_stat), [], 2);

    % Dominanzwerte
    D = zeros(p, q, n);
    for k = 1:n
        D(:,:,k) = abs((C_(:,k) * B_(k,:) * u_j0) ./ (mu * A_(k,k)));
    end
    S = reshape(sum(D, [1, 2]), 1, n);
    M = reshape(max(D, [],[1,2]), 1, n);

    % nicht dominante EW
    [la2, index2] = mink(S, n-m); 

    % reduziertes System in Modalkoordinaten
    [la1, index1] = maxk(S,m);
    A1_ = A_(index1, index1);
    A2_ = A_(index2, index2);
    B1_ = B_(index1, :);
    B2_ = B_(index2, :);

    %Ausgangsgleichung
    C1_ = C_(:, index1);
    C2_ = C_(:, index2);

    %Verfahren nach Guth
    L = - A_ \ B_;
    L1 = L(index1, :);
    L2 = L(index2, :);
    E_Guth = L2 * pinv(L1);
    C_Guth = C1_ + C2_ * E_Guth;
    C1_ = C_Guth;
    D1_ = zeros(q, p);
end

% %% f) Sprungantwort
% figure
% title('Sprungantworten reduzierte Systeme und Originalsystem')
% step(ss(A1_, B1_, Cr_, zeros(q,p)))
% hold on
% step(ss(A1_, B1_, C_Guth, zeros(q,p)))
% step(ss(A_, B_, C_, zeros(q,p)))
% legend('reduziertes System', 'reduziertes System Guth', 'Originalsystem')