function [A1hat, B1hat, C1hat, D1hat] = Krylov(A, B, C)
q = size(C,1);

%     span(A^-1*B, A^-1*A^-1*B, A^-1*A^-1*A^-1*B) 
%     V = orth([A\B, inv(A)*inv(A)*B, inv(A)*inv(A)*inv(A)*B]);

    % Arnoldi Verfahren
    V = A\B/norm(A\B);
    for i = 2:q
        v = A \ V(:, end);
        for j = 1:i-1
           h = v' * V(:,j);
           v = v - h * V(:,j);
        end
        if v == 0
            break
        else
            V = [V, v / norm(v)];    
        end
    end
    W = V;
    
    A1hat = W' * A * V;
    B1hat = W' * B;
    C1hat = C * V;
    D1hat = zeros(size(C, 1), size(B,2));
end