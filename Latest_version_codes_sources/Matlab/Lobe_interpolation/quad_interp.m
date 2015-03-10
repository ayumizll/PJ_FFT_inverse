function y = quad_interp(xp,yp,xq)
Np = length(xp);
Nq = length(xq);
u = (1:Np)';

y = zeros(Nq,1);

A = xq(:,ones(1,Np));
A = A';

B = xp(:,ones(1,Nq));

C = A-B;
C_abs = abs(A-B);

[~,C_tab_min] = min(C_abs,[],1);

% Evaluate difference where the minimum is reached
idx = sub2ind(size(C), C_tab_min, 1:Nq);
min_values = C(idx);

% Keep the closest value towards left
C_tab_min(min_values < 0) = C_tab_min(min_values < 0) - 1; 

C_tab_min(C_tab_min >= (Np-1)) = Np - 2; 

C_tab_min(C_tab_min <= 0) = 1;

[a,b,c] = quad_coeffs(xp(C_tab_min),xp(C_tab_min+1),xp(C_tab_min+2),yp(C_tab_min),yp(C_tab_min+1),yp(C_tab_min+2));

y = a.*xq.*xq + b.*xq + c;

y(xq < min(xp)) = 0;
y(xq > max(xp)) = 0;
end

