function lam = lambert(A0,A1,tau)

a1  = eig(A1);
a1  = a1(:);
a0  = diag(A0);
for k = -1:10
    for i = 1:numel(a1)
        lam(i,k+2) = lambertw(k,tau*a1(i)*exp(-a0(i)*tau))/tau+a0(i);
    end
end
lam = [lam conj(lam)].';
