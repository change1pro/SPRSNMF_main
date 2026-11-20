function F = update_f_splagnmf(X, B, W, Ds, d, F, alpha, beta, innerF, epsv)
Q = diag(1 ./ (4 * (sum(F.^2, 2) + epsv).^(3/4)));
BtB = B' * B;
DXtB = (d .* X') * B;

for t = 1:innerF
    Gp = (d .* F) * BtB + 2 * alpha * Ds * F + beta * Q * F;
    Gm = DXtB + 2 * alpha * W * F;

    Psi = F' * (Gm - Gp);
    Psip = (abs(Psi) + Psi) / 2;
    Psim = (abs(Psi) - Psi) / 2;

    numerator = max(Gm + F * Psim, epsv);
    denominator = max(Gp + F * Psip, epsv);
    F = F .* sqrt(numerator ./ max(denominator, epsv));
    F = max(F, epsv);
    F(~isfinite(F)) = epsv;
    F = F ./ max(sqrt(sum(F.^2, 1)), epsv);

    Q = diag(1 ./ (4 * (sum(F.^2, 2) + epsv).^(3/4)));
end
end
