function B = update_B(X, F, d, B, innerIter, epsv)
%UPDATE_B Multiplicative update for basis matrix B.

DF = d .* F;
for t = 1:innerIter
    numerator = X * DF;
    denominator = B * (F' * DF);
    B = B .* sqrt(numerator ./ max(denominator, epsv));
    B = max(B, epsv);
end
end
