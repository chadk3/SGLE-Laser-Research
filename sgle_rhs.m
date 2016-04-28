
% rhs for the Sinusoidal Ginzburg Landau Equation

function rhs = sgle_rhs(z, psit, dummy, kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau)

psi = ifft(psit);

a=trapz(t,abs(psi).^2);

% gz is the saturating gain
gz = 2*g0/(1 + (1/e0 * a));

% In is the total power of the field at the beginning of the fiber section
In = abs(psi).^2;

w = B * In * sin(2*(alpha1 - alphap));

Q = 1/2 * (exp(-1i * K) * (cos(2 * alpha2 - 2* alpha3 - alphap) + 1i * cos(2 * alpha3 - alphap)) * (1i * cos(2 * alpha1 - alphap - w) - cos(alphap - w)) + exp(1i * K) * (sin(2 * alpha2 - 2 * alpha3 - alphap) - 1i * sin(2 * alpha3 - alphap)) * (sin(alphap - w) - 1i * sin(2 * alpha1 - alphap - w)));
LQ = log(Q);

% Original rhs
rhs = -1i*(D/2 - 1i * gz * tau)*(kt.^2).*psit - 1i* fft( 1i*gz *psi + 1i*psi.* log(Q) - psi.*(abs(psi).^2)-1i*gamma*psi);

%rhs = -1i*(D/2 - 1i* gz * tau) * (kt.^2).*psit + gz*psit - gamma*psit - 1i* fft(1i*psi.*log(Q) - psi.*(abs(psi).^2));

%rhs from the rewritten code
%rhs = (kt.^2).* (-1i * D/2 - gz * tau)* psit + (gz - gamma) * psit + fft(1i * (psi.* abs(psi).^2) + psi.* LQ);

%rhs = (-1i * D/2 - gz * tau + gz - gamma) * (kt.^2).*psit + fft(1i * (psi.* abs(psi).^2) + psi.* LQ);

