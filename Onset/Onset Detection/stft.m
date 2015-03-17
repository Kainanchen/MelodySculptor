function X = stft(x, R, Nfft)
x = x(:).'; 
n = (1:R) - 0.5;
window  = cos(pi*n/R-pi/2);
x = [zeros(1,R) x zeros(1,R)];

Nx = length(x);
Nc = ceil(2*Nx/R)-1;       
L = R/2 * (Nc + 1);
if Nx < L
    x = [x zeros(1,L-Nx)]; 
end
X = zeros(R,Nc);
i = 0;
for k = 1:Nc
    X(:,k) = window .* x(i + (1:R));  
    i = i + R/2;
end
X = fft(X,Nfft);   

