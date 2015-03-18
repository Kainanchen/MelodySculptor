function X = stft(x, R, Nfft)
x = x(:).'; 
%% Half Cosine window
%n = (1:R) - 0.5;
%window  = cos(pi*n/R-pi/2);

%% Hanning window
window = hann(R)';
%0.5+0.5*cos(2*pi*n/R);

%% Hamming window
%window = hamming(R)';
%0.54+0.46*cos(2*pi*n/R);

%% Blackman window
%window = blackman(R)';
%0.42-0.5*cos(2*pi*n/R)+0.08*cos(4*pi*n/R);

%% Blackman-Harris window
%n=1:R;
%window= blackmanharris(R)';
%1/R*(0.35875*cos(2*pi*n/R)+0.48829*cos(4*pi*n/R)+0.14128*cos(6*pi*n/R)+0.01168*cos(8*pi*n/R));


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

