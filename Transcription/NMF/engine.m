clear
close all

[source,fs] = audioread('../../Audio/SheIsMySin.wav');      % Read source
x = source(:,1);                                % Select a channel

L = size(x,1);                                  % Audio length
blockLen = fs*5;                                % Block length = 5s
block = x(1 : blockLen);

N = 1024;                                       % N of FFT
N2 = ceil((N+1)/2);                             % Half N
H = floor(N/2);                                 % Hop size
FH = 2^(nextpow2(0.2*fs/H)-1);                  % Frame Hop size
F = ceil((blockLen-N)/H);                       % No. of frames

X = spectrogram(block,hann(N),N-H,N,fs);        % STFT
M = abs(X);                                     % Spectrum Magnitude
P = unwrap(angle(X));                           % Spectrum Phase
