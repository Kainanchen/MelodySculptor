clear
close all

[source,fs] = audioread('SheIsMySin.wav');      % Read source
x = source(:,1);                                % Select a channel
L = size(x,1);                                  % Audio length
N = 1024;                                       % N of FFT
N2 = ceil((N+1)/2);                             % Half N
H = floor(N/2);                                 % Hop size
F = ceil((L-N)/H);                              % No. of frames
X = spectrogram(x,hann(N),N-H,N,fs);            % STFT
M = abs(X);                                     % Spectrum Magnitude
P = unwrap(angle((X)));                         % Spectrum Phase
SX = M./repmat(sum(M),N2,1);
WX = (SX(:,1:end-2)+SX(:,2:end-1)+SX(:,3:end))/3;
%% Energy
DSE = sum(M(:,2:end)-M(:,1:end-1),1);           % DSE
TM = 0.2*sum(M(:,1:end-1),1);                   % Magnitude threshold
MarkE = DSE>TM;                                 % Mark

%% Phase
% PDA = P(:,3:end)-2*P(:,2:end-1)+P(:,1:end-2);
% PD = abs(sum(PDA.*WX(:,2:end-1),1));
% TP = 35*pi;
% MarkP = PD(2:end)>TP;


%% Complex domain
delta=1;
lamda=2;
EstP=princarg(2*P(:,2:end-1)-P(:,1:end-2));
EstSig=M(:,2:end-1).*exp(1j*EstP);
ComDist=sum(sqrt((real(EstSig)-real(X(:,3:end))).^2+(imag(EstSig)-imag(X(:,3:end))).^2).*WX.^0.5,1);

threshold=zeros(1,F-H-1);
NoT=length(threshold);
for i =1:NoT;
    threshold(i)=delta+lamda*median(ComDist(i:i+H-1)); %moving median solvable for matrix?
end
MarkC=ComDist(H/2:end-H/2)>threshold(1:end);
MarkC=[zeros(1,H/2),MarkC,zeros(1,H/2)];

% stem(MarkP)
% hold on
% stem(MarkE(2:end-1),'r')
% hold on
stem(MarkC,'b')
