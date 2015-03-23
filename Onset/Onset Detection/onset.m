clear
close all
source=audioread('sheismysin.wav');
 %source=audioread('darkchestofwonders.wav');
x = source(:,1)'; 
%y = source(:,2)'; second track
N = length(x);

%% Compute STFT

R = 512;                
Nfft = 1024;       
X = stft(x, R, Nfft); 
    
%% Energy
m=length(X(1,:));
dse=zeros(1,m);
marke=zeros(1,m);
marked=zeros(1,m);
for i=2:m
    
dse(i)=sum(abs(X(:,i)-abs(X(:,i-1))));
delta=0.3*sum(abs(X(:,i-1)));  %assumed threshold
if dse(i)>delta
    marke(i)=1;
end
end


%% Phase
dsp=zeros(1,m);
markp=zeros(1,m);
markpd=zeros(1,m);
deltap=0.2*pi;  

for i=1:m-1
    uX(:,i)=abs(X(:,i))/(sum(abs(X(:,i)))); % weighting matrix
if i<=3
    dsp(i)=0;
else
PhaseArray=unwrap([angle(X(:,i-2)),angle(X(:,i-1)),angle(X(:,i))]);
    pf=PhaseArray(:,3)-2*PhaseArray(:,2)+PhaseArray(:,1);
    dsp(i)=sum(abs(princarg(pf'*uX(:,i))));
    if dsp(i)>deltap
        markp(i)=1;
    end
end
display(i);
end

%% Complex domain
disc=zeros(m);
markc=zeros(m);
H=2*ceil(m/100); % H introduces delay
for i=1:m-1
    uX(:,i)=abs(X(:,i))/(sum(abs(X(:,i))));
    if i<=3
        dsc(i)=0;
    else

        PhaseArray=unwrap([angle(X(:,i-1)),angle(X(:,i))]);
        ep=princarg(2*PhaseArray(:,2)-PhaseArray(:,1));
        amp=abs(X(:,i-1));
        estsig=amp'*exp(j*ep);
        disc(i)=sum(sqrt(real((estsig-X(:,i))).^2+(imag(estsig-X(:,i))).^2)); 
    end
        
deltac=5; % random constant?
lamda=1;
if i<H % assuming the first non-zero index is smaller than half width of median filter
    markc(i)=0;
else
    eta=median(disc(i-H/2:i+H/2));
    if disc(i-H/2)>=deltac+lamda*eta
        markc(i-H/2)=1;
    else
        markc(i-H/2)=0;
    end
end
display(i);
end;


        
        
    
