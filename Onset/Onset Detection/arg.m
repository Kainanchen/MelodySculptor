function y = arg(x)
if real(x)~=0
y=atan(imag(x)/real(x));
else
    y=0;
end