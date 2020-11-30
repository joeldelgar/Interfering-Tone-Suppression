function [y] = Filtre_to_Interferent(x,fs,r)

L=length(x);
N=2^nextpow2(L);

X=fft(x,N);
[val,pos]=max(abs(X(1:N/2)));
omega=(2*pi*pos/N);
Num=poly([exp(1i*omega),exp(-1i*omega)]);
Den=poly([r*exp(1i*omega),r*exp(-1i*omega)]);

h=impz(Num, Den, 100);

%y = conv(x,h);
Y=X.*fft(h,N);
y=ifft(Y,N);
end

