%% Задаём поле Z=0
%Парамеры источника:
global Lx Ly D f c
Lx =10;Ly=10; f =10^6; c =1500000;

%Расчёт шага сетки,k,lymda:
global k lymda step
lymda =c/f; k = 2*pi/lymda; step = lymda/6;

%Задаём поле источника
global N M
N = Lx/step;
M = Ly/step;
R = N/2;
IST = zeros(2*N,2*M);
I1 = 1:2*N; 
I2 = 1:2*M;
x = I1-N;                 % mask x-coordinates 
y = I1-M;                 % mask y-coordinates
[X,Y] = meshgrid(x,y);    % create 2-D mask grid
A = (X.^2 + Y.^2 <= R^2); % circular aperture of radius R
IST(A) = 1; 
%pcolor(IST);


%% УГЛОВОЙ СПЕКТР
global z 
z = 20;
hh = (pi/Lx)^2;
%Расчёт углового спектра(БПФ):
%F =(fft2(IST));
F = fftshift(fft2(IST));
%Домножение на член распространения:

for n = 1:2*N
    for m = 1:2*M
%         kx = (n-N)*pi/Lx;
%         ky = (m-M)*pi/Ly;
        S(n,m) = F(n,m)*exp((-1i*z/(2*k))*((n-N)^2+(m-M)^2)*hh);     
    end
end
%   pcolor(abs(S));
%Обратное БПФ:     
 Res =ifft2(S)/4/pi^2;
% pcolor(abs(Res));

%% Обратная задача
Res_0 = conj(Res);
Fobr = fftshift(fft2(Res_0));
%Домножение на член распространения:

for n = 1:2*N
    for m = 1:2*M
        Sobr(n,m) = Fobr(n,m)*exp((-1i*z/(2*k))*((n-N)^2+(m-M)^2)*hh);     
    end
end

%Обратное БПФ:     
Res_obr =ifft2(Sobr)/4/pi^2;
pcolor(abs(Res_obr));


