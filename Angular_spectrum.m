%% Задаём поле Z = 0
%Парамеры источника:
global Lx Ly D f c
Lx =20; Ly=20; f =10^6; c =1000000;


%Расчёт шага сетки,k,lymda:
global k lymda step
lymda =c/f; k = 2*pi/lymda; step = lymda/2;


%Задаём поле источника
global N M
N = Lx/step;
M = Ly/step;
R = N/2;
IST = zeros(2*N+1,2*M+1);
I1 = 1:2*N+1; 
I2 = 1:2*M+1;
x1 = -Lx:step:Lx;
x2 = -Ly:step:Ly;
x = I1-N-1;                 % x-координаты
y = I1-M-1;                 % y-координаты
[X,Y] = meshgrid(x,y);      % задаём плоскую сетку
A = (X.^2 + Y.^2 <= R^2);   % апертура радиуса R
IST(A) = 1; 


%Начальное распределение
tiledlayout(1,1)
nexttile
pcolor(IST); 
title('Поршневой источник')

%% УГЛОВОЙ СПЕКТР
global z 
z = 10;
hh = (pi/Lx)^2;
%Расчёт углового спектра(БПФ):
F = fftshift(fft2(IST));


%Домножение на член распространения:
for n = 1:(2*N+1)
    for m = 1:(2*M+1)
        kx = (n-N-1)*pi/Lx;
        ky = (m-M-1)*pi/Ly;
        if ((kx^2+ky^2) > k^2 ) 
            kz = k;
            F(n,m) = 0;
        else
            kz = sqrt(k^2-kx^2-ky^2);
        end
        prom = 1i*z*(kz)*0,1;
       S(n,m) = F(n,m)*exp(prom);

    end
end


% Обратное БПФ:     
 Res =ifft2(S);



%% ОБРАТНАЯ ЗАДАЧА
% Res_0 = conj(Res);
Res_0 = Res;
Fobr = fftshift(fft2(Res_0));
% pcolor (abs(Fobr));


%Домножение на член распространения:
for n = 1:2*N+1
    for m = 1:2*M+1
        kx = (n-N-1)*pi/Lx;
        ky = (m-M-1)*pi/Ly;
        if ((kx^2+ky^2) > k^2 ) 
            kz = k;
%             Fobr(n,m) = 0;
        else
            kz = sqrt(k^2-kx^2-ky^2);
        end
        prom2 = -1i*(kz)*z*0,1;
        Sobr(n,m) = Fobr(n,m)*exp(prom2);
    end
end


%Обратное БПФ:     
Res_obr =(ifft2(Sobr));
nexttile
pcolor(abs(Res_obr));
title('Восстановленный источник')



