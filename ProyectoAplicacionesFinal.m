%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aplicaciones de Fisica en Biologia y Medicina I
% Simulacion realista del potencial de
% membrana de una neurona
% Ever Ortega Calderon (everoc.2706@gmail.com) 
% Sebastian Jimenez Carranza (sebas9743@gmail.com)
% Wagner Bermudez Ordonez (webo2900@hotmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametros de la membrana
% Corriente maxima en uA
I_max = 0.9;

% Conductancias maximas en ms/cm2, se eligio 1=K, 2=Na, 3=L
g(1)=36; g(2)=120; g(3)=0.3;
% Area del compartimento en cm2
Area = 1.7795e-3;

% Capacidad en uF/cm2
C=1;
% Potenciales de equilibrio en mV, se eligio 1=n, 2=m, 3=h
E(1)=-12; E(2)=115; E(3)=10.613;
% Condiciones iniciales: V0 en mV, x: 1=n0, 2=m0, 3=h0
I_ext=0.0;
V0=-10;
x=zeros(1,3); x(1) = .318; x(2) = .053; x(3)=.6;

t_inicial=0;
% Pasos en el tiempo 
espaciado_t=0.01;
% Durara 100 ms, para que sea bastante
t_final = 100;
% Implementacion del metodo de Euler para solucionar el sistema
V = V0;
for t=0:espaciado_t:t_final
% Estimulo inicia en t=0ms y termina en 50ms, cuando I_max != 0
% Se divide en dos partes
if t==0; I_ext=I_max/Area; end
if t==50; I_ext=0.0; end
% Se calculan los parametros alfa y beta
alfa(1)=0.01*(-V+10)/((exp((-V+10)/10)-1));
alfa(2)=0.1*(-V+25)/((exp((-V+25)/10)-1));
alfa(3)=0.07*exp((-V)/20);
beta(1)=0.125*exp((-V)/80);
beta(2)=4*exp((-V)/18);
beta(3)=1/(exp((-V+30)/10)+1);
% Se actualiza el vector de parametros n, m y h
x = x.*(1-espaciado_t.*(alfa+beta)) + alfa.*espaciado_t;
% Calcular las conductancias calculadas
gnmh(1)=g(1)*x(1)^4;
gnmh(2)=g(2)*x(2)^3*x(3);
gnmh(3)=g(3);
% Calculo la suma de corrientes finales
I=gnmh.*(V-E);
% Se actualiza el potencial de membrana con el metodo de Euler
V=V+espaciado_t/C*(I_ext-sum(I));
% Se guardan las variables para graficar 
t_inicial=t_inicial+1;
x_plot(t_inicial)=t;
y_plot(t_inicial)=V;
m_plot(t_inicial)=x(2);
n_plot(t_inicial)=x(1);
h_plot(t_inicial)=x(3);
end

% Graficacion del potencial
figure(1);hold off;plot(x_plot,y_plot);
xlabel("Tiempo (ms)"); ylabel("Potencial de membrana (mV)");

% Graficacion de los parametros
figure(2);hold off;plot(x_plot,m_plot,"b");
hold on;plot(x_plot,n_plot,"black");
plot(x_plot,h_plot,"r");xlabel("Tiempo"); ylabel("Conductancias");
legend("m","n","h")