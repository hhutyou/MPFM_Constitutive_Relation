%% This code is developed to give the analytical solution of our 
%% micromechanics-based phase-field model, see You et al., 2021 (CMAME). 
 % Stress state: Conventional triaxial compressive stress
 % Material: Luc du Bonnet granite, see Martin, 1997
clear all
% close all
clc
syms d lambda_p sigma1 sigma2 X g l A k
f1=2*(1-d)*(1/2*X*lambda_p^2)-g/l*d*k*(1-(1-d)^2)^2==0; % Damage evolution law, i.e., Eq. (55) in You et al., 2021, CMAME
d2lambda_p=solve(f1,lambda_p); % lambda_p = f(d)
fg=str2func(['@(d,X,g,l,k)',vectorize(d2lambda_p(1))]); % Construct function handle from function "d2lambda_p".
% The following is the yield function under Conventional triaxial compression, i.e., Eq. (56) in You et al., 2021, CMAME
% f2 = sigma1-(sqrt(6)+2*A)/(sqrt(6)-A)*sigma2+3*sqrt(X*g/l*d*(1-d)^3/k)/(sqrt(6)-A)==0; 
f2 = sqrt(2/3*(sigma1-sigma2)^2)+A/3*(sigma1+2*sigma2)-sqrt(X*g/l/k*d*(1-d)^3)==0; 
d2sigma1=solve(f2,sigma1); % sigma1 = f(d)
f=str2func(['@(d,sigma2,A,X,g,l,k)',vectorize(d2sigma1(1))]); % Construct function handle from function "d2sigma1".
%% Initiation of material parameters. 
i=0;
E=63000;v=0.15;k0=E/3/(1-2*v);u0=E/2/(1+v);
sigma2=-10;A=1.7953;k=10;
X=A^2*k0+2*u0;
l=1; % The parameter l takes no effect on the constitutive behavior. Take a random value.
g=(48.8*16/3/sqrt(3))^2*l*k/X;
%% Load by increasing d
for d=0.0001:0.001:1
    i=i+1;
    s1(i)=f(d,sigma2,A,X,g,l,k); % Axial stress
    s2=sigma2; % Confining pressure
%     s2=s1(i)/k;
%     N=1/2/sqrt((2*(s1(i)-s2)^2)/6)*((2*s1(i)-2*s2)/3)+n/3
    N=-2/sqrt(6)+A/3; % Axial plastic flow direction, see Eq.(40) and Eq.(51)
    e1(i)=fg(d,X,g,l,k)*N+s1(i)/E-v/E*(2*s2); % Axial strain
    N3=1/sqrt(6)+A/3; % Normal plastic flow direction, see Eq.(40) and Eq.(51)
    e3(i)=fg(d,X,g,l,k)*N3+s2/E-v/E*(2*s1(i)); % Normal stain, see Eq.(40) and Eq.(51)
%     fg(d,H,g,l,ka,yita)
    y(i)=d;
    x(i)=e1(i);
end
%% Experimental data
% hd=parula(5);
% % Load data
% load('p_4.txt');  load('p_10.txt'); load('p_15.txt'); load('p_20.txt'); load('p_30.txt');
% data1=p_4;
% num1=size(data1,1);
% s3=zeros(num1,1);s1=zeros(num1,1);
% s3(:,1)=-data1(:,1)/100;
% s1(:,1)=-data1(:,2);%%data transmission accomplished
% plot(s3,s1,'Marker','o','MarkerEdgeColor',hd(1,:),'LineStyle','none');hold on
%%
color=parula(5);
figure(1)
plot([0,e1]*100,[0,s1],'color',color(1,:),'DisplayName','$Pc=-10$ MPa');hold on
%  ylabel
ylabel({'Axial stress $\sigma_1$ [MPa]'},'Interpreter','latex');
%  xlabel
xlabel({'Axial strain $\varepsilon_1$ [\%]'},'Interpreter','latex');
xlim([-1 0])
set(gca,'XDir','reverse','YDir','reverse')
% legend
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex');
% figure(2)
% plot([0,x]*100,[0,y]);hold on
% xlim([-1 0])
% set(gca,'XDir','reverse')
