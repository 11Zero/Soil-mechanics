function [H,F,E,D] = ouput_HFED(alpha,h,x1)
%% ����HFED����
% alphaΪ�Ƕ�ֵ
m = 1/tan(alpha/180*pi);

% syms y;
% x = x1+0.5*m*h;
f = @(x)([x1+0.5*m*h-x(1),(x(1)-x1-0.5*m*h).^2+(x(2)-0.5*h).^2-1/16*(1+m.^2)*h.^2']);
% [x,y]=
x = fsolve(f,[0,0]);
y = x(2);
y=find(y>0);
x = x(1);
H = [x,y];


% x = x1+0.5*m*h;
% f = @(y)((x-x1-0.5*m*h).^2+(y-0.5*h).^2-25/16*(1+m.^2)*h.^2);
f = @(x)([x1+0.5*m*h-x(1),(x(1)-x1-0.5*m*h).^2+(x(2)-0.5*h).^2-25/16*(1+m.^2)*h.^2']);
x = fsolve(f,[0,0]);
y = x(2);
y=find(y>0);
x = x(1);
F = [x,y];

f = @(x)([-m*x(1)+m*x1+0.5*h*(1+m^2)-x(2),(x(1)-x1-0.5*m*h)^2+(x(2)-0.5*h)^2-1/16*(1+m^2)*h^2])
x = fsolve(f,[0,0]);
y = x(2);
y=find(y>0);
x = x(1);
E = [x,y];


f = @(x)([-m*x(1)+m*x1+0.5*h*(1+m^2)-x(2),(x(1)-x1-0.5*m*h)^2+(x(2)-0.5*h)^2-25/16*(1+m^2)*h^2])
x = fsolve(f,[0,0]);
y = x(2);
y=find(y>0);
x = x(1);
D = [x,y];


function  result = search_best_Fs_trapz(HFDE,o,m,h,x1,gamma_w,gamma__,C__,phi__)
%% ��������������ֵ���ַ�(trapz)�����������Ż�Fsֵ
% HFDE�����ĸ�������
% o��ʾԲ��o
% mΪABб�ʵĵ���
% h��֪
% x1ΪA�������
% gamma_wΪ��w
% gamma__��'
% C__ΪC'
% phi__��'

x0 = o(1);
y0 = o(2);

R = sqrt((x0-x1)^2+y0^2)
x2 = x1+m*h;
y2 = h;
x3 = x0+sqrt(R^2+(y0-h)^2);
y3 = h;

Fs0 = 1;
count = 0;
max_count = 1000;
% ���õ�����������1000
while(true)

	step = 0.1;
	x = x1:step:x3;

	ma = sqrt(R^2-(x-x0)^2)/R+(x-x0)*tan(phi__)/R/Fs0;

	if(x>=x1 && x<x1+m*h)
		Fx=1/m*(x-x1);
	else
		Fx=h;
	end

	Gx = y0-sqrt(R^2-(x-x0)^2);

	Fs = step*R*trapz(C__/ma+tan(phi__)*gamma__*(Fx-Gx)-gamma_w*(Fx-Gx)*tan(phi__)*sqrt(R^2-(x-x0)^2)/R)/trapz((x-x0)*phi__*(Fx-Gx));
	count = count+1;
	if(abs(Fs-Fs0)<1e-3)
        disp('���������ﵽ��������ֹ');
		break;
    end
    
    if(count>max_count)
        disp('����'+max_count+'����δ�ﵽ����������������ֹ');
		break;
    end
    
    Fs0 = Fs;
    
end

result = Fs;


function  result = search_best_Fs_quad(HFDE,o,m,h,x1,gamma_w,gamma__,C__,phi__)
%% ��������simpleson���ַ������������Ż�Fsֵ
% HFDE�����ĸ�������
% o��ʾԲ��o
% mΪABб�ʵĵ���
% h��֪
% x1ΪA�������
% gamma_wΪ��w
% gamma__��'
% C__ΪC'
% phi__��'

x0 = o(1);
y0 = o(2);

R = sqrt((x0-x1)^2+y0^2)
x2 = x1+m*h;
y2 = h;
x3 = x0+sqrt(R^2+(y0-h)^2);
y3 = h;

Fs0 = 1;
count = 0;
max_count = 1000;
% ���õ�����������1000
while(true)

    if(x>=x1 && x<x1+m*h)
        Fs_up = @(x)R*C__/(sqrt(R^2-(x-x0)^2)/R+(x-x0)*tan(phi__)/R/Fs0)+tan(phi__)*gamma__*((1/m*(x-x1))-(y0-sqrt(R^2-(x-x0)^2)))-gamma_w*((1/m*(x-x1))-(y0-sqrt(R^2-(x-x0)^2)))*tan(phi__)*sqrt(R^2-(x-x0)^2)/R;
        Fs_down = @(x)(x-x0)*phi__*((1/m*(x-x1))-(y0-sqrt(R^2-(x-x0)^2)));
    else
        Fs_up = @(x)R*C__/(sqrt(R^2-(x-x0)^2)/R+(x-x0)*tan(phi__)/R/Fs0)+tan(phi__)*gamma__*(h-(y0-sqrt(R^2-(x-x0)^2)))-gamma_w*(h-(y0-sqrt(R^2-(x-x0)^2)))*tan(phi__)*sqrt(R^2-(x-x0)^2)/R;
        Fs_down = @(x)(x-x0)*phi__*(h-(y0-sqrt(R^2-(x-x0)^2)));
    end
    
    Fs = (quad(Fs_up,x1,x1+m*h)+quad(Fs_up,x1+m*h,x3))/(quad(Fs_down,x1,x1+m*h)/quad(Fs_down,x1+m*h,x3));
	
    count = count+1;
	
    if(abs(Fs-Fs0)<1e-3)
        disp('���������ﵽ��������ֹ');
		break;
    end
    
    if(count>max_count)
        disp('����'+max_count+'����δ�ﵽ����������������ֹ');
		break;
    end
    
    Fs0 = Fs;
    
end

result = Fs;


function result = get_reliability_args(phi__v,gamma__v,C__v,sreach_fun_index,HFDE,o,m,h,x1,gamma_w)
%% ���ɿ��Ⱥ�������
% gamma__vΪ��'��ȡֵ������������
% phi__vΪ��'��ȡֵ������������
% C__vΪC'��ȡֵ������������
% sreach_fun_indexΪ���Fs�ĺ���ѡ��
v_len = numel(gamma__v);
Fs_v = zeros(v_len,1);

if(sreach_fun_index==1)
    for i=1:v_len
       Fs_v(i) =  search_best_Fs_trapz(HFDE,o,m,h,x1,gamma_w,gamma__v(i),C__v(i),phi__v(i));
    end
elseif (sreach_fun_index==2)
    for i=1:v_len
       Fs_v(i) =  search_best_Fs_quad(HFDE,o,m,h,x1,gamma_w,gamma__v(i),C__v(i),phi__v(i));
    end
end

% a*X=b,X=a\b
a = [ones(v_len,1) phi__v gamma__v C__v phi__v.^2 gamma__v.^2 C__v.^2];
b = Fs_v; 
X = a\b;
    
result = X;


function result = get_Pf(phi__v,gamma__v,C__v,g_args)
%% ����Pf
% g_args�ں�g(x)�ĸ�������ʽϵ��
mean = mean(phi__v,gamma__v,C__v);
std = std(phi__v,gamma__v,C__v);
x0 = mean;
beta_val0 = 1;
count = 0;
max_count = 1000;
while(true)
    g_divide_x = g_args(2:4)+2.*g_args(5:7).*x0;
    alpha = std.*g_divide_x./sqrt(sum(std.*g_divide_x).^2,2);

    gx = @(beta)g_args(1)+g_args(2)*(mean(1)-beta*alpha(1)*std(1))+g_args(3)*(mean(2)-beta*alpha(2)*std(2))+g_args(4)*(mean(3)-beta*alpha(3)*std(3))...
        +g_args(5)*(mean(1)-beta*alpha(1)*std(1))^2+g_args(6)*(mean(2)-beta*alpha(2)*std(2))^2+g_args(7)*(mean(3)-beta*alpha(3)*std(3))^2;

    beta_val = fzero(gx,beta_val0);
    count = count+1;
    if(abs*(beta_val-beta_val0)<1e-3)
       disp('���������ﵽ��������ֹ'); 
       break;
    end
    
    if(count>max_count)
        disp('����'+max_count+'����δ�ﵽ����������������ֹ');
		break;
    end
    
    x0 = mean-beta_val*alpha.*std;
    beta_val0 = beta_val;
end
result = 1 - normcdf(beta_val0,0,1);


function result = calcu_Ex_cost(x1,x2,x3,alpha,lambda,l,b,n,h,B,Rx_str,Gx_str,S1,S2,S3,S4,Pf)
%% �ӵ������������ɱ�����
% Rx_str,Gx_str�����������ַ�����ʽ���룬��ʽΪf(x)
T = atan(h*tan(alpha)/(h-(n-1)*b*tan(alpha)));
diet_v = n^2*l^2/2*cot(T)+n*(n-1)*l*b/2;
Rx = eval(['@(x)',Rx_str]);
Gx = eval(['@(x)',Gx_str]);
d = x1+n*m*l+(n-1)*b;
a = -d;
% 1�ӵ����ڳɱ�
C1 = (quad(Rx,a,d)-diet_v)*S1;

% 2�ӵ�֧�����̳ɱ�
C2 = ((n-1)*l/sin(T)+(n-2)*b)*S2+(l/sin(T)+b)*S3+B*S4;

% 3�ӵ�����ʧЧ�������޸��ɱ�
C3 = Pf*(diet_v+(x3-x2)*h-quad(Gx,x1,x3))*(S1+S5)+Pf*C2;

% 4�ӵ�����ʧЧ�����Ķ�����ʧ
C4 = lambda*C3;

% 5�ӵ������������̬��ʧ
C5 = (S6+S7)*(quad(Rx,a,d)+13);

C = C1+C2+C3+C4+C5;

result = C;
