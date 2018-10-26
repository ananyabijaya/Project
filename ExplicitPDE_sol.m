%% Numerical modelling and simulation project
clc
clear
%Given
L=10;
b=.4;
d=.4;
A=b*d;
rho=7860;
E=200*10^9;
F=@(t) t*0
W=10^5;
I=b*d^3/12;
t_steps=100000;
T_end=5;
%Keep x_steeps even number.
x_steps=46;
x=linspace(0,L,x_steps+1)
%Deflection across the length of the beam is given by w(x,t)
w=zeros(t_steps+2,x_steps+3);

%Response at initial time t=1
w_init=@(x) -(W.*x/(12*E*I)).*(3*L^2/4-x.^2);
x_half=x(1:(1+size(x,2))/2);
w_half=w_init(x_half);
%using symmetry to calculate response from center to the end of the beam
w_N_half=fliplr(w_half);
w_N_half=w_N_half(2:end-1);
w(2,2:end-1)=[0 w_half(2:end) w_N_half 0];

delta_x=L/x_steps;
delta_t=(T_end-1)/t_steps;

m=2 %First time step
%Initialing dummy nodes before x=1 and after x=L
w(m,1)=2*w(m,2)-w(m,3);
w(m,end)=2*w(m,end-1)-w(m,end-2);
    for i=3:x_steps+1
        w(m-1,i)=w(m,i);
        w(m+1,i)=delta_t.^2*F(m)/(rho*A)-(E*I*delta_t^2/(delta_x^4*rho*A))*(w(m,i-2)-4*w(m,i-1)+6*w(m,i)-4*w(m,i+1)+w(m,i+2))+2*w(m,i)-w(m-1,i);
        i=i+1;
    end
 

for m=3:t_steps-1
    w(m,1)=2*w(m,2)-w(m,3);
    w(m,end)=2*w(m,end-1)-w(m,end-2);
    for i=3:x_steps+1
        w(m+1,i)=delta_t.^2*F(m)/(rho*A)-(E*I*delta_t^2/(delta_x^4*rho*A))*(w(m,i-2)-4*w(m,i-1)+6*w(m,i)-4*w(m,i+1)+w(m,i+2))+2*w(m,i)-w(m-1,i);
        i=i+1;
    end
    m=m+1;
end

figure; hAxes = gca;
set(gca, 'FontSize', 16)
hold( hAxes, 'on' )
i=2
while i~=100000
   plot( x, w(i,2:end-1),'linewidth',1);
   title('Deflection response of the beam')
   xlabel('Length (m)')
   ylabel('Deflection') 
   
   if i==2
   i=5000
   else
       i=i+5000
    end
end

