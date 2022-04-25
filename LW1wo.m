clc
clear all
close all

 L = 1;%[0,1]
 delta_x = 0.001*L ;%mesh makani(mitavand riztar ham shavad.)

 
 a = 0.5 ;%Defined, but can be changed
  
 x = 0.0 : delta_x : L ;
  X = length(x) ;

delta_t=[0.001,0.0005,0.0001,0.00005,0.00001];  %%riz kardan mesh(zamani)
col=['r','k','c','g','m'];
error = zeros(5,X) ;
Eu=zeros(5,X);
Eanalytical=zeros(5,X);

for i = 1 : length(delta_t)
    
  %dt and dx 
  x = 0.0 : delta_x : L ;
  X = length(x) ;
  dt = delta_t(i) ;
  total_time = 100 ;
  
  %pre defining matrix
  u = zeros(X,1) ;
  analytical = zeros(X,1) ;
%   error = zeros(5,X) ;
    
 
  %nu=alfa*dt/dx
  nu = [a*delta_t(1)/delta_x, a*delta_t(2)/delta_x, a*delta_t(3)/delta_x, a*delta_t(4)/delta_x, a*delta_t(5)/delta_x] ; 
  
  for ix = 1 : X
      if (x(ix)<0.25)
          u(ix) = 1.0 ;
      else 
          u(ix) = 0.0 ;
      end
      
     
  end
  
  
  for it = 1 : total_time
u(1) = 1 ;
      for ix = 2 : X-1
        u(ix) = u(ix)-(a*0.5*dt/delta_x)*((u(ix+1)-u(ix-1))) + 0.5*a^2.0*dt^2.0*((u(ix+1)-2.0*u(ix)+u(ix-1)))/(delta_x^2) ;
        Eu(i,ix)= u(ix);
      end  
     % u0 = u ;
    
      
        for ix = 1 : X
      if ( x(ix)<0.25)
          analytical(ix) = 1.0 ;
      else 
          analytical(ix) = 0.0 ;
      end
      Eanalytical(i,ix)= analytical(ix);
  end
      
  end
  
    hold on
    scatter (x,u,col(i),'linewidth',1.2)
    axis([0 1 0 1])
    
    
    xlabel 0<x<1
    ylabel u
    title 'Lax-Wendroff:t=100'
    
    for j=1:length(delta_t)
        for ix = 1:X
  error(j,ix) = abs(Eu(j,ix)-Eanalytical(j,ix));
        end
    end
end

plot(x,analytical,'linewidth',1.5)
xlim([0 1])
ylim([-0.4 1.4])
legend('\nu = 0.5','\nu = 0.25','\nu = 0.05','\nu = 0.025','\nu = 0.005', 'analytical')

%


noo1 = [0.5,0.25,0.05,0.025,0.005] ;



%


figure(3)
plot(x,error,'linewidth',1.5)
legend('\nu = 0.5','\nu = 0.25','\nu = 0.05','\nu = 0.025','\nu = 0.005', 'analytical')
title 'abs(Error)'
xlim([0 1])
ylim([-0.4 1.4])
  xlabel 0<x<1
  ylabel |error|
  
  
  noo1 = [0.5,0.25,0.05,0.025,0.005] ;
for p = 1 : length(noo1)
  for ix = 1 : X
            k(ix) = 2.0*pi*ix/(2.0*L) ;
            beta(ix) = k(ix)*delta_x ;
%             amplitude(ix) = noo1(p)^2*cos(beta(ix))/2+(1-noo1(p)^2) ;
%             amplitudey(ix)=-noo1(p)/2*sin(beta(ix))
            amplitude(ix) =0.5* noo1(p)^2*(cos(beta(ix)))+1-0.5* noo1(p)^2
            amplitudey(ix)=0.5* noo1(p)*sin(beta(ix)) ;
%             
      end
    
    figure(2)
    hold on
    grid on
    plot (amplitude,amplitudey,col(p),'linewidth',1.5)
    legend('\nu = 0.5','\nu = 0.25','\nu = 0.05','\nu = 0.025','\nu = 0.005')
   
    xlabel 'real part'
    ylabel 'imaginary part'
    title 'G=(1+\nu^2/2cos(\beta)-\nu^2/2-i \nu/2 sin(\beta))'
end
