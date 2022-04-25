clc
clear all
close all

 L = 1;%[0,1]
 delta_x = 0.01*L ;%mesh makani(mitavand riztar ham shavad.)
 a = 0.5 ;%Defined, but can be changed
  
 x = 0.0 : delta_x : L ;
  X = length(x) ;


error = zeros(5,X) ;
Eu=zeros(5,X);
Eanalytical=zeros(5,X);
 
 
 
%%%Initial condition A

delta_t=[0.00025,0.0002,0.0001,0.00005,0.00001];  %%riz kardan mesh(zamani)
col=['r','k','c','g','m']

for i = 1 : length(delta_t)
    
     %dt and dx 
  x = 0.0 : delta_x : L ;
  X = length(x) ;
  dt = delta_t(i) ;
  total_time = 100 ;
  
  %pre defining matrix
  u = zeros(X,1) ;
  analytical = zeros(X,1) ;
  error = zeros(X,1) ;
    
 
  
  
  
  %nu=alfa*dt/dx
  nu = [a*delta_t(1)/delta_x, a*delta_t(2)/delta_x, a*delta_t(3)/delta_x, a*delta_t(4)/delta_x, a*delta_t(5)/delta_x] ; 
  
  for ix = 1 : X
      if (x(ix)<0.25)
          u(ix) = 1.0 ;
      else 
          u(ix) = 0.0 ;
      end
  end
  
  R = 0.0*u ;
  for it = 1 : total_time
      for ix = 2 : X-1
          u_x = (u(ix)-u(ix-1))/(delta_x) ;  %discritization: u is "positive"
          R(ix) = -a*u_x ;%a=0.5
      end
      u = u + dt*R ;
      u
      
      Eu(i,ix)= u(ix);
     
      
      %Eu(i,ix)= u(ix);
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
    axis([0 1 -0.2 1.2])
    
    
    xlabel 0<x<1
    ylabel u
    title 'Upwind:t=100'
   for j=1:length(delta_t)
        for ix = 1:X
  error(j,ix) = abs(Eu(j,ix)-Eanalytical(j,ix));
        end
    end
end
plot(x,analytical,'linewidth',1.5)
legend('\nu = 0.0125','\nu = 0.01','\nu = 0.005','\nu = 0.0025','\nu = 0.0005', 'Analytical')

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

            amplitude(ix) =1- noo1(p)+noo1(p)*(cos(beta(ix)));
            amplitudey(ix)= noo1(p)*sin(beta(ix)) ;
%             
      end
    
    figure(2)
    hold on
    grid on
    plot (amplitude,amplitudey,col(p),'linewidth',1.5)
    legend('\nu = 0.5','\nu = 0.25','\nu = 0.05','\nu = 0.025','\nu = 0.005')
   
    xlabel 'real part'
    ylabel 'imaginary part'
    title 'G=(1-\nu+\nu (cos(\beta)+i sin(\beta))'
end
