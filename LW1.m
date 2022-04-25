clc
clear all
close all

 L = 1;%[0,1]
 delta_x = 0.001*L ;%mesh makani(mitavand riztar ham shavad.)

 
 a = 0.5 ;%Defined, but can be changed
  
 x = 0.0 : delta_x : L ;
  X = length(x) ;

delta_t=[0.001,0.0005,0.0001,0.00005,0.00001];  %%riz kardan mesh(zamani)
col=['r','k','c','g','m']
error = zeros(X,1) ;
for i = 1 : length(delta_t)
    
     %dt and dx 
  x = 0.0 : delta_x : L ;
  X = length(x) ;
  dt = delta_t(i) ;
  total_time = 100 ;
  
  %pre defining matrix
  u = zeros(X,1) ;
  analytical = zeros(X,1) ;
  error = zeros(5,X) ;
    
 
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
%       for ix = 2 : X-1
%         u(ix) = 0.5*((u(ix+1) + u(ix-1)) - (a*dt/delta_x)*(u(ix+1) - u(ix-1))) ;  %discritization
%          
%       end
u(1) = 1 ;
      for ix = 2 : X-1
        u(ix) = u(ix)-(a*0.5*dt/delta_x)*((u(ix+1)-u(ix-1))) + 0.5*a^2.0*dt^2.0*((u(ix+1)-2.0*u(ix)+u(ix-1)))/(delta_x^2) ;
      end  
     % u0 = u ;
    
      
        for ix = 1 : X
      if ( x(ix)<0.25)
          analytical(ix) = 1.0 ;
      else 
          analytical(ix) = 0.0 ;
      end
  end
      
  end
  
    hold on
    plot (x,u,col(i),'linewidth',1.2)
    axis([0 1 0 1])
    
    
    xlabel 0<x<1
    ylabel u
    title 'Lax-Wendroff'
  error(i,ix) = abs(u(ix)-analytical(ix))
end
plot(x,analytical,'linewidth',1.5)
xlim([0 1])
ylim([-0.4 1.4])
legend('\nu = 0.5','\nu = 0.25','\nu = 0.05','\nu = 0.025','\nu = 0.005', 'analytical')
