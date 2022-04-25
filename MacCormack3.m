clc
clear all
close all

 L = 1;%[0,1]
 delta_x = 0.001*L ;%mesh makani(mitavand riztar ham shavad.)
 %%%har ghadr zarib bala riztar shod javab behtar shod!!!!
 
 a = 0.5 ;%Defined, but can be changed
  


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
      if (x(ix)>=0.0 && x(ix)<0.2)
          u(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          u(ix) = 1.0 ;
      else
          u(ix) = 0.0 ;      
      end
  end 
  
  for it = 1 : total_time
      %game aval
      for ix = 2 : X-1
         us(ix) = u(ix) - a*dt*(u(ix+1)-u(ix))/delta_x ; 
      end
      %game dovom
     for ix = 2 : X-1
         u(ix) = 0.5*((u(ix)+us(ix)) - a*dt*(us(ix)-us(ix-1))/delta_x) ; 
     end

 
        for ix = 1 : X
      if (x(ix)>=0.0 && x(ix)<0.2)
          analytical(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          analytical(ix) = 1.0 ;
      else
          analytical(ix) = 0.0 ;      
      end
  end 
      
  end
  
    hold on
    scatter (x,u,col(i),'linewidth',1.2)
    axis([0 1 0 1])
    
    
    xlabel 0<x<1
    ylabel u
    title 'MacCormack'
  
end
plot(x,analytical,'linewidth',1.5)
xlim([0 1])
ylim([-0.4 1.4])
legend('\nu = 0.125','\nu = 0.1','\nu = 0.05','\nu = 0.025','\nu = 0.005', 'analytical')
%    for ix = 1 : X
%       error(ix) = abs(u(ix)-analytical(ix)) ;
%      fprintf('\n x: %f  Numerical u: %f  Analytical u: %f Error: %f \n', x(ix), u(ix), analytical(ix), error(ix))
%    end