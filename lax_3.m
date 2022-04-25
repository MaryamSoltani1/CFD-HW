clc
clear all
close all

 L = 1;%[0,1]
 delta_x = 0.001*L ;%mesh makani(mitavand riztar ham shavad.)
 %%%har ghadr zarib bala riztar shod javab behtar shod!!!!
 
 a = 0.5 ;%Defined, but can be changed
  
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
      if ( x(ix)<0.2)
          u(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          u(ix) = 1.0 ;
      else
          u(ix) = 0.0 ;      
      end
  end 
    
  %%%%%%noo=1 javab kooooo?
  for it = 1 : total_time
      for ix = 2 : X-1
        u(ix) = 0.5*((u(ix+1) + u(ix-1)) - (a*dt/delta_x)*(u(ix+1) - u(ix-1))) ;  %discritization
         
      end
      
       for ix = 1 : X
      if ( x(ix)<0.2)
          analytical(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          analytical(ix) = 1.0 ;
      else
          analytical(ix) = 0.0 ;      
      end
  end
  end
  
    hold on
    plot (x,u,col(i),'linewidth',1.2)
    axis([0 1 -0.2 1.2])
    
    
    xlabel 0<x<1
    ylabel u
    title 'Lax'
  
end
plot(x,analytical,'linewidth',1.5)
legend('\nu = 0.125','\nu = 0.1','\nu = 0.05','\nu = 0.025','\nu = 0.005', 'analytical')
noo1 = [1.25,1,0.5,0.25,0.05] ;

 for p = 1 : length(noo1)
  for ix = 1 : X
            k(ix) = 2*pi*ix/(2*L) ;
            beta(ix) = k(ix)*delta_x ; %k_m*h
            amplitude(ix) = (cos(beta(ix))*cos(beta(ix)) + noo1(p)^2*sin(beta(ix))*sin(beta(ix)))^0.5 ;
             phi(ix) = atan(-noo1(p)*tan(beta(ix))) ;
      end
    
    figure(2)
    hold on
    grid on
   
    plot (phi/pi,-amplitude,col(p),'linewidth',1.2)
    legend('\nu = 1.25','\nu = 1','\nu = 0.5','\nu = 0.25','\nu = 0.05')
   
    xlabel 'relative phase:\phi/\pi'
    ylabel amplitude
 end




%    for ix = 1 : X
%       error(ix) = abs(u(ix)-analytical(ix)) ;
%      fprintf('\n x: %f  Numerical u: %f  Analytical u: %f Error: %f \n', x(ix), u(ix), analytical(ix), error(ix))
%    end