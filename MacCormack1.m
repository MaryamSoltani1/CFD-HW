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
      if (x(ix)<0.25)
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
      if ( x(ix)<0.25)
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
    title 'MacCormack:t=100'
  
end
plot(x,analytical,'linewidth',1.5)
xlim([0.05 1])
ylim([-0.4 1.4])
legend('\nu = 0.125','\nu = 0.1','\nu = 0.05','\nu = 0.025','\nu = 0.005', 'analytical')

noo1 = [0.25,0.1,0.05,0.025,0.005] ;
for p = 1 : length(noo1)
  for ix = 1 : X
            k(ix) = 2.0*pi*ix/(2.0*L) ;
            beta(ix) = k(ix)*delta_x ;
            amplitude(ix) = noo1(p)^2*cos(beta(ix))/2+(1-noo1(p)^2) ;
            amplitudey(ix)=-i*sin(beta(ix))
%             phi(ix) = atan(-CFL*tan(beta(ix))) ;
      end
    
    figure(2)
    hold on
    grid on
    plot (amplitude,amplitudey,'linewidth',1.2)
    xlabel 'relative phase'
    ylabel amplitude
end