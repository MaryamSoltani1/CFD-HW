delta_t=[0.00025,0.0002,0.0001,0.00005,0.00001] ;
col=['r','k','c','g','m'];
 L = 1 ;%x is in: [0,1]
  dx = 0.01*L ; %game makani: 0.01
  a = 0.5 ;

for i = 1 : length(delta_t)
    
  dt = delta_t(i) ;
  nt = 100 ;%%%% ino age ziad konam kheili ajib mishe!!!
    
  
  x = 0.0 : dx : L ;
  nx = length(x) ;
  
  u = zeros(nx,1) ;
  analytical = zeros(nx,1) ;
  error = zeros(nx,1) ;
  
  nu = [a*delta_t(1)/delta_x, a*delta_t(2)/delta_x, a*delta_t(3)/delta_x, a*delta_t(4)/delta_x, a*delta_t(5)/delta_x] ;  %c*delta_t/delta_x
  
  
  for ix = 1 : nx
      if (x(ix)>=0.0 && x(ix)<0.2)
          u(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          u(ix) = 1.0 ;
      else
          u(ix) = 0.0 ;      
      end
  end
  
  R = 0.0*u ;
  for it = 1 : nt
      for ix = 2 : nx-1
          u_x = (u(ix)-u(ix-1))/(dx) ; % discritization: u is positive. 
          R(ix) = -a*u_x ;
      end
      u = u + dt*R ; 
      
  for ix = 1 : nx
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
    xlabel 0<x<1
    ylabel u
    title 'Upwind:t=100'
  
end
plot(x,analytical,'b','linewidth',1.5)
legend('\nu = 0.0125','\nu = 0.01','\nu = 0.005','\nu = 0.0025','\nu = 0.0005', 'analytical') 

   for ix = 1 : nx
       error(ix) = abs(u(ix)-analytical(ix)) ;
     fprintf('\n x: %f  Numerical u: %f  Analytical u: %f Error: %f \n', x(ix), u(ix), analytical(ix), error(ix))
   end
