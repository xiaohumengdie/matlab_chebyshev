%PCSI Algorithm
function [x_pcsi,m] = Pcsi_iter(A,b,x0,M,tol,MaxEigs,MinEigs,iterNum)
  
  csalpha = 2.0/(MaxEigs-MinEigs);
  csbeta = (MaxEigs+MinEigs)/(MaxEigs-MinEigs);
  csy = csbeta/csalpha;csomga = 2.0/csy;
  x_pcsi=x0;
  r=b-A*x_pcsi;
  r=M*r;
  q = (1.0/csy)*r;
  x_pcsi=x_pcsi+q;
  r=b-A*x_pcsi;
  for m = 1:iterNum
    csomga = 1.0/(csy-csomga/(4.0*csalpha*csalpha));
    %disp(csomga);
    r = M*r;
    q = csomga*r+(csy*csomga-1.0)*q;
    x_pcsi = x_pcsi + q;
    r=b-A*x_pcsi;
    if(norm(r,2)<tol)
      break
    end
  end
end
