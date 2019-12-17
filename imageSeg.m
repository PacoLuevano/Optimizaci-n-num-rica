function [z_med,iter] = imageSeg(w,K,lambda,M,N1,N2,mu)
  rho = 1.9;
  mu=1;
  gamma = 0.01;
  tau = 0.01;
  
  tol = 1e-5;
  maxIter = 20;
  
  z = zeros(M,N1,N2);
  q = zeros(M,1);
  u = zeros(M,N1,N2);
  v = zeros(M,N1,N2,2);
  
  z1 = ones(M,N1,N2);
  
  z_med = zeros(M,N1,N2);
  q_med = zeros(M,1);
  u_med = zeros(M,N1,N2);
  v_med = zeros(M,N1,N2,2);
  
  sigma_u = gamma/tau/(1+1/mu);
  sigma_v = (1-gamma)/tau/8;
  
  proj_simplex_array = @(y) max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);
  proj_l1ball_vector = @(y) max(abs(y)-max(max((cumsum(sort(abs(y),1,'descend'),1)-(K*sqrt(mu*N1*N2)))./(1:size(y,1))'),0),0).*sign(y);
  opD = @(x) cat(4,cat(2,diff(x,1,2),zeros(size(x,1),1,size(x,3))),cat(3,diff(x,1,3),zeros(size(x,1),size(x,2))));
  %opDadj = @(v) - cat(3,cat(1,v(1,:,:,1),diff(v,1,1)(:,:,:,1))) - cat(2,v(:,1,:,2),diff(v,1,2)(:,:,:,2));
  opS = @(u) sum(sum(u,3),2);
  opSadj = @(q) cat(3,q*ones(1,N1),q*ones(1,N1));
  l1_inf_norm = @(z) abs(sum(max(max(z,[],3),[],2),1));
  
  iter = 1;
  while (l1_inf_norm(z-z1)>tol && iter<maxIter)
  par = z - tau*(u + w + opDadj(v));
  z_med = proj_simplex_array(par(:,:,1));
  for i = 2:N2
    z_med=cat(3,z_med,proj_simplex_array(par(:,:,i)));
  end
  
  par2 = q + ((tau/sqrt(mu*N1*N2))*opS(u));
  q_med = projC(par2, K*sqrt(mu*N1*N2));
  
  sadj = (2*q_med-q)*ones(1,N1);
  for i = 2:N2
    sadj = cat(3,sadj,(2*q_med-q)*ones(1,N1));
  end
  par3 = u + sigma_u*(2*z_med - z - sadj/sqrt(mu*N1*N2));
  u_med = zeros(M,N1,N2);
  for i = 1:N1
    for j = 1:N2
      u_med(:,i,j)= lsqnonneg(eye(M),par3(:,i,j));
    end
  end
  
  par4 = v + sigma_v*opD(2*z_med - z);
  v_med = zeros(M,N1,N2,2);
  for i = 1:M
    for j = 1:N1
      for k = 1:N2
        v_med(i,j,k,:) = projball_l2([par4(i,j,k,1);par4(i,j,k,2)], 0.5*lambda);
      end
    end
  end
  z1=z;
  z = z + rho*(z_med-z);
  q = q + rho*(q_med-q);
  u = u + rho*(u_med-u);
  v = v + rho*(v_med-v);
  iter = iter+1
  end
end

  
  
  
  
  
  
  
  
  
  
  
  
  
  