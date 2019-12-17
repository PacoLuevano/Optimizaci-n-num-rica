function proj = projC(x,K)
  if sum(x)<=K
    proj=x;
  else
    A = [K; zeros(size(x,1)-1,1)];
    B = [0; K ; zeros(size(x,1)-2,1)];
    proj = A + (dot(x-A,B-A)/dot(B-A,B-A))*(B-A);
  end
end