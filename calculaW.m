function w = calculaW(pal,y)
  
  N1 = size(y,1);  %height
  N2 = size(y,2);  %width
  
  
  M = size(pal,1);
  
  w = zeros(M,N1,N2);
  
  for i = 1:M
    i
    for j = 1:N1
      %j
      for k = 1:N2
        %k
        aux = double(squeeze(y(j,k,:))) - pal(i,:)';
        w(i,j,k) = 0.5*dot(aux,aux);
      end
    end
  end
end

       