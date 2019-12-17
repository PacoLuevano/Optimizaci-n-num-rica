function img = imagen(z,a)
  for j = 1:size(z,2)
    for k = 1: size(z,3)
      img(j,k,:) = zeros(3,1);
      for i = 1:size(a,1)
        mult = z(i,j,k)*a(i,:)';
        img(j,k,1) = img(j,k,1) + mult(1);
        img(j,k,2) = img(j,k,2) + mult(2);
        img(j,k,3) = img(j,k,3) + mult(3);
      end
    end
  end
end  