function v = opDadj(x)
  a = diff(x,1,1);
  b = diff(x,1,2);
  v = - cat(3,cat(1,x(1,:,:,1),a(:,:,:,1))) - cat(2,x(:,1,:,2),b(:,:,:,2));