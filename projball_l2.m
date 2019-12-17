%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_2 (Frobenius) norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : l2norm(x)<=tau },
% where tau>=0 and
% l2norm = @(x) sqrt(sum(x(:).^2)),
% of the N-D array y.

function x = projball_l2 (y, tau)
	tmp = sqrt(sum(y(:).^2));
	if tmp<=tau, x = y;
	else 
		x = y * (tau/tmp);
	end
end