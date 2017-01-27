function l = isconsecutive(v)
%ISCONSECUTIVE Check if the indices in V are consecutive.
%
% L = ISCONSECUTIVE(V) checks if all the indices in V are consecutive. 

if isempty(v)
	l = true;
else
	l = all((v(2:end) - v(1:end-1)) == 1);
end


end

