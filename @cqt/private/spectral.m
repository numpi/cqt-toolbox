function [l, u, w, n] = spectral(vm, vp)
%SPECTRAL 

switch cqtoption('inversion')
	case 'cr'
		if nargout > 2
			[l, u, w, n] = spectral_cr(vm, vp);
		else
			[l, u] = spectral_cr(vm, vp);
		end
	case 'fft'
		[l, u] = spectral_fft(vm, vp);
end


end

