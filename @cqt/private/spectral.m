function [l, u, w, n] = spectral(vm, vp)
%SPECTRAL 

switch cqtoption('inversion')
	case 'cr'
		[l, u] = spectral_cr(vm, vp);		
	case 'fft'
		[l, u] = spectral_fft(vm, vp);
end


end

