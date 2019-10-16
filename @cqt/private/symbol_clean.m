function [am, ap] = symbol_clean(am, ap, nrm)

epsilon = nrm * cqtoption('threshold');
swapped = false;

while ~isempty(am) && minimal_cut(am, ap) < epsilon      
  if length(am) == 1
        if abs(am) > 0 && abs(am) < abs(ap(end))
            epsilon = epsilon - abs(am);
            am = 0.0;
            ap(1) = 0.0;
        else
            bs = min(256, max(16, length(ap) / 16));
            if length(ap) > bs
                % Do it in blocks, that will likely be much more efficient
                
                ee = abs(ap(end:-1:end-bs+1));
                e1 = cumsum(ee) < epsilon;
                
                idx = find(e1 == 0);
                
                if isempty(idx)
                    idx = bs;
                else
                    idx = idx - 1;
                end
                
                epsilon = epsilon - sum(ee(1:idx));
                ap = ap(1:end-idx);
            else
                if abs(ap(end)) < epsilon
                    epsilon = epsilon - abs(ap(end));
                    ap = ap(1:end-1);

                    if isempty(ap)
                        am = [];
                        ap = [];
                    end
                end
            end
        end
  else  
      if length(ap) == 1
          ap = am; am = ap(1); swapped = ~swapped;
      else
        if abs(am(end)) < abs(ap(end))
          epsilon = epsilon - abs(am(end));
          am = am(1:end-1);
        else
          epsilon = epsilon - abs(ap(end));
          ap = ap(1:end-1);
        end
      end
  end
end

if swapped
    tmp = ap; ap = am; am = tmp;
end

end

function r = minimal_cut(am, ap)

    r = 0.0;

    % Distinguish a few cases
    nam = length(am); nap = length(ap);
    
    if nam == 1 && nap == 1
        r = abs(am(1));
    elseif (nam == 1 && am(1) == 0) || (nap == 1 && ap(1) == 0)
        r = max(abs(am(end)), abs(ap(end)));
    elseif nam > 1 || nap > 1
        r = min(abs(am(end)), abs(ap(end)));
    end
    
end
