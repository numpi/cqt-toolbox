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
            if abs(ap(end)) < epsilon
                epsilon = epsilon - abs(ap(end));
                ap = ap(1:end-1);
                
                if isempty(ap)
                    am = [];
                    ap = [];
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
    
    if length(am) == 1
        r = abs(am(1));
    else
        r = abs(am(end));
    end
    
    if length(ap) == 1 && length(am) > 1
        r = r + abs(ap(1));
    else
        r = r + abs(ap(end));
    end
    
    
end
