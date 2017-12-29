function [am, ap] = symbol_clean(am, ap, nrm)

epsilon = nrm * cqtoption('threshold');

[am, ap, ~] = symbol_clean_rec(am, ap, epsilon);

end

function [am, ap, epsilon] = symbol_clean_rec(am, ap, epsilon)

    if isempty(am)
        return;
    end

    if min(abs([ am(end), ap(end) ])) < epsilon        

        if length(am) == 1
            if abs(am) > 0 && abs(am) < abs(ap(end))
                epsilon = epsilon - abs(am);
                am = 0.0;
                ap(1) = 0.0;
                [am, ap, epsilon] = symbol_clean_rec(am, ap, epsilon);
            else
                if abs(ap(end)) < epsilon
                    epsilon = epsilon - abs(ap(end));
                    ap = ap(1:end-1);
                
                    if isempty(ap)
                        am = [];
                    end
                    
                    [am, ap, epsilon] = symbol_clean_rec(am, ap, epsilon);
                end
            end
            
            
            return;
        end
        
        if length(ap) == 1
            [ap, am, epsilon] = symbol_clean_rec(ap, am, epsilon);
            return;
        end                        
            
        if abs(am(end)) < abs(ap(end))
            epsilon = epsilon - abs(am(end));
            am = am(1:end-1);
        else
            epsilon = epsilon - abs(ap(end));
            ap = ap(1:end-1);
        end
        
        [am, ap, epsilon] = symbol_clean_rec(am, ap, epsilon);
    end
end
