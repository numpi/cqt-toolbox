function TestCqtQuadratic

% Quadratic equations

for j = [1 2 4 7 8]
    [Am1, A0, A1] = cqtgallery('jackson', j);    
    l = A0(10, 10);
    
    Am1 = - Am1 / l;
    A0  = -A0 / l;
    A1  = -A1 / l;
    
    G = quadnewt(Am1, A0, A1);
    
    CheckTestResult(norm(Am1 + A0 * G + A1 * G^2), '<', ...
        1e-6, ...
        sprintf('CQT Newton-based quadratic equation solver (Jackson problem %d)', j));    
    
    G = cr(Am1, A0, A1);
    
    CheckTestResult(norm(Am1 + A0 * G + A1 * G^2), '<', ...
        1e-6, ...
        sprintf('CQT CR quadratic equation solver (Jackson problem %d)', j));
end

end

