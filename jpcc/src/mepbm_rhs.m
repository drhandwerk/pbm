function dn = mepbm_rhs(t,n,mech,kargs,maxsize,r)

switch mech
    case '3step_alt'
        if length(kargs) ~= 7
            error('3step_alt requires 7 parameters.');
        end
        dn = rhs_3step_alt(t,n,kargs(1),kargs(2),kargs(3),kargs(4),kargs(5),kargs(6),kargs(7));
    otherwise
        error('Unrecognized mechanism name.')
end

    function dn = rhs_3step_alt(~,n,S,kf,kb,k1,k2,k3,M)
        M = floor(M);
        dn = zeros(maxsize+1,1);
        dn(1) = -k1.*n(1).*n(2).^2 - k2.*n(1).*sum(r(3:M).*n(3:M)'.*(3:M)) - k3.*n(1).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize)) - kf.*n(1).*S^2 + kb.*n(2).*n(maxsize+1); % precursor A
        dn(2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+1) - 2*k1.*n(1).*n(2).^2; % A_solv
        dn(3) = k1.*n(1).*n(2).^2 - 3.*k2.*n(1).*r(3).*n(3); % nucleus
        
        for j = 4:M
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        
        for j = M+1
            dn(j) = k2.*n(1).*r(j-1).*n(j-1).*(j-1) - k3.*n(1).*r(j).*n(j).*j;
        end
        
        for j = M+2:maxsize
            dn(j) = k3.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        dn(maxsize+1) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+1) + k1.*n(1).*n(2).^2 + k2.*n(1).*sum(r(3:M).*n(3:M)'.*(3:M)) + k3.*n(1).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize)); % POM
    end
end