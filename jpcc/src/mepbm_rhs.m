function dn = mepbm_rhs(t,n,mech,kargs,maxsize,r)

%% Call appropriate rhs with 'kargs' given 'mech'
switch mech
    % 2-step
    %     A -k1> B
    % A + B -k2> 2B
    case '2step'
        if length(kargs) ~= 3
            error('2step requires 3 parameters.');
        end
        dn = rhs_2step(t,n,kargs(1),kargs(2),kargs(3));
        
    case '2step_alt'
        if length(kargs) ~= 5
            error('2step_alt requires 5 parameters.');
        end
        dn = rhs_2step_alt(t,n,kargs(1),kargs(2),kargs(3),kargs(4),kargs(5));
        
    % 3-step
    %     A -k1> B
    % A + B -k2> 2B
    % A + C -k3> 1.5C
    case '3step'
        if length(kargs) ~= 4
            error('3step requires 4 parameters.');
        end
        dn = rhs_3step(t,n,kargs(1),kargs(2),kargs(3),kargs(4));
        
    case '3step_alt'
        if length(kargs) ~= 7
            error('3step_alt requires 7 parameters.');
        end
        dn = rhs_3step_alt(t,n,kargs(1),kargs(2),kargs(3),kargs(4),kargs(5),kargs(6),kargs(7));
    
    % 4-step
    %     A -k1> B
    % A + B -k2> 2B
    % B + B -k3> C
    % A + C -k4> 1.5C
    case '4step'
        if length(kargs) ~= 5
            error('4step requires 5 parameters.');
        end
        dn = rhs_4step(t,n,kargs(1),kargs(2),kargs(3),kargs(4),kargs(5));
        
    case '4step_alt'
        if length(kargs) ~= 8
            error('4step_alt requires 8 parameters.');
        end
        dn = rhs_4step_alt(t,n,kargs(1),kargs(2),kargs(3),kargs(4),kargs(5),kargs(6),kargs(7),kargs(8));
    
    % 5-step
    %     A -k1> B
    % A + B -k2> 2B
    % B + B -k3> C
    % A + C -k4> 1.5C
    % B + C -k5> C
    case '5step'
        if length(kargs) ~= 6
            error('5_step requires 6 parameters.');
        end
        dn = rhs_4step_alt(t,n,kargs(1),kargs(2),kargs(3),kargs(4),kargs(5),kargs(6));
        
    otherwise
        error('Unrecognized mechanism name.')
end


%% 2-step with 3rd-order nucleation
    function dn = rhs_2step(~,n,k1,k2,k3)
        dn = zeros(maxsize+1,1);
        dn(1) = -k3.*k1.*n(1).^k3 - k2.*n(1).*sum(r(3:maxsize).*n(3:maxsize)'.*(3:maxsize)); % precursor
        dn(2) = 0; % No particles of size 2
        dn(3) = k1.*n(1).^k3 - k3.*k2.*n(1).*r(3).*n(3); % nucleus
        
        for j = 4:maxsize
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        dn(maxsize+1) = k2.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses "off-screen".
    end

%% 2-step alternative nucleation
    function dn = rhs_2step_alt(~,n,S,kf,kb,k1,k2)
        dn = zeros(maxsize+2,1);
        dn(1) = -k1.*n(1).*n(2).^2 - k2.*n(1).*sum(r(3:maxsize).*n(3:maxsize)'.*(3:maxsize)) - kf.*n(1).*S^2 + kb.*n(2).*n(maxsize+2); % precursor A
        dn(2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+2) - 2*k1.*n(1).*n(2).^2; % A_solv
        dn(3) = k1.*n(1).*n(2).^2 - 3.*k2.*n(1).*r(3).*n(3); % nucleus
        
        for j = 4:maxsize
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        dn(maxsize+1) = k2.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'.
        dn(maxsize+2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+2) + k1.*n(1).*n(2).^2 + k2.*n(1).*sum(r(3:maxsize).*n(3:maxsize)'.*(3:maxsize)); % POM
    end

%% 3-step standard 3rd-order nucleation
    function dn = rhs_3step(~,n,k1,k2,k3,N)
        N = floor(N);
        dn = zeros(maxsize+1,1);
        dn(1) = -3.*k1.*n(1).^3 - k2.*n(1).*sum(r(3:N).*n(3:N)'.*(3:N)) - k3.*n(1).*sum(r(N+1:maxsize).*n(N+1:maxsize)'.*(N+1:maxsize));
        dn(2) = 0;
        dn(3) = k1.*n(1).^3 - 3.*k2.*n(1).*r(3).*n(3);
        
        for j = 4:N
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        
        for j = N+1
            dn(j) = k2.*n(1).*r(j-1).*n(j-1).*(j-1) - k3.*n(1).*r(j).*n(j).*j;
        end
        
        for j = N+2:maxsize
            dn(j) = k3.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        dn(maxsize+1) = k3.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'.
    end

%% 3-step alternative termolecular nucleation
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

%% 4-step standard nucleation
    function dn = rhs_4step(~,n,k1,k2,k3,k4,maxagglomsize)
        dn = zeros(maxsize+1,1);
        dn(1) = -3.0.*k1.*n(1).^3 - k2.*n(1).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize)) - k4.*n(1).*sum(r(maxagglomsize+1:maxsize).*n(maxagglomsize+1:maxsize)'.*(maxagglomsize+1:maxsize));
        dn(2) = 0;
        dn(3) = k1.*n(1).^3 - k2.*n(1).*r(3).*3.0.*n(3) - k3.*r(3).*3.*n(3).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize));
        
        % Gain from agglomeration only affects sizes of > 6.
        for j = 4:5
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k3.*r(j).*j.*n(j).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize));
        end
        
        % "small" particles
        for j = 6:maxagglomsize
            h = floor(j/2);
            offset = mod(j,2);
            agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
            agglom_r = h+(0:h-3)+offset;
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k3.*r(j).*j.*n(j).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize)) + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
        end
        % "large" particles that agglom
        for j = maxagglomsize+1 % need logic that handles size that is formed from k2 but loses at k4.
            h = floor(j/2);
            offset = mod(j,2);
            agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
            agglom_r = h+(0:h-3)+offset;
            bb_index = agglom_r <= maxagglomsize;
            agglom_l = agglom_l(bb_index);
            agglom_r = agglom_r(bb_index);
            dn(j) = k2.*n(1).*r(j-1).*n(j-1).*(j-1) - k4.*n(1).*r(j).*n(j).*j + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
        end
        for j = maxagglomsize+2:2*maxagglomsize
            h = floor(j/2);
            offset = mod(j,2);
            agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
            agglom_r = h+(0:h-3)+offset;
            bb_index = agglom_r <= maxagglomsize;
            agglom_l = agglom_l(bb_index);
            agglom_r = agglom_r(bb_index);
            dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
        end
        % "large" particles that don't agglom
        for j = 2*maxagglomsize+1:maxsize
            dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        % off-screen growth
        dn(maxsize+1) = k4.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'.
    end

%% 4-step alternative nucleation
    function dn = rhs_4step_alt(~,n,S,kf,kb,k1,k2,k3,k4,maxagglomsize)
        dn = zeros(maxsize+2,1);
        dn(1) = -3.0.*k1.*n(1).*n(2).^2 - k2.*n(1).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize))...
            - k4.*n(1).*sum(r(maxagglomsize+1:maxsize).*n(maxagglomsize+1:maxsize)'.*(maxagglomsize+1:maxsize)) - kf.*n(1).*S^2 + kb.*n(2).*n(maxsize+2);
        dn(2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+2) - 2*k1.*n(1).*n(2).^2; % A_solv;
        dn(3) = k1.*n(1).*n(2).^2 - 3.*k2.*n(1).*r(3).*n(3) - k3.*r(3).*3.*n(3).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize));
        
        % Gain from agglomeration only affects sizes of > 6.
        for j = 4:5
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k3.*r(j).*j.*n(j).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize));
        end
        
        % "small" particles
        for j = 6:maxagglomsize
            h = floor(j/2);
            offset = mod(j,2);
            agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
            agglom_r = h+(0:h-3)+offset;
            dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k3.*r(j).*j.*n(j).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize)) + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
        end
        % "large" particles that agglom
        for j = maxagglomsize+1 % need logic that handles size that is formed from k2 but loses at k4.
            h = floor(j/2);
            offset = mod(j,2);
            agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
            agglom_r = h+(0:h-3)+offset;
            bb_index = agglom_r <= maxagglomsize;
            agglom_l = agglom_l(bb_index);
            agglom_r = agglom_r(bb_index);
            dn(j) = k2.*n(1).*r(j-1).*n(j-1).*(j-1) - k4.*n(1).*r(j).*n(j).*j + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
        end
        for j = maxagglomsize+2:2*maxagglomsize
            h = floor(j/2);
            offset = mod(j,2);
            agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
            agglom_r = h+(0:h-3)+offset;
            bb_index = agglom_r <= maxagglomsize;
            agglom_l = agglom_l(bb_index);
            agglom_r = agglom_r(bb_index);
            dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
        end
        % "large" particles that don't agglom
        for j = 2*maxagglomsize+1:maxsize
            dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
        end
        % off-screen growth
        dn(maxsize+1) = k4.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'.
        dn(maxsize+2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+2) + k1.*n(1).*n(2).^2 + k2.*n(1).*sum(r(3:maxagglomsize).*n(3:maxagglomsize)'.*(3:maxagglomsize))...
            + k4.*n(1).*sum(r(maxagglomsize+1:maxsize).*n(maxagglomsize+1:maxsize)'.*(maxagglomsize+1:maxsize)); % POM
    end

%% 5-step
%     A -k1> B
% A + B -k2> 2B
% B + B -k3> C
% A + C -k4> 1.5C
% B + C -k5> C

function dn = rhs_5step(~,n,k1,k2,k3,k4,k5,M)
dn = zeros(maxsize+1,1);
dn(1) = -3.0.*k1.*n(1).^3 - k2.*n(1).*sum(r(3:M).*n(3:M)'.*(3:M)) - k4.*n(1).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize));
dn(2) = 0;
dn(3) = k1.*n(1).^3 - k2.*n(1).*r(3).*3.0.*n(3) - k3.*r(3).*3.*n(3).*sum(r(3:M).*n(3:M)'.*(3:M)) - k5.*r(3).*3.*n(3).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize));

% Gain from agglomeration only affects sizes of > 6.
for j = 4:5
    dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k3.*r(j).*j.*n(j).*sum(r(3:M).*n(3:M)'.*(3:M)) - k5.*r(j).*j.*n(j).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize));
end

% B+B only
for j = 6:M
    h = floor(j/2);
    offset = mod(j,2);
    agglom_l = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
    agglom_r = h+(0:h-3)+offset;
    dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k3.*r(j).*j.*n(j).*sum(r(3:M).*n(3:M)'.*(3:M))...
        - k5.*r(j).*j.*n(j).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize)) + k3.*sum(r(agglom_l).*agglom_l.*n(agglom_l)'.*r(agglom_r).*agglom_r.*n(agglom_r)');
end
% B+B and B+C range
for j = M+1 % need logic that handles size that is formed at rate k2 but loses at rate k4.
    h = floor(j/2);
    offset = mod(j,2);
    agglom_l_bb = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
    agglom_r_bb = h+(0:h-3)+offset;
    bb_index = agglom_r_bb <= M;
    agglom_l_bb = agglom_l_bb(bb_index);
    agglom_r_bb = agglom_r_bb(bb_index);
    dn(j) = k2.*n(1).*r(j-1).*n(j-1).*(j-1) - k4.*n(1).*r(j).*n(j).*j + k3.*sum(r(agglom_l_bb).*agglom_l_bb.*n(agglom_l_bb)'.*r(agglom_r_bb).*agglom_r_bb.*n(agglom_r_bb)')...
        - k5.*r(j).*j.*n(j).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize));
end

% the first that can form from B+C is M+4 since 3 is the smallest B
% and M+1 is the smallest C. Therefore, These two and the previous one do not have
% gain from B+C
for j = M+2:M+3
    h = floor(j/2);
    offset = mod(j,2);
    agglom_l_bb = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
    agglom_r_bb = h+(0:h-3)+offset;
    bb_index = agglom_r_bb <= M;
    agglom_l_bb = agglom_l_bb(bb_index);
    agglom_r_bb = agglom_r_bb(bb_index);
    agglom_l_bc = 3:(j-(M+1));
    agglom_r_bc = (j-3):-1:M+1;
    dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) + k3.*sum(r(agglom_l_bb).*agglom_l_bb.*n(agglom_l_bb)'.*r(agglom_r_bb).*agglom_r_bb.*n(agglom_r_bb)')...
        - k5.*r(j).*j.*n(j).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize));
end

for j = M+4:2*M
    h = floor(j/2);
    offset = mod(j,2);
    agglom_l_bb = h-(0:h-3); % agglom_l + agglom_r = [j,...,j]
    agglom_r_bb = h+(0:h-3)+offset;
    bb_index = agglom_r_bb <= M;
    agglom_l_bb = agglom_l_bb(bb_index);
    agglom_r_bb = agglom_r_bb(bb_index);
    agglom_l_bc = 3:(j-(M+1));
    agglom_r_bc = (j-3):-1:M+1;
    dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) + k3.*sum(r(agglom_l_bb).*agglom_l_bb.*n(agglom_l_bb)'.*r(agglom_r_bb).*agglom_r_bb.*n(agglom_r_bb)')...
        - k5.*r(j).*j.*n(j).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize)) + k5.*sum(r(agglom_l_bc).*agglom_l_bc.*n(agglom_l_bc)'.*r(agglom_r_bc).*agglom_r_bc.*n(agglom_r_bc)');
end
% B+C only
for j = 2*M+1:maxsize
    agglom_l_bc = 3:(j-(M+1));
    agglom_r_bc = (j-3):-1:M+1;
    dn(j) = k4.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j) - k5.*r(j).*j.*n(j).*sum(r(M+1:maxsize).*n(M+1:maxsize)'.*(M+1:maxsize))...
        + k5.*sum(r(agglom_l_bc).*agglom_l_bc.*n(agglom_l_bc)'.*r(agglom_r_bc).*agglom_r_bc.*n(agglom_r_bc)');
end

% off-screen growth
dn(maxsize+1) = k4.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'. TODO add off screen B+C terms
end

