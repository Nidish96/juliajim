function [U, dUdlam, Ss, flag, Scall] = CONTINUE(func, u0, lam0, lam1, ds, varargin)
%CONTINUE Conducts the continuation and solves the system
%
% USAGE:
%   [U, dUdlam] = CONTINUE(func, u0, lam0, lam1, ds, varargin);
% INPUTS:
%   func	: Function handle (Jacobian is strictly necessary)
%   u0		: (Nu,1) initial guess for initial point
%   lam0	: Parameter value at initial point
%   lam1	: Parameter value at final point
%   ds		: Arc length parameter at initial point
%   Copt 	: (optional) Options structure
% OUTPUTS:
%   U		: (Nu+1, Np) Solution vector
%   dUdlam	: (Nu+1, Np) Solution derivative wrt parameter
%
%
% Solver List: 
%   orthogonal
%   arclength             : Measures arclength based on norm(u./Dscale)
%                           If one has many unknowns, the relative
%                           weighting of 'lam' can become very small. 
%                           In these cases, use 'K0NormalizedArcLength'
%                           instead.
%   K0NormalizedArcLength : Allows for a general stiffness matrix and 
%                           measures vector length as u(1:end-1)'*K0*u(1:end-1).
%                           Default options use Dscale to construct a
%                           diagonal K0.
%                           A diagonal K0 differs from 'arclength' only in
%                           that the relative weights of the unknowns can
%                           be more appropriately controlled. 
%                           Allows variable weighting between lam and norm
%                           of u(1:end-1). However, giving small or 0 
%                           weight to lam may result in NaN. Also, if the
%                           curve has dL/ds = 0 or approx 0, the predictor
%                           will be poorly conditioned.

%% Default options
    Copt = struct('Nmax', 100, 'dsmax', ds*5, 'dsmin', ds/5, ...
                  'angopt', pi/6, 'startdir', sign(lam1-lam0),...
                  'Display', true, 'itDisplay', false, 'lsrch', 0, ...
                  'nev', 1, 'adj', 1, 'itopt', 6, ...
                  'adapt', 1, 'arclengthparm', 'orthogonal', 'ITMAX', 20, ...
                  'parmjac', true, ...
                  'Dscale', ones(length(u0)+1,1), 'DynDscale', 0, ...
                  'opts', struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', ...
                                 1e-6, 'utol', 1e-6, ...
                                 'Display', false, 'lsrch', 0), ...
                  'solverchoice', 1, ...
                  'callbackfun', [], ...
                  'arcsettings', struct());
    if nargin==6
        nflds = fieldnames(varargin{1});
        for i=1:length(nflds)
            if ~strcmp(nflds{i}, 'opts')
                Copt.(nflds{i}) = varargin{1}.(nflds{i});
            else
                nnflds = fieldnames(varargin{1}.(nflds{i}));
                for j=1:length(nnflds)
                    Copt.opts.(nnflds{j}) = varargin{1}.(nflds{i}).(nnflds{j});
                end
            end
        end
    end
    Copt.opts.lsrch = Copt.lsrch;
    Copt.opts.Display = Copt.itDisplay;
    Copt.opts.ITMAX = Copt.ITMAX;  
    
    %% Allocations
    U = zeros(length(u0)+1, Copt.Nmax);
    dUdlam = zeros(length(u0)+1, Copt.Nmax);

    J0 = zeros(length(u0)+1);  
    
    ds0 = ds;
    
    %% Callback Function
    callback = ~isempty(Copt.callbackfun);
    if callback
        Scall = cell(Copt.Nmax, 1);
    else
        Scall = cell(0,1);
    end
    
    %% Correct initial solution
    Copt.opts.Dscale = Copt.Dscale(1:end-1);
    tm = Copt.opts.Display;
    Copt.opts.Display = true;
    switch Copt.solverchoice
      case 3
        oopts = optimset('Jacobian', 'on', 'Display', 'iter')
        [U(1:end-1,1), ~, eflag] = fsolve(@(u) func([u; lam0]), u0, oopts);
      otherwise
        [U(1:end-1,1), ~, eflag] = NSOLVE(@(u) func([u; lam0]), u0, Copt.opts);
    end
    Copt.opts.Display = tm;
    U(end, 1) = lam0;
    if eflag <= 0
        U = zeros(size(u0,1)+1, 0);
        dUdlam = U;
        Ss = [];
        flag = -1;
        fprintf('Initial point non-convergent!');
        return;
    elseif Copt.Display
        disp('Initial Point Converged')
        fprintf('Starting Continuation from %f to %f\n', lam0, lam1);
    end
    
    if Copt.DynDscale
        % mnis = find(abs(U(:,1))>eps & abs(U(:,1))<100*Copt.Dscale);
        Copt.opts.Dscale = max(abs(U(:, 1)), Copt.Dscale);
        % Copt.opts.Dscale(mnis) = abs(U(mnis,1));
    else
        Copt.opts.Dscale = Copt.Dscale;
    end
    
    %% If Parametric jacobian is not supplied
    if ~Copt.parmjac
        ofunc = func;
        %     func = @(ul) REVRES(ul, ofunc, Copt.dsmin/1000);
        func = @(ul) REVRES(ul, ofunc, (lam1-lam0)/1e4);
    end
    
    %% Extract gradient information 
    [~, J0(1:end-1,1:end-1), J0(1:end-1,end)] = func(U(:, 1)); % (This could be eliminated by saving the Jacobian from NSOLVE)
    J0 = J0.*Copt.opts.Dscale';
    dUdlam(:, 1) = [-J0(1:end-1,1:end-1)\J0(1:end-1,end); 1];  % Scaled Tangent
    
    %% Evaluate Callback if necessary
    if callback
        Scall{1} = Copt.callbackfun(U(:,1), J0, dUdlam);
    end
    
    %% Initial tangent (scaled space)
    dxn = Copt.startdir;
    % z = dUdlam(1:end-1, 1);
    al = dxn/vecnorm(dUdlam(:,1));
    
    alp = al;
    Ss = zeros(Copt.Nmax,2);
    
    %% BEGIN CONTINUATION
    lam  = lam0;  % current parameter
    lamp = lam0;  % previous parameter
    n    = 1;
    
    tries = 0;
    
    oopts = optimset('Jacobian', 'on', 'Display', 'off');  % options for fsolve if NSOLVE fails 

    dUn = dUdlam(:,1);  % Scaled tangent
    dUdlam(:, 1) = Copt.opts.Dscale.*dUdlam(:, 1);  % Physical tangent
    
    Ss(1, :) = [0 al];
    u0 = U(:, 1) + ds*al*(Copt.opts.Dscale.*dUn);
    while ( (lam-lam1)*(lamp-lam1) >= 0 && n<Copt.Nmax )
        if Copt.DynDscale
            Dscold = Copt.opts.Dscale;

            % Old, working, but bad.
            % % mnis = find(abs(U(:,n))>eps & abs(U(:,n))<1e-2*Copt.Dscale);
            % Copt.opts.Dscale = max(abs(U(:, n)), Copt.Dscale);
            % % Copt.opts.Dscale(mnis) = abs(U(mnis,n));

            % New
            ratio = abs(U(:,n))./Copt.opts.Dscale;
            cl_rat = min(max(ratio, 0.5), 2.0);

            Dal = 1.0;
            Copt.opts.Dscale = (1-Dal)*Copt.opts.Dscale + Dal*cl_rat.*Copt.opts.Dscale;

            dUn = (Dscold.*dUn)./Copt.opts.Dscale;
            al = sign(al)/vecnorm(dUn);
            alp = al;
        end

        switch Copt.solverchoice
          case 1
            [U(:, n+1), ~, eflag, its, J0, ~] = ...
                ELIMSEQSOLVE(func, u0, U(:, n), al*Copt.opts.Dscale.*dUn, ds, ...
                             Copt.arclengthparm, Copt.opts);
          case 2
            [U(:, n+1), ~, eflag, its, J0] = ...
                NSOLVE(@(u) EXRES(func, u, U(:, n), al*Copt.opts.Dscale.*dUn, ds, ...
                                  Copt.arclengthparm, Copt.opts.Dscale, Copt.arcsettings), ...
                       u0, Copt.opts);
          case 3
            [U(:, n+1), ~, eflag, op, J0] = ...
                fsolve(@(u) EXRES(func, u, U(:, n), al*Copt.opts.Dscale.*dUn, ds, ...
                                  Copt.arclengthparm, Copt.opts.Dscale, Copt.arcsettings), ...
                       u0, oopts);
            its = op.iterations;                  
          otherwise
            error('Unknown Solver choice')
        end
        
        if eflag<=0
            fprintf('Diverged %d - ', tries);
            switch tries
              case {0}
                if ds>Copt.dsmin && Copt.adapt==1
                    fprintf('reducing step-size\n');
                    ds = ds/2;
                    tries = tries-1;
                elseif ds==Copt.dsmin
                    fprintf('increasing step-size\n');
                    ds = ds*2;
                elseif ds~=ds0 && Copt.adapt==1
                    fprintf('resetting step-size\n')
                    ds = ds0;
                end
                u0 = U(:, n) + ds*al*(Copt.opts.Dscale.*dUn);
              case {1}
                fprintf('retrying with dsmax\n')
                  ds = Copt.dsmax;
                  u0 = U(:, n) + ds*al*(Copt.opts.Dscale.*dUn);
              case {2}
                fprintf('retrying with Branch starter\n')
                u0 = U(:, 1);
              case {3}
                fprintf('retrying with zeros\n')
                u0 = U(:, 1)*0;
              otherwise
                fprintf('Stopping\n');
                break;
            end
            tries = tries+1;
            continue;
        else
            %         [~, dR0, U(end, n+1), J0] = ELRES(func, U(1:end-1, n+1), U(:, n), al*dUn, ds, Copt.arclengthparm);
            tries = 0;
        end
        
        J0 = J0.*Copt.opts.Dscale';

        % dUdlam(:, n+1) = [-J0(1:end-1,1:end-1)\J0(1:end-1,end); 1];  % scaled tangent

        [~, kk] = sort(abs(diff(U(:, n:n+1), 1,2)), 'descend');
        for ik=1:length(kk)
            tmp = -J0(1:end-1, setdiff(1:end, kk(ik)))\J0(1:end-1,kk(ik));
            if all(isfinite(tmp))
                dUdlam([setdiff(1:end, kk(ik)) kk(ik)], n+1) = [tmp;1];
                dUdlam(:, n+1) = dUdlam(:, n+1)*sign(dUdlam(end, n+1));
                break;
            end
        end
        
        %% Callback Function
        if callback
            Scall{n+1} = Copt.callbackfun(U(:,n+1), J0, dUdlam);
        end
        
        %% Step Update
        lamp = lam;
        lam  = U(end, n+1);
        % tangent and predictor
        % z = dUdlam(1:end-1, n+1);
        
        dxn = sign(dxn*dUn'*dUdlam(:, n+1));
        %         dxn = sign((U(end,n+1)-U(end,n))*dUdlam(:, n)'*dUdlam(:, n+1));
        % al = dxn/sqrt(1+z'*z);
        al = dxn/vecnorm(dUdlam(:,n+1));

        % theta = acos(alp*al*dUn'*[z; 1]);
        theta = acos(alp*al*dUn'*dUdlam(:,n+1));

        mnis = abs(diff(U(:,n:n+1),1,2));
        mthr = graythresh(mnis)/2;

        alp_ = 1/vecnorm((mnis>mthr).*dUn)*sign(alp);
        al_ = 1/vecnorm((mnis>mthr).*dUdlam(:,n+1))*sign(al);

        % theta = acos(alp_*al_*dUn'*((mnis>mthr).*[z;1]));            

        dUn = dUdlam(:, n+1);
        dUdlam(:, n+1) = Copt.opts.Dscale.*dUdlam(:, n+1);  % tangent in physical space
        Ss(n+1,:) = [Ss(n)+ds al];
        
        %% Step size adaptation
        if Copt.Display
            fprintf('%d %f %f %e %d (%d)\n', n+1, U(end,n+1), ds, theta, its, eflag);
        end
        if n>1 && Copt.adapt==1
            if abs(theta)>1.1*Copt.angopt % angle check
                % if ds>Copt.dsmin
                    ds = max(Copt.dsmin, ds/2);
                %     continue;
                % end
            elseif abs(theta)<0.9*Copt.angopt  % angle and convergence check
                ds = min(max(ds*sqrt(Copt.itopt/its), Copt.dsmin), Copt.dsmax);
            end
        end
        alp = al;    
        n = n+1;

        u0 = U(:, n) + ds*al*(Copt.opts.Dscale.*dUn);
    end
    U = U(:, 1:n);
    dUdlam = dUdlam(:, 1:n);
    Ss = Ss(1:n, :);
    if callback
        Scall = Scall(1:n);
    end
    
    if (lam-lam1)*(lamp-lam1) < 0
        disp('Continuation completed successfully');
        flag = 1;
    else 
        disp('Premature Termination');
        flag = -1;
    end
end

%%
function [R, dRdU, dRdl] = REVRES(ul, func, dl)
    [R, dRdU] = func(ul);
    ulp = ul;
    ulp(end) = ulp(end)+dl;
    
    Rp = func(ulp);
    dRdl = (Rp-R)/dl;
end

%%
function [ul, Re, eflag, its_tot, Je, ulp0] = ELIMSEQSOLVE(func, ul, ul0, ulp0, ds, parm, opts)
%ELIMSEQSOLVE is the residual function along with the continuation constraint
%(with elimination)
%  USAGE:
%   [Re, Je] = ELIMSOLVE(func, ul, ul0, ulp0, ds, parm, opts)
%  INPUTS:
%   ul   : Solution vector (with param) - (Nd, 1)
%   ul0 : previous point (with param) - (Nd+1, 1)
%   ulp0: previous tgt (with param) - (Nd+1, 1)
%   ds  : step
%   parm: parameterization
%  OUTPUTS:
%   Re  : (Nd, 1)
%   Je  : (Nd, Nd)
%   dRe : (Nd+1, Nd+1)

    l0 = ul0(end);
    
    if strcmp(parm, 'none') || strcmp(parm, 'None')  % No Continuation, just step along
        ul(end) = l0 + sign(ul(end)-l0)*ds;
        [ul(1:end-1), ~, eflag, its_tot] = NSOLVE(@(u) func([u; ul(end)]), ul(1:end-1), opts);
        dcR = ulp0';
    else
        its_tot = 0;
        ITMAX = opts.ITMAX;
        eflag = -1;
        while its_tot<ITMAX
            % Newton update on "u" with "l" fixed
            [R, dRdU, dRdL] = func(ul); 
            dRdU = dRdU.*opts.Dscale(1:end-1)';  % scale the Jacobian
            u1 = ul(1:end-1) - opts.Dscale(1:end-1).*(dRdU\R);
            dudL = -opts.Dscale(1:end-1).*(dRdU\dRdL);  % gradient of u wrt lambda

            % Evaluate Constraint
            du = [u1; ul(end)]-ul0;
            switch parm
              case {'orthogonal' , 'Orthogonal'}
                cs = (ulp0./opts.Dscale.^2)'*du-ds;
                dcdL = (ulp0./opts.Dscale.^2)'*[dudL; 1.0];
                dcR = (ulp0./opts.Dscale.^2)';

                l1 = ul(end) - (cs/dcdL);
                %               l = ul(end);

                du = [du(1:end-1); (l1-l0)];

                coefu = (ds-ulp0(end)*(l1-l0)/opts.Dscale(end)^2)/((ulp0(1:end-1)./opts.Dscale(1:end-1).^2)'*du(1:end-1));
                coefl = 1;
              case {'arclength' , 'Arclength'}
                error('This has to be checked for consistency')
                cs = du'*du-ds^2;
                dcdL = 2*du'*[dudL; 1.0];
                dcR = 2*du';

                l1 = ul(end) - opts.Dscale(end)*(cs/dcdL);
                %               l = ul(end);              
                du = [du(1:end-1); (l1-l0)];

                coef = sqrt(ds^2/(du'*du));
                % Determine sign

                [~, mi] = min([vecnorm((ul0./opts.Dscale-coef*du)-[u1;l1]), vecnorm((ul0./opts.Dscale+coef*du)-[u1;l1])]);
                %               coefu = sqrt((ds^2-(l1-l0)^2)/(du(1:end-1)'*du(1:end-1)));

                coefu = (-1)^mi*coef;
                coefl = 1;
            end

            % Save previous point
            lp = ul(end);

            % Update
            % 1. Scale du: Doesn't work
            %       ul = ul + coef*du;

            % 2. Update lam & u through lam
            %       ul(end) = ul(end) - cs/dcdL;
            %       ul(1:end-1) = u1 + dudL*(ul(end)-lp);

            % 3. Update lam & let u update over next iteration
            %       ul(end) = ul(end) - cs/dcdL;
            %       ul(1:end-1) = u1;

            % 4. Update lam & scale du such that the next point is on manifold
            %       ul(end) = ul(end) - cs/dcdL;
            %       ul(1:end-1) = ul0(1:end-1) + du(1:end-1)*coefu;

            % 5. Scale both
            ul(end) = ul(end) - cs/dcdL;
            ul(1:end-1) = ul0(1:end-1)+coefu*du(1:end-1);
            %       ul(end) = l0 + coefl*du(end);

            its_tot = its_tot + 1;

            %           fprintf('%d %e %e\n', its_tot, abs(lp-ul(end)), abs(cs))

            % Update
            if abs((ul(end)-lp)/lp)<opts.reletol
                [ul(1:end-1), ~, eflag, its] = NSOLVE(@(u) func([u; ul(end)]), ul(1:end-1), opts);
                its_tot = its_tot+its;
                %           oopts = optimset('Jacobian', 'on', 'Display', 'off');  % options for fsolve if NSOLVE fails 
                %           [ul(1:end-1), ~, eflag, op] = fsolve(@(u) func([u; ul(end)]), ul(1:end-1), oopts);
                %           its = op.iterations;      

                break;
            else
                lp = l1;
            end
        end
    end
    
    [Re, dRdUe, dRdLe] = func(ul);
    Je = [dRdUe dRdLe; dcR];
end
%%
function [ul, Re, eflag, its_tot, Je, ulp0] = ELIMSOLVE(func, ul, ul0, ulp0, ds, parm, opts)
%ELIMSOLVE is the residual function along with the continuation constraint
%(with elimination)
%  USAGE:
%   [Re, Je] = ELIMSOLVE(func, ul, ul0, ulp0, ds, parm, opts)
%  INPUTS:
%   ul   : Solution vector (with param) - (Nd, 1)
%   ul0 : previous point (with param) - (Nd+1, 1)
%   ulp0: previous tgt (with param) - (Nd+1, 1)
%   ds  : step
%   parm: parameterization
%  OUTPUTS:
%   Re  : (Nd, 1)
%   Je  : (Nd, Nd)
%   dRe : (Nd+1, Nd+1)

    l0 = ul0(end);oopts
    lp = ul(end);

    sgn = sign(ulp0(end));

    swch = 0;
    
    its_tot = 0;
    while true
        switch parm
          case {'orthogonal' , 'Orthogonal'}
            l = l0 + (ds-ulp0(1:end-1)'*(ul(1:end-1)-ul0(1:end-1)))/ulp0(end);
            dcR = ulp0';
          case {'arclength' , 'Arclength'}
            sgn = sign(ulp0(end));              
            l = l0 + sgn*sqrt(ds^2 - (ul(1:end-1)-ul0(1:end-1))'*(ul(1:end-1)-ul0(1:end-1)));

            if ~isreal(l)
                l = lp;
            else
                dcR = 2*(ul-ul0)';
            end
          otherwise
            sgn = sign(ulp0(end));
            [Re0, dRdUe0, dRdLe0] = func(ul0);

            b = 0.5;
            c = (1-b)/(dRdLe0'*dRdUe0*dRdLe0);

            du = (ul-ul0(1:end-1));
            l = l0 + sgn*sqrt(ds^2 - c*(du'*dRdUe0*du))/b;
            
            dcR = [2*c*(du'*dRdUe0) 2*b*(l-l0)];
        end

        % Corrector 
        [un, ~, eflag, its] = NSOLVE(@(u) func([u; l]), ul(1:end-1), opts);
        
        if eflag<=0  % Try direct solve
            [unl, ~, eflag, its] = NSOLVE(@(ul) EXRES(func, ul, ul0, ulp0, ds, parm, Copt.Dscale), ul0+ulp0*ds/2, opts);
            if eflag<=0
                if swch==1
                    break;
                else
                    swch = 1;
                    ulp0 = -ulp0;
                    continue;
                end
            end
            swch = 0;
            
            un = unl(1:end-1);
            l = unl(end);
            lp = l;  % no repeats necessary if direct solve works
            sgn = sign(l-l0);
        end

        its_tot = its_tot + its;
        %       its_tot = its_tot + 1;
        
        % Update
        ul = [un; l];
        if abs((l-lp)/lp)<opts.reletol
            break;
        else
            lp = l;
        end
    end
    
    [Re, dRdUe, dRdLe] = func(ul);
    Je = [dRdUe dRdLe; dcR];
end

%%
function [Re, Je, l, dRe] = ELRES(func, u, ul0, ulp0, ds, parm)
%EXRES is the residual function along with the continuation constraint
%(with elimination)
%  USAGE:
%   [Re, Je] = ELRES(func, u, ul0, ulp0, ds, parm)
%  INPUTS:
%   u   : Solution vector (no param) - (Nd, 1)
%   ul0 : previous point (with param) - (Nd+1, 1)
%   ulp0: previous tgt (with param) - (Nd+1, 1)
%   ds  : step
%   parm: parameterization
%  OUTPUTS:
%   Re  : (Nd, 1)
%   Je  : (Nd, Nd)
%   dRe : (Nd+1, Nd+1)

    l0 = ul0(end);
    [Re0, dRdUe0, dRdLe0] = func(ul0);
    switch parm
      case {'orthogonal' , 'Orthogonal'}
        l = l0 + (ds-ulp0(1:end-1)'*(u-ul0(1:end-1)))/ulp0(end);
        dldu = -ulp0(1:end-1);
      case {'arclength' , 'Arclength'}
        sgn = sign(ulp0(end));
        l = l0 + sgn*sqrt(ds^2 - (u-ul0(1:end-1))'*(u-ul0(1:end-1)));
        
        if ~isreal(l)
            disp('hey');
        end
        
        dldu = -(u-ul0(1:end-1))'/(l-l0);
      otherwise
        sgn = sign(ulp0(end));
        b = 0.5;
        c = (1-b)/(dRdLe0'*dRdUe0*dRdLe0);
        
        du = (u-ul0(1:end-1));
        l = l0 + sgn*sqrt(ds^2 - c*(du'*dRdUe0*du))/b;
        
        dldu = -c*du'*dRdUe0/(b*(l-l0));
    end

    [Re, dRdUe, dRdLe] = func([u; l]);
    Je = dRdUe + dRdLe*dldu;
    
    if nargout>=4
        switch parm
          case {'orthogonal' , 'Orthogonal'}
            %       Re = [Re; up0'*(u-u0)-ds];
            dRe = [dRdUe dRdLe; ulp0'];
          case {'arclength' , 'Arclength'}
            %       Re = [Re; (u-u0)'*(u-u0)-ds^2];
            dRe = [dRdUe dRdLe; 2*([u; l]-ul0)'];
        end
    end
end

%%
function [Re, Je] = EXRES(func, u, u0, up0, ds, parm, Dscale, varargin)
%EXRES is the residual function along with the continuation constraint
%(without elimination)

    [Re, dRdUe, dRdLe] = func(u);
    ups = up0./Dscale; 
    switch parm
      case {'orthogonal' , 'Orthogonal'}
        %       Re = [Re; (up0./Dscale.^2)'*(u-u0)-ds];
        %       Je = [dRdUe dRdLe; (up0./Dscale.^2)'];
        Re = [Re; ups'*((u-u0)./Dscale)-ds];
        Je = [dRdUe dRdLe; (ups./Dscale)'];
      case {'arclength' , 'Arclength'}
        Re = [Re; ((u-u0)./Dscale.^2)'*(u-u0)-ds^2];
        Je = [dRdUe dRdLe; 2*((u-u0)./Dscale.^2)'];
      otherwise        
        
        [Re0, dRdUe0, dRdLe0] = func(u0);
        
        b = 0.5;
        c = (1-b)/(dRdLe0'*dRdUe0*dRdLe0);
        
        du = (u(1:end-1)-u0(1:end-1));
        dL = (u(end)-u0(end));
        Re = [Re; 
              c*(du'*dRdUe0*du)+b*dL^2-ds^2];
        Je = [dRdUe dRdLe;
              c*(2*du'*dRdUe0), 2*b*dL];
    end
end
