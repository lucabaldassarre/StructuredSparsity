% FUNCTION: [fsol, output] = basisPursuitSolver(Aopr, ATopr, b, W, lbx, ...
%                                               ubx, param, normA2, x0)
%
% PURPOSE: Solve the basis pursuit problem of the form:
%
%                       min_x |x|_1
%                       s.t.  A*x == b, lbx <= x <= ubx.
%
% INFORMATION:
%    By Quoc Tran Dinh, Laboratory for Informations and Inference Systems,
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 24.01.2014.
%    Last modified: 24.01.2014.
%    Contact: quoc.trandinh@epfl.ch
%   
function [xsol, output] = tvNormProxSolver(Aopr, ATopr, b, W, lbx, ...
                                           ubx, param, x0)
                                       
% Evaluate the norm of A'*A.
param    = lassoParamSettings(param); 
normA2   = lassoNormAtAeval('operator', Aopr, param.PwMaxIters, ...
                             param.PwRelTol, ATopr, length(x0));
                         
% Initialize the parameters.                                      
LipsA2   = normA2;
rho      = min( max(0.1/normA2, 1.0e-6), 1.0 ); 
tau      = 0.5*(sqrt(5) - 1);
Lg_ub    = LipsA2 + 1;
beta     = Lg_ub*tau^2/(rho*(1 - tau));
norm_b   = norm(b, 2);
irhoL    = 1.0/(rho*LipsA2);
gamFact  = 0.00;
incrGam  = 1.06;
decrGam  = 1.01;
maxGamma = 1.0e+8;
minGamma = 1.0e-8;

% Initialize the outputs.
output   = []; 
cntA     = 0;
cntAt    = 0; 
fx_val   = nan;
%x_zeros  = zeros(size(b));
radius   = 10;

% Define the center points.
x_cen    = x0;
y_cen    = zeros(size(b));
Axs      = Aopr(x_cen);
r_cen    = Axs - b;

% Define the primal point xs and rs.
xs       = x_cen;
if ~isempty(lbx), xs = max(lbx, xs); end
if ~isempty(ubx), xs = min(xs, ubx); end
rs       = r_cen;

% Compute the primal initial point xb.
%xb       = sign(x_cen).*max(abs(x_cen) - irhoL.*W, 0);
m1       = sqrt(length(x0));
n1       = m1;
xb       = chambolle_prox_TV_stop(reshape(x_cen, [m1, n1]), 'lambda', irhoL, 'maxiter', 10);
xb       = xb(:);

% Compute the primal residual rb.
nrs      = norm(rs, 2);
rb       = radius*min(1.0, 1.0/nrs).*rs;
%rb       = x_zeros;


% Update the dual variable yb.
Axb      = Aopr(xb);
cntA     = cntA + 1;
frstFeas = Axb - rb - b;
yb       = y_cen + (1.0/beta)*frstFeas;
sndFeas  = Axs - rs - b;

% The main loop of the algorithm.
for iter = 1:param.MaxIters
    
    % STEP 1: Evaluate the first dual variable ys.
    frstFeas  = Axb - rb - b;
    ybs       = y_cen + (1.0/beta)*frstFeas;
    
    % STEP 2: Compute the initermediate point yh.
    yh        = (1.0 - tau)*yb + tau*ybs;
    
    % STEP 3: Evaluate the objective value.
    if param.isFxEval
        fx_val = norm(W.*xb, 1);
    end
    
    % STEP 4: Compute the primal solution xs.
    iLips     = 1.0/LipsA2;
    irhoL     = 1.0/(rho*LipsA2);
    ss        = xs - ATopr(irhoL*yh + iLips*sndFeas);
    cntAt     = cntAt + 1;
    %xs_next   = sign(ss).*max(abs(ss) - irhoL.*W, 0);
    xs_next   = chambolle_prox_TV_stop(reshape(ss, [m1, n1]), 'lambda', irhoL, 'maxiter', 10);
    xs_next   = xs_next(:);
    if ~isempty(lbx), xs_next = max(lbx, xs_next); end
    if ~isempty(ubx), xs_next = min(xs_next, ubx); end
    Axs_next  = Aopr(xs_next);
    cntA      = cntA + 1;
    
    % STEP 5: Compute the primal residual rs.
    %rs_next   = x_zeros;
    irho      = 1.0/rho;
    %rcs       = Axs - b + irho*yh;
    rcs       = Axs_next - b + irho*yh;
    nrs       = norm(rcs(:), 2);
    radius    = max(0.95*radius, param.RelTolX);
    rs_next   = radius*min(1.0, 1.0/nrs).*rcs;
    
    % STEP 6: Update the dual variable.
    step      = rho;
    sndFeas   = Axs_next - rs_next - b;
    yb_next   = yb + step*sndFeas;
    %step      = rho/Lg_ub;
    %yb_next   = yh + step*sndFeas;
    
    % STEP 7: Update xb_next and rb_next.
    xb_next   = (1 - tau)*xb + tau*xs_next;
    rb_next   = (1 - tau)*rb + tau*rs_next;
    
    % STEP 8: Evaluate the primal and dual feasibility.
    abs_pfeas = norm(sndFeas, 2);
    abs_dfeas = rho*norm(Axs_next - Axs, 2);
    cntAt     = cntAt + 1;
    
    rel_pfeas = abs_pfeas/max(1.0, norm_b);
    rel_dfeas = abs_dfeas/max(1.0, norm(yb, 2));
    rel_schg  = norm(xb_next - xb)/max(norm(xb, 2), 1.0);
    
    % STEP 9a: Update the smoothness parameter beta.
    beta      = (1 - tau)*beta;
    
    % STEP 9b: Update the smoothness parameter rho.
    updateSmoothnessParameter();
    
    % STEP 10: Save the history and print the iterations.
    if param.saveHistMode > 0, saveHistoryList;                 end
    if param.Verbosity >= 2,   printIteration(param.PrintStep); end
    
    % STEP 11: Check the stopping criterion.
    if rel_pfeas <= param.RelTolX && rel_schg <= param.RelTolX
        performFinalPhase(1);
        break;
    end
    
    % STEP 12: Update the parameter tau.
    omeg  = rho_next/rho*tau^2;
    tau   = 0.5*(sqrt(omeg^2 + 4*omeg) - omeg);
    rho   = rho_next;
    
    % STEP 13: Assign to the next iteration.
    rs    = rs_next;
    xs    = xs_next;
    yb    = yb_next;
    xb    = xb_next;
    rb    = rb_next;
    Axs   = Axs_next;
    Axb   = (1.0 - tau)*Axb + tau*Axs;
    
end
% End of the main loop.

% Perform the final phase.
performFinalPhase(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              THE IMPLEMENTATION OF THE NEST FUNCTIONS
%%%              ****************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -----------------------------------------------------------------------
%%% Function: updateSmoothnessParameter()
%%% Purpose:  Update the smoothness parameter.
%%% -----------------------------------------------------------------------
    function updateSmoothnessParameter()
        
        if abs_pfeas >= gamFact*abs_dfeas
            rho_next = rho*incrGam;
        elseif abs_dfeas >= gamFact*abs_pfeas
            rho_next = rho/decrGam;
        else
            rho_next = rho;
        end
        rho_next = min(max(rho_next, minGamma), maxGamma);
        
    end
%%% -----------------------------------------------------------------------
%%% Function: printIteration(step)
%%% Purpose:  Print the output at each iteration.
%%% -----------------------------------------------------------------------
    function printIteration(step)
        
        % Print the header.
        if mod(iter, 10*step) == 1
            fprintf('%s\n', repmat('-', 1, 75)); 
            fprintf([' Iter| RelPGap | RelDGap | RelSchg |', ...
                     '  Beta  |  Rho   |   Tau  |    F(x)\n']);
            fprintf('%s\n', repmat('-', 1, 75));
        end

        % Print the values.
        if mod(iter, step) == 0
            fprintf(['%5d| %3.2e| %3.2e| %3.2e| %3.1e| %3.1e|', ...
                     ' %3.1e| %3.5e \n'], iter, rel_pfeas, rel_dfeas, ...
                     rel_schg, beta, rho, tau, fx_val);
        end
        
    end
%%% -----------------------------------------------------------------------
%%% Function: saveHistoryList
%%% Purpose:  Save the history information if required.
%%% -----------------------------------------------------------------------
    function saveHistoryList()
        
       output.hist.rel_schg( iter, 1) = rel_schg;
       output.hist.rel_pfeas(iter, 1) = rel_pfeas;
       output.hist.rel_dfeas(iter, 1) = rel_dfeas;
       output.hist.abs_pfeas(iter, 1) = abs_pfeas;
       output.hist.abs_dfeas(iter, 1) = abs_dfeas;
       output.hist.fx_val(  iter, 1)  = fx_val;
       
       if param.saveHistMode > 1
           output.hist.beta(iter, 1)  = beta;
           output.hist.rho( iter, 1)  = rho;
           output.hist.tau( iter, 1)  = tau;
       end
       
       if param.saveHistMode > 2
           output.hist.xb{iter}       = xb;
           output.hist.rb{iter}       = rb;
           output.hist.yb{iter}       = yb;
           output.hist.xs{iter}       = xs;
           output.hist.xs{iter}       = rs;
       end
       
    end
%%% -----------------------------------------------------------------------
%%% Function: performFinalPhase(modes)
%%% Purpose:  Perform the final phase after terminating the algorithm.
%%% -----------------------------------------------------------------------
    function performFinalPhase(modes)
        
        if modes == 1
            output.status = 'Convergence achieved';
            output.msg = ['Feasibility gap and diffences of solutions', ...
                          ' are below the given tolerance and the ',...
                          'search direction does not change signficantly'];
        elseif modes == 2
            output.status = 'Convergence achieved';
            output.msg = ['Feasibility gap is below the given tolerance',...
                          ' and the objective does not change signficantly'];
        elseif modes == 0 && iter >= param.MaxIters
            output.status = 'Have not converged yet';
            output.msg = ['Exceed the maximum number of iterations. ', ...
                          'Increase MaxIters if required!'];
        end
        
        % Get the other outputs.
        output.iter     = iter;
        output.rel_pgap = rel_pfeas;
        output.rel_dgap = rel_dfeas;
        output.rel_schg = rel_schg;
        output.fx_val   = fx_val;
        output.cntA     = cntA;
        output.cntAt    = cntAt;
        output.auxi.xs  = xs;
        output.auxi.rb  = rb_next;
        output.auxi.yb  = yb_next;
        
        % Get the final solution.
        xsol = xb_next;
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF THE IMPLEMENTATION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

