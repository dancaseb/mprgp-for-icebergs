% TODO's
% look at the code and identify costly operations
% implement logic for an upper bound

function [u, info] = mprgp_solver(A, b, c, opts)
    % MPRGP_SOLVER  Solve bound-constrained QP  min 1/2 x' A x - b' x  s.t. x >= c
    % Implements MPRGP (Modified Proportioning with Reduced Gradient Projection).
    %
    % Usage:
    %   opts.precond = 'none' | 'jacobi'
    %   opts.epsr = 1e-8
    %   opts.maxit = 1e5
    %   opts.Gamma = 1.0
    %   opts.verbose = true/false
    %   opts.adapt = true/false  (adaptive expansion step size)
    %   opts.bound = 'upper' | 'lower'  (currently only 'lower' implemented)
    %
    % Returns:
    %   u    - solution vector
    %   info - struct with fields: ncg, ne, np, iters, converged, final_norm_gp, runtime
    %
    % NOTES / COST remarks (inline):
    %  - COSTLY operations:
    %    * A*p  (sparse matrix-vector multiply)  <-- dominant per-iteration cost
    %    * Ap computations / matvecs inside loop
    %    * forming ichol(L) factor (one-time, can be costly)
    %    * extracting diag(A) cheap, but pinv/full is NOT used here
    %  - This implementation performs face preconditioning:
    %    z_all = M^{-1} * g  (apply preconditioner to whole g)
    %    z = z_all .* J      (only keep free-set components)
    %
    % Implementation choices:
    %  - Use normest(A) to estimate ||A|| instead of full dense norm
    %  - Robust computation of feasible step alpha_f
    %  - Safe guards for divide-by-zero in CG steps

    t0_all = tic();

    n = size(A, 1);
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'precond'), opts.precond = 'jacobi'; end
    if ~isfield(opts, 'epsr'), opts.epsr = 1e-8; end
    if ~isfield(opts, 'maxit'), opts.maxit = 1e5; end
    if ~isfield(opts, 'Gamma'), opts.Gamma = 1.0; end

    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'adapt'), opts.adapt = false; end
    if ~isfield(opts, 'bound'), opts.bound = 'lower'; end

    if opts.bound == 'upper'
        bs = -1; %bound sign
    else
        bs = 1;
    end

    disp('bs'); disp(bs);
    % bs=1;

    % column vectors
    b = b(:);
    c = c(:);

    % Initial feasible point: start at the bounds
    % if bs == 1
    %     u = max(zeros(n, 1), c);
    % else
    %     u = max(zeros(n, 1), c);
    % end
    u = c;

    lAl = normest(A); % estimate of ||A||_2, cheaper than full norm
    if lAl == 0, lAl = 1; end
    alpha = 1 / lAl; % recommended expansion step size

    % preconditining logic
    precond_type = lower(opts.precond);

    switch precond_type
        case 'none'
            precond_apply = @(g, J) (g .* J); % identity (no preconditioning), then mask free set
            if opts.verbose, disp('Preconditioner: NONE'); end

        case 'jacobi'
            D = diag(A);
            D(abs(D) < eps) = 1; % avoid division by zero
            precond_apply = @(g, J) ((g ./ D) .* J); % z = (diag(A)^{-1} g) on free set
            if opts.verbose, disp('Preconditioner: Jacobi (diagonal)'); end

        otherwise
            error('Unknown preconditioner: %s', opts.precond);
    end

    % initialization
    g = A * u - b; % gradient (Cost: matvec operation)
    disp(' u '); disp(u(1:10));
    disp('A * u'); disp((A * u)(1:10));
    disp('g')
    disp(g(1:10));
    J = (bs .* u > bs .* c); % logical free-set indicator
    % disp('J:');disp(J);
    disp('size J:'); disp(size(J));
    disp('num free:'); disp(sum(J));
    disp('g norm:'); disp(norm(g));
    gf = J .* g; % free gradient

    if bs == 1
        gc = min((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
    end

    if bs == -1
        gc = max((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
    end

    if bs == 1
        gr = min(lAl * (J .* (u - c)), gf); % reduced free gradient (O(n)), phi tilde in paper
    else
        gr = max(lAl * (J .* (u - c)), gf); % reduced free gradient (O(n)), phi tilde in paper
    end
    % probably replace this with a max? in upper bound case. Not sure.
    gp = gf + gc; % projected gradient
    disp('gf, gc, gp norms:'); disp([norm(gf), norm(gc), norm(gp)]);

    z = precond_apply(g, J);
    p = z; % search direction

    % counters & info
    ncg = 0; ne = 0; np = 0;
    iters = 0;
    converged = false;

    % main MPRGP loop
    while norm(gp) > opts.epsr && iters < opts.maxit
        disp(norm(gp));
        iters = iters + 1;

        if (gc' * gc) <= (opts.Gamma^2) * (gr' * gf)
            Ap = A * p; % Cost: matvec operation
            rtp = z' * g;
            pAp = p' * Ap;

            if abs(pAp) < eps
                % direction degenerate; break to avoid division by zero
                if opts.verbose, warning('p''*A*p nearly zero, stopping'); end
                break;
            end

            acg = rtp / pAp; % CG step size
            yy = u - acg * p; % full CG step (not necessarily feasible)

            if all(bs .* yy >= bs .* c)
                % Full CG step accepted
                u = yy;
                g = g - acg * Ap; % update gradient
                J = (bs .* u > bs .* c);
                z = precond_apply(g, J);
                beta = (z' * Ap) / pAp;
                p = z - beta * p;
                % update gradients
                gf = J .* g;

                if bs == 1
                    gc = min((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
                end

                if bs == -1
                    gc = max((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
                end

                if bs == 1
                    gr = min(lAl * (J .* (u - c)), gf);
                else
                    gr = max(lAl * (J .* (u - c)), gf);
                end
                gp = gf + gc;
                ncg = ncg + 1;
            else
                % expansion step -> compute feasible step alpha_f robustly
                pos = (bs .* p > 0) & J; % only indices where p would decrease u towards bound (positive p)

                if any(pos)
                    a_vals = (u(pos) - c(pos)) ./ p(pos);
                    a = min(a_vals); % a_f in paper
                    if ~isfinite(a) || a < 0, a = 0; end
                else
                    a = 0;
                end

                if bs == 1
                    u = max(u - a * p, c); % enforce feasibility
                else
                    u = min(u - a * p, c); % enforce feasibility
                end
                J = (bs .* u > bs .* c);
                g = g - a * Ap; % update gradient
                % recompute preconditioned residual
                z = precond_apply(g, J);
                % Expansion step: projection along -alpha*gf (expansion/projection move)
                % Use gf as J.*g
                % adaptive step size
                if opts.adapt
                    alpha = (gr' * g) / (gr' * (A * gr));

                    if alpha <= 0 || alpha > 1 / lAl
                        alpha = 1 / lAl;
                    end

                end

                if bs == 1
                    u = max(u - alpha * (J .* g), c);
                else
                    u = min(u - alpha * (J .* g), c);
                end
                J = (bs .* u > bs .* c);
                g = A * u - b; % Cost (matvec)
                z = precond_apply(g, J);
                p = z;
                gf = J .* g;

                if bs == 1
                    gc = min((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
                end

                if bs == -1
                    gc = max((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
                end

                if bs ==1
                    gr = min(lAl * (J .* (u - c)), gf);
                else
                    gr = max(lAl * (J .* (u - c)), gf);
                end
                gp = gf + gc;
                ne = ne + 1;
            end

        else
            % proportioning step: move along chopped gradient gc
            Ap = A * gc; % COSTLY: matvec
            denom = gc' * Ap;

            if abs(denom) < eps
                if opts.verbose, warning('denominator in proportioning nearly zero, stopping'); end
                break;
            end

            acg = (gc' * g) / denom;
            u = u - acg * gc;
            if bs == 1
                u = max(u, c); % enforce feasibility numerically
            else
                u = min(u, c); % enforce feasibility numerically
            end
            J = (bs .* u > bs .* c);
            g = g - acg * Ap;
            z = precond_apply(g, J);
            p = z;
            gf = J .* g;

            if bs == 1
                gc = min((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
            end

            if bs == -1
                gc = max((~J) .* g, 0); % chopped gradient (O(n)), beta in paper
            end

            if bs == 1
                gr = min(lAl * (J .* (u - c)), gf);
            else
                gr = max(lAl * (J .* (u - c)), gf);
            end

            gp = gf + gc;
            np = np + 1;
        end

        if iters == opts.maxit
            if opts.verbose, warning('Reached maxit'); end
        end

        % disp(norm(gp));
    end

    info.ncgs = ncg;
    info.nes = ne;
    info.nps = np;
    info.iters = iters;
    info.converged = (norm(gp) <= opts.epsr);
    info.final_norm_gp = norm(gp);
    info.runtime = toc(t0_all);

    if opts.verbose
        fprintf('MPRGP finished: iters=%d, ncg=%d, ne=%d, np=%d, norm_gp=%.3e, time=%.3fs\n', ...
            info.iters, info.ncgs, info.nes, info.nps, info.final_norm_gp, info.runtime);
    end

end
