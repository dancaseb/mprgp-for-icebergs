function [u, info] = mprgp_solver(A, b, c, opts)
    % MPRGP_SOLVER  Solve bound-constrained QP  min 1/2 x' A x - b' x  s.t. x >= c
    % Implements MPRGP (Modified Proportioning with Reduced Gradient Projection).
    %
    % Usage:
    %   opts.precond = 'none' | 'jacobi' | 'ichol'   (default 'jacobi')
    %   opts.epsr = 1e-8
    %   opts.maxit = 1e5
    %   opts.Gamma = 1.0
    %   opts.verbose = true/false
    %
    % Returns:
    %   u    - solution vector (feasible: u >= c)
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
    %
    % References: MPRGP algorithm (see provided PDF/chapter) and face preconditioning.

    t0_all = tic();

    n = size(A, 1);
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'precond'), opts.precond = 'jacobi'; end
    if ~isfield(opts, 'epsr'), opts.epsr = 1e-8; end
    if ~isfield(opts, 'maxit'), opts.maxit = 1e5; end
    if ~isfield(opts, 'Gamma'), opts.Gamma = 1.0; end

    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'adapt'), opts.adapt = false; end

    % Ensure column vectors
    b = b(:);
    c = c(:);

    % Initial feasible point: projection of zero onto bounds
    % u = max(zeros(n, 1), c);
    u = c;  % Start at the lower bound

    % Estimate matrix norm (COST: normest is cheap-ish, avoids full -> dense)
    lAl = normest(A); % estimate of ||A||_2
    if lAl == 0, lAl = 1; end
    alpha = 1 / lAl; % recommended expansion step size

    % Prepare preconditioner (one-time cost)
    precond_type = lower(opts.precond);

    switch precond_type
        case 'none'
            precond_apply = @(g, J) (g .* J); % identity (no preconditioning), then mask free set
            if opts.verbose, disp('Preconditioner: NONE'); end

        case 'jacobi'
            D = diag(A); % COST: extracting diagonal O(n)
            D(abs(D) < eps) = 1; % avoid division by zero
            precond_apply = @(g, J) ((g ./ D) .* J); % z = (diag(A)^{-1} g) on free set
            if opts.verbose, disp('Preconditioner: Jacobi (diagonal)'); end

        % case 'ichol'
        %     % Try incomplete Cholesky; fall back to Jacobi if it fails
        %     try
        %         % COST: ichol factorization is one-time and can be expensive for big A.
        %         % Tweak options if needed (drop tolerance, diagonal compensation)
        %         % Note: ichol expects symmetric positive definite-ish A
        %         opts_ichol.type = 'ict';
        %         opts_ichol.droptol = 1e-3; % tweakable
        %         L = ichol(A, opts_ichol); % COSTLY: one-time factorization (sparse)
        %         precond_apply = @(g, J) ((L' \ (L \ g)) .* J); % z = M^{-1} g on free set
        %         if opts.verbose, disp('Preconditioner: ICHOL (succeeded)'); end
        %     catch ME
        %         warning('ichol failed (%s). Falling back to Jacobi preconditioner.', ME.message);
        %         D = diag(A); D(abs(D) < eps) = 1;
        %         precond_apply = @(g, J) ((g ./ D) .* J);
        %         if opts.verbose, disp('Preconditioner: Jacobi (fallback)'); end
        %     end

        otherwise
            error('Unknown preconditioner: %s', opts.precond);
    end

    % Initial gradient and projected quantities
    g = A * u - b; % gradient (COST: one matvec if u not zero; here A*u cheap for initial u)
    J = (u > c); % logical free-set indicator
    gf = J .* g; % free gradient (O(n))
    gc = min((~J) .* g, 0); % chopped gradient (O(n))

    % gr = min(lAl * (J .* (u - c)), gf); % reduced free gradient (O(n)), phi tilde in paper
    gr = compute_reduced_gradient(gf, J, u, c, lAl, n);

    gp = gf + gc; % projected gradient

    % Preconditioned residual and CG direction
    z = precond_apply(g, J); % preconditioner apply (cost depends on type)
    p = z;

    % counters & info
    ncg = 0; ne = 0; np = 0;
    iters = 0;
    converged = false;

    % Main MPRGP loop
    while norm(gp) > opts.epsr && iters < opts.maxit
        iters = iters + 1;

        % Strict proportionality test (COST: dot prods O(n))
        if (gc' * gc) <= (opts.Gamma^2) * (gr' * gf)
            % Proportional branch: trial preconditioned CG step
            Ap = A * p; % COSTLY: matrix-vector product (dominant per-iteration cost)
            rtp = z' * g; % preconditioned inner product (O(n))
            pAp = p' * Ap; % inner product (O(n))

            if abs(pAp) < eps
                % direction degenerate; break to avoid division by zero
                if opts.verbose, warning('p''*A*p nearly zero, stopping'); end
                break;
            end

            acg = rtp / pAp;
            yy = u - acg * p;

            if all(yy >= c)
                % Full CG step accepted
                u = yy; % O(n)
                g = g - acg * Ap; % O(n)
                % Update preconditioned residual and direction
                J = (u > c);
                z = precond_apply(g, J); % preconditioner apply (COST depends on type)
                % Beta update uses z (preconditioned residual)
                beta = (z' * Ap) / pAp; % O(n)
                p = z - beta * p; % O(n)
                % Update gradients
                gf = J .* g; gc = min((~J) .* g, 0); 

                % gr = min(lAl * (J .* (u - c)), gf);
                gr = compute_reduced_gradient(gf, J, u, c, lAl, n);
                
                gp = gf + gc;
                ncg = ncg + 1;
            else
                % Partial (feasible) CG step -> compute feasible step alpha_f robustly
                pos = (p > 0) & J; % only indices where p would decrease u towards bound (positive p)

                if any(pos)
                    a_vals = (u(pos) - c(pos)) ./ p(pos);
                    a = min(a_vals);
                    if ~isfinite(a) || a < 0, a = 0; end
                else
                    a = 0;
                end

                u = max(u - a * p, c); % enforce feasibility
                J = (u > c);
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
                    %  TODO should lAl be updated here?

                end

                u = max(u - alpha * (J .* g), c); % COST: O(n)
                J = (u > c);
                g = A * u - b; % full matvec to refresh gradient (COSTLY)
                z = precond_apply(g, J);
                p = z;
                gf = J .* g; gc = min((~J) .* g, 0); gr = min(lAl * (J .* (u - c)), gf);
                gp = gf + gc;
                ne = ne + 1;
            end

        else
            % Proportioning step: move along chopped gradient gc
            Ap = A * gc; % COSTLY: matvec
            denom = gc' * Ap;

            if abs(denom) < eps
                if opts.verbose, warning('denominator in proportioning nearly zero, stopping'); end
                break;
            end

            acg = (gc' * g) / denom;
            u = u - acg * gc;
            u = max(u, c); % enforce feasibility numerically
            J = (u > c);
            g = g - acg * Ap;
            z = precond_apply(g, J);
            p = z;
            gf = J .* g; gc = min((~J) .* g, 0); gr = min(lAl * (J .* (u - c)), gf);
            gp = gf + gc;
            np = np + 1;
        end

        if iters == opts.maxit
            if opts.verbose, warning('Reached maxit'); end
        end

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

function gr = compute_reduced_gradient(gf, J, u, c, alpha, n)
    % Compute reduced free gradient: gr = min(lAl * (J .* (u - c)), gf)
    
    % --- reduced free gradient (phi_tilde) -- compute exactly elementwise on free set
    phi = gf;  % free gradient (zero on active)
    phi_tilde = zeros(n,1);
    if any(J)
        % formula: phi_tilde_i = min( (u_i - c_i)/alpha, phi_i) for i in free set
        phi_tilde(J) = min( (u(J) - c(J)) / alpha, phi(J) );
    end
    gr = phi_tilde;
end