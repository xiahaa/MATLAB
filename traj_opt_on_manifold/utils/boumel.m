function varargout = boumel(varargin)
    Rreg = varargin{1};

    n = 3;
    N = size(Rreg,3);
    p = Rreg;

    % For each control point, pick a weight (positive number). A larger value
    % means the regression curve will pass closer to that control point.
    w = ones(N, 1);

    %% Define parameters of the discrete regression curve
    % The curve has Nd points on SO(n)
    if nargin == 1
        Nd = N;
    else
        Nd = varargin{2};
    end
    % Each control point attracts one particular point of the regression curve.
    % Specifically, control point k (in 1:N) attracts curve point s(k).
    % The vector s of length N usually satsifies:
    % s(1) = 1, s(end) = Nd and s(k+1) > s(k).
    s = round(linspace(1, Nd, N));

    % Time interval between two discretization points of the regression curve.
    % This is only used to fix a scaling. It is useful in particular so that
    % other parameter values such as w, lambda and mu (see below) have the same
    % sense even when the discretization parameter Nd is changed.
    delta_tau = 1/(Nd-1);

    % Weight of the velocity regularization term (nonnegative). The larger it
    % is, the more velocity along the discrete curve is penalized. A large
    % value usually results in a shorter curve.
    lambda = 0;

    % Weight of the acceleration regularization term (nonnegative). The larger
    % it is, the more acceleration along the discrete curve is penalized. A
    % large value usually results is a 'straighter' curve (closer to a
    % geodesic.)
    mu = 1;

    %% Pack all data defining the regression problem in a problem structure.
    problem.n = n;
    problem.N = N;
    problem.Nd = Nd;
    problem.p = p;
    problem.s = s;
    problem.w = w;
    problem.delta_tau = delta_tau;
    problem.lambda = lambda;
    problem.mu = mu;

    %% Call the optimization procedure to compute the regression curve.

    % Compute an initial guess for the curve. If this step is omitted, digress
    % (below) will compute one itself. X0 is a 3D matrix of size n x n x Nd,
    % such that each slice X0(:, :, k) is a rotation matrix.
    %
    if nargin == 1
        X0 = Rreg;
    else
        X0 = initguess(problem);
    end
    % Run the optimization procedure to compute X1, the discrete regression
    % curve. X1 is a 3D matrix of size n x n x Nd with each slice a rotation
    % matrix. The second output, info, is a struct-array containing information
    % returned by the optimization algorithm. The third output, optim_problem,
    % is the Manopt optimization problem structure used to produce X1. It can
    % be used to run another algorithm, e.g., for research purposes.
    %
    [X1, info, optim_problem] = digress(problem, X0);

    %% Plot optimization information
    figure(3);
    semilogy([info.time], [info.gradnorm], 'k.-');
    title('Gradient norm');
    xlabel('Computation time [s]');
    pbaspect([1.6, 1, 1]);


    %% Plot speed and acceleration of X0 and X1

    [speed0, acc0] = compute_profiles(problem, X0);
    [speed1, acc1] = compute_profiles(problem, X1);

    % Passage time of each point on the discrete curves.
    time = problem.delta_tau*( 0 : (problem.Nd-1) );

    figure(4);

    subplot(2, 1, 1);
    plot(time, speed0, time, speed1);
    title('Speed of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Speed');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

    subplot(2, 1, 2);
    plot(time, acc0, time, acc1);
    title('Acceleration of initial curve and optimized curve');
    xlabel('Time');
    ylabel('Acceleration');
    legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
    pbaspect([1.6, 1, 1]);

    varargout{1} = X1;
    if nargout == 2
        varargout{2} = info.cost;
    end
    
%     ylim([0, 20]);
end