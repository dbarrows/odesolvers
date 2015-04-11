% Author: Dexter Barrows

% Implementation of 4th order Rosenbrock method for integration of a system of ODEs

% f 		- function handle for the rhs of the DE
% jac 		- function handle for the jacobian of the rhs of the DE
% tspan 	- vector with start end end times for the integration
% y0 		- the initial condition
% h0 		- initial step size
% tol 		- relative tolerance to use

function [tvals Y errvals]  = rosenbrock(f,jac,tspan,y0,h0,tol)


	% Constants
	% ------------------------------------------------

	gam = 1/2;

	a21 = 2;
	a31 = 48/25; 	a32 = 6/25;

	c21 = -8;
	c31 = 372/25; 	c32 = 12/5;
	c41 = -112/125; c42 = -54/125; c43 = -2/5;

	b1 = 19/9; b2 = 1/2; b3 = 25/108; b4 = 125/108;

	e1 = 17/54; e2 = 7/36; e3 = 0; e4 = 125/108;

	c1x = 1/2; c2x = -3/2; c3x = 121/50; c4x = 29/250;

	a2x = 1; a3x = 3/5;

	% ------------------------------------------------

	% Factors
	p_rej 	= 1/4;
	p_acc 	= 1/3;
	S 		= 0.99999;

	% times
	t0 	= tspan(1);
	T 	= tspan(2);

	% factors
	fac_min = 0.5;
	fac_max = 1.5;

	atol = 1e-13;

	% make y0 a row vector if it is passed in as a column
	if size(y0, 1) ~= 1
   		y0 = y0'; 
	end

	h = h0;
	y = y0;
	t = t0;

	% get size of system, number of steps required
	N = length(y0);
	nsteps = floor((T-t0)/h0) + 1;

	% allocate storage for results, populate initial conditions
	Y 			= zeros(nsteps,N);
	tvals 		= zeros(nsteps,1);
	errvals 	= zeros(nsteps,1);
	Y(1,:) 		= y0;
	tvals(1) 	= t0;

	i = 2; % starting on second step

	while t < T

		A = 1/(gam*h)*eye(N) - jac(t,y);

		b = f(t,y);
		g1 = A\b';
		g1 = g1';

		b = f(t + a2x*h, y + a21*g1) + c21*g1/h;
		g2 = A\b';
		g2 = g2';

		dy_common = f(t + a3x*h,y + a31*g1 + a32*g2);

		b = dy_common + (c31*g1 + c32*g2)/h;
		g3 = A\b';
		g3 = g3';

		b = dy_common + (c41*g1 + c42*g2 + c43*g3)/h;
		g4 = A\b';
		g4 = g4';

		ynext = y + b1*g1 + b2*g2 + b3*g3 + b4*g4;

		err = norm(e1*g1 + e2*g2 + e4*g4);
		err_tar = atol + norm(y) .* tol;
		
		% If error is withing desired parameters, advance solution.
		% Otherwise try again with reduced step size.
		if err > err_tar
			rej_flag = 1;
		else
			rej_flag = 0;
		end

		if rej_flag
			fac = S * (err_tar/err)^p_rej;
			fac = max(min(fac,fac_max),fac_min);
			h = h * fac;
		else
			% advance solution
			y = ynext;
			t = t + h;

			% record system state
			Y(i,:) 		= y;
			tvals(i) 	= t;
			errvals(i) 	= err;

			% increment
			i = i + 1;

			fac = S * (err_tar/err)^p_acc;
			fac = max(min(fac,fac_max),fac_min);
			h = h * fac;
		end

	end

	Y = Y(1:(i-2),:);
	tvals = tvals(1:(i-2));
	errvals = errvals(1:(i-2));


end



















