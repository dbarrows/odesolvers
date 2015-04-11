% Author: Dexter Barrows

% Implementation of modified Rosenbrock 2(3) method for integration of a system of ODEs

% f 		- function handle for the rhs of the DE
% jac 		- function handle for the jacobian of the rhs of the DE
% tspan 	- vector with start end end times for the integration
% y0 		- the initial condition
% h0 		- initial step size
% tol 		- relative tolerance to use

function [tvals Y errvals]  = rosen23mod(f,jac,tspan,y0,h0,tol)


	% Parameters
	% ------------------------------------------------

	d = 1/(2 + 2^(1/2));

	e32 = 6 + 2^(1/2);

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

	rej_flag = 1;
	while t < T

		W = eye(N) - h*d*jac(t,y);
		Winv = inv(W);

		if rej_flag 		% last setp was rejected -> F0 must be recalculated
			F0 = f(t,y);					% row
		else 				% last step was successful -> we can take advantage of the fact that F0 = F2
			F0 = F2;
		end
		k1 = Winv * F0';				% column
		F1 = f(t + h/2, y + h*k1'/2);	% row
		k2 = Winv * (F1'-k1) + k1; 		% column

		ynext = y + h*k2';

		F2 = f(t+h,ynext); 				% row
		k3 = Winv * (F2' - e32*(k2-F1') - 2*(k1-F0'));

		err = norm(h*(k1 - 2*k2 + k3)/6);
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



















