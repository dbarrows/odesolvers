% Author: Dexter Barrows

% Implementation of Runge-Kutta-Fehlberg 4(5) method for integration of a system of ODEs

% f 		- function handle for the rhs of the DE
% tspan 	- vector with start end end times for the integration
% y0 		- the initial condition
% h0 		- initial step size
% tol 		- relative tolerance to use

function [tvals Y errvals]  = rkf45(f,tspan,y0,h0,tol)

	% Parameters (Fehlberg set)
	% --------------------------
	a1 = 1/4; 			a2 = 3/8; 			a3 = 12/13; 		a4 = 1; 			a5 = 1/2;

	b11 = 1/4;
	b21 = 3/32; 		b22 = 9/32;
	b31 = 1932/2197; 	b32 = -7200/2197; 	b33 = 7296/2197;
	b41 = 439/216 ; 	b42 = -8; 			b43 = 3680/513; 	b44 = -845/4104;
	b51 = -8/27; 		b52 = 2; 			b53 = -3544/2565; 	b54 = 1859/4104;	b55 = -11/40;

	c1 = 25/216; 		c2 = 0; 			c3 = 1408/2565; 	c4 = 2197/4104; 	c5 = -1/5;

	cs1 = 16/135; 		cs2 = 0; 			cs3 = 6656/12825; 	cs4 = 28561/56430; 	cs5 = -9/50; 	cs6 = 2/55;
	% --------------------------

	% Factors
	p_rej 	= 1/5;
	p_acc 	= 1/4;
	S 		= 0.9;
	fac_min = 0.5;
	fac_max = 1.5;

	% default absolute tolerance
	atol = 1e-13;

	t0 	= tspan(1);
	T 	= tspan(2);

	% make y0 a row vector if it is passed in as a column
	if size(y0, 1) ~= 1
   		y0 = y0'; 
	end

	% get size of system, number of steps required
	N = length(y0);
	nsteps = floor((T-t0)/h0) + 1;

	% allocate storage for results, populate initial conditions
	Y 			= zeros(nsteps,N);
	tvals 		= zeros(nsteps,1);
	errvals     = zeros(nsteps,1);
	Y(1,:) 		= y0;
	tvals(1) 	= t0;
	errvals(1)  = 0;

	% allocate storage for intermediate results
	%k = zeros(N,6);

	t = t0;
	h = h0;
	y = y0;
	i = 2;

	while t < T

		% update solution in time

		k1 = f(t, y);
		k2 = f(t + a1*h, y + k1*h*b11);
		k3 = f(t + a2*h, y + h*( k1*b21 + k2*b22 ) );
		k4 = f(t + a3*h, y + h*( k1*b31 + k2*b32 + k3*b33 ) );
		k5 = f(t + a4*h, y + h*( k1*b41 + k2*b42 + k3*b43 + k4*b44 ) );
		k6 = f(t + a5*h, y + h*( k1*b51 + k2*b52 + k3*b53 + k4*b54 + k5*b55 ) );

		yn 		= y + h*( k1*c1 + k2*c2 + k3*c3 + k4*c4 + k5*c5 );
		ynext 	= y + h*( k1*cs1 + k2*cs2 + k3*cs3 + k4*cs4 + k5*cs5  + k6*cs6 );

		err = norm(yn - ynext);
		
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
	
	tvals = tvals(1:i-2);
	errvals = errvals(1:i-2);
	Y = Y(1:(i-2),:);

end











