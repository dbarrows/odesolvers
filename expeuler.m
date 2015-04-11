% Author: Dexter Barrows
% Implementation of the Explicit Euler method for integration of a system of ODEs

% f 		- function handle for the rhs of the DE
% tspan 	- vector with start end end times for the integration
% y0 		- the initial condition
% h 		- step size

function [tvals Y] = expeuler(f,tspan,y0,h)

	t0 	= tspan(1);
	T 	= tspan(2);

	% make y0 a column vector if it is passed in as a row
	if size(y0, 1) ~= 1
   		y0 = y0'; 
	end

	% get size of system, number of steps required
	N = length(y0);
	nsteps = floor((T-t0)/h) + 1;

	% allocate storage for results, populate initial conditions
	Y 			= zeros(nsteps,N);
	tvals 		= zeros(nsteps,1);
	Y(1,:) 		= y0;
	tvals(1) 	= t0;

	t = t0;
	y = y0;
	i = 2;
	while t < T
		% update solution in time
		y 			= y + h * f(t,y);
		t 			= t + h;

		% record system state
		Y(i,:) 		= y;
		tvals(i) 	= t;

		% increment
		i 			= i + 1;
	end

	% truncate extra term
	Y = Y(1:i-2,:);
	tvals = tvals(1:i-2);

end