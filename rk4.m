% Author: Dexter Barrows

% Implementation of Runge_Kutta 4 method for integration of a system of ODEs

% f 		- function handle for the rhs of the DE
% tspan 	- vector with start end end times for the integration
% y0 		- the initial condition
% h 		- step size

function [tvals Y]  = rk4(f,tspan,y0,h)

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
		k1 = h * f(t, y);
		k2 = h * f(t + h/2, y + k1/2);
		k3 = h * f(t + h/2, y + k2/2);
		k4 = h * f(t, y + k3);
		y  = y + k1/6 + k2/3 + k3/3 + k4/6;

		t  = t + h;

		% record system state
		Y(i,:) 		= y;
		tvals(i) 	= t;

		% increment
		i 			= i + 1;
	end

	Y = Y(1:i-2,:);
	tvals = tvals(1:i-2);

end