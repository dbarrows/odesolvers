% Author: Dexter Barrows

% This function shows how to call the ODE solvers

function main

	% set up integration of Van der Pol system

	T = 1000;
	tspan = [0 T];
	h = 1e-3;
	tol = 1e-4;
	in_con = [2 -2/3];

	% integrate with various methods and time

	tic
	[t Y] = expeuler(@vdpr,tspan,in_con,h);
	expeuler_time = toc;
	expeuler_steps = length(t);

	tic
	[t1 Y1] = rk4(@vdpr,tspan,in_con,h);
	rk4_time = toc;
	rk4_steps = length(t1);

	tic
	[t2 Y2] = rkf45(@vdpr,tspan,in_con,h,tol);
	rkf45_time = toc;
	rkf45_steps = length(t2);

	tic
	[t3 Y3] = rosenbrock(@vdpr,@vdpjac,tspan,in_con,h,tol);
	rosenbrock_time = toc;
	rosenbrock_steps = length(t3);

	tic
	[t4 Y4] = rosen23mod(@vdpr,@vdpjac,tspan,in_con,h,tol);
	rosen23mod_time = toc;
	rosen23mod_steps = length(t4);
	
	options = odeset('Abstol',1e-13,'Reltol',tol);
	tic
	[t5 Y5] = ode23s(@vdpc,tspan,in_con,options);
	ode23s_time = toc;
	ode23s_steps = length(t5);

	% plot the approximation for the first species

	figure
	plot(t,Y(:,1),t1,Y1(:,1),'-',t2,Y2(:,1),'.-',t3,Y3(:,1),'o-',t4,Y4(:,1),'*-', t5,Y5(:,1),'+-');
	lh = legend('Explicit Euler','RK4','RKF45','Rosenbrock4','Rosen23mod', 'ODE23s');
	set(lh,'FontSize',12);
	xlabel('t','FontSize',12);
	ylabel('Species 1','FontSize',12);
	title('Integration of Van der Pol oscillator', 'FontSize', 12);
	box on;
	grid on;
	axis square;

	% print out various statistics
	
	fprintf('\nVan der Pol oscillator\n\n');

	fprintf('Times\n-----\n');
	fprintf('Explicit Euler  %0.10f\n',expeuler_time);
	fprintf('RK4             %0.10f\n',rk4_time);
	fprintf('RKF45           %0.10f\n',rkf45_time);
	fprintf('Rosenbrock      %0.10f\n',rosenbrock_time);
	fprintf('Rosen23mod      %0.10f\n',rosen23mod_time);
	fprintf('Ode23s          %0.10f\n',ode23s_time);
	fprintf('\n');

	fprintf('Steps\n-----\n');
	fprintf('Explicit Euler   %d\n',rk4_steps);
	fprintf('RK4              %d\n',rk4_steps);
	fprintf('RKF45            %d\n',rkf45_steps);
	fprintf('Rosenbrock       %d\n',rosenbrock_steps);
	fprintf('Rosen23mod       %d\n',rosen23mod_steps);
	fprintf('Ode23s           %d\n',ode23s_steps);
	fprintf('\n');


	% set up integration of Gear's system

	T = 200;
	tspan = [0 T];
	h = 1e-3;
	tol = 1e-4;
	in_con = [1 0];

	% integrate with various methods and time
	
	tic
	[t Y] = expeuler(@stiffr,tspan,in_con,h);
	expeuler_time = toc;
	expeuler_steps = length(t);

	tic
	[t1 Y1] = rk4(@stiffr,tspan,in_con,h);
	rk4_time = toc;
	rk4_steps = length(t1);

	tic
	[t2 Y2] = rkf45(@stiffr,tspan,in_con,h,tol);
	rkf45_time = toc;
	rkf45_steps = length(t2);

	tic
	[t3 Y3] = rosenbrock(@stiffr,@stiffjac,tspan,in_con,h,tol);
	rosenbrock_time = toc;
	rosenbrock_steps = length(t3);

	tic
	[t4 Y4] = rosen23mod(@stiffr,@stiffjac,tspan,in_con,h,tol);
	rosen23mod_time = toc;
	rosen23mod_steps = length(t4);
	
	options = odeset('Abstol',1e-13,'Reltol',tol);
	tic
	[t5 Y5] = ode23s(@stiffc,tspan,in_con,options);
	ode23s_time = toc;
	ode23s_steps = length(t5);

	% get exact solution

	tt = linspace(0,T/100,1000);
	yy = stiffex(tt);

	% plot the approximations for the first species against the exact solution

	figure
	hold all
	plot(tt,yy(1,:));
	plot(t,Y(:,1),t1,Y1(:,1),'-',t2,Y2(:,1),'.-',t3,Y3(:,1),'o-',t4,Y4(:,1),'*-', t5,Y5(:,1));
	hold off
	lh = legend('Exact solution','Explicit Euler','RK4','RKF45','Rosenbrock4','Rosen23mod', 'ODE23s');
	set(lh,'FontSize',12);
	xlabel('t','FontSize',12);
	ylabel('Species 1','FontSize',12);
	title('Integration of Gear''s system', 'FontSize', 12);
	axis([-0.1 2 ylim+[-0.1 0.1]]);
	box on;
	grid on;
	axis square;

	% print out various statistics
	
	fprintf('\nGear''s system\n\n')

	fprintf('Times\n-----\n');
	fprintf('Explicit Euler  %0.10f\n',expeuler_time);
	fprintf('RK4             %0.10f\n',rk4_time);
	fprintf('RKF45           %0.10f\n',rkf45_time);
	fprintf('Rosenbrock      %0.10f\n',rosenbrock_time);
	fprintf('Rosen23mod      %0.10f\n',rosen23mod_time);
	fprintf('Ode23s          %0.10f\n',ode23s_time);
	fprintf('\n');

	fprintf('Steps\n-----\n');
	fprintf('Explicit Euler   %d\n',rk4_steps);
	fprintf('RK4              %d\n',rk4_steps);
	fprintf('RKF45            %d\n',rkf45_steps);
	fprintf('Rosenbrock       %d\n',rosenbrock_steps);
	fprintf('Rosen23mod       %d\n',rosen23mod_steps);
	fprintf('Ode23s           %d\n',ode23s_steps);
	fprintf('\n');

	% get E_2 norm measurements for each method

	Error = zeros(6,1);
	Yexact = stiffex(t);
	Error(1) = norm(Y(:,1) - Yexact(1,:)') ./ norm(Yexact);
	Yexact = stiffex(t1);
	Error(2) = norm(Y1(:,1) - Yexact(1,:)') ./ norm(Yexact);
	Yexact = stiffex(t2);
	Error(3) = norm(Y2(:,1) - Yexact(1,:)') ./ norm(Yexact);
	Yexact = stiffex(t3);
	Error(4) = norm(Y3(:,1) - Yexact(1,:)') ./ norm(Yexact);
	Yexact = stiffex(t4);
	Error(5) = norm(Y4(:,1) - Yexact(1,:)') ./ norm(Yexact);
	Yexact = stiffex(t5);
	Error(6) = norm(Y5(:,1) - Yexact(1,:)') ./ norm(Yexact);

	% print out error statistics

	fprintf('Relative errors\n---------------\n');
	fprintf('Explicit Euler  %0.10f\n',Error(1));
	fprintf('RK4             %0.10f\n',Error(2));
	fprintf('RKF45           %0.10f\n',Error(3));
	fprintf('Rosenbrock      %0.10f\n',Error(4));
	fprintf('Rosen23mod      %0.10f\n',Error(5));
	fprintf('Ode23s          %0.10f\n',Error(6));
	fprintf('\n');

end

% Four functions comprise Gear's equation here
% the row one is for my implementations
% the column one is needed for ode23s
% the Jacobian is needed for my Rosenbrock method implementations
% the last on is the exact solution

function dy = stiffr(t,y)
	dy = zeros(1,2);
	dy(1) = 998*y(1) + 1998*y(2);
	dy(2) = -999*y(1) - 1999*y(2);
end

function dy = stiffc(t,y)
	dy = zeros(2,1);
	dy(1) = 998*y(1) + 1998*y(2);
	dy(2) = -999*y(1) - 1999*y(2);
end

function jac = stiffjac(t,y)
	jac = [998 1998;
			-999 -1999];
end

function y = stiffex(t)
	N = length(t);
	y = zeros(2,N);
	y(1,1:N) = 2*exp(-t) - exp(-1000*t);
	y(2,1:N) = -exp(-t) + exp(-1000*t);
end

% Three functions comprise Van der Pol's equation
% the row one is for my implementations
% the column one is needed for ode23s
% the Jacobian is needed for my Rosenbrock method implementations

function dy = vdpr(t,y)
	mu = 500;
	dy = zeros(1,2);
	dy(1) = mu*(y(1) - y(1)^3 / 3 - y(2));
	dy(2) = y(1) / mu;
end

function dy = vdpc(t,y)
	mu = 500;
	dy = zeros(2,1);
	dy(1) = mu*(y(1) - y(1)^3 / 3 - y(2));
	dy(2) = y(1) / mu;
end

function jac = vdpjac(t,y)
	mu = 500;
	jac = [mu*(1-y(1)^2) 	-mu;...
			1/mu 			0];
end






