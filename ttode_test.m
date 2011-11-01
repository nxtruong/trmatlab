function ttode_test(x0, tfinal, T, a, b)
%[tout, yout, te, ie, ye] = ode45(@func, [0, tfinal], x0, odeset('Events', @events,'MaxStep',T/50));
%[tout, yout, te, ie, ye] = ode45(@func, [0, tfinal], x0, odeset('Events',
%@events));
[tout, yout] = ttode(@ode45, struct('ttevents', {{[T, a*T, b*T]}}),...
    @func, [0, tfinal], x0);
plot(yout(:,1), yout(:,2));

    function dy = func(t,y)
        tt = mod(t, T)/T;
        if tt >= a && tt < b
            dy = [-0.2, 0; 0.1, -0.3] * y;
        else
            dy = [-0.1, -0.05; 0, 0.1] * y;
        end
    end
end
