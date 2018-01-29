function ttcvode_test(x0, tfinal, T, a, b)
tic;
[tout, yout] = ttcvode({[T, a*T, b*T]}, @func, [0, tfinal], x0);
toc

plot(yout(1, :), yout(2, :));

    function [dy, flag] = func(t,y)
        flag = 0;
        tt = mod(t, T)/T;
        if tt >= a && tt < b
            dy = [-0.2, 0; 0.1, -0.3] * y;
        else
            dy = [-0.1, -0.05; 0, 0.1] * y;
        end
    end
end
