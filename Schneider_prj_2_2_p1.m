
%% input and orginal data
    x = [1.0; 2.0; 3.0];
    y = [22.14; 21.70; 22.08];

%% Vandermonde matrix 
    A = [ones(size(x)), x, x.^2];

    %disp(A);

    Beta = pinv(A) * y;

    % {disp(Beta(1));
    % disp(Beta(2));
    % disp(Beta(3));}

    xplt = (0.0:0.01:5.0);

    % disp(xplt);
    
    %%Regression Model yplt == y_hat of x plot 
    yplt = Beta(1) + Beta(2) * xplt + Beta(3) * xplt.^2;


    %%optmiaztion 

    x_opt = (-Beta(2) / (2 * Beta(3)));
    y_min = Beta(1) + Beta(2) * x_opt + Beta(3) * x_opt ^ 2;
    
    %% First print statments
    fprintf("\n\tRear Stagger settings and lap time.\n");
    fprintf("\tRear Stagger (RS optimized): %.3f inches \t \n", x_opt);
    fprintf("\tLap time for optimized RS %.3f seconds \n", y_min);

    %%PLOTS
    
    plot(xplt,yplt,"b-");
    hold on;
    plot(x, y, 'r*');
    plot(x_opt, y_min, "go");
    xlabel("Rear Stagger (inches)");
    ylabel("Lap Time (seconds)");
    legend(["Rear Stagger Vs Lap time ", "Known Data", " Optimal Rear Stagger and Time"])
    grid on
    hold off;
    