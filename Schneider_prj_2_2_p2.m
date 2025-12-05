%% Input data
    x = [1.0; 2.0; 3.0; 1.5; 2.5; 3.5; 2.25];
    y = [22.14; 21.70; 22.08; 21.96; 21.79; 22.27; 21.68];

%% Vandermonde matrix 
    A = [ones(size(x)), x, x.^2];

    %disp(A)

    Beta = pinv(A) * y;

    %disp (Beta);
   

%% regression fit math
    y_hat = Beta(1) + Beta(2) * x + Beta(3) * x.^2;
    y_avg = mean(y);
    n = length(y);
    p = 2;
    %%=========
    ss_res = sum((y - y_hat).^ 2);
    ss_t = sum((y - y_avg) .^2);
    %%=========
    R_sqrd = 1 - (ss_res / ss_t);
    R_sqrd_adj = 1 - ((1-R_sqrd) * (n-1) / (n - p - 1));

%% Regression Fit model
    
    fprintf("\n\tRegression fit Model.\n");
    fprintf("\tR-squared: %.4f\n", R_sqrd);
    fprintf("\tAdjusted R-squared: %.4f\n", R_sqrd_adj);

%% Plot math
    xplt = (0.0:0.01:5.0);
    % disp(xplt);
%%Regression Model yplt == y_hat of x plot 
    yplt = Beta(1) + Beta(2) * xplt + Beta(3) * xplt.^2;
   

%%optimization

    x_opt = (-Beta(2) / (2 * Beta(3)));
    y_min = Beta(1) + Beta(2) * x_opt + Beta(3) * x_opt ^ 2;
    
%% optimized print statments
    fprintf("\n\tRear Stagger settings and lap time.\n");
    fprintf("\tRear Stagger (RS optimized): %.3f inches \t \n", x_opt);
    fprintf("\tLap time for optimized RS %.3f seconds \n", y_min);

%%PLOTS
    
    plot(xplt,yplt,"b-");
    hold on;
    % bounding x and y for 3 data points
    xredplt = x([2, 4, 7]); % limited data bounds for x
    yredplt =y([2, 4, 7]); % limited data bounds for y
    %plot(x, y, 'r*');
    plot(xredplt, yredplt, 'r*');
    plot(x_opt, y_min, "go");
    xlabel("Rear Stagger (inches)");
    ylabel("Lap Time (seconds)");
    legend(["Rear Stagger Vs Lap time ", "Known Data", " Optimal Rear Stagger and Time"])
    grid on
    hold off;