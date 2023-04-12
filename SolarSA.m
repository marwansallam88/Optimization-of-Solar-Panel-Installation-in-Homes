%cleaning up the working environment before starting
clc, clear, close all
tic

%Initiialize Panel Databank
Panels.DB = xlsread('optimizationdatabase','F3:I52');%(x,1)Power %(x,2)Area %(x,3)TempCof %(x,4)Cost

%Initiialize input parameters
%Assume all panels are square
amb_Temp = 35;
area_Ava = 100;
budget = 20000;
power_Req = 4500;
allowable_Loss_Percent = 15;
allowable_Waste_Percent = 0.1;
desiredNum = 5;

%Initialize the optimization problem parameters:
Panels.num = 50; %number of panel types
SA.Max_num_iter = 1000;  %maximum number of iterations
SA.T_init = 500;  %initial temperature used in SA
SA.T_final = 0.1;  %define the final temperature to be reached
SA.T_current = SA.T_init;  %save initial temperature as current temperature
SA.Beta = (SA.T_init - SA.T_final)/SA.Max_num_iter;  %calculate the linear cooling parameter (beta)
SA.alpha = 0.99;  %for geometric cooling schedule
SA.Linear_Cooling = 0; %Flag to select linear cooling, if (0) then geometric cooling
SA.Feasible = 0; %Flag to decide whether to check feasibility
SA.Chk_Probability = 1;  %Flag to decide whether to check probability or to bypass it

while(SA.Feasible == 0) %loop until sol is feasible
    %Initialize a random solution:
    SA.Sol_init = zeros(1,Panels.num); %initialize an initial arithmatic sol empty array
    SA.Sol_initB = zeros(1,Panels.num); %initialize an initial binary sol empty array
    cell = randi([1,10]);
    SA.Sol_initB(1,cell) = 1;
    cell = randi([11,20]);
    SA.Sol_initB(1,cell) = 1;
    cell = randi([21,30]);
    SA.Sol_initB(1,cell) = 1;
    cell = randi([31,40]);
    SA.Sol_initB(1,cell) = 1;
    cell = randi([41,50]);
    SA.Sol_initB(1,cell) = 1;

    for i = 1:size(SA.Sol_init,2)  %for all the elements in the array
        val = randi(50);
        SA.Sol_init(1,i) = val;  %update value of element in solution array
    end

    SA.Sol_current = zeros(1,Panels.num); %initialize a current sol empty array
    SA.Sol_currentB = zeros(1,Panels.num);
    SA.Sol_current = SA.Sol_init;  %save initial solution as the current solution
    SA.Sol_currentB = SA.Sol_initB;
    clear i val;  %clear used temporary variables

    %Calculate Cost of initial Solution:
    SA.Sol_cost = 0;
    total_Cost = 0;
    total_Power = 0;
    total_Cof= 0;
    total_Num = 0;
    total_Area = 0;

    for i = 1:size(SA.Sol_current,2)  %for all the elements in the array
        total_Cost = total_Cost + Panels.DB(i,4)*SA.Sol_current(1,i)*SA.Sol_currentB(1,i); %total cost calculation
        total_Power = total_Power + Panels.DB(i,1)*SA.Sol_current(1,i)*SA.Sol_currentB(1,i)*0.9; %total power calculation
        total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(i,3)*SA.Sol_currentB(1,i)*SA.Sol_current(1,i); %total temp coeffiecient calculation
        total_Num = total_Num + SA.Sol_current(1,i)*SA.Sol_currentB(1,i); %total number of panels in current solution calculation
        total_Area = total_Area + SA.Sol_current(1,i)*Panels.DB(i,2)*SA.Sol_currentB(1,i); %total area calculation
    end

    area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
    loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
    if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
        SA.Feasible = 1; %if feasible set flag to 1
        Sol_cost = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2; %calculating cost of initial sol
    end
end
SA.Sol_cost = Sol_cost;
Best_Sol = Sol_cost;
Best_SolA = SA.Sol_current;
Best_SolB = SA.Sol_currentB;

%Plot the cost function figure
Figures.Main_fig = figure;  %create new figure
SA.Cost_fcn_curve = zeros(1,SA.Max_num_iter+1);
SA.Cost_fcn_curve(1,1) = SA.Sol_cost;  %save the cost of initial solution
figure(Figures.Main_fig);  %access the main figure
subplot(1,3,1);  %in the 1st location in the figure
Figures.Cost = plot(0,SA.Cost_fcn_curve(1,1),'m*','LineWidth',1.6); %plot the initial cost
xlim([0,SA.Max_num_iter]); %define limits in X-direction
grid on; %initialize the grid to be on
axis square; %make the axes look like square
xlabel('Number of Iterations');
ylabel('Cost function');
title('The Cost function');  %give title to figure
legend('Cost function');

%Plot the Temperature function figure
SA.Temp_fcn_curve = zeros(1,SA.Max_num_iter+1);
SA.Temp_fcn_curve(1,1) = SA.T_init;  %save the temp of initial solution
figure(Figures.Main_fig);  %access the main figure
subplot(1,3,2);  %in the 2nd location in the figure
Figures.Temp = plot(0,SA.Temp_fcn_curve(1,1),'b*','LineWidth',1.6); %plot the initial temperature
xlim([0,SA.Max_num_iter]); %define limits in X-direction
ylim([0,SA.T_init]); %define limits in Y-direction
grid on; %initialize the grid to be on
axis square; %make the axes look like square
xlabel('Number of Iterations');
ylabel('Temp function (Celcius)');
title('The Temperature function'); %give title to figure
legend('Temp function (C)');

clear total_Area total_Num total_Cof total_Power total_Cost area_Waste loss i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the SA optimization Loop here   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:SA.Max_num_iter  %for the maximum number of iterations
    if iter == 1 %if this is the first iteration
        Sol_old = SA.Sol_init;  %read the initial solution
        Sol_oldB = SA.Sol_initB;
        Sol_old_cost = SA.Sol_cost;  %read the initial solution cost
    end

    %Randomize a new solution:
    %==========================
    SA.Feasible = 0;
    while(SA.Feasible==0)
        Scaler = 0.5; %number [0-1] to define scale of exploration around current solution
        Sol_current = zeros(1,Panels.num);  %initialize an empty array
        Sol_currentB = Sol_oldB;
        counter = 0;
        for i = 1:size(Sol_current,2)  %for all the elements in the array
            val = Sol_old(1,i) + round(randi([-1,1]) * Scaler * rand * 10);
            if(counter<1)
                ranCell2 = randi(Panels.num);
                saveB1 = Sol_currentB(1,ranCell2);
                ranCell =  randi(Panels.num);
                saveB2 = Sol_currentB(1,ranCell);
                Sol_currentB(1,ranCell2) = saveB2;
                Sol_currentB(1,ranCell) = saveB1;
                counter = counter +1;
            end
            if val<0
                val = mod(abs(val),50);
            end
            if val>50
                val=mod(val,50);
            end
            Sol_current(1,i) = val;
        end
        clear i val counter;  %clear used temporary variables

        %Calculate cost of new solution:
        %===================================
        Sol_cost = 0;
        total_Cost = 0;
        total_Power = 0;
        total_Cof= 0;
        total_Num = 0;
        total_Area = 0;
        for i = 1:Panels.num  %for all the elements in the array
            total_Cost = total_Cost + Panels.DB(i,4)*Sol_current(1,i)*Sol_currentB(1,i); %total cost calculation
            total_Power = total_Power + Panels.DB(i,1)*Sol_current(1,i)*Sol_currentB(1,i)*0.9; %total power calculation
            total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(i,3)*Sol_currentB(1,i)*Sol_current(1,i); %total temp coeffiecient calculation
            total_Num = total_Num + Sol_current(1,i)*Sol_currentB(1,i); %total number of panels in current solution calculation
            total_Area = total_Area + Sol_current(1,i)*Panels.DB(i,2)*Sol_currentB(1,i); %total area calculation
        end
        area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
        loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
        if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
            SA.Feasible = 1; %if feasible set flag to 1
            Sol_cost = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2; %calculating cost of sol
        end
    end

    %Check if new solution is better than old one:
    %===============================================
    cost_diff = Sol_cost - Sol_old_cost;  %difference between two costs
    if cost_diff < 0  %if new solution is better
        better_sol_flag = 1;  %raise flag
    else
        better_sol_flag = 0;  %dont raise flag
    end

    %Check what to do in each case:
    %=================================
    if better_sol_flag == 1  %if better solution is evaluated
        SA.Sol_current = Sol_current;  %save the current solution
        SA.Sol_currentB = Sol_currentB;
        SA.Sol_cost = Sol_cost;  %save the cost of current solution

    elseif better_sol_flag == 0  %if worse solution is evaluated
        if SA.Chk_Probability == 1
            T_current = SA.T_current;  %read the current temperature
            P = exp((-1*abs(cost_diff))/T_current);  %calculate the probability to accept
        else
            P = 0;  %OVERWRITE Prob part
        end
        rand_num = rand;  %generate a random number [0-1]
        %Accept or reject worse solution based on probability
        if P > rand_num  %we accept the worse solution
            SA.Sol_current = Sol_current;  %save the current solution
            SA.Sol_currentB = Sol_currentB;
            SA.Sol_cost = Sol_cost;  %save the cost of current solution
        else  %we reject the worse solution
            SA.Sol_current = SA.Sol_current;  %save the old solution
            SA.Sol_currentB = SA.Sol_currentB;  %save the old solution
            SA.Sol_cost = SA.Sol_cost;  %save the cost of old solution
        end
    end

    SA.Cost_fcn_curve(1,iter+1) = SA.Sol_cost;  %save the cost of current solution
    set(Figures.Cost,'XData',0:1:iter,'YData',SA.Cost_fcn_curve(1,1:iter+1));
    subplot(1,3,1); title(['Current Cost Value: ' num2str(SA.Cost_fcn_curve(1,iter+1))]);
    SA.Temp_fcn_curve(1,iter+1) = SA.T_current;  %save the current temperature of solution
    set(Figures.Temp,'XData',0:1:iter,'YData',SA.Temp_fcn_curve(1,1:iter+1));
    pause(0.05);

    %Update the current temperature:
    %=================================
    if SA.Linear_Cooling == 1  %for linear cooling schedule
        SA.T_current = SA.T_init - (SA.Beta * iter);  %update the temperature based on the linear cooling schedule
    else  %then we'll use geometric cooling
        SA.T_current = SA.T_init * (SA.alpha ^ iter);  %update the temperature based on the geometric cooling schedule
    end

    %Update variables for next iteration:
    Sol_old = SA.Sol_current;  %update the current solution
    Sol_oldB = SA.Sol_currentB;
    Sol_old_cost = SA.Sol_cost;  %update the current cost
    if Best_Sol > SA.Sol_cost
        Best_Sol = SA.Sol_cost;
        Best_SolA = SA.Sol_current;
        Best_SolB = SA.Sol_currentB;
        Best_cost = total_Cost;
        Best_power = total_Power;
        Best_area = total_Area;
        Best_area_waste = area_Waste;
        Best_loss = loss;
    end

    [~,Panels1.DB] = xlsread('optimizationdatabase2.xlsx','B3:B52'); %Panel Name and manufacturer
    PanelCount = zeros(1,desiredNum);
    Paneli = zeros(1,desiredNum);
    counter2 = 0;
    for b = 1:Panels.num
        if Best_SolB(1,b)==1
            counter2 = counter2 + 1;
            PanelCount(1,counter2) = Best_SolA(1,b);
            Paneli(1,counter2) = b;
        end
    end
    
    %Plotting pie chart
    labels = {num2str(PanelCount(1,1)), num2str(PanelCount(1,2)), num2str(PanelCount(1,3)), num2str(PanelCount(1,4)), num2str(PanelCount(1,5))};
    subplot(1,3,3);
    pie(PanelCount, labels);
    labels1 = {string(Panels1.DB(Paneli(1),1)),string(Panels1.DB(Paneli(2),1)),string(Panels1.DB(Paneli(3),1)),string(Panels1.DB(Paneli(4),1)),string(Panels1.DB(Paneli(5),1))};
    legend(labels1,'Location','southoutside','Orientation','vertical');
end

disp(['Best cost is:    ' num2str(Best_cost)]);
disp(['Best power is:   ' num2str(Best_power)]);
disp(['Best area is:   ' num2str(Best_area)]);
disp(['Best area wasted is:   ' num2str(Best_area_waste)]);
disp(['Best power loss is:   ' num2str(Best_loss)]);
disp(['Best cost function is: ' num2str(Best_Sol)]);
timElapsed = toc;
disp(['Time of execution: ' num2str(timElapsed)]);



