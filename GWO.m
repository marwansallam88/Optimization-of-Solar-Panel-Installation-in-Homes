%cleaning up the working environment before starting
clc, clear, close all

tic
%Initiialize Panel Databank
Panels.DB = xlsread('optimizationdatabase','F3:I52');%(x,1)Power %(x,2)Area %(x,3)TempCof %(x,4)Cost

%Initiialize input parameters
%Assume all panels are square
amb_Temp = 35;
area_Ava = 70;
budget = 20000;
power_Req = 4500;
allowable_Loss_Percent = 15;
allowable_Waste_Percent = 0.1;
desiredNum = 5;

%Initialize the optimization problem parameters:
Panels.num = 50; %number of panel types
GW.Num_Generations = 1000;  %maximum number of generations to be tested
GW.Pop_size = 50;  %the population size
GW.Feasible = 0;
stopping = 0;

GW.Pop_init = zeros(GW.Pop_size,Panels.num);  %initialize an empty matrix
GW.Pop_initB = zeros(GW.Pop_size,Panels.num);
GW.Pop_cost = zeros(GW.Pop_size,1);  %initialize an empty array of zeros

for i = 1:GW.Pop_size  %for all the population members
    GW.Feasible = 0;
    while(GW.Feasible == 0) %loop until sol is feasible
        GW.Pop_init(i,:) = zeros(1,Panels.num);  %initialize an empty matrix
        GW.Pop_initB(i,:) = zeros(1,Panels.num);
        GW.Pop_cost(i,1) = 0;
        cell1 = randi([1,10]);
        GW.Pop_initB(i,cell1) = 1;
        cell2 = randi([11,20]);
        GW.Pop_initB(i,cell2) = 1;
        cell3 = randi([21,30]);
        GW.Pop_initB(i,cell3) = 1;
        cell4 = randi([31,40]);
        GW.Pop_initB(i,cell4) = 1;
        cell5 = randi([41,50]);
        GW.Pop_initB(i,cell5) = 1;
        for j = 1:Panels.num  %for all the elements in the array
            val = randi(50);
            GW.Pop_init(i,j) = val;  %update value of element in solution array
        end

        %Calculate Cost of initial Solution:
        total_Cost = 0;
        total_Power = 0;
        total_Cof= 0;
        total_Num = 0;
        total_Area = 0;
        for j = 1:Panels.num  %for all the elements in the array
            total_Cost = total_Cost + Panels.DB(j,4)*GW.Pop_init(i,j)*GW.Pop_initB(i,j); %total cost calculation
            total_Power = total_Power + Panels.DB(j,1)*GW.Pop_init(i,j)*GW.Pop_initB(i,j)*0.9; %total energy calculation
            total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*GW.Pop_init(i,j)*GW.Pop_initB(i,j); %total temp coeffiecient calculation
            total_Num = total_Num + GW.Pop_init(i,j)*GW.Pop_initB(i,j); %total number of panels in current solution calculation
            total_Area = total_Area + GW.Pop_init(i,j)*Panels.DB(j,2)*GW.Pop_initB(i,j); %Total area calculation
        end
        area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
        loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
        if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
            GW.Feasible = 1; %if feasible set flag to 1
            %GW.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste; %calculating cost of initial sol
            GW.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2; %calculating cost of initial sol
        end
    end
end

GW.Alpha = zeros(1,Panels.num);
GW.AlphaB = zeros(1,Panels.num);
GW.Beta = zeros(1,Panels.num);
GW.BetaB = zeros(1,Panels.num);
GW.Delta = zeros(1,Panels.num);
GW.DeltaB = zeros(1,Panels.num);
GW.Pop_current = GW.Pop_init;  %save initial population as the current population
GW.Pop_currentB = GW.Pop_initB;
clear i j k m cell val total_Area total_Num total_Cof total_Power total_Cost;  %clear used temporary variables

%Plot the cost function figure
Figures.Main_fig = figure;  %create new figure
Best_Cost_Array = zeros(1,GW.Num_Generations+1);
Best_Cost_Array(1,1) = min(GW.Pop_cost);  %save the cost of initial solution
figure(Figures.Main_fig);  %access the main figure
subplot(1,2,1);  %in the 1st location in the figure
Figures.Cost = plot(0,Best_Cost_Array(1,1),'m*','LineWidth',1.6); %plot the initial cost
xlim([0,GW.Num_Generations]); %define limits in X-direction
grid on; %initialize the grid to be on
axis square; %make the axes look like square
xlabel('Number of Generations');
ylabel('Cost function');
title('The Cost function');  %give title to figure
legend('Cost function');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the grey wolf optimization Loop here   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for gen = 1:GW.Num_Generations  %for the maximum number of generations
    if stopping == 300
        break;
    end

    %Update a,A,C
    GW.a = 2 - (2/GW.Num_Generations)*gen;

    %Initialize a new generation matrix to be filled by the newly created members:
    GW_new = zeros(GW.Pop_size,Panels.num);  %initialize an empty matrix
    GW_newB = zeros(GW.Pop_size,Panels.num);  %initialize an empty matrix
    GW_new_cost = zeros(GW.Pop_size,1);  %initialize an empty array of zeros
    Best_SolA = zeros(1,Panels.num);
    Best_SolB = zeros(1,Panels.num);

    [~, member_index] = sort(GW.Pop_cost);
    GW.Alpha = GW.Pop_current(member_index(1),:);
    GW.AlphaB = GW.Pop_currentB(member_index(1),:);
    GW.Beta = GW.Pop_current(member_index(2),:);
    GW.BetaB = GW.Pop_currentB(member_index(2),:);
    GW.Delta = GW.Pop_current(member_index(3),:);
    GW.DeltaB = GW.Pop_currentB(member_index(3),:);
    D_alpha = 0;
    D_beta = 0;
    D_delta = 0;
    x1 = 0;
    x2 = 0;
    x3 = 0;
    xAvg = zeros(1,Panels.num);
    omegaB = zeros(1,Panels.num);

    GW_new(1,:) = GW.Alpha;
    GW_new(2,:) = GW.Beta;
    GW_new(3,:) = GW.Delta;
    GW_newB(1,:) = GW.AlphaB;
    GW_newB(2,:) = GW.BetaB;
    GW_newB(3,:) = GW.DeltaB;
    GW_new_cost(1,1) = GW.Pop_cost(member_index(1),1);
    GW_new_cost(2,1) = GW.Pop_cost(member_index(2),1);
    GW_new_cost(3,1) = GW.Pop_cost(member_index(3),1);

    for i = 4:GW.Pop_size  %for the required number of omegas
        GW.Feasible = 0;
        while(GW.Feasible == 0)
            GW.A = 2*rand()*GW.a - GW.a;
            GW.C = 2*rand();
            GW.Feasible = 0;
            if(abs(GW.A) < 1)
                for j = 1:Panels.num
                    D_alpha = abs((GW.C*GW.Alpha(1,j)) - GW.Pop_current(member_index(i),j));
                    D_beta = abs((GW.C*GW.Beta(1,j)) - GW.Pop_current(member_index(i),j));
                    D_delta = abs((GW.C*GW.Delta(1,j)) - GW.Pop_current(member_index(i),j));
                    x1 = abs(GW.Alpha(1,j) - (abs(GW.A)*D_alpha));
                    if(x1 > 50)
                        x1 = mod(x1,50);
                    end
                    x2 = abs(GW.Beta(1,j) - (abs(GW.A)*D_beta));
                    if(x2 > 50)
                        x2 = mod(x2,50);
                    end
                    x3 = abs(GW.Delta(1,j) - (abs(GW.A)*D_delta));
                    if(x3 > 50)
                        x3 = mod(x3,50) ;
                    end
                    xAvg(1,j) = round((x1 + x2 + x3)/3);
                end
            elseif(abs(GW.A) >= 1)
                for j = 1:Panels.num
                    D_alpha = abs((GW.C*GW.Alpha(1,j)) - GW.Pop_current(member_index(i),j));
                    D_beta = abs((GW.C*GW.Beta(1,j)) - GW.Pop_current(member_index(i),j));
                    D_delta = abs((GW.C*GW.Delta(1,j)) - GW.Pop_current(member_index(i),j));
                    x1 = GW.Alpha(1,j) + (abs(GW.A)*D_alpha);
                    if(x1 > 50)
                        x1 = mod(x1,50);
                    end
                    x2 = GW.Beta(1,j) + (abs(GW.A)*D_beta);
                    if(x2 > 50)
                        x2 = mod(x2,50);
                    end
                    x3 = GW.Delta(1,j) + (abs(GW.A)*D_delta);
                    if(x3 > 50)
                        x3 = mod(x3,50) ;
                    end
                    xAvg(1,j) = round((x1 + x2 + x3)/3);
                end
            end
            omegaB = GW.Pop_currentB(member_index(i),:);
            for counter = 1:10
                ranCell2 = 0;
                ranCell = 0;
                while ranCell2 == ranCell
                    ranCell2 = randi([1,Panels.num]);
                    ranCell =  randi([1,Panels.num]);
                end
                saveB1 = omegaB(1,ranCell2);
                saveB2 = omegaB(1,ranCell);
                omegaB(1,ranCell2) = saveB2;
                omegaB(1,ranCell) = saveB1;
            end
            GW_newB(i,:) = omegaB;
            total_Cost = 0;
            total_Power = 0;
            total_Cof= 0;
            total_Num = 0;
            total_Area = 0;
            for j = 1:Panels.num  %for all the elements in the
                total_Cost = total_Cost + Panels.DB(j,4)*xAvg(1,j)*GW_newB(i,j); %total cost calculation
                total_Power = total_Power + Panels.DB(j,1)*xAvg(1,j)*GW_newB(i,j)*0.9; %total energy calculation
                total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*xAvg(1,j)*GW_newB(i,j); %total temp coeffiecient calculation
                total_Num = total_Num + xAvg(1,j)*GW_newB(i,j); %total number of panels in current solution calculation
                total_Area = total_Area + xAvg(1,j)*Panels.DB(j,2)*GW_newB(i,j); %Total area calculation
            end
            area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
            loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
            if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
                GW.Feasible = 1;
                GW_new(i,:)= xAvg;
                %GW_new_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste;
                GW_new_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
            end
        end
    end
    GW.Pop_current = GW_new;
    GW.Pop_currentB = GW_newB;
    GW.Pop_cost = GW_new_cost;
    [~, member_index] = sort(GW_new_cost);
    Best_SolB(1,:) = GW_newB(member_index(1),:);
    Best_SolA(1,:) = GW_new(member_index(1),:);
    bestSol = min(GW_new_cost);
    Best_Cost_Array(1,gen+1) = bestSol;
    if Best_Cost_Array(1,gen) == bestSol
        stopping = stopping + 1;
    else
        stopping = 0;
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

    %Plotting function curve
    set(Figures.Cost,'XData',0:1:gen,'YData',Best_Cost_Array(1,1:gen+1));
    subplot(1,2,1); title(['Current Cost Value: ' num2str(Best_Cost_Array(1,gen+1))]);
    pause(0.05);

    %Plotting pie chart
    labels = {num2str(PanelCount(1,1)), num2str(PanelCount(1,2)), num2str(PanelCount(1,3)), num2str(PanelCount(1,4)), num2str(PanelCount(1,5))};
    subplot(1,2,2);
    pie(PanelCount, labels);
    labels1 = {string(Panels1.DB(Paneli(1),1)),string(Panels1.DB(Paneli(2),1)),string(Panels1.DB(Paneli(3),1)),string(Panels1.DB(Paneli(4),1)),string(Panels1.DB(Paneli(5),1))};
    legend(labels1,'Location','southoutside','Orientation','vertical');
end

[~, member_index1] = sort(GW_new_cost);
finalSol = zeros(1,Panels.num);
for j = 1:Panels.num
    finalSol(1,j) = GW_new(member_index1(1),j)*GW_newB(member_index1(1),j);
end
total_Cost = 0;
total_Power = 0;
total_Cof= 0;
total_Num = 0;
total_Area = 0;
for j = 1:Panels.num
    total_Cost = total_Cost + Panels.DB(j,4)*finalSol(1,j); %total cost calculation
    total_Power = total_Power + Panels.DB(j,1)*finalSol(1,j)*0.9; %total energy calculation
    total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*finalSol(1,j); %total temp coeffiecient calculation
    total_Num = total_Num + finalSol(1,j); %total number of panels in current solution calculation
    total_Area = total_Area + finalSol(1,j)*Panels.DB(j,2); %Total area calculation
end
area_Waste = (area_Ava-total_Area)/area_Ava;
loss = -(total_Cof/total_Num);
disp(['Best cost is:    ' num2str(total_Cost)]);
disp(['Best power is:   ' num2str(total_Power)]);
disp(['Best area is:   ' num2str(total_Area)]);
disp(['Best area wasted is:   ' num2str(area_Waste)]);
disp(['Best power loss is:   ' num2str(loss)]);
disp(['Best cost function is: ' num2str(bestSol)]);
timElapsed = toc;
disp(['Time of execution: ' num2str(timElapsed)]);