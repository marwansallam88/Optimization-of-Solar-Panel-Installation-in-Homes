%cleaning up the working environment before starting
clc, clear, close all

tic
%Initiialize Panel Databank
Panels.DB = xlsread('optimizationdatabase','F3:I52');

%Initiialize input parameters
%Assume all panels are square
amb_Temp = 35;
area_Ava = 70;
budget = 50000;
power_Req = 3000;
allowable_Loss_Percent = 15;
allowable_Waste_Percent = 0.5;
desiredNum = 5;

%Initialize the optimization problem parameters:
Panels.num = 50; %number of panel types
MVO.Num_Generations = 1000;  %maximum number of generations to be tested
MVO.Pop_size = 10;  %the population size
MVO.Feasible = 0;
MVO.original_pos = zeros(1,Panels.num);
MVO.WEP_Max = 1;
MVO.WEP_Min = 0.1;
stopping = 0;

MVO.Pop_init = zeros(MVO.Pop_size,Panels.num);  %initialize an empty matrix
MVO.Pop_initB = zeros(MVO.Pop_size,Panels.num);
MVO.Pop_cost = zeros(MVO.Pop_size,1);  %initialize an empty array of zeros

for i = 1:MVO.Pop_size  %for all the population members
    MVO.Feasible = 0;
    while(MVO.Feasible == 0) %loop until sol is feasible
        MVO.Pop_init(i,:) = zeros(1,Panels.num);  %initialize an empty matrix
        MVO.Pop_initB(i,:) = zeros(1,Panels.num);
        MVO.Pop_cost(i,1) = 0;
        cell1 = randi([1,10]);
        MVO.Pop_initB(i,cell1) = 1;
        cell2 = randi([11,20]);
        MVO.Pop_initB(i,cell2) = 1;
        cell3 = randi([21,30]);
        MVO.Pop_initB(i,cell3) = 1;
        cell4 = randi([31,40]);
        MVO.Pop_initB(i,cell4) = 1;
        cell5 = randi([41,50]);
        MVO.Pop_initB(i,cell5) = 1;
        for j = 1:Panels.num  %for all the elements in the array
            val = randi(50);
            MVO.Pop_init(i,j) = val;  %update value of element in solution array
        end

        %Calculate Cost of initial Solution:
        total_Cost = 0;
        total_Power = 0;
        total_Cof= 0;
        total_Num = 0;
        total_Area = 0;
        for j = 1:Panels.num  %for all the elements in the array
            total_Cost = total_Cost + Panels.DB(j,4)*MVO.Pop_init(i,j)*MVO.Pop_initB(i,j); %total cost calculation
            total_Power = total_Power + Panels.DB(j,1)*MVO.Pop_init(i,j)*MVO.Pop_initB(i,j)*0.9; %total energy calculation
            total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*MVO.Pop_init(i,j)*MVO.Pop_initB(i,j); %total temp coeffiecient calculation
            total_Num = total_Num + MVO.Pop_init(i,j)*MVO.Pop_initB(i,j); %total number of panels in current solution calculation
            total_Area = total_Area + MVO.Pop_init(i,j)*Panels.DB(j,2)*MVO.Pop_initB(i,j); %Total area calculation
        end
        area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
        loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
        if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
            MVO.Feasible = 1; %if feasible set flag to 1
            %MVO.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste; %calculating cost of initial sol
            MVO.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2; %calculating cost of initial sol
        end
    end
end

MVO.Pop_current = MVO.Pop_init;  %save initial population as the current population
MVO.Pop_currentB = MVO.Pop_initB;
clear i j k m cell val total_Area total_Num total_Cof total_Power total_Cost;  %clear used temporary variables

%Plot the cost function figure
Figures.Main_fig = figure;  %create new figure
Best_Cost_Array = zeros(1,MVO.Num_Generations+1);
Best_Cost_Array(1,1) = min(MVO.Pop_cost);  %save the cost of initial solution
figure(Figures.Main_fig);  %access the main figure
subplot(1,2,1);  %in the 1st location in the figure
Figures.Cost = plot(0,Best_Cost_Array(1,1),'m*','LineWidth',1.6); %plot the initial cost
xlim([0,MVO.Num_Generations]); %define limits in X-direction
grid on; %initialize the grid to be on
axis square; %make the axes look like square
xlabel('Number of Generations');
ylabel('Cost function');
title('The Cost function');  %give title to figure
legend('Cost function');

MVO_new_cost = MVO.Pop_cost;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the Multi-Verse Optimizer Loop here   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gen = 1:MVO.Num_Generations  %for the maximum number of generations
    if stopping == 200
        break;
    end
    [~, member_index] = sort(MVO_new_cost);
    bestSol = MVO_new_cost(member_index(1),1);
    Best_Cost_Array(1,gen+1) = bestSol;

    %Initialize a new generation matrix to be filled by the newly created members:
    MVO_new = zeros(MVO.Pop_size,Panels.num);  %initialize an empty matrix
    MVO_newB = zeros(MVO.Pop_size,Panels.num);  %initialize an empty matrix
    MVO_new_cost = zeros(MVO.Pop_size,1);  %initialize an empty array of zeros

    WEP = MVO.WEP_Min + gen*((MVO.WEP_Max - MVO.WEP_Min)/MVO.Num_Generations);
    TDR = 1-((gen)^(1/6)/(MVO.Num_Generations)^(1/6));

    if gen == 1 %if this is the initial generation
        MVO.Pop_current = MVO.Pop_init;  %read the initial population
        MVO.Pop_currentB = MVO.Pop_initB;  %read the initial population
        MVO.Pop_cost = MVO_new_cost;  %Read the initial population cost
        Best_SolA = MVO.Pop_init(member_index(1),:);
        Best_SolB = MVO.Pop_initB(member_index(1),:);
    end

    for i = 1:MVO.Pop_size
        MVO.Feasible = 0;
        while(MVO.Feasible == 0)
            universeB = MVO.Pop_currentB(i,:);
            MVO_new(i,:) = zeros(1,Panels.num);
            MVO_newB(i,:) = zeros(1,Panels.num);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Randomizing Binary 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             for counter = 1:10
            %                 ranCell2 = 0;
            %                 ranCell = 0;
            %                 while ranCell2 == ranCell
            %                     ranCell2 = randi([1,Panels.num]);
            %                     ranCell =  randi([1,Panels.num]);
            %                 end
            %                 saveB1 = universeB(1,ranCell2);
            %                 saveB2 = universeB(1,ranCell);
            %                 universeB(1,ranCell2) = saveB2;
            %                 universeB(1,ranCell) = saveB1;
            %             end
            %            MVO_newB(i,:) = universeB;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Randomizing Binary 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1:Panels.num
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Randomizing Binary 2
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ranCell =  randi([1,Panels.num]);
                saveB1 = universeB(1,j);
                saveB2 = universeB(1,ranCell);
                universeB(1,j) = saveB2;
                universeB(1,ranCell) = saveB1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Randomizing Binary 2
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                r1 = rand();
                fitness_sum = 0;
                for k = 1:MVO.Pop_size
                    fitness_sum = fitness_sum + MVO_new_cost(k);
                end
                norm_const = 1/sqrt(fitness_sum);
                if(norm_const*MVO_new_cost(i)>r1)
                    white_hole = randi([1,MVO.Pop_size]);
                    MVO_new(i,j) = MVO.Pop_current(member_index(white_hole),j);
                end
                r2 = rand();
                if(r2<WEP)
                    r3 = rand();
                    r4 = rand();
                    ub = max(max(MVO.Pop_current));
                    lb = min(min(MVO.Pop_current));
                    if(r3<0.5)
                        MVO_new(i,j) = round(MVO.Pop_current(member_index(1),j) + TDR*((ub-lb)*r4 + lb));

                    else
                        MVO_new(i,j) = round(MVO.Pop_current(member_index(1),j) - TDR*((ub-lb)*r4 + lb));
                    end
                    if(MVO_new(i,j)<0 || MVO_new(i,j)>50)
                        MVO_new(i,j) = round(abs(mod(MVO_new(i,j),50)));
                    end
                end
            end
            MVO_newB(i,:) = universeB;
            total_Cost = 0;
            total_Power = 0;
            total_Cof= 0;
            total_Num = 0;
            total_Area = 0;
            for j = 1:Panels.num  %for all the elements in the
                total_Cost = total_Cost + Panels.DB(j,4)*MVO_new(i,j)*MVO_newB(i,j); %total cost calculation
                total_Power = total_Power + Panels.DB(j,1)*MVO_new(i,j)*MVO_newB(i,j)*0.9; %total energy calculation
                total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*MVO_new(i,j)*MVO_newB(i,j); %total temp coeffiecient calculation
                total_Num = total_Num + MVO_new(i,j)*MVO_newB(i,j); %total number of panels in current solution calculation
                total_Area = total_Area + MVO_new(i,j)*Panels.DB(j,2)*MVO_newB(i,j); %Total area calculation
            end
            area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
            loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
            if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
                MVO.Feasible = 1;
                %MVO_new_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste;
                MVO_new_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
                MVO.Pop_current(i,:) = MVO_new(i,:);
                MVO.Pop_currentB(i,:) = MVO_newB(i,:);
            end
        end
    end

    [~, member_index] = sort(MVO_new_cost);
    bestSol = min(MVO_new_cost);
    if(bestSol<Best_Cost_Array(1,gen))
        Best_Cost_Array(1,gen+1) = bestSol;
        Best_SolB(1,:) = MVO_newB(member_index(1),:);
        Best_SolA(1,:) = MVO_new(member_index(1),:);
    else
        Best_Cost_Array(1,gen+1) = Best_Cost_Array(1,gen);
    end
    if Best_Cost_Array(1,gen+1) == Best_Cost_Array(1,gen)
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


finalSol = zeros(1,Panels.num);
for j = 1:Panels.num
    finalSol(1,j) = Best_SolB(1,j)*Best_SolA(1,j);
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
bestSol= total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
disp(['Best cost is:    ' num2str(total_Cost)] + " $");
disp(['Best power is:   ' num2str(total_Power)] + " W");
disp(['Best area is:   ' num2str(total_Area)] + " m^2");
disp(['Best area wasted is:   ' num2str(area_Waste*100)] + " %");
disp(['Best power loss is:   ' num2str(loss)] + " %");
disp(['Best cost function is: ' num2str(bestSol)]);
timElapsed = toc;
disp(['Time of execution: ' num2str(timElapsed)] + " s");







