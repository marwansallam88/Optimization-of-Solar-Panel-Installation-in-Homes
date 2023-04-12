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
GA.Num_Generations = 1000;  %maximum number of generations to be tested
GA.Pop_size = 10;  %the population size
GA.Elite_ratio = 0.1;  %percentage of elitism
GA.CrossOver_ratio = 0.6; %percentage of cross over processes
GA.Mutation_ratio = 1 - GA.Elite_ratio - GA.CrossOver_ratio; %the rest of mutation ratio
GA.Alpha = 0.6;  %alpha used for cross over process to generate new children
GA.Noise_Scale = 0.5;  %used for mutation to add noise around specific genes, number from [0 to 1]
GA.Feasible = 0;
stopping = 0;

%Initialize a random population:
GA.Pop_init = zeros(GA.Pop_size,Panels.num);  %initialize an empty matrix
GA.Pop_initB = zeros(GA.Pop_size,Panels.num);
GA.Pop_cost = zeros(GA.Pop_size,1);  %initialize an empty array of zeros
Best_SolA = zeros(1,Panels.num);
Best_SolB = zeros(1,Panels.num);

for i = 1:GA.Pop_size  %for all the population members
    GA.Feasible = 0;
    while(GA.Feasible == 0) %loop until sol is feasible
        GA.Pop_init(i,:) = zeros(1,Panels.num);  %initialize an empty matrix
        GA.Pop_initB(i,:) = zeros(1,Panels.num);
        GA.Pop_cost(i,1) = 0;
        cell1 = randi([1,10]);
        GA.Pop_initB(i,cell1) = 1;
        cell2 = randi([11,20]);
        GA.Pop_initB(i,cell2) = 1;
        cell3 = randi([21,30]);
        GA.Pop_initB(i,cell3) = 1;
        cell4 = randi([31,40]);
        GA.Pop_initB(i,cell4) = 1;
        cell5 = randi([41,50]);
        GA.Pop_initB(i,cell5) = 1;
        for j = 1:Panels.num  %for all the elements in the array
            val = randi(50);
            GA.Pop_init(i,j) = val;  %update value of element in solution array
        end

        %Calculate Cost of initial Solution:
        total_Cost = 0;
        total_Power = 0;
        total_Cof= 0;
        total_Num = 0;
        total_Area = 0;
        for j = 1:Panels.num  %for all the elements in the array
            total_Cost = total_Cost + Panels.DB(j,4)*GA.Pop_init(i,j)*GA.Pop_initB(i,j); %total cost calculation
            total_Power = total_Power + Panels.DB(j,1)*GA.Pop_init(i,j)*GA.Pop_initB(i,j)*0.9; %total energy calculation
            total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*GA.Pop_init(i,j)*GA.Pop_initB(i,j); %total temp coeffiecient calculation
            total_Num = total_Num + GA.Pop_init(i,j)*GA.Pop_initB(i,j); %total number of panels in current solution calculation
            total_Area = total_Area + GA.Pop_init(i,j)*Panels.DB(j,2)*GA.Pop_initB(i,j); %Total area calculation
        end
        area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
        loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
        if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
            GA.Feasible = 1; %if feasible set flag to 1
            if i==1
                Best_Sol = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
                for j = 1:Panels.num
                    Best_SolA(1,j) = GA.Pop_init(i,j);
                    Best_SolB(1,j) = GA.Pop_initB(i,j);
                end
            end
            GA.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2; %calculating cost of initial sol
            if Best_Sol > GA.Pop_cost(i,1)
                Best_Sol = GA.Pop_cost(i,1);
                for j = 1:Panels.num
                    Best_SolA(1,j) = GA.Pop_init(i,j);
                    Best_SolB(1,j) = GA.Pop_initB(i,j);
                end
            end
        end
    end
end
GA.Pop_current = GA.Pop_init;  %save initial population as the current population
GA.Pop_currentB = GA.Pop_initB;
clear i j k m cell val total_Area total_Num total_Cof total_Power total_Cost;  %clear used temporary variables

%Plot the cost function figure
Figures.Main_fig = figure;  %create new figure
Elite_Cost_Array = zeros(1,GA.Num_Generations+1);
Elite_Cost_Array(1,1) = min(GA.Pop_cost);  %save the cost of initial solution
figure(Figures.Main_fig);  %access the main figure
subplot(1,2,1);  %in the 1st location in the figure
Figures.Cost = plot(0,Elite_Cost_Array(1,1),'m*','LineWidth',1.6); %plot the initial cost
xlim([0,GA.Num_Generations]); %define limits in X-direction
grid on; %initialize the grid to be on
axis square; %make the axes look like square
xlabel('Number of Generations');
ylabel('Cost function');
title('The Cost function');  %give title to figure
legend('Cost function');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the GA optimization Loop here   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for gen = 1:GA.Num_Generations  %for the maximum number of generations
    if stopping == 200
        break;
    end

    if gen == 1 %if this is the initial generation
        Gen_old = GA.Pop_init;  %read the initial population
        Gen_oldB = GA.Pop_initB;  %read the initial population
        Gen_old_cost = GA.Pop_cost;  %Read the initial population cost
    end

    %Initialize a new generation matrix to be filled by the newly created members:
    Gen_new = zeros(GA.Pop_size,Panels.num);  %initialize an empty matrix
    Gen_newB = zeros(GA.Pop_size,Panels.num);  %initialize an empty matrix
    Gen_new_cost = zeros(GA.Pop_size,1);  %initialize an empty array of zeros

    %sort the current population based on cost:
    [~, member_index] = sort(Gen_old_cost);  %sort and save index of each element
    Num_elite = round(GA.Elite_ratio*GA.Pop_size); %calculate required number of elite elements
    for i = 1:Num_elite  %for the required number of elites
        Gen_new(i,:) = Gen_old(member_index(i),:);  %get the best members in prev generatio
        Gen_newB(i,:) = Gen_oldB(member_index(i),:);
        Gen_new_cost(i,1) = Gen_old_cost(member_index(i),1);  %save elites cost
    end

    %calculate the required number of crossover children:
    Num_CO_children = round(GA.CrossOver_ratio*GA.Pop_size);    %calculate required number of children elements
    for i = Num_elite+1 : 2: Num_elite+Num_CO_children  %for the number of following members (with steps of 2)
        GA.Feasible = 0;
        sol1_found = 0;
        sol2_found = 0;
        save1 = zeros(1,Panels.num);
        save1B = zeros(1,Panels.num);
        save2 = zeros(1,Panels.num);
        save2B = zeros(1,Panels.num);
        while(GA.Feasible == 0)
            parent_1_num = randi([1,GA.Pop_size]); %randomly select a member from the generation
            parent_1 = Gen_old(parent_1_num,:);  %read the info of this parent
            parent_1B = Gen_oldB(parent_1_num,:);  %read the info of this parent
            parent_2_num = randi([1,GA.Pop_size]); %randomly select a member from the generation
            while(parent_2_num == parent_1_num)
                parent_2_num = randi([1,GA.Pop_size]);
            end
            parent_2 = Gen_old(parent_2_num,:);  %read the info of this parent
            parent_2B = Gen_oldB(parent_2_num,:);  %read the info of this parent

            %initialize two children arrays empty:
            child_1 = zeros(1,Panels.num);  %initialize an empty array
            child_2 = zeros(1,Panels.num);  %initialize an empty array
            child_1B = zeros(1,Panels.num);  %initialize an empty array
            child_2B = zeros(1,Panels.num);  %initialize an empty array

            for j = 1:Panels.num  %for all the elements in the child array
                child_1(1,j) = round(GA.Alpha*parent_1(1,j) + (1-GA.Alpha)*parent_2(1,j));  %calculate gene for child-1
                child_2(1,j) = round((1-GA.Alpha)*parent_1(1,j) + GA.Alpha*parent_2(1,j));  %calculate gene for child-2
            end

            for j = 1:Panels.num  %for all the elements in the child array
                cross = randi(Panels.num);
                half11 = parent_1B(1:cross);
                half12 = parent_1B(cross+1:Panels.num);
                half21 = parent_2B(1:cross);
                half22 = parent_2B(cross+1:Panels.num);
                child_1B(1,:) = [half11 half22];
                child_2B(1,:) = [half21 half12];
                if(nnz(child_1B)>5)
                    num_of_1s = nnz(child_1B);
                    indicies = find(child_1B,10);
                    for k = 1:(num_of_1s-5)
                        index = indicies(1,k);
                        child_1B(1,index) = 0;
                        while(child_2B(1,index) == 1)
                            index = index+1;
                            if(index>50)
                                index = 0;
                            end
                        end
                        child_2B (1,index) = 1;
                    end
                elseif(nnz(child_2B)>5)
                    num_of_1s = nnz(child_2B);
                    indicies = find(child_2B,10);
                    for k = 1:(num_of_1s-5)
                        index = indicies(1,k);
                        child_2B(1,index) = 0;
                        while(child_1B(1,index) == 1)
                            index = index+1;
                            if(index>50)
                                index = 0;
                            end
                        end
                        child_1B (1,index) = 1;
                    end
                end
            end
            if(sol1_found == 0)
                Gen_new(i,:) = child_1;  %save child 1 to the generation
                Gen_newB(i,:) = child_1B;  %save child 1 to the generation
                total_Cost = 0;
                total_Power = 0;
                total_Cof= 0;
                total_Num = 0;
                total_Area = 0;
                for j = 1:size(GA.Pop_init,2)  %for all the elements in the
                    total_Cost = total_Cost + Panels.DB(j,4)*Gen_new(i,j)*Gen_newB(i,j); %total cost calculation
                    total_Power = total_Power + Panels.DB(j,1)*Gen_new(i,j)*Gen_newB(i,j)*0.9; %total energy calculation
                    total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*Gen_new(i,j)*Gen_newB(i,j); %total temp coeffiecient calculation
                    total_Num = total_Num + Gen_new(i,j)*Gen_newB(i,j); %total number of panels in current solution calculation
                    total_Area = total_Area + Gen_new(i,j)*Panels.DB(j,2)*Gen_newB(i,j); %Total area calculation
                end
                area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
                loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
                if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
                    sol1_found= 1;
                    Gen_new_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
                    if Best_Sol > Gen_new_cost(i,1)
                        Best_Sol = Gen_new_cost(i,1);
                        for j = 1:Panels.num
                            Best_SolA(1,j) = Gen_new(i,j);
                            Best_SolB(1,j) = Gen_newB(i,j);
                        end
                    end
                end

            end
            if(sol2_found == 0)
                Gen_new(i+1,:) = child_2;  %save child 2 to the generation
                Gen_newB(i+1,:) = child_2B;  %save child 2 to the generation
                total_Cost = 0;
                total_Power = 0;
                total_Cof= 0;
                total_Num = 0;
                total_Area = 0;
                for j = 1:size(GA.Pop_init,2)  %for all the elements in the
                    total_Cost = total_Cost + Panels.DB(j,4)*Gen_new(i+1,j)*Gen_newB(i+1,j); %total cost calculation
                    total_Power = total_Power + Panels.DB(j,1)*Gen_new(i+1,j)*Gen_newB(i+1,j)*0.9; %total energy calculation
                    total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*Gen_new(i+1,j)*Gen_newB(i+1,j); %total temp coeffiecient calculation
                    total_Num = total_Num + Gen_new(i+1,j)*Gen_newB(i+1,j); %total number of panels in current solution calculation
                    total_Area = total_Area + Gen_new(i+1,j)*Panels.DB(j,2)*Gen_newB(i+1,j); %Total area calculation
                end
                area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
                loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
                if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
                    sol2_found = 1;
                    Gen_new_cost(i+1,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
                    if Best_Sol > Gen_new_cost(i+1,1)
                        Best_Sol = Gen_new_cost(i+1,1);
                        for j = 1:Panels.num
                            Best_SolA(1,j) = Gen_new(i,j);
                            Best_SolB(1,j) = Gen_newB(i,j);
                        end
                    end
                end
            end
            GA.Feasible = sol1_found && sol2_found;
        end
    end

    %calculate the required number of mutated children:
    Num_Mutants = GA.Pop_size - Num_elite - Num_CO_children;  %the rest of population are mutants
    for i = Num_elite + Num_CO_children + 1 : GA.Pop_size  %for the rest of population size members
        %we need to select the worst members to mutate, so in this code we will use the sorted cost array
        %to select the worst members and add random noise to their values hoping they will improve.
        member = Gen_old(member_index(end-i+1),:);  %select the elements from the end of list
        memberB = Gen_oldB(member_index(end-i+1),:);  %select the elements from the end of list
        
        %Randomize a new member:
        %==========================
        Scaler = GA.Noise_Scale; %number [0-1] to define scale of exploration around current solution
        member_new = zeros(size(member,1), size(member,2));  %initialize an empty array
        member_newB = zeros(size(member,1), size(member,2));  %initialize an empty array
        GA.Feasible = 0;
        while(GA.Feasible == 0)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mutate Binary
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            member_newB = memberB;
            for counter = 1:10
                ranCell2 = 0;
                ranCell = 0;
                while ranCell2 == ranCell
                    ranCell2 = randi([1,Panels.num]);
                    ranCell =  randi([1,Panels.num]);
                end
                saveB1 = member_newB(1,ranCell2);
                saveB2 = member_newB(1,ranCell);
                member_newB(1,ranCell2) = saveB2;
                member_newB(1,ranCell) = saveB1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Mutate Arithmatic
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            member_new = member;
            noise = round(randi([-10,10],[1,size(member_new,2)]));
            member_new= member + noise;
            for j = 1:size(member_new,2)
                while member_new(j) <0
                    member_new(j) = member_new(j) + 50;
                end
                while member_new(j) > 50
                    member_new(j) = member_new(j) - 50;
                end
            end

            total_Cost = 0;
            total_Power = 0;
            total_Cof= 0;
            total_Num = 0;
            total_Area = 0;
            for j = 1:size(GA.Pop_init,2)  %for all the elements in the
                total_Cost = total_Cost + Panels.DB(j,4)*member_new(j)*member_newB(j); %total cost calculation
                total_Power = total_Power + Panels.DB(j,1)*member_new(j)*member_newB(j)*0.9; %total energy calculation
                total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*member_new(j)*member_newB(j); %total temp coeffiecient calculation
                total_Num = total_Num + member_new(j)*member_newB(j); %total number of panels in current solution calculation
                total_Area = total_Area + member_new(j)*Panels.DB(j,2)*member_newB(j); %Total area calculation
            end
            area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
            loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
            if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
                GA.Feasible = 1; %if feasible set flag to 1
                Gen_new_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2;
                if Best_Sol > Gen_new_cost(i,1)
                    Best_Sol = Gen_new_cost(i,1);
                    for j = 1:Panels.num
                        Best_SolA(1,j) = member_new(j);
                        Best_SolB(1,j) = member_newB(j);
                    end
                end
                Gen_new(i,:) = member_new;  %save MUTANT to the generation
                Gen_newB(i,:) = member_newB;  %save MUTANT to the generation
            end
        end
    end
    GA.Pop_current = Gen_new;
    GA.Pop_currentB = Gen_newB;
    GA.Pop_cost = Gen_new_cost;

    %Update variables for next iteration:
    Gen_old = Gen_new;  %update the current generation
    Gen_oldB = Gen_newB;  %update the current generation
    Gen_old_cost = Gen_new_cost;  %update teh current cost
    [~, member_index2] = sort(Gen_new_cost);
    bestSol = min(Gen_new_cost(member_index2(1)));
    Elite_Cost_Array(1,gen+1) = bestSol;
    if Elite_Cost_Array(1,gen) == bestSol
        stopping = stopping + 1;
    else
        stopping = 0;
    end

    %cleaning up used variables so far:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GA.Num_Elites = Num_elite;
    GA.Num_CrossOver_children = Num_CO_children;
    GA.Mutants = Num_Mutants;
    clear child_1 child_2 i j member member_index member_new;
    clear Num_CO_children Num_elite Num_Mutants parent_1 parent_2;
    clear parent_1_num parent_2_num Scaler val;

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
    set(Figures.Cost,'XData',0:1:gen,'YData',Elite_Cost_Array(1,1:gen+1));
    subplot(1,2,1); title(['Current Cost Value: ' num2str(Elite_Cost_Array(1,gen+1))]);
    pause(0.05);
    
    %Plotting pie chart
    labels = {num2str(PanelCount(1,1)), num2str(PanelCount(1,2)), num2str(PanelCount(1,3)), num2str(PanelCount(1,4)), num2str(PanelCount(1,5))};
    subplot(1,2,2);
    pie(PanelCount, labels);
    labels1 = {string(Panels1.DB(Paneli(1),1)),string(Panels1.DB(Paneli(2),1)),string(Panels1.DB(Paneli(3),1)),string(Panels1.DB(Paneli(4),1)),string(Panels1.DB(Paneli(5),1))};
    legend(labels1,'Location','southoutside','Orientation','vertical');
end

[~, member_index1] = sort(Gen_new_cost);
finalSol = zeros(1,Panels.num);
for j = 1:Panels.num
    finalSol(1,j) = Gen_new(member_index1(1),j)*Gen_newB(member_index1(1),j);
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
disp(['Best cost function is: ' num2str(Best_Sol)]);
timElapsed = toc;
disp(['Time of execution: ' num2str(timElapsed)]);









