clear
% Read the excel from iso_corr; Please feed in EVEN CHAIN fatty acids!
input_path = '';

% inital the name for output file
input_path_size = size(input_path);
input_pathLlen = input_path_size(1,2);
filename = strcat(input_path(1:input_pathLlen-5),'_MID.xlsx');

%% Change tracer enrichment 
% output : sheet1 measured MID data of each compound
% output : sheet2 C enrichment
% output : sheet3 simulated MID
% output : sheet4 fraction abundance of AcCoA
% output : sheet5 DNL (de novo fatty acids synthesis) and SSR of the
% simulation

% Find the indexes which helps locate compounds of interest
[RAW] = readcell(input_path,'sheet', 'enrichment');
[absolute] = readcell(input_path,'sheet', 'absolte');
[total_ion] = readcell(input_path,'sheet', 'total ion');
compound_list = cell2mat(unique(RAW(2:end,1),'stable'));
[r_c,c_c] = size(compound_list);
[r_raw,c_raw] = size(RAW);




new_raw(1,:) = RAW(1,:);
for i = 1:r_c
fa_species = compound_list(i,:);
fat_length = regexp(fa_species ,'\d*','Match');
fat_length = str2double(fat_length(1,1));
index1 = find(strcmp(RAW(:,1),fa_species ));
index2 = find(strcmp(RAW(:,11),strcat('C12 PARENT')));
position = intersect(index1,index2);
new_raw (end+1,:) =  RAW(position,:);
for l = 1: fat_length
index1 = find(strcmp(RAW(:,1),fa_species ));
index2 = find(strcmp(RAW(:,11),strcat('C13-label-',int2str(l))));
position = intersect(index1,index2);
if isempty( position)
fittting_row = num2cell(zeros(1,c_raw));
fittting_row(1,[1 2  11]) = [{fa_species} {fa_species}  {strcat('C13-label-',int2str(l))}];
new_raw (end+1,:) =  fittting_row ;
else
new_raw (end+1,:) =  RAW(position,:);
end
end
end
RAW = new_raw;

%% define sim_MID sheet
RAW2 = RAW;
starting_column = 16; % the column index of your data

%% def carbon enrichment sheet 
%% M+n sum(n*n_enrichment)/#total carbon

enrich_sheet(1,:) = [RAW(1,2) RAW(1,starting_column:end)];


%% define fit result sheet (building block contribution)
fit_sheet(1,:) = [RAW(1,2) {'Fraction'} RAW(1,starting_column:end)];
fraction_columns = [{'M+0'} {'M+1'} {'M+2'} {'DNL'}]';

%% define DNL_ssr sheet
ssr_sheet(1,:) = [RAW(1,2) {'Type'} RAW(1,starting_column:end)];

for i = 1:r_c
index = find(strcmp(RAW(:,1),compound_list(i,:)));
compound_position(i,1) = index(1);
compound_position(i,2) = index(end);
MID = cell2mat(RAW(compound_position(i,1):compound_position(i,2),starting_column:end));
[r_m,c_m] = size(MID);
[sim_mid,ssr_fit,fit] = MID_simulation(MID);
% SIM_mid sheet
RAW2(compound_position(i,1):compound_position(i,2),starting_column:end) = num2cell(sim_mid);



% Fit result sheet
fit_sheet((i-1)*4+2:i*4+1,2) = fraction_columns;
compound_type = repmat({compound_list(i,:)},4,1);
fit_sheet((i-1)*4+2:i*4+1,1) = compound_type;
fit_sheet((i-1)*4+2:i*4+1,3:end) =  num2cell(fit(1:4,:));



% Enrichment sheet
carbon_list = 0:r_m -1;
enrichment = carbon_list * MID/(r_m -1);%./tracer_enrichment; %%% change the tracer enrichment here!
enrich_sheet(i+1,:) = [compound_list(i,:)  num2cell(enrichment)];
% DNL and SSR sheet
ssr_sheet((i-1)*2+2,:) = [compound_list(i,:) {'DNL'} num2cell(fit(4,:))];
ssr_sheet((i-1)*2+3,:) = [compound_list(i,:) {'SSR'} num2cell(ssr_fit)];
end


cellArray = fit_sheet(:,1:2);
rowToFind = {'C16:0', 'M+2'};
rowIndex = find(strcmp(cellArray(:, 1), rowToFind{1}) & strcmp(cellArray(:, 2), rowToFind{2}));

%% Change tracer enrichment 
%%tracer_enrichment =  0.25*ones(1,c_raw -15);
disp('Enrichment norm to C16:0 M+2 AcCoA, check the numbers!')
tracer_enrichment =  fit_sheet(rowIndex ,3:end);
cell2mat(tracer_enrichment)
enrich_result_raw = cell2mat(enrich_sheet(2:end,2:end));
enrich_sheet(2:end,2:end) = num2cell(enrich_result_raw./cell2mat(tracer_enrichment));

% Enrichment average and stand
enrich_std_mean = [enrich_sheet(:,1) cal_mean_std(enrich_sheet(:,2:end))];


RAW(cellfun(@(x) any(ismissing(x)), RAW)) = {'NaN'};
RAW2(cellfun(@(x) any(ismissing(x)), RAW2)) = {'NaN'};

absolute(cellfun(@(x) any(ismissing(x)), absolute)) = {'NaN'};
total_ion(cellfun(@(x) any(ismissing(x)), total_ion)) = {'NaN'};
writecell(RAW,filename,'Sheet','MID')
writecell(absolute,filename,'Sheet','absolute')
writecell(total_ion,filename,'Sheet','total ion')
writecell(enrich_sheet,filename,'Sheet','carbon_enrichment')
writecell(enrich_std_mean,filename,'Sheet','enrich_mean')
writecell(RAW2,filename,'Sheet','Sim_MID')
writecell(fit_sheet,filename,'Sheet','fraction contribution')
%% fraction contribution mean
writecell([fit_sheet(:,1:2) cal_mean_std(fit_sheet(:,3:end))],filename,'Sheet','frac contribution_mean')

%% dnl_ssr mean
writecell(ssr_sheet,filename,'Sheet','DNL_SSR')
writecell([ssr_sheet(:,1:2) cal_mean_std(ssr_sheet(:,3:end))],filename,'Sheet','DNL_SSR_mean')

%% MID, sim MID mean
writecell([RAW(:,[1 11 12]) cal_mean_std(RAW(:,starting_column :end))],filename,'Sheet','MID_mean')
writecell([RAW2(:,[1 11 12]) cal_mean_std(RAW2(:,starting_column :end))],filename,'Sheet','MID_mean')

disp('Job done :P')

disp('Please verify your replicates here !')
replicates = find_replicates(fit_sheet(1,3:end))





function [sim_mid,ssr_fit,x] = MID_simulation(MID)
[r,c] = size(MID);
MID = MID./sum(MID);
for i = 1:c
    data = MID(:,i);
    weight = eye(r);
     % convolution
    mid1 = @(x) multi_conv(x,r/2-1);
    % simulation
    sim1 = @(x) mid1(x(1:3))*x(4) + (1-x(4))*[1;zeros([r-1,1])]; 
    % calculate SSR
    ssr1 = @(x) (sim1(x)-data)'*weight*(sim1(x)-data);
    % linear programming
    x0 = [0.95;0.02;0.03;0.1]; %%% You can change the initial condition here
    Aeq = [1 1 1 0];
    beq = 1;
    % out put
    fit1 = fmincon(ssr1,x0,[],[],Aeq,beq,[0 0 0 0],[1 1 1 1]);
    x(:,i) = fit1;
    sim_mid(:,i) = sim1(fit1);
    ssr_fit(:,i) = ssr1(fit1);
end
end

function conv_re = multi_conv(x,layers)
conv_re  = x;
for i = 1:layers
    conv_re  = conv(conv_re ,x);
end
end

function replicates = find_replicates(lst)
    % Takes a cell array of strings and returns a struct where each field is a prefix that has replicates,
    % and the corresponding value is a cell array of strings representing the replicates for that prefix.
    
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Iterate over the list
    for i = 1:numel(lst)
        % Extract the prefix and suffix from the string
        parts = strsplit(lst{i}, '_');
        prefix = strjoin(parts(1:end-1), '_');
        
        % Add the item to the corresponding group based on its prefix
        if isKey(groups, prefix)
            groups(prefix) = [groups(prefix), lst{i}];
        else
            groups(prefix) = {lst{i}};
        end
    end
    
    % Identify prefixes with replicates
    replicates = struct();
    keys = groups.keys();
    for i = 1:numel(keys)
        prefix = keys{i};
        group = groups(prefix);
        if numel(group) > 1
            replicates.(prefix) = group;
        end
    end
end


function mean_std_cell = cal_mean_std(my_cell)


my_list = my_cell(1, :);
% Call the find_replicates function to get the replicates struct
replicates = find_replicates(my_list);
keys = fieldnames(replicates);

% Calculate the average of each row's elements based on the group they belong to
%num_rows = size(my_cell, 1) - 1;
averages = num2cell(zeros(size(my_cell, 1) , numel(keys)));
stds = num2cell(zeros(size(my_cell, 1) , numel(keys)));
for i = 1:numel(keys)
    prefix = keys{i};
    group = replicates.(prefix);
    group_elements = cellfun(@(x) find(strcmp(my_list, x)), group);
    averages(2:end,i) = num2cell(mean(cell2mat(my_cell(2:end,group_elements)),2));
    averages(1,i) = {strcat(prefix,'_mean')};
    stds(2:end,i) = num2cell(std(cell2mat(my_cell(2:end,group_elements)),0,2));
    stds(1,i) = {strcat(prefix,'_std')};
end
mean_std_cell = [averages stds ];
end

