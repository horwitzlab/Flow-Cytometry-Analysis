% Flow Cytometry Analysis
% Based on FACSworkup.m
% GDLH 12/20/19
clear all
% Local directory
% folder_name = 'C:\Users\Matt\OneDrive\College\Summer 2021 Horwitz Lab\2020 04 19 Chihiro screen AAV1\2020 04 19 Chihiro screen AAV1\CSV 04222020';

% Tower directory
folder_name = 'Z:\Viral Vectors\AAV Development\AAV Clinical';
old_dir = pwd;
disp(['Currently in ',old_dir])
cd(folder_name);
[fname_nc,pathname] = uigetfile('*.*','Which file for negative control?','MultiSelect','off');
if fname_nc == 0
   error('No negative control file selected. Aborting.'); 
end
Data = importdata([pathname,filesep,fname_nc]); 
cd(old_dir);

% Gating the events
FSCA = Data.data(:,strcmp(Data.textdata,'"FSC-A"')); % Forward scatter absorbance
FSCH = Data.data(:,strcmp(Data.textdata,'"FSC-H"')); % Forward scatter height
FSCW = Data.data(:,strcmp(Data.textdata,'"FSC-W"')); % Forward scatter width
SSCA = Data.data(:,strcmp(Data.textdata,'"SSC-A"')); % Side scatter absorbance
GFP = Data.data(:,strcmp(Data.textdata,'"GFP-A"')); % GFP fluorescence
DAPI = Data.data(:,strcmp(Data.textdata,'"DAPI-A"')); % DAPI fluorescence
MCHERRY = Data.data(:,strcmp(Data.textdata,'"mCherry-A"')); % mCherry fluorescence
TIME = Data.data(:,strcmp(Data.textdata,'"TIME"')); % Time


%% Gate 1
h_fig = figure;
hold on; box off
plot(FSCA,SSCA,'.k','MarkerSize',1)
set(gca,'TickDir','out')
title('Draw a Polygon around the Cluster of Interest')
xlabel('Forward Scatter')
ylabel('Side Scatter');
btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
        'Position', [20 20 60 20],...
        'Callback', 'uiresume');
h = impoly; % Allows the user to interactively draw a polygon on the data.

% Pausing the analysis until the gate is completed and the user decides to
% continue
uiwait

% Calling the vertices of the polygon 'Gate1'
Gate1 = getPosition(h);

% Defining the x and y components of the Gate
xv = Gate1(:,1);
yv = Gate1(:,2);

% inpolygon takes inputs of the x and y data as well as the x and y
% positions of the gate to output a logical vector which indicates whether a
% particular point was included or excluded from the center of the polygon
% defined by the gate.
Lgate1 = inpolygon(FSCA,SSCA,xv,yv);

%% Gate 2
% Plotting a new subspace of the reduced data set
figure(h_fig);
clf;
MyClr = [0,0.8,0.6];
hold on; box off
plot(FSCH(Lgate1),FSCW(Lgate1),'.k','MarkerSize',1.5)
set(gca,'TickDir','out')
title('Draw a Polygon around the Cluster of Interest')
xlabel('Forward Scatter Height')
ylabel('Forward Scatter Width')

% Creating a button that will resume the analysis after being pressed
btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
        'Position', [20 20 60 20],...
        'Callback', 'uiresume');
    
% Drawing a new polygon
i = imrect;

% Pausing the analysis until the gate is completed and the user decides to
% continue
uiwait
% Establishing the coordinates for gate 2 given from impoly
Gate2 = getPosition(i);
Gate2 = [Gate2(1) Gate2(2)
         Gate2(1)+Gate2(3) Gate2(2)
         Gate2(1)+Gate2(3) Gate2(2)+Gate2(4)
         Gate2(1) Gate2(2)+Gate2(4)];
     
% Lgate_both = inpolygon(FSCH(Lgate1),FSCW(Lgate1),Gate2(:,1),Gate2(:,2));
Lgate_both = Lgate1 & inpolygon(FSCH,FSCW,Gate2(:,1),Gate2(:,2));


%%
% Which reporter signal
reporteridx = listdlg('PromptString','Which reporter?',...
                      'SelectionMode','single',...
                      'ListString',{'GFP','mCherry'});
if reporteridx == 1
    reporter_data = GFP;
    reporter_label = 'GFP';
else
    reporter_data = MCHERRY;
    reporter_label = 'mCherry';
end
%% 
% Plotting the reporter signal for the negative control
figure(h_fig);
clf;

% L_sensible_value = reporter_data > 0;
% bins = linspace(min(log10(reporter_data(L_sensible_value & Lgate_both)))-.5, max(log10(reporter_data(L_sensible_value & Lgate_both)))+.5,100);
% [n1,x] = hist(log10(reporter_data(L_sensible_value & Lgate_both)),bins);
if any(reporter_data < 0)
    error('cannot have negative values');
end

log_gated_reporter_data = log10(reporter_data(Lgate_both));

bins = linspace(min(log_gated_reporter_data)-.5, max(log_gated_reporter_data)+.5,100);
histogram(log_gated_reporter_data,bins);

set(gca,'Xlim',[bins(1) bins(end)],'Yscale','linear');

% Determine threshold; defaults to threshold = 0.01, or the 1% highest
% fluorescence reads
threshold = inputdlg('Enter a threshold', 'Setting Threshold', [1 20], {'0.01'});
threshold = str2double(threshold{1});
threshold = prctile(log_gated_reporter_data, ((1 - threshold) * 100));

neg_cntrl_counts_and_total = [sum(log_gated_reporter_data > threshold) sum(Lgate_both)];
xlabel([reporter_label,' signal']);
ylabel('counts');
title(['Negative control: ',fname_nc]);
%%
% Loading experimental data
cd(folder_name);
[file_names,pathname] = uigetfile('*.*','Which experimental files?','MultiSelect','on');
cd(old_dir);

condition_names = cell(1, length(file_names));
if ~iscell(file_names)
    file_names = {file_names};
end

% for i = 1:length(file_names)
%    condition_names(i) = inputdlg(['Enter condition for ',file_names{i}]);
% end

ratios = zeros(length(file_names), 1);

% Assign serum dilution ratios from file name
% Requires that dilution be written as '1 to XXXX' or '1toXXXX' somwhere
% in the file name
for i = 1:length(file_names)
    curr_file = file_names{i}(~isspace(file_names{i}));
    
    % Find start of dilution value
    ratio_start = strfind(curr_file, '1to') + length('1to');
    curr_ratio = curr_file(ratio_start:end); % 'XXXX...' (dilution plus trailing text)
    ratio_end = 0;
    
    % Find start of trailing text after last digit of dilution value
    for j = 1:length(curr_ratio)
        if isnan(str2double(curr_ratio(j)))
            ratio_end = j - 1;
            break
        end
    end
    
    curr_ratio = curr_ratio(1:ratio_end); % 'XXXX' (still a string)
    ratios(i) = 1 / str2double(curr_ratio); % XXXX (numeric)
    
    % Make sure all cells are numerics and the cell is not empty
    if isempty(ratios(i))
        ratios(i) = 0;
    elseif isnan(ratios(i))
        ratios(i) = 0;
    end
end

% Graph fluorescence histograms and threshold for each run
figure(h_fig);
set(gcf,'Position',[128    50   440   850])
data = [];
for i = 1:length(file_names)
    curr_file  = file_names{i};
    Data = importdata([pathname,filesep,curr_file]);
    FSCA = Data.data(:,strcmp(Data.textdata,'"FSC-A"'));
    FSCH = Data.data(:,strcmp(Data.textdata,'"FSC-H"'));
    FSCW = Data.data(:,strcmp(Data.textdata,'"FSC-W"'));
    SSCA = Data.data(:,strcmp(Data.textdata,'"SSC-A"'));
    GFP = Data.data(:,strcmp(Data.textdata,'"GFP-A"'));
    MCHERRY = Data.data(:,strcmp(Data.textdata,'"mCherry-A"'));

    Lgate1 = inpolygon(FSCA,SSCA,xv,yv);
    Lgate_both = Lgate1 & inpolygon(FSCH,FSCW,Gate2(:,1),Gate2(:,2));
    subplot(length(file_names),1,i); hold on;
    if reporteridx == 1
        reporter_data = GFP;
    else
        reporter_data = MCHERRY;        
    end
    
    log_gated_reporter_data = log10(reporter_data(Lgate_both));
    
    if ~isempty(reporter_data)
        histogram(log_gated_reporter_data,bins);
        [n2,x] = hist(log_gated_reporter_data,bins);
        
        plot(threshold, max(n2),'kv');
        set(gca,'Xlim',[bins(1) bins(end)]);
        % data holds serum concentration, counts of fluorescence above
        % threshold, and total counts in columns 1, 2, and 3 respectively
        data = [data; ratios(i) sum(log_gated_reporter_data > threshold) sum(Lgate_both)];
        title(ratios(i));
        ylabel('counts');
        set(gca,'Yscale','log');
    end
end

% ratios = cat(1, [0], ratios);

xlabel([reporter_label,' signal']);
fit_fig = figure; axes; hold on;
plot(data(:,1), data(:,2)./data(:,3),'o','MarkerFace','black','DisplayName', 'Data');
ylabel('Proportion of Transformed Cells');
xlabel('Serum Concentration Ratio');
% set(gca,'Xtick',1:length(file_names)+1,'XtickLabel',cat(2,'NC',ratios));
%%
% Linear regression to exponential curve
response_observations = data(:,2) ./ data(:,3);
log_response_observations = log(response_observations);
design_matrix = ones(height(data), 2);
design_matrix(:,2) = data(:,1);

% exponential of the form y = beta(1)*e^(beta(2)*x)
log_betas = regress(log_response_observations, design_matrix);
log_betas(1) = exp(log_betas(1));

regression_x = linspace(0, max(data(:,1)), 100);
regression_y = log_betas(1) .* exp(log_betas(2) .* regression_x);

plot(regression_x, regression_y, 'red', 'DisplayName', 'Regression Fit');

% Linear regression using fminsearch() to find MLE
% exponential of the form y = alpha*e^(beta*x) + delta
% Includes negative control in calculation as (infinity, NC)
x_values = [inf; data(:,1)];
y_observed = [neg_cntrl_counts_and_total(1); data(:,2)];
total_observations = [neg_cntrl_counts_and_total(2); data(:,3)];
delta = neg_cntrl_counts_and_total(1) / neg_cntrl_counts_and_total(2);
log_betas(1) = data(height(data), 2) / data(height(data), 3) - delta; % Improve alpha guess
abd = [log_betas; delta]; % alpha, beta, and delta guesses for fit
abd(2,1) = min([0, abd(2,1)]); % beta should not be positive

[alpha, beta, delta] = exponential_fit(x_values, y_observed, total_observations, abd);
abd(:,1) = [alpha; beta; delta];

y_mlefit = alpha .* exp(beta .* regression_x) + delta;

plot(regression_x, y_mlefit, 'blue', 'DisplayName', 'MLE Fit');
ylim([0 max(y_mlefit)])

% Calculate x and y coordinates of the half-maximal proportion of fluorescing cells
% and associated serum concentration
[D50_x, D50_y] = calculate_D50(alpha, beta, delta);
plot(D50_x, D50_y, 'Marker', 'd', 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor','c', 'DisplayName', 'D 50')

%%
% Bootstrapped Standard Errors
num_sims = 200;
bootstrap_dist = zeros(height(x_values), num_sims); % Bootstrapped values for y_observed
for i = 1:height(x_values)
    bootstrap_dist(i,:) = binornd(total_observations(i), y_observed(i) ./ total_observations(i), ...
        1, num_sims);
end

bootstrap_fits = zeros(num_sims, 3); % Exponential fits from bootstrapped values
bootstrap_D50s = zeros(num_sims, 2); % D50's from bootstrapped exponential fits
for i = 1:num_sims
    [bootstrap_fits(i,1), bootstrap_fits(i,2), bootstrap_fits(i,3)] ...
        = exponential_fit(x_values, bootstrap_dist(:,i), total_observations, abd);
    [bootstrap_D50s(i,1), bootstrap_D50s(i,2)] ...
        = calculate_D50(bootstrap_fits(i,1), bootstrap_fits(i,2), bootstrap_fits(i,3));
end

bstrap_fig = figure;
histogram(bootstrap_D50s(:,1), 50)
title('Bootstrapped D50s')

% Use distribution of bootstrapped D50's as an approximation of standard
% error
SE_estimate = std(bootstrap_D50s);
disp(SE_estimate(1))

figure(fit_fig)
errorbar(D50_x, D50_y, SE_estimate(1), 'horizontal',...
    'DisplayName', 'Bootstrapped Error')
legend

%%
% Subfunction(s)
function [alpha, beta, delta] = exponential_fit(x_values, y_observed, total_observations, params)

    % Constrain parameter values for reasonable fit
    % Lower bound of delta is 0
    lb = [-inf, -inf, 0];
    % Upper bound of delta is the largest GFP+ proportion
    proportions = y_observed ./ total_observations;
    ub = [inf, inf, max(proportions)];

    options = optimoptions(@fmincon, 'Display', 'off');
    % Optimize log likelihood equation (MLE)
    fmin_params = fmincon(@ll, params, [], [], [], [], lb, ub, [], options);
    
    alpha = fmin_params(1);
    beta = fmin_params(2);
    delta = fmin_params(3);
    
    function likelihood = ll(input) 
        
        % Calculates log likelihood for parameters of an exponential fit
        alpha = input(1);
        beta = input(2);
        delta = input(3);
        
        y_prob = alpha .* exp(beta .* x_values) + delta;
        
        % Log likelihood function
        LL = y_observed .* log(y_prob) + (total_observations - y_observed) .* log(1 - y_prob);
        likelihood = -sum(LL);
    end
end

function [D50_x, D50_y] = calculate_D50(alpha, beta, delta)
    % Calculate 1/2 max of fit curve
    % x at midway point between y(0) and asymptote for exponential curve
    D50_x = log(0.5) / beta;
    % y at that value of x
    D50_y = alpha * exp(beta * D50_x) + delta;
end
