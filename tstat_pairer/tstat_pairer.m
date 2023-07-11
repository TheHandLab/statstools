%% given 2 columns of (unpaired but) repeated-measures data, find a combination that matches the reported mean differences, t-test and/or p-value, by stepwise approximations
%% unfortunately, this provides a tool to create real-looking effects from real-looking (but random) data...
%% INPUTS:
%% data unpaired list of datapoints (eg extracted from a graph)
%% se   reported standard error of *difference* between columns - one of 3 possible critical inputs
%% t    reported t-statistic (with N-1 degrees of freedom), for *difference* between columns - one of 3 possible critical inputs; must be POSITIVE
%% p    reported TWO-TAILED p-value for *difference* between columns - one of 3 possible critical inputs
%% s    how many swaps to do? default=10000
function [bestfit,output,diagnostics]=tstat_pairer(data,se,t,p,s)
    if isempty(data)
        error("No data provided");
    end
    if nargin==1
        error("Target argument empty: need SE, t, or p");
    elseif nargin==2
        if isempty(se)
            error("Target argument empty: need SE, t, or p");
        end
    elseif nargin==3
        if isempty(se) && isempty(t)
            error("Target argument empty: need SE, t, or p");
        end
    elseif nargin==4
        if isempty(se) && isempty(t) && isempty(p)
            error("Target argument empty: need SE, t, or p");
        end
        s=10000;
    end
    if ~exist('s','var') 
        s=10000;
    end
    if ~isempty(se) && isempty(t)
        t=(mean(data(:,1))-mean(data(:,2)))./se;                            % get t-value from mean difference and known se of the difference
    end
    if exist('p','var')
        if p<0
            error("P-value (2-tailed) cannot be below zero");
        elseif p>1
            error("P-value (2-tailed) cannot be above one");
        elseif isempty(t)
            % convert p to t
            t=tinv(1-(p./2),size(data,1)-1);                                % get target t-value from p value
        end
    end

    %% initialise the output variables
    bestfit=data;                                                           % the data to analyse
    bestfit(:,2)=data(randperm(size(data,1)),2);                            % randomise condition 2 data
    D=bestfit(:,1)-bestfit(:,2);                                            % current best estimate of differences between conditions
    [~,~,N,~,~,T,~]=describe(D);                                            % current best-fit statistics
    ix=[1:N]';                                                              % index for selecting data
    i=0;                                                                    % keep track of how many swaps so far
    output=nan(4,s);                                                        % for storing the simulation data
    pgap=1;                                                                 % gap between actual and target statistic, as a proportion (can be>1, stops at <.005)
    while i<s && pgap>.0001                                                 % keep checking until maximum iterations, or we're within 0.01% of the target

        %% serial swap approach____________________________________________
        for w=1:N                                                           % for each value of condition 1
            skip=0;                                                         % reset the skip value
            r=0;                                                            % reset condition 2 number

            % find first datapoint in condition 2 which moves the sample closer to the target, and swap it over
            while skip==0 && r<N                                            % until we have found a good swap, or until the end of condition 1
                r=r+1;                                                      % try the next value in condition 2
                if r~=w                                                     % can't swap a number with itself
                    i=i+1;                                                  % keep track of swaps (limit the overall simulation to eg 10000)
                    gap=abs(T-t);                                           % how far are is the current T value from the target t value?
                    pgap=gap./abs(t);                                       % gap as a proportion of the target score
                    ix2=ix;                                                 % new index to swap data around
                    ix2([w,r])=ix([r,w]);                                   % swap the data around
                    D=bestfit(:,1)-bestfit(ix2,2);                          % current best estimate of differences between conditions
                    [~,~,~,~,~,T2,~]=describe(D);                           % current best-fit statistics
                    
                    % if new gap is smaller than old gap, swap in the new order
                    if abs(T2-t)<gap
                        bestfit(:,2)=bestfit(ix2,2);                        % replace the best-fit data with the current pairing
                        skip=1;                                             % toggles the end of this loop
                    else
                        skip=0;                                             % keep trying the next value of condition 2
                    end

                    %% save some summary statistics________________________
                    D=bestfit(:,1)-bestfit(:,2);                            % current best estimate of differences between conditions
                    [M,~,N,SE,~,T,P]=describe(D);                           % current best-fit statistics
                    output(1,i)=M;                                          % current best estimate of mean difference between conditions
                    output(2,i)=SE;                                         % current best estimate of SE between conditions
                    output(3,i)=T;                                          % current best estimate of T-value for difference between conditions
                    output(4,i)=P;                                          % current best estimate of TWO-TAILED P-values for difference between conditions
                end                                                         % end of swap check loop
            end                                                             % end of while loop for checking
        end                                                                 % end of condition 1 loop
    end                                                                     % end of maximum iterations while loop
    output=output(:,1:i);                                                   % reduce output size
    diagnostics.i=i;                                                        % iterations needed to converge
    diagnostics.gap=gap;                                                    % final gap between simulated and target t-value
    diagnostics.pgap=pgap;                                                  % final gap as a proportion of the target t-value
    R=corrcoef(bestfit);                                                    % current best-fit correlation between columns
    diagnostics.r=R(1,2);                                                   % the correlation coefficient

    % reconstruct the order of bestfit data
    for n=1:N
        diagnostics.perm(n)=find(bestfit(n,2)==data(:,2),1);                % find first matching datapoint, save the position
    end
end                                                                         % end of function