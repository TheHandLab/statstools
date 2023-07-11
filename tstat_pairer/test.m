%% load some data to run the tstat_pairer function
data=load("data.txt");                                                      % unpaired raw data, 24x2
reps=10;                                                                    % how many times to repeat this simulation?
perms=nan(size(data,1),reps);                                               % save the solutions here
figure(1);
hold on;
for rep=1:reps                                                              % how many times to run the function?
    t=2.701;                                                                % target t-value
    [bestfit,output,diagnostics]=tstat_pairer(data,[],t);                   % search for a solution
    disp(diagnostics.i);
    
    %% print some data
    plot(output(3,:),'k--');                                                % the simulation results
    plot([diagnostics.i,diagnostics.i],[t-0.05,t+0.05],'r-');               % show the endpoint
    perms(:,rep)=diagnostics.perm;                                          % to compare the solutions
end

%% format the figure
plot([1,diagnostics.i],[t,t],'r-');
xlabel('Iteration');
ylabel('Simulated t-value');
title('A discrete approach to finding matches between unpaired data');
a=axis;
text(a(2).*.90,t.*1.01,'target t-score','Color','r');