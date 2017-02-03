%% ETMHIST
%  loads a text file with the layer thicknesses in units of lambda
%  and makes a histogram of the reflectivities after varying the
%  parameters
%
%  Rana  Dec 16, 2009
%

z = load('etm_layers_091216_172635.txt');

%% Run the trials
N = 10000;   % number of trials
ts = zeros(N,2);  % empty array with 2 columns for the 2 wavelengths

for jj = 1:N
    
    % perturb the layer thicknesses with a 
    x0 = z + 0.002 * randn(size(z));
    
    % calculate the error and the Transmission
    [y, T] = optETM(x0,0);

    ts(jj,:) = T;
    
end

%% Do plots and histograms
subplot(211)
[n532, x532] = hist(ts(:,1),100);
bar(x532,n532,'g')
set(gca,'Yscale','log')
xlabel('T @ 532 nm')
text(0.055,300,['Mean = ' num2str(mean(ts(:,1)))])
text(0.055,180,['Std Dev = ' num2str(std(ts(:,1)))])

subplot(212)
[n1064, x1064] = hist(ts(:,2),100);
bar(x1064,n1064,'r')
set(gca,'Yscale','log')
xlabel('T @ 1064 nm')
text(20.5e-6,300,['Mean = ' num2str(mean(ts(:,2)))])
text(20.5e-6,180,['Std Dev = ' num2str(std(ts(:,2)))])

gak = [x532' n532'];
save histdata1.txt gak -ascii -double -tabs

gak = [x1064' n1064'];
save histdata2.txt gak -ascii -double -tabs