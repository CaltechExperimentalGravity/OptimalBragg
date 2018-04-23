% run doAlGas.m many times to tune the hyper-parameters

selfies   = [0.5 1 1.5];
socialies = [0.5 1 1.5];
N = 3; % how many times to repeat with each setting

scores = zeros(length(selfies), length(socialies), N, 5);

for ii = 1:length(selfies)
    for jj = 1:length(socialies)
        for kk = 1:N
            selfie = selfies(ii);
            socialie = socialies(jj);
            doAlGaAs
            scores(ii, jj, kk, :) = [TNout.yy t_toc];
        end
    end
end

tnowstr = datestr(now, 'yymmdd_HHMM');

%% SAVE layer structure, IFOmodel, and noise
savename = ['Data/' NUMTOOLS.opt_name '_tune_' tnowstr];
save(savename, 'scores', 'selfies', 'socialies');
