function benchmark(m,p)
% BENCHMARK would test the performance of kernels

ntrail = 10;

p.use_gpu = exist('grandn','builtin');

fprintf('Benchmark for Motion Pattern Separating Model @ %s\n',...
    datestr(now));
% Load data into one animation
fprintf('Loading Data ... ');
anim = readdata(m,p,'chunks', ...
    mod(randperm(p.num_chunks*ntrail),p.num_chunks)+1);
fprintf('Done\n');

signal = norm(anim,'fro'); % Square Root of Signal

% First-Layer Test
z = anim2z(anim,m,p);
recover = z2anim(z,m);
snr = 20*log10(signal/norm(recover-anim,'fro'));
fprintf('SNR(Recover Image from 1st-Layer Codes) >> %.2f dB\n',snr);

% Suppress Display
p.ampmodel.minFunc_ind_Opts.Display = 0;
p.phasetrans.minFunc_ind_Opts.Display = 0;

% Second-Layer Test
[pcode,mcode] = z2pm(z,m,p);
zrec = pm2z(pcode,mcode,m);
snr = 20*log10(norm(z,'fro')/norm(zrec-z,'fro'));
fprintf('SNR(Recover 1st-Layer Codes from 2nd-Layer Codes) >> %.2f dB\n',snr);

% Combine Test
recover = z2anim(zrec,m);
snr = 20*log10(signal/norm(recover-anim,'fro'));
fprintf('SNR(Recover Image from 2nd-Layer Codes) >> %.2f dB\n',snr);

end