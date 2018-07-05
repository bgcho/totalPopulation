% Estended Kalman filter test with synthetic data
calTSFast = freqRespLib('calTSFast');
EKF = freqRespLib('EKF');

% Generate synthetic data (Attenuated)
waterDepth = 120;
r = linspace(0,5000,101);
dr = mean(diff(r));
z1 = 50 - 0.005*r;
H1 = 50*ones(1,length(r));
z2 = 110*ones(1,length(r));
H2 = 10*ones(1,length(r));
ratioLayer = 0.5*ones(1,length(r));
zNb = 10*ones(1,length(r));
nAdB = -((r-2500)./1000).^2 + 10*log10(3);
% nA = 10.^((-((r-2500)./1000).^2 + 10*log10(3))./10);

nssh = matfile('nsshFlPdf.mat');
fishStruct.forkLength = nssh.fl;
fishStruct.fracMajor = 0.33;
fishStruct.wPower = 3.35;
fishStruct.wCoeff = 0.00335;
fishStruct.totalLengthExpFactor = 1.103;
fishStruct.totalLengthOffset = 0.01;
fishStruct.fleshViscosity = 50;
fishStruct.fleshDensity = 1071;
fishStruct.species = 'Herring';
fishStruct.meanForkLength = nssh.meanFl; % 24.2;
fishStruct.stdForkLength = nssh.stdFl; % fishStruct.meanForkLength*0.068;
fishStruct.flpdf = nssh.flPdf; % normpdf(fishStruct.forkLength,fishStruct.meanForkLength,fishStruct.stdForkLength);


freqHz = [850 955 1125 1335 1465 1600];
numFreq = length(freqHz);
c = 1500;
k = 2*pi*freqHz./c;

for ii=1:length(r)
	TS = calTSFast(fishStruct, z1(ii), H1(ii), z2(ii), H2(ii), ratioLayer(ii), zNb(ii), freqHz);
	scsTot(:,ii) = 4*pi./k.^2.*imag(TS.meanSFn);
	ssTrue(:,ii) = TS.meanTS + nAdB(ii);
	% ssTrue(:,ii) = TS.meanTS + 10*log10(nA(ii));
end


nV = 10.^(nAdB./10)./waterDepth;
% nV = nA./waterDepth;
nVMat = repmat(nV,[numFreq, 1]);

cumAttnTmp = 10*log10(exp(1))*dr*cumsum(nVMat.*scsTot,2);
cumAttn = [zeros(numFreq,1) cumAttnTmp(:,1:end-1)];

ssAttn = ssTrue - cumAttn;

nPts = 10;
gnoise = 0+3*randn(numFreq, length(r), nPts);
ssSynthetic = repmat(ssAttn,[1,1,nPts]) + gnoise;


stateInit = [z1(1), H1(1), z2(1), H2(1), ratioLayer(1), zNb(1), nAdB(1)]';
% stateInit = [z1(1), H1(1), z2(1), H2(1), ratioLayer(1), zNb(1), nA(1)]';
numStateMem = length(stateInit);
stateCovInit = zeros(numStateMem, numStateMem);	% 1e-3*diag([1 1 1 1 1 0.1 10], 0)
stateCov = 1e-3*diag([1 1 1 1 1 0.1 10], 0); % zeros(numStateMem, numStateMem); % eye(numStateMem, numStateMem);
% measurementCov = repmat(2e-4*eye(numFreq,numFreq), [nPts, nPts]); % kron(eye(nPts), 1e-3*eye(numFreq,numFreq));	% 2e-4, 5e-4
measurementCov = diag(0.1*rand(1,nPts*numFreq),0);% Use this configuration for 3dB noise
% measurementCov = diag(1e-6*rand(1,nPts*numFreq),0);% Use this configuration for no noise
correctionInit = zeros(numFreq,1);


% 
for ii=1:length(r)
	measurement = squeeze(ssSynthetic(:,ii,:));
	% measurement = measurementTmp(:);
	[stateFinal, stateCovFinal, correctionFinal] = EKF(fishStruct, stateInit, stateCovInit, stateCov, measurement, measurementCov, correctionInit, dr, waterDepth, freqHz);

	stateInit = stateFinal;
	stateCovInit = stateCovFinal;
	correctionInit = correctionFinal;
	stateAll(:,ii) = stateFinal;
	correctionAll(:,ii) = correctionFinal;
	stateCovAll(:,:,ii) = stateCovFinal;
end




% For analysis
TSForAnlysis(1) = calTSFast(fishStruct, z1(1), H1(1), z2(1), H2(1), ratioLayer(1), zNb(1), logspace(1,4,1000));
TSForAnlysis(2) = calTSFast(fishStruct, z1(end), H1(end), z2(end), H2(end), ratioLayer(end), zNb(end), logspace(1,4,1000));



