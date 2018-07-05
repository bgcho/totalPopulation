function fHandle = freqRespLib(fName)
	switch lower(fName)
	case 'avgtl'
		fHandle = @avgTL;
	case 'gwindow'
		fHandle = @gWindow;
	case 'calscs'
		fHandle = @calSCS;
	case 'calscsgeneral'
		fHandle = @calSCSGeneral;
	case 'calpdf'
		fHandle = @calPdf;
	case 'calpdfgradient'
		fHandle = @calPdfGradient;
	case 'calgradient'
		fHandle = @calGradient;
	case 'gradientdescent'
		fHandle = @gradientDescent;
	case 'createslmap'
		fHandle = @createSlMap;
	case 'gettstampfrombimap'
		fHandle = @getTStampFromBimap;
	case 'caltsfast'
		fHandle = @calTSFast;
	case 'caljacobiants'
		fHandle = @calJacobianTS;
	case 'caljacobiants2'
		fHandle = @calJacobianTS2;
	case 'ekf'
		fHandle = @EKF;
	case 'ekf2'
		fHandle = @EKF2;
    case 'avgss'
        fHandle = @avgSS;
    case 'avgareamap'
        fHandle = @avgAreaMap;
	end

	function tlMapAvg = avgTL(tlName, xgrid, ygrid, depthOption)
		for oo=1:length(tlName)			
			tlStruct = matfile(tlName{oo});
			switch lower(depthOption)
			case 'bottom'
				tlMap = tlStruct.TL_twoway_bott_dB_final;
			case 'fish'
				tlMap = tlStruct.TL_twoway_fish_dB_final;
			case '030_140'
				tlMap = tlStruct.TL_twoway_030_140;
			end
			% tlMap = tlStruct.TL_twoway_fish_dB_final;
			% areaMap = tlStruct.area_term_dB;		
			% areaMap = 10*log10(read_arr(areaName{oo}));					
			% ssMapBig = splMap - tlMap - areaMap;
			tlMapBig(:,:,oo) = 10*log10(interp2(tlStruct.datax, (tlStruct.datay)', 10.^(tlMap./10), xgrid, ygrid(:)));
		end
		tlMapAvg = 10*log10(squeeze(nanmean(10.^(tlMapBig./10),3)));		
	end

	function gWindow = gWindow(R,stdWin)
		gWin = 1./(sqrt(2*pi)*stdWin)*exp(-R.^2./(2*stdWin^2));
		normFactor = nansum(gWin(:));
		gWindow = gWin./normFactor;
	end

	function scs = calSCS(fishStruct, zOc, zNb, freqHz)
    % zOc and zNb should be the same size vectores, in meters
    % sigmaBs: numDepth*numLength*numFreq, m^2
    % Typical example of fishStruct:
    % fishStruct.species='Herring';
    % .fleshDensity = 1071;
    % .fleshViscosity = 50;
    % .forkLength = [15:50];
    % .totalLengthOffset = 0.01;
    % .totalLengthExpFactor = 1.103;
    % .fracMajor = 0.33;
    % .wCoeff = 0.00335;
    % .wPower = 3.35;

  	% nssh = matfile('nsshFlPdf.mat');
		% fishStruct.forkLength = nssh.fl;
		% fishStruct.fracMajor = 0.33;
		% fishStruct.wPower = 3.35;
		% fishStruct.wCoeff = 0.00335;
		% fishStruct.totalLengthExpFactor = 1.103;
		% fishStruct.totalLengthOffset = 0.01;
		% fishStruct.fleshViscosity = 50;
		% fishStruct.fleshDensity = 1071;
		% fishStruct.species = 'Herring';
		% fishStruct.meanForkLength = nssh.meanFl; % 24.2;
		% fishStruct.stdForkLength = nssh.stdFl; % fishStruct.meanForkLength*0.068;
		% fishStruct.flpdf = nssh.flPdf; % normpdf(fishStruct.forkLength,fishStruct.meanForkLength,fishStruct.stdForkLength);



    if length(zOc)>1 
      dz = mean(diff(zOc));
    else
      dz = 1;
    end

    if length(zNb)>1 
      dzNb = mean(diff(zNb));
    else
      dzNb = 1;
    end



    % if nargin>4
    %   zpdf = interp1(varargin{1}(:,1), varargin{1}(:,2), zOc(:));
    %   normFactor = sum(zpdf)*dz;
    %   zpdf = [zOc(:) zpdf(:)./normFactor];
    % else
    %   zpdf = ones(length(zOc),1)./(dz*length(zOc));
    %   zpdf = [zOc(:) zpdf(:)];
    % end
    
    name = fishStruct.species;
    dFish = fishStruct.fleshDensity;  %kg/m^3
    visFish = fishStruct.fleshViscosity; % Pa*sec
    lf = fishStruct.forkLength; % cm
    lt = fishStruct.totalLengthOffset + fishStruct.totalLengthExpFactor*lf; % cm
    b = (fishStruct.fracMajor*lt)/2;  % cm
    weight = fishStruct.wCoeff*lf.^(fishStruct.wPower); % g

    numLength = length(lf);
    numDepth = length(zOc);
    numFreq = length(freqHz);
    numNBDepth = length(zNb);

    V0Tmp = 0.05*(zNb(:)/10+1)*weight(:)'; % zNb*lf, cm^3
    V0 = repmat(reshape(V0Tmp,[1,numNBDepth,numLength]), [numDepth,1,1]);
    

    zOcMat = repmat(zOc(:), [1, numNBDepth, numLength]);
    V = V0./(1+zOcMat./10);
    r = (3/(4*pi)*V).^(1/3)/100;  % zOc*zNb*lf, m    

    bMat = repmat(reshape(b(:), [1, 1, numLength]), [numDepth, numNBDepth, 1]);
    a = sqrt(3/(4*pi)*1./(1+zOcMat./10).*(V0./bMat)); % zOc*zNb*lf, cm



    ecc = a./bMat;  %zOc*zNb*lf

    
    zeta = sqrt(2)*(1-ecc.^2).^(1/4)./ecc.^(1/3)./sqrt(log((1+sqrt(1-ecc.^2))./(1-sqrt(1-ecc.^2))));  % zOc*zNb*lf

    pAmb = 10^5;  % Ambient pressure: 10^5 Pa
    ratioHeat = 1.4;  % Specific heat ratio    

    f0 = 1/(2*pi)*zeta./r.*sqrt(3*ratioHeat*pAmb*(1+zOcMat./10)/dFish); % zOc*zNb*lf, Hz

    rMat = repmat(reshape(r,[numDepth,numNBDepth,numLength,1]), [1,1,1,numFreq]);	% zOc*zNb*lf*freq, m
    freqMat = repmat(reshape(freqHz(:),[1,1,1,numFreq]), [numDepth,numNBDepth,numLength,1]);	% zOc*zNb*lf*freq, Hz
    f0Mat = repmat(reshape(f0,[numDepth,numNBDepth,numLength,1]),[1,1,1,numFreq]);	% zOc*zNb*lf*freq, Hz
    c = 1500; % sound speed in water, m/s
    kMat = 2*pi*freqMat./c;	% zOc*zNb*lf*freq, 1/m
    
    hInv = (2*pi*rMat.*freqMat.^2)./(c*f0Mat) + visFish./(pi*dFish*rMat.^2.*f0Mat);
    fleshDampFactor = 2*visFish./(dFish*rMat.^2*2*pi.*freqMat);    
    % radDampFactor = totDampFactor - fleshDampFactor;   
    radDampFactor = kMat.*rMat;
    % totDampFactor = hInv.*f0Mat./freqMat;
    totDampFactor = radDampFactor + fleshDampFactor;


    scs.sigmaTot = (4*pi*rMat.^2)./(((f0Mat./freqMat).^2-1).^2+((f0Mat./freqMat).*hInv).^2);
    scs.sigmaBs = scs.sigmaTot./(4*pi);
    scs.f0 = f0;

    scs.sFnRad = ((f0Mat./freqMat).^2-1).*kMat.*rMat./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2) ...
              + 1i*radDampFactor.*(kMat.*rMat)./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    % scs.sFnFlesh = 1i*fleshDampFactor.*(kMat.*rMat).^2./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    scs.sFnFlesh = 1i*fleshDampFactor.*(kMat.*rMat)./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    scs.sFn = scs.sFnRad + scs.sFnFlesh;

    pdfMatFl = repmat(reshape(fishStruct.flpdf, [1,1,numLength,1]),[numDepth,numNBDepth,1,numFreq]);
    % pdfMatZ = repmat(reshape(zpdf(:,2), [numDepth,1,1,1]),[1,numNBDepth,numLength,numFreq]);
    dlf = mean(abs(diff(lf)));
    
    scs.sigmaBSMeanOverLf = squeeze(nansum(scs.sigmaBs.*pdfMatFl,3)*dlf);	% numDepth*numNBDepth*numFreq

    if isfield(fishStruct,'zocZnbdPdf')
    	pdfMatZ = repmat(reshape(fishStruct.zOCzNBDepthPdf,[numDepth,numNBDepth,1,1]),[1,1,numLength,numFreq]);
    	scs.sigmaBSMeanOverAll = squeeze(nansum(nansum(nansum(scs.sigmaBS.*pdfMatFl.*pdfMatZ,1),2),3)*dlf);
    end
    
    
    % sFnMat = sqrt(scs.sigmaBs).*freqMat./c;
    % sFnMat = scs.sFnRad;
    % scs.meanSFn = squeeze(nansum(nansum(sFnMat.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.meanSFnFlesh = squeeze(nansum(nansum(scs.sFnFlesh.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.meanSqrSFn = squeeze(nansum(nansum(abs(sFnMat).^2.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.meanTS = 10*log10(squeeze(nansum(nansum(scs.sigmaBs.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz));
    % TSmat = 10*log10(scs.sigmaBs);
    % TSmeandB = squeeze(trapz(trapz(TSmat.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % TSsqrdB = squeeze(trapz(trapz(TSmat.^2.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    
    % scs.meanFleshDampFactor = squeeze(nansum(nansum(fleshDampFactor.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.meanRadDampFactor = squeeze(nansum(nansum(radDampFactor.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.meanTotDampFactor = squeeze(nansum(nansum(totDampFactor.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.stdTS = sqrt(TSsqrdB - TSmeandB.^2);
    scs.freqHz = freqHz;
    scs.depth = zOc;
    scs.nbDepth = zNb;
    scs.lf = lf;

    % k = 2*pi*freqHz/c;
    % kMat = repmat(reshape(k, [1,1,numFreq]),[numDepth,numLength,1]);
    % rMat = repmat(reshape(r, [numDepth,numLength,1]),[1,1,numFreq]);
    % sFnPrSphere = -kMat.*rMat + 1i*(kMat.*rMat).^2;
    % scs.meanSFnPrSphere = squeeze(nansum(nansum(sFnPrSphere.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % scs.meanTSPrSphere = 10*log10(squeeze(nansum(nansum(abs(sFnPrSphere./kMat).^2.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz));
    % scs.meanRadius = (squeeze(nansum(nansum((r.^3).*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz)).^(1/3);
    % scs.radiusEq = r;
    % scs.meanVolume = squeeze(nansum(nansum(V.*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz);
    % scs.stdVolume = sqrt(nansum(nansum((V-scs.meanVolume).^2.*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz);
    % scs.stdRadius = sqrt(nansum(nansum((r-scs.meanRadius).^2.*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz);
    % scs.totDampFactor = totDampFactor;

        
  end

  function scs = calSCSGeneral(fishStruct, zOc, zNb, freqHz, varargin)
    % zOc and zNb should be the same size vectores, in meters
    % sigmaBs: numDepth*numLength*numFreq, m^2
    % Typical example of fishStruct:
    % fishStruct.species='Herring';
    % .fleshDensity = 1071;
    % .fleshViscosity = 50;
    % .forkLength = [15:50];
    % .totalLengthOffset = 0.01;
    % .totalLengthExpFactor = 1.103;
    % .fracMajor = 0.33;
    % .wCoeff = 0.00335;
    % .wPower = 3.35;

    if length(zOc)>1 
      dz = mean(diff(zOc));
    else
      dz = 1;
    end

    if nargin>4
      zpdf = interp1(varargin{1}(:,1), varargin{1}(:,2), zOc(:), 'linear',0);
      normFactor = sum(zpdf)*dz;
      zpdf = [zOc(:) zpdf(:)./normFactor];
    else
      zpdf = ones(length(zOc),1)./(dz*length(zOc));
      zpdf = [zOc(:) zpdf(:)];
    end
    
    name = fishStruct.species;
    dFish = fishStruct.fleshDensity;  %kg/m^3
    visFish = fishStruct.fleshViscosity; % Pa*sec
    lf = fishStruct.forkLength; % cm
    lt = fishStruct.totalLengthOffset + fishStruct.totalLengthExpFactor*lf; % cm
    b = (fishStruct.fracMajor*lt)/2;  % cm
    weight = fishStruct.wCoeff*lf.^(fishStruct.wPower); % g
    V0 = 0.05*(zNb(:)/10+1)*weight(:)'; % zNb*lf, cm^3

    numLength = length(lf);
    numDepth = length(zOc);
    numFreq = length(freqHz);

    zOcMat = repmat(zOc(:), [1,numLength]);
    V = V0./(1+zOcMat./10);
    r = (3/(4*pi)*V).^(1/3)/100;  % zNb*lf, m    

    bMat = repmat(b(:)', [numDepth,1]);
    a = sqrt(3/(4*pi)*1./(1+zOcMat./10).*(V0./bMat)); % zNb*lf, cm
    ecc = a./bMat;  %zNb*lf
    zeta = sqrt(2)*(1-ecc.^2).^(1/4)./ecc.^(1/3)./sqrt(log((1+sqrt(1-ecc.^2))./(1-sqrt(1-ecc.^2))));  % zNb*lf

    pAmb = 10^5;  % Ambient pressure: 10^5 Pa
    ratioHeat = 1.4;  % Specific heat ratio    

    f0 = 1/(2*pi)*zeta./r.*sqrt(3*ratioHeat*pAmb*(1+zOcMat./10)/dFish); % zNb*lf, Hz

    rMat = repmat(reshape(r,[numDepth,numLength,1]), [1,1,numFreq]);
    freqMat = repmat(reshape(freqHz(:),[1,1,numFreq]), [numDepth,numLength,1]);
    f0Mat = repmat(reshape(f0,[numDepth,numLength,1]),[1,1,numFreq]);
    c = 1500; % sound speed in water, m/s
    kMat = freqMat./c;

    hInv = (2*pi*rMat.*freqMat.^2)./(c*f0Mat) + visFish./(pi*dFish*rMat.^2.*f0Mat);
    fleshDampFactor = 2*visFish./(dFish*rMat.^2*2*pi.*freqMat);    
    % radDampFactor = totDampFactor - fleshDampFactor;   
    radDampFactor = kMat.*rMat;
    % totDampFactor = hInv.*f0Mat./freqMat;
    totDampFactor = radDampFactor + fleshDampFactor;


    scs.sigmaTot = (4*pi*rMat.^2)./(((f0Mat./freqMat).^2-1).^2+((f0Mat./freqMat).*hInv).^2);
    scs.sigmaBs = scs.sigmaTot./(4*pi);
    scs.f0 = f0;

    scs.sFnRad = ((f0Mat./freqMat).^2-1).*kMat.*rMat./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2) ...
              + 1i*radDampFactor.*(kMat.*rMat)./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    % scs.sFnFlesh = 1i*fleshDampFactor.*(kMat.*rMat).^2./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    scs.sFnFlesh = 1i*fleshDampFactor.*(kMat.*rMat)./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    scs.sFn = scs.sFnRad + scs.sFnFlesh;

    pdfMatFl = repmat(reshape(fishStruct.flpdf, [1,numLength,1]),[numDepth,1,numFreq]);
    pdfMatZ = repmat(reshape(zpdf(:,2), [numDepth,1,1]),[1,numLength,numFreq]);
    
    dlf = mean(abs(diff(lf)));
    % sFnMat = sqrt(scs.sigmaBs).*freqMat./c;
    sFnMat = scs.sFnRad;
    scs.meanSFn = squeeze(nansum(nansum(sFnMat.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.meanSFnFlesh = squeeze(nansum(nansum(scs.sFnFlesh.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.meanSqrSFn = squeeze(nansum(nansum(abs(sFnMat).^2.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.meanTS = 10*log10(squeeze(nansum(nansum(scs.sigmaBs.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz));
    TSmat = 10*log10(scs.sigmaBs);
    TSmeandB = squeeze(trapz(trapz(TSmat.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    % TSsqrdB = squeeze(trapz(trapz(TSmat.^2.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    TSsqrdB = squeeze(trapz(squeeze((trapz(TSmat.*pdfMatFl,2)*dlf).^2).*repmat(zpdf(:,2),[1,numFreq]),1)*dz);
    TSsqrdB = TSsqrdB(:);
  
    
    scs.meanFleshDampFactor = squeeze(nansum(nansum(fleshDampFactor.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.meanRadDampFactor = squeeze(nansum(nansum(radDampFactor.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.meanTotDampFactor = squeeze(nansum(nansum(totDampFactor.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.stdTS = sqrt(TSsqrdB - TSmeandB.^2);
    scs.freqHz = freqHz;
    scs.depth = zOc;
    scs.nbDepth = zNb;

    k = 2*pi*freqHz/c;
    kMat = repmat(reshape(k, [1,1,numFreq]),[numDepth,numLength,1]);
    rMat = repmat(reshape(r, [numDepth,numLength,1]),[1,1,numFreq]);
    sFnPrSphere = -kMat.*rMat + 1i*(kMat.*rMat).^2;
    scs.meanSFnPrSphere = squeeze(nansum(nansum(sFnPrSphere.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz);
    scs.meanTSPrSphere = 10*log10(squeeze(nansum(nansum(abs(sFnPrSphere./kMat).^2.*pdfMatFl.*pdfMatZ,1),2)*dlf*dz));
    scs.meanRadius = (squeeze(nansum(nansum((r.^3).*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz)).^(1/3);
    scs.radiusEq = r;
    scs.meanVolume = squeeze(nansum(nansum(V.*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz);
    scs.stdVolume = sqrt(nansum(nansum((V-scs.meanVolume).^2.*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz);
    scs.stdRadius = sqrt(nansum(nansum((r-scs.meanRadius).^2.*squeeze(pdfMatFl(:,:,1).*pdfMatZ(:,:,1)),1),2)*dlf*dz);
    scs.totDampFactor = totDampFactor;

        
  end



  function [zOc,zNb,zPdfMat] = calPdf(z, zOcBar, zOcStd, zNbBar, zNbStd, option)
		zOc = z;
		zNb = z;
		numDepth = length(z);
		dz = mean(abs(diff(z)));
		[zOcMat, zNbMat] = ndgrid(zOc, zNb);

		zOcPdf = normpdf(z, zOcBar, zOcStd);
		zOcPdf = zOcPdf./trapz(zOcPdf*dz);
		zOcPdfMat = repmat(zOcPdf(:), [1,numDepth]);

		zNbPdf = normpdf(z, zNbBar, zNbStd);
		zNbPdfMatTmp = repmat(zNbPdf(:)', [numDepth,1]);
		switch lower(option)
		case 'negtivelybuoyant'
			zNbPdfMatTmp(zNbMat>zOcMat) = 0;
		case 'noconstraint'
		end

		zNbPdfMat = zNbPdfMatTmp./repmat(trapz(zNbPdfMatTmp*dz,2),[1,numDepth]);

		zPdfMat = zOcPdfMat.*zNbPdfMat;	%zOc*zNb
		
	end

	
	function [zOc, zNb, dP_dZOcBar, dP_dZOcStd, dP_dZNbBar, dP_dZNbStd] = calPdfGradient(z, zOcBar, zOcStd, zNbBar, zNbStd, option)
		zOc = z;
		zNb = z;
		numDepth = length(z);
		dz = mean(abs(diff(z)));
		[zOcMat, zNbMat] = ndgrid(zOc, zNb);

		% P(zOc;zOcBar,zOcStd)
		PzOcVec = normpdf(z, zOcBar, zOcStd);
		% PzOcVec = PzOcVec./trapz(PzOcVec*dz);
		PzOcTmp = repmat(PzOcVec(:), [1,numDepth]);
		normFactorPzOc = repmat(trapz(PzOcTmp,1)*dz,[numDepth,1]);
		PzOc = PzOcTmp./normFactorPzOc;


		% P(zNb|zOc;zOcBar,zOcStd)
		PzNbVec = normpdf(z, zNbBar, zNbStd);
		PzNbTmp = repmat(PzNbVec(:)', [numDepth,1]);
		switch lower(option)
		case 'negtivelybuoyant'
			PzNbTmp(zNbMat>zOcMat) = 0;
		case 'noconstraint'
		end
		normFactorPzNb = repmat(trapz(PzNbTmp*dz,2),[1,numDepth]);
		PzNb = PzNbTmp./normFactorPzNb;

		% dP(zOc;zOcBar,zOcStd)/dzOcBar, dP(zOc;zOcBar,zOcStd)/dzOcStd
		dNzOc_dZOcBar = (zOcMat-zOcBar)./(sqrt(2*pi)*zOcStd^3).*exp(-(zOcMat-zOcBar).^2./(2*zOcStd^2));
		dNzOc_dZOcStd = ((zOcMat-zOcBar).^2 - zOcStd^2)./(sqrt(2*pi)*zOcStd^4).*exp(-(zOcMat-zOcBar).^2./(2*zOcStd^2));

		normFactor_dNzOc_dZOcBar = repmat(trapz(dNzOc_dZOcBar,1)*dz,[numDepth,1]);
		normFactor_dNzOc_dZOcStd = repmat(trapz(dNzOc_dZOcStd,1)*dz,[numDepth,1]);

		dPzOc_dZOcBar = dNzOc_dZOcBar./normFactorPzOc - PzOc./normFactorPzOc.*normFactor_dNzOc_dZOcBar;
		dPzOc_dZOcStd = dNzOc_dZOcStd./normFactorPzOc - PzOc./normFactorPzOc.*normFactor_dNzOc_dZOcStd;

		% dP(zNb|zOc;zNbBar,zNbStd)/dzNbBar, dP(zNb|zOc;zNbBar,zNbStd)/dzNbStd
		dNzNb_dZNbBar = (zNbMat-zNbBar)./(sqrt(2*pi)*zNbStd^3).*exp(-(zNbMat-zNbBar).^2./(2*zNbStd^2));
		switch lower(option)
		case 'negtivelybuoyant'
			dNzNb_dZNbBar(zNbMat>zOcMat) = 0;
		case 'noconstraint'
		end

		dNzNb_dZNbStd = ((zNbMat-zNbBar).^2 - zNbStd^2)./(sqrt(2*pi)*zNbStd^4).*exp(-(zNbMat-zNbBar).^2./(2*zNbStd^2));
		switch lower(option)
		case 'negtivelybuoyant'
			dNzNb_dZNbStd(zNbMat>zOcMat) = 0;
		case 'noconstraint'
		end
		

		normFactor_dNzNb_dZNbBar = repmat(trapz(dNzNb_dZNbBar,2)*dz,[1,numDepth]);
		normFactor_dNzNb_dZNbStd = repmat(trapz(dNzNb_dZNbStd,2)*dz,[1,numDepth]);

		dPzNb_dZNbBar = dNzNb_dZNbBar./normFactorPzNb - PzNb./normFactorPzNb.*normFactor_dNzNb_dZNbBar;
		dPzNb_dZNbStd = dNzNb_dZNbStd./normFactorPzNb - PzNb./normFactorPzNb.*normFactor_dNzNb_dZNbStd;

		% dP(zOc,zNb;zOcBar,zOcStd,zNbBar,zNbStd)/dzOcBar, dP(zOc,zNb;zOcBar,zOcStd,zNbBar,zNbStd)/dzOcStd
		dP_dZOcBar = dPzOc_dZOcBar.*PzNb;
		% dP_dZOcStd = dNzOc_dZOcStd.*PzNb;
		dP_dZOcStd = dPzOc_dZOcStd.*PzNb;

		% dP(zOc,zNb;zOcBar,zOcStd,zNbBar,zNbStd)/dzNbBar, dP(zOc,zNb;zOcBar,zOcStd,zNbBar,zNbStd)/dzNbStd
		dP_dZNbBar = PzOc.*dPzNb_dZNbBar;
		dP_dZNbStd = PzOc.*dPzNb_dZNbStd;
	end

	

	function [J, dJ_dTheta] = calGradient(ssData, windowVec, zOcBar, zOcStd, zNbBar, zNbStd, nAdB, z, fishStruct, freqHz, option)
		% ssData: numPixel*numFreq(5 or 6)
		% windowVec: numPixel*1(5 or 6)
		% theta = [zOcBar, zOcStd, zNbBar, zNbStd, nAdB]';
		zOc = z;
		zNb = z;
		dzOc = mean(abs(diff(zOc)));
		dzNb = mean(abs(diff(zNb)));

		numZOc = length(zOc);
		numZNb = length(zNb);
		numFreq = length(freqHz);
		numPixel = size(ssData,1);

		scs = calSCS(fishStruct, zOc, zNb, freqHz);
		[~, ~, P] = calPdf(z, zOcBar, zOcStd, zNbBar, zNbStd, option);
		[~, ~, dP_dZOcBar, dP_dZOcStd, dP_dZNbBar, dP_dZNbStd] = calPdfGradient(z, zOcBar, zOcStd, zNbBar, zNbStd, option);

		% numFreq*1
		scsBack = squeeze(trapz(trapz(scs.sigmaBSMeanOverLf.*repmat(P,[1,1,numFreq]),1),2))*dzOc*dzNb;
		TS = 10*log10(scsBack);
		dTS_dTheta1 = 10*log10(exp(1))./scsBack.*squeeze(trapz(trapz(scs.sigmaBSMeanOverLf.*repmat(dP_dZOcBar,[1,1,numFreq]),1),2))*dzOc*dzNb;
		dTS_dTheta2 = 10*log10(exp(1))./scsBack.*squeeze(trapz(trapz(scs.sigmaBSMeanOverLf.*repmat(dP_dZOcStd,[1,1,numFreq]),1),2))*dzOc*dzNb;
		dTS_dTheta3 = 10*log10(exp(1))./scsBack.*squeeze(trapz(trapz(scs.sigmaBSMeanOverLf.*repmat(dP_dZNbBar,[1,1,numFreq]),1),2))*dzOc*dzNb;
		dTS_dTheta4 = 10*log10(exp(1))./scsBack.*squeeze(trapz(trapz(scs.sigmaBSMeanOverLf.*repmat(dP_dZNbStd,[1,1,numFreq]),1),2))*dzOc*dzNb;

		TSMat = repmat(TS(:)', [numPixel,1]);
		dTS_dThetaMat1 = repmat(dTS_dTheta1(:)', [numPixel,1]);
		dTS_dThetaMat2 = repmat(dTS_dTheta2(:)', [numPixel,1]);
		dTS_dThetaMat3 = repmat(dTS_dTheta3(:)', [numPixel,1]);
		dTS_dThetaMat4 = repmat(dTS_dTheta4(:)', [numPixel,1]);
		% nAdBMat = repmat(nAdBMat, [numPixel,numFreq]);
		windowMat = repmat(windowVec(:), [1,numFreq]);

		% numFreq*1
		ssStd = nanstd(ssData,0,1);
		ssStdMat = repmat(ssStd(:)', [numPixel,1]);

		dJ_dTheta1 = -nansum(nansum(windowMat./ssStdMat.^2.*(ssData - nAdB - TSMat).*dTS_dThetaMat1,2),1);
		dJ_dTheta2 = -nansum(nansum(windowMat./ssStdMat.^2.*(ssData - nAdB - TSMat).*dTS_dThetaMat2,2),1);
		dJ_dTheta3 = -nansum(nansum(windowMat./ssStdMat.^2.*(ssData - nAdB - TSMat).*dTS_dThetaMat3,2),1);
		dJ_dTheta4 = -nansum(nansum(windowMat./ssStdMat.^2.*(ssData - nAdB - TSMat).*dTS_dThetaMat4,2),1);
		dJ_dTheta5 = -nansum(nansum(windowMat./ssStdMat.^2.*(ssData - nAdB - TSMat),2),1);		

		dJ_dTheta = [dJ_dTheta1 ; dJ_dTheta2 ; dJ_dTheta3 ; dJ_dTheta4 ; dJ_dTheta5];
		J = nansum(nansum(windowMat./ssStdMat.^2.*(ssData - nAdB - TSMat).^2,2),1);		

	end

	function [theta,costVal] = gradientDescent(ssData, windowVec, zOcBar, zOcStd, zNbBar, zNbStd, nAdB, z, fishStruct, freqHz, learnRate, numIter, option, threshold, doPlot)

		theta = [zOcBar zOcStd zNbBar zNbStd nAdB]';
		scs = calSCS(fishStruct, z, z, freqHz);
		numFreq = length(freqHz);
		dz = mean(abs(diff(z)));

		if doPlot
			h = figure;
		end

		for ii=1:numIter
			[J(ii), dJ_dTheta(:,ii)] = calGradient(ssData, windowVec, theta(1), theta(2), theta(3), theta(4), theta(5), z, fishStruct, freqHz, option);
			theta = theta - learnRate*dJ_dTheta(:,ii);
			display([sprintf('%d th iteration: ', ii), num2str(theta')])
			if doPlot
				figure(h)
				subplot(2,1,1); hold on
				plot(ii,J(ii),'ro')
				set(gca,'YScale','log')
				

				[~,~,probDensity] = calPdf(z, theta(1), theta(2), theta(3), theta(4), option);
				TS = squeeze(10*log10(trapz(trapz(scs.sigmaBSMeanOverLf.*repmat(probDensity,[1,1,numFreq]),1),2)*dz^2));
				subplot(2,1,2);
				errorbar(freqHz, 10*log10(nanmean(10.^(ssData./10))), nanstd(ssData),'linewidth',2, 'color','blue')
				hold on
				plot(freqHz, TS + theta(5), '-^','color','red','linewidth',2)
				hold off

				drawnow
			end
			if ii>1
				JDiff = J(ii) - J(ii-1);
			else
				JDiff = J(ii);
			end

			if abs(JDiff/J(ii)) < threshold
				break
			end

		end

		costVal = J(end);

	end

	% function background = getBackgroundLevelNearSource(frequency,ssMapBulk,xgrid,ygrid,tStamp,numMix)
	function [xgrid,ygrid,tStamp,slDiffMap] = createSlMap(frequency, xgrid, ygrid, tStamp, bimapFile);
		% frequency = 955;
		freqVec = [850 955 1125 1335 1465 1600];
		freqIndx = find(freqVec==frequency);
		% ssFileName = fullfile(pwd,'..','ssMapsFullArray',sprintf('%04d',frequency),sprintf('ssMapBulk%dPingAvg.mat',nAvg));
		% load(ssFileName,'xgrid','ygrid','tStamp');
		% ssGridFileName = fullfile(pwd,'..','ssMapsFullArray','0955WideArea','gridSS.mat');
		% load(ssGridFileName,'xgrid','ygrid')
		% splDir = fullfile(pwd,'..','..','..','..','geoclutter','home','bgcho','splMapFullArray','feb20',sprintf('ah_t6a_wt1s_%d',frequency));
		% fileAll = dir(splDir);
		% slFileName = fullfile(pwd,'..','slMaps','slAzDependence.mat');
		slFileName = fullfile('/', 'scratch','bgcho','OpticalFlow','slMaps','slAzDependence.mat');
		load(slFileName,'azDeg','slHalf');
		azDeg = azDeg(:)';
		azDegLong = fliplr(180-[-fliplr(azDeg(2:end-1)) azDeg]);
		slHalfLong = fliplr([fliplr(slHalf(freqIndx,2:end-1)) slHalf(freqIndx,:)]);

		
		% bimapCount = 1;
		% for ii=1:length(fileAll)
		% 	if ~isempty(regexp(fileAll(ii).name,'bimap_user.txt'))			
		% 		bimapFile{bimapCount} = fileAll(ii).name;
		% 		YY = str2num(fileAll(ii).name(20:23));
		% 		JDN = str2num(fileAll(ii).name(26:28));
		% 		Greg = JDN2Greg(YY,JDN);
		% 		yy = Greg(1);
		% 		mm = Greg(2);
		% 		dd = Greg(3);
		% 		HH = str2num(fileAll(ii).name(30:31));
		% 		MM = str2num(fileAll(ii).name(32:33));
		% 		SS = str2num(fileAll(ii).name(34:35));
		% 		tStampBimap(bimapCount) = datenum([yy mm dd HH MM SS]);
		% 		bimapCount = bimapCount+1;
		% 	end
		% end
		% validMask = ismember(tStampBimap,tStamp);
		% bimapFile = bimapFile(validMask);
		% tStamp = tStampBimap;
		numPing = length(tStamp);

		for ii=1:numPing
			% fid = fopen(fullfile(splDir,bimapFile{ii}));
			fid = fopen(bimapFile{ii});
			while ~feof(fid)
				entry = fgetl(fid);			
				if ~isempty(regexp(entry,'BIMAP Source UTM X'))
					delim = regexp(entry,'=');
					srcUTMPos(ii,1) = str2num(entry(delim+1:end-1));
				elseif ~isempty(regexp(entry,'BIMAP Source UTM Y'))
					delim = regexp(entry,'=');
					srcUTMPos(ii,2) = str2num(entry(delim+1:end-1));
				elseif ~isempty(regexp(entry,'BIMAP Source Latitude \(Y\)'))
					delim = regexp(entry,'=');
					srcDegPos(ii,1) = str2num(entry(delim+1:end-1));
				elseif ~isempty(regexp(entry,'BIMAP Source Longitude \(X\)'))
					delim = regexp(entry,'=');
					srcDegPos(ii,2) = str2num(entry(delim+1:end-1));
				elseif ~isempty(regexp(entry,'BIMAP Array Heading \(True North\)'))				
					delim = regexp(entry,'=');
					arrayHeading(ii) = str2num(entry(delim+1:end-1));
				end
			end
			fclose(fid);
		end

		arrayHeading = 90-arrayHeading;
		ahMask = arrayHeading<=-180;
		arrayHeading(ahMask) = arrayHeading(ahMask)+360;

		
		slDiffLong = slHalfLong-max(slHalfLong);

		h = figure; 
		for ii=1:numPing
			[X,Y] = meshgrid(xgrid-srcUTMPos(ii,1),ygrid-srcUTMPos(ii,2));
			% TH = 90-180/pi*atan2(Y-srcUTMPos(ii,2), X-srcUTMPos(ii,1));
			TH = 180/pi*atan2(Y,X);
			THFromAH = TH-arrayHeading(ii);
			THFromAHVec = THFromAH(:);
			thMask = THFromAHVec<0;
			THFromAHVec(thMask) = 360+THFromAHVec(thMask);
			slDiffVec = interp1(azDegLong,slDiffLong,THFromAHVec,'linear','extrap');
			slDiffMap(:,:,ii) = reshape(slDiffVec,size(X));

			figure(h); clf(h)
			imagesc(squeeze(slDiffMap(:,:,ii)))
			axis ij equal tight	
			colorbar
		end

		close(h)


	end

	function tStamp = getTStampFromBimap(bimapFileName)

		% bimapCount = 1;


		for ii=1:length(bimapFileName)
			[fPath,fNameTmp,fExt] = fileparts(bimapFileName{ii});
			fName = [fNameTmp fExt];
			% if ~isempty(regexp(fileAll(ii).name,'bimap_user.txt'))			
			% bimapFile{bimapCount} = fileAll(ii).name;
			YY = str2num(fName(20:23));
			JDN = str2num(fName(26:28));
			Greg = JDN2Greg(YY,JDN);
			yy = Greg(1);
			mm = Greg(2);
			dd = Greg(3);
			HH = str2num(fName(30:31));
			MM = str2num(fName(32:33));
			SS = str2num(fName(34:35));
			tStamp(ii) = datenum([yy mm dd HH MM SS]);
			% bimapCount = bimapCount+1;
			% end
		end
	end



	function TS = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer, zNb, freqHz)
		
		name = fishStruct.species;
    dFish = fishStruct.fleshDensity;  %kg/m^3
    visFish = fishStruct.fleshViscosity; % Pa*sec
    lf = fishStruct.forkLength; % cm
    lt = fishStruct.totalLengthOffset + fishStruct.totalLengthExpFactor*lf; % cm
    b = (fishStruct.fracMajor*lt)/2;  % cm
    weight = fishStruct.wCoeff*lf.^(fishStruct.wPower); % g

    numLength = length(lf);
    V0Tmp = 0.05*(zNb/10+1)*weight(:)'; % 1*lf, cm^3

    numDepth = 20;
    zOc(1,:) = linspace(z1-H1/2, z1+H1/2, numDepth);
    zOc(2,:) = linspace(z2-H2/2, z2+H2/2, numDepth);

    V0 = repmat(reshape(V0Tmp,[1,numLength,1]), [numDepth,1,2]);	% numDepth*numLength*2, cm^3

    zOcMat(:,:,1) = repmat(zOc(1,:)', [1, numLength]);	% numDepth*numLength*2, m  
    zOcMat(:,:,2) = repmat(zOc(2,:)', [1, numLength]);
    V = V0./(1+zOcMat./10);
    r = (3/(4*pi)*V).^(1/3)/100;  % numDepth*numLength*2, m  

    bMat = repmat(reshape(b(:), [1, numLength, 1]), [numDepth, 1, 2]);
    a = sqrt(3/(4*pi)*1./(1+zOcMat./10).*(V0./bMat)); % numDepth*numLength*2, cm

    ecc = a./bMat;  % numDepth*numLength*2

    zeta = sqrt(2)*(1-ecc.^2).^(1/4)./ecc.^(1/3)./sqrt(log((1+sqrt(1-ecc.^2))./(1-sqrt(1-ecc.^2))));  % numDepth*numLength*2

    pAmb = 10^5;  % Ambient pressure: 10^5 Pa
    ratioHeat = 1.4;  % Specific heat ratio    

    f0 = 1/(2*pi)*zeta./r.*sqrt(3*ratioHeat*pAmb*(1+zOcMat./10)/dFish); % numDepth*numLength*2, Hz

    numFreq = length(freqHz);
    rMat = repmat(r, [1,1,1,numFreq]);	% numDepth*numLength*2*numFreq, m
    freqMat = repmat(reshape(freqHz(:),[1,1,1,numFreq]), [numDepth,numLength,2,1]);	% numDepth*numLength*2*numFreq, Hz
    f0Mat = repmat(reshape(f0,[numDepth,numLength,2,1]),[1,1,1,numFreq]);	% numDepth*numLength*2*numFreq, Hz
    c = 1500; % sound speed in water, m/s
    kMat = 2*pi*freqMat./c;	% numDepth*numLength*2*numFreq, 1/m
    
    hInv = (2*pi*rMat.*freqMat.^2)./(c*f0Mat) + visFish./(pi*dFish*rMat.^2.*f0Mat);
    fleshDampFactor = 2*visFish./(dFish*rMat.^2*2*pi.*freqMat);    
    % radDampFactor = totDampFactor - fleshDampFactor;   
    radDampFactor = kMat.*rMat;
    % totDampFactor = hInv.*f0Mat./freqMat;
    totDampFactor = radDampFactor + fleshDampFactor;

    TS.sigmaTot = (4*pi*rMat.^2)./(((f0Mat./freqMat).^2-1).^2+((f0Mat./freqMat).*hInv).^2);
    TS.sigmaBs = TS.sigmaTot./(4*pi);
    TS.f0 = f0;

    TS.sFnRad = ((f0Mat./freqMat).^2-1).*kMat.*rMat./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2) ...
              + 1i*radDampFactor.*(kMat.*rMat)./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    % scs.sFnFlesh = 1i*fleshDampFactor.*(kMat.*rMat).^2./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    TS.sFnFlesh = 1i*fleshDampFactor.*(kMat.*rMat)./(((f0Mat./freqMat).^2-1).^2 + totDampFactor.^2);
    TS.sFn = TS.sFnRad + TS.sFnFlesh;

    pdfMatFl = repmat(reshape(fishStruct.flpdf, [1,numLength,1,1]),[numDepth,1,2,numFreq]);	% numDepth*numLength*2*numFreq
    % pdfMatZ = repmat(reshape(zpdf(:,2), [numDepth,1,1,1]),[1,numNBDepth,numLength,numFreq]);
    dlf = mean(abs(diff(lf)));
    
    TS.sigmaBSMeanOverLf = squeeze(nansum(TS.sigmaBs.*pdfMatFl,2)*dlf);	% numDepth*2*numFreq

    TS.sigmaBSMeanOverLfDepthEachLayer = squeeze(nanmean(TS.sigmaBSMeanOverLf,1));	% 2*numFreq
    TS.sigmaBSMeanOverLfDepth = ratioLayer*TS.sigmaBSMeanOverLfDepthEachLayer(1,:) + (1-ratioLayer)*TS.sigmaBSMeanOverLfDepthEachLayer(2,:);
    TS.meanTS = 10*log10(TS.sigmaBSMeanOverLfDepth);

    TS.sFnMeanOverLf = squeeze(nansum(TS.sFn.*pdfMatFl,2)*dlf);	% 
    TS.sFnMeanOverLfDepthEachLayer = squeeze(nanmean(TS.sFnMeanOverLf,1));
    TS.meanSFn = ratioLayer*TS.sFnMeanOverLfDepthEachLayer(1,:) + (1-ratioLayer)*TS.sFnMeanOverLfDepthEachLayer(2,:);
	end

	function tsJacobian = calJacobianTS(fishStruct, z1, H1, z2, H2, ratioLayer, zNb, freqHz)
		% tsJacobian: numFreq*6
		% c = 1500;
		% k = 2*pi*freqHz./c;
		% nV = 10.^(nAdB./10)./waterDepth;
		dz = 1;
		dh = 1;
		dratio = 0.05;
		dznb = 1;

		TS2 = calTSFast(fishStruct, z1+dz, H1, z2, H2, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1-dz, H1, z2, H2, ratioLayer, zNb, freqHz);
		tsJacobian(:,1) = (TS2.meanTS - TS1.meanTS)./(2*dz);
		% alphaJacobian(:,1) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dz);

		TS2 = calTSFast(fishStruct, z1, H1+dh, z2, H2, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1-dh, z2, H2, ratioLayer, zNb, freqHz);
		tsJacobian(:,2) = (TS2.meanTS - TS1.meanTS)./(2*dh);
		% alphaJacobian(:,2) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dh);

		TS2 = calTSFast(fishStruct, z1, H1, z2+dz, H2, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2-dz, H2, ratioLayer, zNb, freqHz);
		tsJacobian(:,3) = (TS2.meanTS - TS1.meanTS)./(2*dz);
		% alphaJacobian(:,3) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dz);

		TS2 = calTSFast(fishStruct, z1, H1, z2, H2+dh, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2, H2-dh, ratioLayer, zNb, freqHz);
		tsJacobian(:,4) = (TS2.meanTS - TS1.meanTS)./(2*dh);
		% alphaJacobian(:,4) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dh);

		TS2 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer+dratio, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer-dratio, zNb, freqHz);
		tsJacobian(:,5) = (TS2.meanTS - TS1.meanTS)./(2*dratio);
		% alphaJacobian(:,5) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dratio);

		TS2 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer, zNb+dznb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer, zNb-dznb, freqHz);
		tsJacobian(:,6) = (TS2.meanTS - TS1.meanTS)./(2*dznb);
		% alphaJacobian(:,6) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dznb);
	end

	function [tsJacobian, alphaJacobian] = calJacobianTS2(fishStruct, z1, H1, z2, H2, ratioLayer, zNb, nAdB, waterDepth, freqHz)
		% tsJacobian: numFreq*6
		c = 1500;
		k = 2*pi*freqHz./c;
		nV = 10.^(nAdB./10)./waterDepth;
		dz = 1;
		dh = 1;
		dratio = 0.05;
		dznb = 1;

		TS2 = calTSFast(fishStruct, z1+dz, H1, z2, H2, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1-dz, H1, z2, H2, ratioLayer, zNb, freqHz);
		tsJacobian(:,1) = (TS2.meanTS - TS1.meanTS)./(2*dz);
		alphaJacobian(:,1) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dz);

		TS2 = calTSFast(fishStruct, z1, H1+dh, z2, H2, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1-dh, z2, H2, ratioLayer, zNb, freqHz);
		tsJacobian(:,2) = (TS2.meanTS - TS1.meanTS)./(2*dh);
		alphaJacobian(:,2) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dh);

		TS2 = calTSFast(fishStruct, z1, H1, z2+dz, H2, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2-dz, H2, ratioLayer, zNb, freqHz);
		tsJacobian(:,3) = (TS2.meanTS - TS1.meanTS)./(2*dz);
		alphaJacobian(:,3) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dz);

		TS2 = calTSFast(fishStruct, z1, H1, z2, H2+dh, ratioLayer, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2, H2-dh, ratioLayer, zNb, freqHz);
		tsJacobian(:,4) = (TS2.meanTS - TS1.meanTS)./(2*dh);
		alphaJacobian(:,4) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dh);

		TS2 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer+dratio, zNb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer-dratio, zNb, freqHz);
		tsJacobian(:,5) = (TS2.meanTS - TS1.meanTS)./(2*dratio);
		alphaJacobian(:,5) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dratio);

		TS2 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer, zNb+dznb, freqHz);
		TS1 = calTSFast(fishStruct, z1, H1, z2, H2, ratioLayer, zNb-dznb, freqHz);
		tsJacobian(:,6) = (TS2.meanTS - TS1.meanTS)./(2*dznb);
		alphaJacobian(:,6) = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS2.meanSFn - TS1.meanSFn)./(2*dznb);
	end

	function [stateFinal, stateCovFinal, correctionFinal] = EKF(fishStruct, stateInit, stateCovInit, stateCov, measurement, measurementCov, correctionInit, dr, waterDepth, freqHz)
		% state: [z1, H1, z2, H1, ratioLayer, zNb, nAdB]', 7*1 column vector
		% stateCovaraince: 7*1 matrix
		% measurement: numFreq*numPts matrix
		% correction: numFreq*1 row vector
		% freqHz = [850 955 1125 1335 1465 1600];
		stateInter = stateInit;
		stateCovInter = stateCovInit + stateCov;

		numPts = size(measurement,2);
		numFreq = length(freqHz);
		numStateMem = length(stateInit);
		measurementCorrected = measurement + repmat(correctionInit, [1,numPts]);

		
		
		TS = calTSFast(fishStruct, stateInit(1), stateInit(2), stateInit(3), stateInit(4), stateInit(5), stateInit(6), freqHz);
		
		residual = measurementCorrected(:) - repmat(TS.meanTS', [numPts,1]) - stateInit(7);
		% residual = residual(:);
		TSJac = calJacobianTS(fishStruct, stateInit(1), stateInit(2), stateInit(3), stateInit(4), stateInit(5), stateInit(6), freqHz);
		H = [TSJac ones(numFreq,1)];
		% H = [TSJac 10*log10(exp(1))*ones(numFreq,1)];
		HMat = repmat(H, [numPts, 1]);
		
		resCov = HMat*stateCovInter*HMat' + measurementCov;
		
		KalmanGain = stateCovInter*HMat'*inv(resCov);

		% keyboard

		% Update
		stateFinal = stateInit + KalmanGain*residual;
		I = eye(numStateMem, numStateMem);
		stateCovFinal = (I - KalmanGain*HMat)*stateCovInter;

		TSFinal = calTSFast(fishStruct, stateFinal(1), stateFinal(2), stateFinal(3), stateFinal(4), stateFinal(5), stateFinal(6), freqHz);
		c = 1500;
		k = 2*pi*freqHz./c;
		scsTot = 4*pi./(k.^2).*imag(TSFinal.meanSFn);
		nV = 10.^(stateFinal(7)./10)/waterDepth;
		% nV = stateFinal(7)/waterDepth;
		
		correctionFinal = correctionInit + 10*log10(exp(1))*(nV*dr*scsTot');
		
	end

	function [stateFinal, stateCovFinal] = EKF2(fishStruct, stateInit, stateCovInit, stateCov, measurement, measurementCov, dr, waterDepth, freqHz)
		% state: [z1, H1, z2, H1, ratioLayer, zNb, nAdB, DL_850 DL_955 DL_1125 DL_1335 DL_1465 DL_1600]', 13*1 column vector
		% stateCovaraince: 13*13 matrix
		% measurement: numFreq*numPts matrix
		% correction: numFreq*1 row vector
		% freqHz = [850 955 1125 1335 1465 1600];

		c = 1500;
		k = 2*pi*freqHz./c;
		nV = 10.^(stateInit(7)./10)/waterDepth;

		numPts = size(measurement,2);
		numFreq = length(freqHz);
		numStateMem = length(stateInit);
		% measurementCorrected = measurement;

		
		
		TS = calTSFast(fishStruct, stateInit(1), stateInit(2), stateInit(3), stateInit(4), stateInit(5), stateInit(6), freqHz);
		attnFactor = 10*log10(exp(1))*nV*4*pi./k.^2.*imag(TS.meanSFn);
		dAttndnAdB = 10*log10(exp(1))*4*pi./k.^2.*imag(TS.meanSFn)*log(10)/(10*waterDepth)*10^(stateInit(7)/10);

		
		
		% residual = residual(:);
		[TSJac, attnJac] = calJacobianTS2(fishStruct, stateInit(1), stateInit(2), stateInit(3), stateInit(4), stateInit(5), stateInit(6), stateInit(7), waterDepth, freqHz);
		FBlock1 = diag(ones(1,7),0);
		FBlock2 = zeros(7,6);
		FBlock3 = [attnJac dAttndnAdB(:)]*dr;
		FBlock4 = diag(ones(1,6),0);
		
		FMat = [FBlock1 FBlock2 ; FBlock3 FBlock4];

		H = [TSJac ones(numFreq,1) -diag(ones(1,6),0)];
		% H = [TSJac 10*log10(exp(1))*ones(numFreq,1)];
		HMat = repmat(H, [numPts, 1]);

		% Predict
		stateInter = stateInit + [zeros(7,1) ; attnFactor(:)*dr];
		stateCovInter = FMat*stateCovInit*FMat' + stateCov;

		% Update		
		residual = measurement(:) - (repmat(TS.meanTS', [numPts,1]) + stateInit(7) - repmat(stateInit(8:end), [numPts,1]) - repmat(attnFactor', [numPts,1])*dr);
		resCov = HMat*stateCovInter*HMat' + measurementCov;
		KalmanGain = stateCovInter*HMat'*inv(resCov);

		stateFinal = stateInit + KalmanGain*residual;
		I = eye(numStateMem, numStateMem);
		stateCovFinal = (I - KalmanGain*HMat)*stateCovInter;

		% TSFinal = calTSFast(fishStruct, stateFinal(1), stateFinal(2), stateFinal(3), stateFinal(4), stateFinal(5), stateFinal(6), freqHz);
		
		% scsTot = 4*pi./(k.^2).*imag(TSFinal.meanSFn);
		% nV = 10.^(stateFinal(7)./10)/waterDepth;
		% nV = stateFinal(7)/waterDepth;
		
		% correctionFinal = correctionInit + 10*log10(exp(1))*(nV*dr*scsTot');
		
	end

	function ssMapAvg = avgSS(splName, slMap, tlName, areaName, xgrid, ygrid, depthOption)
		for oo=1:length(splName)			
			splMap = read_arr(splName{oo});
			tlStruct = matfile(tlName{oo});
			switch lower(depthOption)
			case 'bottom'
				tlMap = tlStruct.TL_twoway_bott_dB_final;
			case 'fish'
				tlMap = tlStruct.TL_twoway_fish_dB_final;
			end
			% tlMap = tlStruct.TL_twoway_fish_dB_final;
			% areaMap = tlStruct.area_term_dB;		
			areaMap = 10*log10(read_arr(areaName{oo}));					
			ssMapBig = splMap - tlMap - areaMap;
			ssMap(:,:,oo) = 10*log10(interp2(tlStruct.datax, (tlStruct.datay)', 10.^(ssMapBig./10), xgrid, ygrid(:)))-slMap(:,:,oo);
		end
		ssMapAvg = 10*log10(squeeze(nanmean(10.^(ssMap./10),3)));		
    end
    
    function areaMapAvg = avgAreaMap(areaName, mName, xgrid, ygrid)
        for oo=1:length(areaName)			
            [xgridBig, ygridBig] = readMFile(mName{oo});                        
			areaMapBig = read_arr(areaName{oo});					
			areaMap(:,:,oo) = interp2(xgridBig, ygridBig', areaMapBig, xgrid, ygrid(:));
        end		
		areaMapAvg = squeeze(nanmean(areaMap,3));
    end

    function [xgrid, ygrid] = readMFile(mFileName)
        fid = fopen(mFileName);
        while(~feof(fid))
            tline = fgetl(fid);
            if ~isempty(regexp(tline, 'grid_inc =', 'once'))
                delim = regexp(tline, '=');
                grid_inc = str2double(tline(delim+2:end-1));
            elseif ~isempty(regexp(tline, 'grid_xmin = ', 'once'))
                delim = regexp(tline, '=');
                grid_xmin = str2double(tline(delim+2:end-1));
            elseif ~isempty(regexp(tline, 'grid_xmax = ', 'once'))
                delim = regexp(tline, '=');
                grid_xmax = str2double(tline(delim+2:end-1));
            elseif ~isempty(regexp(tline, 'grid_ymin = ', 'once'))
                delim = regexp(tline, '=');
                grid_ymin = str2double(tline(delim+2:end-1));
            elseif ~isempty(regexp(tline, 'grid_ymax = ', 'once'))
                delim = regexp(tline, '=');
                grid_ymax = str2double(tline(delim+2:end-1));
            end                    
        end
        xgrid = [grid_xmin:grid_inc:grid_xmax];
        ygrid = [grid_ymax:-grid_inc:grid_ymin];
    end
        
end


