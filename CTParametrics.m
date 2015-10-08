function CTParametrics
    clear all;
    close all;
    fclose all;
    warning off;

    
    
    %     CTInpfn = 'C:\Documents and Settings\maqsood.yaqub\My Documents\VUmc\code\Matlab\progs\CTPerfusieprogramma\Data\Parametrics\PeerAP.txt';
    %     ROIdata = load(fname, '-ascii');
    [fnameb, pathb] = uigetfile('*.txt',['Select Input TAC file']);
    cd(pathb);
    fnTACInput = [pathb fnameb];
    dataInput=importdata(fnTACInput,' ',3);
    yInput=dataInput.data;

    [fnameb, pathb] = uigetfile('*.txt',['Select Tumor TAC file']);
    fnTACTumour = [pathb fnameb];
    idname = fnameb;
    
    fnCSV = [fnTACTumour(1:end-4) '_resultsNew.csv'];
    fnPS = [fnTACTumour(1:end-4) '_resultsNew'];
    fnParametric = [fnTACTumour(1:end-4) '_Perfusion'];
    
    dataTumour=importdata(fnTACTumour,' ',3);
    yTumour=dataTumour.data;
    
    ROIdata = [yInput(:,1) yInput(:,4)-yInput(1,4) yTumour(:,4)-yTumour(1,4)];
    
    choiceParM = 0;
%    choiceParM = menu('Include parametric analysis ?', 'Yes', 'No');
    
    
    if choiceParM == 1
        [fnameb, pathb] = uigetfile('*.v',['Select ECAT file with dynamic CT data']);
        CTDynfn = [pathb fnameb];
        
        midposx = input('Enter position of Tumor (x in vinci)');
        midposy = input('Enter position of Tumor (y in vinci)');
        [mhI,shI,data] = readECAT7(CTDynfn);
        nrframesfixed = mhI.num_frames;
        Iwidth = 100;
        vol = zeros(Iwidth,Iwidth,shI{1}.z_dimension,nrframesfixed);
%    midposx = 300;
%    midposy = 525;
        posxl = midposx - Iwidth/2;
        posyl = midposy - Iwidth/2;
        for i=1:nrframesfixed
            tmpImg=single(data{i})*mhI.ecat_calibration_factor*shI{i}.scale_factor;
            vol(:,:,:,i) = tmpImg(posxl:posxl+Iwidth-1,posyl:posyl+Iwidth-1,:);
        end;
        clear tmpImg;
        clear data;
    else
        mhI = [];
        shI = [];
        vol = 0;
        CTDynfn ='';
        Iwidth = 0;
        posxl = 0;
        posyl = 0;
    end

    [spl] = initPSfile(0, fnPS, ['CTP analysis']);
    dataR = [];

    [headernm, dataR, spl] = analysisCTTACs(ROIdata, fnPS, fnCSV, spl, idname, dataR,choiceParM, mhI,shI,vol, CTDynfn, Iwidth, posxl, posyl,fnParametric);
    writecsv2(0,fnCSV,headernm,dataR, idname,',');
%    close all
    close 1 2 5
end

% ---------------------------------------------
function [hdr, dataR, spl] = analysisCTTACs(data, fnPS, fnCSV, spl, idname, dataR,choiceParM,mhI,shI,dataI, CTDynfn, Iwidth, posxl, posyl, fnParametric)
    [spl] = chapterPSfile(0, fnPS, idname);
    %data = load(fname, '-ascii');
    
    print('-dpsc', '-append', [fnPS '.ps']);
    spl=1; fig=fig_settings(1);


    % [time arterial tumor], alle tissue enhanced curves, i.e. non-contrast
    % image substracted

    % define t, arterial and tissue curve
    stp = data(2,1) - data(1,1);
    t_old = data(:,1);
    a_old = data(:,2);
    c_old = data(:,3);
    spl = addsubplotTAC(t_old, a_old, spl, 'Original arterial TAC', 'r.', fnPS);
    spl = addsubplotTAC(t_old, c_old, spl, 'Original Tumor TAC', 'b.', fnPS);

    % Correct fit and correct for delay
    t = data(:,1);
    a = data(:,2);
    c = data(:,3);
    percUp = 0.15;
    percDown = 0.15;
    [down50, up20] = excluderecirculation(a, percUp, percDown);
    [anew, A_art, B_art, D_art, toArt] = fitGamma(a(up20:down50), t(up20:down50));
    C_art = 1 / ((B_art^(A_art+1)) * (gamma(A_art+1)));
    afit = D_art*C_art*exp(-(t-toArt)/(B_art)).*(t-toArt).^A_art;
    afit(t<=toArt) = 0;
    spl = addsubplotTACs(t, a, afit, spl, 'Fitted arterial TAC', 'r.', fnPS, 'b--');

    choiceDebug=menu('TAC needs better fit?','Yeah','Noo');
    while choiceDebug==1
        percUp = input('Enter percUp for arterial TAC: ');
        percDown = input('Enter percDown: ');
        [down50, up20] = excluderecirculation(a, percUp, percDown);
        [anew, A_art, B_art, D_art, toArt] = fitGamma(a(up20:down50), t(up20:down50));
        C_art = 1 / ((B_art^(A_art+1)) * (gamma(A_art+1)));
        afit = D_art*C_art*exp(-(t-toArt)/(B_art)).*(t-toArt).^A_art;
        afit(t<=toArt) = 0;
        spl = addsubplotTACs(t, a, afit, spl, 'Fitted arterial TAC', 'r.', fnPS, 'b--');
        choiceDebug=menu('again?','Yes','Go away');
    end
    
    percUp = 0.15;
    percDown = 0.50;
    [down50, up20] = excluderecirculation(c, percUp, percDown);
    [anew, A_tum, B_tum, D_tum, toTum] = fitGamma(c(up20:down50), t(up20:down50));
    C_tum = 1 / ((B_tum^(A_tum+1)) * (gamma(A_tum+1)));
    cfit = D_tum*C_tum*exp(-(t-toTum)/(B_tum)).*(t-toTum).^A_tum;
    cfit(t<=toTum) = 0;
    spl = addsubplotTACs(t, c, cfit, spl, 'Fitted Tumor TAC', 'r.', fnPS, 'b--');

    choiceDebug=menu('TAC needs better fit?','Yeah','Noo');
    while choiceDebug==1
        percUp = input('Enter percUp for arterial TAC: ');
        percDown = input('Enter percDown: ');
        [down50, up20] = excluderecirculation(c, percUp, percDown);
        [anew, A_tum, B_tum, D_tum, toTum] = fitGamma(c(up20:down50), t(up20:down50));
        C_tum = 1 / ((B_tum^(A_tum+1)) * (gamma(A_tum+1)));
        cfit = D_tum*C_tum*exp(-(t-toTum)/(B_tum)).*(t-toTum).^A_tum;
        cfit(t<=toTum) = 0;
        spl = addsubplotTACs(t, c, cfit, spl, 'Fitted Tumor TAC', 'r.', fnPS, 'b--');
        choiceDebug=menu('again?','Yes','Go away');
    end

    if choiceParM == 1
        DoAndSaveParametric(afit, t, up20, down50, A_tum, B_tum, D_tum, toTum, dataR,mhI,shI,dataI, CTDynfn, Iwidth, posxl, posyl, A_art, B_art, D_art, toArt, fnParametric);
    end

    delayf = toTum - toArt;
    C_art = 1 / ((B_art^(A_art+1)) * (gamma(A_art+1)));
    afit = D_art*C_art*exp(-(t-toArt-delayf)/(B_art)).*(t-toArt-delayf).^A_art;
    afit(t<=toTum) = 0;
    spl = addsubplotTAC(t, afit, spl, 'Delay corrected and fitted arterial TAC', 'r.', fnPS);
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MOMENTS
    AUC_tissue = trapz(cfit)*stp;
    AUC_artery = trapz(afit)*stp;
    Hct_artery = 0.43;
    r = 0.7; % 0.8 for kids
    k = 0.05; % estimated, is the ratio of standard deviation of tissue transit times to the mean transit time 
    %cog_tissueA=integrate(t(length(t)),t,t.*cfit)/integrate(t(length(t)),t,cfit);
    cog_tissue = trapz(t.*cfit) / trapz(cfit);
    cog_artery = trapz(t.*afit) / trapz(afit);
    m = 0; % difference arrival time true artery vs measured artery
    n = (1-Hct_artery) * (1+k*k) / (2 * (1-r*Hct_artery)*(1-m));
    BF_moments = n * (AUC_tissue / AUC_artery) / (cog_tissue - cog_artery);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DECONV, Not working correct Yet
    
    noconf = 0;

    Amat = zeros(length(t),length(t));

    for ii=1:length(t)
        Amat(1:ii) = afit(ii:-1:1);
    end
        
    cfft = cfit;
    afft = afit;
    tfft = t;
%     if length(cfit) > 32
%         cfft = cfit(1:32);
%         afft = afit(1:32);
%         tfft = t(1:32);
%     else
%         if length(cfit) > 16
%             cfft = cfit(1:16);
%             afft = afit(1:16);
%             tfft = t(1:16);
%         else
%             if length(cfit) > 8
%                 cfft = cfit(1:8);
%                 afft = afit(1:8);
%                 tfft = t(1:8);
%             else
%                 if length(cfit) > 4
%                     cfft = cfit(1:4);
%                     afft = afit(1:4);
%                     tfft = t(1:4);
%                 end
%             end
%         end
%     end
%     
    RMultCBF = ifft(fft(cfft)./fft(afft));
    [val, imaxRMultCBF] = max(RMultCBF); 
    
    
    
    
    Xdat_CONV = RMultCBF(imaxRMultCBF:end);
    Ydat_CONV = tfft(imaxRMultCBF:end);
    
    
    [fitsCONV, MTT, CBF] = fitexp(Ydat_CONV, Xdat_CONV, tfft(imaxRMultCBF:end)); 
    
    
    %cn = [cfit' zeros(1,length(afit')-1)];
    %P_IRF = customDeconvolution(cfit', afit');
    spl = addsubplotTACs(tfft(imaxRMultCBF:end), RMultCBF(imaxRMultCBF:end), fitsCONV, spl, 'Deconvolution Rt*CBF fit', 'r.', fnPS, 'b--');
    P_V_deconv = CBF;
    
    
    %%[testc] = conv(P_IRF, afit);
    %%figure; plot([1:length(testc)], testc, 'r.', [1:length(cfit)], cfit, 'b.');
    %P_V_deconv = max(P_IRF);
    %MTT_deconv = trapz(P_IRF/P_V_deconv)*stp;
    %Vb_deconv =  trapz(P_IRF)*stp;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MULLANI-GOULD
    [val, itprime] = max(cfit);  
    P_V_mullanigould = cfit(itprime) / (trapz(afit(1:itprime))*stp);
    A_int_MG = integrate(t,t,afit)';
    spl = addsubplotTAC(t, A_int_MG, spl, 'Integral arterial', 'b.', fnPS);

    % Fick curve, blood flow, perfusion = slope
    spl = addsubplotFick(A_int_MG, cfit, spl, 'Fick, Tumor / Integral arterial', 'b.', fnPS);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SLOPE
    dc = diff(cfit)./diff(t);
    [val, imaxdc] = max(dc);
    [val, imaxa] = max(afit);
    P_V_slope = dc(imaxdc)/afit(imaxa);
    spl = addsubplotTAC(t(1:end-1), dc, spl, 'Differential Tumor', 'b.', fnPS);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use delay corrected and fitten curves for SKM
    % Patlack and Single Tissue Plasma Input Model don't use fitted curves
    % Single Tissue Plasma Input Model has its own delay correction
    %c_old = cfit;
    %a_old = afit;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % patlack
    Xdat_pat = integrate(t_old,t_old,a_old)'./(a_old + eps);
    Ydat_pat = c_old./(a_old + eps);
    lng = length(t_old);
    [fits, intercept_pat, slope_pat] = fitLin(Ydat_pat(round(lng*0.3):end), Xdat_pat(round(lng*0.3):end), Xdat_pat);
    spl = addsubplotPatlackFit(Xdat_pat, Ydat_pat, fits, spl, 'Patlack', 'r.', fnPS, 'b--');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single tissue model
    Xdat_stm = t_old;
    Ydat_stm = c_old;
    ca_stm = a_old;
    [fits, K1, k2, Vb, del] = fitSingleTissue(Ydat_stm, Xdat_stm, ca_stm);
    spl = addsubplotTACs(t_old, c_old, fits, spl, 'Fitted TAC Single Tissue Model', 'r.', fnPS, 'b--');


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results
    hdr = ['"Name","CT methods fit to_art (min)","CT methods fit to_Tum (min)","CT methods delay_corr (min)","Moments BF","Deconv PV","Mullani_Gould PV","Slope PV","Patlack K","Patlack rBV","STM K1","STM k2","STM Vb","STM Vd", "STM delay correction (min)",'];
    dataR = [dataR; toArt toTum toTum - toArt BF_moments P_V_deconv P_V_mullanigould  P_V_slope slope_pat intercept_pat K1 k2 Vb K1/k2 del];

    print('-dpsc', '-append', [fnPS '.ps']);
    spl=0; fig=fig_settings(1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DoAndSaveParametric(afit, t, up20, down50, Ain, Bin, Din, toin, dataR,mhI,shI,vol, CTDynfn, Iwidth, posxl, posyl, A_art, B_art, D_art, toArt, fnParametric)
    % initials
%    dimx = shI{1}.x_dimension;
%    dimy = shI{1}.y_dimension;
%    dimz = shI{1}.z_dimension;
    hdr = loadhdr('M:\DAC\SW\EXPLORATORY\MATLAB\CTP\tmp\tmp.hdr');

    dimx = Iwidth;
    dimy = Iwidth;
    dimz = shI{1}.z_dimension;

    lowerthrs = -300;
    upperthres = 500;

    lnT = length(vol(1,1,1,:));
    AVG = mean(vol, 4);
    Par = zeros(dimx, dimy, dimz);

    
    for i=1:lnT
        vol(:,:,:,i) = vol(:,:,:,i) - vol(:,:,:,1);
    end

    
    [val, imaxa] = max(afit);
                        
    h=waitbar(0.000001,'Reading dicom data');
    valC = 0;
    for ix=1:2:dimx
        for iy=1:2:dimy
            for iz=1:dimz
                valC = valC + 1;
                perc = valC*4/(double(dimx*dimy*dimz));
                str = ['Analyzing ' num2str(valC*4) ' of ' num2str(dimx*dimy*dimz) ''];
                waitbar(perc,h,str);
%                if AVG(ix,iy,iz) > lowerthrs & AVG(ix,iy,iz) < upperthres
                    tmp1 = vol(ix,iy,iz,:);
                    tmp2 = vol(ix+1,iy,iz,:);
                    tmp3 = vol(ix,iy+1,iz,:);
                    tmp4 = vol(ix+1,iy+1,iz,:);
                    c = (reshape(tmp1,lnT,1) + reshape(tmp2,lnT,1)+reshape(tmp3,lnT,1)+reshape(tmp4,lnT,1))*0.25;

%                    [down50, up20] = excluderecirculation(c, percUp, percDown);
                    [anew, A_tum, B_tum, D_tum] = fitGammaQ(c(up20:down50), t(up20:down50), Ain, Bin, Din, toin);
                    C_tum = 1 / ((B_tum^(A_tum+1)) * (gamma(A_tum+1)));
                    cfit = D_tum*C_tum*exp(-(t-toin)/(B_tum)).*(t-toin).^A_tum;
                    cfit(t<=toin) = 0;
                    c = cfit;

                    %delayf = toTum - toArt;
                    %C_art = 1 / ((B_art^(A_art+1)) * (gamma(A_art+1)));
                    %afit = D_art*C_art*exp(-(t-toArt-delayf)/(B_art)).*(t-toArt-delayf).^A_art;
                    %afit(t<=toTum) = 0;
                    %afit = afit;
                    
                    
                    % slope
                    dc = diff(cfit)./diff(t);
                    [val, imaxdc] = max(dc);

                    perfusion = dc(imaxdc)/afit(imaxa);
                    %

                    Par(ix  , iy  , iz) = perfusion;
                    Par(ix+1, iy  , iz) = perfusion;
                    Par(ix  , iy+1, iz) = perfusion;
                    Par(ix+1, iy+1, iz) = perfusion;
 %              end
            end
        end
    end

    % loop over data
    %% get a tac
    %% fit
    %[result, A, B, D, to] = fitGammaQ(Ydat, Xdat, Ain, Bin, Din, toin);
    %% apply delay
    %% calc parameters
    % end loop

    
    
    dimx = shI{1}.x_dimension;
    dimy = shI{1}.y_dimension;
    dimz = shI{1}.z_dimension;
    spx = shI{1}.x_pixel_size*10;
    spy = shI{1}.y_pixel_size*10;
    spz = shI{1}.z_pixel_size*10;

    Par2 = zeros(dimx,dimy,dimz);
    Par2(posxl:posxl+Iwidth-1,posyl:posyl+Iwidth-1,:) = Par(:,:,:);
    Par2 = flipdim(Par2,3);
    Par2 = flipdim(Par2,2);
    
    basefn = [fnParametric '_slopeMethod_'];
    
    WriteNormalisedVolume(Par2, basefn, hdr, dimx, dimy, dimz, spx, spy, spz);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IRF = customDeconvolution(Output, Input)
    %[IRF] = deconvreg(Output,Input,[],[1e-9 100]); % need to look at
    FF_Output = fft(Output, 64);
    FF_Input = fft(Input, 64);
    FF_dev = FF_Output ./ FF_Input;
    IRF = ifft(FF_dev);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [down50, up20] = excluderecirculation(cfit, percUp, percDown)
    [maxV,maxI]=max(cfit);
    down50 = maxI;
    for i=maxI:length(cfit)
        if (cfit(i) < cfit(down50)) && (cfit(i)>maxV*percDown)
            down50 = i;
        end
        if cfit(i)<maxV*percDown
            break;
        end
    end

    up20 = 1;
    for i=1:maxI
        if (cfit(i)>maxV*percUp)
            up20 = i;
        end
        if cfit(i)>maxV*percUp
            break;
        end
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fits, K1, k2, Vb, del] = fitSingleTissue(Ydat, Xdat, ca) 
    % K1 = parn(1);
    % k2 = parn(2);
    % Vb = parn(3);

    stpar = [1  1 0.05 0];
    lB = [1e-50 1e-50 0 -0.5];
    uB = [1e50 1e50 1 0.5];
    dependent = Ydat;
    independent = Xdat;

    Data = [independent dependent];
    Input = [independent ca ca];
    weights = [];
    %weights = dependent.*dependent;
    %weights = weights./max(weights);

    [par,rssbest,sdpar]=modFminsIt(['norm(fit_mod_stn_del(parIn,P1,P2,P3,P4))'],stpar,[],[], lB, uB, 5, 1, Data, Input, weights,0);
    [f,fits] = fit_mod_stn_del(par,Data,Input,[],0);
    K1 = par(1);
    k2 = par(2);
    Vb = par(3);
    del = par(4);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fits, MTT, CBF] = fitexp(Ydat, Xdat, fullX) 
    stpar = [1  1]; % cbf and MTT
    lB = [0.000001 0.000001];
    uB = [10000 10000];
    dependent = Ydat;
    independent = Xdat;

    Data = [independent dependent];
    Input = independent;
    weights = [];
    %weights = dependent.*dependent;
    %weights = weights./max(weights);

    [par,rssbest,sdpar]=modFminsIt(['norm(fit_exp(parIn,P1,P2,P3,P4))'],stpar,[],[], lB, uB, 1, 1, Data, Input, weights,0);
    MTT = par(2);
    CBF = par(1);
    fits = CBF*exp(-fullX/MTT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fits, intercept, slope] = fitLin(Ydat, Xdat, fullX) 
    stpar = [1  0];
    lB = [-10000 -10000];
    uB = [10000 10000];
    dependent = Ydat;
    independent = Xdat;

    Data = [independent dependent];
    Input = independent;
    weights = [];
    %weights = dependent.*dependent;
    %weights = weights./max(weights);

    [par,rssbest,sdpar]=modFminsIt(['norm(fit_lin(parIn,P1,P2,P3,P4))'],stpar,[],[], lB, uB, 1, 1, Data, Input, weights,0);
    intercept = par(2);
    slope = par(1);
    fits = slope*fullX + intercept;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, A, B, D, to, tob, toe] = fitGamma(Ydat, Xdat)

    tob = Xdat(1)-5/60;
    toe = Xdat(1) - eps;
    tos = (toe+tob)/2;

    %pars = [A B D t0]
    %stpar = [75 0.005 1 0];
    %lB = [0 0.001 0.000001 tob];
    %uB = [150 1 1000 toe];

    stpar = [2.7 0.08 100 0];
    lB = [0 0.001 0.000001 tob];
    uB = [35 2 100000 toe];



    dependent = Ydat;
    independent = Xdat;

    Data = [independent dependent];
    Input = independent;
    weights = [];
    %weights = dependent.*dependent;
    %weights = weights./max(weights);

    [par,rssbest,sdpar]=modFminsIt(['norm(fit_gammavariate(parIn,P1,P2,P3,P4))'],stpar,[],[], lB, uB, 0, 1, Data, Input, weights,0);

    t = Xdat;
    A = par(1);
    B = par(2);
    D =  par(3);
    to = par(4);
    C = 1 / ((B^(A+1)) * (gamma(A+1)));
    result = D*C*exp(-(t-to)/(B)).*(t-to).^A;

    %plot(Xdat,Ydat,'bo',Xdat,result,'r.-');
    %drawnow;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, A, B, D] = fitGammaQ(Ydat, Xdat, Ain, Bin, Din, toin)
%    tob = Xdat(1)-5/60;
%    toe = Xdat(1) - eps;
%    tos = (toe+tob)/2;

%    stpar = [Ain Bin Din toin];
%    lB = [0 0.001 0.000001 tob];
%    uB = [35 2 1000 toe];
    
    stpar = [Ain Bin Din];
    lB = [0 0.001 0.000001];
    uB = [35 2 1000];
    
    dependent = Ydat;
    independent = Xdat;

    Data = [independent dependent];
    Input = independent;
    weights = [];
    %weights = dependent.*dependent;
    %weights = weights./max(weights);

    [par,rssbest,sdpar]=modFminsItQ(['norm(fit_gammavariateQ(parIn,P1,P2,P3,P4))'],stpar,[],[], lB, uB, 0, 1, Data, Input, weights,0, toin);

    t = Xdat;
    A = par(1);
    B = par(2);
    D =  par(3);
    C = 1 / ((B^(A+1)) * (gamma(A+1)));
    result = D*C*exp(-(t-toin)/(B)).*(t-toin).^A;

    %plot(Xdat,Ydat,'bo',Xdat,result,'r.-');
    %drawnow;

end


function result = integrate(tm,tbl,cbl)
    %
    % INTEGRATE (AAL - 89/05/25) is used for integrating the curve (tbl,cbl) at
    % times tm:
    %
    p=[1 0];
    result=convexvarc(tm,p',tbl,cbl);
end

%----------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = loadhdr(fn); % from spm2
    fid = fopen(fn,'r','ieee-le');
    [F,count] = fread(fid, 1, 'long');
    fclose(fid);

    if F == 348
       fid = fopen(fn,'r','ieee-le');
    else
       fid = fopen(fn,'r','ieee-be');
       [F,count] = fread(fid, 1, 'long');
       fclose(fid);
       if F == 348
          %fprintf('\n reading ieee-be');
          fid = fopen(fn,'r','ieee-be');
       else
          fprintf('\nHDR-file is not Analyze format!!\n\n');
          fid = -1;      
       end
    end

    if fid~=-1
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       fseek(fid,0,'bof');
       hk.sizeof_hdr 		= fread(fid,1,'int32'); %4
       hk.data_type  		= mysetstr(fread(fid,10,'uchar'))'; %14
       hk.db_name    		= mysetstr(fread(fid,18,'uchar'))'; %32
       hk.extents    		= fread(fid,1,'int32'); %36
       hk.session_error	= fread(fid,1,'int16'); %38
       hk.regular			= mysetstr(fread(fid,1,'uchar'))'; %39
       hk.hkey_un0			= mysetstr(fread(fid,1,'uchar'))'; %40
       dime.dim		= fread(fid,8,'int16')'; %56
       dime.vox_units	= mysetstr(fread(fid,4,'uchar'))'; %60
       dime.cal_units	= mysetstr(fread(fid,8,'uchar'))'; %68
       dime.unused1	= fread(fid,1,'int16'); %70
       dime.datatype	= fread(fid,1,'int16'); %72
       dime.bitpix		= fread(fid,1,'int16'); %74
       dime.dim_un0	= fread(fid,1,'int16'); %76
       dime.pixdim		= fread(fid,8,'float')'; %108
       dime.vox_offset	= fread(fid,1,'float'); %112
       dime.funused1	= fread(fid,1,'float'); %116
       dime.funused2	= fread(fid,1,'float'); %120
       dime.funused3	= fread(fid,1,'float'); %124
       dime.cal_max	= fread(fid,1,'float'); %128
       dime.cal_min	= fread(fid,1,'float'); %132
       dime.compressed	= fread(fid,1,'int32'); %136
       dime.verified	= fread(fid,1,'int32'); %140
       dime.glmax		= fread(fid,1,'int32'); %144
       dime.glmin		= fread(fid,1,'int32'); %148
       hist.descrip	= mysetstr(fread(fid,80,'uchar'))'; %228
       hist.aux_file	= mysetstr(fread(fid,24,'uchar'))'; %252
       hist.orient		= fread(fid,1,'uchar'); %253
       hist.origin		= fread(fid,5,'int16')';%263
       hist.generated	= mysetstr(fread(fid,10,'uchar'))'; %273
       hist.scannum	= mysetstr(fread(fid,10,'uchar'))'; %283
       hist.patient_id	= mysetstr(fread(fid,10,'uchar'))'; %293
       hist.exp_date	= mysetstr(fread(fid,10,'uchar'))'; %303
       hist.exp_time	= mysetstr(fread(fid,10,'uchar'))'; %313
       hist.hist_un0	= mysetstr(fread(fid,3,'uchar'))'; %316
       hist.views		= fread(fid,1,'int32'); %320
       hist.vols_added	= fread(fid,1,'int32'); %324
       hist.start_field= fread(fid,1,'int32'); %328
       hist.field_skip	= fread(fid,1,'int32'); %332
       hist.omax		= fread(fid,1,'int32'); %336
       hist.omin		= fread(fid,1,'int32'); %340
       hist.smax		= fread(fid,1,'int32'); %344
       hist.smin		= fread(fid,1,'int32'); %348

       hdr.hk   = hk;
       hdr.dime = dime;
       hdr.hist = hist;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteNormalisedVolume(frameA, basefn, hdr, xdim, ydim, zdim, spx, spy, spz)
    IMGFileOut = [basefn '.img'];
    HDRFileOut = [basefn '.hdr'];

    
    if exist(IMGFileOut)==2
        delete IMGFileOut
    end
    if exist(HDRFileOut)==2
        delete HDRFileOut
    end


    ENDIAN = 'ieee-le'; % for files to be created
    Data_type = 'float';
    IOWriteAnalyzeImg(IMGFileOut,frameA, Data_type, ENDIAN, 1);

    minfr = 0;
    maxfr = 0;
    if max(frameA) > maxfr
        maxfr = max(frameA);
    end
    if min(frameA) < minfr
        minfr = min(frameA);
    end

    hdr.dime.dim(1) = 4;
    hdr.dime.dim(2) = xdim;
    hdr.dime.dim(3) = ydim;
    hdr.dime.dim(4) = zdim;
    hdr.dime.pixdim(2) = spx;
    hdr.dime.pixdim(3) = spy;
    hdr.dime.pixdim(4) = spz;
            
    hdr.dime.datatype = 16;
    hdr.dime.bitpix = 32;
    hdr.dime.dim(5) = 1;
    hdr.dime.glmax		= ceil(maxfr); %144
    hdr.dime.glmin		= floor(minfr); %148
    hdr.dime.cal_max	= maxfr; %128
    hdr.dime.cal_min	= minfr; %132
    hdr.dime.funused1 = 1;
    savehdr(HDRFileOut, hdr, ENDIAN); % saving default hdr from file 1, adjusted for nrframes
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IOWriteAnalyzeImg(IMGFile, dataArray, typeData, ENDIAN, index)
    fprintf('\n writing ieee-le');
    if index == 1
        FileID = fopen(IMGFile,'w',ENDIAN);
    else
        FileID = fopen(IMGFile,'a',ENDIAN);
    end

    if FileID ~= -1
        fwrite(FileID, dataArray, typeData);
    end
    fclose(FileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savehdr(fn, hdr, ENDIAN); % from spm2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(fn,'w',ENDIAN);

    %hk.sizeof_hdr 		= fread(fid,1,'int32'); %4
    cc = fwrite(fid,hdr.hk.sizeof_hdr,'int32');
    %hk.data_type  		= mysetstr(fread(fid,10,'uchar'))'; %14
    cc = fwrite(fid,hdr.hk.data_type,'uchar');
    %hk.db_name    		= mysetstr(fread(fid,18,'uchar'))'; %32
    cc = fwrite(fid,hdr.hk.db_name,'uchar');
    %hk.extents    		= fread(fid,1,'int32'); %36
    cc = fwrite(fid,hdr.hk.extents,'int32');
    %hk.session_error	= fread(fid,1,'int16'); %38
    cc = fwrite(fid,hdr.hk.session_error,'int16');
    %hk.regular			= mysetstr(fread(fid,1,'uchar'))'; %39
    cc = fwrite(fid,hdr.hk.regular,'uchar');
    %hk.hkey_un0			= mysetstr(fread(fid,1,'uchar'))'; %40
    cc = fwrite(fid,hdr.hk.hkey_un0,'uchar');
    %dime.dim		= fread(fid,8,'int16')'; %56
    cc = fwrite(fid,hdr.dime.dim,'int16');
    %dime.vox_units	= mysetstr(fread(fid,4,'uchar'))'; %60
    cc = fwrite(fid,hdr.dime.vox_units,'uchar');
    %dime.cal_units	= mysetstr(fread(fid,8,'uchar'))'; %68
    cc = fwrite(fid,hdr.dime.cal_units,'uchar');
    %dime.unused1	= fread(fid,1,'int16'); %70
    cc = fwrite(fid,hdr.dime.unused1,'int16');
    %dime.datatype	= fread(fid,1,'int16'); %72
    cc = fwrite(fid,hdr.dime.datatype,'int16');
    %dime.bitpix		= fread(fid,1,'int16'); %74
    cc = fwrite(fid,hdr.dime.bitpix,'int16');
    %dime.dim_un0	= fread(fid,1,'int16'); %76
    cc = fwrite(fid,hdr.dime.dim_un0,'int16');
    %dime.pixdim		= fread(fid,8,'float')'; %108
    cc = fwrite(fid,hdr.dime.pixdim,'float');
    %dime.vox_offset	= fread(fid,1,'float'); %112
    cc = fwrite(fid,hdr.dime.vox_offset,'float');
    %dime.funused1	= fread(fid,1,'float'); %116
    cc = fwrite(fid,hdr.dime.funused1,'float');
    %dime.funused2	= fread(fid,1,'float'); %120
    cc = fwrite(fid,hdr.dime.funused2,'float');
    %dime.funused3	= fread(fid,1,'float'); %124
    cc = fwrite(fid,hdr.dime.funused3,'float');
    %dime.cal_max	= fread(fid,1,'float'); %128
    cc = fwrite(fid,hdr.dime.cal_max,'float');
    %dime.cal_min	= fread(fid,1,'float'); %132
    cc = fwrite(fid,hdr.dime.cal_min,'float');
    %dime.compressed	= fread(fid,1,'int32'); %136
    cc = fwrite(fid,hdr.dime.compressed,'int32');
    %dime.verified	= fread(fid,1,'int32'); %140
    cc = fwrite(fid,hdr.dime.verified,'int32');
    %dime.glmax		= fread(fid,1,'int32'); %144
    cc = fwrite(fid,hdr.dime.glmax,'int32');
    %dime.glmin		= fread(fid,1,'int32'); %148
    cc = fwrite(fid,hdr.dime.glmin,'int32');
    %hist.descrip	= mysetstr(fread(fid,80,'uchar'))'; %228
    cc = fwrite(fid,hdr.hist.descrip,'uchar');
    %hist.aux_file	= mysetstr(fread(fid,24,'uchar'))'; %252
    cc = fwrite(fid,hdr.hist.aux_file,'uchar');
    %hist.orient		= fread(fid,1,'uchar'); %253
    cc = fwrite(fid,hdr.hist.orient,'uchar');
    %hist.origin		= fread(fid,5,'int16')';%263
    cc = fwrite(fid,hdr.hist.origin,'int16');
    %hist.generated	= mysetstr(fread(fid,10,'uchar'))'; %273
    cc = fwrite(fid,hdr.hist.generated,'uchar');
    %hist.scannum	= mysetstr(fread(fid,10,'uchar'))'; %283
    cc = fwrite(fid,hdr.hist.scannum,'uchar');
    %hist.patient_id	= mysetstr(fread(fid,10,'uchar'))'; %293
    cc = fwrite(fid,hdr.hist.patient_id,'uchar');
    %hist.exp_date	= mysetstr(fread(fid,10,'uchar'))'; %303
    cc = fwrite(fid,hdr.hist.exp_date,'uchar');
    %hist.exp_time	= mysetstr(fread(fid,10,'uchar'))'; %313
    cc = fwrite(fid,hdr.hist.exp_time,'uchar');
    %hist.hist_un0	= mysetstr(fread(fid,3,'uchar'))'; %316
    cc = fwrite(fid,hdr.hist.hist_un0,'uchar');
    %hist.views		= fread(fid,1,'int32'); %320
    cc = fwrite(fid,hdr.hist.views,'int32');
    %hist.vols_added	= fread(fid,1,'int32'); %324
    cc = fwrite(fid,hdr.hist.vols_added,'int32');
    %hist.start_field= fread(fid,1,'int32'); %328
    cc = fwrite(fid,hdr.hist.start_field,'int32');
    %hist.field_skip	= fread(fid,1,'int32'); %332
    cc = fwrite(fid,hdr.hist.field_skip,'int32');
    %hist.omax		= fread(fid,1,'int32'); %336
    cc = fwrite(fid,hdr.hist.omax,'int32');
    %hist.omin		= fread(fid,1,'int32'); %340
    cc = fwrite(fid,hdr.hist.omin,'int32');
    %hist.smax		= fread(fid,1,'int32'); %344
    cc = fwrite(fid,hdr.hist.smax,'int32');
    %hist.smin		= fread(fid,1,'int32'); %348
    cc = fwrite(fid,hdr.hist.smin,'int32');
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mysetstr(in)
    tmp = find(in == 0);
    tmp = min([min(tmp) length(in)]);
    out = setstr([in(1:tmp)' zeros(1,length(in)-(tmp))])';
end

