
function [l, ps, ps_err, pmax,ci_pmax, ek, ci_ek, alpha, ci_a, R2, mused, modelfit]= PvI_ys4(light, pvar, calc_no, mused)

% Analyze PvI curves/relationships. Apply one of 2 possible PvI models
% (See Yang et al., 2020). Keep best fitting model (lowest R2). 
%
%INPUT
% light = total DC light during light curves. 
% pvar = ETR, or carbon uptake rates. must be same length as light
% Set calc_no to 1 to calculate ETR using the kinetic algorithm (Gorbunov
% et al., 2020). calc_no 2 will calculate ETR using the amplitude based
% algorithm. 
% mused - set the PI model to be used. optional parameter, default is best fit. 1 = no photoinhibition model, 2 = photoinhibition 
%
%OUTPUT
% l = light level 
% ps = calculated etr
% ps_err = standard error of ps estimates
% pmax = maximum rate of ETR. ci_pmax = 95% confidence interval
% ek = light saturation level. w 95% confidence interval (ci_ek)
% alpha = initial light limited slope of etr with 95% confidence interval
% (ci_a)
% R2 = model fit R2. 
% mused = model number used (1-4 value)
% parameters = 2-3 parameters of model used 
% uncertainty = 95% confidence interval of each parameter
%
% Y.Sez Updated Oct 2022
% Esat and Emax more specified for each curve. 
% only 2 models instead of 4
% ETR calculated outside this fxn. 

%define a threshold CoV to discard data points
threshold= .1;

%choose how many points to drop from the initial light period at each light

skip_points= 2;

switch calc_no
    case 1
        colored = [0 0 0];
    case 2
        colored = [1 0 0];
end

%Light step ETR values

l = unique(light);
[ps, ps_err,ps_std,count] = deal(nan * ones(size(l)));
if max(l) >= 40
    for i = 1:numel(l)
        light_step = find(light == l(i));
        light_step = light_step(skip_points:end);
        out = isoutlier(pvar(light_step));
        count(i) = sum(~out);
        if count(i) < 3
            break
        end
        ps(i) = nanmedian(pvar(light_step(~out)));
        ps_err(i) = nanstd(pvar(light_step(~out)))./sqrt(sum(~out));
        ps_std(i) = mad(pvar(light_step(~out)));
    end
end

cov = ps_std./ps; % coefficient of variance = std normalized to mean
bad_data = any([isnan(ps), cov > threshold, isinf(ps), count < 3],2);
ps(bad_data) = nan; ps_err(bad_data) = nan;
l(bad_data) = nan;

% curve fit: 2 PI models with photoinhibition terms: use the best fit.
model = {'pmax .* ((alpha .* x) ./(pmax + alpha.* x))';...
    'ps .* ((1 - exp(-alpha .* x./ps)) .* exp(-beta .* x ./ ps))'};
coeff = [{'pmax', 'alpha',nan};{'ps','alpha','beta'}];

[pmax_guess, emaxind] = max(ps); alpha_guess = pmax_guess./l(emaxind); 
xdat= linspace(1, max(light) + 100, 1000);
initial_guess = [pmax_guess alpha_guess nan; 1.1*pmax_guess alpha_guess 0.5*alpha_guess]; % [pmax, alpha, beta]

if numel(light(~bad_data)) < 5 % min points for curve fit w 3 params
    [pmax,ci_pmax,ek, ci_ek, alpha, ci_a, R2,mused]= deal(nan);
    
    % Edit saved and returned ps data
    [ps, ps_err,l] = deal(nan(size(l)));
    modelfit = struct();

    return
else
    
    output = struct();
    [rmse, r2] = deal(ones(size(model)));
    for i = 1:numel(model)
        nonan = isnan(initial_guess(i,:)); %Deal with models that have different number of coefficients
        ft = fittype(model{i},'coeff', coeff(i,~nonan));
        [params,stats]= fit(l(~bad_data),ps(~bad_data),ft,'startpoint',initial_guess(i,~nonan));
        output.(['model_' num2str(i)]).('params') = params;
        output.(['model_' num2str(i)]).('stats') = stats;
        rmse(i) = stats.rmse;
        r2(i) = stats.rsquare;      
    end
    
    %Model 1 results:
    pmax = output.model_1.params.pmax;
    alpha = output.model_1.params.alpha;
            
    %Model 2 results:
    Ps = output.model_2.params.ps;
    alpha1 = output.model_2.params.alpha;
    beta1 = output.model_2.params.beta;
    
    
    pred = [pmax .* (alpha .* xdat./ (pmax + alpha.*xdat));...
        Ps .* ((1 - exp(-alpha1 .* xdat./Ps)) .* exp(-beta1 .* xdat ./ Ps))];
    
if ~exist('mused','var')
    [~, mused] = min(rmse); % else, choose the best model of the 3 options
end
% If there is no photoinhibition, don't use the photoinhibition model
if beta1 <= 0
    mused = 1;
end

    R2 = r2(mused);
    
    ci = confint(output.(['model_' num2str(mused)]).params,0.95);
    ci(isnan(ci)) = 0;
    
    modelfit = output.(['model_' num2str(mused)]).params;
    
    % Calculate Pmax and Ek + uncertainties only for the best model.
    
    switch mused
        case 1
            
            %Model 1 results:
            ek = pmax./alpha;
            ci_pmax = (ci(2,1) - ci(1,1))./2;
            ci_a = (ci(2,2) - ci(1,2))./2;
            ci_ek = sqrt((ci_pmax./pmax)^2 + (ci_a./alpha)^2);
            
        case 2
            
            ek = Ps ./ alpha1 .* log((beta1 + alpha1)./beta1);
            pmax = Ps .* (alpha1./(alpha1 + beta1)) .* (beta1./(alpha1 + beta1)).^(beta1./alpha1);
            ci_ps = (ci(2,1) - ci(1,1))./2;
            ci_alpha = (ci(2,2) - ci(1,2))./2;
            ci_beta = (ci(2,3) - ci(1,3))/2;
            alpha = alpha1; ci_a = ci_alpha;
            
            % Propagating error (Need to review this)
            
            syms P_S ALPH BET
            f = P_S./ALPH .* log((ALPH + BET)/BET);
            g = P_S .* (ALPH./(ALPH + BET)) .* (BET./(ALPH + BET))^(BET/ALPH);
            
            dek_dps = vpa(subs(diff(f,P_S),[P_S,ALPH,BET],[Ps,alpha1,beta1]));
            dek_dalpha = vpa(subs(diff(f,ALPH),[P_S,ALPH,BET],[Ps,alpha1,beta1]));
            dek_dbeta = vpa(subs(diff(f,BET),[P_S,ALPH,BET],[Ps,alpha1,beta1]));
            
            dpmax_dps = vpa(subs(diff(g,P_S),[P_S,ALPH,BET],[Ps,alpha1,beta1]));
            dpmax_dalpha = vpa(subs(diff(g,ALPH),[P_S,ALPH,BET],[Ps,alpha1,beta1]));
            dpmax_dbeta = vpa(subs(diff(g,BET),[P_S,ALPH,BET],[Ps,alpha1,beta1]));
            
            ci_ek = sqrt(double((ci_ps.* dek_dps)^2 + (ci_alpha.* dek_dalpha)^2 + (ci_beta.* dek_dbeta)^2));
            ci_pmax = sqrt(double(ci_ps^2 .* dpmax_dps^2 + ci_alpha^2 .* dpmax_dalpha^2 + ci_beta^2 .* dpmax_dbeta^2));

        
    end
    
    
    % ------------- Plotting-------------------------------------------
    errorbar(l(~bad_data),ps(~bad_data),ps_err(~bad_data),'o','MarkerFaceCol',colored,'linestyle','n'); hold on
    plot(l(bad_data),ps(bad_data),'x','markeredgecol',[0.5 0.5 0.5],'markersize',10, 'markerfacecol',[0.5 0.5 0.5],'linestyle','n');
   % plot(light,pvar,'.','markeredgecol',[.7 .7 .7]),  %plots all of the individual data points for each light step
    %c = {'r--', 'b--', 'c--','g--'};
    %for i = 1:numel(model)
        plot(xdat,pred(mused,:),'color',colored,'linewidth',0.5);    %fitted curves
    %end
    xlabel('PAR (\muE m^{-2} s^{-1})')
    ylabel('ETR (mol e mol RCII^{-1} s^{-1})')
%     xlim([0 max(l(~bad_data))])
%     ylim([0 max(ps(~bad_data) + ps_err(~bad_data))])

    set(gca,'xminorti','on','yminorti','on','xtickmode','auto','ytickmode','auto','box','on',...
        'xticklabelmode','auto','yticklabelmode','auto', 'FontSize', 20);
    
    %add results to plot window
%    text(.6,.15,['P_{max} =  ', num2str(round(pmax)), '\pm', num2str(round(ci_pmax))], 'unit','norm','fontsize',16);
%    text(.6,.24,strcat('e_k =  ','  ',num2str(round(ek)), '\pm', num2str(round(ci_ek))),'unit','norm','fontsize',16);

    text(.6, .33 + (calc_no - 1) *.1, ['model ', num2str(mused), ' used'],'unit','norm','fontsize',16,'color',colored);
    if r2(mused)<.9
        text(.6,.9,'Low fit quality','unit','norm')
    end
    if any(isinf(ps))
        text(.6, .8, 'ETRk = inf', 'unit','norm')
    end
    set(gca,'xminorti','on','yminorti','on','xtickmode','auto','ytickmode','auto','box','on');
end

% Edit saved and returned ps data
ps(bad_data) = nan;
ps_err(bad_data) = nan;
l(bad_data) = nan;

return