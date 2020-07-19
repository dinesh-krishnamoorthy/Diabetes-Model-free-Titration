clear
clc
par.Ts = 1/24; % 1 hour

% Virtual patient parameters
par.tau1 = 0.5;
par.p2 = 15.8;     % mean value from Kanderian et al.
par.pGEZI = 3.31;  % mean value from Kanderian et al.

[f,sys] = MVP_T2D(par);
rng(3)
nMC = 100;
CV =  abs(2 + 1.*randn(nMC,1));

rng(2)
for iMC = 1:nMC
    
    % Virtual patient initial state
    mc.SI(iMC) = 1.0 + (2.4-1.0)*rand;
    mc.EGP(iMC) = 350 + (380-350)*rand;
    mc.B(iMC) = 1.1 + (1.5-1.1)*rand;
    d_in = [mc.SI(iMC);mc.EGP(iMC);mc.B(iMC)];
    
    u0  = 5*24; % Initial insulin
    xfs = [0;0;15.2;12];
    
    % setup
    dT = 6/(60*24); % 6 min
    Ns = 60/6; % 1hour/6min = 10 samples
    k = 1;
    
    % RLS initialization
    lambda = 0.01;
    P = eye(2);
    theta  = [0;0];
    u(1) = u0;
    
    method = 1;
    
    adhere_1 = 1;
    sim.u(k) = u0;
    u_last = u0;
    for sim_k = 1:60*24
        
        if  rem(sim_k,24)==0
            k = k+1;
            adhere =1; % round(0.35 + 0.65*rand);
            
            if adhere
                % Read SMBG data
                SMBG(k) = xfs(4);
                % Update Dose
                switch method
                    case 1
                        [u_in,theta,P,g(k)] = dose_guidance(u(k-1),SMBG(k),k,theta,P,lambda,adhere_1);
                    case 2
                        u_in = dose_202(u(k-1),SMBG(k),k);
                    case 3
                        u_in = dose_stepwise(u(k-1),SMBG(max(1,k-2:k)),k);
                end
                u(k) = u_in;
                sim.u(k) = u_in;
                adhere_1 = adhere;
                
            else
                SMBG(k) = NaN;  %fsolve(@(x) cost(x,theta'*[u_last;1]),7);
                u(k) = u(k-1);
                u_in = 0;
                adhere_1 = adhere;
                sim.u(k) = u_in;
                
            end
        else
            u_in=0;
        end
        
        % -----------  Simulator -----------
        
        for i = 1:Ns
            % Stocahstic simulation - Euler Maruyama
            xfs = xfs + full(f(xfs,u_in,d_in))*dT + CV(iMC)*sqrt(dT)*randn(4,1);
        end
    end
    
    mc.SMBG(:,iMC) = SMBG';
    mc.Ib(:,iMC) = sim.u./24;
end

save('mc','mc')

%%
x1 = mc.SMBG(1:end,:);
x2 = mc.SMBG([10:end],:);

days = (1:k-1);
figure(22)
clf
subplot(2,2,[1:2])
hold all
fill([0 k-1 k-1 0],[4 4 6 6],[0.8 0.8 0.8],'facealpha',.35)
plot(days,mc.SMBG(2:end,:),'r','linewidth',0.5)
ylim([3,17])
title ('SMBG data for 50 virtual patients','Interpreter','Latex')
ylabel('SMBG [mmol/L]','Interpreter','Latex')
xlabel('days','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 12;
box on
grid on

glyc_var = sum((SMBG-5).^2 + 8.*(min(0,SMBG-5)).^2);
avg_gl = mean(SMBG);
avg_30 = [mean(SMBG(1:30));mean(SMBG(31:end))];

subplot(223)
hold all
fill([4 6 6 4],[0 0 2e4 2e4],[0.8 0.8 0.8],'facealpha',.35)
histogram(x2,'Normalization','pdf','FaceColor','b')
title ('Historgram after 10 days','Interpreter','Latex')
ylim([0,1.2])
xlim([3,8])
box on
xlabel('SMBG [mmol/L]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 12;
box on
grid on

subplot(224)
hold all
fill([4 6 6 4],[0 0 2e4 2e4],[0.8 0.8 0.8],'facealpha',.35)

for i =1:nMC
    [N,edges] = histcounts(mc.SMBG([14:end],i),10, 'Normalization','cdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N,'color',[0.6,0.6,1])
end
ylim([0,1])
xlim([3,8])

[N,edges] = histcounts(x2,50, 'Normalization','cdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N,'color',[0,0,1],'linewidth',2)
box on
title ('Cumulative distribution after 10 days','Interpreter','Latex')
xlabel('SMBG [mmol/L]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 12;
box on
grid on

x1 = x1(:);
n1 = numel(x1);

iTIR = find(x1>=3.9 & x1<=10);
iTBR = find(x1>=3.0 &  x1<3.9);
iTBR2 = find(x1<3.0);
iTAR = find(x1>10 & x1 <=13.9);
iTAR2 = find(x1>13.9);
iTIR2 = find(x1>=3.9 & x1<=6);

stats.TIR = numel(iTIR)/n1*100;
stats.TBR = numel(iTBR)/n1*100;
stats.TBR2 = numel(iTBR2)/n1*100;
stats.TAR = numel(iTAR)/n1*100;
stats.TAR2 = numel(iTAR2)/n1*100;
stats.TIR2 = numel(iTIR2)/n1*100;

% after 2 weeks
x2 = x2(:);
n2 = numel(x2);

iTIR = find(x2>=3.9 & x2<=10);
iTBR = find(x2>=3.0 &  x2<3.9);
iTBR2 = find(x2<3.0);
iTAR = find(x2>10 & x2 <=13.9);
iTAR2 = find(x2>13.9);
iTIR2 = find(x2>=3.9 & x2<=6);

stats2.TIR = numel(iTIR)/n2*100;
stats2.TBR = numel(iTBR)/n2*100;
stats2.TBR2 = numel(iTBR2)/n2*100;
stats2.TAR = numel(iTAR)/n2*100;
stats2.TAR2 = numel(iTAR2)/n2*100;
stats2.TIR2 = numel(iTIR2)/n2*100;

%% Functions 
function res = cost(x,Jhat)
res = -Jhat + (x-5)^2 + 8*(min(0,x-5))^2;
end

function [u,theta,P,gradient] = dose_guidance(u,SMBG,day,theta,P,lambda,adhere)
delta = 0.5*24;

if ~adhere
    u_1 = 0;
else
    u_1 = u;
end
J = (SMBG-5)^2 + 8*(min(0,SMBG-5))^2;

% RLS gradient estimation
phi = [u_1;1];
K = (lambda + phi'*P*phi)\(P*phi);
P = (P - (lambda + phi'*P*phi)\(P*phi*phi'*P))/lambda;
theta = theta + K*(J-phi'*theta);
gradient = theta(1);

if SMBG >7  && gradient >0.04
    Ki = -100;  % Do no decrease insulin if SMBG >7
else if SMBG < 4.2 && gradient<-0.04
        Ki = -100; % Do no increase insulin if SMBG <4
    else
        Ki = 1500;
    end
end


if SMBG < 6 && SMBG > 4
    roc = 2;
else
    roc = 8;
end

u = u - min(roc*24,max(-roc*24,Ki*gradient));

if rem(day,2)==0
    u_in = u + delta;
else
    u_in = u - delta;
end

end


function u = dose_202(u,SMBG,day)

if  rem(day,7)==0
    if SMBG <4
        u = u -2*24;
    else if SMBG>5
            u = u + 2*24;
        else
            u = u;
        end
    end
end

end

function u = dose_stepwise(u,SMBG_3,day)

if  SMBG_3(end)<4
    u = u -2*24;
else  if rem(day,7)==0
        
        bg = mean(SMBG_3);
        
        if bg <=3.1
            u = u -4*24;
        else if bg>=3.1 && bg <= 4
                u = u -2*24;
            else if bg>4 && bg<=5
                    u = u;
                else if bg>5 && bg <=7
                        u = u + 2*24;
                    else if bg>7 && bg <=8
                            u = u + 4*24;
                        else if bg>8 && bg <=9
                                u = u + 6*24;
                            else if bg >9
                                    u = u + 8*24;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


end


