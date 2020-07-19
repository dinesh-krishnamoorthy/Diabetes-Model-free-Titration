clear
clc
par.Ts = 1/24; % 1 hour

% Virtual patient parameters
par.tau1 = 0.5;
par.p2 = 15.8;     % mean value from Kanderian et al.
par.pGEZI = 3.31;  % mean value from Kanderian et al.

[f,sys] = MVP_T2D(par);

% Virtual patient initial state
d_in = [1.8;368;1.27];
u0  = 5*24; % Initial insulin
xfs = [0;0;15.2;12];

% setup
dT = 6/(60*24); % 6 min
Ns = 60/6; % 1hour/6min = 10 samples
k = 1;
init_insulin = 5; % day to start insulin titration

% RLS initialization
lambda = 0.01;
P = eye(2);
theta  = [0;0];
u(1) = u0;

method = 1;

rng(2)
adhere_1 = 1;
sim.u(k) = u0;
u_last = u0;
for sim_k = 1:1*90*24
    
    if  rem(sim_k,24)==0
        k = k+1;
        adhere =round(0.35 + 0.65*rand); % non-adherence to treatment regimen
       
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
            SMBG(k) = NaN;  
            u(k) = u(k-1);
            u_in = 0;
            adhere_1 = adhere;
            sim.u(k) = u_in;
            
        end
    else
        u_in=0;
    end
    
    % -----------  Simulator -----------
    if sim_k>60*24
        d_in(1) = 1.8*(1+0.3); % 30% increase in insulin sensitivity
    else
        d_in(1) = 1.8;
    end
    
    for i = 1:Ns
        % Stocahstic simulation - Euler Maruyama
        dw = randn(4,1);
        xfs = xfs + full(f(xfs,u_in,d_in))*dT + 2*sqrt(dT)*dw;
    end
    sim.G(sim_k) = xfs(4);
    sim.t(sim_k) = sim_k/24;
    
end

%%
days = (0:k-1);
figure(22)

set(gcf, 'Position',  [200, 200, 500, 400])
subplot(211)
hold all
fill([0 k-1 k-1 0],[4 4 6 6],[0.8 0.8 0.8],'facealpha',.35)
plot(days,SMBG,'rx','linewidth',1)
ylim([3,15])
ylabel('SMBG [mmol/L]','Interpreter','latex')
xlabel('days','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 12;
box on
grid on

subplot(212)
hold all
plot(days,sim.u./24,'ko','linewidth',1)
ylabel('Insulin [U]','Interpreter','latex')
xlabel('days','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 12;
box on
grid on

glyc_var = sum((SMBG-5).^2 + 8.*(min(0,SMBG-5)).^2);
avg_gl = mean(SMBG)
avg_30 = [mean(SMBG(1:30));mean(SMBG(31:end))];

%% ----------------- Functions --------------------

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


