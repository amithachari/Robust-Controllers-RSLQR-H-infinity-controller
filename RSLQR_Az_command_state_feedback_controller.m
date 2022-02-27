clear all
close all

format short e

disp('****************************** Program Start ****************');
plot_file_name = 'stab_analysis1.ppt';
save_plots = 0; % Flag to bypass saving

V = 886.78;%fps
rtd  = 180/pi;
dt = 0.002;
t = 0:dt:2;
w = logspace(-1,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);

%Plant Model
Ma = 47.711;
Za_V = -1.3046;
Zd_V = -0.2142;
Md = -104.83;
Za = V*Za_V;
Zd = V*Zd_V;
w_act = 2.*pi*11;
z_act = 0.707;
grav = 32.174;

Ap =[Za_V   1.  Zd_V  0.;
     Ma     0.  Md  0;
     0.     0.  0.  1.;
     0.     0.  -w_act*w_act -2*z_act*w_act]
 
Bp =[0.; 0.; 0.; w_act*w_act]

Cp=[Za 0. Zd 0.;
    eye(4)]

Dp=[ 0*Cp*Bp]

pout(1,:) = 'Az (fps2)';
pout(2,:) = 'AOA (rad)';
pout(3,:) = 'q   (rps)';
pout(4,:) = 'de  (rad)';
pout(5,:) = 'dedt(rps)';

[nCp, nAp] = size(Cp);
[~, nBp] = size(Bp);

%Wiggle Model
Aw = [0. Za   0.;
      0. Za_V 1.;
      0. Ma   0.]
Bw = [Zd;
      Zd_V;
      Md]
  
%LQR Setup
Q = 0.*Aw;
R = 1;
xeig = [];
LQReig = [];

qq=logspace(-6,0,100);

xopenloop = eig(Aw);

az_st = 0.* ones(numel(qq), numel(t));
del_st = 0.* ones(numel(qq), numel(t));
deldot_st = 0.* ones(numel(qq), numel(t));
T_st = 0.* ones(numel(qq), numel(w));
S_st = 0.* ones(numel(qq), numel(w));

npts = numel(qq);
ii = 1;
i_cc = 0;

while ii < npts,
    i_cc = ii -1;
    Q(1,1) = qq(ii);
    [Kx_lqr, ~, ~] = lqr(Aw,Bw,Q,R);
    
    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = 0.;
    
    Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
    Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc)
        (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
    Bcl = [       Bp*Z*Dc2
        (Bc2+Bc1*Dp*Z*Dc2)];
    Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
    Dcl =(Dp*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);
    
    if max(real(eig(Acl))) > 0,
      disp('Closed-Loop System is Unstable');
      disp(['  Most unstable eig = ', num2str(max(real(eig(Acl))))]);
      return
    end;
    
    xx = eig(Acl);
    xeig = [xeig;xx];
    xx = eig(Aw - Bw*Kx_lqr);
    LQReig = [LQReig;xx];
    
    % Frequency Domain Analysis
    % SS model of loop gain at the plant input Lu
    A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
    B_Lu = [ Bp; Bc1*Dp];
    C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
    D_Lu = -[ Dc1*Dp];
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
    
    Lu = freqresp(sys_Lu,w);
    
    % SS model of loop gain at the plant input Ly
    Aout = [ Ap Bp*Cc; 0.*Bc1*Cp Ac];
    Bout = [ Bp*Dc1; Bc1];
    Cout = -[ Cp Dp*Cc];%change sign for loop gain
    Dout = -[ Dp*Dc1];
    sys_Ly = ss(Aout,Bout,Cout,Dout);
 %Analysis at Plant Input
    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    sr = sigma(sys_Lu,w,3); % Stability Robustness
    sru_min = min(abs(sr));
    rd = sigma(sys_Lu,w,2); % Return Difference
    rdu_min = min(abs(rd));
    
 %Analysis at Plant Output
    T = freqresp(sys_cl,w); % Complementary Sensitivity
    S = 1 - T; % Sensitivity
    T_st(ii,:) = 20*log10(abs(squeeze(T(1,1,:))));
    S_st(ii,:) = 20*log10(abs(squeeze(S(1,1,:))));
    Tmax = max(T_st(ii,:)); % Inf Norm of T in dB
    Smax = max(S_st(ii,:)); % Inf Norm of S in dB
    
    for i=1:numel(w),
        s = sqrt(-1)*w(i);
        GG = Cp*inv(s*eye(size(Ap))-Ap)*Bp+Dp;
        KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1;
        Lu_HS(i)  = -KK*GG;
        RDu_HS(i)  = 1. + Lu_HS(i);
        SRu_HS(i) = 1. + 1./Lu_HS(i);
        Lin = C_Lu*inv(s*eye(size(A_Lu))-A_Lu)*B_Lu+D_Lu;
        Lout = Cout*inv(s*eye(size(Aout)) - Aout)*Bout + Dout;
        for jj = 1:nBp,
            Fyjj = eye(size(Lin));
            Fyjj(jj,jj) = 0.;
            Tujj(jj,i) = inv(eye(size(Lin)) + Lin*Fyjj)*Lin;
        end
        for jj = 1:nCp,
            Fyjj = eye(nCp);
            Fyjj(jj,jj) = 0.;
            Tyjj = inv(eye(nCp) + Lout*Fyjj)*Lout;
            Tyj(jj,i) = Tyjj(jj,jj);
        end
    end
    for jj=1:nCp,
        v_min(jj) = min(min(abs(1 + Tyj(jj,:))));
    end
    
    if rdu_min <= 0.3,
        disp('Exit loop on rdmin')
        ip = ii -1;
        ii = npts;
    end
    
    % Time Domain Analysis
    y = step(sys_cl,t);
    az = y(:,1); % acceleration (fps2)
    aze = abs(ones(size(az))-az); % error for az
    taur = 0.; taus= 0.; % rise time and settling time
    fv = aze(numel(aze)); % final value of the error
    e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taur = crosst(e_n1,t); % rise time
    e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taus = crosst(e_n1,t); % settling time
    azmin = abs(min(az))*100; % undershoot
    azmax = (abs(max(az))-1)*100; % overshoot
    dmax = max(abs(y(:,4)))*rtd*grav; % compute in per g commanded
    ddmax = max(abs(y(:,5)))*rtd*grav;
    metric=[qq(ii) rdu_min sru_min wc taur taus azmin azmax dmax ddmax Tmax Smax];
    data(ii,:) = metric;
    az_st(ii,:) = az';
    del_st(ii,:) = rtd*y(:,4);
    deldot_st(ii,:) = rtd*y(:,5);
    ii = ii+1;

end
    
figure('Name','Accel Time Histories');
plot(t,az_st);
title('Acceleration Step Response');
xlabel('Time (sec)');
ylabel('Az (fps2)');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Surface Time Histories');
plot(t,del_st);
title('Surface Deflection to Step Command');
xlabel('Time (sec)');
ylabel('Del (deg)');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Surface Rate Time Histories');
plot(t,deldot_st);
title('Surface Rate tp Step Command');
xlabel('Time (sec)');
ylabel('Deldot (deg)');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name', 'Comp Sens T');
semilogx(w, T_st);
title('Comp Sens T');
xlabel('Frequency (rps)');
ylabel('Mag (dB)');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name', 'Sens S');
semilogx(w, S_st);
title('Sens S');
xlabel('Frequency (rps)');
ylabel('Mag (dB)');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end

%Root Locus Plot
figure('Name', 'LQR Root Locus');
plot(real(xopenloop),imag(xopenloop), 'rd', real(xeig),imag(xeig),'x','Linewidth',2);
title('RSLQR Root Locus');
xlabel('Real');
ylabel('Imaginary');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end


%extract the data from the metric matrix
qm = data(1:i_cc,1); %LQR Penalty
rdm = data(1:i_cc,2); %min Singular Value of RD
srm = data(1:i_cc,3); %min Singular value of SRM
wcm = data(1:i_cc,4); %LGCF
trm = data(1:i_cc,5); %Rise Time
tsm = data(1:i_cc,6); %Settling Time
usm = data(1:i_cc,7); %Undershoot
osm = data(1:i_cc,8); %Overshoot
fam = data(1:i_cc,9); %Max control position from step response
frm = data(1:i_cc,10); %Max fin rate from step response
Tinf = data(1:i_cc,11); %Inf Norm of Comp Sens T
Sinf = data(1:i_cc,12); %Inf Norm of Sens S


err = abs(rdm - srm);
[err_min, i_min] = min(err);
ip = i_min - 3

%Plot minimum values of the return difference and stability roboustness
%matrix
figure('Name','LQR Design Chart - RDM, SRM vs wc')
semilogx(wcm,rdm,'b-', wcm,srm,'r-','LineWidth',2)
hold on
semilogx(wcm(ip),rdm(ip),'bo', wcm(ip),srm(ip),'ro','LineWidth',2)
xlabel('wc (rps)');
ylabel('min(I + L) min(I + inv(L))');
legend('min RDM Lu', 'min SRM Lu', 'Location', 'Best');
if(save_plots == 1) saveppt2(plot_file_name); end
hold off
grid

%Plot of the q's vs bandwidth
figure('Name','LQR Design Chart - Penalty vs wc')
loglog(wcm,qm,wcm(ip),qm(ip),'bo','LineWidth',2);grid
text(.5,.5,num2str(qm(ip)));
xlabel('wc (rps)');
ylabel('LQR Penalty q11');
if(save_plots == 1) saveppt2(plot_file_name); end

%Plot of q's vs bandwidth
figure('Name','LQR Design Chart - Tinf, Sinf vs wc')
semilogx(wcm,Tinf,'r', wcm,Sinf,'b','LineWidth',2)
hold on
semilogx(wcm(ip),Tinf(ip),'ro', wcm(ip),Sinf(ip),'bo','LineWidth',2)
legend('inf norm Ty', 'inf norm Sy', 'Location', 'Best');
xlabel('wc (rps)');
ylabel('Inf Norm (dB)');
if(save_plots == 1) saveppt2(plot_file_name); end
hold off
grid

%Plot of rise time and settling time
figure('Name','LQR Design Chart - Tr, Ts vs wc')
semilogx(wcm,trm,'b-', wcm,tsm,'r-','LineWidth',2)
hold on
semilogx(wcm(ip),trm(ip),'bo', wcm(ip),tsm(ip),'ro','LineWidth',2)
xlabel('wc (rps)');
ylabel('Tr Ts (sec)');
if(save_plots == 1) saveppt2(plot_file_name); end
hold off
grid

%Plot of overshoot
figure('Name','LQR Design Chart - %OS vs wc')
semilogx(wcm,osm, wcm(ip),osm(ip),'bo','LineWidth',2); grid
xlabel('wc (rps)');
ylabel('% OS');
if(save_plots == 1) saveppt2(plot_file_name); end

%Plot of undershoot and overshoot
figure('Name','LQR Design Chart - %US and %OS vs wc')
semilogx(wcm,usm,'b', wcm,osm,'r','LineWidth',2)
hold on
semilogx(wcm(ip),usm(ip),'bo', wcm(ip),osm(ip),'ro','LineWidth',2)
xlabel('wc (rps)');
ylabel('%US %OS');
if(save_plots == 1) saveppt2(plot_file_name); end
hold off
grid

%Plot of undershoot
figure('Name','LQR Design Chart - %US vs wc')
semilogx(wcm,usm, wcm(ip),usm(ip),'bo','LineWidth',2); grid
xlabel('wc (rps)');
ylabel('% US');
if(save_plots == 1) saveppt2(plot_file_name); end

%Plot of Elevon Rate
figure('Name','LQR Design Chart - Deldot_max vs wc')
semilogx(wcm,frm, wcm(ip),frm(ip),'bo','LineWidth',2); grid
xlabel('wc (rps)');
ylabel('Elevon Rate deg/s/g');
if(save_plots == 1) saveppt2(plot_file_name); end

%Plot of Elevon
figure('Name','LQR Design Chart - Del_max vs wc')
semilogx(wcm,fam, wcm(ip),fam(ip),'bo','LineWidth',2); grid
xlabel('wc (rps)');
ylabel('Elevon deg/g');
if(save_plots == 1) saveppt2(plot_file_name); end


%Plot of the elevon rate
figure('Name','LQR Design Chart - Dele and Deledot vs wc')
semilogx(wcm,frm,'b', wcm,fam,'r','LineWidth',2)
hold on
semilogx(wcm(ip),frm(ip),'bo', wcm(ip),fam(ip),'ro','LineWidth',2)
xlabel('wc (rps)');
ylabel('Dele (deg) and Deledot (dps)');
if(save_plots == 1) saveppt2(plot_file_name); end
hold off
grid

%%Plot of Elevon
figure('Name','LQR Design Chart - Tr vs %US')
semilogx(usm,trm, usm(ip),trm(ip),'bo','LineWidth',2); grid
xlabel('Tr (sec)');
ylabel('% US');
if(save_plots == 1) saveppt2(plot_file_name); end

% linestyle = {'b-','r-','k-','g-','m-','c-'};
% linestyle2 = {'bo','ro','ko','go','mo','co'};
% figure('Name','SISO Min RDM All Input/Output')
% semilogx(wcm,data(1:i_cc,12+1), linestyle{1},'LineWidth',2); grid
% hold on
% for jj = 2:nCp,
%     semilogx(wcm, data(1:i_cc,12+jj), linestyle{1},'LineWidth',2); grid
% end
% semilogx(wcm,rdm,'b--','LineWidth',2); grid
% semilogx(wcm(ip),data(ip,12+1), linestyle2{1},'LineWidth',2); grid
% 
% for jj = 2:nCp,
%     semilogx(wcm(ip), data(ip,12+jj), linestyle2{jj},'LineWidth',2); grid
% end

semilogx(wcm(ip), rdm(ip),'bo', 'LineWidth',2); grid
hold off
xlabel('wc (rps)');
ylabel('Min RDM');
title('SISO Min RDM SISO Min RDM All Input/Output')
legend(pout(1,:),pout(2,:),pout(3,:),pout(4,:),pout(5,:),'RDM Lu')
grid
if(save_plots == 1) saveppt2(plot_file_name); end


Q(1,1) = qq(ip)
[Kx_lqr, Sw, Ew] = lqr(Aw, Bw,Q,R);
Acl = Aw - Bw*Kx_lqr;
[re,im] = nyquist(Aw, Bw, Kx_lqr, 0,1,w);
Lu = (re + sqrt(-1)*im);
SR_lqr = ones(size(Lu)) + ones(size(Lu))./Lu;
Lu_lqr = Lu;


%%%%
% Frequency Domain Analysis
    % SS model of loop gain at the plant input Lu
    
    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = 0.;
    
    A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
    B_Lu = [ Bp; Bc1*Dp];
    C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
    D_Lu = -[ Dc1*Dp];
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
    
    Lu = freqresp(sys_Lu,w);
    
    % SS model of loop gain at the plant input Ly
    Aout = [ Ap Bp*Cc; 0.*Bc1*Cp Ac];
    Bout = [ Bp*Dc1; Bc1];
    Cout = -[ Cp Dp*Cc];%change sign for loop gain
    Dout = -[ Dp*Dc1];
    sys_Ly = ss(Aout,Bout,Cout,Dout);
 %Analysis at Plant Input
    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    sr = sigma(sys_Lu,w,3); % Stability Robustness
    sru_min = min(abs(sr));
    rd = sigma(sys_Lu,w,2); % Return Difference
    rdu_min = min(abs(rd));
    
 %Analysis at Plant Output
    T = freqresp(sys_cl,w); % Complementary Sensitivity
    S = 1 - T; % Sensitivity
    T_st(ii,:) = 20*log10(abs(squeeze(T(1,1,:))));
    S_st(ii,:) = 20*log10(abs(squeeze(S(1,1,:))));
    Tmax = max(T_st(ii,:)); % Inf Norm of T in dB
    Smax = max(S_st(ii,:)); % Inf Norm of S in dB
    
    for i=1:numel(w),
        s = sqrt(-1)*w(i);
        GG = Cp*inv(s*eye(size(Ap))-Ap)*Bp+Dp;
        KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1;
        Lu_HS(i)  = -KK*GG;
        RDu_HS(i)  = 1. + Lu_HS(i);
        SRu_HS(i) = 1. + 1./Lu_HS(i);
        Lin = C_Lu*inv(s*eye(size(A_Lu))-A_Lu)*B_Lu+D_Lu;
        Lout = Cout*inv(s*eye(size(Aout)) - Aout)*Bout + Dout;
        for jj = 1:nBp,
            Fyjj = eye(size(Lin));
            Fyjj(jj,jj) = 0.;
            Tujj(jj,i) = inv(eye(size(Lin)) + Lin*Fyjj)*Lin;
        end
        for jj = 1:nCp,
            Fyjj = eye(nCp);
            Fyjj(jj,jj) = 0.;
            Tyjj = inv(eye(nCp) + Lout*Fyjj)*Lout;
            Tyj(jj,i) = Tyjj(jj,jj);
        end
    end
    for jj=1:nCp,
        v_min(jj) = min(min(abs(1 + Tyj(jj,:))));
    end
    
    if rdu_min <= 0.3,
        disp('Exit loop on rdmin')
        ip = ii -1;
        ii = npts;
    end


%%%%%

rho_Acl = max(abs(eig(Acl)));
accel_zero = tzero(Ap,Bp,Cp(1,:),Dp(1,:));

%Root Locus Plot
figure('Name', 'LQR Root Locus');
plot(real(xopenloop),imag(xopenloop), 'rX', real(LQReig(1:3*ip)),imag(LQReig(1:3*ip)),'x',real(accel_zero), imag(accel_zero),'rO',rho_Acl*xx1,rho_Acl*yy1,'g-','Linewidth',2);
axis([-30 10 -20 20]);
title('RSLQR Root Locus');
xlabel('Real');
ylabel('Imaginary');
legend('Open Loop Poles', 'Closed Loop Poles', 'Accel Zero', 'Max Spectral Radius', 'Location', 'Best');
grid;
axis;
if(save_plots == 1) saveppt2(plot_file_name); end

    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = 0.;
    
    Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
    Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc)
        (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
    Bcl = [       Bp*Z*Dc2
        (Bc2+Bc1*Dp*Z*Dc2)];
    Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
    Dcl =(Dp*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);

% Compute LGCF
bode_options                = bodeoptions;
% bode_options.FreqUnits      = 'Hz'; 
bode_options.PhaseWrapping  = 'on';
[mag,pha] = bode(sys_Lu,w);

[GM, PM_deg, wc_GM, wc_Pm] = margin(sys_Lu)

disp('Classical Margins')
allmargin(sys_Lu)

disp('  ')
disp('SV Margins')
RDu_nGM = 1/(1+rdu_min);
RDu_pGM = 1/(1-rdu_min);
RDu_Pha = 2*asin(rdu_min/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi ;
disp('RDu_nGM RDu_pGM RDu_Pha')
disp([num2str(RDu_nGM) ' ' num2str(RDu_pGM) ' ' num2str(RDu_Pha)])
disp([num2str(RDu_nGM_dB) ' ' num2str(RDu_pGM_dB) ' ' num2str(RDu_Pha_deg)])
SRu_nGM = 1-sru_min;
SRu_pGM = 1+sru_min;
SRu_Pha = 2*asin(sru_min/2);
SRu_nGM_dB = 20*log10(SRu_nGM);
SRu_pGM_dB = 20*log10(SRu_pGM);
SRu_Pha_deg = 180*SRu_Pha/pi ;
disp('SRu_nGM SRu_pGM SRu_Pha')
disp([num2str(SRu_nGM) ' ' num2str(SRu_pGM) ' ' num2str(SRu_Pha)])
disp([num2str(SRu_nGM_dB) ' ' num2str(SRu_pGM_dB) ' ' num2str(SRu_Pha_deg)])
disp('  ')


figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu)),imag(squeeze(Lu)),'k',real(Lu_HS),imag(Lu_HS),'r--',...
real(Tujj),imag(Tujj),'g--','LineWidth',2);grid
axis([-2 2 -2 2]);
legend('Unit Circle','Loop','Short Form','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

% Compute LGCF
bode_options                = bodeoptions;
% bode_options.FreqUnits      = 'Hz'; 
bode_options.PhaseWrapping  = 'on';
[mag,pha] = bode(sys_Lu,w);

% bode(sys_Lu)
grid
figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(Lu_HS)),'k',w,20*log10(squeeze(mag)),'r--','LineWidth',2);grid
legend('Loop','Short Form','Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode Magnitude at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Bode Phase at Plant Input'),
semilogx(w,squeeze(pha),'k','LineWidth',2);grid
% legend('Loop','Short Form','Location','Best');
xlabel('Frequency (rps)')
ylabel('Phase (deg)')
title('Bode Phase at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(RDu_HS)),'r--','LineWidth',2);grid
    RDu_min = min(min(abs(RDu_HS)));
    legend([num2str(jj) ' min(I+Lu) = ' num2str(RDu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end



figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(SRu_HS)),'r--','LineWidth',2);grid
    SRu_min = min(min(abs(SRu_HS)));
    legend([num2str(jj) ' min(I+invLu) = ' num2str(SRu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end


%Uncertainity Model
s = tf('s');
delta = ((w_act*w_act)/(s*s + 2*s*z_act*w_act + w_act*w_act)) - 1;
sv_delta = sigma(delta,w);

figure('Name','Robustness Check using Small Gain Theorem'),
semilogx(w,20*log10(abs(SRu_HS)),'r','LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(sv_delta)),'k','LineWidth',2);grid
SRu_min = min(min(abs(SRu_HS)));
legend([num2str(jj) ' min(I+invLu) = ' num2str(SRu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end
grid


[cl_EigVec,cl_eigVal] = eig(Acl)
