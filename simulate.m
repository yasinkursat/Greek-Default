clc
clear all
load data_sim.txt
load def_per.txt
load def_per_ss.txt
load param.txt
delta = 0.155;
indicator_CCD = 0;
r=0.04;
coupon = 1 - exp(-r)*(1-delta);
delta_LC = delta;
coupon_LC = delta_LC;
per_num = param(1);  %HOW MANY PERIODS IN EACH SAMPLE
n = param(2);        %HOW MANY SAMPLES
y = zeros(per_num, n);
b_LC = zeros(per_num, n);
b = zeros(per_num, n);
q_LC = zeros(per_num, n);
q = zeros(per_num, n);
c = zeros(per_num, n);
tb = zeros(per_num, n);
d = zeros(per_num, n);
excl = zeros(per_num, n);
duration = zeros(per_num,n);
i_g_sim = zeros(per_num, n);
i_premium_sim = zeros(per_num, n);
def_prob = zeros(per_num,n);
pi = zeros(per_num,n);
q_LC_covenant = zeros(per_num, n);
if indicator_CCD > .5
    for i=1:n
       y(:,i) = data_sim((i-1)*per_num+1:i*per_num,1); 
       b_LC(:,i) = data_sim((i-1)*per_num+1:i*per_num,2); 
       b(:,i) = data_sim((i-1)*per_num+1:i*per_num,3); 
       q(:,i) = data_sim((i-1)*per_num+1:i*per_num,4);
       q_LC(:,i) = data_sim((i-1)*per_num+1:i*per_num,5);
       c(:,i) = data_sim((i-1)*per_num+1:i*per_num,6); 
       tb(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
       d(:,i) = data_sim((i-1)*per_num+1:i*per_num,8)-1;
       excl(:,i) = data_sim((i-1)*per_num+1:i*per_num,9);
       i_premium_sim(:,i) = data_sim((i-1)*per_num+1:i*per_num,10);
       duration(:,i) = (((1-delta).*q(:,i) +coupon)./(q(:,i)))./(((1-delta).*q(:,i) +coupon)./(q(:,i))-1+delta);
       duration_LC(:,i) = (((1-delta_LC).*q_LC(:,i)+ coupon_LC)./((1+pi(:,i)).*coupon_LC));
       pi(:,i) = data_sim((i-1)*per_num+1:i*per_num,15);
    end
    q_LC_covenant = q_LC;
else
    for i=1:n
       y(:,i) = data_sim((i-1)*per_num+1:i*per_num,1); 
       b_LC(:,i) = data_sim((i-1)*per_num+1:i*per_num,2); 
       b(:,i) = data_sim((i-1)*per_num+1:i*per_num,3); 
       q(:,i) = data_sim((i-1)*per_num+1:i*per_num,4);
       q_LC(:,i) = data_sim((i-1)*per_num+1:i*per_num,5);
       c(:,i) = data_sim((i-1)*per_num+1:i*per_num,6); 
       tb(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
       d(:,i) = data_sim((i-1)*per_num+1:i*per_num,8)-1;
       excl(:,i) = data_sim((i-1)*per_num+1:i*per_num,9);
       i_premium_sim(:,i) = data_sim((i-1)*per_num+1:i*per_num,10);
       duration(:,i) = (((1-delta).*q(:,i) +coupon)./(q(:,i)))./(((1-delta).*q(:,i) +coupon)./(q(:,i))-1+delta);
       duration_LC(:,i) = (((1-delta_LC).*q_LC(:,i)+ coupon_LC)./((1+pi(:,i)).*coupon_LC));
       pi(:,i) = data_sim((i-1)*per_num+1:i*per_num,15);
       q_LC_covenant(:,i) = data_sim((i-1)*per_num+1:i*per_num,16);
    end   
end
bnext = b(2:per_num,:);
%CREATE OUTPUT SERIES
num = 20; %HOW MANY PERIODS IN EACH SUBSAMPLE
indices = find(sum(d(per_num - num - 5 +1:per_num,:))<1); 


num_observations = length(indices);
last_y = y(per_num - num+1:per_num, indices);  
last_c = c(per_num - num+1:per_num, indices);  
last_tb = tb(per_num - num+1:per_num, indices);  
last_d = d(per_num - num+1:per_num, indices);
last_b = b(per_num - num+1:per_num, indices);
last_bnext= last_b(2:num,:);
last_b_LC = b_LC(per_num - num+1:per_num, indices);
last_b_LCnext= last_b_LC(2:num,:);
last_i_premium = i_premium_sim(per_num - num+1:per_num, indices);
last_i_g = i_g_sim(per_num - num+1:per_num, indices);
last_def_prob = def_prob(per_num - num+1:per_num, indices);
last_spread = (log(coupon./q(per_num - num+1:per_num,indices) - delta + 1)-r);
last_spread_LC = (log(coupon./q_LC(per_num - num+1:per_num,indices) - delta + 1)-r);

last_spread_LC_covenant = (log(coupon./q_LC(per_num - num+1:per_num,indices) - delta + 1)-log(coupon./q_LC_covenant(per_num - num+1:per_num,indices) - delta + 1));
last_pi = pi(per_num - num+1:per_num, indices);
last_duration = duration(per_num - num+1:per_num, indices);
last_totb = last_b + last_b_LC;
last_totby = (b(per_num - num+1:per_num, indices)+b_LC(per_num - num+1:per_num, indices))./exp(y(per_num - num+1:per_num, indices));
last_LC_share = b_LC(per_num - num+1:per_num, indices)./(b_LC(per_num - num+1:per_num, indices)+b(per_num - num+1:per_num, indices));
last_duration_LC = duration_LC(per_num - num+1:per_num, indices);
last_bFXy = b(per_num - num+1:per_num, indices)./exp(y(per_num - num+1:per_num, indices));
last_bLCy = b_LC(per_num - num+1:per_num, indices)./exp(y(per_num - num+1:per_num, indices));
last_LC_sharey = last_LC_share./exp(y(per_num - num+1:per_num, indices));
issuance_FX = last_bnext - (1-delta).*last_b(1:num-1,:);
issuance_LC = last_b_LCnext - (1-delta).*last_b_LC(1:num-1,:);

%CREATE VECTORS OF STD OF RETURNS AND CORRELATION BETWEEN RETURNS, OUTPUT
%AND TB. NEED TO TRIM OBSERVATIONS WHILE THE COUNTRY IS IN DEFAULT.

lambda = 100;

y_trend = hpfilter(last_y,lambda);            
c_trend = hpfilter(last_c, lambda);          
tb_trend = hpfilter(last_tb, lambda);         
spread_trend = hpfilter(last_spread, lambda); 
pi_trend = hpfilter(last_pi, lambda);
spread_LC_trend = hpfilter(last_spread_LC_covenant, lambda); 
b_trend = hpfilter(last_b, lambda);
b_LC_trend = hpfilter(last_b_LC, lambda);
issuance_FX_trend = hpfilter(issuance_FX, lambda);
issuance_LC_trend = hpfilter(issuance_LC, lambda);
b_next_trend = hpfilter(last_bnext, lambda);
bLC_next_trend = hpfilter(last_b_LCnext, lambda);
bFXy_trend = hpfilter(last_bFXy,lambda);
bLCy_trend = hpfilter(last_bLCy,lambda);
totby_trend = hpfilter(last_totby, lambda);
LC_sharey_trend = hpfilter(last_LC_sharey, lambda);
LC_share_trend = hpfilter(last_LC_share, lambda);

y_dev = last_y - y_trend';
c_dev = last_c - c_trend';
tb_dev = last_tb - tb_trend';
spread_dev = last_spread - spread_trend';
spread_LC_dev = last_spread_LC_covenant - spread_LC_trend';

pi_dev = last_pi - pi_trend';
b_dev = last_b - b_trend';
b_LC_dev = last_b_LC - b_LC_trend';
issuance_FX_dev = issuance_FX - issuance_FX_trend';
issuance_LC_dev = issuance_LC - issuance_LC_trend';
bFX_next_dev = last_bnext - b_next_trend';
bLC_next_dev = last_b_LCnext - bLC_next_trend';
bFXy_dev = last_bFXy - bFXy_trend';
bLCy_dev = last_bLCy - bLCy_trend';
totby_dev = last_totby - totby_trend';
LC_sharey_dev = last_LC_sharey - LC_sharey_trend';
LC_share_dev = last_LC_share - LC_share_trend';

for i=1:num_observations
    matrix = corrcoef([y_dev(:,i) c_dev(:,i) tb_dev(:,i) spread_dev(:,i) bFXy_dev(:,i) bLCy_dev(:,i) spread_LC_dev(:,i) LC_sharey_dev(:,i) totby_dev(:,i)]);
    matrix_corr(i,:) = matrix(1,2:9);

    %COMPUTE CORRELATION BETWEEN FILTERED SPREAD AND TRADE BALANCE
    matrix1 = corrcoef([spread_dev(:,i) spread_LC_dev(:,i)]);
    matrix_corr1(i) = matrix1(1,2);
    
    matrix_pi = corrcoef([pi_dev(:,i) y_dev(:,i) spread_LC_dev spread_dev(:,i) b_dev(:,i) b_LC_dev(:,i) bFXy_dev(:,i) bLCy_dev(:,i) c_dev(:,i) LC_sharey_dev(:,i) totby_dev(:,i)]);
    matrix_corr_pi(i,:) = matrix_pi(1,2:11);
    
    matrix_issuance = corrcoef([pi_dev(2:num,i) issuance_FX_dev(1:num-1,i) issuance_LC_dev(1:num-1,i)]);
    matrix_corr_issuance(i,:) = matrix_issuance(1,2:3);
    
    matrix_issuance_y = corrcoef([y_dev(2:num,i) issuance_FX_dev(1:num-1,i) issuance_LC_dev(1:num-1,i)]);
    matrix_corr_issuance_y(i,:) = matrix_issuance_y(1,2:3);
    
    matrix_nextb = corrcoef([pi_dev(1:num-1,i) bFX_next_dev(:,i) bLC_next_dev(:,i)]);
    matrix_corr_nextb(i,:) = matrix_nextb(1,2:3);
    
    matrix_LCshare = corrcoef([LC_share_dev(1:num,i) spread_LC_dev(:,i) spread_dev(:,i) y_dev(:,i) pi_dev(:,i)]);
    matrix_corr_LCshare(i,:) = matrix_LCshare(1,2:5);
    
%     matrix_issuance_not_detrended = corrcoef([last_pi(2:num,i) issuance_FX(1:num-1,i) issuance_LC(1:num-1,i)]);
%     matrix_corr_issuance_not_detrended(i,:) = matrix_issuance_not_detrended(1,2:3);
    
end
fprintf(['Statistics when default-free ',num2str(num_observations),' samples of ',num2str(num),' periods are considered \n'])

fprintf('corr(c, y)   = %6.2f \n',mean(matrix_corr(:,1)))
fprintf('corr(tb,y)   = %6.2f \n',mean(matrix_corr(:,2)))
fprintf('corr(Rs FX,y)   = %6.2f \n',mean(matrix_corr(:,3)))
fprintf('corr(bFXy,y)   = %6.2f \n',mean(matrix_corr(:,4)))
fprintf('corr(bLCy,y)   = %6.2f \n',mean(matrix_corr(:,5)))
fprintf('corr(Rs LC,y)   = %6.2f \n',mean(matrix_corr(:,6)))
fprintf('corr(LC share,y)   = %6.2f \n',mean(matrix_corr(:,7)))
fprintf('tot debt,y)   = %6.2f \n',mean(matrix_corr(:,8)))
fprintf('corr(RsFX,RsLC)   = %6.2f \n',mean(matrix_corr1(1,:)))
fprintf('corr(pi, y)   = %6.2f \n',mean(matrix_corr_pi(:,1)))
fprintf('corr(pi, RsLC)   = %6.2f \n',mean(matrix_corr_pi(:,2)))
fprintf('corr(pi, RsFX)   = %6.2f \n',mean(matrix_corr_pi(:,3)))
fprintf('corr(bFX, pi)   = %6.2f \n',mean(matrix_corr_pi(:,4)))
fprintf('corr(bLC, pi)   = %6.2f \n',mean(matrix_corr_pi(:,5)))
fprintf('corr(bFXy, pi)   = %6.2f \n',mean(matrix_corr_pi(:,6)))
fprintf('corr(bLCy, pi)   = %6.2f \n',mean(matrix_corr_pi(:,7)))
fprintf('corr(c, pi)   = %6.2f \n',mean(matrix_corr_pi(:,8)))
fprintf('corr(LC share/y, pi)   = %6.2f \n',mean(matrix_corr_pi(:,9)))
fprintf('corr(tot debt/y, pi)   = %6.2f \n',mean(matrix_corr_pi(:,10)))
fprintf('corr(bFX issuance, pi)   = %6.2f \n',mean(matrix_corr_issuance(:,1)))
fprintf('corr(bLC issuance, pi)   = %6.2f \n',mean(matrix_corr_issuance(:,2)))
fprintf('corr(bFX issuance, y)   = %6.2f \n',mean(matrix_corr_issuance_y(:,1)))
fprintf('corr(bLC issuance, y)   = %6.2f \n',mean(matrix_corr_issuance_y(:,2)))
fprintf('corr(bFX next, pi)   = %6.2f \n',mean(matrix_corr_nextb(:,1)))
fprintf('corr(bLC next, pi)   = %6.2f \n',mean(matrix_corr_nextb(:,2)))
fprintf('corr(LC share, Rs LC)   = %6.2f \n',mean(matrix_corr_LCshare(:,1)))
fprintf('corr(LC share, Rs)   = %6.2f \n',mean(matrix_corr_LCshare(:,2)))
fprintf('corr(LC share, y)   = %6.2f \n',mean(matrix_corr_LCshare(:,3)))
fprintf('corr(LC share, pi)   = %6.2f \n',mean(matrix_corr_LCshare(:,4)))
fprintf('std of output  = %6.2f \n',100*mean(std(y_dev)))
fprintf('std cons / std y   = %6.2f \n',mean(std(c_dev))/mean(std(y_dev)))
fprintf('std of TB/Y    = %6.2f \n',100*mean(std(tb_dev)))
fprintf('std of R_s FX     = %6.2f \n',100*mean(std(spread_dev)))
fprintf('std of R_s LC     = %6.2f \n',100*mean(std(spread_LC_dev)))
fprintf('std of pi     = %6.2f \n',100*mean(std(pi_dev)))
fprintf('std of bFXy     = %6.2f \n',100*mean(std(bFXy_dev)))
fprintf('std of bLCy     = %6.2f \n',100*mean(std(bLCy_dev)))
fprintf('std of totby     = %6.2f \n',100*mean(std(totby_dev)))

ind = find(data_sim(:,9)<2); % find periods in which the government is not in default
fprintf('Annual mean default rate (all periods) = %6.2f \n', 100*sum(data_sim(:,8)-1)/length(ind))
fprintf('Mean FX debt/y    = %6.2f \n',100*mean(mean(last_b./(exp(last_y)))))
fprintf('Mean LC debt/y    = %6.2f \n',100*mean(mean(last_b_LC./(exp(last_y)))))
fprintf('Total debt/y    = %6.2f \n',100*mean(mean(last_b./(exp(last_y))))+100*mean(mean(last_b_LC./(exp(last_y)))))

fprintf('E(R_s)         = %6.2f \n',100*mean(mean(last_spread)))
fprintf('E(R_s LC covenant)         = %6.2f \n',100*mean(mean(last_spread_LC_covenant)))
fprintf('inflation         = %6.2f \n',100*mean(mean(last_pi)))
fprintf('share of purchase back for FX debt      = %6.2f \n',100*(length(find(last_bnext - (1 - delta)*last_b(1:num-1,:)<-0.01)))/(length(last_bnext)*(num-1)))
fprintf('share of purchase back for LC debt     = %6.2f \n',100*(length(find(last_b_LCnext - (1 - delta)*last_b_LC(1:num-1,:)./(1+last_pi(1:num-1,:))<-0.01)))/(length(last_b_LCnext)*(num-1)))
fprintf('Duration     = %f5 \n',mean(mean(last_duration)))

if indicator_CCD < .5
    column = 3;
    m = length(last_y(1,:)); %150;
    figure
    plot(exp(last_y(:,1:m)), last_pi(:,1:m)*100, 'LineStyle','none', 'Marker','o','Color','b','Markersize',4);
    hold on
    plot(exp(last_y(:,column)), last_pi(:,column)*100, 'LineStyle','none', 'markerfacecolor','k','Marker','d','Color','k','Markersize',7);
    xlabel('Income','FontSize',16, 'FontWeight','bold','Interpreter','latex')
    ylabel('Inflation (annual, \%)','FontSize',16, 'FontWeight','bold','Interpreter','latex')
    hold on
    x = xlim;
    plot([x(1) x(2)],[7.17 7.17], 'LineStyle','--', 'Color','k','LineWidth',2,'Markersize',10)
    dim = [.15 .33 .0 .0];
    str = 'mean inflation';
    annotation('textbox', dim, 'String', str,'FitBoxtoText','on','EdgeColor','none')
    set(gca,'TickLabelInterpreter','latex','Fontsize',16)
    hold off
end
