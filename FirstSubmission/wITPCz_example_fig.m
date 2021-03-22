

%% Randomly generate values
randvects = rand(50,1)*2*pi-pi;
vectormod = (randvects + randn(size(randvects))).^2;
vectormod2 = vectormod-min(vectormod)+1; % make sure no negative values



%% ITPC
% subplot(2,2,1)
figure
% polar([zeros(size(randvects)) randvects]',[zeros(size(randvects)) ones(size(randvects))]','k')
polarplot([zeros(size(randvects)) randvects]',[zeros(size(randvects)) ones(size(randvects))]','k')
ax = gca;
ax.GridLineStyle = '--';
title([ 'ITPC = ' num2str(round(abs(mean(exp(1i*randvects))),4)) ])

%% wITPC
% subplot(2,2,3)
figure
% polar([zeros(size(randvects)) randvects]',[zeros(size(randvects)) vectormod]','k')
polarplot([zeros(size(randvects)) randvects]',[zeros(size(randvects)) vectormod]','k')
ax = gca;

witpc_fake = abs(mean(vectormod.*exp(1i*randvects)));

% wITPCz
perm_witpc = zeros(1,1000);
for i=1:1000
    perm_witpc(i) = abs(mean(vectormod(randperm(length(vectormod))).*exp(1i*randvects)));
end

witpc_z_fake = (witpc_fake-mean(perm_witpc))./std(perm_witpc);
title([ 'wITPCz = ' num2str(round(witpc_z_fake,4)) ])

%% Histogram of null-hypothesis WITPC
% subplot(2,2,[2,4])
figure
[y,x]=hist(perm_witpc,50);
h=bar(x,y,'histc','EdgeColor',[1 1 1]);
set(h,'linestyle','-')
hold on
plot([witpc_fake witpc_fake],get(gca,'ylim')/2,'m')
title('Histogram of H_0 wITPC values')
ylabel('Count'); xlabel('wITPC');
% title([ 'wITPC = ' num2str(witpc) ])


clear x y h i 


clear witpc_fake witpc_z_fake perm_witpc randvects vectormod vectormod2





