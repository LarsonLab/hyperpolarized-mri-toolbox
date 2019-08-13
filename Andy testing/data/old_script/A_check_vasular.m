%% save
filename='JW318';
lac_v=lac_v';
pyr_v=pyr_v';
save(['./Andy testing/test_vascular/',filename(1:5)],'lac_c','lac_t','lac_v','pyr_c','pyr_t','pyr_v');
%% load
clearvars -except filename
%%
% filename='JW318';
load(['./Andy testing/test_vascular/',filename(1:5),'.mat'])
%% plot
figure
subplot(221);
	yyaxis left;plot(lac_v);hold on;plot(lac_t);
	yyaxis right;plot(lac_v./lac_t);
subplot(222);
	yyaxis left;plot(lac_v);hold on;plot(lac_c);
	yyaxis right;plot(lac_v./lac_c);
subplot(223);
	yyaxis left;plot(pyr_v);hold on;plot(pyr_t);
	yyaxis right;plot(pyr_v./pyr_t);
subplot(224);
	yyaxis left;plot(pyr_v);hold on;plot(pyr_c);
	yyaxis right;plot(pyr_v./pyr_c);