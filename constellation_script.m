%% read in data
data=xlsread('Stars.xlsx');
RA_all = data(:,2);
HR_all = data(:,1);
DEC_all = data(:,3);
appmag_all = data(:,4);
real_partition = data(:,5);
N_all = length(HR_all);
stars_lines=xlsread('stars_linesonly.xlsx');
HR_lines = stars_lines(:,1);
N_lines = length(HR_lines);
RA_lines = zeros(N_lines,1);
DEC_lines = zeros(N_lines,1);
mag_lines = zeros(N_lines,1);
partition_lines = zeros(N_lines,1);
for i = 1:N_lines
    for j = 1:N_all
        if HR_lines(i) == HR_all(j)
            RA_lines(i) = RA_all(j);
            DEC_lines(i) = DEC_all(j);
            mag_lines(i) =appmag_all(j);
            partition_lines(i) = real_partition(j);
        end
    end
end
shiftedappmags = xlsread('Star_mags_shifted.xlsx','F2:F9082');

%% open or make finalscalednew If you make it, the numbers will be slightly different
figure; imshow('saccadegraph.jpg');
% extract saccade probabilities from graph by clicking on 11 points on the black line of the graph 
[x,y]=ginput(11);
% click on bottom and top of y axis
[xyaxis,yyaxis]=ginput(2)
% calculate information about the y axis on graph to put x and y into
% context
range_yaxis=range(yyaxis);
yaxisinterval=(range_yaxis)/22
max_yyaxis=max(yyaxis);
% scale y coordinates
ycoorscaled = zeros(11,1)
ycoorscaled(1,1) = (max_yyaxis)-y(1,1);
ycoorscaled(2,1) = (max_yyaxis)-y(2,1);
ycoorscaled(3,1) = (max_yyaxis)-y(3,1);
ycoorscaled(4,1) = (max_yyaxis)-y(4,1);
ycoorscaled(5,1) = (max_yyaxis)-y(5,1);
ycoorscaled(6,1) = (max_yyaxis)-y(6,1);
ycoorscaled(7,1) = (max_yyaxis)-y(7,1);
ycoorscaled(8,1) = (max_yyaxis)-y(8,1);
ycoorscaled(9,1) = (max_yyaxis)-y(9,1);
ycoorscaled(10,1) = (max_yyaxis)-y(10,1);
ycoorscaled(11,1) = (max_yyaxis)-y(11,1);
finalscalednew=ycoorscaled./yaxisinterval;
% make coordinates add to 100
add = sum(finalscalednew);
diff = (100-add)./11
finalscalednew = finalscalednew+diff

%% make constellations_indices based on real_partition
constellation_indices = zeros(88,2);
constellation_indices(1,1) = 1;
constellation_indices(88,2) = 9081;
for i = 1:88
    count = 0;
    for j = constellation_indices(i,1):9081
        if real_partition(j) == i;
            count = count+1;
        else
            constellation_indices(i,2)=constellation_indices(i,1)+count-1;
            constellation_indices(i+1,1) = j;
            break
        end
    end
end
%% make bycon, where 1=in line-drawing and 0=out of line-drawing. this is
% used for bycon sequence t/p-value analysis
inorout_lines = zeros(9081,1);
for i = 1:9081
    for j = 1:696
        if HR_all(i) == HR_lines(j)
            inorout_lines(i) = 1;
        end
    end
end
% by con cells
bycon_1or0 = {1,88};
for i = 1:88
    temp = inorout_lines(constellation_indices(i,1):constellation_indices(i,2));
    bycon_1or0{i} = temp;
end 
% by con vector
bycon_1or0_vector = zeros(9081,1);
index = 1;
for i = 1:88
    temp = bycon_1or0{i};
    for j = 1:length(temp)
        bycon_1or0_vector(index) = temp(j);
        index = index+1;
    end
end

%% create the distance matrix
distance_real = zeros(N_all,N_all);
for i = 1:N_all
    for j = (i+1):(N_all)
        distance_real(i,j) = acosd(((sind(DEC_all(i,1)))*(sind(DEC_all(j,1))))+((cosd(DEC_all(i,1)))*(cosd(DEC_all(j,1)))*(cosd((RA_all(i,1))-(RA_all(j,1)))))); 
    end
end
distance_real = distance_real'+distance_real;
%% make null voronoi model
voronoi_partition = zeros(9081,1000);
voronoi_a = zeros(1000,1);
for h=1:1000
    RA_voronoi = 360*rand(88,1);
    DEC_voronoi = 2*rand(88,1)-1;
    DEC_voronoi = asind(DEC_voronoi);
    for i = 1:9081
        min_dis = acosd(((sind(DEC_voronoi(1)))*(sind(DEC_all(i))))+((cosd(DEC_voronoi(1)))*(cosd(DEC_all(i)))*(cosd((RA_voronoi(1))-(RA_all(i))))));
        for j = 1:88
            distance = acosd(((sind(DEC_voronoi(j)))*(sind(DEC_all(i))))+((cosd(DEC_voronoi(j)))*(cosd(DEC_all(i)))*(cosd((RA_voronoi(j))-(RA_all(i))))));
            if distance < min_dis
                min_dis = distance;
                voronoi_partition(i,h) = j;
            elseif distance == min_dis && j ==1
                voronoi_partition(i,h) = 1;
            end
        end
    end
    voronoi_a(h) = find_a(voronoi_partition(:,h),real_partition);
end
mean_v = mean(voronoi_a);
std_v = std(voronoi_a);
%% make real transitions probability matricies
%% level 1: distance model with 1/distance matrix (with smallest distance capped at 1)
capped_distance = distance_real;
for i = 1:N_all
    for j = 1:N_all
        if capped_distance(i,j) < 1.0
            capped_distance(i,j) = 1;
        end
    end
end
oneoverdistance = 1./capped_distance;
for i = 1:9081
    oneoverdistance(i,i) = 0;
end
for i = 1:9081
    oneoverdistance(i,:) = oneoverdistance(i,:)./sum(oneoverdistance(i,:));
end
%% level 2: saccade model with saccade probabilities and distance
saccade_trans= zeros(N_all,N_all); % replace all distance values with corresponding saccade probability
for i=(1:N_all)
    for j = 1:N_all
        if (distance_real(i,j) >= 0.0 && distance_real(i,j) < 2.0)
            saccade_trans(i,j) = finalscalednew(1);
        elseif (distance_real(i,j) >= 2.0 && distance_real(i,j) < 4.0)
            saccade_trans(i,j) = finalscalednew(2);
        elseif (distance_real(i,j) >= 4.0 && distance_real(i,j) < 6.0)
            saccade_trans(i,j) = finalscalednew(3);
        elseif (distance_real(i,j) >= 6.0 && distance_real(i,j) < 8.0)
            saccade_trans(i,j) = finalscalednew(4);
        elseif (distance_real(i,j) >= 8.0 && distance_real(i,j) < 10.0)
            saccade_trans(i,j) = finalscalednew(5);
        elseif (distance_real(i,j) >= 10.0 && distance_real(i,j) < 12.0)
            saccade_trans(i,j) = finalscalednew(6);
        elseif (distance_real(i,j) >= 12.0 && distance_real(i,j) < 14.0)
            saccade_trans(i,j) = finalscalednew(7);
        elseif (distance_real(i,j) >= 14.0 && distance_real(i,j) < 16.0)
            saccade_trans(i,j) = finalscalednew(8);
        elseif (distance_real(i,j) >= 16.0 && distance_real(i,j) < 18.0)
            saccade_trans(i,j) = finalscalednew(9);
        elseif (distance_real(i,j) >= 18.0 && distance_real(i,j) < 20.0)
            saccade_trans(i,j) = finalscalednew(10);
        elseif (distance_real(i,j) >= 20.0) 
            saccade_trans(i,j) = finalscalednew(11);
        else
            saccade_trans(i,j)=0;
        end
    end
end
for i = 1:N_all
    saccade_trans(i,i) = 0; % make diagonal=0
end
% divide each probability by the frequency of that same probability within the row
saccade_trans_prop=zeros(N_all,N_all);
for i=1:N_all
    saccade_sums_byrow = zeros(11,1);
    saccade_sums_byrow(1) = (sum(saccade_trans(i,:)==finalscalednew(1)));
    saccade_sums_byrow(2) = (sum(saccade_trans(i,:)==finalscalednew(2)));
    saccade_sums_byrow(3) = (sum(saccade_trans(i,:)==finalscalednew(3)));
    saccade_sums_byrow(4) = (sum(saccade_trans(i,:)==finalscalednew(4)));
    saccade_sums_byrow(5) = (sum(saccade_trans(i,:)==finalscalednew(5)));
    saccade_sums_byrow(6) = (sum(saccade_trans(i,:)==finalscalednew(6)));
    saccade_sums_byrow(7) = (sum(saccade_trans(i,:)==finalscalednew(7)));
    saccade_sums_byrow(8) = (sum(saccade_trans(i,:)==finalscalednew(8)));
    saccade_sums_byrow(9) = (sum(saccade_trans(i,:)==finalscalednew(9)));
    saccade_sums_byrow(10) = (sum(saccade_trans(i,:)==finalscalednew(10)));
    saccade_sums_byrow(11) = (sum(saccade_trans(i,:)==finalscalednew(11)));
    for j = 1:N
        if saccade_trans(i,j) == finalscalednew(1)
    saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(1);
        elseif saccade_trans(i,j) == finalscalednew(2)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(2);
        elseif saccade_trans(i,j) == finalscalednew(3)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(3);
        elseif saccade_trans(i,j) == finalscalednew(4)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(4);
        elseif saccade_trans(i,j) == finalscalednew(5)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(5);
        elseif saccade_trans(i,j) == finalscalednew(6)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(6);
        elseif saccade_trans(i,j) == finalscalednew(7)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(7);
        elseif saccade_trans(i,j) == finalscalednew(8)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(8);
        elseif saccade_trans(i,j) == finalscalednew(9)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(9);
        elseif saccade_trans(i,j) == finalscalednew(10)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(10);
        elseif saccade_trans(i,j) == finalscalednew(11)
            saccade_trans_prop(i,j)=saccade_trans(i,j)./saccade_sums_byrow(11);
        else
            saccade_trans_prop(i,j)=0;
        end
    end
end
for i=1:N_all % make each row sum to 1
    A=sum(saccade_trans_prop(i,:));
    saccade_trans_prop(i,:)=saccade_trans_prop(i,:)./A;
end
%% level 3: magnitude model with apparent magnitudes, saccade probabilities, and distance
r = 10;
x = 0.05;
for i = 1:9081
    count = 0;
    local_stars = [];
    bright = [];
    for j = 1:9081
        if distance_real(i,j) < r
            count = count+1;
            local_stars(count) = shiftedappmags(j);
        end
    end
    A = ceil(x.*length(local_stars));
    bright = mink(local_stars, A);
    m_avg = mean(bright);
    saccade_mags_trans(:,i) = (10^(-0.4*(shiftedappmags(i) - m_avg)))*saccade_trans_prop(:,i);
end
for i=1:N_all
    saccade_mags_trans(i,:)=(saccade_mags_trans(i,:))./(sum(saccade_mags_trans(i,:)));
end

%% clustering by constellation are on real transition probability matricies
% level 1 - distance model
ood_ahat = human_expectations(oneoverdistance, 0.3);
% make matrix symmetrical for spectral clustering
ood_ahat = (ood_ahat+ood_ahat')./2;
ood_partition = zeros(9081,100);
ood_a = zeros(100,1);
ood_zscore = zeros(100,1);
for i = 1:100
    ood_partition(:,i) = spectralcluster(ood_ahat, 88, 'Distance','precomputed', 'ClusterMethod', 'kmedoids');
    ood_a(i) = find_a(ood_partition(:,i), real_partition);
    ood_zscore(i) = (ood_a(i) - mean_v)./std_v;
end
% level 2 - saccade model
sacc_ahat = human_expectations(saccade_trans_prop, 0.3);
sacc_ahat = (sacc_ahat+sacc_ahat')./2;
sacc_partition = zeros(9081,100);
sacc_a = zeros(100,1);
sacc_zscore = zeros(100,1);
for i = 1:100
    sacc_partition(:,i) = spectralcluster(sacc_ahat, 88, 'Distance','precomputed', 'ClusterMethod', 'kmedoids');
    sacc_a(i) = find_a(sacc_partition(:,i), real_partition);
    sacc_zscore(i) = (sacc_a(i) - mean_v)./std_v;
end
% level 3 - magnitude model
mags_ahat = human_expectations(saccade_mags_trans, 0.3);
mags_ahat = (mags_ahat+mags_ahat')./2;
mags_partition = zeros(9081,100);
mags_a = zeros(100,1);
mags_zscore = zeros(100,1);
for i = 1:100
    mags_partition(:,i) = spectralcluster(mags_ahat, 88, 'Distance','precomputed', 'ClusterMethod', 'kmedoids');
    mags_a(i) = find_a(mags_partition(:,i), real_partition);
    mags_zscore(i) = (mags_a(i) - mean_v)./std_v;
end

%% clustering by constelation line-drawings on real transition probability matricies
%% level 1 real t- and p-vals
% 100 different sequences from ahat real distance model
sequ_by_con_ood = zeros(1500,88,100);
for i = 1:100
    for j = 1:88
        sequ_by_con_ood(:,j,i)=generate_sequence(ood_ahat(constellation_indices(j,1):constellation_indices(j,2),constellation_indices(j,1):constellation_indices(j,2)),1500);
    end
end
% collect 100x88 p-values and t-values
t_vals_ood = zeros(88,100);
p_vals_ood = zeros(88,100);
for i = 1:100
    freqs_ood = {1,88};
    for j = 1:88
        temp = zeros(constellation_indices(j,2)-constellation_indices(j,1)+1);
        for k = 1:length(temp)
            for l = 1:1500
                if sequ_by_con_ood(l,j,i) == k
                    temp(k) = temp(k) +1;
                end
            end
        end
        freqs_ood{j} = temp;
    end
    for m = 1:88
        X = freqs_ood{m}(find(bycon_1or0{m}>0));
        Y = freqs_ood{m}(find(bycon_1or0{m}<1));
        [h p cstat tstat] = ttest2(X, Y);
        t_vals_ood(m,i) = tstat.tstat;
        p_vals_ood(m,i) = p;
    end
end
% get mean p-value of each constellation
mean_p_vals_ood = zeros(88,1);
for i = 1:88
    mean_p_vals_ood(i) = mean(p_vals_ood(i,:));
end
% average to get 1 t-value per constellation
av_tvals_ood = zeros(88,1);
for i = 1:88
    av_tvals_ood(i) = mean(t_vals_ood(i,:));
end
%% level 2 real t- and p-vals
% 100 different sequences from ahat real saccade model
sequ_by_con_sacc = zeros(1500,88,100);
for i = 1:100
    for j = 1:88
        sequ_by_con_sacc(:,j,i)=generate_sequence(sacc_ahat(constellation_indices(j,1):constellation_indices(j,2),constellation_indices(j,1):constellation_indices(j,2)),1500);
    end
end
% collect 100x88 p-values and t-values
t_vals_sacc = zeros(88,100);
p_vals_sacc = zeros(88,100);
for i = 1:100
    freqs_sacc = {1,88};
    for j = 1:88
        temp = zeros(constellation_indices(j,2)-constellation_indices(j,1)+1);
        for k = 1:length(temp)
            for l = 1:1500
                if sequ_by_con_sacc(l,j,i) == k
                    temp(k) = temp(k) +1;
                end
            end
        end
        freqs_sacc{j} = temp;
    end
    for m = 1:88
        X = freqs_sacc{m}(find(bycon_1or0{m}>0));
        Y = freqs_sacc{m}(find(bycon_1or0{m}<1));
        [h p cstat tstat] = ttest2(X, Y);
        t_vals_sacc(m,i) = tstat.tstat;
        p_vals_sacc(m,i) = p;
    end
end
mean_p_vals_sacc = zeros(88,1);
for i = 1:88
    mean_p_vals_sacc(i) = mean(p_vals_sacc(i,:));
end
av_tvals_sacc = zeros(88,1);
for i = 1:88
    av_tvals_sacc(i) = mean(t_vals_sacc(i,:));
end
%% level 3 real t- and p-vals
% 100 different sequences from ahat real magnitude model
sequ_by_con_mags = zeros(1500,88,100);
for i = 1:100
    for j = 1:88
        sequ_by_con_mags(:,j,i)=generate_sequence(mags_ahat(constellation_indices(j,1):constellation_indices(j,2),constellation_indices(j,1):constellation_indices(j,2)),1500);
    end
end
% collect 100x88 p-values and t-values
t_vals_mags = zeros(88,100);
p_vals_mags = zeros(88,100);
for i = 1:100
    freqs_mags = {1,88};
    for j = 1:88
        temp = zeros(constellation_indices(j,2)-constellation_indices(j,1)+1);
        for k = 1:length(temp)
            for l = 1:1500
                if sequ_by_con_mags(l,j,i) == k
                    temp(k) = temp(k) +1;
                end
            end
        end
        freqs_mags{j} = temp;
    end
    for m = 1:88
        X = freqs_mags{m}(find(bycon_1or0{m}>0));
        Y = freqs_mags{m}(find(bycon_1or0{m}<1));
        [h p cstat tstat] = ttest2(X, Y);
        t_vals_mags(m,i) = tstat.tstat;
        p_vals_mags(m,i) = p;
    end
end
mean_p_vals_mags = zeros(88,1);
for i = 1:88
    mean_p_vals_mags(i) = mean(p_vals_mags(i,:));
end
av_tvals_mags = zeros(88,1);
for i = 1:88
    av_tvals_mags(i) = mean(t_vals_mags(i,:));
end

%% make null transition probability matricies
%% make null level 1 -distance model- and find 100 zscores
% permuted RA and real DEC mode
ood_perm_partition = zeros(9081,100);
ood_perm_a = zeros(100,1);
ood_perm_zscore = zeros(100,1);
for i = 1:100
    RA_perm = RA_all(randperm(N_all));
    distanceperm = zeros(N_all,N_all);
    for j = 1:N_all
        for k = (j+1):(N_all)
            distanceperm(j,k) = acosd(((sind(DEC_all(j,1)))*(sind(DEC_all(k,1))))+((cosd(DEC_all(j,1)))*(cosd(DEC_all(k,1)))*(cosd((RA_perm(j,1))-(RA_perm(k,1))))));
        end
    end
    distanceperm = distanceperm' + distanceperm;
    for n = 1:N_all
        for o = 1:N_all
            if distanceperm(n,o) < 1.0
                distanceperm(n,o) = 1; % cap min distance at 1 degree
            end
        end
    end
    ood_perm = 1./distanceperm;
    for m = 1:N_all
        ood_perm(m,m) = 0; % make diagonal=0
    end
    for l = 1:9081
        ood_perm(l,:) = ood_perm(l,:)./sum(ood_perm(l,:)); % normalize
    end
    ood_rand_ahat = human_expectations(ood_perm, 0.3);
    ood_rand_ahat = (ood_rand_ahat + ood_rand_ahat')./2;
    ood_perm_partition(:,i) = spectralcluster(ood_rand_ahat, 88, 'Distance','precomputed', 'ClusterMethod', 'kmedoids');
    ood_perm_a(i) = find_a(ood_perm_partition(:,i), real_partition);
    ood_perm_zscore(i) = (ood_perm_a(i) - mean_v)./std_v;
end
% random RA and random DEC model
ood_rand_partition = zeros(9081,100);
ood_rand_a = zeros(100,1);
ood_rand_zscore = zeros(100,1);
for i = 1:100
    RA_rand = 360*rand(9081,1);
    sin_DEC_rand = 2*rand(9081,1)-1;
    DEC_rand = asind(sin_DEC_rand);
    distance_rand = zeros(N_all,N_all);
    for j = 1:N_all
        for k = (j+1):(N_all)
        distance_rand(j,k) = acosd(((sind(DEC_rand(j,1)))*(sind(DEC_rand(k,1))))+((cosd(DEC_rand(j,1)))*(cosd(DEC_rand(k,1)))*(cosd((RA_rand(j,1))-(RA_rand(k,1)))))); 
        end
    end
    distancerand = distancerand' + distancerand;
    for n = 1:N_all
        for o = 1:N_all
            if distancerand(n,o) < 1.0
                distancerand(n,o) = 1; % cap min distance at 1 degree
            end
        end
    end
    ood_rand = 1./distancerand;
    for m = 1:N_all
        ood_rand(m,m) = 0; % make diagonal=0
    end
    for l = 1:9081
        ood_rand(l,:) = ood_rand(l,:)./sum(ood_rand(l,:)); % normalize
    end
    ood_rand_ahat = human_expectations(ood_rand, 0.3);
    ood_rand_ahat = (ood_rand_ahat + ood_rand_ahat')./2;
    ood_rand_partition(:,i) = spectralcluster(ood_rand_ahat, 88, 'Distance','precomputed', 'ClusterMethod', 'kmedoids');
    ood_rand_a(i) = find_a(ood_rand_partition(:,i), real_partition);
    ood_rand_zscore(i) = (ood_rand_a(i) - mean_v)./std_v;
end
%% make null level 2 -saccade model- and find 100 zscores
unisacc= zeros(N_all,N_all);
for i=(1:N_all)
    for j = 1:N_all
        if (distance_real(i,j) <= 180 && distance_real(i,j) > (10.*(180./11)))
            unisacc(i,j) = finalscalednew(1);
        elseif (distance_real(i,j) <= (10.*(180./11)) && distance_real(i,j) > (9.*(180./11)))
            unisacc(i,j) = finalscalednew(2);
        elseif (distance_real(i,j) <= (9.*(180./11)) && distance_real(i,j) > (8.*(180./11)))
            unisacc(i,j) = finalscalednew(3);
        elseif (distance_real(i,j) <= (8.*(180./11)) && distance_real(i,j) > (7.*(180./11)))
            unisacc(i,j) = finalscalednew(4);
        elseif (distance_real(i,j) <= (7.*(180./11)) && distance_real(i,j) > (6.*(180./11)))
            unisacc(i,j) = finalscalednew(5);
        elseif (distance_real(i,j) <= (6.*(180./11)) && distance_real(i,j) > (5.*(180./11)))
            unisacc(i,j) = finalscalednew(6);
        elseif (distance_real(i,j) <= (5.*(180./11)) && distance_real(i,j) > (4.*(180./11)))
            unisacc(i,j) = finalscalednew(7);
        elseif (distance_real(i,j) <= (4.*(180./11)) && distance_real(i,j) > (3.*(180./11)))
            unisacc(i,j) = finalscalednew(8);
        elseif (distance_real(i,j) <= (3.*(180./11)) && distance_real(i,j) > (2.*(180./11)))
            unisacc(i,j) = finalscalednew(9);
        elseif (distance_real(i,j) <= (2.*(180./11)) && distance_real(i,j) > (1.*(180./11)))
            unisacc(i,j) = finalscalednew(10);
        elseif (distance_real(i,j) <= (1.*(180./11))) 
            unisacc(i,j) = finalscalednew(11);
        else
            unisacc(i,j)=0;
        end
    end
end
for i = 1:N_all
    unisacc(i,i) = 0;
end
unisacc_prop = zeros(9081,9081);
for i=1:N_all
    saccade_sums_byrow = zeros(11,1);
    saccade_sums_byrow(1) = (sum(unisacc(i,:)==finalscalednew(1)));
    saccade_sums_byrow(2) = (sum(unisacc(i,:)==finalscalednew(2)));
    saccade_sums_byrow(3) = (sum(unisacc(i,:)==finalscalednew(3)));
    saccade_sums_byrow(4) = (sum(unisacc(i,:)==finalscalednew(4)));
    saccade_sums_byrow(5) = (sum(unisacc(i,:)==finalscalednew(5)));
    saccade_sums_byrow(6) = (sum(unisacc(i,:)==finalscalednew(6)));
    saccade_sums_byrow(7) = (sum(unisacc(i,:)==finalscalednew(7)));
    saccade_sums_byrow(8) = (sum(unisacc(i,:)==finalscalednew(8)));
    saccade_sums_byrow(9) = (sum(unisacc(i,:)==finalscalednew(9)));
    saccade_sums_byrow(10) = (sum(unisacc(i,:)==finalscalednew(10)));
    saccade_sums_byrow(11) = (sum(unisacc(i,:)==finalscalednew(11)));
    for j = 1:N
        if unisacc(i,j) == finalscalednew(1)
    unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(1);
        elseif unisacc(i,j) == finalscalednew(2)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(2);
        elseif unisacc(i,j) == finalscalednew(3)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(3);
        elseif unisacc(i,j) == finalscalednew(4)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(4);
        elseif unisacc(i,j) == finalscalednew(5)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(5);
        elseif unisacc(i,j) == finalscalednew(6)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(6);
        elseif unisacc(i,j) == finalscalednew(7)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(7);
        elseif unisacc(i,j) == finalscalednew(8)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(8);
        elseif unisacc(i,j) == finalscalednew(9)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(9);
        elseif unisacc(i,j) == finalscalednew(10)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(10);
        elseif unisacc(i,j) == finalscalednew(11)
            unisacc_prop(i,j)=unisacc(i,j)./saccade_sums_byrow(11);
        else
            unisacc_prop(i,j)=0;
        end
    end
end
unisacc_partition = zeros(9081,100);
unisacc_a = zeros(100,1);
unisacc_zscore = zeros(100,1);
for h=1:100
    unisacc = unisacc_prop;
    for i = 1:9081 % add noise
        for j = 1:9081
            a = randn*2.5;
            a = a./100;
            unisacc(i,j) = unisacc(i,j)+(a*unisacc(i,j));
        end
    end
    for k = 1:9081
        unisacc(k,:) = unisacc(k,:)./sum(unisacc(k,:));
    end
    unisacc_ahat = human_expectations(unisacc,0.3);
    unisacc_ahat = (unisacc_ahat+unisacc_ahat')./2;
    unisacc_partition(:,h) = spectralcluster(unisacc_ahat, 88, 'Distance','precomputed', 'ClusterMethod', 'kmedoids');
    unisacc_a(h) = find_a(unisacc_partition(:,h), real_partition);
    unisacc_zscore(h) = (unisacc_a(h) - mean_v)./std_v;
end
%% make null level 3 -magnitude model- and find 88x100 t-values
% random magnitude model
max_mags = max(shiftedappmags);
min_mags = min(shiftedappmags);
range_mags = range(shiftedappmags);
% make 100 vectors of random magnitude values
rand_mags = zeros(9081,100);
for i = 1:100
    rand_mags(:,i)= (range_mags*rand(1,N))+min_mags;
end
% make 100 sequences, one on each model
rand_mags_sequ = zeros(1500,88,100);
for i = 1:100
    rand_mags_trans = zeros(N_all,N_all);
    for h = 1:N_all
        count = 0;
        local_stars = [];
        bright = [];
        for j = 1:N_all
            if distance_real(h,j) < r
                count = count+1;
                local_stars(count) = rand_mags(j,i);
            end
        end
        A = ceil((x/100).*length(local_stars));
        bright = mink(local_stars, A);
        m_avg = mean(bright);
        rand_mags_trans(:,h) = (10^(-0.4*(rand_mags(h,i) - m_avg)))*saccade_trans_prop(:,h);
    end
    for k = 1:N_all
        A=sum(rand_mags_trans(k,:));
        rand_mags_trans(k,:)=rand_mags_trans(k,:)./A;
    end
    rand_mags_trans_ahat = human_expectations(rand_mags_trans, 0.3);
    for l = 1:88
        rand_mags_sequ(:,l,i)=generate_sequence(rand_mags_trans_ahat(constellation_indicies(l,1):constellation_indicies(l,2),constellation_indicies(l,1):constellation_indicies(l,2)),1500);
    end
end
% find 100x88 t-values and p values
t_vals_rand_mags = zeros(88,100);
p_vals_rand_mags = zeros(88,100);
for i = 1:100
    freqs_randmags = {1,88};
    for j = 1:88
        temp = zeros(constellation_indicies(j,2)-constellation_indicies(j,1)+1);
        for k = 1:length(temp)
            for l = 1:1500
                if rand_mags_sequ(l,j,i) == k
                    temp(k) = temp(k) +1;
                end
            end
        end
        freqs_randmags{j} = temp;
    end
    for m = 1:88
        X = freqs_randmags{m}(find(bycon_1or0{m}>0));
        Y = freqs_randmags{m}(find(bycon_1or0{m}<1));
        [h p cstat tstat] = ttest2(X, Y);
        t_vals_rand_mags(m,i) = tstat.tstat;
        p_vals_rand_mags(m,i) = p;
    end
end
% average p-values to get 1 per constellation
mean_p_vals_rand_mags = zeros(88,1);
for i = 1:88
    mean_p_vals_rand_mags(i) = mean(p_vals_rand_mags(i,:));
end
% average t-values to get 1 per constellation
av_tvals_rand_mags = zeros(100,1);
for i = 1:100
    av_tvals_rand_mags(i) = mean(t_vals_rand_mags(:,i));
end

% repulsive magnitude model
% make vector of repulsive magnitudes
repulsive_mags = zeros(9081,1);
for i=1:9081
    repulsive_mags(i) = max_mags - shiftedappmags(i) + min_mags;
end
% add repmags to real saccade matrix
rep_mags_trans = zeros(N_all,N_all);
for j=1:N_all
    count = 0;
    local_stars = [];
    bright = [];
    for k = 1:N_all
        if distance_real(j,k) < r
            count = count+1;
            local_stars(count) = repulsive_mags(k);
        end
    end
    A = ceil((x/100).*length(local_stars));
    bright = mink(local_stars, A);
    m_avg = mean(bright);
    rep_mags_trans(:,j) = (10^(-0.4*(repulsive_mags(j) - m_avg)))*saccade_trans_prop(:,j);
end
% add noise and run 100x88 sequences
rep_mags_sequ = zeros(1500,88,100);
for i = 1:100
    rep_mags_trans_noise = rep_mags_trans;
    for j = 1:9081
       for k = 1:9081
           a = randn*2.5;
           a = a/100;
           rep_mags_trans_noise(i,j) = rep_mags_trans_noise(i,j)+rep_mags_trans_noise(i,j)*a;
       end
    end
    for m=1:N
        rep_mags_trans_noise(m,:) = rep_mags_trans_noise(m,:)./sum(rep_mags_trans_noise(m,:));
    end
    rep_mags_trans_noise_ahat = human_expectations(rep_mags_trans_noise, 0.3);
    for l = 1:88
        rep_mags_sequ(:,l,i)=generate_sequence(rep_mags_trans_noise_ahat(constellation_indicies(l,1):constellation_indicies(l,2),constellation_indicies(l,1):constellation_indicies(l,2)),1500);
    end
end
% find 100x88 t-values and p values
t_vals_rep_mags = zeros(88,100);
p_vals_rep_mags = zeros(88,100);
for i = 1:100
    freqs_repmags = {1,88};
    for j = 1:88
        temp = zeros(constellation_indicies(j,2)-constellation_indicies(j,1)+1);
        for k = 1:length(temp)
            for l = 1:1500
                if rep_mags_sequ(l,j,i) == k
                    temp(k) = temp(k) +1;
                end
            end
        end
        freqs_repmags{j} = temp;
    end
    for m = 1:88
        X = freqs_repmags{m}(find(bycon_1or0{m}>0));
        Y = freqs_repmags{m}(find(bycon_1or0{m}<1));
        [h p cstat tstat] = ttest2(X, Y);
        t_vals_rep_mags(m,i) = tstat.tstat;
        p_vals_rep_mags(m,i) = p;
    end
end
% average p-values to get 1 per constellation
mean_p_vals_rep_mags = zeros(88,1);
for i = 1:88
    mean_p_vals_rep_mags(i) = mean(p_vals_rep_mags(i,:));
end
% average t-values to get 1 per constellation
av_tvals_rep_mags = zeros(88,1);
for i = 1:88
    av_tvals_rep_mags(i) = mean(t_vals_rep_mags(i,:));
end

%% whole sky edge frequency analysis 
% make 100 whole-sky sequences on A of magnitude model
sky_sequences_a = zeros(132000,100);
for i = 1:100
    sky_sequences_a(:,i) = generate_sequence(saccade_mags_trans, 132000);
    i
end
% make an edge frequency matrix for A
for i = 1:100
    a = zeros(9081,9081);
    for j = 1:132000-1
        start_star = sky_sequences_a(j,i);
        end_star = sky_sequences_a(j+1,i);
        a(start_star, end_star) = a(start_star, end_star)+1;
    end
    sky_edge_freqs_a{i} = a;
end
% calculate ratio of transitions between 2 line-drawing stars/other 
% transitions for A
ratio_a = zeros(100,1);
for i = 1:100
    sum_in = 0;
    sum_out = 0;
    for j = 1:N_all
        for k = 1:N_all
            if inorout_lines(j) == 1 && inorout_lines(k) == 1
                sum_in = sum_in + sky_edge_freqs_a{i}(j,k);
            else
                sum_out = sum_out + sky_edge_freqs_a{i}(j,k);
            end
        end
    end
    sum_in_avg = sum_in./((N_lines*N_lines)-696);
    sum_out_avg = sum_out./((N_all*N_all) - ((N_lines*N_lines)-696));
    ratio_a(i) = sum_in_avg./sum_out_avg;
end

% make 100 whole-sky sequences on A_hat of mangitude model
saccade_mags_ahat = human_expectations(saccade_mags_trans, 0.3);
for i = 1:9081
    saccade_mags_ahat(i,i) = 0;
end
sky_sequences_ahat = zeros(132000,100);
for i = 1:100
    sky_sequences_ahat(:,i) = generate_sequence(saccade_mags_ahat, 132000);
    i
end
% make an edge frequency matrix for A_hat
for i = 1:100
    a = zeros(9081,9081);
    for j = 1:132000-1
        start_star = sky_sequences_ahat(j,i);
        end_star = sky_sequences_ahat(j+1,i);
        a(start_star, end_star) = a(start_star, end_star)+1;
    end
    sky_edge_freqs_ahat{i} = a;
end
% calculate ratio of transitions between 2 line-drawing stars/other 
% transitions for Ahat 
ratio_ahat = zeros(100,1);
for i = 1:100
    sum_in = 0;
    sum_out = 0;
    for j = 1:N_all
        for k = 1:N_all
            if inorout_lines(j) == 1 && inorout_lines(k) == 1
                sum_in = sum_in + sky_edge_freqs_ahat{i}(j,k);
            else
                sum_out = sum_out + sky_edge_freqs_ahat{i}(j,k);
            end
        end
    end
    sum_in_avg = sum_in./((N_lines*N_lines)-696);
    sum_out_avg = sum_out./((N_all*N_all) - ((N_lines*N_lines))-696);
    ratio_ahat(i) = sum_in_avg./sum_out_avg;
end
%% whole sky sequence null ratios
% for A
null_ratio_a = zeros(100,1);
for i = 1:100
    inorout_lines_null = inorout_lines(randperm(9081));
    sum_in = 0;
    sum_out = 0;
    for j = 1:N_all
        for k = 1:N_all
            if inorout_lines_null(j) == 1 && inorout_lines_null(k) == 1
                sum_in = sum_in + sky_edge_freqs_a{i}(j,k);
            else
                sum_out = sum_out + sky_edge_freqs_a{i}(j,k);
            end
        end
    end
    sum_in_avg = sum_in./((N_lines*N_lines)-696);
    sum_out_avg = sum_out./((N_all*N_all) - ((N_lines*N_lines)-696));
    null_ratio_a(i) = sum_in_avg./sum_out_avg;
end
mean_null_ratio1_a = mean(null_ratio_a);
std_null_ratio1_a = std(null_ratio_a);
% for Ahat
null_ratio_ahat = zeros(100,1);
for i = 1:100
    inorout_lines_null = inorout_lines(randperm(9081));
    sum_in = 0;
    sum_out = 0;
    for j = 1:N_all
        for k = 1:N_all
            if inorout_lines_null(j) == 1 && inorout_lines_null(k) == 1
                sum_in = sum_in + sky_edge_freqs_ahat{i}(j,k);
            else
                sum_out = sum_out + sky_edge_freqs_ahat{i}(j,k);
            end
        end
    end
    sum_in_avg = sum_in./((N_lines*N_lines)-696);
    sum_out_avg = sum_out./((N_all*N_all) - ((N_lines*N_lines)-696));
    null_ratio_ahat(i) = sum_in_avg./sum_out_avg;
end
mean_null_ratio1_ahat = mean(null_ratio_ahat);
std_null_ratio1_ahat = std(null_ratio_ahat);

%% Make boxplots for figures
% figure 3 - comparing real models with null models
% level 1 real vs null zscores
ood_zscores = [ood_zscore, ood_perm_zscore ood_rand_zscore];
boxplot(ood_zscores, {'real', 'permuted', 'random'}, 'PlotStyle','compact')
title('100 real vs null distance model z-scores')
% level 2 real vs null zscores
sacc_zscores = [sacc_zscore repsac_zscore];
boxplot(sacc_zscores, {'real', 'uniform'}, 'PlotStyle','compact')
title('100 real vs null saccade model zscores')
% level 3 real vs null t-values
mags_tvals = [av_tvals_mags av_tvals_rand_mags av_tvals_rep_mags];
boxplot(mags_tvals, {'real', 'random', 'repulisve'}, 'PlotStyle','compact')
title('100 real vs null magnitude model tvals')

% figure 4 - comparing levels of the real model
% comparing zscores
real_zscores = [ood_zscore sacc_zscore mags_zscore];
boxplot(real_zscores, {'level1', 'level2', 'level3'}, 'PlotStyle','compact')
title('100 real zscores for each model')
% comparing tvalues
real_tvals = [av_tvals_ood av_tvals_sacc av_tvals_mags];
boxplot(real_tvals, {'level1', 'level2', 'level3'}, 'PlotStyle','compact')
title('100 real tvals for each model')

% figure 5 - comparing ratios for whole sky sequences
all_ratios = [ratio_a, ratio_ahat, null_ratio_a, null_ratio_ahat];
boxplot(all_ratios, {'real on A', 'null on A', 'real on Ahat', 'null on ahat'}, 'PlotStyle', 'compact')
title('ratio of stars in line-drawings to stars out of line-drawings')