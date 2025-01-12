clear c1New c1
get_GPIO=readtable(uigetfile('.csv')) % get the GPIO file for identification of start of the shock protocol. 
 c1=get_GPIO(find(ismember(get_GPIO{:,2},'GPIO-1')>0),:) % identify the start of the experiment
c1New= removevars(c1,{'ChannelName'});
c1New=table2array(c1New);
clear c1
c1=c1New;

get_valley=find(diff(c1(:,2))>35000); % identify the shock periods based on the GPIO value. 
get_peak=get_valley+1; %check for 33 if it is same for every session
gpioTime_shock= [c1(get_peak,1)];

%ManualTime=[730 780 820 880 910 950 990 1050 1100 1160] % This is the shock time in seconds from the start of the experiment
 %ManualTime=[300 350 390 420 480 530 560 600 660 690]
%gpioTime_shock= [c1(get_peak,1)]+1;

  calcium_data99 = readtable(uigetfile('.csv'),'Range', '2:500000'); % read calcium signals
  calcium_data=table2array(calcium_data99);
  ca_vals=num2cell(calcium_data);

gpioTime_shock=[300 350 390 420 480 530 560 600 660 690 730 780 820 880 910 950 990 1050 1100 1160]% This is the shock time in seconds from the start of the experiment

            for k=2:size(calcium_data,2) 
               % calcium_epochs(k).Vals=zeros(100,11)
                
  for TrialNo=1:length(gpioTime_shock) % get trial based response

  ca_Trialtimes = find(calcium_data(:,1)>gpioTime_shock(TrialNo)-5 & calcium_data(:,1)<gpioTime_shock(TrialNo)+7);
    
   calcium_epochs(k).Vals{TrialNo,1}=calcium_data(ca_Trialtimes,1);
   calcium_epochs(k).Vals{TrialNo,2}=calcium_data(ca_Trialtimes,k);
  end
            end


            % deconvolution of calcium signals to obtain spikes. Using 2 sd as it gives
% good response in the ACC as seen in navigation experiments, also can
% change to 3 to see and compare results. Take the best parameter. 
for k=2:size(calcium_data,2)
pres_data=calcium_data(:,k);
std_dev=std(pres_data);
STD=2;
std_dF=std(calcium_data(:,1:end));
std_dF(1)=NaN;
[pks,locs] = findpeaks(calcium_data(:,k),1:size(calcium_data,1), 'MinPeakHeight', std_dF(:,k)*STD,'MinPeakProminence', (std_dF(:,k)));
pres_spk=calcium_data(locs,1);
Spk_data(k).Spikes=pres_spk ; % creating spike timestamps for each identified cell
end

%% Calculating the spike based response in each trial 

for k = 2:length(Spk_data)
sum1=0; % preallocating the variable
sum2=0;

for TrialNo=1:10
        postShockTimes=find(Spk_data(k).Spikes>gpioTime_shock(TrialNo)+2 & Spk_data(k).Spikes<gpioTime_shock(TrialNo)+5); % postshock spike times
        from2toSpike=find(Spk_data(k).Spikes>gpioTime_shock(TrialNo) & Spk_data(k).Spikes<gpioTime_shock(TrialNo)+2); % shockspike times - make a temporary spike file that gets timestamps within the range
        if isempty(from2toSpike)==0
Spk_data(k).Shock_spk{TrialNo}=Spk_data(k).Spikes(from2toSpike);
        else
           Spk_data(k).Shock_spk{TrialNo}=[]; % creating an empty set if there is no spikes
        end
 if isempty(postShockTimes)==0

Spk_data(k).PostShock_spk{TrialNo}=Spk_data(k).Spikes(postShockTimes);
 else
     Spk_data(k).PostShock_spk{TrialNo}=[]; % creating an empty set if there is no spikes
 end

end

for ii=1:length(Spk_data(k).Shock_spk)
    sum1=sum1+length(Spk_data(k).Shock_spk{ii}) % identify the number of shock spikes 

end

for ii=1:length(Spk_data(k).PostShock_spk)
    sum2=sum2+length(Spk_data(k).PostShock_spk{ii}) % identify the number of post shock spikes 
end

Spk_data(k).AllShockSpkNo=sum1; % total number of shock spikes over all trials
Spk_data(k).LengthSpk=length(Spk_data(k).Spikes); %total number of spikes
Spk_data(k).meanSpk=length(Spk_data(k).Spikes)/1200; % mean firing rate
Spk_data(k).meanShockSpk=sum1/20; % mean firing rate during shock
Spk_data(k).ExpectedshockSpk=(length(Spk_data(k).Spikes)/1200)*20; % expected spike number during shock 
Spk_data(k).ObservedShockSpk=sum1; % observed spike number during shock 

% Calculation of the above features, but in the post shock period here. 
Spk_data(k).AllPostShockSpkNo=sum2;
%Spk_data(k).LengthSpk=length(Spk_data(k).Spikes);
%Spk_data(k).meanSpk=length(Spk_data(k).Spikes)/1200;
Spk_data(k).meanPostShockSpk=sum2/20;
Spk_data(k).ExpectedPostshockSpk=(length(Spk_data(k).Spikes)/1200)*50;
Spk_data(k).ObservedPostShockSpk=sum2;

% chi-square test to study if shok or postshock response is by chance. 
clear n1 n2 N1 N2
       n1 = sum1; N1 = 20;
       n2 = length(Spk_data(k).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

       chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)

       %x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       %x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       %[tbl,chi2stat,pval] = crosstab(x1,x2)

Spk_data(k).Chi_pvalue=pval;

if Spk_data(k).meanShockSpk < Spk_data(k).ObservedShockSpk
Spk_data(k).ShockChange{1}='Increased';
else
    Spk_data(k).ShockChange{1}='Decreased';

end

if pval<0.05
     Spk_data(k).ShockStatus{1}='Significant';
else
     Spk_data(k).ShockStatus{1} = ' Insignificant';
end

       clear n1 n2 N1 N2
       n1 = sum2; N1 = 50;
       n2 = length(Spk_data(k).Spikes); N2 =1200;

       p0 = (n1+n2) / (N1+N2)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

       chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)

Spk_data(k).PostChi_pvalue=pval;

if Spk_data(k).meanPostShockSpk < Spk_data(k).ObservedPostShockSpk
    Spk_data(k).PostShockChange{1}='Increased';
else
    Spk_data(k).PostShockChange{1}='Decreased';
end

if pval<0.05
     Spk_data(k).PostShockStatus{1}='Significant';
else
     Spk_data(k).PostShockStatus{1} = ' Insignificant';
end

end
%%
for kk=2:size(calcium_data,2)

pres_data=calcium_data(:,kk);
std_dev=std(pres_data);
STD=2;
std_dF=std(calcium_data(:,1:end));
std_dF(1)=NaN;
[pks,locs] = findpeaks(calcium_data(:,kk),1:size(calcium_data,1), 'MinPeakHeight', std_dF(:,kk)*STD,'MinPeakProminence', (std_dF(:,kk)));
pres_spk=calcium_data(locs,1)

shuffling_sec=20:1:1200; % shuffling the data randomly by time shifting from 20th second to 1200 seconds

for ii=1:1000
    increment=datasample(shuffling_sec,1,'Replace',true);      % picking UNIQUE random shuffling increment from suffling matrix
    all_inc(kk)=increment;
    pres_spk1=pres_spk+increment; % add the time increment leading to different time in spikes. 

    for isp=1:length(pres_spk1)
        if pres_spk1(isp)>1200
            pres_spk1(isp)=pres_spk1(isp)-1200; % create time shifted spikes
        end
    end  
    Spk_data1(ii).Spikes=sort(pres_spk1) ; % sort time shifted spikes
    sum1=0;
end
Shuffle_Spike(kk).Data=Spk_data1; 
clear Spk_data1
end 

   %% for each shuffled spike session, analyze the shock and postshock response and get the statistics as same as in the original data. 
for k=2:length(Shuffle_Spike)
for ii=1:1000
sum1=0;
sum2=0;
for TrialNo=1:10

        from2toSpike=find(Shuffle_Spike(k).Data(ii).Spikes>gpioTime_shock(TrialNo) & Shuffle_Spike(k).Data(ii).Spikes<gpioTime_shock(TrialNo)+2); % make a temporary spike file that gets timestamps within the range
        fromShocktoPost=find(Shuffle_Spike(k).Data(ii).Spikes>gpioTime_shock(TrialNo)+2 & Shuffle_Spike(k).Data(ii).Spikes<gpioTime_shock(TrialNo)+7);

        if isempty(from2toSpike)==0
Shuffle_Spike(k).Data(ii).Shock_spk{TrialNo}=Shuffle_Spike(k).Data(ii).Spikes(from2toSpike);
        else 
          Shuffle_Spike(k).Data(ii).Shock_spk{TrialNo}=[];  

        end
%post shock definition
%         fromShocktoPost=find(Shuffle_Spike(k).Data(ii).Spikes>gpioTime_shock(TrialNo)+2 & Shuffle_Spike(k).Data(ii).Spikes<gpioTime_shock(TrialNo)+7);
% 
if isempty(fromShocktoPost)==0
Shuffle_Spike(k).Data(ii).PostShock_spk{TrialNo}=Shuffle_Spike(k).Data(ii).Spikes( fromShocktoPost);
else 
    Shuffle_Spike(k).Data(ii).PostShock_spk{TrialNo}=[];   
end
%


end

for isi=1:length(Shuffle_Spike(k).Data(ii).Shock_spk)
    sum1=sum1+length(Shuffle_Spike(k).Data(ii).Shock_spk{isi});

end

for isi=1:length(Shuffle_Spike(k).Data(ii).PostShock_spk)
    sum2=sum2+length(Shuffle_Spike(k).Data(ii).PostShock_spk{isi});

end


Shuffle_Spike(k).Data(ii).AllShockSpkNo=sum1;
Shuffle_Spike(k).Data(ii).LengthSpk=length(Shuffle_Spike(k).Data(ii).Spikes);
Shuffle_Spike(k).Data(ii).meanSpk=length(Shuffle_Spike(k).Data(ii).Spikes)/1200;
Shuffle_Spike(k).Data(ii).meanShockSpk=sum1/20;
Shuffle_Spike(k).Data(ii).ExpectedshockSpk=(length(Shuffle_Spike(k).Data(ii).Spikes)/1200)*20;
Shuffle_Spike(k).Data(ii).ObservedShockSpk=sum1;

Shuffle_Spike(k).Data(ii).PostAllShockSpkNo=sum2;
%Shuffle_Spike(k).Data(ii).LengthSpk=length(Shuffle_Spike(k).Data(ii).Spikes);
%Shuffle_Spike(k).Data(ii).meanSpk=length(Shuffle_Spike(k).Data(ii).Spikes)/1200;
Shuffle_Spike(k).Data(ii).meanPostShockSpk=sum2/20;
Shuffle_Spike(k).Data(ii).ExpectedPostshockSpk=(length(Shuffle_Spike(k).Data(ii).Spikes)/1200)*50;
Shuffle_Spike(k).Data(ii).ObservedPostShockSpk=sum2;


clear n1 n2 N1 N2
       n1 = sum1; N1 = 20;
       n2 = length(Shuffle_Spike(k).Data(ii).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)

       
       %x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       %x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       %[tbl,chi2stat,pval] = crosstab(x1,x2)
Shuffle_Spike(k).Data(ii).Chi_pvalue=pval;

if Shuffle_Spike(k).Data(ii).meanShockSpk < Shuffle_Spike(k).Data(ii).ObservedShockSpk
Shuffle_Spike(k).Data(ii).ShockChange{1}='Increased';
else
   Shuffle_Spike(k).Data(ii).ShockChange{1}='Decreased';

end

if pval<0.05
    Shuffle_Spike(k).Data(ii).ShockStatus{1}='Significant';
else
     Shuffle_Spike(k).Data(ii).ShockStatus{1} = ' Insignificant';
end


clear n1 n2 N1 N2
       n1 = sum2; N1 = 50;
       n2 = length(Shuffle_Spike(k).Data(ii).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)
Shuffle_Spike(k).Data(ii).PostChi_pvalue=pval;


if Shuffle_Spike(k).Data(ii).meanPostShockSpk < Shuffle_Spike(k).Data(ii).ObservedPostShockSpk
Shuffle_Spike(k).Data(ii).PostShockChange{1}='Increased';
else
   Shuffle_Spike(k).Data(ii).PostShockChange{1}='Decreased';

end

if pval<0.05
    Shuffle_Spike(k).Data(ii).PostShockStatus{1}='Significant';
else
     Shuffle_Spike(k).Data(ii).PostShockStatus{1} = ' Insignificant';
end




end
end

%% find the shock shuffle stats
for k=2:length(Spk_data)
    Spk_data(k).meanShockSpk;          
   for ii=1:1000
       
     Shuffle_meanShockSpk(ii)=Shuffle_Spike(k).Data(ii).meanShockSpk;
   end

Spk_data(k).ShockShuffleParms=Shuffle_meanShockSpk(ii);
figure
histogram(Shuffle_meanShockSpk,0:0.04:0.4)
xline(prctile(Shuffle_meanShockSpk,95),'k')
xline(  Spk_data(k).meanShockSpk,'r')
if Spk_data(k).meanShockSpk > prctile(Shuffle_meanShockSpk,95)

Spk_data(k).ShockShuffleStats='Significant';
title(strcat('Cell#',num2str(k),'Significant'))
else
    Spk_data(k).ShockShuffleStats='Not';
    title(strcat('Cell#',num2str(k),'Not'))
end

end
close all

%% do the same for the second condition
for k = 2:length(Spk_data)
sum1=0;
sum2=0;
for TrialNo=11:20
         postShockTimes=find(Spk_data(k).Spikes>gpioTime_shock(TrialNo)+2 & Spk_data(k).Spikes<gpioTime_shock(TrialNo)+5); 

        from2toSpike=find(Spk_data(k).Spikes>gpioTime_shock(TrialNo) & Spk_data(k).Spikes<gpioTime_shock(TrialNo)+2); % make a temporary spike file that gets timestamps within the range
        if isempty(from2toSpike)==0
Spk_data(k).Shock_spk{TrialNo}=Spk_data(k).Spikes(from2toSpike);

        else 
            Spk_data(k).Shock_spk{TrialNo}=[];
        end
 if isempty(postShockTimes)==0

Spk_data(k).PostShock_spk{TrialNo}=Spk_data(k).Spikes(postShockTimes);

 else
     Spk_data(k).PostShock_spk{TrialNo}=[];
 end
 end

for ii=11:length(Spk_data(k).Shock_spk)
    sum1=sum1+length(Spk_data(k).Shock_spk{ii})
end

for ii=11:length(Spk_data(k).PostShock_spk)
    sum2=sum2+length(Spk_data(k).PostShock_spk{ii})
end

Spk_data(k).AllShockSpkNo(2)=sum1;
Spk_data(k).LengthSpk(2)=length(Spk_data(k).Spikes);
Spk_data(k).meanSpk(2)=length(Spk_data(k).Spikes)/1200;
Spk_data(k).meanShockSpk(2)=sum1/20;
Spk_data(k).ExpectedshockSpk(2)=(length(Spk_data(k).Spikes)/1200)*20;
Spk_data(k).ObservedShockSpk(2)=sum1;

Spk_data(k).AllPostShockSpkNo(2)=sum2;
%Spk_data(k).LengthSpk(2)=length(Spk_data(k).Spikes);
%Spk_data(k).meanSpk(2)=length(Spk_data(k).Spikes)/1200;
Spk_data(k).meanPostShockSpk(2)=sum2/20;
Spk_data(k).ExpectedPostshockSpk(2)=(length(Spk_data(k).Spikes)/1200)*50;
Spk_data(k).ObservedPostShockSpk(2)=sum2;

clear n1 n2 N1 N2
       n1 = sum1; N1 = 20;
       n2 = length(Spk_data(k).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)

       
       %x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       %x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       %[tbl,chi2stat,pval] = crosstab(x1,x2)

Spk_data(k).Chi_pvalue(2)=pval;

if Spk_data(k).meanShockSpk(2) < Spk_data(k).ObservedShockSpk(2)
Spk_data(k).ShockChange{2}='Increased';
else
    Spk_data(k).ShockChange{2}='Decreased';

end

if pval<0.05
     Spk_data(k).ShockStatus{2}='Significant';
else
     Spk_data(k).ShockStatus{2} = ' Insignificant';
end


clear n1 n2 N1 N2
       n1 = sum2; N1 = 50;
       n2 = length(Spk_data(k).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)

       
       %x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       %x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       %[tbl,chi2stat,pval] = crosstab(x1,x2)

Spk_data(k).PostChi_pvalue(2)=pval;

if Spk_data(k).meanPostShockSpk(2) < Spk_data(k).ObservedPostShockSpk(2)
Spk_data(k).PostShockChange{2}='Increased';
else
    Spk_data(k).PostShockChange{2}='Decreased';

end

if pval<0.05
     Spk_data(k).PostShockStatus{2}='Significant';
else
     Spk_data(k).PostShockStatus{2} = ' Insignificant';
end
end

%% Creating the shuffle spikes for the second condition
for k=2:length(Shuffle_Spike)

for ii=1:1000
sum1=0;
for TrialNo=11:20

        from2toSpike=find(Shuffle_Spike(k).Data(ii).Spikes>gpioTime_shock(TrialNo) & Shuffle_Spike(k).Data(ii).Spikes<gpioTime_shock(TrialNo)+2); % make a temporary spike file that gets timestamps within the range
          fromShocktoPost=find(Shuffle_Spike(k).Data(ii).Spikes>gpioTime_shock(TrialNo)+2 & Shuffle_Spike(k).Data(ii).Spikes<gpioTime_shock(TrialNo)+7);
        
        
        if isempty(from2toSpike)==0
Shuffle_Spike(k).Data(ii).Shock_spk2{TrialNo}=Shuffle_Spike(k).Data(ii).Spikes(from2toSpike);
        else 
          Shuffle_Spike(k).Data(ii).Shock_spk2{TrialNo}=[];  

        end

 if isempty(fromShocktoPost)==0
Shuffle_Spike(k).Data(ii).PostShock_spk2{TrialNo}=Shuffle_Spike(k).Data(ii).Spikes( fromShocktoPost);
        else 
          Shuffle_Spike(k).Data(ii).PostShock_spk2{TrialNo}=[];  

        end
end

for isi=11:length(Shuffle_Spike(k).Data(ii).Shock_spk2)
    sum1=sum1+length(Shuffle_Spike(k).Data(ii).Shock_spk2{isi});
end

for isi=11:length(Shuffle_Spike(k).Data(ii).PostShock_spk2)
    sum2=sum2+length(Shuffle_Spike(k).Data(ii).PostShock_spk2{isi});
end

Shuffle_Spike(k).Data(ii).AllShockSpkNo2=sum1;
Shuffle_Spike(k).Data(ii).LengthSpk2=length(Shuffle_Spike(k).Data(ii).Spikes);
Shuffle_Spike(k).Data(ii).meanSpk2=length(Shuffle_Spike(k).Data(ii).Spikes)/1200;
Shuffle_Spike(k).Data(ii).meanShockSpk2=sum1/20;
Shuffle_Spike(k).Data(ii).ExpectedshockSpk2=(length(Shuffle_Spike(k).Data(ii).Spikes)/1200)*20;
Shuffle_Spike(k).Data(ii).ObservedShockSpk2=sum1;

Shuffle_Spike(k).Data(ii).AllPostShockSpkNo2=sum2;
%Shuffle_Spike(k).Data(ii).LengthSpk2=length(Shuffle_Spike(k).Data(ii).Spikes);
%Shuffle_Spike(k).Data(ii).meanSpk2=length(Shuffle_Spike(k).Data(ii).Spikes)/1200;
Shuffle_Spike(k).Data(ii).meanPostShockSpk2=sum1/20;
Shuffle_Spike(k).Data(ii).ExpectedPostshockSpk2=(length(Shuffle_Spike(k).Data(ii).Spikes)/1200)*50;
Shuffle_Spike(k).Data(ii).ObservedPostShockSpk2=sum2;

clear n1 n2 N1 N2
       n1 = sum1; N1 = 20;
       n2 = length(Shuffle_Spike(k).Data(ii).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)

       %x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       %x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       %[tbl,chi2stat,pval] = crosstab(x1,x2)
Shuffle_Spike(k).Data(ii).Chi_pvalue2=pval;

if Shuffle_Spike(k).Data(ii).meanShockSpk2 < Shuffle_Spike(k).Data(ii).ObservedShockSpk2
Shuffle_Spike(k).Data(ii).ShockChange{2}='Increased';
else
   Shuffle_Spike(k).Data(ii).ShockChange{2}='Decreased';

end

if pval<0.05
    Shuffle_Spike(k).Data(ii).ShockStatus{2}='Significant';
else
     Shuffle_Spike(k).Data(ii).ShockStatus{2} = ' Insignificant';
end
clear n1 n2 N1 N2
       n1 = sum2; N1 = 50;
       n2 = length(Shuffle_Spike(k).Data(ii).Spikes); N2 =1200;

 p0 = (n1+n2) / (N1+N2)
 n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)
       pval = 1 - chi2cdf(chi2stat,1)
       %x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       %x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       %[tbl,chi2stat,pval] = crosstab(x1,x2)
Shuffle_Spike(k).Data(ii).PostChi_pvalue2=pval;
end
end


%%
%% Identify the mean shock spike for the second condition
for k=2:length(Spk_data)
    Spk_data(k).meanShockSpk;          
   for ii=1:1000
       
     Shuffle_meanShockSpk(ii)=Shuffle_Spike(k).Data(ii).meanShockSpk2;
   end

Spk_data(k).ShockShuffleParms=Shuffle_meanShockSpk(ii);
end

