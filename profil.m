%{
This code takes xyz input files for visualizing profiles taken from
bathymetry. Also, this code calculates Dean profile for average and each
individual profiles.

C. Özsoy (2020)

%}

clc
clear
close all
tic

%% Inputs

proNameFormat='profile%d.xyz';                                              % Enter profile name format (for reading)
proNum=3;                                                                   % Enter profile number
closDepth=0;                                                                % Enter depth for maximum calculation limit of Dean profile (for not limiting give 0)
manualLim=1;                                                                % (1):Axis limits are manual, (0):Axis limits are automatic
lDistLim=-100;                                                              % Enter lower distance limit
uDistLim=2200;                                                              % Enter upper distance limit
lineWidth=2;                                                                % Enter linewidth
fontName='Calibri';                                                         % Enter fontname
fontSize=25;                                                                % Enter fontsize
language=2;                                                                 % Enter language

%% Profile Drawer

for i=1:proNum
    proName=sprintf(proNameFormat,i);
    temp=dlmread(proName);
    x(:,i)=temp(:,1);
    y(:,i)=temp(:,2);
    z(:,i)=temp(:,3);
    clear temp;
    [index,~]=find(min(abs(z(:,i)))==abs(z(:,i)));
    xRel(:,i)=x(index,i)-x(:,i);
    yRel(:,i)=y(index,i)-y(:,i);
    zRel(:,i)=z(:,i)-z(index,i);
    dist(:,i)=sqrt(xRel(:,i).^2+yRel(:,i).^2);
    [index2,~]=find(zRel(:,i)>0);
    dist(index2,i)=-dist(index2,i);
    figure('windowState','maximized');
    hold all
    plot(dist(:,i),zRel(:,i),'LineWidth',lineWidth);
    if language==1
        xlabel('Distance to the Shore (m)');
        ylabel('Depth (m)');
        legendText{i}=sprintf('Profile %d',i);
    elseif language==2
        xlabel('Kýyýya Dik Uzaklýk (m)');
        ylabel('Derinlik (m)');
        legendText{i}=sprintf('Profil %d',i);
    end
    if manualLim==1
        xlim([lDistLim uDistLim]);
    end
    legend(legendText{i});
    ax=gca;
    grid on;
    ax.FontName=fontName;
    ax.FontSize=fontSize;
    hold off;
end

%% Average Profile Calculator/Drawer

distTarget=(min(min(dist)):0.1:max(max(dist)))';                            % Defining a distance range
%distTarget=(-200:0.1:10000)';                                               % Defining a distance range
zTarget=zeros(numel(distTarget),proNum);                            

for i=1:proNum
    zTarget(:,i)=interp1(dist(:,i),zRel(:,i),distTarget);                   % Interpolating depth values at designated distance range
end

zAvg=sum(zTarget,2)/proNum;
indexRemove=find(isnan(zAvg));
distTarget(indexRemove,:)=[];                                               % Removing NAN values
zAvg=rmmissing(zAvg);                                                       % Removing NAN values
[indexAvg,~]=find(min(abs(zAvg(:,1)))==abs(zAvg(:,1)));
distAvgRel(:,1)=distTarget(indexAvg,1)-distTarget(:,1);
zAvgRel(:,1)=zAvg(:,1)-zAvg(indexAvg,1);
distAvgRel=-distAvgRel;
figure('windowstate','maximized');
hold all;

for i=1:proNum
    plot(dist(:,i),zRel(:,i),'LineWidth',lineWidth);
end

plot(distAvgRel,zAvgRel,'--','Color','k','LineWidth',lineWidth+1);
if language==1
    xlabel('Distance to the Shore (m)');
    ylabel('Depth (m)');
    legendText{proNum+1}=sprintf('Average Profile');
elseif language==2
    xlabel('Kýyýya Dik Uzaklýk (m)');
    ylabel('Derinlik (m)');
    legendText{proNum+1}=sprintf('Ortalama Profil');
end
if manualLim==1
    xlim([lDistLim uDistLim]);
end
legend(legendText);
ax=gca;
grid on;
ax.FontName=fontName;
ax.FontSize=fontSize;
hold off;

%% Dean Profile Calculator/Drawer

if closDepth==0
    closDepth=9999;
end

ft=fittype('a*x^(2/3)','independent','x','dependent','y');                  % Defining fit parameters
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Display='Off';
opts.StartPoint=0.001;
A=zeros(1,proNum+1);

for i=1:proNum                                                              % Dean profile for profiles
    tempZ=zRel(:,i);
    %index3=find(tempZ<=0 & abs(tempZ)<=closDepth);
    index3=find(tempZ<=0 & abs(tempZ)>closDepth);
    tempX=dist(:,i);
    if isempty(index3)
        coefZ=tempZ;
        coefX=tempX;
    else
        coefZ=tempZ(1:index3(1));
        coefX=tempX(1:index3(1));
    end
    coefZ=coefZ(coefZ<=0);
    %coefZ=tempZ(index3);
    coefZ=-coefZ;
    coefX=coefX(coefX>=0);
    [fitRes,~]=fit(coefX,coefZ,ft,opts);
    A(i)=coeffvalues(fitRes);
end

%index4=find(zAvgRel<=0 & abs(zAvgRel)<=closDepth);                          % Dean profile for average profile
index4=find(zAvgRel<=0 & abs(zAvgRel)>closDepth);
if isempty(index4)
    zDean=zAvgRel;
    deanDist=distAvgRel;
else
    zDean=zAvgRel(1:index4(1));
    deanDist=distAvgRel(1:index4(1));
end
zDean=zDean(zDean<=0);
%zDean=zAvgRel(index4);
zDean=-zDean;
deanDist=deanDist(deanDist>=0);
[fitResults,gof]=fit(deanDist,zDean,ft,opts);
A(proNum+1)=coeffvalues(fitResults);

figure('windowstate','maximized');                                          % Plotting figure
hold all;
plot(deanDist,-zDean,'--','Color','k','LineWidth',lineWidth);
if language==1
    xlabel('Distance to the Shore (m)');
    ylabel('Depth (m)');
    legendTextDean={'Average Profile','Dean Profile'};
elseif language==2
    xlabel('Kýyýya Dik Uzaklýk (m)');
    ylabel('Derinlik (m)');
    legendTextDean={'Ortalama Profil','Dean Profili'};
end
plot(deanDist,-deanDist.^(2/3)*A(proNum+1),'b','LineWidth',lineWidth);
if manualLim==1
    xlim([lDistLim uDistLim]);
end
legend(legendTextDean);
ax=gca;
grid on;
ax.FontName=fontName;
ax.FontSize=fontSize;
hold off;

toc;

