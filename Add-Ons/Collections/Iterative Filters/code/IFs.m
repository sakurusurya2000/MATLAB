function IMF = IFs(f,options)

%
%  function IMF = IFs(f,options)
%
% It generates the decomposition of the signal f :
%
%  f = IMF(1,:) + IMF(2,:) + ... + IMF(size(IMF, 1), :)
%
% where the last row in the matrix IMF is the trend and the other rows
% are actual IMFs
%
%                                Inputs
%
%   f         Signal to be decomposed
%
%   options    Structure, generated using function decompSettings_IF, containing
%              all the parameters needed in the various algorithms
%
%                               Output
%
%   IMF       Matrices containg in row i the i-th IMF. The last row
%              contains the remainder-trend.
%
%   See also SETTINGS_IFs, GETMASK_V1, MAXMINS, PLOT_IMF_V8.
%
% Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%      Signal Decomposition and Instantaneous Frequency analysis'. 
%      http://arxiv.org/abs/1411.6051
%

%% deal with the input

if nargin < 1,  help IFs; return; end
if nargin < 2, options = Settings_IFs; end

FigCol = 'ckmygr'; % Plot Colors

N = length(f);
if size(f,1)>size(f,2)
    f = f.';
end
if size(f,1)>1
    disp('Wrong dataset, the signal must be a single row vector')
end
IMF =[];

nameFile=sprintf('%1.0d',sum(round(clock*1000)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Iterative Filtering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
load('prefixed_double_filter','MM');

maxmins=Maxmins(f,options.IFs.extensionType);

k = length(maxmins);
diffMaxmins=diff(maxmins);

if options.plots>1
    figMask=figure;
    title(['Mask length IMF_{' num2str(size(IMF,1)+1) '}'])
    plot([f f],'b','LineWidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    hold on
    plot([maxmins N+maxmins],[f(maxmins) f(maxmins)],'kx')
    legend('Signal','Maxes and mins')
    set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %pause(1)
    if figMask > 10
        close all
    end
end

while size(IMF,1) < options.IFs.NIMFs && k>=options.IFs.ExtPoints && toc<options.maxTime
    if options.verbose>0
        fprintf('\n IMF # %1.0d   -   # Extreme points %5.0d\n',size(IMF,1)+1,k)
    end
    SD=1;
    SDlog=[];
    havelog=[];
    h=f;
    
    if isa(options.IFs.alpha,'char')
        if strcmp(options.IFs.alpha,'ave') % Using an average mask length
            m = 2*round(N/k*options.IFs.Xi);
        else
            disp(' Value of alpha not recognized')
            return
        end
    else % using a fixed value alpha
        m = 2*round(options.IFs.Xi*(max(diffMaxmins)*options.IFs.alpha+min(diffMaxmins)*(1-options.IFs.alpha)));
    end
    inStepN=0;
    if options.verbose>0
        fprintf('\n  step #            SD             Mask length \n\n')
    end
    
    a = get_mask_v1(MM,m);
    if options.plots>0 %&& rem(inStepN,5)==0
        figN=figure;
        set(figN,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    while SD>options.IFs.delta && inStepN < options.IFs.MaxInner
        inStepN=inStepN+1;
        
        %             if m > N
        %                 ff =[h];
        %                 for i=1:floor(m/N)
        %                     ff=[h ff h];
        %                 end
        %                 len=rem(m,N);
        %                 ff=[h(N-len+1:N) ff h(1:len)];
        %             else
        %                 ff=[h(N-m+1:N) h h(1:m)];
        %             end
        if strcmp(options.IFs.extensionType,'p')
            if m > N
                ff = [h(N-rem(m,N)+1:N), h, h, h, h(1:rem(m,N))];
            else
                ff = [h(N-m+1:N), h, h(1:m)]; %this could also just be ff = [h, h, h]
            end
        elseif strcmp(options.IFs.extensionType,'r')
            if m > N
                ff = [h(N-rem(m,N)+1:N), fliplr(h), h, fliplr(h), h(1:rem(m,N))];
            else
                ff = [h(m:-1:1), h, h(N:-1:N-m+1)];
            end
        elseif strcmp(options.IFs.extensionType,'c')
            ff = [h(1)*ones(1,m), h, h(end)*ones(1,m)];
        end
        
        
        
        h_ave=zeros(1,N);
        for i=1:N
            h_ave(i)=a*transpose(ff(i:i+2*m));
        end
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        SD=norm(h_ave)^2/norm(h)^2;
        SDlog=[SDlog SD];
        havelog=[havelog norm(h_ave)];
        if options.verbose>0
            fprintf('    %2.0d      %1.14f          %2.0d\n',inStepN,SD,m)
        end
        
        %%%%%%%%%%%%%%%%%% generating f_n %%%%%%%%%%%%%%%%%%
        
        h=h-h_ave;
        
        if options.plots>0 %&& rem(inStepN,5)==0
            if size(IMF,1)>0
                plot_imf_v8([IMF(end,:);h;f-h],1:length(h),sprintf('IMF%1.0d  Step # %5.0d',size(IMF,1)+1,inStepN),[],figN);
            else
                plot_imf_v8([h;f-h],1:length(h),sprintf('IMF%1.0d  Step # %5.0d',size(IMF,1)+1,inStepN),[],figN);
            end
            if options.saveplots>0
                saveas(figN,[nameFile '_IMF' num2str(size(IMF,1)+1) '_fig_' num2str(inStepN)], 'fig')
                saveas(figN,[nameFile '_IMF' num2str(size(IMF,1)+1) '_fig_' num2str(inStepN)], 'epsc')
                saveas(figN,[nameFile '_IMF' num2str(size(IMF,1)+1) '_fig_' num2str(inStepN)], 'png')
            end
            %pause(0.01)
        end
    end    
    if inStepN >= options.IFs.MaxInner
        disp('Max # of inner steps reached')
        return
    end
    if options.saveInter==1
        save([nameFile '_intermediate_results.mat'])
    end
    if options.plots>0 %&& rem(inStepN,5)==0
        close(figN)
    end
    IMF = [IMF; h];
    f=f-h;
    
    maxmins = Maxmins(f,options.IFs.extensionType);
    k = length(maxmins);
    diffMaxmins=diff(maxmins);
    
    if options.plots>1
        figMask=figure;
        title(['Mask length IMF_{' num2str(size(IMF,1)+1) '}'])
        set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        plot([f f],'b','LineWidth',2)
        hold on
        plot([maxmins N+maxmins],[f(maxmins) f(maxmins)],'kx')
        legend('Signal','Maxes and mins')
        set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        %pause(1)
        if figMask > 5
            close all
        end
    end
    
    if options.plots>0
        figN=plot_imf_v8([IMF(end,:);f]);
        for i=1:length(figN)
            set(figN(i),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        end
    end
    
    
end %end of while

IMF = [IMF; f];

if options.plots>=1
    figN=plot_imf_v8(IMF,1:size(IMF,2));
    %pause(0.5)
    for i=1:length(figN)
        set(figN(i),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        if options.saveplots>0
            saveas(figN(i),[nameFile '_IMFs'], 'fig')
            saveas(figN(i),[nameFile '_IMFs'], 'epsc')
            saveas(figN(i),[nameFile '_IMFs'], 'png')
        end
    end
end

if options.saveEnd == 1
    save([ nameFile '_IFs.mat'])
end

end


%% Auxiliar functions


function a=get_mask_v1(y,k)
% 
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
            a=[];
            return
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>10^-14
                fprintf('\n\n Warning!\n\n')
                fprintf(' Area under the mask equals %2.20f\n',norm(a,1))
                fprintf(' it should be equal to 1\n We rescale it using its norm 1\n\n')
        a=a/norm(a,1);
    end
end

end




function maxmins = Maxmins(f,extensionType)


% Identify the maxima and minima of a signal f

if nargin == 1, extensionType = 'c'; end
N = length(f);
maxmins=zeros(1,N);
df = diff(f);


h = 1;
cIn=0;
if strcmp(extensionType,'p') && df(1) == 0 && df(end) == 0
    while df(h)==0
        cIn=cIn+1;
        h=h+1;
    end
end

c = 0;
cmaxmins=0;
for i=h:N-2
    if   df(i)*df(i+1) <= 0
        if df(i+1) == 0
            if c == 0
                posc = i;
            end
            c = c + 1;
        else
            if c > 0
                cmaxmins=cmaxmins+1;
                maxmins(cmaxmins)=posc+floor((c-1)/2)+1;
                c = 0;
            else
                cmaxmins=cmaxmins+1;
                maxmins(cmaxmins)=i+1;
            end
        end
    end
end
if c > 0
    cmaxmins=cmaxmins+1;
    maxmins(cmaxmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    if maxmins(cmaxmins)==0
        maxmins(cmaxmins)=N;
    end
end

maxmins=maxmins(1:cmaxmins);

if strcmp(extensionType,'p') % we deal with a periodic signal
    if isempty(maxmins)
        maxmins = 1;
    else
        if maxmins(1)~=1 && maxmins(end)~=N
            if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
                maxmins=[1 maxmins];
            end
        end
    end
elseif strcmp(extensionType,'c')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
elseif strcmp(extensionType,'r')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
end

end