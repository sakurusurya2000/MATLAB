% Time varying filtering based EMD
function imf = tvf_emd(x)
%input
% x: input signal x(t)
%output
% imf: resulting intrinsic mode functions

 %%%%%%%%%%% preprocessing
 
 THRESH_BWR=0.1; % inst bandwidth  threshold (stopping criterion)
 MAX_IMF=50; % max output imfs
 BSP_ORDER=26; % b-spline order
 end_flag=0;
 imf=zeros(MAX_IMF,numel(x));
if size(x,1) > 1 % convert to row vector
  x = x.';
end
temp_x=x;

t=1:numel(x); 
for nimf=1:MAX_IMF
    
    [indmin_x, indmax_x] = extr(temp_x); % is there enough extrema to continue
    if nimf==MAX_IMF
        imf(nimf,:)=temp_x;
        nimf=nimf+1;
        break;
    end
    
    if(numel([indmin_x, indmax_x])<4)
        imf(nimf,:)=temp_x;
        if numel(find(temp_x~=0))>0
            nimf=nimf+1;
        end
        end_flag=1;
    end
    if end_flag==1
        break;
    end
    
    num_padding = round(length(temp_x)*0.5);%  padding number 
    
    y = temp_x;

    flag_stopiter=0;
    for iter=1:100
        y = [fliplr(y(2:2+num_padding-1)) y fliplr(y(end-num_padding:end-1))]; % padding to deal with boundary effect (symmetric)
        %y=extendsignal(y,num_padding); % padding to deal with boundary effect (match slope)
               
         tt=1:numel(y);
         ind_remov_pad=num_padding+1:numel(y)-num_padding;
         

        %%%%%%%%%%%%% extract inst amp and inst freq using hilbert transform    
        [indmin_y, indmax_y] = extr(y);
        indexC_y=sort([indmin_y, indmax_y]);
        [instAmp0,instFreq0] = INST_FREQ_local(y);

        %%%%%%%%%%%%%% divide y(t) into two parts and compute bis freq (cutoff freq)
        [a1 f1 a2 f2 bis_freq instBWR avgFreq]=devide_y(y,instAmp0,instFreq0);
      
       instBWR2=instBWR;
       for j=1:2:numel(indexC_y)-2
           ind=indexC_y(j):indexC_y(j+2);
           instBWR2(ind)=mean(instBWR(ind));      
       end
        
        bis_freq(instBWR2<THRESH_BWR)=1e-12;
        bis_freq(bis_freq>0.5)=0.45; bis_freq(bis_freq<=0)=1e-12; 
          
        %%%%%%%%%%%% deal with mode mixing
        bis_freq = anti_modemixing(y,bis_freq,ind_remov_pad,num_padding);
        bis_freq=bis_freq(ind_remov_pad);
        bis_freq = [fliplr(bis_freq(2:2+num_padding-1)) bis_freq fliplr(bis_freq(end-num_padding:end-1))];
        
        bis_freq = anti_modemixing(y,bis_freq,ind_remov_pad,num_padding);
        bis_freq=bis_freq(ind_remov_pad);
        bis_freq = [fliplr(bis_freq(2:2+num_padding-1)) bis_freq fliplr(bis_freq(end-num_padding:end-1))];
       
       %%%%%%%%%%%%%% stopping criteria
        temp_instBWR=instBWR2(ind_remov_pad);
        ind_start=round(numel(temp_instBWR)*0.05);% incase boundary effect
        ind_end=round(numel(temp_instBWR)*0.95);
        

        if (iter>=2 && mean(temp_instBWR(ind_start:ind_end))<THRESH_BWR+THRESH_BWR/4*iter) || iter>=6 || (nimf>1 && mean(temp_instBWR(ind_start:ind_end))<THRESH_BWR+THRESH_BWR/4*iter)
            flag_stopiter=1; 
        end
        if numel(find(temp_instBWR(ind_start:ind_end)>THRESH_BWR))/numel(instBWR2(ind_remov_pad))<0.2 % most of the signal are local narrowband
           flag_stopiter=1;
        end

        %%%%%%%%%%%%%% obtain local mean using time varying filtering
        phi=zeros(1,numel(bis_freq)); % build knots via h(t)
        for i=1:numel(bis_freq)-1
            phi(i+1)=phi(i)+2*pi*bis_freq(i);
        end
         [indmin_knot, indmax_knot] = extr(cos(phi)); 
    
        indexC_knot=sort([indmin_knot, indmax_knot]);
        if numel(indexC_knot)>2
            pp_spline = splinefit(1:length(y),y,indexC_knot,BSP_ORDER,'p'); % TVF filtering with the extrema of h(t) as knots
            localmean = ppval(pp_spline,1:length(y));
        else
            flag_stopiter=1;  %not enough knots to perform filtering
        end       
        if (max(abs(y(ind_remov_pad)-localmean(ind_remov_pad)))/min(abs(localmean(ind_remov_pad)))<1e-3) % prevent signal of very small amplitude to be extracted  
            flag_stopiter=1;
        end
    
        temp_residual=y-localmean;
        temp_residual=temp_residual(ind_remov_pad);
        temp_residual=temp_residual(round(numel(temp_residual)*0.1):end-round(numel(temp_residual)*0.1));
        localmean2=localmean(ind_remov_pad);
        localmean2=localmean2(round(numel(localmean2)*0.1):end-round(numel(localmean2)*0.1));
        if abs(max(localmean2))/abs(max(instAmp0(ind_remov_pad)))<3.5e-2 || abs(max(temp_residual))/abs(max(instAmp0(ind_remov_pad)))<1e-2
            flag_stopiter=1;
        end
        
        %%%%%%%%  check stopping criteria
        if flag_stopiter
            imf(nimf,:)=y(ind_remov_pad);
            temp_x=temp_x-y(ind_remov_pad);
            break;
        end
        
        y=y-localmean;      
        y=y(ind_remov_pad);
       
    end
end
imf(nimf:MAX_IMF,:)=[];
end

function output_cutoff = anti_modemixing(y,bis_freq,ind_remov_pad,num_padding)
        org_bis_freq=bis_freq;
        output_cutoff=bis_freq;
        flag_intermitt=0;
        t=1:numel(bis_freq);
        intermitt=[];
        %%%%% locate the maxima of the input signal
        [indmin_y, indmax_y] = extr(y);
        indexC_y=sort([indmin_y, indmax_y]);
        zero_span=[];
        %%%% preprocessing
         for i=2:numel(indmax_y)-1
            time_span=indmax_y(i-1):indmax_y(i+1);
            if (max(bis_freq(time_span))-min(bis_freq(time_span)))/min(bis_freq(time_span))>0.25
                zero_span=[zero_span time_span];
            end
        end
        bis_freq(zero_span)=0;
        %%%% find out all intermittences
        diff_bis_freq=zeros(size(bis_freq));
        for i=1:numel(indmax_y)-1
            time_span=indmax_y(i):indmax_y(i+1);
            if (max(bis_freq(time_span))-min(bis_freq(time_span)))/min(bis_freq(time_span))>0.25
                intermitt=[intermitt indmax_y(i)];
                diff_bis_freq(indmax_y(i))=bis_freq(indmax_y(i+1))-bis_freq(indmax_y(i));
            end
        end
        ind_remov_pad([1:round(0.1*numel(ind_remov_pad)),round(0.9*numel(ind_remov_pad)):end])=[];
        inters=intersect(ind_remov_pad,intermitt);
        if numel(inters) >0
            flag_intermitt=1;
        end
        %%%% find out floors  
       % plot(t(intermitt),bis_freq(intermitt),'r.',t,bis_freq,'b')
        for i=2:numel(intermitt)-1  
            u1=intermitt(i-1);
            u2=intermitt(i);
            u3=intermitt(i+1); % check the derivative of cutoff frequency
            if diff_bis_freq(u2)>0  % rising edge
                bis_freq(u1:u2)=0;
            end
            if diff_bis_freq(u2)<0 % falling edge
                bis_freq(u2:u3)=0;
            end
        end
      
        temp_bis_freq=bis_freq;
        temp_bis_freq(temp_bis_freq<1e-9)=0; %cutoff freq of very small value is considered to be zero (floor value)  
        temp_bis_freq=temp_bis_freq(ind_remov_pad);
        temp_bis_freq = [fliplr(temp_bis_freq(2:2+num_padding-1)) temp_bis_freq fliplr(temp_bis_freq(end-num_padding:end-1))];
        flip_bis_freq=fliplr(bis_freq);
        if numel(find(temp_bis_freq>1e-9))>0 && numel(find(flip_bis_freq>1e-9,1,'first'))>0
             temp_bis_freq(1)=bis_freq(find(bis_freq>1e-9,1,'first'));temp_bis_freq(end)=flip_bis_freq(find(flip_bis_freq>1e-9,1,'first')); %padding boundary for interpolation
        else
            temp_bis_freq(1)=bis_freq(1); temp_bis_freq(end)=bis_freq(end);
        end
       
        bis_freq=temp_bis_freq;
        %%%% interpolate between peaks
        if numel(t(bis_freq~=0))<2
            return;
        end
        bis_freq = interp1(t(bis_freq~=0),bis_freq(bis_freq~=0) ,t,'pchip');
        
         flip_bis_freq=fliplr(org_bis_freq);
        if numel(find(org_bis_freq>1e-9))>0 && numel(find(flip_bis_freq>1e-9,1,'first'))>0
             org_bis_freq(1)=org_bis_freq(find(org_bis_freq>1e-9,1,'first'));org_bis_freq(end)=flip_bis_freq(find(flip_bis_freq>1e-9,1,'first')); %padding boundary for interpolation
        end
        org_bis_freq(org_bis_freq<1e-9)=0;
        org_bis_freq(1)=bis_freq(1);org_bis_freq(end)=bis_freq(end);
        org_bis_freq = interp1(t(org_bis_freq~=0),org_bis_freq(org_bis_freq~=0) ,t,'pchip');
     
        if flag_intermitt && max(temp_bis_freq(ind_remov_pad))>1e-9
            output_cutoff=bis_freq;
        else
            output_cutoff=org_bis_freq;
        end
        
        output_cutoff(output_cutoff>0.45)=0.45;
        output_cutoff(output_cutoff<0)=0;
end



% devide y(t) into two sub-signals a1(t)exp(2*pi*f1(t)) and  a2(t)exp(2*pi*f2(t))
function [a1 f1 a2 f2 bis_freq ratio_bw avgFreq] = devide_y(y,instAmp0,instFreq0)
%input 
%  y: input signal y(t)
%  instAmp0: inst amp of y(t)
%  instFreq0: inst freq of y(t)
%output
%  a1,f1: inst amp and inst amp of the first sub-signal
%  a2,f2: inst amp and inst amp of the second sub-signal
% bis_freq: bisecting frequency = [f1+f2]/2
%  ratio_bw: inst bandwidth ratio (stopping criterion)

    %%%%%%%%%%% preprocessing
    bis_freq=zeros(size(instAmp0));
     tt=1:length(instAmp0);
     squar_instAmp0=instAmp0.^2;
     
     [indmin_y, indmax_y] = extr(y);
     
     
    %%%%%%%%%%% obtain a1(t) and a2(t) 
    [indmin_amp0, indmax_amp0] = extr(squar_instAmp0);

    if numel(indmin_amp0) <2 || numel(indmax_amp0)<2
        a1=zeros(size(instAmp0));
        a2=a1;
        f1=instFreq0;
        f2=instFreq0;
        ratio_bw=a1;
        bis_freq=zeros(size(instAmp0));
        avgFreq=zeros(size(instAmp0));
        return;
    end
    %envpmax_instAmp = pchip(indmax_amp0,instAmp0(indmax_amp0),1:length(instAmp0));
    %envpmin_instAmp = pchip(indmin_amp0,instAmp0(indmin_amp0),1:length(instAmp0));
    envpmax_instAmp = interp1(indmax_amp0,instAmp0(indmax_amp0) ,tt,'pchip');
    envpmin_instAmp = interp1(indmin_amp0,instAmp0(indmin_amp0) ,tt,'pchip');

    a1 = (envpmax_instAmp+envpmin_instAmp)/2;  
    a2 = (envpmax_instAmp-envpmin_instAmp)/2;
    [indmin_a2, indmax_a2] = extr(a2);
   
    %%%%%%%%%%% obtain phi1(t) and phi2(t)
    instAmpinstAmp2=instFreq0.*instAmp0.^2;
    instAmp_tmax = pchip(indmax_amp0,instAmpinstAmp2(indmax_amp0),1:length(instAmpinstAmp2));
    instAmp_tmin = pchip(indmin_amp0,instAmpinstAmp2(indmin_amp0),1:length(instAmpinstAmp2));
    f1=zeros(size(instFreq0));
    f2=zeros(size(instFreq0));
    for i=1:numel(instFreq0)
        A = [a1(i).^2+a1(i).*a2(i) a2(i).^2+a1(i).*a2(i);
           a1(i).^2-a1(i).*a2(i) a2(i).^2-a1(i).*a2(i)];
        B = [instAmp_tmax(i)
            instAmp_tmin(i)];
        C=A\B;
        f1(i)=C(1);
        f2(i)=C(2);
    end
  
    
    bis_freq=(instAmp_tmax-instAmp_tmin)./(4*a1.*a2);
    if numel(indmax_a2)>3
        bis_freq = interp1(indmax_a2,bis_freq(indmax_a2) ,tt,'pchip');
    end
  
    avgFreq=(instAmp_tmax+instAmp_tmin)./(2*(a1.^2+a2.^2));
     %%%%%%%%%%% obtain instantaneous bandwidth
    cos_diffphi=(instAmp0.^2-a1.^2-a2.^2)./(2*a1.*a2);
    cos_diffphi(cos_diffphi>1.2)=1;
    cos_diffphi(cos_diffphi<-1.2)=-1;
    [instAmp1,instFreq_diff_phi] = INST_FREQ_local(cos_diffphi);

    
    diff_a1 = (a1(3:end) - a1(1:end-2))/2;
    diff_a1 = [diff_a1(1) diff_a1 diff_a1(end)];
    diff_a2 = (a2(3:end) - a2(1:end-2))/2;
    diff_a2 = [diff_a2(1) diff_a2 diff_a2(end)];
    
    instBW=((diff_a1.^2+diff_a2.^2)./(a1.^2+a2.^2)+(a1.^2).*(a2.^2).*(instFreq_diff_phi.^2)./((a1.^2+a2.^2).^2)).^0.5;
    ratio_bw=abs(instBW./avgFreq);
    ratio_bw(a2./a1<5e-3)=0;
    ratio_bw(avgFreq<1e-7)=0;
    ratio_bw(ratio_bw>1)=1;
       
    ff1=(instFreq_diff_phi+2*bis_freq)/2;  %in case a1=a2
    ff2=(2*bis_freq-instFreq_diff_phi)/2;
    f1(abs((a1-a2)./a1)<0.05)=ff1(abs((a1-a2)./a1)<0.05);
    f2(abs((a1-a2)./a1)<0.05)=ff2(abs((a1-a2)./a1)<0.05);

    temp_instAmp0=instAmp0;
    for j=1:numel(indmax_y)-1
      ind=indmax_y(j):indmax_y(j+1);
      temp_instAmp0(ind)=mean(instAmp0(ind));
    end
    ratio_bw(abs(temp_instAmp0)./max(abs(y))<5e-2)=0;
    f1(abs(temp_instAmp0)./max(abs(y))<4e-2)=1/numel(y)/1000; % very small freq
    f2(abs(temp_instAmp0)./max(abs(y))<4e-2)=1/numel(y)/1000;
    bis_freq(bis_freq>0.5)=0.5;
    bis_freq(bis_freq<0)=0;
end






function padded = extendsignal(y,nPadding)
%PADSIGNAL Pads a signal to decrease border effects
%   INPUTS:
%       y - A real vector that will be extended 
%       nPadding - number of points that are going to be added to one
%       side of the signal.
%   OUTPUTS:
%       padded - A version of y that has been extended
    leftPadded = 0;
    rightPadded = 0;
    nAddedEnd = 0;
    nAddedStart = 0;
    padded = y;
while(leftPadded ==0 || rightPadded==0)    
    nPoints = numel(y); 
    endSlope = y(nPoints)-y(nPoints-1);
    endInds = FastCrossing(y,y(nPoints));% find points index has the same level 
    nAddEnd = nPoints-1;
    endInds(endInds == numel(y)) = [];
    nAddEnd = MatchSlope(y,endInds,endSlope)+1;
    if isempty(nAddEnd)      
        if(rightPadded == 0)
            if(nPoints-1+nAddedEnd>=nPadding) % no points has the same level as end point, consider a symmetric padding
                rightPadded = 1;
            end 
            nAddEnd2 = max(1,nPoints-1-nPadding+1);
            padded = [padded y(end-1:-1:nAddEnd2) ];
            nAddedEnd = length(y(end-1:-1:nAddEnd2))+nAddedEnd;
        end
    end

    startSlope = y(2)-y(1);

    startInds = FastCrossing(y,y(1));% find points index has the same level 
    nAddStart = 2;
    startInds(startInds+1 == 2) = [];
    nAddStart = MatchSlope(y,startInds+1,startSlope)-2;
    if isempty(nAddStart) 
        if(leftPadded == 0)
            if(nPoints-1+nAddedStart>=nPadding) % no points has the same level as 1st point, consider a symmetric padding
                leftPadded = 1;
            end     
            nAddStart2 = min(nPoints,nPadding+1);
            padded = [y(nAddStart2:-1:2) padded];
            nAddedStart = nAddedStart+length(y(nAddStart2:-1:2));
        end      
    end
    if(leftPadded)
        nAddStart2 = nAddStart;
    else
        nAddStart2 = 1;
    end
    if(nAddStart+nAddedStart>=nPadding)
        nAddStart2 = nAddStart+nAddedStart-nPadding+1;
        leftPadded = 1;
    end
    if(rightPadded)
        nAddEnd2 = nAddEnd;
    else      
        nAddEnd2 = nPoints;
    end
    if(nPoints-nAddEnd+1+nAddedEnd>=nPadding)
        nAddEnd2 = nAddEnd+nPadding-1-nAddedEnd;
        rightPadded = 1;
    end
    padded = [y(nAddStart2:1:nAddStart) padded y(nAddEnd:1:nAddEnd2)];
    nAddedEnd = nAddedEnd+nAddEnd2-nAddEnd+1;
    nAddedStart = nAddedStart+nAddStart-nAddStart2+1;
end
end

% Calculate when S crosses level
function ind = FastCrossing(s,level)
    s   = s - level;
    ind = zerocross(s);
end

function outInd = MatchSlope(y,inds,slope)
    nPoints = numel(y); 
    inds((y(inds)-y(inds-1))*slope<0)=[];
    inds(inds>nPoints-5)=[];
    inds(inds<5)=[];
    slopes = y(inds)-y(inds-1);
    [slopeDiff,ind] = min(abs(slopes-slope));
    if ~isempty(slopeDiff)    
        outInd = inds(ind);
    else
        outInd = [];
    end
end


%extracts the indices of extrema
% Writtezn by Gabriel Rilling
function [indmin, indmax, indzer] = extr(x,t)

if(nargin==1)
  t=1:length(x);
end

m = length(x);

if nargout > 2
  x1=x(1:m-1);
  x2=x(2:m);
  indzer = find(x1.*x2<0);

  if any(x == 0)
    iz = find( x==0 );
    indz = [];
    if any(diff(iz)==1)
      zer = x == 0;
      dz = diff([0 zer 0]);
      debz = find(dz == 1);
      finz = find(dz == -1)-1;
      indz = round((debz+finz)/2);
    else
      indz = iz;
    end
    indzer = sort([indzer indz]);
  end
end

d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;


% when two or more successive points have the same value we consider only one extremum in the middle of the constant area
% (only works if the signal is uniformly sampled)

if any(d==0)

  imax = [];
  imin = [];

  bad = (d==0);
  dd = diff([0 bad 0]);
  debs = find(dd == 1);
  fins = find(dd == -1);
  if debs(1) == 1
    if length(debs) > 1
      debs = debs(2:end);
      fins = fins(2:end);
    else
      debs = [];
      fins = [];
    end
  end
  if length(debs) > 0
    if fins(end) == m
      if length(debs) > 1
        debs = debs(1:(end-1));
        fins = fins(1:(end-1));

      else
        debs = [];
        fins = [];
      end
    end
  end
  lc = length(debs);
  if lc > 0
    for k = 1:lc
      if d(debs(k)-1) > 0
        if d(fins(k)) < 0
          imax = [imax round((fins(k)+debs(k))/2)];
        end
      else
        if d(fins(k)) > 0
          imin = [imin round((fins(k)+debs(k))/2)];
        end
      end
    end
  end

  if length(imax) > 0
    indmax = sort([indmax imax]);
  end

  if length(imin) > 0
    indmin = sort([indmin imin]);
  end

end
end
