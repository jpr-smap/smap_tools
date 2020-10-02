function thr_ref=PR_quick(tableref1,varargin);

if( nargin>1 )
%     tableref2=varargin{1};
    prec_thrs=varargin{1};
else
%     tableref2=tableref1;
    prec_thrs=[0.99 0.95 0.9 0.85];

end;

tableref=tableref1;

val_min=min([tableref.peak_part_ctrl; tableref.peak_part]);
val_max=max([tableref.peak_part_ctrl; tableref.peak_part]);
    
xx=linspace(val_min,val_max,5e3);
peaks=hist(tableref.peak_part,xx);
peaks_ctrl=hist(tableref.peak_part_ctrl,xx);

% peaks_s=sort(tableref1.peak_part,'ascend');
% peaks_s_ctrl=sort(tableref1.peak_part_ctrl,'ascend');
% nAbove=zeros(1,length(peaks_s));
% for i=1:length(peaks_s)
%     nAbove(i)=length(find(peaks_s>peaks_s_ctrl(i)));
% end;
% nAbove_ctrl=fliplr([1:length(peaks_s_ctrl)]);
cs=sum(peaks)-cumsum(peaks);
cs_ctrl=sum(peaks_ctrl)-cumsum(peaks_ctrl);

prec=cs./(cs+cs_ctrl);
recall=cs./max(cs);
plot(recall,prec);
thr_ref=[];
for i=1:length(prec_thrs)
    try
        thr_ref(i)=xx(find(prec>prec_thrs(i),1,'first'));
    catch
        fprintf('no thr. found for precision %3.2f\n',prec_thrs(i));
    end;
end;
if( nargin>1 )
    hold on;
    plot(cs./max(cs),cs./(cs+cs_ctrl));
end;
grid on;
xlabel('est. recall'); ylabel('est. precision');

% clf
% plot(xx,peak_old_ctrl);
% hold on;
% plot(xx,cs_old)
% plot(xx,cs_new);
% figure
% plot(xx,cs_old./(cs_old+cs_old_ctrl));
% hold on;
% plot(xx,cs_new./(cs_new+cs_mew_ctrl));
% plot(xx,cs_new./(cs_new+cs_new_ctrl));

% val_min=min([tableref1.peak_part_ctrl; tableref1.peak_part; tableref2.peak_part_ctrl; tableref2.peak_part]);
% val_max=max([tableref1.peak_part_ctrl; tableref1.peak_part; tableref2.peak_part_ctrl; tableref2.peak_part]);

%[peak_old,sI]=sort(tableref1.peak_part,'ascend');
% peak_new_ctrl=hist(tableref2.peak_part_ctrl,xx);
% peak_new=hist(tableref2.peak_part,xx);
% cs_new=sum(peak_new)-cumsum(peak_new);
% cs_new_ctrl=sum(peak_new_ctrl)-cumsum(peak_new_ctrl);
% plot(prec,cs_old./(cs_old+cs_old_ctrl));
% prec=cs_old./max(cs_old);
