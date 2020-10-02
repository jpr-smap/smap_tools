function d=pairwiseQD(q1,varargin);
% 031417: no longer calls external function

if( (size(q1,1)~=4) | (isempty(strfind(class(q1),'double'))) )
    if( (size(q1,1)==3) & isa(q1,'quaternion')==0 )
        %q1=normalizeRM(q1);
        q1=squeeze(double(quaternion.rotationmatrix(q1)));
    elseif( strcmp(class(q1),'quaternion') )
        q1=squeeze(double(q1));
    else
        fprintf('format not recognized...\n');
    end;
end;
if( nargin>1 )
    q2=varargin{1};
    if( (size(q2,1)~=4) | (isempty(strfind(class(q2),'double'))) )
        if( (size(q2,1)==3) & isa(q2,'quaternion')==0 )
            %q2=normalizeRM(q2);
            q2=squeeze(double(quaternion.rotationmatrix(q2)));
        elseif( strcmp(class(q2),'quaternion') )
            q2=squeeze(double(q2));
        else
            fprintf('format not recognized...\n');
        end;
    end;
else
    q2=q1;
end;

pfact=2.*180./pi;

d=zeros(size(q1,2),size(q2,2),'double');
for i=1:size(q1,2)
%     d(i,:)=smap.measureQD(q1(i),q2); 
%     d(i,:)=2.*real(acos(abs(squeeze(sum(bsxfun(@times,q1(:,i),q2)))))).*180./pi;
    
    dq=arrayfun(@real,arrayfun(@acos,complex(arrayfun(@abs,sum(bsxfun(@times,q1(:,i),q2))))));
    d(i,:)=bsxfun(@times,pfact,dq);

end;

if( nargin>1 )
else
%     dummy=ones(size(d,1),size(d,2));
%     dummy_nan=nan(size(d,1),size(d,2));
%     dummy=dummy.*diag(nan(size(d,1),1),0);
%     dummy=dummy.*tril(dummy_nan);
%     dummy(~isnan(dummy))=1;
%     d=d.*dummy;
%     %d(find(d==0))=nan;
%     %d=round((d));
    d=(1-(triu(d).*tril(d))).*triu(d);
end;



