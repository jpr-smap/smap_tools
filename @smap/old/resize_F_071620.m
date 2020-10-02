function outref=resize_F(inref,sf,varargin);
% Fourier crop or pad to a new size; use 'newSize' to scale boundaries, too
% skip varargin to keep the original size;
% 
% outref=resize_F(inref,sf,varargin);

method='fixedSize';
if( nargin>2 )
    method=varargin{1};
end;

sfAmp=abs(1-sf);
if( sfAmp>0.2 )
    sfType='large';
else
    sfType='small';
end;

if( sf>1 ) % Fourier pad (upsample/zoom in):
    os=[size(inref,1) size(inref,2) size(inref,3)];
    bgVal=mode(inref(:));
    finalSize=floor([size(inref,1) size(inref,2)].*sf);
    inref_F=fftshift(fftn(ifftshift(inref)));
    if( size(inref,3)==1 )
        disp('downsample');
        
%         outref=resize_2d(inref_F,finalSize,'up');

        inref_F_padded=smap.cropOrPad(inref_F,[finalSize(1),finalSize(end)],0);
        outref=real(fftshift(ifftn(ifftshift(inref_F_padded))));

        outref=outref.*(sf).^2;
%        
        switch method
            case 'fixedSize'
                outref=smap.cutj(outref,[os(1),os(2)]);
            case 'newSize'
                outref=outref;
            otherwise
                outref=outref;
        end;
        if( strcmp(sfType,'small')==1 )
%             tempCC=smap.ccf(inref,outref); % % XX
%             peakVal=max(tempCC(:));
%             if( peakVal<100 )
%                 disp('centering...');
%             end;
%             ctr=1;
%             while( peakVal<100 )
%                 fs=os(1)+ctr;
%                 tempRs=resize_2d(inref_F,fs,'up');
%                 tempRs=tempRs.*(sf).^2;
%                 tempRs=smap.cutj(tempRs,[os(1),os(2)]);
%                 tempCC=smap.ccf(tempRs,outref);
%                 peakVal=max(tempCC(:)); % % XX
%                 ctr=ctr+1;
%             end;
%             ind=find(tempCC==peakVal,1,'first');
%             [x,y]=ind2sub(size(inref),ind);
%             cp=smap.getcp(inref); % % XX
%             nxny=smap.maxInterpF(tempCC,16,20); % % XX
%             if( max(abs(nxny)>0 ) )
%                 outref=smap.applyPhaseShifts(outref,-fliplr(nxny)); % % XX
%             end;
        end;
    else
        inref_F_padded=smap.cropOrPad(inref_F,[finalSize,finalSize,finalSize],bgVal);
        outref=real(fftshift(ifftn(ifftshift(inref_F_padded))));
%         outref=smap.iftj(inref_F_padded);
%        
        outref=outref.*(1./sf).^3;
%        
        if( strcmp(method,'fixedSize') )
            if( size(outref,1)<os(1) )
                outref=smap.cropOrPad(outref,os,bgVal);
            elseif( size(outref,1)>os(1) )
                outref=smap.cropOrPad(outref,os);
            end;

        end;

    end;
    
else % Fourier crop (downsample/zoom out):
    if( size(inref,3)==1 )
        disp('upsample');
        
        os=[size(inref,1) size(inref,2)];
        dsSize=[];
        
        finalSize=[]; dsSize=[];
        
        for k=1:2
            finalSize(k)=size(inref,k);
            dsSize(k)=finalSize(k).*sf;
            while( mod(dsSize(k),1)>1e-2 )%1e-3 ) % % 062117
                finalSize(k)=finalSize(k)+1;
                dsSize(k)=finalSize(k).*sf;
            end;
        end;
        dsSize
%         finalSize
        finalSize=floor([size(inref,1) size(inref,2)].*sf);

%         inref=smap.cropOrPad(inref,finalSize,mode(inref(:)));
        inref_F=fftshift(fftn(ifftshift(inref)));        
%         outref=resize_2d(inref_F,dsSize,'down');
        
%         inref_F_cut=double(smap.cropOrPad(inref_F,os));%[finalSize(1),finalSize(end)]));
        inref_F_cut=double(smap.cropOrPad(inref_F,[finalSize(1),finalSize(end)]));
        outref=real(fftshift(ifftn(ifftshift(inref_F_cut))));

                
        switch method
            case 'fixedSize'
                outref=smap.cropOrPad(smap.cutj(outref,[os(1),os(2)]),[os(1),os(2)],mean(outref(:)));
            otherwise
                outref=outref;
        end;
%        
        outref=outref.*(sf.^2);
%

        if( strcmp(sfType,'tiny')==1 )
            tempCC=smap.ccf(inref,outref); % % XX
            peakVal=max(tempCC(:));
            if( peakVal<100 )
                disp('centering...');
            end;
            ctr=1;
            while( peakVal<100 )
                fs=os(1)-ctr;
                tempRs=resize_2d(inref_F,fs,'down');
%                 tempRs=tempRs.*(sf).^2;
                tempRs=smap.extendj(tempRs,[os(1),os(2)],tempRs(1,1));
                tempCC=smap.ccf(tempRs,outref); % % XX
                peakVal=max(tempCC(:));
                ctr=ctr+1;
            end;
            ind=find(tempCC==peakVal,1,'first');
            [x,y]=ind2sub(size(inref),ind);
            cp=smap.getcp(inref);
            nxny=smap.maxInterpF(tempCC,16,20);
            if( max(abs(nxny)>0 ) )
                outref=smap.applyPhaseShifts(outref,-fliplr(nxny)); % % XX
            end;
        end;
        
        
    else
        
        os=[size(inref,1) size(inref,2) size(inref,3)];
        bgVal=mode(inref(:));
%         finalSize=max([size(inref,1) size(inref,2) size(inref,3)]);
        
        finalSize=[]; dsSize=[];
        for k=1:3
            finalSize(k)=size(inref,k);
            dsSize(k)=finalSize(k).*sf;
            while( mod(dsSize(k),1)>1e-2 )%1e-3 ) % % 062117
                finalSize(k)=finalSize(k)+1;
                dsSize(k)=finalSize(k).*sf;
            end;
        end;

        disp(['will pad original volume to ' num2str(finalSize) ' voxels...']);
        disp(['extending to full size...']);
%         inref=smap.extendj(inref,[finalSize,finalSize,finalSize],inref(1,1,1));
        inref=smap.cropOrPad(inref,finalSize,bgVal);
        
        ns=[size(inref,1).*sf  size(inref,2).*sf size(inref,3).*sf];
        
        disp('forward FFT...');
        inref=fftshift(fftn(ifftshift(inref)));
        disp('cropping...');
        inref=single(smap.cropOrPad(inref,dsSize));%[ns(2),ns(1),ns(3)]));
        disp('reverse FFT...');
        outref=real(fftshift(ifftn(ifftshift(inref))));
%        
        outref=outref.*(sf.^3);
        if( strcmp(method,'fixedSize') )
            if( size(outref,1)<os(1) )
                outref=smap.cropOrPad(outref,os,bgVal);
            elseif( size(outref,1)>os(1) )
                outref=smap.cropOrPad(outref,os);
            end;
        end;
%

    end;
    
end;

    function outref=resize_2d(inref_F,finalSize,resampleDir);
        switch resampleDir
            case 'up'
                inref_F_padded=smap.cropOrPad(inref_F,[finalSize(1),finalSize(end)],0);
                outref=real(fftshift(ifftn(ifftshift(inref_F_padded))));
%                 outref=smap.iftj(inref_F_padded);
            case 'down'
                inref_F_cut=double(smap.cropOrPad(inref_F,[finalSize(1),finalSize(end)]));
                outref=real(fftshift(ifftn(ifftshift(inref_F_cut))));
%                 outref=smap.iftj(inref_F_cut);
        end;
%     end;
    
%     function outref=smap.getcp(inref);
%         outref=[floor(size(inref,1)./2)+1 floor(size(inref,2)./2)+1];
%     end;
    