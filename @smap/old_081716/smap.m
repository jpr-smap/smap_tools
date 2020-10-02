classdef smap < handle
    
    %% properties:
    properties
        computerID
        baseDir
        smapID
        info
        im
        status
        search
    end;
    %     events
    %       stepDone
    %     end
    
    %% ordinary methods:
    methods
        %% constructor:
        function obj=smap(tempID,varargin)
            if( nargin>1 )
                fn=varargin{1};
            end;
            obj.computerID=smap.whoami;
            switch obj.computerID
                case '.B945E8C52A98'
                    %                     disp('helium');
                    obj.baseDir='/Volumes/denk';
                case '.DD81502405EA.E8EFCF40A15A'
                    %                     disp('neon');
                    obj.baseDir='/tier2/denk';
            end;
            
            masterList=[obj.baseDir '/masterList.mat'];
            load(masterList,'s');
            filePath={}; ID={};
            [filePath{1:length(s),1}] = deal(s.filePath);
            [ID{1:length(s),1}] = deal(s.ID);
            try
                idx = find(~cellfun('isempty',strfind({s.ID},tempID))==1);
                disp('loading file...');
                load([obj.baseDir s(idx).filePath],'t');
                obj=t;
                obj.computerID=smap.whoami; % in case loading the instance from a different computer
                switch obj.computerID
                    case '.B945E8C52A98'
                        %                                                 disp('helium');
                        obj.baseDir='/Volumes/denk';
                    case '.DD81502405EA.E8EFCF40A15A'
                        %                                                 disp('neon');
                        obj.baseDir='/tier2/denk';
                end;
            catch
                disp('does not exist...creating new smap instance');
                obj.smapID=tempID;
                disp(['skipping read of' fn]);
                obj.im.image=[];%tr(fn);
                obj.info.Npix=[size(obj.im.image,1) size(obj.im.image,2)];
                obj.info.gainRef=[];%'/images/gainRef/K2-0001 1 Gain Ref. x1.m3.kv[300].dm4';
                obj.info.rawFile=strrep(fn,obj.baseDir,'');
                obj.status='';
                obj.search=[];
                t=obj;
                sl=length(s)+1;
                s(sl).ID=tempID;
                s(sl).filePath=['/obj/' tempID '.mat']
                
                [obj.baseDir s(sl).filePath]
                save([obj.baseDir s(sl).filePath],'t');
                save(masterList,'s');
            end;
        end;
        
        %% gain-correction:
        function obj=gainCorr(obj)
            
            basedir=obj.baseDir;
            gainref_fn=[basedir obj.info.gainRef];
            inref_fn=[basedir obj.info.rawFile];
            outref_fn=[basedir '/images/preproc/' obj.smapID '_gc.mrc'];
            
            disp('loading gain ref...');
            [m,sx,units]=ReadDMFile(gainref_fn);
            disp('loading raw image file...');
            tempOrig=tr(inref_fn);
            tempOrig=tempOrig(:,:,2:end);
            shortEdge=min([size(tempOrig,1) size(tempOrig,2)]);
            
            tempOrig=smap.cutj(tempOrig,[shortEdge,shortEdge,size(tempOrig,3)]);
            
            temp=zeros(size(tempOrig),'single');
            m=smap.cutj(m,[shortEdge,shortEdge]);
            
            gcAll=zeros(size(tempOrig),'single');
            
            for j=1:size(tempOrig,3)
                temp(:,:,j)=tempOrig(:,:,j)';
            end;
            
            tempDs=[];
            for j=1:size(temp,3)
                gcAll(:,:,j)=temp(:,:,j).*m;
                tempDs(:,:,j)=smap.resize_F(gcAll(:,:,j),0.125,'newSize');
                tempDiff=abs(diff(tempDs(:,:,j),1)); % diff along dim 2 if im not rotated/flipped
                nZeros(j)=length(find(tempDiff==0));
            end;
            
            badFrames=find(nZeros>100);
            if( length(badFrames)>0 )
                disp(['found ' num2str(length(badFrames)) 'bad frames...']);
            end;
            
            goodFrames=setdiff(1:size(gcAll,3),badFrames)
            %             goodFrames=1:size(gcAll,3);
            gcAll=gcAll(:,:,goodFrames);
            %             tempDsToWrite=sum(tempDs(:,:,goodFrames),3);
            temp2=sum(gcAll,3);
            
            clear temp;
            
            gs=sum(gcAll,3);
            [nx_hot,ny_hot]=find(nm(gs)>7);
            nx=nx_hot; ny=ny_hot;
            allPix=[nx ny];
            disp(['replacing hot pixels... (found ' num2str(length(nx_hot)) ')']);
            
            imOut=zeros(shortEdge,shortEdge,size(gcAll,3),'single');
            for j=1:size(gcAll,3)
                movref=gcAll(:,:,j);
                outIm=movref;
                for i=1:size(allPix,1)
                    rs=[(allPix(i,1)-1):(allPix(i,1)+1)];
                    cs=[(allPix(i,2)-1):(allPix(i,2)+1)];
                    if(min(rs)<1)
                        rs=rs-min(rs)+1;
                    end;
                    if(min(cs)<1)
                        cs=cs-min(cs)+1;
                    end;
                    if( max(rs)>size(movref,1) )
                        rs=rs-(max(rs)-size(movref,1));
                    end;
                    if( max(cs)>size(movref,2) )
                        cs=cs-(max(cs)-size(movref,2));
                    end;
                    tempPatch=movref(rs,cs);
                    mv=mean(tempPatch(:));
                    outIm(allPix(i,1),allPix(i,2))=mv;
                end;
                disp(j);
                imOut(:,:,j)=outIm;
            end;
            disp('Saving gain-corrected mov..');
            obj.im.gaincorr=sum(imOut,3);
            mw(imOut,outref_fn,32);
            obj.im.temp=[];
            
            obj.status='gain-corrected';
            obj.info.gainCorr_fn=['/images/preproc/' obj.smapID '_gc.mrc'];
            [~]=updateObj(obj);
            disp('done');
        end;
        
        %% motion-correction:
        function obj=motionCorr(obj,varargin);
            
            fprintf('Getting GPU...\n');
            gdev=gpuDevice(1);
            
            warning off;
            
            basedir=obj.baseDir;
            inref_fn=[basedir obj.info.gainCorr_fn];
            outref_fn=strrep(inref_fn,'gc.mrc','mc.mrc');
            outref_sum_fn=strrep(inref_fn,'gc.mrc','mc_fullSum.mrc');
            outref_temp_fn=strrep(inref_fn,'gc.mrc','mc_temp.mrc');
            fnOut=strrep(inref_fn,'gc.mrc','shifts.mat');
            
            binVal=4;
            maxShift=20/binVal; % search region
            
            % % read in file, get properties and cut to a square:
            fprintf('Reading gain-corrected movie...\n');
            movref=mr(inref_fn);
            nFrames=size(movref,3);
            temp=smap.resize_F(movref(:,:,1),1./binVal,'newSize');
            temp=zeros(size(temp,1),size(temp,2),size(movref,3));
            if( binVal>1 )
                for j=1:nFrames
                    temp(:,:,j)=smap.resize_F(movref(:,:,j),1./binVal,'newSize');
                end;
                movref=single(temp);
            end;
            
            minDim=min([size(movref,1) size(movref,2)]);
%             movref=smap.cutj(movref,[minDim,minDim,size(movref,3)]);
            temp=movref; movref=zeros(minDim,minDim,nFrames,'single');
            for j=1:nFrames
                movref(:,:,j)=nm(cutj(temp(:,:,j),[minDim,minDim]));
            end;
            
            %             tmLP=variableCosMask(minDim,[0.44 0.46],aPerPix);
            
            % % get initial reference frame (standard: 10-30 for 40 Hz):
            fprintf('Finding reference frame...\n');
            searchFrames=[10 size(movref,3)-10];
            % searchFrames=[1 9];
            testP=zeros(1,nFrames); Y=[]; temp=[];
            for j=searchFrames(1):searchFrames(2)
                temp=ccf(movref(:,:,j),movref(:,:,j+1));
                testP(j)=max(temp(:));
            end;
            refNum=find(testP==max(testP))+1
            
            doneFrames=zeros(1,nFrames);
            shifts=zeros(nFrames,2);
            refImOrig=movref(:,:,refNum);
            peaks=zeros(1,nFrames);
            peaksBest=peaks;
            nDone=1;
            
            refIm=refImOrig;
            
            inputMov=movref;
            movref=gpuArray(movref); wait(gdev);
            movNew=movref;
            
            fprintf('Starting motion-correction...\n');
            movref_F=gpuArray(complex(zeros(size(movref),'single'),zeros(size(movref),'single')));
            for j=1:nFrames
                movref_F(:,:,j)=smap.ftj(gpuArray(movref(:,:,j)));
            end;
            ctr=1; pbd=[];
            shiftsAll=zeros(nFrames,2);
            pbd=[1e7 1e8];
            shiftsByRound=[];
            annealingMode=0;
            while( ((ctr<100)&(abs(pbd(end)-pbd(1:end-1))>0))|(ctr<=nFrames) )
                if( ctr>=nFrames )
                    annealingMode=1;
                end;
                
                for j=1:nFrames
                    refIm_F(:,:,j)=smap.ftj(nm(sum(movNew(:,:,setdiff(1:nFrames,j)),3)));
                end;
                
                shifts=zeros(nFrames,2);
                peakVals=zeros(1,nFrames);
                for j=1:nFrames
                    cc=gather(smap.iftj(refIm_F(:,:,j).*conj(movref_F(:,:,j)))); wait(gdev);
                    [shifts(j,:),peakVals(j)]=smap.maxInterpF(cc,maxShift,10*binVal);
%                     if( mod(j,10)==0 )
%                         fprintf('%d...',j);
%                     end;
                end;
                
                % apply and store:
                peaksNew=peakVals;
                
                pb=[peaksBest; peaksNew];
                peaksBest=max(peaksBest,peaksNew);
                improvedPeaks=find(pb(2,:)>pb(1,:));
                pbd(ctr)=sum(abs(pb(2,:)-pb(1,:)));
                if( nDone==1 )
                    [X,Y]=sort(peaksNew,'descend');
                    newFrames=Y(1:2); % ref frame and next-best frame
                    shiftsAll(newFrames,:)=shifts(newFrames,:);%.*binVal;
                else
                    peaksEx=(1-doneFrames).*peaksNew;
                    newFrames=find(peaksEx==max(peaksEx));
                    dpbd=[0 diff(pbd)];
%                     if( dpbd(end)<dpbd(end-1) )
%                         improvedPeaks=[];
%                     end;
                end;
                if( ~isempty( newFrames ) )
                    doneFrames(newFrames)=1;
                    goodFrames=find(doneFrames==1);
                    shiftsAll(newFrames,:)=shifts(newFrames,:);%.*binVal;
                    shiftsToUpdate=intersect(improvedPeaks,goodFrames);
                    shiftsAll(shiftsToUpdate,:)=shifts(shiftsToUpdate,:);
                    %                     if( isempty(shiftsToUpdate)==1 )
                    %                         temp=smap.applyPhaseShifts(movref(:,:,newFrames),-fliplr(shiftsAll(newFrames,:)));
                    %                         movNew(:,:,newFrames)=temp;
                    %                     else
                    %                         temp=smap.applyPhaseShifts(movref(:,:,shiftsToUpdate),-fliplr(shiftsAll(shiftsToUpdate,:)));
                    %                         for j=1:size(shiftsToUpdate,2)
                    %                             movNew(:,:,shiftsToUpdate(j))=temp(:,:,j);
                    %                         end;
                    %                     end;
                    movNew=gpuArray(smap.applyPhaseShifts(inputMov,-fliplr(shiftsAll))); wait(gdev);
                end;
                
                figure(101); subplot(1,2,1);
                plot(shiftsAll,'x-'); 
                xlim([-0.5 nFrames+0.5]); ylim([-maxShift maxShift]); pause(0.05);
                
                shiftsByRound{ctr}=shiftsAll;
                
                nDone=length(find(doneFrames==1));
                fprintf('%d %f\n',[nDone pbd(ctr)]);
                mw(single(refIm),outref_temp_fn,32);
                save(fnOut,'shiftsByRound');
                ctr=ctr+1;
%                 toc;
            end;
            fprintf('...done.\n');
            
            shiftsFinal=shiftsAll-repmat(shiftsAll(refNum,:),size(shifts,1),1)
            
%             movref=tr(inref_fn);
%             nFrames=size(movref,3);
%             temp=smap.resize_F(movref(:,:,1),1./binVal,'newSize');
%             temp=zeros(size(temp,1),size(temp,2),size(movref,3));
%             for j=1:nFrames
%                 temp(:,:,j)=smap.resize_F(movref(:,:,j),1./binVal,'newSize');
%             end;
%             movref=temp;
            
            movref=mr(inref_fn);
            ftu=1:size(movref,3)
            movShifted=single(smap.applyPhaseShifts(movref(:,:,ftu),-fliplr(shiftsFinal).*binVal));
            
            minDim=min([size(movShifted,1) size(movShifted,2)]);
            temp=movShifted; movShifted=zeros(minDim,minDim,nFrames,'single');
            for j=1:nFrames
                movShifted(:,:,j)=cutj(temp(:,:,j),[minDim,minDim]);
            end;
            
            movShifted=movShifted(:,:,ftu);
            temp=movShifted; clear movShifted;
            for j=1:nFrames
                disp(j);
                movShifted(:,:,j)=resize_F(temp(:,:,j),0.5,'newSize');
            end;
            
            fullSum=sum(movShifted,3);
            fprintf('Writing motion-corrected files...\n');
            mw(single(fullSum),outref_sum_fn,32);
            mw(single(movShifted),outref_fn,32);
            
            obj.im.fullSum=fullSum;
            obj.im.motionTracking.shiftsAll=shiftsAll;
            obj.im.motionTracking.shiftsFinal=shiftsFinal;
            obj.im.motionTracking.shiftsByRound=shiftsByRound;
            obj.status='motion-corrected';
            obj.info.motionCorr_fn=strrep(outref_fn,obj.baseDir,'');
            [~]=updateObj(obj);
            
            
            disp('done');
            
            
        end;
        
        %% update object:
        function obj=updateObj(obj)
            disp(['updating instance...']);
            t=obj;
            fnOut=[obj.baseDir '/obj/' obj.smapID '.mat'];
            try
                movefile(fnOut,[obj.baseDir '/obj/bak/']);
            catch
                disp('writing new object...');
            end;
            save(fnOut,'t');
        end;
        
        %% dummy:
        function obj=dummy(obj);
            
        end;
        
    end;
    

    
    %% static methods
    
    methods (Static)
        %% whoami:
        function outref=whoami(varargin)
            sid = '';
            ni = java.net.NetworkInterface.getNetworkInterfaces;
            while ni.hasMoreElements
                addr = ni.nextElement.getHardwareAddress;
                if ~isempty(addr)
                    addrStr = dec2hex(int16(addr)+128);
                    sid = [sid, '.', reshape(addrStr,1,2*length(addr))];
                end
            end
            outref=sid;
        end;
        
        %% ftj:
        function outref=ftj(inref)%,varargin)
            Npix=prod(size(inref));
            inref_F=fftshift(fftn(ifftshift(inref)))./sqrt(Npix);
            outref=inref_F;
        end;
        
        %% iftj:
        function outref=iftj(inref)
            Npix=prod(size(inref));
            inref(find(isnan(inref)==1))=0;
            outref=fftshift(real(ifftn(ifftshift(inref)))).*sqrt(Npix);
        end;
        
        %% cutj:
        function outref=cutj(inref,varargin);
            nDims=length(size(inref));
            for i=1:nDims
                oldDim(i)=size(inref,i);
                newDim(i)=varargin{1}(i);
                dimL(i)=(newDim(i)>oldDim(i));
            end;
            
            centerPixOld=smap.getcp(inref);
            for i=1:nDims
                if( dimL(i)==0 )
                    oddOldFlag=mod(oldDim(i),2);
                    if( oddOldFlag==1 )
                        halfOldDim(i)=ceil(oldDim(i)./2);
%                         centerPixOld(i)=halfOldDim(i);
                    else
                        halfOldDim(i)=oldDim(i)./2;
%                         centerPixOld(i)=halfOldDim(i)+1;
                    end;
                    edges{i}=ceil([centerPixOld(i)-newDim(i)./2 centerPixOld(i)+(newDim(i)./2)-1]);
                else
                    edges{i}=[1 oldDim(i)];
                end;
            end;
            %disp(edges{1}); disp(edges{2});
            if( nDims<3 )
                outIm=inref(edges{1}(1):edges{1}(2),edges{2}(1):edges{2}(2));
            else
                for j=1:oldDim(3)
                    outIm(:,:,j)=inref(edges{1}(1):edges{1}(2),edges{2}(1):edges{2}(2),j);
                end;
                outIm=outIm(:,:,(edges{3}(1):edges{3}(2)));
            end;
            outref=outIm;
        end;
        
        %% extendj:
        function outref=extendj(inref,varargin);
            nDims=length(size(inref));
            for i=1:nDims
                oldDim(i)=size(inref,i);
                newDim(i)=varargin{1}(i);
            end;
            padVal=varargin{end};
            
            for i=1:nDims
                oddOldFlag=mod(oldDim(i),2);
                if( oddOldFlag==1 )
                    halfOldDim(i)=ceil(oldDim(i)./2);
                    halfNewDim(i)=floor(newDim(i)./2);
                    centerPixOld(i)=halfOldDim(i);
                    centerPixNew(i)=halfNewDim(i)+1;
                else
                    halfOldDim(i)=oldDim(i)./2;
                    halfNewDim(i)=newDim(i)./2;
                    centerPixOld(i)=halfOldDim(i)+1;
                    centerPixNew(i)=floor(halfNewDim(i)+1);
                end;
                edges{i}=[centerPixOld(i)-newDim(i)./2 centerPixOld(i)+(newDim(i)./2)-1];
            end;
            diffSize=newDim(1)-oldDim(1);
            
            if( nDims<3 )
                outref=ones(newDim(1),newDim(2),'single').*padVal;
                outref((centerPixNew(1)-oldDim(1)./2):(centerPixNew(1)+oldDim(1)./2)-1,(centerPixNew(2)-oldDim(2)./2):(centerPixNew(2)+oldDim(2)./2)-1)=inref;
            else
                outref=ones(newDim(1),newDim(2),newDim(3),'single').*padVal;
                for j=1:oldDim(3)
                    temp=ones(newDim(1),newDim(2)).*padVal;
                    temp(1:oldDim(1),1:oldDim(2))=inref(:,:,j);
                    outref(:,:,j)=temp;
                end;
                outref=circshift(outref,[centerPixNew(1)-centerPixOld(1),centerPixNew(2)-centerPixOld(2),centerPixNew(3)-centerPixOld(3)]);
            end;
        end;
        
        %% resize_F:
        function outref=resize_F(inref,sf,varargin);
            % outref=resize_F(inref,sf,varargin);
            % skip varargin to keep the original size;
            % use 'newSize' for varargin to scale the image boundaries, too
            
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
                finalSize=floor(min([size(inref,1) size(inref,2)]).*sf);
                inref_F=smap.ftj(inref);
                if( size(inref,3)==1 )
                    outref=resize_2d(inref_F,finalSize,'up');
                    outref=outref.*(sf).^2;
                    switch method
                        case 'fixedSize'
                            outref=smap.cutj(outref,[os(1),os(2)]);
                        case 'newSize'
                            outref=outref;
                        otherwise
                            outref=outref;
                    end;
                    if( strcmp(sfType,'small')==1 )
                        tempCC=ccf(inref,outref); % % XX
                        peakVal=max(tempCC(:));
                        if( peakVal<100 )
                            disp('centering...');
                        end;
                        ctr=1;
                        while( peakVal<100 )
                            fs=os(1)+ctr;
                            tempRs=resize_2d(inref_F,fs,'up');
                            tempRs=tempRs.*(sf).^2;
                            tempRs=smap.cutj(tempRs,[os(1),os(2)]);
                            tempCC=ccf(tempRs,outref);
                            peakVal=max(tempCC(:)); % % XX
                            ctr=ctr+1;
                        end;
                        ind=find(tempCC==peakVal,1,'first');
                        [x,y]=ind2sub(size(inref),ind);
                        cp=getcp(inref); % % XX
                        nxny=maxInterpF(tempCC,16,20); % % XX
                        if( max(abs(nxny)>0 ) )
                            outref=smap.applyPhaseShifts(outref,-fliplr(nxny)); % % XX
                        end;
                    end;
                else
                    inref_F_padded=smap.extendj(inref_F,[finalSize,finalSize,finalSize],0);
                    outref=smap.iftj(inref_F_padded);
                    outref=outref.*(1./sf).^3;
                end;
                
            else % Fourier crop (downsample/zoom out):
                if( size(inref,3)==1 )
                    
                    os=[size(inref,1) size(inref,2)];
                    finalSize=min([size(inref,1) size(inref,2)]);
                    dsSize=finalSize.*sf;
                    while( mod(dsSize,1)>1e-3 )
                        finalSize=finalSize+1;
                        dsSize=finalSize.*sf;
                    end;
                    finalSize=round(finalSize);
                    %         disp(['will pad original image to ' num2str(finalSize) ' pixels...']);
                    %         disp(['key to continue...']);
                    %         pause;
                    inref=single(smap.extendj(inref,[finalSize,finalSize],inref(1,1)));
                    
                    ns=round([size(inref,1).*sf size(inref,2).*sf]);
                    
                    inref_F=smap.ftj(inref);
                    
                    outref=resize_2d(inref_F,ns(1),'down');
                    
                    switch method
                        case 'fixedSize'
                            outref=smap.extendj(smap.cutj(outref,[os(1),os(2)]),[os(1),os(2)],mean(outref(:)));
                        otherwise
                            outref=outref;
                    end;
                    
                    outref=outref.*(sf.^2);
                    if( strcmp(sfType,'small')==1 )
                        tempCC=ccf(inref,outref); % % XX
                        peakVal=max(tempCC(:));
                        if( peakVal<100 )
                            disp('centering...');
                        end;
                        ctr=1;
                        while( peakVal<100 )
                            fs=os(1)-ctr;
                            tempRs=resize_2d(inref_F,fs,'down');
                            tempRs=tempRs.*(sf).^2;
                            tempRs=smap.extendj(tempRs,[os(1),os(2)],tempRs(1,1));
                            tempCC=ccf(tempRs,outref); % % XX
                            peakVal=max(tempCC(:));
                            ctr=ctr+1;
                        end;
                        ind=find(tempCC==peakVal,1,'first');
                        [x,y]=ind2sub(size(inref),ind);
                        cp=smap.getcp(inref);
                        nxny=maxInterpF(tempCC,16,20);
                        if( max(abs(nxny)>0 ) )
                            outref=smap.applyPhaseShifts(outref,-fliplr(nxny)); % % XX
                        end;
                    end;
                    
                    
                else
                    
                    os=[size(inref,1) size(inref,2) size(inref,3)];
                    finalSize=min([size(inref,1) size(inref,2) size(inref,3)]);
                    %     dsSize=finalSize/2;
                    dsSize=finalSize.*sf;
                    while( mod(dsSize,1)>1e-3 )
                        finalSize=finalSize+1;
                        dsSize=finalSize.*sf;
                    end;
                    finalSize=round(finalSize);
                    disp(['will pad original volume to ' num2str(finalSize) ' voxels...']);
                    %         disp(['key to continue...']);
                    %         pause;
                    inref=single(smap.extendj(inref,[finalSize,finalSize,finalSize],inref(1,1,1)));
                    
                    ns=[size(inref,1).*sf size(inref,2).*sf size(inref,3).*sf];
                    
                    inrefF=smap.ftj(inref);
                    inrefF_cut=single(smap.cutj(inrefF,[ns(2),ns(1),ns(3)]));
                    outref=smap.iftj(inrefF_cut);
                    
                    outref=outref.*(sf.^3);
                    
                end;
                
            end;
            
            function outref=resize_2d(inref_F,finalSize,resampleDir);
                switch resampleDir
                    case 'up'
                        inref_F_padded=smap.extendj(inref_F,[finalSize,finalSize],0);
                        outref=smap.iftj(inref_F_padded);
                    case 'down'
                        inref_F_cut=double(smap.cutj(inref_F,[finalSize,finalSize]));
                        outref=smap.iftj(inref_F_cut);
                end;
            end;
            
            function outref=getcp(inref);
                outref=[floor(size(inref,1)./2)+1 floor(size(inref,2)./2)+1];
            end;
            
        end;
        
        %% applyPhaseShifts:
        function outref=applyPhaseShifts(inref,shifts,varargin);
            % shifts is a MxN vector with offsets in units of pixels
            gpuFlag=0;
            if( isa(inref,'gpuArray')==1 )
                inref=gather(inref); gpuFlag=1; %wait(gdev);
            end;
            shifts_int=fix(shifts);
            shifts_subpix=shifts-fix(shifts);
            
            if( nargin<3 )
                cropFlag=0;
            else
                cropFlag=varargin{1};
            end;
            
            if( (size(shifts,1).*size(shifts,2))==3 )
                volFlag=1;
                disp('applying 3D phase shifts to volume...');
            else
                volFlag=0;
            end;
            
            if( volFlag==1 )
                movr=size(inref,1); movc=size(inref,2); movp=size(inref,3);
                xs=size(inref,1); ys=size(inref,2); zs=size(inref,3);
                
                inref=circshift(inref,[shifts_int(1) shifts_int(2) shifts_int(3)]);
                if( max(abs(shifts_subpix))>0 )
                    [xd,yd,zd]=meshgrid(((0:xs-1)-xs/2)/xs,((0:ys-1)-ys/2)/ys,((0:zs-1)-zs/2)/zs);
                    di_f=smap.ftj(inref);
                    clear inref;
                    dphs=yd.*(-shifts_subpix(1))+xd.*(-shifts_subpix(2))+zd.*(-shifts_subpix(3));
                    clear xd yd zd;
                    dphs=exp(1i.*2.*pi.*dphs);
                    d_done=di_f.*dphs;
                    mrDone=smap.iftj(d_done);
                else
                    mrDone=inref;
                end;
                outref=mrDone;
            else
                
                movNP=inref;
                nFrames=size(movNP,3);
                
                xs=size(inref,1); ys=size(inref,2);
                [xd,yd]=meshgrid(((0:xs-1)-xs/2)/xs,((0:ys-1)-ys/2)/ys);
                movShifted=zeros(xs,ys,nFrames);
                
                for i=1:nFrames
                    movNP(:,:,i)=circshift(movNP(:,:,i),[shifts_int(i,2) shifts_int(i,1)]);
                    if( max(abs(shifts_subpix(i,:)))>0 )
                        di_f=smap.ftj(movNP(:,:,i));
                        dphs=xd'.*(-shifts_subpix(i,2))+yd'.*(-shifts_subpix(i,1));
                        dphs=exp(1i.*2.*pi.*dphs);
                        d_done=di_f.*dphs;
                        movShifted(:,:,i)=smap.iftj(d_done);
                    else
                        movShifted(:,:,i)=movNP(:,:,i);
                    end;
                    if( cropFlag==1 )
                        dummy=circshift(movNP_dummy,ceil([shifts(i,2),shifts(i,1)]));
                        movShifted(:,:,i)=movShifted(:,:,i).*dummy;
                    end;
                end;
                
                if( cropFlag==1 )
                    movSum=sum(movShifted,3);
                    badCols=find((nansum(movSum,1))==0);
                    badRows=find((nansum(movSum,2))==0);
                    goodRows=setdiff(1:size(movSum,1),badRows);
                    goodCols=setdiff(1:size(movSum,2),badCols);
                    movShifted=movShifted(goodRows,goodCols,:);
                end;
                
                outref=movShifted;
            end;
            
            if( gpuFlag==1 )
                outref=gpuArray(outref); %wait(gdev);
            end;
            
        end;
        
        %% cropPatchFromImage3:
        function outref=cropPatchFromImage3(inref,halfDim,rowColInds);
            centerPixel(1)=floor((size(inref,1)/2)+1);
            centerPixel(2)=floor((size(inref,2)/2)+1);
            shifts=fliplr(centerPixel)-fliplr(rowColInds);
            movref_shifted=smap.applyPhaseShifts(inref,shifts);
            outref=smap.cutj(movref_shifted,[halfDim.*2 halfDim.*2]);
        end;
        
        %% maxInterpF:
        function [nxny,peakVal,errorFlags]=maxInterpF(imref,varargin);
        % [nxny,peakVal,errorFlags]=maxInterpF(imref,varargin);
        % varargin can include: [halfWidth],[centerPixel]

            rtu=min(size(imref));
            iFactor=4;
            if( nargin>1 )
                rtu=varargin{1};
                if( nargin>2 )
                    iFactor=varargin{2};
                end;
            end;
            cp_init=floor((size(imref,1)/2)+1);
            cc=smap.cropPatchFromImage3(imref,rtu,[cp_init cp_init]);
            cp=floor((size(cc,1)/2)+1);
            [nx,ny]=find(cc==max(cc(:)),1,'first');
            shifts_fullpix=-([cp cp]-[nx ny]);
            ccPatch=smap.cropPatchFromImage3(cc,rtu,[nx ny]);
            ccPatch_i=smap.resize_F(ccPatch,iFactor,'newSize');
            cp_i=floor((size(ccPatch_i,1)/2)+1);
            nanmask=smap.rrj(ccPatch_i);
            nanmask(find(nanmask>0.5))=nan;
            [nx_i,ny_i]=find(ccPatch_i==nanmax(ccPatch_i(:)),1,'first');
            shifts_subpix=-([cp_i cp_i]-[nx_i ny_i])./iFactor;
            shifts=shifts_fullpix+shifts_subpix;
            shiftsToUse=-fliplr(shifts);
            nxny=-shifts;
            peakVal=max(ccPatch_i(:));
        end;
        
        %% getcp:
        function outref=getcp(inref);
            nDims=length(size(inref));
            outref=floor(size(inref)./2)+1;
%             outref=[floor(size(inref,1)./2)+1 floor(size(inref,2)./2)+1];
        end;
        
        %% rrj:
        function outref=rrj(inref);
            % inref is any matrix that is the size of the k-space
            
            dims=[size(inref,1) size(inref,2) size(inref,3)];
            Npix=size(inref,1);
            cp=smap.getcp(inref);
            cp=cp(1);
            x=[1:Npix]-cp;
            
            X=zeros(dims(1),dims(2),dims(3));
            for i=1:dims(3)
                X(:,:,i)=repmat(x,dims(1),1);
            end;
            
            if( dims(3)>1 )
                Y=permute(X,[2 1 3]);
                Z=permute(X,[3 1 2]);
                R=sqrt(X.^2+Y.^2+Z.^2);
                outref=R./(2.*R(1,cp,cp));%Npix;
            else
                Y=permute(X,[2 1]);
                R=sqrt(X.^2+Y.^2);
                outref=R./(2.*R(cp,1));
            end;
            
        end;
        
        %% searchForPDB:
        function outref=searchForPDB(varargin);
            outref=[];
            % open a web browser to rcsb to search for PDBs
            url='http://www.rcsb.org/pdb/home/home.do';
            loc=[];
            [stat,h,url]=web(url,'-new');
            while( length(loc)==0 )
                pause(0.1);
                loc=get(h,'CurrentLocation');
            end;
            while(isempty(loc)==0)
                locLast=loc;
                loc=get(h,'CurrentLocation');
            end;
            url=locLast;
            disp(url);
            assignin('base','url',url);
            fileParts=regexp(url,'structureId=','split');
            if( length(fileParts)==2 )
                newPDB=fileParts{2};
                disp(['new PDB is ' newPDB]);
                
                %
                downloadDir='~/Documents/MATLAB/';
                structureDir='~/matching/temp/';
                %
                
                % % make a getComputerProps method
                cd(downloadDir);
                z=dir([downloadDir '*.gz']);
                for j=1:length(z)
                    tempName=char(regexp(z(j).name,newPDB,'match','ignorecase'));
                    if( length(tempName)>0 )
                        zipName=z(j).name;
%                         zipName=[tempName '.gz'];
                        disp(['unzipping ' downloadDir zipName '...']);
                        words=['gzip -d ' downloadDir zipName];
                        system(words);
                    end;
                end;
                
                z=dir([downloadDir '*.pdb*']); ctr=1;
                for j=1:length(z)
                    tempName=char(regexp(z(j).name,newPDB,'match','ignorecase'));
                    
                    if( length(tempName)>0 )
                        if(ctr>1)
                            disp(['Enter to move ' tempName ' into structure dir ' structureDir upper(tempName) '/ ...']);
                        else
                            disp(['Moving ' tempName ' into structure dir ' structureDir upper(tempName) '/ ...']);
                        end;
                        if( exist([structureDir upper(tempName) '/'],'dir')==7 )
                            
                        else
                            disp(['making new structure dir for ' upper(tempName) '...']);
                            mkdir([structureDir upper(tempName) '/']);
                        end;
                        newStructure=[structureDir upper(tempName) '/' upper(tempName) '.pdb'];
                        movefile([downloadDir z(j).name],newStructure);
                        outref=newStructure;
                        ctr=ctr+1;
                    end;
                end;
                
%                 disp('no PDB downloaded');
                
            else
                disp('last page was not a pdb download page');
            end;
        end;
    
        %% calculate the electrostatic potential from a new PDB file:
        function varargout=PDB2EP(structureName,nmPerVoxel);
            projName='trash.mrc';

            structureDir=['~/matching/temp/'];
            
            
            pdbLine=['pdb_file_in = ' structureName '.pdb'];
            voxelLine=['voxel_size = ' num2str(nmPerVoxel)];
            
            df=70/1e3;
            logLine=['log_file = ' structureName '_' datestr(now,30) '.log'];
            dfLine=strcat('defocus_nominal = ',[' ' sprintf('%4.3f',df)]);
            imFileLine=['image_file_out = ' projName ];
            
            
            of={
                '=== simulation === '
                'generate_micrographs = yes'
                logLine
                ''
                '=== sample ==='
                'diameter = 300 # in nm; diameter of hole in carbon'
                'thickness_center = 0'
                'thickness_edge = 0'
                ''
                '=== particle arb ==='
                'source = pdb'
                pdbLine
                'add_hydrogen = no'
                voxelLine
                'map_file_re_out = arb_real_map.mrc'
                'map_file_im_out = arb_im_map.mrc'
                ''
                '=== particleset ==='
                'particle_type = arb'
                'particle_coords = file'
                'coord_file_in = arb_location.txt'
                ''
                '=== geometry ==='
                'gen_tilt_data = no'
                'tilt_mode = single_particle'
                'geom_file_in = arb_rotations.txt'
                'ntilts = 1'
                'theta_start = 0 '
                'theta_incr = 0 '
                'geom_errors = none'
                ''
                '=== electronbeam ==='
                'acc_voltage = 300 # in kV'
                'energy_spread = 0.0 # in eV'
                'gen_dose = yes'
                'dose_per_im = 1000 # e/nm^2'
                ''
                '=== optics ==='
                'magnification = 103626 # set to give final pixel size of 0.965A/pix (= 5e4 A/pix / (0.965A/pix))'
                'cs = 2.7 # in mm'
                'cc = 0 # in mm'
                'aperture = 300 # in mm; aperture in focal plane (back focal plane?)'
                'focal_length = 3.5 # in mm; focal length of ?primary lens? (objective, presumably)'
                %         'cond_ap_angle = 0.05 # in millirad (aperture angle of condenser lens)'
                'cond_ap_angle = 0.00 # in millirad (aperture angle of condenser lens)'
                'gen_defocus = yes'
                dfLine
                ''
                '=== detector ==='
                'det_pix_x = 512'
                'det_pix_y = 512'
                'padding = 150'
                'pixel_size = 5'
                'gain = 1'
                'use_quantization = no'
                'dqe = 1'
                'mtf_a = 0.5'
                'mtf_b = 0.5'
                'mtf_c = 0 '
                'mtf_alpha = 0'
                'mtf_beta = 0'
                imFileLine
                };
            
            
            
            fnOut=[structureDir structureName '/PDB2EP.txt'];
            cd([structureDir structureName]);
            fid=fopen(fnOut,'w');
            for i=1:length(of)
                fprintf(fid,'%s\n',of{i});
            end;
            fclose(fid);
            
            fnOut='arb_location.txt';
            fid=fopen(fnOut,'w');
            fprintf(fid,'1  6\n');
            fprintf(fid,'#\t x\t y\t z\t phi\t theta\t psi\n');
            fprintf(fid,'\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\n',zeros(6,1));
            fclose(fid);
            
            fnOut='arb_rotations.txt';
            fid=fopen(fnOut,'w');
            fprintf(fid,'1  3\n');
            fprintf(fid,'#\t phi\t theta\t psi\n');
            fprintf(fid,'\t %7.4f\t %7.4f\t %7.4f\n',zeros(6,1));
            fclose(fid);

            disp(['Using TEM-simulator to calculate EP for ' structureName '...']);
            
            words=['TEM-simulator PDB2EP.txt'];
            system(words);

            disp('done');
                        
        end;
            
        
        %% EP to scattering potential:
        function outref=EP2SP(structureName,edgeSize,varargin);
            % % modify this to use imaginary part too
            
%             structureName='2R7Pb0';
%             pixelSize=0.5;

            if( nargin>2 )
                pixelSize=varargin{1};
            end;
            
            structureDir=['~/matching/temp/' structureName];
            fname=[structureDir '/arb_real_map.mrc'];
            
            [pv,aPerPix_pv]=mr(fname);
            fprintf('Read in %s with %7.4f nm voxel size\n',fname,aPerPix_pv);
            
            nPix_vol=edgeSize;
            aPerPix_final=pixelSize;
            
            dx=aPerPix_final./1e10;
            
            cc=smap.def_consts();
                        
            % convert into a cube:
            mpv=min(pv(:));
            pvMinDim=min(size(pv));
            pvMaxDim=max(size(pv));
            if( pvMinDim<pvMaxDim )
                disp('making a cube...');
                pv_cube=single(smap.extendj(pv,[pvMaxDim,pvMaxDim,pvMaxDim],mpv));
            else
                pv_cube=pv;
            end;
            clear pv;
            
            % subtract the background (water) potential:
            bgVal=pv_cube(1,1,1);
            pv_cube=pv_cube-bgVal;
            
            % 120615:
            disp('resizing...');
            SPV=smap.resize_F(pv_cube,aPerPix_pv/aPerPix_final);
            
            disp('padding...');
            SPV_out=single(smap.extendj(SPV,[edgeSize,edgeSize,edgeSize],0));
            
            disp('computing phase shifts...');
            SPV_out_2=SPV_out.*cc.IC.*dx./(2.*cc.k);
            
            % %
            % % this seems like a mistake because it reverses the handedness (052416):
            % %
            % % and as of 080215 we are reversing the slice order, as before Rullgard:
            temp=zeros(size(SPV_out_2));
            for i=1:size(SPV_out_2,3)
                temp(:,:,i)=SPV_out_2(:,:,size(SPV_out_2,3)+1-i);
            end;
            SPV_out_3=temp;
            
            SPV_out_3=smap.applyPhaseShifts(SPV_out_3,[0 0 1]);
            
            % % write to files for matching:
            % %
            % % move this into a separate method
            % % 
            SPV=single(SPV_out_3);
            if( exist('spv_ref','var')==1 )
                disp('registering with reference volume...');
                spv_reg=reg2vols(SPV,spv_ref,200);
                for i=1:3
                    imsc(squeeze(sum(spv_ref,i))); pause; imsc(squeeze(sum(spv_reg,i))); pause;
                end;
                SPV=spv_reg;
            end;
            
            SPVName=[structureDir '/SPV.mrc'];
            mw(single(SPV),SPVName,aPerPix_final);
            disp('done');
            
            
        end;
            
            
        %% definte constants (things that are true about the world):
        function consts=def_consts();
            consts=[];
            consts.V=300e3;
            consts.m_e=1.587.*(9.109e-31); % at 300 kV
            % consts.m_e=9.10938188e-31;
            consts.h=6.636e-34;
            consts.hbar=(6.636e-34)./(2.*pi);
            consts.q_e=1.602e-19;
            consts.wl=0.00197e-9; % wavelength at 300 kV
            consts.IC=(2.*consts.m_e./(consts.hbar.^2)).*consts.q_e;
            consts.k=2.*pi./consts.wl;
            consts.Cs=2.7e-3;
            consts.Cc=2.7e-3;
            consts.a_i=0.05e-3; % Illumination aperture [rad]
            consts.dE=0.7;
            consts.c = 2.99792458e8;
        end;
        
        

    end; % methods (Static)
    
    
end %classdef
