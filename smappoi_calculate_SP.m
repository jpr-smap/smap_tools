function varargout=smappoi_calculate_SP(paramsFile,jobref);
% compile with:
% /misc/local/matlab-2019a/bin/mcc -mv -I ~/smappoi_min/v2.1/src/ -R -singleCompThread smappoi_calculate_SP.m

clearvars -global 

jobNum=str2num(jobref);

[~,this_server]=system('uname -n');

%
% set up GPUs:
% gtu=mod(jobNum-1,8)+1;
words=['nvidia-smi -L'];
[~,resp]=system(words);
temp=splitlines(resp);
inds_gpus=find(~cellfun(@isempty,temp));
temp=temp(inds_gpus)
n_gpus=length(inds_gpus);
gtu=mod(jobNum-1,n_gpus)+1
gpu_name=temp{gtu}


fprintf('getting gpu # %s...',num2str(gtu));
try
    gdev=gpuDevice(gtu);
catch
    fprintf('Failed to get gpu # %s\n',num2str(gtu));
    fidFail=fopen([pwd '/fail_' smap.zp(jobNum,4) '.txt'],'w');
    fprintf(fidFail,'%s\n',datestr(now,31));
    try
        fprintf(fidFail,'%s\n',gpu_name);
    end;
%     fclose(fidFail);
%     if( debug_flag )
%         gdev=gpuDevice(1);
%     else
        return;
%     end;
end;


if( ~isempty(paramsFile) )
    try
        params=smap.readParamsFile(paramsFile);
    catch
        fprintf('Problem reading parameter file %s...\n',char(paramsFile));
        return;
    end;
end;

[~,fn,ext]=fileparts(params.PDBFile);
fn_len=min([length(fn) 4]);
scratchDir=fullfile(params.outputDir,['scratch_' fn(1:fn_len)]);
if( exist(params.outputDir,'dir')~=7 )
    fprintf('Making new project directory... [%s]\n',params.outputDir);
    mkdir(params.outputDir);
end;
if( exist(scratchDir,'dir')~=7 )
    fprintf('Making new scratch directory... [%s]\n',scratchDir);
    mkdir(scratchDir);
end;

disp(datestr(now));
fnLog=[scratchDir '/output_' smap.zp(jobNum,4) '.txt'];
fprintf('Making new logfile... [%s]\n',fnLog);
fidLog=fopen(fnLog,'w');
fprintf(fidLog,'%s\n',datestr(now,31));

fprintf(fidLog,'job %i of %i...\n',jobNum,params.nCores);
disp(params);

aPerPix=params.aPerPix;
baseDir=[pwd '/'];
disp(datestr(now));
fileBase='EP_';
fNumPadded=smap.zp(num2str(jobNum),4);

searchDir=[scratchDir '/'];
outputDir=[params.outputDir '/'];

inc=1;
cc=smap.def_consts();
baseDir=[pwd '/'];
fn_pdb=params.PDBFile;

[~,fn,ext]=fileparts(fn_pdb);
switch ext
    case {'.pdb','.PDB','.pdb1'}
        smap.read_pdb_file(fn_pdb);
    case {'.cif'}
        params_cif.CIFFile=fn_pdb;
        if( ~isempty(params.chains) )
            params_cif.chains=regexp(params.chains,' ','split');
        else
            params_cif.chains=[];
        end
        smap.read_cif_file(params_cif);
    otherwise
        fprintf('File format for %s not recognized\n',char(fn_pdb));
        fprintf(fidLog,'File format for %s not recognized\n',char(fn_pdb));
        fclose(fidLog);
        return
end;

fn_pdb=[fn ext];

if( ~isempty(params.bArb) )
    if( params.bArb == -1 )
        params.bArb=[];
    else
        fprintf(fidLog,'using imposed b-factor %5.3f for all atoms...\n',params.bArb);
        bFactor=ones(1,length(bFactor)).*params.bArb;
    end;
end;

padVal_A=5.0;

% move COM to (0,0,0):
xyz_orig=xyz;
xyz=xyz_orig-repmat((nanmean(xyz_orig,2)),1,size(xyz_orig,2));

if( isfield(params,'MW_target') )
    if( ~isempty(params.MW_target) )
        MW_target=params.MW_target
        atom_rad=sqrt(sum(xyz.^2,1));
        [atom_rad_s,sI]=sort(atom_rad,'ascend');
        dummy=[1:size(atom_rad,2)].*13.824./1000;
        ind=find(dummy<MW_target,1,'last');
        inds=sI(1:ind);
        atom_rad_s(ind)
        atomList=atomList(inds);
        bFactor=bFactor(inds);
        chainIDs=chainIDs(inds);
        atomNums=atomNums(inds);
        xyz=xyz(:,inds);
        % %
    end;
end;
xyz=xyz-repmat((nanmean(xyz,2)),1,size(xyz,2));


mm=[min(xyz,[],2) max(xyz,[],2)];
edgeSize_A=(2.*max(abs(mm(:))))+padVal_A;
edgeSize=ceil(edgeSize_A./aPerPix);
edgeSize=max([edgeSize 211]);

temp=single(smap.rrj(zeros(edgeSize,edgeSize,edgeSize)));
D=temp.*(size(temp,1)-1); % real-space values here
cp=floor(size(D,1)./2)+1;

tempVec=D(cp,:,cp);
tempVec(1:cp)=-tempVec(1:cp);
[X0,Y0,Z0]=meshgrid(tempVec,tempVec,tempVec);
[X0,Y0,Z0]=deal(X0.*aPerPix,Y0.*aPerPix,Z0.*aPerPix);

X=X0;
Y=Y0;
Z=Z0;
D=sqrt(X.^2+Y.^2+Z.^2);

nAtoms=length(atomList);
itu=jobNum:params.nCores:nAtoms;
al=char;
atu=atomList(itu);
btu=bFactor(itu);
xyz=xyz(:,itu);
for i=1:length(atu)
    temp=char(atu{i});
    al(i)=temp(1);
end;
u=unique(al);
clear uI; clear V_s;

fprintf(fidLog,'calculating potential for %d atoms...\n',length(itu));

V=zeros(size(X0));

temp=unique(sort(abs(X0(:)),'ascend'));
dx=temp(2);
atomRad=5;

m_e=cc.m_e;
h=cc.h;
q_e=cc.q_e;
c_v=cc.c_v;
IC=cc.IC;

    %         if( ~isempty(params.units) )
    %             temp=mod(params.units-1,8)+1;
    %         else
    %             temp=mod([0:(params.nCores-1)],8)+1
    %         end;
    %
    %         gtu=temp(jobNum);
    %         fprintf(fidLog,'getting gpu # %s...',num2str(gtu));
    %         tic;
    %         try
    %             gdev=gpuDevice(gtu);
    %             fprintf(fidLog,'%f seconds\n',toc);
    %             %         reset(gdev)
    %         catch
    %             fprintf(fidLog,'Failed to get gpu # %s\n',num2str(gtu));
    %             fidFail=fopen([scratchDir '/fail_' fNumPadded '.txt'],'w');
    %             fprintf(fidFail,'%s\n',datestr(now,31));
    %             exit;
    %         end;
    
    % list vars for xfer:
    gpuVars={'X','Y','Z','X0','Y0','Z0','D','V','inds','atomRad', ...
        'a','b_init','b','ctr','m_e','h','q_e','c_v','IC'}
    
    for j=1:length(gpuVars)
        if( exist(gpuVars{j},'var')==0 )
            eval([gpuVars{j} '=[];']);
        end;
        eval([gpuVars{j} '=gpuArray(' gpuVars{j} ');']); wait(gdev);
    end;

update_interval=200;

for ctr=1:length(itu)
    X=X0-xyz(1,ctr);
    Y=Y0-xyz(2,ctr);
    Z=Z0-xyz(3,ctr);
    D=sqrt(X.^2+Y.^2+Z.^2);
    inds=find(D<=atomRad); % indices for which to calculate potential
    
    el=al(ctr);
    [a,b_init]=smap.parameterizeSF(el);
    b=b_init+btu(ctr)+16.*(aPerPix.^2);
    lead_term=((16.*pi.^(5/2).*(h./(2.*pi)).^2)./(m_e.*q_e)).*1e20; % does not use relativistic electron mass
    sum_term=gpuArray.zeros(size(inds,1),5);
    for i=1:5
        sum_term(:,i)=(a(i)./(b(i).^(3./2))).*exp(-((4.*pi.^2).*(D(inds).^2))./b(i));
    end;
    sum_term=sum(sum_term,2);
    V(inds)=V(inds)+lead_term.*sum_term;
    
    if( mod(ctr,update_interval)==0 )
        fprintf('%d/%d\n',ctr,length(itu));
        fprintf(fidLog,'%d/%d\n',ctr,length(itu));
    end;
end;

for i=1:length(gpuVars)
    eval([gpuVars{i} '=gather(' gpuVars{i} ');']); wait(gdev);
end;

if( ~isempty( params.modelName ) )
    fn=strtrim(char(params.modelName));
else
    [~,fn]=fileparts(params.PDBFile);
end;

fn_out=[searchDir fn '_EP_' fNumPadded '.mrc'];
fn_EP=[outputDir fn '_EP.mrc'];
fn_SP=strrep(fn_EP,'_EP.mrc','_SP.mrc');

smap.mw(single(V),fn_out,aPerPix);

fprintf(['done with potential calculation at ' datestr(now,31) '\n']);
interval=0.01;
nextFile=1;


if( jobNum==params.nCores )
    nFilesExpected=params.nCores;%*2;
    fn_pdb=params.PDBFile;
    
    while 1
        fprintf(fidLog,'Looking for %i %s_EP files...\n',nFilesExpected,fn);
        
        fnList={};
        numFound=[]; fileTypesFound={};
        A=dir([searchDir fn '_EP_*.mrc']);
        ctr=1;
        for i=1:length(A)
            tempNum=regexp(A(i).name,[fn '_EP_(\d{4,4})'],'tokens');
            if( length(tempNum)>0 )
                numFound(ctr)=str2num(char(tempNum{1}));
                fnList{ctr}=[searchDir A(i).name];
                fprintf(fidLog,'%s\n',fnList{ctr});
                ctr=ctr+1;
            end;
        end;
        if( ctr>nFilesExpected )
            break;
        end;
        pause(1);
        
    end;
    
    fprintf(fidLog,'reading files...\n');
    for i=1:length(fnList)
        try
            if( i==1 )
                temp=smap.mr(fnList{i});
                EPV=zeros(size(temp,1),size(temp,2),size(temp,3));
            end;
            temp=smap.mr(fnList{i});
            EPV=EPV+temp;
        catch
            fprintf(fidLog,'problem reading file %i\n',i)
        end;
    end;
    
    %pause;
    
    wPot=4.871; % Rullgard's calculation
    out=EPV+ones(size(EPV)).*wPot;
    EPV=out;
    
    fprintf(fidLog,'writing combined electrostatic potential...\n');
    smap.mw(single(EPV),fn_EP,params.aPerPix);
    fprintf(fidLog,'done combining electrostatic potentials at %s\n',datestr(now,31));
    
    cc=smap.def_consts();
    
    SPV=EPV;
    try
        pd=smap.particleDiameter(SPV);%,0.05)
        edgeSize=size(smap.resizeForFFT(zeros(round(2.5*pd),round(2.5*pd)),'crop'),1);
    catch
        disp('problem estimating particle diameter...');
        pd=size(SPV,1)
        edgeSize=2.5*pd;
    end;
    
    if( edgeSize>params.edge_max )
        edgeSize=params.edge_max;
        fprintf(fidLog,'using edge_max for padding (%i^3)...\n',params.edge_max);
    end;
    edgeSize=max([edgeSize max(size(SPV))]);
    fprintf(fidLog,'padding to %i^3 volume...\n',edgeSize);
    SPV_out=single(smap.cropOrPad(SPV,[edgeSize,edgeSize,edgeSize],wPot));
    
    fprintf(fidLog,'computing phase shifts...\n');
    dx=params.aPerPix./1e10;
    
    lambda = h/sqrt(q_e*params.V_acc*m_e*(q_e/m_e*params.V_acc/c_v^2 + 2 ));
    k=2.*pi./lambda;
    SPV_out_2=SPV_out.*IC.*dx./(2.*k);
    
    %         %%
    fprintf(fidLog,'centering scattering potential...\n');
    
    inref=SPV_out_2;
    
    bg_val=mode(inref(:));
    inref=inref-bg_val;
    edge_size=size(inref,1);
    cp=floor(edge_size./2)+1;
    dummy=[-(cp-1):(cp-(2-mod(edge_size,2)))];
    [X,Y,Z]=meshgrid(dummy,dummy,dummy);
    
    M=sum(inref(:));
    COM_x=dot(X(:),inref(:))/M;
    COM_y=dot(Y(:),inref(:))/M;
    COM_z=dot(Z(:),inref(:))/M;
    disp([COM_x COM_y COM_z]);
    
    outref=smap.applyPhaseShifts(inref,-[COM_y COM_x COM_z]);
    
    M=sum(outref(:));
    COM_x=dot(X(:),outref(:))/M;
    COM_y=dot(Y(:),outref(:))/M;
    COM_z=dot(Z(:),outref(:))/M;
    
    disp([COM_x COM_y COM_z])
    
    outref=outref+bg_val;
    
    SPV_out_2=outref;
    
    %         %%
    
    fprintf(fidLog,'writing scattering potential...\n');
    
    smap.mw(single(real(SPV_out_2)),fn_SP,params.aPerPix);
    
    fprintf(fidLog,'deleting scratch files...\n');
    for i=1:length(fnList)
        delete(fnList{i});
    end;
    
    
    fprintf(fidLog,'Done combining intermediate output at %s\n',datestr(now,31));
    fprintf(fidLog,'Consolidating log files...\n');
    
    fclose(fidLog);
    
    fnFinal=strrep(fn_SP,'_SP.mrc','.log');
    fidFinal=fopen(fnFinal,'w');
    
    tline={''};
    for i=1:jobNum
        fn=[searchDir 'output_' smap.zp(i,4) '.txt'];
        try
            fid=fopen(fn,'r');
            while 1
                temp=fgets(fid);
                disp(temp);
                if ~ischar(temp)
                    break;
                else
                    fprintf(fidFinal,temp);
                end;
            end
            fclose(fid);
            delete(fn);
        catch
            tline{end+1}=['could not open ' fn];
            if( exist('fid','var') )
                if( fid>0 )
                    fclose(fid);
                end;
            end;
        end;
    end;
    
    
    fprintf(fidFinal,'*****\nFinished at %s\n*****',datestr(now,31));
    fclose(fidFinal);
    
    try
        rmdir(scratchDir);
        if( exist(fullfile([queueDir '/done/']),'dir')<7 )
            mkdir(fullfile([queueDir '/done/']));
        end;
        movefile(paramsFile,fullfile([queueDir '/done/']));
        
    catch
        
    end
    
end;





%%

