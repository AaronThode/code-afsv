classdef MatrixND
    %%%Methods%%%%
    %  out=extract_slice(obj,fixed_label,fixed_val);  %fixed val can be a
    %   scalar or two-element vector (for extracing a range of slices)
    %  out=sum_slice(obj,sum_label,bins_sum); %Sums over one dimension
    %  obj=append_time(obj,obj_new);
    %  obj=append_time_with_gaps(obj,obj_new)
    %  obj=plus(obj,obj_new)
    %  obj=histogram(obj,sum_label)
    %  obj=transpose_2D(obj)
    %  hh=plot1Ddistribution(obj,varargin)
    %  out=get_Tabs_bounds(obj)
    %  out=get_grid(obj,name)
    %  out=accumulate_slice(obj,sum_label,Nbins_to_combine)
    %  matrixx=image_2D_slice(obj,varargin)
    
    properties
        N
        bin_grid
        labels
    end
    
    properties ( Hidden=true)
        permitted_labels
        plot_labels
    end
    
    methods
        function td = MatrixND(N,bin_grid,labels,permitted_labels,plot_labels)
            %function td = MatrixND(N,bin_grid,labels,permitted_labels,plot_labels)
            if nargin > 0
                Ndim=length(bin_grid);
                
                if isempty(N)
                    dimms=[];
                    for Id=1:Ndim
                        dimms=[dimms length(bin_grid{Id})];
                    end
                    td.N=zeros(dimms,'single');
                else
                    keep_going=true;
                    
                    if Ndim==1 &size(N,2)>1 %flip row vector
                        N=N.';
                    end
                    
                    for Id=1:Ndim
                        keep_going=keep_going&(size(N,Id)==length(bin_grid{Id}));
                    end
                    
                    if keep_going&strcmpi(class(N),'single')
                        td.N=N;
                    else
                        disp('MatrixND: bin_grid and N not compatible, or N not ''single'' type');
                        td=[];
                    end
                    
                end
                
                if length(bin_grid)==Ndim
                    td.bin_grid=bin_grid;
                else
                    disp('Constructor MatrixND: mismatch in bin_grid dimension');
                    td=[];
                end
                
                if ~exist('permitted_labels','var')
                    td.permitted_labels={'Azimuth','Elevation','ItoE', 'KEtoPE','IntensityPhase','PSD','Frequency','Time'};
                else
                    td.permitted_labels=permitted_labels;
                end
                
                
                
                if ~exist('plot_labels','var')
                    td.plot_labels={'Azimuth (deg)','Elevation (deg)','Transport velocity', ...
                        'Kinetic/Potential ratio','Intensity Phase','PSD (dB re 1uPa^2/Hz)','Frequency (Hz)','Time'};
                else
                    td.plot_labels=plot_labels;
                end
                pass_me=all(contains(labels,td.permitted_labels,'IgnoreCase',true));
                if pass_me&length(labels)==Ndim
                    td.labels=labels;
                else
                    disp('Constructor MatrixND: mismatch in labels dimension');
                    td=[];
                end
                
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%extract_slice.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out=extract_slice(obj,fixed_label,fixed_val)
            % out=extract_slice(fixed_label,fixed_val)
            %    fixed_label: cell array of strings of dimension you want to fix
            %                   Can also be a single string, if only one
            %                   fixed value
            %    fixed_val: fixed value of fixed_dimension to take the slice.
            %           if two values, return slices with values between
            %           bounds
            
            out=[];
            if ~iscell(fixed_label)&ischar(fixed_label)
               % disp('converting fixed_label char to cell')
                temp{1}=fixed_label;
                fixed_label=temp;
            end
            
            if ~iscell(fixed_val)&isnumeric(fixed_val)
                %disp('converting fixed_val char to cell')
                temp{1}=fixed_val;
                fixed_val=temp;
            end
            pass_me=all(contains(fixed_label,obj.permitted_labels,'IgnoreCase',true));
            if ~pass_me
                disp('Labels not permitted');
                return
            end
            
            Iindex=zeros(size(fixed_label));
            for I=1:length(fixed_label)
                Iindex(I)=find(contains(obj.labels,fixed_label{I},'IgnoreCase',true));
                if isempty(Iindex(I))
                    disp('Requested label does not exist!');
                    return
                end
            end
            
            my_grid=obj.bin_grid(Iindex);
            
            %%%Indicies of dimensions that will remain
            Iindex_out=find(~contains(obj.labels,fixed_label,'IgnoreCase',true));
            
             %%% Keep fixed_label variable in output if fixed_val is a range
            
            for I=1:length(fixed_val)
                if length(fixed_val{I})>1
                    Iindex_out=sort([Iindex_out Iindex(I)]);
                end
            end
            
            %%%Remove variable from output if only a single bin
            for I=1:length(Iindex_out)
                if length(obj.bin_grid{Iindex_out(I)})==1
                    Iindex_out(I)=-1;
                end
            end
            Iindex_out(Iindex_out<0)=[];
            
            out_label=obj.labels(Iindex_out);
            out_grid=obj.bin_grid(Iindex_out);
            
            % Islice=zeros(size(fixed_val));
            for I=1:length(fixed_label)
                if length(fixed_val{I})==1
                    [~,Islice{I}]=min(abs(my_grid{I}-fixed_val{I}));
                    %fprintf('extract_slice:Requested Value: %6.2f, returned value: %6.2f\n',fixed_val{I},my_grid{I}(Islice{I}));
                elseif length(fixed_val{I})==2
                    Islice{I}=find((my_grid{I}>=fixed_val{I}(1))&(my_grid{I}<fixed_val{I}(2)));
                    out_grid{Iindex(I)}=out_grid{Iindex(I)}(Islice{I});  %restrict grid as well
                else
                    Islice{I}=fixed_val{I};
                    out_grid{Iindex(I)}=out_grid{Iindex(I)}(Islice{I});  %restrict grid as well
                end
            end
            
            %Iindex is the dimension, Islice is the element in that dimension
            %subsasgn, subsref
            %https://www.mathworks.com/matlabcentral/answers/344423-index-slice-of-nd-array-of-unknown-dimension
            %  Above lists how to remove a dimension of an array
            
            S.type = '()';
            S.subs = repmat({':'},1,ndims(obj.N));
            
            for I=1:length(Iindex)
                S.subs{Iindex(I)} = Islice{I}; % Specifiy index to extract
                
            end
            out = MatrixND(squeeze(subsref(obj.N,S)),out_grid,out_label,obj.permitted_labels,obj.plot_labels); %Subsref is generic indexing reference command
            
        end
        
        %%%%%%%%%%%%%%%%%%sum_slice.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out=sum_slice(obj,sum_label,bins_sum)
            % out=sum_slice(sum_label,bins_sum)
            %  sum_label: string of dimension you want to sum over
            % bins_sum: optional argument to sum over only bins_sum bins
            %           
            
            out=[];
            
            if nargin>2
               if bins_sum<1
                   error('sum_slice: bins_sum must be greater than zero');
               end
            end
            if ~iscell(sum_label)&ischar(sum_label)
               % disp('converting char to cell')
                temp{1}=sum_label;
                sum_label=temp;
            end
            
            pass_me=all(contains(sum_label,obj.permitted_labels,'IgnoreCase',true));
            if ~pass_me
                disp('Requested label not permitted');
                return
            end
            Iindex=find(contains(obj.labels,sum_label,'IgnoreCase',true));  %Highest dimension listed first
            Iindex_out=find(~contains(obj.labels,sum_label,'IgnoreCase',true));
             
            
            if nargin>2  %resampling
                
                   
                for I=1:length(Iindex)
                    sizz=size(obj.N);
               
                    obj.N=cumsum(obj.N,Iindex(I));
                    
                    S.type = '()';
                    S.subs = repmat({':'},1,ndims(obj.N));
                    
                    ngrid_new=sizz(Iindex(I));
                    S.subs{Iindex(I)} = bins_sum:bins_sum:ngrid_new; % Specifiy index to extract
                    
                    test2=subsref(obj.N,S);
                    
                    %%%%create layer of zeros and prepend
                    sizz(Iindex(I))=1;
                    obj.N=diff(cat(Iindex(I),zeros(sizz),test2),1,Iindex(I));
                    %  obj.N=diff([0 test(bins_sum:bins_sum:end)]); 1-D equivalent
                    obj.bin_grid{Iindex(I)}=obj.bin_grid{Iindex(I)}(bins_sum:bins_sum:ngrid_new);
                    if length(bins_sum:bins_sum:ngrid_new)==1
                        out_label=obj.labels(Iindex_out);
                        out_grid=obj.bin_grid(Iindex_out);
                    else
                        out_label=obj.labels;
                        out_grid=obj.bin_grid;
                    end
                end
                
               
            else  %% sum and remove dimension
                out_label=obj.labels(Iindex_out);
                out_grid=obj.bin_grid(Iindex_out);
                
                for I=1:length(Iindex)
                    obj.N=sum(obj.N,Iindex(I));  %Iindex does not need to be sorted if squeeze is performed once.
                end
            end
            
        
        
        out = MatrixND(squeeze(obj.N),out_grid,out_label,obj.permitted_labels,obj.plot_labels);
        
        
    end %sum_slice
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%append_time.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=append_time(obj,obj_new)
            %obj=append_time(obj_new)
            Index1=(contains(obj.labels,'Time','IgnoreCase',true));
            Index2=(contains(obj_new.labels,'Time','IgnoreCase',true));
            
            if ~any(Index1&Index2)
                disp('Two matricies don''t use same dimension for time; abort');
                return
            end
            Index1=find(Index1);
            obj.N=cat(Index1,obj.N,obj_new.N);
            obj.bin_grid{Index1}=[obj.bin_grid{Index1} obj_new.bin_grid{Index1}];
        end
        
        
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%append_time_with_gaps.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        function obj=append_time_with_gaps(obj,obj_new)
            %obj=append_time(obj_new)
            Index1=(contains(obj.labels,'Time','IgnoreCase',true));
            Index2=(contains(obj_new.labels,'Time','IgnoreCase',true));
            
            if ~any(Index1&Index2)
                disp('Two matricies don''t use same dimension for time; abort');
                return
            end
            Index1=find(Index1);  %%%Index number of bin_grid with 'Time'
            tabs=obj.bin_grid{Index1};
            dt=tabs(2)-tabs(1);
            
            tabs_new=obj_new.bin_grid{Index1};
            if dt~=tabs_new(2)-tabs_new(1)
                 disp('Time interval of new object not the same; abort');
                return
            end
            
            tabs_middle=tabs(end):dt:tabs_new(1);
            Nother=setdiff(1:ndims(obj.N),Index1);
            
            if length(Nother)>1
                
                disp('append_time_with_gaps requires 2D matrix; abort');
                return
            end
            
            if Nother==1
                obj_mid.N=zeros(length(obj.bin_grid{Nother}),length(tabs_middle));
            else
                obj_mid.N=zeros(length(tabs_middle),length(obj.bin_grid{Nother}));
            end
            
            obj.N=cat(Index1,obj.N,obj_mid.N);
            obj.N=cat(Index1,obj.N,obj_new.N);
            
            obj.bin_grid{Index1}=[obj.bin_grid{Index1} tabs_middle];
            
            obj.bin_grid{Index1}=[obj.bin_grid{Index1} obj_new.bin_grid{Index1}];
        end
        
        %%%%%%%%%%%%%%%%%%plus.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=plus(obj,obj_new)
            if ndims(obj.N)==ndims(obj_new.N)
                obj.N=obj.N+obj_new.N;
            end
        end
        
        %%%%%%%%%%%%%%%%%%histogram.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=histogram(obj,sum_label)
            %histogram(sum_label);
            other_label=~contains(obj.labels,sum_label,'IgnoreCase',true);
            obj=sum_slice(obj,obj.labels(other_label));
        end
        
        %%%%%%%%%%%%%%%%%%transpose_2D.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=transpose_2D(obj)
            %transpose_2D transposes a 2D matrix for alternate plotting
            %choices
            if length(obj.labels)~=2
                disp('MatrixND: transpose_2D: must be 2D matrix');
                return
            end
            obj.N=obj.N';
            obj.bin_grid=obj.bin_grid([2 1]);
            obj.labels=obj.labels([2 1]);
            
        end
        
        %%%%%%%%%%%%%%%%%%plot1Ddistribution.m%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function hh=plot1Ddistribution(obj,varargin)
            % function plot1Ddistribution(obj,plot_chc)
       
            %%%Check that I am 2D
            if ~isvector(obj.N)
                disp('MatrixND: bar: must be 1D matrix');
                return
            end
            
            if nargin==1
                plot_chc='bar';
            else
                plot_chc=varargin{1};
            end
            
             dx=obj.bin_grid{1}(2)-obj.bin_grid{1}(1);
            
           
            if nargin<3  %no normalization
                NN=obj.N;
                ylab='Count';
            elseif contains(varargin{2},'pdf')
                NN=obj.N/(dx*sum(obj.N));
               
                ylab='Probability Density';
            elseif contains(varargin{2},'cdf')
                NN=cumsum(obj.N)/sum(obj.N);
                ylab='Cumulative Probability';
                
            else
                NN=obj.N/max(obj.N);
                ylab='Relative Probability';
            end
                
            
            if strcmpi(plot_chc,'bar')
               hh = bar(obj.bin_grid{1},NN);grid on
            elseif strcmpi(plot_chc,'plot')
               hh = plot(obj.bin_grid{1},NN,'linewidth',3);grid on
            end
            xlabel(obj.labels{1});ylabel(ylab);
            set(gca,'fontweight','bold','fontsize',14);
        end
        
        %%%%%%%%%%%%%%%%%%get_Tabs_bounds%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out=get_Tabs_bounds(obj)
            Ilab=find(contains(obj.labels,'Time','IgnoreCase',true));
            out=[min(obj.bin_grid{Ilab}) max(obj.bin_grid{Ilab})];
        end
        
          %%%%%%%%%%%%%%%%%%get_grid%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out=get_grid(obj,name)
            Ilab=find(contains(obj.labels,name,'IgnoreCase',true));
            out=obj.bin_grid{Ilab};
        end
        
        %%%%%%%%%%%accumulate_slice%%%%%%%%%%%%%%%%%%%%%%
        %  consolidate into fewer bins 
        function obj=accumulate_slice(obj,sum_label,Nbins_to_combine)
            %out=accumulate_slice(obj,sum_label,Nbins_to_combine)
            
            if ~iscell(sum_label)&ischar(sum_label)
                disp('accumulate_slice:converting char to cell')
                temp{1}=sum_label;
                sum_label=temp;
            end
            
            pass_me=all(contains(sum_label,obj.permitted_labels,'IgnoreCase',true));
            if ~pass_me
                disp('accumulate_slice: Labels not permitted');
                return
            end
            
            %%%Check that I am 2D
            if length(obj.labels)>2
                disp('MatrixND: accumulate_slice: must be 1 or 2D matrix');
                return
            end
            
            Iindex=find(contains(obj.labels,sum_label,'IgnoreCase',true));  %Highest dimension listed first
            Iother=find(~contains(obj.labels,sum_label,'IgnoreCase',true));  %for loop
           
            Nbin=length(obj.bin_grid{Iindex});
            bin_edges=obj.bin_grid{Iindex};
            dX=bin_edges(2)-bin_edges(1);
            bin_edges=[bin_edges-0.5*dX bin_edges(end)+0.5*dX];
            
            
            edges_A=bin_edges(unique([1:Nbins_to_combine:length(bin_edges) length(bin_edges)]));
            edges_A_mid=edges_A(1:(end-1))+0.5*dX*Nbins_to_combine;
            obj.bin_grid{Iindex}=edges_A_mid;
            
            Nbins_out=floor(Nbin./Nbins_to_combine);
            subs=ones(Nbins_to_combine,1)*(1:Nbins_out);
            subs=subs(:);
            subs=[subs; (Nbins_out+1)*ones(Nbin-length(subs),1)];
            Nbins_out=Nbins_out+1;
            
            if isempty(Iother) %1-D arrayu
                A=accumarray(subs(:),obj.N);
                obj.N=A;
            else
                Nother=length(obj.bin_grid{Iother});
                for I=1:Nother
                    
                    S.type = '()';
                    S.subs = repmat({':'},1,ndims(obj.N));
                    S.subs{Iother} = I; % Specifiy index to extract
                    
                    A=accumarray(subs(:),  subsref(obj.N,S));
                    
                    S.subs{Iindex}=1:length(A);
                    obj.N=subsasgn(obj.N,S,A);
                end
                %%%Trim array
                S.type = '()';
                S.subs = repmat({':'},1,ndims(obj.N));
                S.subs{Iindex}=1:length(A);
                obj.N= subsref(obj.N,S);
            
            end
        
            
        end
        
        %%%%%%%%%%%%%%%%%%image_2D_slice%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function matrixx=image_2D_slice(obj,varargin)
            %image_2D_slice(scale_str,is_log,plot_chc,azimuth_from,contours)
            % scale_str: 'raw', joint, conditional, maxnorm, 'conditionalcumulative'
            %           raw: plot raw counts
            %           joint; plot p(X and Y)
            %           conditional; plot(Y|X)  (matrix/sum(matrix));
            %           MaxnormPerX: plot(X and Y) divided by maximum value of
            %                   Y along X
            %           MaxnormPerY:  plot(X and Y) divided by maxicumm
            %           value of X along Y
            %           conditional cumulative:  plot cumulative
            %           distribution of p(Y|X);
            % is_log, if true plot on log scale
            %plot_chc: 'contourf','image'
            %aziuth_from:  if 'true' (the default) show azimuth from which sound is
            %       coming from, not to
            %contours: vector of contour intervals
            
            is_log=false;
            plot_chc='image';
            azimuth_from=true;
            contourss=[];
            if nargin==1
                scale_str='raw';     
            elseif nargin==2
                scale_str=varargin{1};
            elseif nargin==3
                scale_str=varargin{1};
                is_log=varargin{2};
            elseif nargin==4
                scale_str=varargin{1};
                is_log=varargin{2};
                plot_chc=varargin{3};
            elseif nargin==5
                scale_str=varargin{1};
                is_log=varargin{2};
                plot_chc=varargin{3};
                azimuth_from=varargin{4};
            elseif nargin==6
                scale_str=varargin{1};
                is_log=varargin{2};
                plot_chc=varargin{3};
                azimuth_from=varargin{4};
                contourss=varargin{5};
            end
            
            %%%Check that I am 2D
            if length(obj.labels)~=2
                disp('MatrixND: image_2D_slice: must be 2D matrix');
                return
            end
            
            index1=contains(obj.permitted_labels,obj.labels{1},'IgnoreCase',true);
            plot_label{1}=obj.plot_labels{index1};
            index2=contains(obj.permitted_labels,obj.labels{2},'IgnoreCase',true);
            plot_label{2}=obj.plot_labels{index2};
            
            matrixx=double(obj.N);
            Ntotal=double(sum(sum(obj.N)));
            
            xplot=obj.bin_grid{2};
            yplot=obj.bin_grid{1};
            if strcmpi(obj.labels{2},'Frequency') && max(xplot>10000)
                xplot=xplot/1000;
                plot_label{2}='Frequency (kHz)';
            elseif strcmpi(obj.labels{1},'Frequency') && max(yplot>10000)
                yplot=yplot/1000;
                plot_label{1}='Frequency (kHz)';
            end
            %%%Change azimuth to directon of propagation toward...
            Itick=(strcmpi(obj.labels,'Azimuth'));
            if any(Itick)
                if ~azimuth_from
                    
                    %if strcmpi(obj.labels{1},'Azimuth')
                    [obj.bin_grid{Itick},Isort]=sort(bnorm(180+obj.bin_grid{Itick}));
                    if Itick(1)
                        matrixx=matrixx(Isort,:);
                    else
                        matrixx=matrixx(:,Isort);
                    end
                    plot_label{Itick}='Azimuth toward (deg)';
                else
                    plot_label{Itick}='Azimuth from (deg)';
                end
                
            end
            %if ~any(contains(obj.labels,'Frequency'))&~any(contains(obj.labels,'Time'))
            %   matrixx=matrixx./Ntotal;
            %end
            
            tit_str='raw';
            if strcmpi(scale_str,'joint')
                matrixx=matrixx./(Ntotal);
            elseif strcmpi(scale_str,'conditional') %% shows p(label(1)|label(2))
                px=histogram(obj,obj.labels(2));
                px_norm=(double(px.N.')./double(sum(px.N)));
                matrixx=matrixx./(Ntotal);
                
                matrixx=matrixx./px_norm;
                tit_str='conditional';
                disp('conditional plot');
            elseif strcmpi(scale_str,'conditionalcumulative')
                px=histogram(obj,obj.labels(2));
                px_norm=(double(px.N.')./sum(double(px.N)));
                matrixx=matrixx./(Ntotal);
                
                matrixx=cumsum(matrixx./px_norm);
                tit_str='conditionalcumulative';
                disp('conditionalcumulative plot');
                  
            elseif contains(scale_str,'MaxnormPerX','IgnoreCase',true)
                matrixx=matrixx./max(matrixx);
                tit_str='normalized';
                disp('normalized plot');
            elseif contains(scale_str,'MaxnormPerY','IgnoreCase',true)
                matrixx=matrixx';
                matrixx=(matrixx./max(matrixx))';
                tit_str='normalized';
                disp('normalized plot');
            end
            
            if is_log
                matrixx=log10(double(matrixx));
            end
            
            
            switch plot_chc
                case 'image'
                    
                    imagesc(xplot, yplot,matrixx);
                case 'pcolor'
                    
                    pcolor(xplot, yplot,matrixx);
                case 'contourf'
                    
                    if isempty(contourss)
                        %contourrs=10:2:30;
                        contourss=0:0.1:1;
                        
                    end
                    [C,H]=contourf(xplot,yplot, matrixx,contourss);axis ij
                    clabel(C,H)
                    
                    
            end
            xlabel(plot_label{2});ylabel(plot_label{1});
            title(sprintf('%s,%i samples',tit_str,Ntotal));
            strr='yx';
            
            Itick=(strcmpi(obj.labels,'Azimuth'));
            if any(Itick)
                set(gca,[strr(Itick) 'tick'],0:30:360);
            end
            
            Itick=(strcmpi(obj.labels,'PSD'));
            if any(Itick)
                set(gca,[strr(Itick) 'tick'],10*floor(min(obj.bin_grid{Itick})/10):10:10*ceil(max(obj.bin_grid{Itick})/10));
            end
            
            
            Itick=strcmpi(obj.labels,'Time');
            if any(Itick)
                datetick(strr(Itick),'HH:MM');
            end
            
            axis xy
            set(gca,'fontweight','bold','fontsize',14);
            colorbar('East','color','w')
            set(gca,'gridcolor','w');grid on
            
        end
        
      
    end
end