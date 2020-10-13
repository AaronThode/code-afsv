classdef MatrixND
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
                    td.permitted_labels={'Azimuth','ItoE', 'KEtoPE','IntensityPhase','PSD','Frequency','Time'};
                else
                    td.permitted_labels=permitted_labels;
                end
                
                
                
                if ~exist('plot_labels','var')
                    td.plot_labels={'Azimuth (deg)','Transport velocity', ...
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
                disp('converting fixed_label char to cell')
                temp{1}=fixed_label;
                fixed_label=temp;
            end
            
            if ~iscell(fixed_val)&isnumeric(fixed_val)
                disp('converting fixed_val char to cell')
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
                    fprintf('extract_slice:Requested Value: %6.2f, returned value: %6.2f\n',fixed_val{I},my_grid{I}(Islice{I}));
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
        function out=sum_slice(obj,sum_label)
            % out=sum_slice(sum_label)
            %  sum_label: string of dimension you want to sum over
            
            out=[];
            
            if ~iscell(sum_label)&ischar(sum_label)
                disp('converting char to cell')
                temp{1}=sum_label;
                sum_label=temp;
            end
            
            pass_me=all(contains(sum_label,obj.permitted_labels,'IgnoreCase',true));
            if ~pass_me
                disp('Labels not permitted');
                return
            end
            Iindex=find(contains(obj.labels,sum_label,'IgnoreCase',true));  %Highest dimension listed first
            
            Iindex_out=find(~contains(obj.labels,sum_label,'IgnoreCase',true));
            out_label=obj.labels(Iindex_out);
            out_grid=obj.bin_grid(Iindex_out);
            
            
            for I=1:length(Iindex)
                obj.N=sum(obj.N,Iindex(I));  %Iinded does not need to be sorted if squeeze is performed once.
            end
            out = MatrixND(squeeze(obj.N),out_grid,out_label,obj.permitted_labels,obj.plot_labels);
            
        end %sum_slice
        
        %%%%%%%%%%%%%%%%%%append_time.m%%%%%%%%%%%%%%%%%%%%%%
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
        function plot1Ddistribution(obj,varargin)
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
            
            if nargin<3  %no normalization
                NN=obj.N;
                ylab='Count';
            elseif contains(varargin{2},'pdf')
                NN=obj.N/sum(obj.N);
                ylab='Probability';
            else
                NN=obj.N/max(obj.N);
                ylab='Relative Probability';
            end
                
            
            if strcmpi(plot_chc,'bar')
                bar(obj.bin_grid{1},NN);grid on
            elseif strcmpi(plot_chc,'plot')
               plot(obj.bin_grid{1},NN);grid on
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
        
        %%%%%%%%%%%accumulate_slice%%%%%%%%%%%%%%%%%%%%%%
        %  consolidate into fewer bins 
        function obj=accumulate_slice(obj,sum_label,N_accum)
            %out=accumulate_slice(obj,sum_label,Nbins_to_sum)
            
            
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
            
            
            edges_A=bin_edges(unique([1:N_accum:length(bin_edges) length(bin_edges)]));
            edges_A_mid=edges_A(1:(end-1))+0.5*dX*N_accum;
            obj.bin_grid{Iindex}=edges_A_mid;
            
            Nbins_out=floor(Nbin./N_accum);
            subs=ones(N_accum,1)*(1:Nbins_out);
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
        function image_2D_slice(obj,varargin)
            %image_2D_slice(scale_str,is_log,plot_chc)
            % scale_str: 'raw',joint,conditional,norm, 'conditionalcumulative'
            
            plot_chc='image';
            if nargin==1
                scale_str='raw';
                is_log=false;
                
            elseif nargin==2
                is_log=false;
                scale_str=varargin{1};
                
            elseif nargin==3
                is_log=varargin{2};
                scale_str=varargin{1};
            elseif nargin>=4
                is_log=varargin{2};
                scale_str=varargin{1};
                plot_chc=varargin{3};
            else
                
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
                
                %%%Change azimuth to directon of propagation toward...
                Itick=(strcmpi(obj.labels,'Azimuth'));
                if any(Itick)
                    %if strcmpi(obj.labels{1},'Azimuth')
                    [obj.bin_grid{Itick},Isort]=sort(bnorm(180+obj.bin_grid{Itick}));
                    if Itick(1)
                        matrixx=matrixx(Isort,:);
                    else
                        matrixx=matrixx(:,Isort);
                    end
                    plot_label{Itick}='Azimuth toward (deg)';
                end
                matrixx=cumsum(matrixx./px_norm);
                tit_str='conditionalcumulative';
                disp('conditionalcumulative plot');
                
                
            elseif strfind(scale_str,'norm')
                matrixx=matrixx./max(matrixx);
                tit_str='normalized';
                disp('normalized plot');
            end
            
            if is_log
                matrixx=log10(double(matrixx));
            end
            
%             if strcmpi(scale_str,'conditionalcumulative')
%                 plot(obj.bin_grid{1},matrixx);grid on
%                 if strcmpi(obj.labels{1},'Azimuth')
%                     set(gca,'xtick',0:30:360);
%                     %set(gca,'xtick',obj.bin_grid{1}(sortt));
%                     xtickss=get(gca,'xticklabel');
%                     for I=1:length(xtickss)
%                        mytick{I}=int2str(bnorm(180+str2num(xtickss{I})));
%                     end
%                     set(gca,'xticklabel',mytick);
%                 end
%                 
%                 ylim([0 1]);
%                 for I=1:length(obj.bin_grid{2})
%                     legstr{I}=num2str(obj.bin_grid{2}(I));
%                 end
%                 legend(legstr)
%                 
%                 xlabel(plot_label{1});ylabel('Cumulative Fraction');
%                 axis xy
%                 set(gca,'fontweight','bold','fontsize',14);
%                 return
%             end
            
            switch plot_chc
                case 'image'
                    imagesc(obj.bin_grid{2}, obj.bin_grid{1},matrixx);
                case 'contourf'
                    
                    if nargin<5
                       %contourrs=10:2:30;
                        contourrs=0:0.1:1;
                    else
                        contourrs=varargin{end};
                    end
                    [C,H]=contourf(obj.bin_grid{2},obj.bin_grid{1}, matrixx,contourrs);axis ij
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
                datetick(strr(Itick),'HH');
            end
            
            axis xy
            set(gca,'fontweight','bold','fontsize',14);
            colorbar('East')
            
        end
        
      
    end
end