function plot_directional_metric(TT,FF,output_array,handles,param,PdB,use_wavelets)
%function plot_directional_metric(TT,FF,output_array,handles,param)
%  TT, FF, time (sec) and frequency (Hz) associated with output_array
%  handles: handles to figure
%  param:  climm, alg

climm=eval(param.climm);
if strcmpi(handles.display_view,'Elevation')& max(climm)==360
      climm=[-90 90];
end
alg_mult=eval(param.alg);
if ~exist('PdB','var')
    PdB=[];
end

if strcmpi(handles.display_view,'AdditiveBeamforming')
    
    hh=imagesc(TT,FF/1000,10*log10(abs(output_array)));
    format_spectrogram_image(true);
     titstr='Additive Beamforming';
elseif strcmpi(handles.display_view,'Azimuth')|strcmpi(handles.display_view,'Elevation')

    %%%convert to mask if desired
    if eval(param.mask)
        crange=eval(param.climm);
        temp=zeros(size(output_array));
        Igood=find(output_array>=crange(1)&output_array<=crange(2));
        temp(Igood)=1;
        output_array=temp;

    end
    if ~use_wavelets
        hh=imagesc(TT,FF/1000,output_array);
    else
        %         dbcont
        %         hh=pcolor(TT,log2(FF/1000),output_array);
        %         hh.EdgeColor='none';
        %         ax=gca;
        %         ytick=round(pow2(ax.YTick),3);
        %         ax.YTickLabel=ytick;
        
        hh=surface(TT,(FF/1000),output_array);
        hh.EdgeColor='none';
        ax=gca;
        ax.XLabel.String='Time';
        ax.YLabel.String='Frequency';
        
    end
    colorbar

    titstr=handles.display_view;
    try
        if get(handles.checkbox_grayscale,'Value')==1
            colormap(flipud(gray));
            clim(climm);
        elseif eval(param.mask)==1
            colormap(flipud(gray))
        else
            if strcmpi(handles.display_view,'Elevation')
                colormap(jet)
            else
                colormap(hsv);
            end
            clim(climm);
        end
    catch
        colormap(hsv);
    end
    
    
    
    %%%Use alpha adjustment to display transport ratio as well
    
    %          set(hh,'AlphaData',intensity./energy_density);
    %          set(hh,'AlphaDataMapping','scaled')
    %          alim([0 1]);
    %
    
    %%%Setting transparency to SPL
    if handles.checkbox_transparency.Value&~isempty(PdB)&&all(size(PdB)==size(output_array))
        set(hh,'AlphaDataMapping','scaled')
        set(hh,'AlphaData',PdB);
        %set(hh,'FaceAlpha','flat');
        alim(str2double(handles.edit_mindB.String) + [0 str2double(handles.edit_dBspread.String)]);
        
    end
    
else
    %%%Effective velocity j/Sp
    if strcmpi(handles.display_view,'ItoERatio')
        
        hh=imagesc(TT,FF/1000,output_array);caxis([0 1])
        hbar=colorbar;
        titstr='Intensity/Energy Ratio';
    elseif strcmpi(handles.display_view,'KEtoPERatio')
        
        %imagesc(TT,FF/1000,10*log10(output_array));
        imagesc(TT,FF/1000,(output_array));
        hbar=colorbar;
        titstr='Kinetic/Potential Ratio';
        %caxis([0 90])
        %caxis([-20 20])
        caxis([-10 10]);
        caxis([-3 3]);
        hbar.Ticks=[-10 -6 -3 0 3 6 10];
    elseif strcmpi(handles.display_view,'IntensityPhase')
        
        hh=imagesc(TT,FF/1000,output_array);
        hbar=colorbar;
        titstr='IntensityPhase';
        caxis([0 90])
        
        %%%Setting transparency to SPL
        if handles.checkbox_transparency.Value&~isempty(PdB)&&all(size(PdB)==size(output_array))
            set(hh,'AlphaDataMapping','scaled')
            set(hh,'AlphaData',PdB);
            %set(hh,'FaceAlpha','flat');
            alim(str2double(handles.edit_mindB.String) + [0 str2double(handles.edit_dBspread.String)]);
            
        end
    elseif strcmpi(lower(handles.display_view),lower('Polarization'))
        
        %choice=menu('Choice for Polarization focus:',param.polar_desc(:));
        %Icol=ceil(choice/3);Irow=choice-3*(Icol-1);
       
        Ndims=size(output_array,3);
        choice_stokes=menu('Focal plot for which Stokes parameter?','S0','S1','S2','S3','DegreePolarization');
        stokes_str={'S0','S1','S2','S3','DegreePolarization'};
        stokes_col=[1 2 3 3 4];
        Icol=stokes_col(choice_stokes);
       
        if Ndims==1
            Irow=1;
        else
            choice_dimension=menu('Focal plot for which polarization plane?','xy','yz','zx');
            dim_str={'xy','yz','zx'};
            Irow=find(contains(param.polar_desc(:,1),dim_str{choice_dimension}));
        end
        fprintf('Selected focal plot to use %s\n',param.polar_desc{Irow,Icol});

        if Icol<=2
            hh=imagesc(TT,FF/1000,10*log10(squeeze(output_array(:,:,Irow,Icol))));
        elseif Icol==4
            hh=imagesc(TT,FF/1000,(squeeze(output_array(:,:,Irow,4))));clim([0 1])
        elseif choice_stokes==3
            hh=imagesc(TT,FF/1000,real(squeeze(output_array(:,:,Irow,3))));clim([-1 1])
        elseif choice_stokes==4
            hh=imagesc(TT,FF/1000,imag(squeeze(output_array(:,:,Irow,3))));clim([-1 1])
        end
        
        hbar=colorbar;
        titstr=sprintf('Polarization: %s %s',dim_str{choice_dimension},stokes_str{choice_stokes});
        myfig=gcf;

        %%%Setting transparency to SPL
        if handles.checkbox_transparency.Value&~isempty(PdB)&&all(size(PdB)==size(output_array(:,:,1,1)))
            set(hh,'AlphaDataMapping','scaled')
            set(hh,'AlphaData',PdB);
            %set(hh,'FaceAlpha','flat');
            alim(str2double(handles.edit_mindB.String) + [0 str2double(handles.edit_dBspread.String)]);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%Plot full polarization values %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure;
        for JJ=1:Ndims
            h1(JJ)=subplot(Ndims,1,JJ);
            imagesc(TT,FF/1000,imag(squeeze(output_array(:,:,JJ,3))));clim([-1 1]);colorbar
            format_spectrogram_image(false)
            title(sprintf('S3_%s',dim_str{JJ}),'Interpreter','none');
           
           
        end
        linkaxes(h1,'xy')
        
        figure
        for JJ=1:Ndims
            h2(JJ)=subplot(Ndims,1,JJ);
            imagesc(TT,FF/1000,(squeeze(output_array(:,:,JJ,4))));clim([0 1]);colorbar
            format_spectrogram_image(false)
            title(sprintf('DegreePolarization_%s',dim_str{JJ}),'Interpreter','none');
            
        end
        linkaxes(h2,'xy')
        
        figure(myfig)

        keyboard
    elseif strcmpi(handles.display_view,'PhaseSpeed')
        
        %%%Phase speed directly
        hh=imagesc(TT,FF/1000,output_array);
        hbar=colorbar;
        titstr='Phase Speed';
        caxis([1000 3000])
        %
        %         hh=imagesc(TT,FF/1000,real(acosd(1500./output_array)));
        %         titstr='Elevation angle';
        %         caxis([0 90])
        
        %%%Setting transparency to SPL
        if handles.checkbox_transparency.Value&~isempty(PdB)&&all(size(PdB)==size(output_array))
            
            set(hh,'AlphaData',PdB);
            set(hh,'AlphaDataMapping','scaled')
            alim(str2double(handles.edit_mindB.String) + [0 str2double(handles.edit_dBspread.String)]);
            
        end
    end
    
    
    
    %%%Set colorscale
    if get(handles.checkbox_grayscale,'Value')==1
        colormap(flipud(gray));
    else
        colormap(jet);
    end
    
end

grid on
axis('xy')

try
    fmax=str2double(get(handles.edit_fmax,'String'));
    fmin=str2double(get(handles.edit_fmin,'String'));
    if fmax==0
        ylim([0 handles.Fs/2000]);
        set(handles.edit_fmax,'String',num2str(Fs/2000));
    else
        if ~use_wavelets
            ylim([fmin fmax]);
        end
    end
    %ylim([0 1]);axis('xy')
    
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
    
    set(handles.text_filename,'String',sprintf('%s, %s ... ', ...
        get(handles.text_filename,'String'),titstr));
catch
    disp('catch in plot_directional_metric');
end

    function format_spectrogram_image(adjust_caxis)
        grid on
        axis('xy')
        fmax=str2double(get(handles.edit_fmax,'String'));
        fmin=str2double(get(handles.edit_fmin,'String'));
        if fmax==0
            ylim([0 Fs/2000]);
            set(handles.edit_fmax,'String',num2str(Fs/2000));
        else
            ylim([fmin fmax]);
        end
        %ylim([0 1]);axis('xy')

        if adjust_caxis
            climm(1)=str2double(get(handles.edit_mindB,'String'));
            climm(2)=climm(1)+str2double(get(handles.edit_dBspread,'String'));
            %%If switching from correlogram, reset to suggested values
            if climm(1)==0&&climm(2)<1
                climm=[40 70];
                set(handles.edit_mindB,'String',num2str(climm(1)));
                set(handles.edit_dBspread,'String',num2str(climm(2)));
                climm(2)=sum(climm);
            end

            caxis(climm);
        end

        if get(handles.checkbox_grayscale,'Value')==1
            colormap(flipud(gray));
        else
            colormap(jet);
        end
        colorbar;
        % set(gcf,'pos',[30   322  1229   426])
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (sec)');ylabel('Frequency (kHz)');
        xlim([0 str2num(handles.edit_winlen.String)]);
        if ~strcmp(handles.display_view,'Spectrogram')
            title(get(handles.text_filename,'String'));
        end
    end  %%%format_spectrogram_image

end