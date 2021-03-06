function plot_directional_metric(TT,FF,output_array,handles,param,PdB,use_wavelets)
%function plot_directional_metric(TT,FF,output_array,handles,param)
%  TT, FF, time (sec) and frequency (Hz) associated with output_array
%  handles: handles to figure
%  param:  climm, alg

climm=eval(param.climm);
alg_mult=eval(param.alg);
if ~exist('PdB','var')
    PdB=[];
end
if strcmpi(handles.display_view,'Directionality')
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
    titstr=' Azimuth';
    try
        if get(handles.checkbox_grayscale,'Value')==1
            colormap(flipud(gray));
        else
            colormap(hsv);
        end
    catch
        colormap(hsv);
    end
    caxis(climm);
    
    
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
        
        imagesc(TT,FF/1000,output_array);caxis([0 1])
        titstr='Intensity/Energy Ratio';
    elseif strcmpi(handles.display_view,'KEtoPERatio')
        
        %imagesc(TT,FF/1000,10*log10(output_array));
        imagesc(TT,FF/1000,(output_array));
        
        titstr='Kinetic/Potential Ratio';
        %caxis([0 90])
        caxis([-20 20])
    elseif strcmpi(handles.display_view,'IntensityPhase')
        
        imagesc(TT,FF/1000,output_array);
        titstr='Kinetic/Potential Ratio';
        caxis([0 90])
    elseif strcmpi(handles.display_view,'PhaseSpeed')
        
        %%%Phase speed directly
        hh=imagesc(TT,FF/1000,output_array);
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
    
    colorbar
    
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
    
    
    colorbar;
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
    
    set(handles.text_filename,'String',sprintf('%s, %s ... ', ...
        get(handles.text_filename,'String'),titstr));
catch
    disp('catch in plot_directional_metric');
end