%%%%%%%%%plot_filtered_spectrograms%%%%%%%%%%%%%%%%%
function plot_feature_filtered_spectrogram_and_morph(inputt,filename,param,feature_names,Ichcx,Ichcy,coords,plot_label)
for KK=1:length(inputt)
    Ix=give_N_closest_fits(coords(1),inputt{KK}.(feature_names{Ichcx}),100);
    Iy=give_N_closest_fits(coords(2),inputt{KK}.(feature_names{Ichcy})(Ix),10);
    inputt{KK}.Igood=Ix(Iy);
    figure
    %disp('Warning: starting at index 6');
    for L=1:length(Iy),
        close all;

        subplot(3,1,1)
        tlen=inputt{KK}.duration_debug(inputt{KK}.Igood(L));
        ctime=inputt{KK}.ctime_debug(inputt{KK}.Igood(L));
        [y,t,head]=readfile(filename{KK},ctime,tlen,1,'ctime','calibrate');
        spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');

        subplot(3,1,2)
        imagesc(full(inputt{KK}.Image{inputt{KK}.Igood(L)}))
        title('stored image');
       

        yes=1;
        while (yes==1);
            %%Set processor you want
            %Ithresh=menu('Which threshold?','otsu','nonlinear peak');
            Ithresh=1;
            if Ithresh==1,
                param.morph.threshold_chc='local_peaks';
            else
                param.morph.threshold_chc='nonlinear';
            end
            
           %param.morph.equalization(:,1)=auto_station(Istation).param.equalization_freq;
            param.morph.equalization(:,1)=inputt{KK}.param.equalization_freq;
            param.morph.equalization(:,2)=inputt{KK}.equalization(:,inputt{KK}.Igood(L));

            [stats,labeled_image]=extract_image_features(y,ctime,param,0);

            subplot(3,1,3);
            imagesc(sum(labeled_image,3));
            stored_str=sprintf('Station %i,%s, %s call %i of %i, %s: %6.2f, %s: %6.2f',KK, ctime2str(inputt{KK}.ctime_debug(inputt{KK}.Igood(L))), ...
                plot_label, L,length(Iy), ...
                feature_names{Ichcx},inputt{KK}.(feature_names{Ichcx})(inputt{KK}.Igood(L)),...
                feature_names{Ichcy},inputt{KK}.(feature_names{Ichcy})(inputt{KK}.Igood(L)));
            title(stored_str);
            disp(stored_str);

            [feature_namex,index.x]=modify_values(feature_names{Ichcx});
            [feature_namey,index.y]=modify_values(feature_names{Ichcy});


            for II=1:length(stats),


                disp(sprintf('New morph value of %s is %6.2f, for %s is %6.2f',feature_namex,stats(II).(feature_namex)(index.x), ...
                    feature_namey,stats(II).(feature_namey)(index.y)));
            end
            yes=menu('Redo?','Yes','No');
        end
    end
end

    function [nameout,index_out]=modify_values(input_name)
        index_out=1;
        nameout=input_name;

        switch input_name
            case 'BoundingBox'
                nameout='duration';
            case 'Centroid'
                index_out=2;
                %nameout='Centroid Frequency';
        end

    end
end %function plot_feature_...


function Ix=give_N_closest_fits(xpick,xdata,N)
fitt=(abs(xpick-xdata));
[junk,Isort]=sort(fitt,'ascend');
Imax=min([length(Isort) N]);
Ix=Isort(1:Imax);
end