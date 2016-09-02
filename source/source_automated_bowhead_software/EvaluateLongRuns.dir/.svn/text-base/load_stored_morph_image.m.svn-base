%%%load_stored_morph_image.m%%%
function [stored_image,new_image,stats]=load_stored_morph_image(param,station,station_index,automated_dir,auto_template_name,filename,Idebug)

            stored_image=[];new_image=[];stats=[];
            
            if ~isempty(station_index), %morph result does exist

                tlen=param.energy.bufferTime+mean(station.duration_debug(station_index));
                start_ctime=min(station.ctime_debug(station_index));
                index_offset_t=auto_station.ctime(station_index)-start_ctime-bufferTime;

                %%Download original image

                fname=dir([automated_dir '/' auto_template_name '*morph*']);
                auto=load([ automated_dir '/' fname.name]);
                %end
                I_image=min(station.indicies(1,station_index));
                stored_image=full(auto.best_calls.features{I_image}(1).final_image);

               
            elseif ~exist('start_ctime');
                disp('adjust_morph: start_ctime does not exist...');
                keyboard
            else
                index_offset_t=param.energy.bufferTime;

            end

            
           
        end
        [y,t,head]=readfile(filename,start_ctime,tlen,1,'ctime','calibrate');
        
        param_out=param;

        change_flag=1;
        while change_flag==1,
            [stats,new_image]=extract_image_features(y,start_ctime,param_out,1);
            figure(5);
            if ~isempty(new_image),
                if ~isempty(stats),
                    FF=stats(1).dF*(0:(size(new_image,1)-1));
                    TT=stats(1).dT*(0:(size(new_image,2)-1));
                else
                    FF=(0:(size(new_image,1)-1));
                    TT=(0:(size(new_image,2)-1));

                end

                subplot(2,1,1)
                plot_morph_image(TT,FF,station_index,index_offset_t,new_image);
                title(sprintf('direct from sio, start time (including buffer): %s:%i', ...
                    datestr(datenum(1970,1,1,0,0,start_ctime),31),round(1000*(start_ctime-floor(start_ctime)))));

                if ~isempty(stored_image)
                    subplot(2,1,2)
                    plot_morph_image(TT,FF,station_index,index_offset_t,stored_image);
                    title(sprintf('Stored image at %ith index in station, %ith index in best_calls',station_index,I_image));

                else
                    title('Stored image does not exist');
                    subplot(2,1,2);
                    cla;
                end
                hold off
            else
                disp('No morph returned');

            end
            if exist('show_spec')&&show_spec>0,
                figure(6);
                spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
                caxis([-60 80])
                pause;
                close(6)
            end

            [tmp,change_flag]=alter_parameters(param_out.morph);
            param_out.morph=tmp;

        end