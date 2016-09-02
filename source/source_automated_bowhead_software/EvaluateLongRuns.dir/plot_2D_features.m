 function [Ichcx,Ichcy]=plot_2D_features(whale,other,feature_names)

        Ichcx=menu('Select x parameter to plot:',feature_names);
        Ichcy=menu('Select y parameter to plot:',feature_names);
        %Ichcz=menu('Select z parameter to plot:',feature_names)

        chc=[Ichcx Ichcy ];

        xlimm=[Inf -Inf];
        ylimm=[Inf -Inf];
        for K=1:length(whale),
            xlimm(1)=min([xlimm(1) whale{K}.(feature_names{Ichcx})]);
            xlimm(2)=max([xlimm(2) whale{K}.(feature_names{Ichcx})]);

            ylimm(1)=min([ylimm(1) whale{K}.(feature_names{Ichcy})]);
            ylimm(2)=max([ylimm(2) whale{K}.(feature_names{Ichcy})]);

            xlimm(1)=min([xlimm(1) other{K}.(feature_names{Ichcx})]);
            xlimm(2)=max([xlimm(2) other{K}.(feature_names{Ichcx})]);

            ylimm(1)=min([ylimm(1) other{K}.(feature_names{Ichcy})]);
            ylimm(2)=max([ylimm(2) other{K}.(feature_names{Ichcy})]);


        end
        Nx=linspace(xlimm(1),xlimm(2),100);
        Ny=linspace(ylimm(1),ylimm(2),91);


        for K=1:length(whale),

            [xgrid,ygrid,count]=contour_2D_histogram(whale{K}.(feature_names{Ichcx}),whale{K}.(feature_names{Ichcy}),Nx,Ny);
            if K==1,
                total_whale_count=count;
            else
                total_whale_count=total_whale_count+count;
            end

            [xgrid,ygrid,count]=contour_2D_histogram(other{K}.(feature_names{Ichcx}),other{K}.(feature_names{Ichcy}),Nx,Ny);
            if K==1,
                total_false_count=count;
            else
                total_false_count=total_false_count+count;
            end



        end

        subplot(2,1,1);
        imagesc(xgrid,ygrid,log10(total_whale_count));colorbar;
        set(gca,'fontweight','bold','fontsize',14);
        xlabel(feature_names{Ichcx});ylabel(feature_names{Ichcy});grid on;
        %adjust_feature_plot_limits(feature_names,Ichcx,Ichcy);
        subplot(2,1,2);
        imagesc(xgrid,ygrid,log10(total_false_count));colorbar;

        set(gca,'fontweight','bold','fontsize',14);
        xlabel(feature_names{Ichcx});ylabel(feature_names{Ichcy});grid on;