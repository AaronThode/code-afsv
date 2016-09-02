%%%%%%%%initialize_location_object.m%%%%%%%%
%function locations=initialize_location_object(station,Istation,Nstations,options),
%
% initialize_location_object,
% choice of basic fields or complete.
% output: locations{I}.dt (Nstations x1)
%                     .indicies (station index , Nstations x1)
%                     .ctime (Nstationsx1)
%                     .SNR (Nstationsx1);
function locations=initialize_location_object(station,Istation,Nstations,options)

Nobjects=length(station.ctime_min);
%options='basics';
locations=cell(1,Nobjects);


if strcmp(options,'basic'),
    for I=1:Nobjects,
        if rem(I,1000)==0,disp(sprintf('Initialization of %ith anchor station %i percent complete',Istation, round(100*I/Nobjects)));end
        locations{I}.dt=NaN*ones(Nstations,1);
        locations{I}.dt(Istation)=0;

        locations{I}.station_indicies=zeros(Nstations,1);
        locations{I}.station_indicies(Istation)=I;

%         for J=1:length(fnames)
%             locations{I}.(fnames{J})=zeros(Nstations,1);
%             locations{I}.(fnames{J})(Istation)=station.(fnames{J})(I);
%         end
         locations{I}.ctime_min=zeros(Nstations,1);
        locations{I}.ctime_min(Istation)=station.ctime_min(I);

        locations{I}.SEL=zeros(Nstations,1);
        locations{I}.SEL(Istation)=station.SEL(I);
    end
end

