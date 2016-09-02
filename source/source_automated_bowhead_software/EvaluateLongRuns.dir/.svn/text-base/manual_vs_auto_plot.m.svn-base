load location_test

run_options.min_stations=2;  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria

run_options.Ncalls_to_sample=200;  %Number of snip samples to read into memory at once
run_options.force_energy_detector_run=1;  %If 1, always force JAVA cell CFAR to run
run_options.load_single_DASAR_results=0;  %If one, load stage 3 output (morph processing) from file
run_options.dt_slop=0.5;  %How much tolerance to give time delays between DASARs
run_options.bearing_alg='sel_ratio';
run_options.plot_locations=0;  %If one, plot locations..
run_options.localization_alg='Huber'; %HuberCalibratedKappa, repeated_median, MedianHuber
run_options.filter_chc='min_weight'; %dt_error, min_weight
param.morph.threshold_chc='otsu';

manual_kappa=manual.individual{1}.kappa(Icall_manual,:)';
auto_kappa=auto.locations{1}{Icall_auto}.kappa;

locations=compute_position(locations,goodFile,param,Icase,Isite,run_options);
                       