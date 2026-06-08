function results = run_csopf_dc(mpc, mpopt)
    %% SET OPTIONS
    if nargin < 2 || isempty(mpopt)
        mpopt = mpoption;
    end
    %% INTERNAL NUMBERING
    % Convert to internal indexing (buses, gens, branches)
    % Keep same pattern as MATPOWER's opf.m
    mpc_int = ext2int(mpc, mpopt);

    %% BUILD BASE-CASE DC OPF MODEL
    om = opf_setup(mpc_int, mpopt);   % standard MATPOWER OPF model object

    %% EXTEND WITH CS-SCOPF VARIABLES & CONSTRAINTS
    om = add_csopf_dc(om, mpc_int);

    %% SOLVE
    [results_int, success, raw] = opf_execute(om, mpopt);

    results_int.success = success;
    results_int.raw     = raw;

    %% BACK TO EXTERNAL NUMBERING
    results = int2ext(results_int);
end
