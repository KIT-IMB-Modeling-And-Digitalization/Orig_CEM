from cement_sim import module as m

id = "01"

# === Input for genpartnew ===
# NOTE: use short filenames (no {RUN}) because genpartnew runs with cwd=run_dir
gp_cfg = {
    # Basic controls
    "seed": -3034,                     # Random number seed
    "place_menu": 2,                   # Menu choice to place particles
    "n_size_classes": 16,              # Number of size classes to place
    "dispersion_px": 0,                # Dispersion distance between particles in pixels

    # Calcium sulfate parameters
    "calcium_sulfate_vf": "0.0604",    # Total volume fraction
    "calcium_sulfate_split": ["0.515", "0.041"],  # [hemihydrate, anhydrite] fractions

    # Size classes: list of (count, radius, phase_id)
    "size_classes": [
        (1, 17, 1), (1, 15, 1), (1, 14, 1), (1, 13, 1),
        (2, 12, 1), (2, 11, 1), (4, 10, 1), (5,  9, 1),
        (8,  8, 1), (13, 7, 1), (21, 6, 1), (38, 5, 1),
        (73, 4, 1), (174, 3, 1), (450, 2, 1), (2674, 1, 1)
    ],

    # Flocculation & output menus
    "report_phase_counts_menu": 4,     # Report phase counts
    "flocculate_menu": 3,              # Flocculate particles
    "n_flocs": 1,                      # Number of flocs to create
    "output_menu": 8,                  # Output current microstructure to file

    # Output filenames (short names; genpartnew runs in the sandbox folder)
    # Tip: keep {ID} so you can vary id without editing strings
    "out_image":        "cem140w04floc_{ID}.img",
    "out_particle_ids": "pcem140w04floc_{ID}.img",

    "exit_menu": 1                     # Exit program
}


# === Input for distrib3d ===
# NOTE: use {RUN}/... because distrib3d runs with cwd=_BIN_DIR
d3_cfg = {
    "seed": -99,
    "in_name":  "cem140w04floc_{ID}.img",      # short name; no path
    "filters_root": "cement140",
    "out_name": "cement140w04flocf_{ID}.img",  # short name; no path
    # order: C3S vf, C3S sa, C2S vf, C2S sa, C3A vf, C3A sa, C4AF vf, C4AF sa
    "targets": ["0.7344","0.6869","0.0938","0.1337","0.1311","0.1386","0.0407","0.0408"],
}


# === Input for disrealnew ===
# NOTE: again {RUN}/... because disrealnew runs with cwd=_BIN_DIR
dr_cfg = {
    "seed": -2794,

    # Files (short names; module will prefix with {RUN}/ and format {ID})
    "phase_file": "cement140w04flocf_{ID}.img",
    "part_file":  "pcem140w04floc_{ID}.img",

    # Phase assignments for: C3S C2S C3A C4AF gypsum hemihydrate anhydrite flyash slag
    # (either list of 9 items or a single space-separated string)
    "phase_map": [1, 2, 3, 4, 5, 6, 7, 28, 26],

    # Phase ID for C3A in fly ash particles
    "c3a_fa": 35,

    # One-pixel particles to add: 7 (count, phase_id) pairs
    "one_px_pairs": [
        ("44990","1"), ("5850","2"), ("8692","3"),
        ("2631","4"),  ("1100","5"), ("2062","6"),
        ("839","7"),
    ],
    # The trailing lone count line
    "one_px_extra": "0",

    # Main controls
    "cycles": "1000",
    "sat_flag": "0",
    "max_diff": "500",

    # Nucleation parameters (4 lines in this order)
    "nuc_params": ["0.0001 9000.", "0.01 9000.", "0.00002 10000.", "0.002 2500."],

    # Frequencies: [pore_perc, solids_perc, stats_out, micro_out]
    "freqs": ["50","5","5000","100"],

    # Thermal: [induction_h, T_init_C, T_amb_C, U]
    "thermal": ["0.00","20.0","20.0","0.0"],

    # Activation energies: [cement, pozz, slag]
    "Ea": ["40.0","83.14","80.0"],

    "cycle_to_time": "0.00035",
    "agg_vf": "0.72",

    # Flags (8 lines): isothermal, forbid_CSH, allow_CH_agg, movie_slices,
    #                  diss_bias, pre_resat_cycles, random_CSH, pH_affects
    "flags": ["0","0","1","10","1.0","0","0","1"],
}


# Run the pipeline
m.run_pipeline(id, gp_cfg, d3_cfg, dr_cfg)
