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
distrib3d_input = "\n".join([
    "-99",
    "{RUN}/cem140w04floc_{ID}.img",
    "cement140",
    "{RUN}/cement140w04flocf_{ID}.img",  # fixed name
    "0.7344","0.6869","0.0938","0.1337","0.1311","0.1386","0.0407","0.0408"
])

# === Input for disrealnew ===
# NOTE: again {RUN}/... because disrealnew runs with cwd=_BIN_DIR
disrealnew_input = "\n".join([
    "-2794",
    "{RUN}/cement140w04flocf_{ID}.img",
    "1 2 3 4 5 6 7 28 26",
    "35",
    "{RUN}/pcem140w04floc_{ID}.img",
    "44990","1","5850","2","8692","3","2631","4","1100","5","2062","6","839","7",
    "0",
    "1000","0","500",
    "0.0001 9000.","0.01 9000.","0.00002 10000.","0.002 2500.",
    "50","5","5000","100",
    "0.00","20.0","20.0","0.0",
    "40.0","83.14","80.0",
    "0.00035","0.72",
    "0","0","1","10","1.0","0","0","1"
]) + "\n"   # ensure newline at end

# Run the pipeline
m.run_pipeline(id, gp_cfg, distrib3d_input, disrealnew_input)
