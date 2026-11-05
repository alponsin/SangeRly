# Sangerly offline (Windows) - how to build and use

## What you get
A folder `dist/Sangerly/` containing portable R + all packages + the app + a launcher.
Double-click `Sangerly.vbs`: the app opens in its own window (Edge/Chrome app mode) and works with NO internet.
Only these need internet: BLAST, "Fetch" of a reference by NCBI accession, and "Validate taxon names".
Everything else (assembly, viewer, editing, SNPs, translation, export) is offline.
Tip: save the consensus FASTA in the field and run BLAST later when you have signal (upload it in the BLAST tab).

## Build it (once, on a Windows PC with internet where the app already runs)
1. Put the `desktop/` folder next to `app.R` (so you have app.R, R/, www/, desktop/).
2. In R:  `setwd("D:/alponsin/sangerly"); source("desktop/build_portable.R")`
   It takes a few minutes and ends with "OK: the portable R loads every package."
3. Add ONE line to `app.R`, as the first line inside `server <- function(input, output, session) {`:
   `if (exists("sangerly_register_session")) sangerly_register_session(session)`
   (It does nothing on shinyapps.io; in the desktop version it stops R when you close the window.)
   Do this BEFORE step 2, or re-run step 2 afterwards.
4. Test: turn the internet off and double-click `dist/Sangerly/Sangerly.vbs`.
   If nothing opens, run `Sangerly_debug.bat` (shows the error) or read `%LOCALAPPDATA%\Sangerly\last_error.txt`.

## Make a real Setup.exe (optional)
Install Inno Setup 6.3+ (free), open `dist/installer.iss`, press Compile -> `dist/Output/Sangerly_Setup.exe`.
It installs per user (no administrator rights), adds a Start-menu entry and an optional desktop shortcut.
To give your app an icon, put `Sangerly.ico` in `desktop/` before building.

## Notes
- Size: about 0.7-1 GB as a folder; the Setup.exe is about 250-400 MB. A USB stick is fine.
- 64-bit Windows 10/11 only. Needs Edge or Chrome (Edge is built into Windows) for the app window; otherwise the default browser is used.
- Windows may show a SmartScreen warning for an unsigned Setup.exe: "More info" -> "Run anyway".
- Pressing F5 in the app window is fine; closing the window stops the app after about 5 seconds.
- If the first start offline is slow (a few seconds), it is the theme trying to reach a web font; it still works.
- Updating the app later: re-run step 2 and copy the new folder / installer.
