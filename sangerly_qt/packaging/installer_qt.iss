; Inno Setup 6.3+ script for the Qt version. Compile it after packaging/build_exe.py (folder dist_qt\SangeRly).
#define AppName "SangeRly"
#ifndef AppVer
#define AppVer "1.0.0"
#endif
#define Src "..\dist_qt\SangeRly"
; AppPublisher is shown in the installer and in the Windows "unknown publisher" prompt; use the name on the certificate.
#define AppPublisher "SangeRly"

[Setup]
AppId={{B2F5C3E0-6B1A-4E8E-9D7A-5A1C0E5E1A01}
AppName={#AppName}
AppVersion={#AppVer}
AppPublisher={#AppPublisher}
VersionInfoVersion={#AppVer}
DefaultDirName={localappdata}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=..\Output
OutputBaseFilename=SangeRly_Setup
SetupIconFile=Sangerly.ico
UninstallDisplayIcon={app}\SangeRly.exe
Compression=lzma2/max
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

; Code signing. Compile with   ISCC /DSign /Ssg="<command with $f>" installer_qt.iss   (packaging\sign_and_build.ps1 does this for you).
; Without /DSign the installer is built unsigned, exactly as before.
#ifdef Sign
SignTool=sg
SignedUninstaller=yes
#endif

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"

[Files]
Source: "{#Src}\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
Source: "{#Src}\app\Data_example\*"; DestDir: "{userdocs}\Sangerly\Data_example"; Flags: recursesubdirs createallsubdirs ignoreversion uninsneveruninstall skipifsourcedoesntexist

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\SangeRly.exe"; WorkingDir: "{app}"; IconFilename: "{app}\Sangerly.ico"
Name: "{autoprograms}\{#AppName} - example data"; Filename: "{userdocs}\Sangerly"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\SangeRly.exe"; WorkingDir: "{app}"; IconFilename: "{app}\Sangerly.ico"; Tasks: desktopicon

[Run]
Filename: "{app}\SangeRly.exe"; Description: "Start SangeRly"; Flags: postinstall nowait skipifsilent
