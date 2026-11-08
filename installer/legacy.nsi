; Separate native-Tk legacy beta; never replaces the modern application.
Unicode True
!include "MUI2.nsh"
!include "x64.nsh"
!include "LogicLib.nsh"
!include "FileFunc.nsh"
!include "WinVer.nsh"
!include "StrFunc.nsh"
${Using:StrFunc} StrRep
${Using:StrFunc} StrStr
Var NoIntegration
!ifndef PAYLOAD
  !error "Build with build_legacy.py --setup."
!endif
!ifndef EXE_SHA256
  !error "The verified executable SHA256 is required."
!endif
Name "In-Silico Docking Studio Legacy Beta"
OutFile "${OUTPUT}"
InstallDir "$LOCALAPPDATA\Programs\InSilicoDockingStudioLegacy"
InstallDirRegKey HKCU "Software\InSilicoDockingStudioLegacy" "InstallDir"
RequestExecutionLevel user
SetCompressor /SOLID lzma
SetCompressorDictSize 32
BrandingText "Legacy beta · Medical education"
VIProductVersion "1.3.0.1"
VIAddVersionKey "ProductName" "In-Silico Docking Studio Legacy Beta"
VIAddVersionKey "FileDescription" "Legacy beta for Windows 7 SP1, 8 and 8.1 x64"
VIAddVersionKey "ProductVersion" "1.3.0-legacy.1"
VIAddVersionKey "FileVersion" "1.3.0.1"
VIAddVersionKey "LegalCopyright" "Dr Manu Kumar Shetty and Dr Pratibha Nadig, 2026"
!define MUI_ICON "assets\icon.ico"
!define MUI_UNICON "assets\icon.ico"
!define MUI_ABORTWARNING
!define MUI_WELCOMEPAGE_TITLE "Docking Studio — Legacy Beta"
!define MUI_WELCOMEPAGE_TEXT "Separate native desktop edition intended for Windows 7 SP1, Windows 8 and Windows 8.1, 64-bit.$\r$\n$\r$\nIncludes Python and the scientific tools. No WebView2, browser or Anaconda installation is required for docking.$\r$\n$\r$\nThis is an experimental legacy build. It has not been executed on those older Windows systems; test it on your exact teaching computer before use.$\r$\n$\r$\nThis installation is separate from the Windows 10/11 application and keeps projects in a separate user-data folder."
!define MUI_FINISHPAGE_TITLE "Legacy Beta installed"
!define MUI_FINISHPAGE_TEXT "Open the Legacy Beta shortcut. Follow Targets, Drugs, Prepare, Dock and Results.$\r$\n$\r$\nRead LEGACY_WINDOWS.md for prerequisites, verification limits and installation help. Your saved projects are retained on uninstall."
!define MUI_FINISHPAGE_RUN "$INSTDIR\InSilicoDockingStudioLegacy.exe"
!define MUI_FINISHPAGE_RUN_TEXT "Open Docking Studio Legacy Beta"
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\LICENSE"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Function .onInit
  ${IfNot} ${RunningX64}
    MessageBox MB_ICONSTOP "This scientific legacy build needs 64-bit Windows. A 32-bit Windows edition is not included."
    Abort
  ${EndIf}
  ${IfNot} ${AtLeastWin7}
    MessageBox MB_ICONSTOP "This legacy build needs Windows 7 Service Pack 1 or newer. Windows XP and Vista are not supported."
    Abort
  ${EndIf}
  ${If} ${IsWin7}
    ${IfNot} ${AtLeastServicePack} 1
      MessageBox MB_ICONSTOP "Install Windows 7 Service Pack 1 before using this legacy build."
      Abort
    ${EndIf}
  ${EndIf}
  ; CPython 3.8 requires the supported DLL-directory API on Windows 7.
  ; Checking presence is read-only; setup never installs OS updates.
  System::Call 'kernel32::GetModuleHandleW(w "kernel32.dll")p.r0'
  System::Call 'kernel32::GetProcAddress(p r0,m "AddDllDirectory")p.r1'
  ${If} $1 == 0
    MessageBox MB_ICONSTOP "Windows needs the AddDllDirectory update (KB2533623 or a later update containing it). Install the Microsoft update, restart Windows, then run this setup again."
    Abort
  ${EndIf}
  SetShellVarContext current
  StrCpy $NoIntegration "0"
  ${GetParameters} $0
  ClearErrors
  ${GetOptions} $0 "/NOINTEGRATION" $1
  ${IfNot} ${Errors}
    StrCpy $NoIntegration "1"
  ${EndIf}
FunctionEnd

Section "Legacy Beta" MainSection
  IfFileExists "$INSTDIR\InSilicoDockingStudioLegacy.exe" 0 upgrade_unlocked
  ClearErrors
  FileOpen $0 "$INSTDIR\InSilicoDockingStudioLegacy.exe" a
  IfErrors upgrade_locked
  FileClose $0
  Goto upgrade_unlocked
  upgrade_locked:
    IfSilent upgrade_locked_exit
    MessageBox MB_OK|MB_ICONEXCLAMATION "Close the running Legacy Beta app, then run setup again. No application files have been changed."
  upgrade_locked_exit:
    SetErrorLevel 10
    Quit
  upgrade_unlocked:
  SetOutPath "$INSTDIR"
  ClearErrors
  File /r "${PAYLOAD}\*"
  IfErrors copy_failed
  nsExec::ExecToStack '"$SYSDIR\certutil.exe" -hashfile "$INSTDIR\InSilicoDockingStudioLegacy.exe" SHA256'
  Pop $0
  Pop $1
  StrCmp $0 "0" 0 verify_failed
  ; Certutil on older Windows can separate digest bytes with spaces.
  ; Normalize horizontal whitespace and CR, preserving line boundaries.
  ${StrRep} $2 "$1" "$\r" ""
  ${StrRep} $2 "$2" " " ""
  ${StrRep} $2 "$2" "$\t" ""
  StrCpy $2 "$\n$2$\n"
  ; Match one complete 64-hex digest line, never a heading/path substring.
  ${StrStr} $3 "$2" "$\n${EXE_SHA256}$\n"
  StrCmp $3 "" verify_failed
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  WriteINIStr "$INSTDIR\installation.ini" "Setup" "NoIntegration" "$NoIntegration"
  ${If} $NoIntegration == "1"
    Goto installed
  ${EndIf}
  CreateDirectory "$SMPROGRAMS\In-Silico Docking Studio Legacy"
  CreateShortcut "$SMPROGRAMS\In-Silico Docking Studio Legacy\Legacy Beta.lnk" "$INSTDIR\InSilicoDockingStudioLegacy.exe"
  CreateShortcut "$DESKTOP\In-Silico Docking Studio Legacy Beta.lnk" "$INSTDIR\InSilicoDockingStudioLegacy.exe"
  CreateShortcut "$SMPROGRAMS\In-Silico Docking Studio Legacy\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  WriteRegStr HKCU "Software\InSilicoDockingStudioLegacy" "InstallDir" "$INSTDIR"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy" "DisplayName" "In-Silico Docking Studio Legacy Beta"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy" "DisplayVersion" "1.3.0-legacy.1"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy" "Publisher" "Dr Manu Kumar Shetty and Dr Pratibha Nadig"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy" "UninstallString" '$\"$INSTDIR\Uninstall.exe$\"'
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy" "NoModify" 1
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy" "NoRepair" 1
  installed:
    Goto finished
  copy_failed:
    IfSilent copy_failed_exit
    MessageBox MB_OK|MB_ICONSTOP "Some Legacy Beta application files could not be replaced. Close the app and run setup again."
  copy_failed_exit:
    SetErrorLevel 11
    Quit
  verify_failed:
    IfSilent verify_failed_exit
    MessageBox MB_OK|MB_ICONSTOP "The copied Legacy Beta executable could not be verified. Run setup again before opening it."
  verify_failed_exit:
    SetErrorLevel 12
    Quit
  finished:
SectionEnd

Section "Uninstall"
  SetShellVarContext current
  ReadINIStr $NoIntegration "$INSTDIR\installation.ini" "Setup" "NoIntegration"
  ${If} $NoIntegration != "1"
    Delete "$DESKTOP\In-Silico Docking Studio Legacy Beta.lnk"
    Delete "$SMPROGRAMS\In-Silico Docking Studio Legacy\Legacy Beta.lnk"
    Delete "$SMPROGRAMS\In-Silico Docking Studio Legacy\Uninstall.lnk"
    RMDir "$SMPROGRAMS\In-Silico Docking Studio Legacy"
    DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoDockingStudioLegacy"
    DeleteRegKey HKCU "Software\InSilicoDockingStudioLegacy"
  ${EndIf}
  !include "${UNINSTALL_FILES}"
  Delete "$INSTDIR\installation.ini"
  Delete "$INSTDIR\Uninstall.exe"
  RMDir "$INSTDIR"
SectionEnd
