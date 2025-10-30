import os
import subprocess
from typing import Optional


def _volumes_list() -> set[str]:
    try:
        return set(os.listdir("/Volumes"))
    except Exception:
        return set()


def is_mounted(volume_name: str) -> bool:
    """
    Returns True if a volume with this name appears under /Volumes.
    """
    return volume_name in _volumes_list()


def mount_volume(smb_url: str, volume_name: Optional[str] = None) -> tuple[bool, str]:
    """
    Mount an SMB share via macOS Finder using AppleScript:
      smb_url examples:
        'smb://idnas37.d.uzh.ch/g_psyplafor_methlab$'
        'smb://idnas37.d.uzh.ch/g_psyplafor_methlab_data$'
    If Keychain has credentials, this is seamless (no prompt).
    Returns (success, message).
    """
    # If the caller knows the expected mount name, skip work if already mounted
    if volume_name and is_mounted(volume_name):
        return True, f"'{volume_name}' already mounted."

    cmd = ["osascript", "-e", f'mount volume "{smb_url}"']
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode == 0:
        if volume_name and is_mounted(volume_name):
            return True, f"Mounted '{volume_name}' successfully."
        return True, "Mounted successfully."
    else:
        # Common reasons: missing credentials in Keychain; network ACLs; typos.
        return False, proc.stderr.strip() or proc.stdout.strip()


def ensure_mounted(smb_url: str, volume_name: str) -> None:
    """
    Idempotent convenience: mount if not present; print a clear status.
    """
    if is_mounted(volume_name):
        print(f"[INFO] {volume_name} already mounted.")
        return
    ok, msg = mount_volume(smb_url, volume_name)
    status = "SUCCESS" if ok else "FAIL"
    print(f"[{status}] {volume_name}: {msg}")
