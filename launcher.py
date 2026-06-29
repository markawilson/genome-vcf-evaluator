"""
Desktop launcher for Genome VCF Evaluator.
Starts the Streamlit server and opens the browser automatically.
Used as the entry point for PyInstaller builds.
"""

import os
import sys
import threading
import time
import webbrowser

PORT = 8501


def _open_browser():
    """Wait for the server to start, then open the default browser."""
    time.sleep(3)
    webbrowser.open(f"http://localhost:{PORT}")


def main():
    # When running as a PyInstaller bundle, sys._MEIPASS points to the
    # temporary directory containing the unpacked files.
    if getattr(sys, "frozen", False):
        bundle_dir = sys._MEIPASS
    else:
        bundle_dir = os.path.dirname(os.path.abspath(__file__))

    app_path = os.path.join(bundle_dir, "app.py")

    # Launch browser in background
    threading.Thread(target=_open_browser, daemon=True).start()

    # Start Streamlit programmatically
    from streamlit.web import cli as stcli

    sys.argv = [
        "streamlit", "run", app_path,
        f"--server.port={PORT}",
        "--server.headless=true",
        "--server.fileWatcherType=none",
        "--browser.gatherUsageStats=false",
        "--global.developmentMode=false",
    ]
    stcli.main()


if __name__ == "__main__":
    main()
