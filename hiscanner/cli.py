import os
import sys
import subprocess
from pathlib import Path

def main():
    # Determine the current operating system
    if sys.platform.startswith('linux'):
        binary_name = 'mbicseq'
        binary_path = Path(__file__).parent / 'cli_tools/linux/mbicseq'
    elif sys.platform.startswith('darwin'):
        binary_name = 'mbicseq'
        binary_path = Path(__file__).parent / 'cli_tools/osx/mbicseq'
    elif sys.platform.startswith('win'):
        binary_name = 'mbicseq.exe'
        binary_path = Path(__file__).parent / 'cli_tools/win/mbicseq.exe'
    else:
        raise OSError("Unsupported operating system")

    # Ensure the binary is executable
    binary_path.chmod(0o755)

    # Execute the binary with any arguments passed to the wrapper
    result = subprocess.run([binary_path] + sys.argv[1:], check=True)

    return result.returncode

if __name__ == "__main__":
    main()
