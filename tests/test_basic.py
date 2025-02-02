from click.testing import CliRunner
from hiscanner.cli import cli
from hiscanner.__version__ import __version__

def test_cli_help():
    """Test if --help option works"""
    runner = CliRunner()
    result = runner.invoke(cli, ['--help'])
    assert result.exit_code == 0
    assert 'HiScanner:' in result.output

def test_cli_version():
    """Test if --version option works"""
    runner = CliRunner()
    result = runner.invoke(cli, ['--version'])
    assert result.exit_code == 0
    assert f'HiScanner, version {__version__}' in result.output
