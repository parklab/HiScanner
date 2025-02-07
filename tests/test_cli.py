import pytest
from click.testing import CliRunner
from pathlib import Path
from hiscanner.cli import cli
from hiscanner.__version__ import __version__

@pytest.fixture
def runner():
    return CliRunner()

@pytest.fixture
def temp_dir(tmp_path):
    return tmp_path

def test_cli_version(runner):
    result = runner.invoke(cli, ['--version'])
    assert result.exit_code == 0
    assert f'HiScanner, version {__version__}' in result.output

def test_cli_help(runner):
    result = runner.invoke(cli, ['--help'])
    assert result.exit_code == 0
    assert 'HiScanner:' in result.output
    assert 'Commands:' in result.output

def test_init_command(runner, temp_dir):
    result = runner.invoke(cli, ['init', '--output', str(temp_dir)])
    assert result.exit_code == 0
    assert (temp_dir / 'config.yaml').exists()
    assert (temp_dir / 'cluster.yaml').exists()

def test_run_command_validates_steps(runner):
    result = runner.invoke(cli, ['run', '--step', 'invalid_step'])
    assert result.exit_code != 0
    assert "Error" in result.output

def test_run_command_requires_config(runner):
    result = runner.invoke(cli, ['run', '--step', 'all'])
    assert result.exit_code != 0