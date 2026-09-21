"""CLI registration. Cheap, but it catches a renamed command or a broken import
before a Nextflow process fails mid-run."""
import re

from typer.testing import CliRunner

from kmeans_py.cli import app

runner = CliRunner()

EXPECTED_COMMANDS = ["cluster", "choose-k", "aggregate", "project-profiles"]

_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def plain(result) -> str:
    """Output with colour removed. Rich splits option-like tokens across styled
    spans, so `--jaccard` is not a literal substring of the coloured text. Colour is
    off locally and forced on in CI, so asserting on the raw output passes on a laptop
    and fails only in CI."""
    return _ANSI.sub("", result.output)


def test_help_exits_zero():
    assert runner.invoke(app, ["--help"]).exit_code == 0


def test_all_commands_registered():
    out = plain(runner.invoke(app, ["--help"]))
    for command in EXPECTED_COMMANDS:
        assert command in out


def test_subcommand_help():
    for command in EXPECTED_COMMANDS:
        assert runner.invoke(app, [command, "--help"]).exit_code == 0, command


def test_k_is_required():
    result = runner.invoke(app, ["cluster", "--project", "P", "--matrix", "m.csv"])
    assert result.exit_code != 0
    assert "--k" in plain(result)


def test_aggregate_needs_stats():
    result = runner.invoke(app, ["aggregate", "--labels", "l.txt", "--k", "3"])
    assert result.exit_code != 0


def test_no_p_value_is_offered():
    """The package must not grow an --alpha or a p-value flag back by accident: a
    flat partition has no AU p-value, and the README explains why."""
    for command in EXPECTED_COMMANDS:
        out = plain(runner.invoke(app, [command, "--help"]))
        assert "--alpha" not in out
        assert "--jaccard" not in out
