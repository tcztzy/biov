# Codex and Claude Code plugins

This repository is a marketplace containing one `biov` plugin for each host.
Both use the same repository root, `skills/` tree, reference files and MCP
configuration. No copied skill tree or generated catalog is required.

## Prerequisites

Install the current BioV source as a command available on the host application's
PATH. From this checkout, with uv installed:

```sh
uv tool install .
biov --help
```

If BioV is already installed, use `uv tool install --force .` to replace it with
this checkout. Restart the host application if its PATH predates the installation.
The plugin starts `biov mcp`; it does not fetch a potentially older PyPI release
or install scientific environments during startup. Update the installed BioV
runtime when updating the plugin's source version.

Scientific software, model weights and data are configured separately as described
in the [execution guide](environments.md). A plugin installation makes guidance
discoverable; it does not establish that every listed program is deployed.

## Install from a local checkout

Run these commands from the repository root in the host you use.

Codex:

```sh
codex plugin marketplace add .
codex plugin add biov@biov
```

Claude Code:

```sh
claude plugin marketplace add .
claude plugin install biov@biov
```

Start a new session after installation. Skills are namespaced by the plugin;
ask the agent to use BioV's `scientific-software` or `biological-data` skill.
The host exposes skill descriptions and loads the relevant instruction and
reference files on demand. The MCP server supplies identifier tools, provider
records and [local managed Python analysis](analysis.md), including saved-record
queries and result access. The host's terminal tool runs native commands through
`biov exec`; managed analysis does not dispatch through SSH or LSF.

For local development, Claude Code also supports `claude --plugin-dir .`.
Do not enable this alongside an installed copy of the same plugin.

## Install from GitHub

After these files have been published to GitHub, replace `.` in the marketplace
add command with `tcztzy/biov`. Install the BioV runtime from the corresponding
checkout or release first. The marketplace files are:

- Codex: [`.agents/plugins/marketplace.json`](https://github.com/tcztzy/biov/blob/main/.agents/plugins/marketplace.json)
- Claude Code: [`.claude-plugin/marketplace.json`](https://github.com/tcztzy/biov/blob/main/.claude-plugin/marketplace.json)

The plugin is distributed through Git, separately from the Python wheel and
source distribution. Its root
includes all referenced skills, documentation and packaged resource descriptions.
External datasets, credentials, model weights and environment caches are not
plugin content. The `lab-protocols` skill links to publishers for protocol
discovery; third-party protocol full-text collections are not distributed.
Both plugin manifest versions must equal the package version in
`pyproject.toml`; `tests/test_plugin_distribution.py` enforces that, so bump
all three together when the contents change.

## Validate

```sh
uv run --locked pytest tests/test_plugin_distribution.py -q
claude plugin validate .claude-plugin/plugin.json --strict
claude plugin validate .claude-plugin/marketplace.json --strict
```

The repository test checks both marketplace roots, shared component paths, local
skill references and a real MCP startup with the installed BioV command. It does
not install the plugin into either host or run scientific workflows.

Host formats: [Codex plugin layout](https://github.com/openai/plugins/blob/main/README.md),
[Claude Code marketplace reference](https://code.claude.com/docs/en/plugins/marketplace-reference),
and [Claude Code plugin reference](https://code.claude.com/docs/en/plugins-reference).
