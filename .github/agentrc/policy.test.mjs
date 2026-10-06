import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { mkdtemp, mkdir, writeFile, rm, symlink } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";

const directory = path.dirname(fileURLToPath(import.meta.url));
const cli = path.join(directory, "node_modules/@microsoft/agentrc/dist/index.js");
const policy = path.join(directory, "policy.mjs");
const ids = [
  "lint-config",
  "format-config",
  "typecheck-config",
  "build-script",
  "test-script",
  "lockfile",
];

function scan(repo) {
  const output = execFileSync(
    process.execPath,
    [cli, "readiness", repo, "--json", "--policy", policy],
    {
      cwd: directory,
      env: { PATH: process.env.PATH, HOME: repo },
      encoding: "utf8",
      timeout: 10000,
    },
  );
  const result = JSON.parse(output);
  assert.equal(result.ok, true);
  assert.equal(result.status, "success");
  return Object.fromEntries(
    result.data.criteria.map((criterion) => [criterion.id, criterion.status]),
  );
}

test("pinned CLI detects Rust/Python configuration and missing evidence", async () => {
  const repo = await mkdtemp(path.join(os.tmpdir(), "biov-agentrc-"));
  try {
    await mkdir(path.join(repo, ".github/workflows"), { recursive: true });
    await mkdir(path.join(repo, "tests"));
    await mkdir(path.join(repo, "crates/example/src"), { recursive: true });
    await writeFile(path.join(repo, "Cargo.toml"), '[workspace]\nmembers = ["crates/example"]\n');
    await writeFile(
      path.join(repo, "crates/example/Cargo.toml"),
      '[package]\nname = "example"\nversion = "0.1.0"\n',
    );
    await writeFile(path.join(repo, "crates/example/src/lib.rs"), "#[test]\nfn works() {}\n");
    await writeFile(
      path.join(repo, "pyproject.toml"),
      '[build-system]\nrequires = ["setuptools"]\n[tool.ruff.lint]\n[tool.mypy]\n',
    );
    await writeFile(
      path.join(repo, ".pre-commit-config.yaml"),
      "repos:\n  - hooks:\n      - id: ruff-format\n",
    );
    const commands =
      "cargo clippy --workspace\ncargo fmt --all -- --check\ncargo test --workspace --locked\nuv run --locked pytest tests/ -q\n";
    await writeFile(path.join(repo, ".github/workflows/ci.yml"), commands);
    await writeFile(path.join(repo, "Cargo.lock"), "version = 4\n");
    await writeFile(path.join(repo, "uv.lock"), "version = 1\n");
    await writeFile(path.join(repo, "tests/test_example.py"), "def test_example(): pass\n");
    const passing = scan(repo);
    for (const id of ids) assert.equal(passing[id], "pass", id);

    await rm(path.join(repo, "uv.lock"));
    assert.equal(scan(repo).lockfile, "fail");
    await writeFile(path.join(repo, "uv.lock"), "version = 1\n");
    await rm(path.join(repo, "tests/test_example.py"));
    assert.equal(scan(repo)["test-script"], "fail");
    await writeFile(path.join(repo, "tests/test_example.py"), "def test_example(): pass\n");
    await writeFile(
      path.join(repo, ".github/workflows/ci.yml"),
      commands.replace("uv run --locked pytest tests/ -q\n", ""),
    );
    assert.equal(scan(repo)["test-script"], "fail");
    await writeFile(
      path.join(repo, ".github/workflows/ci.yml"),
      commands
        .split("\n")
        .map((line) => `# ${line}`)
        .join("\n"),
    );
    const commented = scan(repo);
    for (const id of ["lint-config", "format-config", "test-script"])
      assert.equal(commented[id], "fail", id);
    await writeFile(path.join(repo, ".github/workflows/ci.yml"), commands);
    await writeFile(
      path.join(repo, "pyproject.toml"),
      "# [build-system]\n# [tool.ruff.lint]\n# [tool.mypy]\n",
    );
    for (const id of ["build-script", "lint-config", "typecheck-config"])
      assert.equal(scan(repo)[id], "fail", id);
    await rm(path.join(repo, "Cargo.lock"));
    await symlink("uv.lock", path.join(repo, "Cargo.lock"));
    assert.equal(scan(repo).lockfile, "fail");
  } finally {
    await rm(repo, { recursive: true, force: true });
  }
});
