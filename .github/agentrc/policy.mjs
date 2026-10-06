import { lstat, readFile, readdir } from "node:fs/promises";
import path from "node:path";

// Presence checks only. The existing CI runs the real build, tests and linters.
async function read(context, file) {
  const filename = path.join(context.repoPath, file);
  try {
    const stat = await lstat(filename);
    if (!stat.isFile() || stat.size > 1024 * 1024) return "";
    return await readFile(filename, "utf8");
  } catch (error) {
    if (error.code === "ENOENT" || error.code === "ENOTDIR") return "";
    throw error;
  }
}

function result(found, evidence) {
  return {
    status: found ? "pass" : "fail",
    evidence,
    reason: found ? undefined : "Missing Rust/Python configuration or command evidence.",
  };
}

function criterion(id, title, pillar, level, check) {
  const impact = ["lint-config", "build-script", "test-script", "lockfile"].includes(id)
    ? "high"
    : "medium";
  return {
    id,
    title,
    pillar,
    level,
    scope: "repo",
    impact,
    effort: "low",
    check,
  };
}

export default {
  name: "biov-rust-python",
  criteria: {
    add: [
      criterion(
        "lint-config",
        "Rust/Python lint commands configured",
        "style-validation",
        1,
        async (context) =>
          result(
            /^\[tool\.ruff\.lint\]$/mu.test(await read(context, "pyproject.toml")) &&
              /^\s*(?:run:\s*)?cargo clippy --workspace\b/mu.test(
                await read(context, ".github/workflows/ci.yml"),
              ),
            ["pyproject.toml", ".github/workflows/ci.yml"],
          ),
      ),
      criterion(
        "format-config",
        "Rust/Python formatter commands configured",
        "code-quality",
        2,
        async (context) =>
          result(
            /^\s*-?\s*id: ruff-format\s*$/mu.test(await read(context, ".pre-commit-config.yaml")) &&
              /^\s*(?:run:\s*)?cargo fmt --all -- --check\b/mu.test(
                await read(context, ".github/workflows/ci.yml"),
              ),
            [".pre-commit-config.yaml", ".github/workflows/ci.yml"],
          ),
      ),
      criterion(
        "typecheck-config",
        "Rust/Python type-check configuration present",
        "style-validation",
        2,
        async (context) =>
          result(
            /^\[workspace\]$/mu.test(await read(context, "Cargo.toml")) &&
              /^\[tool\.(?:mypy|pyright)\]$/mu.test(await read(context, "pyproject.toml")),
            ["Cargo.toml", "pyproject.toml"],
          ),
      ),
      criterion(
        "build-script",
        "Rust/Python build configuration present",
        "build-system",
        1,
        async (context) =>
          result(
            /^\[workspace\]$/mu.test(await read(context, "Cargo.toml")) &&
              /^\[build-system\]$/mu.test(await read(context, "pyproject.toml")),
            ["Cargo.toml", "pyproject.toml"],
          ),
      ),
      criterion(
        "test-script",
        "Rust/Python test commands and Python tests present",
        "testing",
        1,
        async (context) => {
          const ci = await read(context, ".github/workflows/ci.yml");
          const tests = await readdir(path.join(context.repoPath, "tests"), {
            withFileTypes: true,
          }).catch((error) => {
            if (error.code === "ENOENT" || error.code === "ENOTDIR") return [];
            throw error;
          });
          return result(
            /^\s*(?:run:\s*)?cargo test --workspace --locked\b/mu.test(ci) &&
              /^\s*(?:run:\s*)?uv run --locked pytest tests\//mu.test(ci) &&
              tests.some((entry) => entry.isFile() && /^test_.+\.py$/u.test(entry.name)),
            [".github/workflows/ci.yml", "tests/test_*.py"],
          );
        },
      ),
      criterion(
        "lockfile",
        "Cargo and Python lockfiles present",
        "dev-environment",
        1,
        async (context) =>
          result(
            Boolean(await read(context, "Cargo.lock")) && Boolean(await read(context, "uv.lock")),
            ["Cargo.lock", "uv.lock"],
          ),
      ),
    ],
  },
};
