#![cfg(all(target_os = "linux", target_arch = "x86_64"))]
use std::{fs, process::Command};
#[path = "../../biov-tools/tests/support/mod.rs"]
mod support;

#[test]
fn native_cli_setup_inspection_literal_argv_and_status() {
    let dir = tempfile::tempdir().unwrap();
    let pixi = support::manager(dir.path());
    let root = dir.path().join("environments");
    let run = |words: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_biov-rs"))
            .args(words)
            .env("BIOV_ENVIRONMENT_ROOT", &root)
            .env("BIOV_PIXI_BIN", &pixi)
            .output()
            .unwrap()
    };
    let before = run(&["tools", "inspect", "goatools"]);
    assert!(before.status.success());
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&before.stdout).unwrap()["status"],
        "unavailable"
    );
    assert!(!root.exists());
    let unavailable_cwd = dir.path().join("unavailable");
    let bad_cwd = run(&[
        "tools",
        "exec",
        "--cwd",
        unavailable_cwd.to_str().unwrap(),
        "goatools",
    ]);
    assert_eq!(bad_cwd.status.code(), Some(2));
    assert!(!root.exists());
    let missing = run(&["tools", "exec", "--no-install", "goatools", "--help"]);
    assert_eq!(missing.status.code(), Some(2));
    assert!(!root.exists());
    let first = run(&[
        "tools",
        "exec",
        "--cwd",
        dir.path().to_str().unwrap(),
        "goatools",
        "--",
        "",
        "O'Connor",
        "$(touch must-not-exist)",
        "--environment-root",
        "native-option",
    ]);
    assert_eq!(
        first.status.code(),
        Some(37),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert_eq!(first.stdout, b"native stdout\n");
    assert!(String::from_utf8(first.stderr)
        .unwrap()
        .contains("native stderr"));
    let native: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.path().join("native.json")).unwrap()).unwrap();
    assert_eq!(native["argv"][1], "");
    assert_eq!(native["argv"][2], "O'Connor");
    assert_eq!(native["argv"][3], "$(touch must-not-exist)");
    assert_eq!(native["argv"][4], "--environment-root");
    let second = run(&["tools", "exec", "--no-install", "goatools", "--help"]);
    assert_eq!(second.status.code(), Some(37));
    let log = fs::read_to_string(dir.path().join("argv.jsonl")).unwrap();
    assert_eq!(
        log.lines()
            .filter(|line| line.starts_with("[\"install\""))
            .count(),
        1
    );
    assert!(!dir.path().join("must-not-exist").exists());
}

#[test]
fn parser_errors_and_help_have_no_execution_side_effects() {
    // Synchronize fixture materialization before either test starts a child.
    support::prepare();
    for args in [
        vec!["tools"],
        vec!["tools", "unknown"],
        vec!["tools", "exec", "--pixi"],
        vec!["tools", "exec", "--no-install", "--no-install", "samtools"],
        vec!["tools", "setup", "goatools", "--help"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_biov-rs"))
            .args(args)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
    }
    let output = Command::new(env!("CARGO_BIN_EXE_biov-rs"))
        .args(["tools", "--help"])
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8(output.stderr).unwrap().contains("setup"));
}

#[test]
fn direct_wrapper_interrupt_forwards_reaps_and_retains_workspace_lock() {
    use fs4::fs_std::FileExt;
    use std::{
        os::unix::fs::PermissionsExt,
        process::{Child, Stdio},
        time::{Duration, Instant},
    };
    struct Cleanup {
        child: Child,
        pids: std::path::PathBuf,
        done: bool,
    }
    impl Drop for Cleanup {
        fn drop(&mut self) {
            if !self.done {
                if let Ok(bytes) = fs::read(&self.pids) {
                    let ids: Vec<u32> = serde_json::from_slice(&bytes).unwrap();
                    let _ = Command::new("kill")
                        .args(["-KILL", "--", &format!("-{}", ids[0])])
                        .status();
                }
                let _ = self.child.kill();
                let _ = self.child.wait();
            }
        }
    }
    for (signal, code) in [("TERM", 143), ("INT", 130), ("HUP", 129)] {
        let dir = tempfile::tempdir().unwrap();
        let pixi = support::manager(dir.path());
        let root = dir.path().join("environments");
        let output = Command::new(env!("CARGO_BIN_EXE_biov-rs"))
            .args(["tools", "setup", "samtools"])
            .env("BIOV_ENVIRONMENT_ROOT", &root)
            .env("BIOV_PIXI_BIN", &pixi)
            .output()
            .unwrap();
        assert!(output.status.success());
        let receipt: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
        let executable =
            std::path::Path::new(receipt["prefix"].as_str().unwrap()).join("bin/samtools");
        fs::write(&executable, r#"#!/usr/bin/env python3
import json, os, pathlib, signal, subprocess, sys, time
home = pathlib.Path.cwd()
child = subprocess.Popen([sys.executable, '-c', "import pathlib,time; pathlib.Path('descendant-ready').write_text('ready'); time.sleep(60)"])
def interrupted(number, frame):
    child.wait(timeout=3)
    (home / 'interrupted').write_text(str(number))
    time.sleep(0.4)
    sys.exit(128 + number)
for number in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
    signal.signal(number, interrupted)
while not (home / 'descendant-ready').exists():
    time.sleep(0.005)
(home / 'pids.json').write_text(json.dumps([os.getpid(), child.pid]))
while True:
    time.sleep(1)
"#).unwrap();
        fs::set_permissions(&executable, fs::Permissions::from_mode(0o755)).unwrap();
        let child = Command::new(env!("CARGO_BIN_EXE_biov-rs"))
            .args([
                "tools",
                "exec",
                "--no-install",
                "--cwd",
                dir.path().to_str().unwrap(),
                "samtools",
            ])
            .env("BIOV_ENVIRONMENT_ROOT", &root)
            .env("BIOV_PIXI_BIN", &pixi)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let mut cleanup = Cleanup {
            child,
            pids: dir.path().join("pids.json"),
            done: false,
        };
        let deadline = Instant::now() + Duration::from_secs(8);
        while !cleanup.pids.exists() {
            assert!(Instant::now() < deadline, "native child did not start");
            std::thread::sleep(Duration::from_millis(5));
        }
        let ids: Vec<u32> = serde_json::from_slice(&fs::read(&cleanup.pids).unwrap()).unwrap();
        assert!(Command::new("kill")
            .args(["-s", signal, &cleanup.child.id().to_string()])
            .status()
            .unwrap()
            .success());
        while !dir.path().join("interrupted").exists() {
            assert!(
                cleanup.child.try_wait().unwrap().is_none(),
                "wrapper exited before forwarding/reaping"
            );
            assert!(
                Instant::now() < deadline,
                "interrupt did not reach child and descendant"
            );
            std::thread::sleep(Duration::from_millis(5));
        }
        let lock = fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(
                root.join("native-tool-locks")
                    .join(format!("{}.lock", biov_tools::workspace_identity())),
            )
            .unwrap();
        assert!(
            !FileExt::try_lock_exclusive(&lock).unwrap(),
            "workspace lock released before native child exited"
        );
        let status = loop {
            if let Some(status) = cleanup.child.try_wait().unwrap() {
                break status;
            }
            assert!(Instant::now() < deadline, "wrapper did not reap and finish");
            std::thread::sleep(Duration::from_millis(5));
        };
        assert_eq!(status.code(), Some(code));
        assert!(FileExt::try_lock_exclusive(&lock).unwrap());
        for pid in ids {
            assert!(
                !std::path::Path::new(&format!("/proc/{pid}")).exists(),
                "native process survived wrapper"
            );
        }
        cleanup.done = true;
    }
}

#[test]
fn interactive_stdin_and_interrupt_delivery_are_preserved() {
    let dir = tempfile::tempdir().unwrap();
    let pixi = support::manager(dir.path());
    let output = Command::new("python3")
        .args(["-c", r#"
import fcntl, json, os, pathlib, pty, signal, subprocess, sys, termios, time
binary, home, manager = sys.argv[1:]
home = pathlib.Path(home)
env = dict(os.environ, BIOV_ENVIRONMENT_ROOT=str(home / 'envs'), BIOV_PIXI_BIN=manager)
receipt = json.loads(subprocess.run([binary, 'tools', 'setup', 'samtools'], env=env, capture_output=True, check=True).stderr)
native = pathlib.Path(receipt['prefix']) / 'bin/samtools'
native.write_text('''#!/usr/bin/env python3
import json, os, pathlib, signal, sys, time
home = pathlib.Path.cwd()
count = 0
def interrupted(number, frame):
    global count
    count += 1
    (home / 'count').write_text(str(count))
signal.signal(signal.SIGINT, interrupted)
line = sys.stdin.readline()
(home / 'ready').write_text(json.dumps([os.getpid(), os.getpgrp(), line]))
while count == 0:
    time.sleep(.005)
time.sleep(.2)
sys.exit(33)
''')
native.chmod(0o755)
for delivery in ('direct', 'terminal'):
    for name in ('ready', 'count'):
        (home / name).unlink(missing_ok=True)
    master, slave = pty.openpty()
    def terminal_session():
        os.setsid()
        fcntl.ioctl(0, termios.TIOCSCTTY, 0)
    wrapper = subprocess.Popen([binary, 'tools', 'exec', '--no-install', '--cwd', str(home), 'samtools'], env=env, stdin=slave, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, preexec_fn=terminal_session)
    os.close(slave)
    native_pid = None
    try:
        os.write(master, b'input from terminal\n')
        deadline = time.monotonic() + 8
        while not (home / 'ready').exists():
            assert time.monotonic() < deadline, 'native terminal input stalled'
            time.sleep(.005)
        native_pid, group, line = json.loads((home / 'ready').read_text())
        assert group == wrapper.pid, 'interactive foreground group changed'
        assert line == 'input from terminal\n'
        if delivery == 'direct':
            os.kill(wrapper.pid, signal.SIGINT)
        else:
            os.write(master, b'\x03')
        assert wrapper.wait(timeout=8) == 33
        assert (home / 'count').read_text() == '1', 'interrupt delivered twice'
        assert not pathlib.Path(f'/proc/{native_pid}').exists(), 'native child not reaped'
    finally:
        if wrapper.poll() is None:
            wrapper.kill()
            wrapper.wait()
        if native_pid is not None and pathlib.Path(f'/proc/{native_pid}').exists():
            os.kill(native_pid, signal.SIGKILL)
        os.close(master)
"#, env!("CARGO_BIN_EXE_biov-rs"), dir.path().to_str().unwrap(), pixi.to_str().unwrap()])
        .output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
