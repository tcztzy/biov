//! Keep the operation lock alive while forwarding interrupts and reaping the tool.
use std::process::{Command, ExitStatus};

#[cfg(target_os = "linux")]
pub(super) fn run(command: &mut Command) -> Result<ExitStatus, String> {
    use signal_hook::{
        consts::signal::{SIGHUP, SIGINT, SIGTERM},
        low_level::unregister,
    };
    use std::{
        io::IsTerminal,
        os::unix::process::CommandExt,
        sync::{
            atomic::{AtomicBool, Ordering},
            Arc, Mutex, OnceLock,
        },
        time::Duration,
    };

    // Unregistering signal-hook actions does not restore the OS default handler.
    // Keep one conditional default action per signal, enabled outside execution.
    struct SignalState {
        idle: Arc<AtomicBool>,
        active: Mutex<usize>,
    }
    static STATE: OnceLock<Result<SignalState, String>> = OnceLock::new();
    let state = STATE
        .get_or_init(|| {
            let idle = Arc::new(AtomicBool::new(true));
            for signal in [SIGINT, SIGTERM, SIGHUP] {
                signal_hook::flag::register_conditional_default(signal, Arc::clone(&idle))
                    .map_err(|e| {
                        format!("cannot register native execution default handler: {e}")
                    })?;
            }
            Ok(SignalState {
                idle,
                active: Mutex::new(0),
            })
        })
        .as_ref()
        .map_err(Clone::clone)?;
    struct Active<'a>(&'a SignalState);
    impl Drop for Active<'_> {
        fn drop(&mut self) {
            let mut active = self.0.active.lock().unwrap_or_else(|e| e.into_inner());
            *active -= 1;
            self.0.idle.store(*active == 0, Ordering::SeqCst);
        }
    }

    struct Handlers(Vec<signal_hook::SigId>);
    impl Drop for Handlers {
        fn drop(&mut self) {
            for id in self.0.drain(..) {
                unregister(id);
            }
        }
    }
    let mut handlers = Handlers(Vec::new());
    let mut pending = Vec::new();
    let separate_group = !std::io::stdin().is_terminal();
    let child_started = Arc::new(AtomicBool::new(false));
    for signal in [SIGINT, SIGTERM, SIGHUP] {
        let flag = Arc::new(AtomicBool::new(false));
        let pending_flag = Arc::clone(&flag);
        let established_child = Arc::clone(&child_started);
        handlers.0.push(
            // SAFETY: the handler only reads siginfo and stores an atomic flag.
            // Terminal-generated signals already reach an interactive child in
            // the shared foreground group; do not deliver them a second time.
            unsafe {
                signal_hook_registry::register_sigaction(signal, move |info| {
                    if separate_group
                        || info.si_code != libc::SI_KERNEL
                        || !established_child.load(Ordering::SeqCst)
                    {
                        pending_flag.store(true, Ordering::SeqCst);
                    }
                })
            }
            .map_err(|e| format!("cannot register native execution signal handler: {e}"))?,
        );
        pending.push((signal, flag));
    }
    {
        let mut active = state.active.lock().unwrap_or_else(|e| e.into_inner());
        *active += 1;
        state.idle.store(false, Ordering::SeqCst);
    }
    let _active = Active(state);
    // A separate group lets direct wrapper interrupts reach pipeline descendants.
    // Interactive children retain the terminal's foreground group so stdin and
    // the shell's normal job control keep working.
    if separate_group {
        command.process_group(0);
    }
    let mut child = command
        .spawn()
        .map_err(|e| format!("native executable launch failed: {e}"))?;
    child_started.store(true, Ordering::SeqCst);
    let target = if separate_group {
        -(child.id() as i32)
    } else {
        child.id() as i32
    };
    loop {
        // Check exit first: never send a pending signal to an already reaped PID.
        if let Some(status) = child
            .try_wait()
            .map_err(|e| format!("native executable wait failed: {e}"))?
        {
            return Ok(status);
        }
        for (signal, flag) in &pending {
            if flag.swap(false, Ordering::SeqCst) {
                // SAFETY: target is the live child PID or its process group,
                // created by this invocation; signal is a fixed Unix constant.
                // ESRCH is harmless when the child exits between polling and kill.
                unsafe {
                    libc::kill(target, *signal);
                }
            }
        }
        std::thread::sleep(Duration::from_millis(10));
    }
}

#[cfg(not(target_os = "linux"))]
pub(super) fn run(command: &mut Command) -> Result<ExitStatus, String> {
    command
        .status()
        .map_err(|e| format!("native executable launch failed: {e}"))
}
