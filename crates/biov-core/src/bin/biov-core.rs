//! Development binary for the first core slice; does not replace `biov` yet.
use biov_core::sequence::{normalize, reverse_complement, Kind};
use std::{env, process::ExitCode};
const USAGE: &str = "Usage: biov-core <normalize|reverse-complement> <dna|rna|protein> <sequence>";
fn run(args: &[String]) -> Result<String, String> {
    if args == ["--help"] || args == ["-h"] {
        return Ok(USAGE.to_owned());
    }
    if args.len() != 3 {
        return Err(USAGE.to_owned());
    }
    let kind: Kind = args[1]
        .parse()
        .map_err(|e: biov_core::sequence::SequenceError| e.to_string())?;
    match args[0].as_str() {
        "normalize" => normalize(&args[2], kind),
        "reverse-complement" => reverse_complement(&args[2], kind),
        _ => return Err(USAGE.to_owned()),
    }
    .map_err(|error| error.to_string())
}
fn main() -> ExitCode {
    match run(&env::args().skip(1).collect::<Vec<_>>()) {
        Ok(output) => {
            println!("{output}");
            ExitCode::SUCCESS
        }
        Err(message) => {
            eprintln!("{message}");
            ExitCode::from(2)
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dispatch() {
        assert_eq!(
            run(&["reverse-complement".into(), "rna".into(), "aug".into()]).unwrap(),
            "CAU"
        );
        assert!(run(&[]).is_err());
        assert!(run(&["normalize".into(), "dna".into(), "U".into()]).is_err());
    }
}
