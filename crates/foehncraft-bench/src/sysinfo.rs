//! Host description for reports.

/// CPU model name, best effort.
pub fn cpu_model() -> String {
    #[cfg(target_os = "linux")]
    if let Some(v) = std::fs::read_to_string("/proc/cpuinfo").ok().and_then(|s| {
        s.lines()
            .find(|l| l.starts_with("model name"))
            .and_then(|l| l.split_once(':'))
            .map(|(_, v)| v.trim().to_string())
    }) {
        return v;
    }
    #[cfg(target_os = "macos")]
    if let Ok(o) = std::process::Command::new("sysctl")
        .args(["-n", "machdep.cpu.brand_string"])
        .output()
    {
        let s = String::from_utf8_lossy(&o.stdout).trim().to_string();
        if !s.is_empty() {
            return s;
        }
    }
    #[cfg(target_os = "windows")]
    if let Ok(s) = std::env::var("PROCESSOR_IDENTIFIER") {
        return s;
    }
    format!("{} CPU", std::env::consts::ARCH)
}

/// Operating system and architecture.
pub fn platform() -> String {
    format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH)
}
